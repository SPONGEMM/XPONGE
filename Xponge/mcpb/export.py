"""Export helpers for MCPB workflows."""

from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

from ._compat import atom_element, atom_to_residue_map, atom_type_name, save_pdb, save_sponge_input


def audit_sponge_ready(molecule, *, prefix: str = "mcpb_audit") -> tuple[bool, str | None]:
    try:
        with tempfile.TemporaryDirectory(prefix="xponge_mcpb_audit_") as tmpdir:
            save_sponge_input(molecule, prefix, dirname=tmpdir)
        return True, None
    except Exception as exc:
        return False, f"{type(exc).__name__}: {exc}"


def save_pdb_with_connect(molecule, filename) -> str:
    return save_pdb(molecule, filename)


def _charge_override_records(result) -> list[dict[str, object]]:
    mapping = atom_to_residue_map(result.molecule)
    records: list[dict[str, object]] = []
    for atom_id in result.updated_charge_atoms:
        atom = result.molecule.atoms[int(atom_id)]
        residue = result.molecule.residues[mapping[int(atom_id)]]
        records.append(
            {
                "atom_id": int(atom_id),
                "residue_id": int(mapping[int(atom_id)]),
                "residue_name": str(residue.name),
                "atom_name": str(atom.name),
                "element": atom_element(atom),
                "type": atom_type_name(atom),
                "charge": float(atom.charge),
            }
        )
    return records


def _json_ready(value):
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _json_ready(val) for key, val in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_ready(item) for item in value]
    return repr(value)


def _materialize_optional_frcmod(result, directory: Path, prefix: str) -> str | None:
    if not result.frcmod_path:
        return None
    source = Path(result.frcmod_path)
    target = directory / f"{prefix}.frcmod"
    if source.resolve() != target.resolve():
        shutil.copyfile(source, target)
    else:
        target = source
    return str(target)


def write_mcpb_artifacts(
    result,
    directory,
    *,
    prefix: str = "mcpb",
    write_sponge_input_files: bool = True,
    write_local_models: bool = True,
) -> dict[str, object]:
    outdir = Path(directory)
    outdir.mkdir(parents=True, exist_ok=True)

    pdb_path = save_pdb_with_connect(result.molecule, outdir / f"{prefix}.pdb")
    frcmod_path = _materialize_optional_frcmod(result, outdir, prefix)

    charge_override_records = _charge_override_records(result)
    charge_override_path = outdir / f"{prefix}.charge_overrides.json"
    charge_override_path.write_text(
        json.dumps(charge_override_records, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    small_model_pdb = None
    large_model_pdb = None
    if write_local_models:
        if result.small_model is not None:
            small_model_pdb = save_pdb_with_connect(result.small_model.molecule, outdir / f"{prefix}.small_model.pdb")
        if result.large_model is not None:
            large_model_pdb = save_pdb_with_connect(result.large_model.molecule, outdir / f"{prefix}.large_model.pdb")

    sponge_files: list[str] = []
    if write_sponge_input_files:
        save_sponge_input(result.molecule, prefix, dirname=str(outdir))
        sponge_files = sorted(path.name for path in outdir.glob(f"{prefix}_*.txt"))

    manifest = {
        "prefix": prefix,
        "pdb_path": str(pdb_path),
        "frcmod_path": frcmod_path,
        "charge_override_path": str(charge_override_path),
        "small_model_pdb": small_model_pdb,
        "large_model_pdb": large_model_pdb,
        "sponge_files": sponge_files,
        "sponge_ready": bool(result.sponge_ready),
        "pending_requirements": list(result.pending_requirements),
        "warnings": list(result.warnings),
        "metadata": _json_ready(dict(result.metadata)),
        "connect_records": [list(record) for record in result.connect_records],
        "updated_charge_atoms": list(result.updated_charge_atoms),
        "registered_metal_templates": list(result.registered_metal_templates),
    }
    manifest_path = outdir / f"{prefix}.mcpb.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    manifest["manifest_path"] = str(manifest_path)
    return manifest
