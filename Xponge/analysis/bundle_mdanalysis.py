"""MDAnalysis topology and trajectory readers for CIF and SPONGE H5MD."""

from __future__ import annotations

import hashlib
import json
import os
import re
import warnings
from typing import Any

import h5py
import numpy as np

import MDAnalysis as mda
from MDAnalysis.coordinates import H5MD, base
from MDAnalysis.core import topologyattrs
from MDAnalysis.core.topology import Topology
from MDAnalysis.topology.base import TopologyReaderBase
from MDAnalysis.units import get_conversion_factor
from MDAnalysis.lib.util import openany

from ..helper import guess_element_from_mass
from ..io_bundle.errors import (
    AmbiguousH5MDLayoutError,
    BundleSchemaError,
    BundleTopologyError,
    BundleTrajectoryError,
    BundleUnitError,
    IncompleteBundleError,
    UnverifiedBundlePairError,
)


_INPUT_SCHEMA_VERSION = "sponge.input.v2"
_OUTPUT_SCHEMA_VERSION = "sponge.output.v2"
_TOPOLOGY_SCHEMA = "sponge.topology.h5"
_TRAJECTORY_SCHEMA = "sponge.output.h5md"
_LEGACY_WALKER_RE = re.compile(r"^trajectory(?P<index>\d+)$")


def _mokda_mapping_digest(mapping_rows: list[dict[str, Any]]) -> str:
    """Return Mokda's canonical digest for ordered simulation atom mappings."""

    fields = (
        "simulation_index",
        "external_id",
        "canonical_atom_id",
        "simulation_residue_index",
        "simulation_residue_id",
    )
    payload = {
        "atom_mapping": [
            {field: row[field] for field in fields} for row in mapping_rows
        ]
    }
    canonical = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def _decode_text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("utf-8")
    if isinstance(value, np.ndarray) and value.shape == ():
        return _decode_text(value.item())
    return str(value)


def _read_text(handle: h5py.File, path: str) -> str:
    if path not in handle:
        raise BundleSchemaError(f"{handle.filename}:{path} is missing")
    return _decode_text(handle[path][()])


def _optional_text(handle: h5py.File, path: str) -> str | None:
    return _decode_text(handle[path][()]) if path in handle else None


def _read_string_array(handle: h5py.File, path: str) -> list[str]:
    if path not in handle:
        return []
    return [_decode_text(value) for value in np.asarray(handle[path][...]).reshape(-1)]


def _open_hint_target(thing: Any):
    if isinstance(thing, h5py.File):
        return thing, False
    if not isinstance(thing, (str, os.PathLike)):
        return None, False
    try:
        return h5py.File(os.fspath(thing), "r"), True
    except (OSError, TypeError, ValueError):
        return None, False


def _validate_schema(
    handle: h5py.File,
    *,
    name_path: str,
    version_path: str,
    expected_name: str,
    strict: bool,
) -> None:
    if name_path not in handle or version_path not in handle:
        if strict:
            raise BundleSchemaError(
                f"{handle.filename} is missing schema metadata at {name_path} or {version_path}"
            )
        return
    name = _read_text(handle, name_path)
    version = _read_text(handle, version_path)
    if name != expected_name:
        raise BundleSchemaError(
            f"{handle.filename}:{name_path} is {name!r}, expected {expected_name!r}"
        )
    expected_version = (
        _OUTPUT_SCHEMA_VERSION if expected_name == _TRAJECTORY_SCHEMA else _INPUT_SCHEMA_VERSION
    )
    if version != expected_version:
        raise BundleSchemaError(
            f"{handle.filename}:{version_path} is {version!r}, expected {expected_version!r}"
        )


def _array(handle: h5py.File, path: str, *, dtype=None) -> np.ndarray:
    if path not in handle:
        raise BundleTopologyError(f"{handle.filename}:{path} is missing")
    return np.asarray(handle[path][...], dtype=dtype)


def _optional_array(handle: h5py.File, path: str, *, dtype=None) -> np.ndarray | None:
    if path not in handle:
        return None
    return np.asarray(handle[path][...], dtype=dtype)


def _require_vector(values: np.ndarray, size: int, *, filename: str, path: str) -> np.ndarray:
    if values.shape != (size,):
        raise BundleTopologyError(
            f"{filename}:{path} has shape {values.shape}, expected ({size},)"
        )
    return values


def _charge_values(dataset: h5py.Dataset, *, strict: bool) -> np.ndarray:
    values = np.asarray(dataset[...], dtype=np.float64)
    unit_value = dataset.attrs.get("unit")
    if unit_value is None:
        if strict:
            warnings.warn(
                f"{dataset.file.filename}:{dataset.name} has no unit; interpreting schema v1 values as Amber charge",
                RuntimeWarning,
                stacklevel=3,
            )
        unit = "Amber"
    else:
        unit = _decode_text(unit_value)
    aliases = {"electron": "e", "elementary_charge": "e", "SPONGE": "Amber"}
    unit = aliases.get(unit, unit)
    try:
        factor = get_conversion_factor("charge", unit, "e")
    except (KeyError, ValueError) as exc:
        raise BundleUnitError(
            f"{dataset.file.filename}:{dataset.name} has unsupported charge unit {unit!r}"
        ) from exc
    return values * factor


def _deduplicate_connections(values: np.ndarray, width: int, *, normalize_reverse: bool) -> list[tuple[int, ...]]:
    if values.size == 0:
        return []
    if values.ndim != 2 or values.shape[1] != width:
        raise BundleTopologyError(
            f"connectivity table has shape {values.shape}, expected (*, {width})"
        )
    result: list[tuple[int, ...]] = []
    seen: set[tuple[int, ...]] = set()
    for row in values:
        item = tuple(int(index) for index in row)
        key = min(item, item[::-1]) if normalize_reverse else item
        if key not in seen:
            seen.add(key)
            result.append(item)
    return result


class CIFTopologyParser(TopologyReaderBase):
    """Read mmCIF ``_atom_site`` metadata into an MDAnalysis topology.

    The parser follows Xponge's mmCIF conventions: author atom/residue/chain
    identifiers take precedence, label identifiers fill missing author fields,
    and alternate location ``A`` is selected by default.  Explicit PDB segment
    IDs are retained as MDAnalysis ``segids``; otherwise ``label_asym_id`` is
    used as the segment ID while the author chain remains available as
    ``chainID``.
    """

    format = "XPONGE_CIF"

    @staticmethod
    def _format_hint(thing: Any) -> bool:
        if not isinstance(thing, (str, os.PathLike)):
            return False
        filename = os.fspath(thing)
        lower_filename = filename.lower()
        if lower_filename.endswith(
            (".cif", ".mmcif", ".mcif", ".cif.gz", ".mmcif.gz", ".mcif.gz")
        ):
            return True
        try:
            with openany(filename, "rt") as handle:
                return "_atom_site." in handle.read(65536).lower()
        except (OSError, TypeError, ValueError):
            return False

    def parse(self, **kwargs) -> Topology:
        # Reuse Xponge's CIF tokenizer and missing-value conventions so quoted
        # values, loop rows, and custom Mokda categories are interpreted alike.
        from ..load import (
            _mmcif_clean,
            _mmcif_first,
            _mmcif_parse,
            _mmcif_parse_float,
            _mmcif_parse_int,
            _mmcif_rows,
        )

        with openany(self.filename, "rt") as handle:
            data = _mmcif_parse(handle)
        atom_rows = _mmcif_rows(data, "atom_site")
        if not atom_rows:
            raise ValueError("mmCIF file does not contain _atom_site records")

        def first(row: dict[str, str], *tags: str) -> str:
            return _mmcif_first(row, list(tags))

        def parse_formal_charge(value: str) -> int | None:
            cleaned = _mmcif_clean(value)
            charge = _mmcif_parse_int(cleaned)
            if charge is None and len(cleaned) > 1 and cleaned[-1] in "+-":
                magnitude = _mmcif_parse_int(cleaned[:-1])
                if magnitude is not None:
                    charge = magnitude if cleaned[-1] == "+" else -magnitude
            return charge

        model_values = {
            first(row, "_atom_site.pdbx_pdb_model_num") or "1" for row in atom_rows
        }
        model_id = kwargs.get("model_id")
        if model_id is None and len(model_values) > 1:
            raise ValueError("mmCIF contains multiple models; pass model_id explicitly")
        selected_model = str(model_id) if model_id is not None else next(iter(model_values))
        if selected_model not in model_values:
            raise ValueError(f"mmCIF does not contain model {selected_model!r}")

        altloc = kwargs.get("altloc", "A")
        selected_rows = []
        for row in atom_rows:
            row_model = first(row, "_atom_site.pdbx_pdb_model_num") or "1"
            if row_model != selected_model:
                continue
            row_altloc = first(row, "_atom_site.label_alt_id")
            if altloc is not None and row_altloc and row_altloc != str(altloc):
                continue
            selected_rows.append(row)
        if not selected_rows:
            raise ValueError("mmCIF model/alternate-location selection produced no atoms")

        atom_ids: list[int] = []
        atom_names: list[str] = []
        atom_types: list[str] = []
        elements: list[str] = []
        residue_names: list[str] = []
        residue_ids: list[int] = []
        residue_numbers: list[int] = []
        residue_index_by_atom: list[int] = []
        residue_segment_index: list[int] = []
        segment_ids: list[str] = []
        chain_ids: list[str] = []
        insertion_codes: list[str] = []
        alternate_locations: list[str] = []
        occupancies: list[float | None] = []
        temperature_factors: list[float | None] = []
        formal_charges: list[int | None] = []
        record_types: list[str] = []
        residue_components: list[str] = []
        residue_insertion_codes: list[str] = []
        residue_atom_names: list[dict[str, int]] = []
        atom_by_site_id: dict[str, int] = {}
        atom_by_label_key: dict[tuple[str, ...], int] = {}
        atom_by_auth_key: dict[tuple[str, ...], int] = {}
        segment_index_by_id: dict[str, int] = {}
        residue_index_by_key: dict[tuple[str, ...], int] = {}
        atom_id_set: set[int] = set()

        def identity_key(asym_id, seq_id, comp_id, atom_id, insertion_code):
            return (
                _mmcif_clean(asym_id),
                _mmcif_clean(seq_id),
                _mmcif_clean(comp_id).upper(),
                _mmcif_clean(atom_id).upper(),
                _mmcif_clean(insertion_code),
            )

        for row in selected_rows:
            auth_asym = first(row, "_atom_site.auth_asym_id")
            label_asym = first(row, "_atom_site.label_asym_id")
            chain_id = auth_asym or label_asym
            segid = (
                first(
                    row,
                    "_atom_site.pdbx_pdb_segment_id",
                    "_atom_site.seg_id",
                    "_atom_site.segid",
                )
                or label_asym
                or auth_asym
                or "SYSTEM"
            )
            if segid not in segment_index_by_id:
                segment_index_by_id[segid] = len(segment_ids)
                segment_ids.append(segid)
            segment_index = segment_index_by_id[segid]

            auth_seq = first(row, "_atom_site.auth_seq_id")
            label_seq = first(row, "_atom_site.label_seq_id")
            auth_resid = _mmcif_parse_int(auth_seq)
            label_resid = _mmcif_parse_int(label_seq)
            if auth_resid is None:
                auth_resid = label_resid
            if label_resid is None:
                label_resid = auth_resid
            auth_comp = first(row, "_atom_site.auth_comp_id")
            label_comp = first(row, "_atom_site.label_comp_id")
            residue_name = auth_comp or label_comp or "SYSTEM"
            insertion_code = first(row, "_atom_site.pdbx_pdb_ins_code")
            residue_key = (
                segid,
                chain_id,
                label_asym,
                auth_seq,
                label_seq,
                insertion_code,
                residue_name,
            )
            if residue_key not in residue_index_by_key:
                residue_index = len(residue_names)
                residue_index_by_key[residue_key] = residue_index
                resid = auth_resid if auth_resid is not None else residue_index + 1
                resnum = label_resid if label_resid is not None else resid
                residue_names.append(residue_name)
                residue_ids.append(resid)
                residue_numbers.append(resnum)
                residue_insertion_codes.append(insertion_code)
                residue_segment_index.append(segment_index)
                residue_components.append(label_comp or auth_comp or residue_name)
                residue_atom_names.append({})
            residue_index = residue_index_by_key[residue_key]
            residue_index_by_atom.append(residue_index)

            atom_index = len(atom_names)
            raw_atom_id = first(row, "_atom_site.id")
            parsed_atom_id = _mmcif_parse_int(raw_atom_id)
            if parsed_atom_id is None or parsed_atom_id in atom_id_set:
                parsed_atom_id = atom_index + 1
                while parsed_atom_id in atom_id_set:
                    parsed_atom_id += 1
            atom_id_set.add(parsed_atom_id)
            atom_ids.append(parsed_atom_id)

            auth_atom = first(row, "_atom_site.auth_atom_id")
            label_atom = first(row, "_atom_site.label_atom_id")
            atom_name = auth_atom or label_atom or f"A{atom_index + 1}"
            atom_names.append(atom_name)
            type_symbol = first(row, "_atom_site.type_symbol")
            element = type_symbol[:1].upper() + type_symbol[1:].lower() if type_symbol else ""
            elements.append(element)
            atom_types.append(element or atom_name)
            chain_ids.append(chain_id)
            insertion_codes.append(insertion_code)
            alternate_locations.append(first(row, "_atom_site.label_alt_id"))
            occupancies.append(
                _mmcif_parse_float(first(row, "_atom_site.occupancy"))
            )
            temperature_factors.append(
                _mmcif_parse_float(first(row, "_atom_site.b_iso_or_equiv"))
            )
            formal_charge = parse_formal_charge(
                first(row, "_atom_site.pdbx_formal_charge")
            )
            formal_charges.append(formal_charge)
            record_types.append(first(row, "_atom_site.group_pdb").upper() or "ATOM")

            label_seq = first(row, "_atom_site.label_seq_id")
            label_comp = first(row, "_atom_site.label_comp_id")
            label_key = identity_key(label_asym, label_seq, label_comp, label_atom, insertion_code)
            auth_key = identity_key(
                auth_asym or label_asym,
                auth_seq or label_seq,
                auth_comp or label_comp,
                auth_atom or label_atom,
                insertion_code,
            )
            if any(label_key):
                atom_by_label_key[label_key] = atom_index
            if any(auth_key):
                atom_by_auth_key[auth_key] = atom_index
            if raw_atom_id:
                atom_by_site_id[raw_atom_id] = atom_index
            for alias in (auth_atom, label_atom):
                if alias:
                    residue_atom_names[residue_index][alias.upper()] = atom_index

        atom_count = len(atom_names)
        residue_count = len(residue_names)
        segment_count = len(segment_ids)

        trajectory_mapping = []
        mapping_rows = _mmcif_rows(data, "mokda_trajectory_mapping")
        if mapping_rows:
            mapping_by_index: dict[int, dict[str, Any]] = {}
            mapping_hashes: set[str] = set()
            canonical_atom_ids: set[int] = set()
            external_ids: set[str] = set()
            simulation_residue_id_by_index: dict[int, str] = {}
            simulation_residue_index_by_id: dict[str, int] = {}
            for row in mapping_rows:
                simulation_index = _mmcif_parse_int(
                    first(row, "_mokda_trajectory_mapping.simulation_index")
                )
                atom_site_id = first(row, "_mokda_trajectory_mapping.atom_site_id")
                canonical_atom_id = _mmcif_parse_int(
                    first(row, "_mokda_trajectory_mapping.canonical_atom_id")
                )
                external_id = first(row, "_mokda_trajectory_mapping.external_id")
                simulation_residue_index = _mmcif_parse_int(
                    first(row, "_mokda_trajectory_mapping.simulation_residue_index")
                )
                simulation_residue_id = first(
                    row, "_mokda_trajectory_mapping.simulation_residue_id"
                )
                mapping_hash = first(row, "_mokda_trajectory_mapping.mapping_hash")
                if (
                    simulation_index is None
                    or canonical_atom_id is None
                    or simulation_residue_index is None
                    or atom_site_id not in atom_by_site_id
                    or not external_id
                    or not simulation_residue_id
                    or canonical_atom_id <= 0
                    or simulation_residue_index < 0
                ):
                    raise BundleTopologyError(
                        "invalid _mokda_trajectory_mapping row"
                    )
                if simulation_index in mapping_by_index:
                    raise BundleTopologyError(
                        "duplicate _mokda_trajectory_mapping.simulation_index"
                    )
                if canonical_atom_id in canonical_atom_ids:
                    raise BundleTopologyError(
                        "duplicate _mokda_trajectory_mapping.canonical_atom_id"
                    )
                if external_id in external_ids:
                    raise BundleTopologyError(
                        "duplicate _mokda_trajectory_mapping.external_id"
                    )
                canonical_atom_ids.add(canonical_atom_id)
                external_ids.add(external_id)
                previous_residue_id = simulation_residue_id_by_index.setdefault(
                    simulation_residue_index, simulation_residue_id
                )
                previous_residue_index = simulation_residue_index_by_id.setdefault(
                    simulation_residue_id, simulation_residue_index
                )
                if (
                    previous_residue_id != simulation_residue_id
                    or previous_residue_index != simulation_residue_index
                ):
                    raise BundleTopologyError(
                        "conflicting _mokda_trajectory_mapping simulation residue identity"
                    )
                mapping_by_index[simulation_index] = {
                    "simulation_index": simulation_index,
                    "atom_site_id": atom_site_id,
                    "canonical_atom_id": canonical_atom_id,
                    "external_id": external_id,
                    "simulation_residue_index": simulation_residue_index,
                    "simulation_residue_id": simulation_residue_id,
                    "mapping_hash": mapping_hash,
                }
                mapping_hashes.add(mapping_hash)

            expected_indices = list(range(len(selected_rows)))
            if sorted(mapping_by_index) != expected_indices:
                raise BundleTopologyError(
                    "_mokda_trajectory_mapping simulation indices must cover "
                    "the selected CIF atoms from 0 to n_atoms - 1"
                )
            trajectory_mapping = [
                mapping_by_index[index] for index in expected_indices
            ]
            mapped_site_ids = [row["atom_site_id"] for row in trajectory_mapping]
            selected_site_ids = [first(row, "_atom_site.id") for row in selected_rows]
            if len(set(selected_site_ids)) != len(selected_site_ids):
                raise BundleTopologyError(
                    "_atom_site.id values must be unique when a Mokda trajectory mapping is present"
                )
            if mapped_site_ids != selected_site_ids:
                raise BundleTopologyError(
                    "_mokda_trajectory_mapping does not match _atom_site row order"
                )
            if len(mapping_hashes) != 1:
                raise BundleTopologyError(
                    "_mokda_trajectory_mapping rows have inconsistent mapping hashes"
                )
            mapping_hash = next(iter(mapping_hashes))
            if (
                len(mapping_hash) != 64
                or mapping_hash != mapping_hash.lower()
                or any(character not in "0123456789abcdef" for character in mapping_hash)
            ):
                raise BundleTopologyError(
                    "_mokda_trajectory_mapping.mapping_hash must be lowercase SHA-256"
                )
            if _mokda_mapping_digest(trajectory_mapping) != mapping_hash:
                raise BundleTopologyError(
                    "_mokda_trajectory_mapping rows do not match mapping_hash"
                )

        partial_charges: list[float | None] = [None] * atom_count
        for row in _mmcif_rows(data, "mokda_charge_state"):
            site_id = first(row, "_mokda_charge_state.atom_site_id")
            atom_index = atom_by_site_id.get(site_id)
            if atom_index is None:
                # A multi-model CIF may contain charge rows for atoms outside
                # the model selected for this topology.
                continue
            partial_charges[atom_index] = _mmcif_parse_float(
                first(row, "_mokda_charge_state.partial_charge")
            )
            if formal_charges[atom_index] is None:
                formal_charges[atom_index] = parse_formal_charge(
                    first(row, "_mokda_charge_state.formal_charge")
                )

        attrs: list[Any] = [
            topologyattrs.Atomids(np.asarray(atom_ids, dtype=np.int64)),
            topologyattrs.Atomnames(np.asarray(atom_names, dtype=object)),
            topologyattrs.Atomtypes(np.asarray(atom_types, dtype=object)),
            topologyattrs.Resids(np.asarray(residue_ids, dtype=np.int64)),
            topologyattrs.Resnums(np.asarray(residue_numbers, dtype=np.int64)),
            topologyattrs.Resnames(np.asarray(residue_names, dtype=object)),
            topologyattrs.Segids(np.asarray(segment_ids, dtype=object)),
            topologyattrs.ChainIDs(np.asarray(chain_ids, dtype=object)),
            topologyattrs.ICodes(np.asarray(residue_insertion_codes, dtype=object)),
            topologyattrs.RecordTypes(np.asarray(record_types, dtype=object)),
        ]
        if any(elements):
            attrs.append(topologyattrs.Elements(np.asarray(elements, dtype=object)))
        if any(alternate_locations):
            attrs.append(topologyattrs.AltLocs(np.asarray(alternate_locations, dtype=object)))
        if any(value is not None for value in occupancies):
            attrs.append(
                topologyattrs.Occupancies(
                    np.asarray([np.nan if value is None else value for value in occupancies], dtype=np.float64)
                )
            )
        if any(value is not None for value in temperature_factors):
            attrs.append(
                topologyattrs.Tempfactors(
                    np.asarray(
                        [np.nan if value is None else value for value in temperature_factors],
                        dtype=np.float64,
                    )
                )
            )
        if any(value is not None for value in formal_charges):
            attrs.append(
                topologyattrs.FormalCharges(
                    np.asarray([0 if value is None else value for value in formal_charges], dtype=np.int64)
                )
            )
        if any(value is not None for value in partial_charges):
            attrs.append(
                topologyattrs.Charges(
                    np.asarray(
                        [np.nan if value is None else value for value in partial_charges],
                        dtype=np.float64,
                    )
                )
            )

        bonds: set[tuple[int, int]] = set()
        bond_orders: dict[tuple[int, int], float] = {}
        bond_types: dict[tuple[int, int], str] = {}

        def parse_bond_order(value: str) -> float | None:
            token = _mmcif_clean(value).upper()
            if token in {"SING", "SINGLE"}:
                return 1.0
            if token in {"DOUB", "DOUBLE"}:
                return 2.0
            if token in {"TRIP", "TRIPLE"}:
                return 3.0
            if token in {"AROM", "AROMATIC", "DELO"}:
                return 1.5
            return _mmcif_parse_float(token)

        def add_bond(
            atom1: int | None,
            atom2: int | None,
            *,
            order: float | None = None,
            bond_type: str = "",
        ) -> None:
            if atom1 is None or atom2 is None or atom1 == atom2:
                return
            bond = tuple(sorted((int(atom1), int(atom2))))
            bonds.add(bond)
            if order is not None:
                bond_orders[bond] = order
            if bond_type:
                bond_types[bond] = bond_type

        def remove_bond(atom1: int | None, atom2: int | None) -> None:
            if atom1 is None or atom2 is None or atom1 == atom2:
                return
            bond = tuple(sorted((int(atom1), int(atom2))))
            bonds.discard(bond)
            bond_orders.pop(bond, None)
            bond_types.pop(bond, None)

        for row in _mmcif_rows(data, "chem_comp_bond"):
            comp_id = first(row, "_chem_comp_bond.comp_id").upper()
            atom1_name = first(row, "_chem_comp_bond.atom_id_1").upper()
            atom2_name = first(row, "_chem_comp_bond.atom_id_2").upper()
            if not comp_id or not atom1_name or not atom2_name:
                continue
            for residue_index, residue_comp in enumerate(residue_components):
                if residue_comp.upper() == comp_id:
                    add_bond(
                        residue_atom_names[residue_index].get(atom1_name),
                        residue_atom_names[residue_index].get(atom2_name),
                        order=parse_bond_order(
                            first(row, "_chem_comp_bond.value_order")
                        ),
                    )

        def resolve_struct_conn_atom(row: dict[str, str], partner: str) -> int | None:
            label_key = identity_key(
                first(row, f"_struct_conn.{partner}_label_asym_id"),
                first(row, f"_struct_conn.{partner}_label_seq_id"),
                first(row, f"_struct_conn.{partner}_label_comp_id"),
                first(row, f"_struct_conn.{partner}_label_atom_id"),
                first(row, f"_struct_conn.pdbx_{partner}_pdb_ins_code"),
            )
            auth_key = identity_key(
                first(row, f"_struct_conn.{partner}_auth_asym_id", f"_struct_conn.{partner}_label_asym_id"),
                first(row, f"_struct_conn.{partner}_auth_seq_id", f"_struct_conn.{partner}_label_seq_id"),
                first(row, f"_struct_conn.{partner}_auth_comp_id", f"_struct_conn.{partner}_label_comp_id"),
                first(row, f"_struct_conn.{partner}_auth_atom_id", f"_struct_conn.{partner}_label_atom_id"),
                first(row, f"_struct_conn.pdbx_{partner}_pdb_ins_code"),
            )
            return atom_by_auth_key.get(auth_key, atom_by_label_key.get(label_key))

        for row in _mmcif_rows(data, "struct_conn"):
            add_bond(
                resolve_struct_conn_atom(row, "ptnr1"),
                resolve_struct_conn_atom(row, "ptnr2"),
                bond_type=first(row, "_struct_conn.conn_type_id"),
            )
        for row in _mmcif_rows(data, "mokda_bond_semantic"):
            add_bond(
                atom_by_site_id.get(first(row, "_mokda_bond_semantic.atom_site_id_1")),
                atom_by_site_id.get(first(row, "_mokda_bond_semantic.atom_site_id_2")),
                order=_mmcif_parse_float(
                    first(row, "_mokda_bond_semantic.order")
                ),
                bond_type=first(row, "_mokda_bond_semantic.bond_type"),
            )
        for row in _mmcif_rows(data, "mokda_edit_operation"):
            operation = first(row, "_mokda_edit_operation.operation_type").lower()
            atom1 = atom_by_site_id.get(first(row, "_mokda_edit_operation.atom_site_id_1"))
            atom2 = atom_by_site_id.get(first(row, "_mokda_edit_operation.atom_site_id_2"))
            if atom1 is None or atom2 is None or atom1 == atom2:
                continue
            bond = tuple(sorted((int(atom1), int(atom2))))
            if operation == "delete_bond":
                remove_bond(atom1, atom2)
            elif operation in {"add_bond", "create_bond", "update_bond"}:
                add_bond(
                    atom1,
                    atom2,
                    order=_mmcif_parse_float(
                        first(row, "_mokda_edit_operation.bond_order")
                    ),
                    bond_type=first(row, "_mokda_edit_operation.bond_type"),
                )
        if bonds:
            ordered_bonds = sorted(bonds)
            bond_kwargs = {}
            if bond_orders:
                bond_kwargs["order"] = np.asarray(
                    [bond_orders.get(bond, np.nan) for bond in ordered_bonds],
                    dtype=np.float64,
                )
            if bond_types:
                bond_kwargs["types"] = np.asarray(
                    [bond_types.get(bond, "") for bond in ordered_bonds],
                    dtype=object,
                )
            attrs.append(topologyattrs.Bonds(ordered_bonds, **bond_kwargs))

        topology = Topology(
            atom_count,
            residue_count,
            segment_count,
            attrs,
            np.asarray(residue_index_by_atom, dtype=np.int32),
            np.asarray(residue_segment_index, dtype=np.int32),
        )
        topology.cif_metadata = {
            # Keep every parsed CIF column, including Mokda-specific categories
            # and auth/label identifiers that have no one-to-one MDA attribute.
            "raw": data,
            "selected_model": selected_model,
            "selected_atom_site_ids": tuple(
                first(row, "_atom_site.id") for row in selected_rows
            ),
            "mokda_trajectory_mapping": trajectory_mapping,
        }
        return topology


class BundleTopologyParser(TopologyReaderBase):
    """Parse ``topology.spgt.h5`` into an MDAnalysis topology."""

    format = "SPONGE_TOPOLOGY_H5"

    @staticmethod
    def _format_hint(thing: Any) -> bool:
        handle, close = _open_hint_target(thing)
        if handle is None:
            return False
        try:
            return _optional_text(handle, "/schema/name") == _TOPOLOGY_SCHEMA
        except (OSError, KeyError, ValueError):
            return False
        finally:
            if close:
                handle.close()

    def parse(self, **kwargs) -> Topology:
        strict = bool(kwargs.get("strict", True))
        with h5py.File(self.filename, "r") as handle:
            _validate_schema(
                handle,
                name_path="/schema/name",
                version_path="/schema/version",
                expected_name=_TOPOLOGY_SCHEMA,
                strict=strict,
            )
            if "/topology/atom_count" not in handle:
                raise BundleTopologyError(
                    f"{handle.filename}:/topology/atom_count is required"
                )
            atom_count = int(np.asarray(handle["/topology/atom_count"][()]).reshape(-1)[0])
            if atom_count <= 0:
                raise BundleTopologyError(
                    f"{handle.filename}:/topology/atom_count must be positive, got {atom_count}"
                )

            atom_resindex, residue_count = self._residue_mapping(handle, atom_count)
            attrs: list[Any] = [
                topologyattrs.Atomids(np.arange(1, atom_count + 1, dtype=np.int64)),
                topologyattrs.Resids(np.arange(1, residue_count + 1, dtype=np.int64)),
                topologyattrs.Resnums(np.arange(1, residue_count + 1, dtype=np.int64)),
                topologyattrs.Segids(np.asarray(["SYSTEM"], dtype=object)),
            ]

            masses = _optional_array(handle, "/atoms/mass", dtype=np.float64)
            if masses is not None:
                masses = _require_vector(
                    masses,
                    atom_count,
                    filename=handle.filename,
                    path="/atoms/mass",
                )
                attrs.append(topologyattrs.Masses(masses))

            names = _read_string_array(handle, "/parameters/xponge/atoms/name")
            type_names = _read_string_array(handle, "/parameters/xponge/atoms/type_name")
            elements: list[str] | None = None
            if masses is not None:
                elements = [guess_element_from_mass(float(mass)) for mass in masses]
                attrs.append(topologyattrs.Elements(elements, guessed=True))
            if names:
                if len(names) != atom_count:
                    raise BundleTopologyError(
                        f"{handle.filename}:/parameters/xponge/atoms/name has {len(names)} values, expected {atom_count}"
                    )
                attrs.append(topologyattrs.Atomnames(names))
            else:
                guessed_names = elements or [f"A{index + 1}" for index in range(atom_count)]
                attrs.append(topologyattrs.Atomnames(guessed_names, guessed=True))
            if type_names:
                if len(type_names) != atom_count:
                    raise BundleTopologyError(
                        f"{handle.filename}:/parameters/xponge/atoms/type_name has {len(type_names)} values, expected {atom_count}"
                    )
                attrs.append(topologyattrs.Atomtypes(type_names))
            else:
                attrs.append(topologyattrs.Atomtypes(names or elements or ["X"] * atom_count, guessed=True))

            if "/atoms/charge" in handle:
                charges = _charge_values(handle["/atoms/charge"], strict=strict)
                _require_vector(
                    charges,
                    atom_count,
                    filename=handle.filename,
                    path="/atoms/charge",
                )
                attrs.append(topologyattrs.Charges(charges))

            residue_names = _read_string_array(handle, "/parameters/xponge/residues/name")
            if residue_names:
                if len(residue_names) != residue_count:
                    raise BundleTopologyError(
                        f"{handle.filename}:/parameters/xponge/residues/name has {len(residue_names)} values, expected {residue_count}"
                    )
            else:
                residue_names = ["SYSTEM"] * residue_count
            attrs.append(topologyattrs.Resnames(residue_names))

            self._append_connectivity(handle, attrs, atom_count)
            residue_segindex = np.zeros(residue_count, dtype=np.int32)
            return Topology(
                atom_count,
                residue_count,
                1,
                attrs,
                atom_resindex,
                residue_segindex,
            )

    @staticmethod
    def _residue_mapping(handle: h5py.File, atom_count: int) -> tuple[np.ndarray, int]:
        residue_index = _optional_array(handle, "/atoms/residue_index", dtype=np.int32)
        offsets = _optional_array(handle, "/residues/atom_offset", dtype=np.int64)
        if residue_index is None and offsets is None:
            return np.zeros(atom_count, dtype=np.int32), 1
        if residue_index is None:
            if offsets.ndim != 1 or offsets.size < 2 or offsets[0] != 0 or offsets[-1] != atom_count:
                raise BundleTopologyError(
                    f"{handle.filename}:/residues/atom_offset is invalid for {atom_count} atoms"
                )
            counts = np.diff(offsets)
            if np.any(counts < 0):
                raise BundleTopologyError(
                    f"{handle.filename}:/residues/atom_offset is not monotonic"
                )
            residue_index = np.repeat(np.arange(counts.size, dtype=np.int32), counts)
            return residue_index, int(counts.size)

        _require_vector(
            residue_index,
            atom_count,
            filename=handle.filename,
            path="/atoms/residue_index",
        )
        if np.any(residue_index < 0):
            raise BundleTopologyError(
                f"{handle.filename}:/atoms/residue_index contains negative indices"
            )
        residue_count = int(residue_index.max()) + 1 if atom_count else 1
        if offsets is not None:
            if offsets.shape != (residue_count + 1,) or offsets[0] != 0 or offsets[-1] != atom_count:
                raise BundleTopologyError(
                    f"{handle.filename}:/residues/atom_offset is inconsistent with residue_index"
                )
            expected = np.repeat(
                np.arange(residue_count, dtype=np.int32),
                np.diff(offsets),
            )
            if not np.array_equal(expected, residue_index):
                raise BundleTopologyError(
                    f"{handle.filename} residue_index and atom_offset describe different mappings"
                )
        return residue_index, residue_count

    @staticmethod
    def _append_connectivity(handle: h5py.File, attrs: list[Any], atom_count: int) -> None:
        definitions = (
            ("/forcefield/bond/atoms", 2, topologyattrs.Bonds, True),
            ("/forcefield/angle/atoms", 3, topologyattrs.Angles, True),
            ("/forcefield/dihedral/atoms", 4, topologyattrs.Dihedrals, True),
            ("/forcefield/improper/atoms", 4, topologyattrs.Impropers, False),
        )
        for path, width, attr_type, normalize_reverse in definitions:
            if path not in handle:
                continue
            values = np.asarray(handle[path][...], dtype=np.int64)
            connections = _deduplicate_connections(
                values,
                width,
                normalize_reverse=normalize_reverse,
            )
            if any(index < 0 or index >= atom_count for item in connections for index in item):
                raise BundleTopologyError(
                    f"{handle.filename}:{path} contains an atom index outside [0, {atom_count})"
                )
            attrs.append(attr_type(connections))


def _creator_name(handle: h5py.File) -> str | None:
    if "/h5md/creator" not in handle:
        return None
    value = handle["/h5md/creator"].attrs.get("name")
    return _decode_text(value) if value is not None else None


def _legacy_groups(handle: h5py.File) -> tuple[str | None, list[tuple[int, str]]]:
    if "/particles" not in handle:
        return None, []
    keys = list(handle["/particles"].keys())
    single = "trajectory" if "trajectory" in keys else None
    numbered: list[tuple[int, str]] = []
    for key in keys:
        match = _LEGACY_WALKER_RE.match(key)
        if match:
            numbered.append((int(match.group("index")), key))
    numbered.sort()
    return single, numbered


def _detect_layout(handle: h5py.File, requested: str = "auto") -> str:
    requested = requested.lower()
    if requested not in {"auto", "bundle", "legacy", "h5md"}:
        raise ValueError(f"layout must be auto, bundle, legacy, or h5md; got {requested!r}")
    if requested != "auto":
        return requested
    schema = _optional_text(handle, "/parameters/sponge/schema/name")
    if schema == _TRAJECTORY_SCHEMA:
        return "bundle"
    single, numbered = _legacy_groups(handle)
    creator = (_creator_name(handle) or "").lower()
    if (single is not None or numbered) and creator and "sponge" not in creator and "xponge" not in creator:
        if "/h5md" in handle:
            return "h5md"
    if single is not None and numbered:
        raise AmbiguousH5MDLayoutError(
            f"{handle.filename} contains both /particles/trajectory and numbered walker groups"
        )
    if single is not None or numbered:
        return "legacy"
    if "/h5md" in handle and "/particles" in handle:
        return "h5md"
    raise BundleTrajectoryError(f"{handle.filename} is not a recognized H5MD trajectory")


def _select_bundle_stream(handle: h5py.File, particle_stream: str | None) -> str:
    if "/particles" not in handle:
        raise BundleTrajectoryError(f"{handle.filename}:/particles is missing")
    groups = list(handle["/particles"].keys())
    if particle_stream is not None:
        if particle_stream not in groups:
            raise BundleTrajectoryError(
                f"{handle.filename} has no particle stream {particle_stream!r}; available: {groups}"
            )
        return particle_stream
    declared = _read_string_array(handle, "/parameters/sponge/output/particle_streams")
    if len(declared) == 1:
        if declared[0] not in groups:
            raise BundleTrajectoryError(
                f"{handle.filename} declares missing particle stream {declared[0]!r}"
            )
        return declared[0]
    if len(declared) > 1:
        raise BundleTrajectoryError(
            f"{handle.filename} declares multiple particle streams {declared}; pass particle_stream"
        )
    if len(groups) == 1:
        return groups[0]
    raise BundleTrajectoryError(
        f"{handle.filename} has multiple particle streams {groups}; pass particle_stream"
    )


def _select_legacy_group(handle: h5py.File, walker: int | None) -> str:
    if walker is None:
        walker = 0
    if not isinstance(walker, int) or isinstance(walker, bool):
        raise TypeError(f"walker must be an int, got {type(walker).__name__}")
    if walker < 0:
        raise ValueError(f"walker must be non-negative, got {walker}")
    single, numbered = _legacy_groups(handle)
    if single is not None and numbered:
        raise AmbiguousH5MDLayoutError(
            f"{handle.filename} contains both single and numbered legacy walker groups"
        )
    if single is not None:
        if walker != 0:
            raise ValueError(f"{handle.filename} contains one walker; walker={walker} is invalid")
        return single
    if not numbered:
        raise BundleTrajectoryError(f"{handle.filename} contains no legacy walker groups")
    indices = [index for index, _ in numbered]
    by_index = dict(numbered)
    if walker not in by_index:
        raise ValueError(
            f"{handle.filename} has walker indices {indices}; walker={walker} is invalid"
        )
    return by_index[walker]


_UNIT_ALIASES = {
    "length": {"Angstrom": "Angstrom", "A": "Angstrom", "angstrom": "Angstrom", "nm": "nm"},
    "time": {"ps": "ps", "fs": "fs", "ns": "ns", "AKMA": "AKMA"},
    "velocity": {
        "Angstrom ps-1": "Angstrom/ps",
        "A ps-1": "Angstrom/ps",
        "Angstrom/ps": "Angstrom/ps",
        "nm ps-1": "nm/ps",
        "nm/ps": "nm/ps",
    },
    "force": {
        "kcal mol-1 Angstrom-1": "kcal/(mol*Angstrom)",
        "kcal/(mol*Angstrom)": "kcal/(mol*Angstrom)",
        "kJ mol-1 Angstrom-1": "kJ/(mol*Angstrom)",
        "kJ/(mol*Angstrom)": "kJ/(mol*Angstrom)",
    },
}

_BASE_UNITS = {
    "length": "Angstrom",
    "time": "ps",
    "velocity": "Angstrom/ps",
    "force": "kJ/(mol*Angstrom)",
}

_UNIT_TYPES = {"length": "length", "time": "time", "velocity": "speed", "force": "force"}


def _dataset_unit(
    dataset: h5py.Dataset | None,
    kind: str,
    *,
    strict: bool,
    legacy: bool,
) -> str | None:
    if dataset is None:
        return None
    value = dataset.attrs.get("unit")
    if value is None:
        if legacy:
            return _BASE_UNITS[kind]
        if strict:
            raise BundleUnitError(f"{dataset.file.filename}:{dataset.name} has no {kind} unit")
        return _BASE_UNITS[kind]
    raw = _decode_text(value)
    try:
        return _UNIT_ALIASES[kind][raw]
    except KeyError as exc:
        raise BundleUnitError(
            f"{dataset.file.filename}:{dataset.name} has unsupported {kind} unit {raw!r}"
        ) from exc


def _unit_factor(kind: str, unit: str | None) -> float:
    if unit is None:
        return 1.0
    try:
        return float(get_conversion_factor(_UNIT_TYPES[kind], unit, _BASE_UNITS[kind]))
    except (KeyError, ValueError) as exc:
        raise BundleUnitError(f"cannot convert {kind} unit {unit!r} to {_BASE_UNITS[kind]!r}") from exc


def _edges_to_dimensions(edges: np.ndarray) -> np.ndarray:
    vectors = np.asarray(edges, dtype=np.float64)
    if vectors.shape == (3,):
        if np.any(vectors <= 0) or not np.all(np.isfinite(vectors)):
            raise BundleTrajectoryError(f"box edges are invalid: {vectors}")
        return np.asarray([*vectors, 90.0, 90.0, 90.0], dtype=np.float32)
    if vectors.shape != (3, 3):
        raise BundleTrajectoryError(f"box edges have shape {vectors.shape}, expected (3,) or (3, 3)")
    lengths = np.linalg.norm(vectors, axis=1)
    if np.any(lengths <= 0) or not np.all(np.isfinite(lengths)):
        raise BundleTrajectoryError("box edges contain a zero or non-finite vector")

    def angle(left: int, right: int) -> float:
        cosine = np.dot(vectors[left], vectors[right]) / (lengths[left] * lengths[right])
        return float(np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0))))

    return np.asarray(
        [*lengths, angle(1, 2), angle(0, 2), angle(0, 1)],
        dtype=np.float32,
    )


class SpongeH5MDReader(base.ReaderBase):
    """Unified MDAnalysis reader for bundled, legacy-walker, and native H5MD."""

    format = "SPONGE_H5MD"

    @staticmethod
    def _format_hint(thing: Any) -> bool:
        handle, close = _open_hint_target(thing)
        if handle is None:
            return False
        try:
            if _optional_text(handle, "/parameters/sponge/schema/name") == _TRAJECTORY_SCHEMA:
                return True
            single, numbered = _legacy_groups(handle)
            if single is None and not numbered:
                return False
            creator = (_creator_name(handle) or "").lower()
            return not creator or "sponge" in creator or "xponge" in creator
        except (OSError, KeyError, ValueError):
            return False
        finally:
            if close:
                handle.close()

    @staticmethod
    def is_h5md(thing: Any) -> bool:
        """Return whether a file has a layout this reader can consume.

        Unlike ``_format_hint``, this includes ordinary native H5MD files.  The
        distinction keeps generic H5MD files out of automatic MDAnalysis
        format guessing while allowing explicit CIF + H5MD loading.
        """

        handle, close = _open_hint_target(thing)
        if handle is None:
            return False
        try:
            _detect_layout(handle)
            return True
        except (AmbiguousH5MDLayoutError, BundleTrajectoryError, OSError, ValueError):
            return False
        finally:
            if close:
                handle.close()

    @classmethod
    def parse_n_atoms(
        cls,
        filename,
        *,
        layout: str = "auto",
        particle_stream: str | None = None,
        walker: int | None = None,
        **kwargs,
    ) -> int:
        with h5py.File(filename, "r") as handle:
            selected_layout = _detect_layout(handle, layout)
            if selected_layout == "bundle":
                group_name = _select_bundle_stream(handle, particle_stream)
            elif selected_layout == "legacy":
                group_name = _select_legacy_group(handle, walker)
            else:
                return int(H5MD.H5MDReader.parse_n_atoms(filename))
            path = f"/particles/{group_name}/position/value"
            if path not in handle:
                raise BundleTrajectoryError(f"{handle.filename}:{path} is missing")
            shape = handle[path].shape
            if len(shape) != 3 or shape[2] != 3:
                raise BundleTrajectoryError(
                    f"{handle.filename}:{path} has shape {shape}, expected (frames, atoms, 3)"
                )
            return int(shape[1])

    def __init__(
        self,
        filename,
        n_atoms: int | None = None,
        *,
        layout: str = "auto",
        particle_stream: str | None = None,
        walker: int | None = None,
        strict: bool = True,
        allow_incomplete: bool = False,
        convert_units: bool = True,
        **kwargs,
    ):
        self._file: h5py.File | None = None
        self._native = None
        self._group_path: str | None = None
        self._frame = -1
        self._n_atoms = 0
        self._n_frames = 0
        self._observables: dict[str, h5py.Dataset] = {}
        self._strict = bool(strict)
        self._allow_incomplete = bool(allow_incomplete)
        self._requested_layout = layout
        self._particle_stream = particle_stream
        self.walker = walker
        super().__init__(
            filename,
            convert_units=convert_units,
            layout=layout,
            particle_stream=particle_stream,
            walker=walker,
            strict=strict,
            allow_incomplete=allow_incomplete,
            **kwargs,
        )
        self.convert_units = bool(convert_units)
        self._file = h5py.File(filename, "r")
        self.layout = _detect_layout(self._file, layout)
        if self.layout == "h5md":
            if particle_stream is not None or walker is not None:
                self.close()
                raise ValueError("particle_stream and walker are not valid for native H5MD layout")
            self.close()
            self._native = H5MD.H5MDReader(filename, convert_units=convert_units, **kwargs)
            self.ts = self._native.ts
            self.units = dict(self._native.units)
            return

        if self.layout == "bundle":
            if walker is not None:
                self.close()
                raise ValueError("walker is only valid for legacy H5MD layout")
            _validate_schema(
                self._file,
                name_path="/parameters/sponge/schema/name",
                version_path="/parameters/sponge/schema/version",
                expected_name=_TRAJECTORY_SCHEMA,
                strict=self._strict,
            )
            group_name = _select_bundle_stream(self._file, particle_stream)
        else:
            if particle_stream is not None:
                self.close()
                raise ValueError("particle_stream is only valid for bundle H5MD layout")
            group_name = _select_legacy_group(self._file, walker)
        self._group_path = f"/particles/{group_name}"
        try:
            self._initialize_particle_group(expected_n_atoms=n_atoms)
        except Exception:
            self.close()
            raise

    @property
    def n_atoms(self) -> int:
        return self._native.n_atoms if self._native is not None else self._n_atoms

    @property
    def n_frames(self) -> int:
        return self._native.n_frames if self._native is not None else self._n_frames

    def _initialize_particle_group(self, *, expected_n_atoms: int | None) -> None:
        group = self._particle_group()
        if "position/value" not in group:
            self.close()
            raise BundleTrajectoryError(f"{self.filename}:{group.name}/position/value is missing")
        position = group["position/value"]
        if len(position.shape) != 3 or position.shape[2] != 3:
            self.close()
            raise BundleTrajectoryError(
                f"{self.filename}:{position.name} has shape {position.shape}, expected (frames, atoms, 3)"
            )
        self._n_atoms = int(position.shape[1])
        if expected_n_atoms is not None and int(expected_n_atoms) != self._n_atoms:
            self.close()
            raise BundleTrajectoryError(
                f"{self.filename}:{position.name} has {self._n_atoms} atoms, topology has {expected_n_atoms}"
            )

        datasets: dict[str, h5py.Dataset] = {"position": position}
        for component in ("velocity", "force"):
            if f"{component}/value" in group:
                datasets[component] = group[f"{component}/value"]
        if "box/edges/value" in group:
            datasets["box"] = group["box/edges/value"]

        counts: dict[str, int] = {}
        for name, dataset in datasets.items():
            if name != "box" and (len(dataset.shape) != 3 or dataset.shape[1:] != (self._n_atoms, 3)):
                self.close()
                raise BundleTrajectoryError(
                    f"{self.filename}:{dataset.name} has shape {dataset.shape}, expected (*, {self._n_atoms}, 3)"
                )
            if name == "box" and (
                len(dataset.shape) not in {2, 3}
                or dataset.shape[1:] not in {(3,), (3, 3)}
            ):
                self.close()
                raise BundleTrajectoryError(
                    f"{self.filename}:{dataset.name} has unsupported box shape {dataset.shape}"
                )
            counts[name] = int(dataset.shape[0])

        for axis in (self._step_dataset(group), self._time_dataset(group)):
            if axis is not None:
                counts[axis.name] = int(axis.shape[0])
        unique_counts = set(counts.values())
        if len(unique_counts) > 1 and not self._allow_incomplete:
            self.close()
            raise BundleTrajectoryError(
                f"{self.filename}:{group.name} has mismatched frame counts {counts}"
            )
        self._n_frames = min(unique_counts) if unique_counts else 0
        if self._n_frames <= 0:
            self.close()
            raise BundleTrajectoryError(f"{self.filename}:{group.name} has no complete frames")
        self._validate_completion()

        legacy = self.layout == "legacy"
        time_dataset = self._time_dataset(group)
        self.units = {
            "length": _dataset_unit(position, "length", strict=self._strict, legacy=legacy),
            "time": _dataset_unit(time_dataset, "time", strict=self._strict, legacy=legacy),
            "velocity": _dataset_unit(datasets.get("velocity"), "velocity", strict=self._strict, legacy=legacy),
            "force": _dataset_unit(datasets.get("force"), "force", strict=self._strict, legacy=legacy),
        }
        self._unit_factors = {
            kind: _unit_factor(kind, unit) for kind, unit in self.units.items()
        }
        self.ts = self._Timestep(
            self._n_atoms,
            positions=True,
            velocities="velocity" in datasets,
            forces="force" in datasets,
            **self._ts_kwargs,
        )
        if time_dataset is not None and self._n_frames > 1:
            dt = float(time_dataset[1]) - float(time_dataset[0])
            if self.convert_units:
                dt *= self._unit_factors["time"]
            self.ts.dt = dt
        self._discover_observables()
        self._frame = -1
        self.ts.frame = -1
        self._read_next_timestep()

    def _particle_group(self) -> h5py.Group:
        if self._file is None or self._group_path is None:
            raise BundleTrajectoryError("trajectory is closed")
        return self._file[self._group_path]

    @staticmethod
    def _step_dataset(group: h5py.Group) -> h5py.Dataset | None:
        for path in ("position/step", "step", "box/edges/step", "velocity/step", "force/step"):
            if path in group:
                return group[path]
        return None

    @staticmethod
    def _time_dataset(group: h5py.Group) -> h5py.Dataset | None:
        for path in ("position/time", "time", "box/edges/time", "velocity/time", "force/time"):
            if path in group:
                return group[path]
        return None

    def _validate_completion(self) -> None:
        if self.layout != "bundle" or self._file is None:
            return
        status = _optional_text(self._file, "/parameters/sponge/output/status")
        if status != "finalized" and not self._allow_incomplete:
            raise IncompleteBundleError(
                f"{self.filename} output status is {status!r}, expected 'finalized'"
            )
        path = "/parameters/sponge/output/frame_count"
        if path in self._file and status == "finalized":
            values = np.asarray(self._file[path][...]).reshape(-1)
            if values.size == 0:
                raise IncompleteBundleError(f"{self.filename}:{path} is empty")
            # SPONGE records completion as an append-only commit journal. The
            # last entry is authoritative; a one-element legacy dataset still
            # works through the same rule.
            declared = int(values[-1])
            if declared != self._n_frames:
                raise IncompleteBundleError(
                    f"{self.filename}:{path} is {declared}, actual complete frames are {self._n_frames}"
                )

    def _discover_observables(self) -> None:
        self._observables = {}
        if self.layout != "bundle" or self._file is None:
            return
        stream = self._group_path.rsplit("/", 1)[-1]
        root_path = f"/observables/{stream}"
        if root_path not in self._file:
            return
        root = self._file[root_path]

        def visitor(name: str, obj: Any) -> None:
            if not isinstance(obj, h5py.Group) or "value" not in obj:
                return
            dataset = obj["value"]
            if dataset.ndim >= 1 and dataset.shape[0] >= self._n_frames:
                self._observables[f"observables/{name}"] = dataset

        root.visititems(visitor)

    def _read_frame(self, frame: int):
        if self._native is not None:
            self.ts = self._native._read_frame(frame)
            return self.ts
        if frame < 0 or frame >= self._n_frames:
            raise IndexError(f"frame {frame} is outside [0, {self._n_frames})")
        group = self._particle_group()
        ts = self.ts
        ts.positions = np.asarray(group["position/value"][frame], dtype=np.float32)
        if self.convert_units:
            ts.positions *= self._unit_factors["length"]
        if ts.has_velocities:
            ts.velocities = np.asarray(group["velocity/value"][frame], dtype=np.float32)
            if self.convert_units:
                ts.velocities *= self._unit_factors["velocity"]
        if ts.has_forces:
            ts.forces = np.asarray(group["force/value"][frame], dtype=np.float32)
            if self.convert_units:
                ts.forces *= self._unit_factors["force"]
        if "box/edges/value" in group:
            edges = np.asarray(group["box/edges/value"][frame], dtype=np.float64)
            if self.convert_units:
                edges *= self._unit_factors["length"]
            ts.dimensions = _edges_to_dimensions(edges)
        else:
            ts.dimensions = None
        time_dataset = self._time_dataset(group)
        if time_dataset is not None:
            time_value = float(time_dataset[frame])
            if self.convert_units:
                time_value *= self._unit_factors["time"]
            ts.time = time_value
        step_dataset = self._step_dataset(group)
        if step_dataset is not None:
            ts.data["step"] = int(step_dataset[frame])
        for key, dataset in self._observables.items():
            ts.data[key] = np.asarray(dataset[frame])
        self._frame = frame
        ts.frame = frame
        return ts

    def _read_next_timestep(self):
        if self._native is not None:
            self.ts = self._native._read_next_timestep()
            return self.ts
        next_frame = self._frame + 1
        if next_frame >= self._n_frames:
            raise EOFError
        return self._read_frame(next_frame)

    def _reopen(self) -> None:
        if self._native is not None:
            self._native._reopen()
            self.ts = self._native.ts
            return
        if self._file is not None:
            self._file.close()
        self._file = h5py.File(self.filename, "r")
        self._frame = -1
        self.ts.frame = -1
        self._discover_observables()

    def close(self) -> None:
        native = getattr(self, "_native", None)
        if native is not None:
            native.close()
            self._native = None
        handle = getattr(self, "_file", None)
        if handle is not None:
            handle.close()
            self._file = None


# Historical public spelling; it remains fully supported.
SPONGEH5MDReader = SpongeH5MDReader


def register_mdanalysis_formats() -> tuple[str, str]:
    """Return the idempotently registered XPONGE MDAnalysis format names."""

    return BundleTopologyParser.format, SpongeH5MDReader.format


def validate_bundle_pair(
    topology,
    trajectory,
    *,
    particle_stream: str | None = None,
    strict: bool = True,
    allow_unverified_pair: bool = False,
    allow_incomplete: bool = False,
) -> str:
    """Validate explicit topology/trajectory paths and return the selected stream."""

    with h5py.File(topology, "r") as top, h5py.File(trajectory, "r") as traj:
        _validate_schema(
            top,
            name_path="/schema/name",
            version_path="/schema/version",
            expected_name=_TOPOLOGY_SCHEMA,
            strict=strict,
        )
        _validate_schema(
            traj,
            name_path="/parameters/sponge/schema/name",
            version_path="/parameters/sponge/schema/version",
            expected_name=_TRAJECTORY_SCHEMA,
            strict=strict,
        )
        stream = _select_bundle_stream(traj, particle_stream)
        position_path = f"/particles/{stream}/position/value"
        if position_path not in traj:
            raise BundleTrajectoryError(f"{traj.filename}:{position_path} is missing")
        shape = traj[position_path].shape
        if len(shape) != 3 or shape[2] != 3:
            raise BundleTrajectoryError(
                f"{traj.filename}:{position_path} has shape {shape}, expected (frames, atoms, 3)"
            )
        if "/topology/atom_count" not in top:
            raise BundleTopologyError(f"{top.filename}:/topology/atom_count is missing")
        atom_count = int(np.asarray(top["/topology/atom_count"][()]).reshape(-1)[0])
        if shape[1] != atom_count:
            raise BundleTrajectoryError(
                f"{traj.filename}:{position_path} has {shape[1]} atoms, {top.filename} declares {atom_count}"
            )

        pairs = (
            (
                "/topology/topology_hash",
                "/parameters/sponge/topology_compatibility/topology_hash",
            ),
            (
                "/topology/atom_order_hash",
                "/parameters/sponge/topology_compatibility/atom_order_hash",
            ),
        )
        missing = [(left, right) for left, right in pairs if left not in top or right not in traj]
        if missing and not allow_unverified_pair:
            raise UnverifiedBundlePairError(
                f"cannot prove topology/trajectory compatibility; missing hash paths: {missing}"
            )
        if missing:
            warnings.warn(
                "topology/trajectory compatibility is only verified by atom count",
                RuntimeWarning,
                stacklevel=2,
            )
        for top_path, traj_path in pairs:
            if top_path not in top or traj_path not in traj:
                continue
            top_hash = _read_text(top, top_path)
            traj_hash = _read_text(traj, traj_path)
            if top_hash != traj_hash:
                raise UnverifiedBundlePairError(
                    f"{traj.filename}:{traj_path} is {traj_hash!r}, expected {top_hash!r} from {top.filename}:{top_path}"
                )

        status = _optional_text(traj, "/parameters/sponge/output/status")
        if status != "finalized" and not allow_incomplete:
            raise IncompleteBundleError(
                f"{traj.filename} output status is {status!r}, expected 'finalized'"
            )
        return stream


def _mokda_mapping_file_for_cif(topology) -> str | None:
    """Find Mokda's standard mapping sidecar beside an ordered topology CIF."""

    if not isinstance(topology, (str, os.PathLike)):
        return None
    filename = os.fspath(topology)
    suffix = "_trajectory_topology.cif"
    basename = os.path.basename(filename)
    if not basename.endswith(suffix):
        return None
    prefix = basename[: -len(suffix)]
    return os.path.join(
        os.path.dirname(filename), f"{prefix}_atom_order_mapping.json"
    )


def _verify_mokda_mapping_file(mapping_file, cif_mapping):
    """Compare Mokda's JSON mapping sidecar with its CIF mapping category."""

    if not cif_mapping:
        raise BundleTopologyError(
            "a Mokda atom-order mapping file is present, but the CIF has no "
            "_mokda_trajectory_mapping category"
        )
    try:
        with open(mapping_file, "rt", encoding="utf-8") as handle:
            document = json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise BundleTopologyError(
            f"could not read Mokda atom-order mapping file {mapping_file!r}: {exc}"
        ) from exc
    if (
        not isinstance(document, dict)
        or document.get("schema") != "sponge-atom-order-mapping"
    ):
        raise BundleTopologyError(
            f"{mapping_file!r} is not a sponge-atom-order-mapping document"
        )
    mapping_atoms = document.get("atoms")
    if not isinstance(mapping_atoms, list):
        raise BundleTopologyError(
            f"{mapping_file!r} does not contain an atoms list"
        )

    fields = (
        "simulation_index",
        "external_id",
        "canonical_atom_id",
        "simulation_residue_index",
        "simulation_residue_id",
    )
    try:
        sidecar_mapping = [
            {field: row[field] for field in fields}
            for row in mapping_atoms
            if isinstance(row, dict)
        ]
    except KeyError as exc:
        raise BundleTopologyError(
            f"{mapping_file!r} has an incomplete atom mapping row"
        ) from exc
    if len(sidecar_mapping) != len(mapping_atoms):
        raise BundleTopologyError(
            f"{mapping_file!r} contains a non-object atom mapping row"
        )
    cif_identity_mapping = [
        {field: row[field] for field in fields} for row in cif_mapping
    ]
    cif_mapping_hash = cif_mapping[0]["mapping_hash"]
    if document.get("mapping_hash") != cif_mapping_hash:
        raise BundleTopologyError(
            "Mokda mapping file hash does not match the CIF trajectory mapping"
        )
    if sidecar_mapping != cif_identity_mapping:
        raise BundleTopologyError(
            "Mokda mapping file atoms do not match the CIF trajectory mapping"
        )
    return {
        "path": os.fspath(mapping_file),
        "verified": True,
        "mapping_hash": cif_mapping_hash,
        "atom_count": len(cif_identity_mapping),
    }


def load_cif_h5md_universe(
    topology,
    trajectory,
    *,
    mapping_file=None,
    model_id: str | int | None = None,
    altloc: str | None = "A",
    particle_stream: str | None = None,
    walker: int | None = None,
    strict: bool = True,
    allow_incomplete: bool = False,
    convert_units: bool = True,
    **universe_kwargs,
):
    """Load an mmCIF topology and an H5MD trajectory into an MDAnalysis Universe.

    CIF topology parsing happens before constructing the Universe so CIF-only
    options such as ``model_id`` and ``altloc`` are not forwarded to the H5MD
    reader. Mokda's standard atom-order mapping sidecar is checked when present.
    """

    if not CIFTopologyParser._format_hint(topology):
        raise ValueError(f"{topology!r} is not a recognized CIF/mmCIF topology path")
    if not SpongeH5MDReader.is_h5md(trajectory):
        raise BundleTrajectoryError(f"{trajectory!r} is not a recognized H5MD trajectory")

    cif_topology = CIFTopologyParser(topology).parse(model_id=model_id, altloc=altloc)
    if mapping_file is None:
        mapping_file = _mokda_mapping_file_for_cif(topology)
    if mapping_file is not None and os.path.isfile(mapping_file):
        cif_topology.cif_metadata["mapping_file_validation"] = (
            _verify_mokda_mapping_file(
                mapping_file,
                cif_topology.cif_metadata["mokda_trajectory_mapping"],
            )
        )
    else:
        cif_topology.cif_metadata["mapping_file_validation"] = {
            "available": False,
            "verified": False,
        }
    trajectory_atom_count = SpongeH5MDReader.parse_n_atoms(
        trajectory,
        particle_stream=particle_stream,
        walker=walker,
    )
    if trajectory_atom_count != cif_topology.n_atoms:
        raise BundleTrajectoryError(
            f"{trajectory!r} has {trajectory_atom_count} atoms, CIF topology has "
            f"{cif_topology.n_atoms}"
        )
    if cif_topology.cif_metadata["mokda_trajectory_mapping"]:
        mapping_file_verified = cif_topology.cif_metadata["mapping_file_validation"][
            "verified"
        ]
        mapping_status = (
            "the CIF simulation mapping matches its companion mapping file"
            if mapping_file_verified
            else "the CIF simulation mapping matches its atom-site row order and digest"
        )
        warning = (
            f"CIF/H5MD compatibility: {mapping_status}, and the H5MD atom count "
            "matches; this H5MD file has no per-atom IDs for an independent "
            "coordinate-order check."
        )
    else:
        warning = (
            "CIF/H5MD compatibility is checked by atom count only; ensure the CIF "
            "atom order matches the trajectory atom order."
        )
    warnings.warn(warning, RuntimeWarning, stacklevel=2)
    universe = mda.Universe(
        cif_topology,
        trajectory,
        format=SpongeH5MDReader,
        particle_stream=particle_stream,
        walker=walker,
        strict=strict,
        allow_incomplete=allow_incomplete,
        convert_units=convert_units,
        **universe_kwargs,
    )
    universe.cif_metadata = cif_topology.cif_metadata
    return universe


def load_bundle_universe(
    topology,
    trajectory,
    *,
    particle_stream: str | None = None,
    strict: bool = True,
    allow_unverified_pair: bool = False,
    allow_incomplete: bool = False,
    convert_units: bool = True,
    **universe_kwargs,
):
    """Validate a bundle pair and load it through registered MDAnalysis formats."""

    stream = validate_bundle_pair(
        topology,
        trajectory,
        particle_stream=particle_stream,
        strict=strict,
        allow_unverified_pair=allow_unverified_pair,
        allow_incomplete=allow_incomplete,
    )
    return mda.Universe(
        topology,
        trajectory,
        topology_format=BundleTopologyParser.format,
        format=SpongeH5MDReader.format,
        particle_stream=stream,
        strict=strict,
        allow_incomplete=allow_incomplete,
        convert_units=convert_units,
        **universe_kwargs,
    )


__all__ = [
    "CIFTopologyParser",
    "BundleTopologyParser",
    "SpongeH5MDReader",
    "SPONGEH5MDReader",
    "register_mdanalysis_formats",
    "validate_bundle_pair",
    "load_cif_h5md_universe",
    "load_bundle_universe",
]
