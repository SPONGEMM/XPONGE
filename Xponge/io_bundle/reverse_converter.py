"""Bundled SPONGE input to direct/legacy input conversion."""

from __future__ import annotations

import json
from pathlib import Path

from .bundle_case import BundleCase, scan_bundle_case
from .bundle_reader import BundleReader
from .contracts import (
    CONTRACTS,
    IOContract,
    contracts_by_legacy_key,
    reversible_contracts,
    validate_contract_registry,
)
from .errors import (
    BundleCapabilityError,
    BundleExportError,
    BundlePathError,
    BundleValidationError,
)
from .exporters import EXPORTERS, ExportContext, validate_exporter_registry
from .legacy_case import render_mdin_without_keys
from .legacy_materializer import LegacyMaterializer, LegacyPayload
from .manifest import ManifestEntry, ReverseConversionManifest


_H5_INPUT_KEYS = {
    "input_h5_topology_path",
    "input_h5_protocol_path",
    "input_h5_restart_path",
    "input_h5_restart_load",
    "input_h5_trajectory_path",
    "input_h5_trajectory_particle_stream",
    "default_in_file_prefix",
}


class BundleToLegacyConverter:
    """Materialize direct/legacy SPONGE input files from bundled artifacts."""

    def __init__(
        self,
        case: BundleCase,
        output_dir: str | Path,
        *,
        prefix: str | None = None,
        strict: bool = True,
        overwrite: bool = False,
    ):
        self.case = case
        self.output_dir = Path(output_dir).resolve()
        self.prefix = prefix or "input"
        self.strict = strict
        self.materializer = LegacyMaterializer(self.output_dir, overwrite=overwrite)
        self.manifest = ReverseConversionManifest(
            bundle_root=str(self.case.root),
            output_root=str(self.output_dir),
            mode=self.case.mode,
        )
        self._contracts_by_key = contracts_by_legacy_key()
        self._materialized_keys: set[str] = set()
        self._materialized_groups: set[str] = set()

    def convert(self, *, dry_run: bool = False) -> ReverseConversionManifest:
        """Run reverse conversion and return its manifest."""

        validate_contract_registry()
        validate_exporter_registry()
        with BundleReader(self.case, strict=self.strict) as reader:
            self.manifest.warnings.extend(reader.warnings)
            sidecars = {
                bundle_file: reader.read_legacy_sidecars(bundle_file)
                for bundle_file in (
                    "topology.spgt.h5",
                    "protocol.spgp.h5",
                    "restart.spgr.h5",
                    "trajectory.spg.h5md",
                )
                if reader.has_bundle_file(bundle_file)
            }
            context = ExportContext(
                mode=self.case.mode,
                prefix=self.prefix,
                particle_stream=self.case.particle_stream,
            )
            contracts = sorted(
                reversible_contracts(self.case.mode),
                key=lambda contract: {
                    "typed_required": 0,
                    "typed_or_sidecar": 1,
                    "embedded_text": 2,
                    "sidecar_only": 3,
                    "scalar": 4,
                }.get(contract.reverse_policy, 5),
            )
            for contract in contracts:
                self._materialize_contract(contract, reader, sidecars, context)
            self._materialize_dynamic_sidecars()
            self._materialize_xponge_metadata(reader)

        mdin_target = self._plan_legacy_mdin()
        self.manifest.generated_mdin = str(mdin_target)
        self.manifest.add(
            ManifestEntry(
                contract_id="run_mdin.legacy_generated",
                status="mdin_binding_generated",
                source_path=str(self.case.mdin_path),
                target_path=str(mdin_target),
                bundle_file="run.mdin",
                direction="input",
                component="run_policy",
                payload_kind="file",
                source_kind="bundled_mdin",
            )
        )

        if not dry_run:
            manifest_target = self.output_dir / "manifest.bundle_to_legacy.json"
            self.materializer.plan_payload(
                LegacyPayload(
                    key="__manifest__",
                    data=json.dumps(self.manifest.to_dict(), indent=2, sort_keys=True) + "\n",
                    source_kind="manifest",
                ),
                manifest_target.name,
            )
            self.materializer.write_all()
        else:
            self.materializer.validate_targets()
        return self.manifest

    def _materialize_contract(
        self,
        contract: IOContract,
        reader: BundleReader,
        sidecars: dict[str, dict[str, Path]],
        context: ExportContext,
    ) -> None:
        if contract.payload_kind != "file" or contract.bundle_file in {"run.mdin", "*.legacy"}:
            return
        key = contract.legacy_keys[0]
        if key in self._materialized_keys:
            return
        if contract.materialization_group in self._materialized_groups:
            return
        if not reader.has_bundle_file(contract.bundle_file):
            return

        required_paths = contract.required_bundle_paths
        if contract.bundle_file == "trajectory.spg.h5md" and context.particle_stream != "all":
            required_paths = tuple(
                path.replace("/particles/all/", f"/particles/{context.particle_stream}/")
                for path in required_paths
            )
        typed_present = bool(required_paths) and all(
            reader.contains(contract.bundle_file, path) for path in required_paths
        )
        payloads: list[LegacyPayload] | None = None
        status = "typed_exported"
        source_path = str(self.case.path_for_bundle_file(contract.bundle_file))

        exporter = EXPORTERS.get(contract.exporter_id or "")
        if typed_present and exporter is not None:
            try:
                payloads = exporter(contract, reader, context)
            except (ValueError, TypeError, IndexError) as exc:
                raise BundleExportError(
                    f"failed to export {contract.contract_id} from "
                    f"{contract.bundle_file}:{contract.bundle_path}"
                ) from exc
        else:
            fallback = self._fallback_payload(contract, reader, sidecars)
            if fallback is not None:
                payloads = [fallback]
                status = {
                    "sidecar": "sidecar_restored",
                    "embedded": "embedded_exported",
                    "compatibility": "compatibility_restored",
                }[fallback.source_kind]
                source_path = fallback.source_path or source_path
            elif typed_present:
                message = f"no reverse exporter or fallback is available for {contract.contract_id}"
                if self.strict:
                    raise BundleCapabilityError(message)
                self.manifest.add(
                    ManifestEntry(
                        contract_id=contract.contract_id,
                        status="unsupported",
                        source_key=key,
                        source_path=source_path,
                        bundle_file=contract.bundle_file,
                        bundle_path=contract.bundle_path,
                        direction="input",
                        component=contract.component,
                        payload_kind=contract.payload_kind,
                        message=message,
                    )
                )
                return
            else:
                return

        if not payloads:
            return
        for payload in payloads:
            payload_contract = self._canonical_contract_for_key(payload.key, contract)
            filename = payload.filename or self._filename_for_contract(payload_contract)
            target = self.materializer.plan_payload(payload, filename)
            self._materialized_keys.add(payload.key)
            entry_contract_id = (
                payload_contract.contract_id
                if payload.key in payload_contract.legacy_keys
                else contract.contract_id
            )
            if payload.key not in payload_contract.legacy_keys:
                entry_contract_id += ".data." + payload.key.removesuffix("_in_file")
            self.manifest.add(
                ManifestEntry(
                    contract_id=entry_contract_id,
                    status=status,
                    source_key=payload.key,
                    source_path=source_path,
                    target_path=str(target),
                    bundle_file=contract.bundle_file,
                    bundle_path=contract.bundle_path,
                    direction="input",
                    component=contract.component,
                    payload_kind="file",
                    override_policy=contract.override_policy,
                    comparison_rule=contract.comparison_rule,
                    source_kind=payload.source_kind,
                )
            )
        if contract.materialization_group is not None:
            self._materialized_groups.add(contract.materialization_group)

    def _fallback_payload(
        self,
        contract: IOContract,
        reader: BundleReader,
        sidecars: dict[str, dict[str, Path]],
    ) -> LegacyPayload | None:
        key = contract.legacy_keys[0]
        if contract.reverse_policy in {"embedded_text", "sidecar_only"} and reader.contains(
            contract.bundle_file, contract.bundle_path
        ):
            return LegacyPayload(
                key=key,
                data=reader.read_text(contract.bundle_file, contract.bundle_path),
                source_kind="embedded",
            )

        sidecar_path = sidecars.get(contract.bundle_file, {}).get(key)
        if sidecar_path is not None:
            data = sidecar_path.read_bytes()
            try:
                rendered: str | bytes = data.decode("utf-8")
                binary = False
            except UnicodeDecodeError:
                rendered = data
                binary = True
            return LegacyPayload(
                key=key,
                data=rendered,
                binary=binary,
                source_kind="sidecar",
                source_path=str(sidecar_path),
            )

        compatibility_root = "/compatibility/legacy_import/" + _safe_h5_name(key)
        raw_text_path = compatibility_root + "/raw_text"
        raw_bytes_path = compatibility_root + "/raw_bytes"
        if reader.contains(contract.bundle_file, raw_text_path):
            return LegacyPayload(
                key=key,
                data=reader.read_text(contract.bundle_file, raw_text_path),
                source_kind="compatibility",
            )
        if reader.contains(contract.bundle_file, raw_bytes_path):
            values = reader.read(contract.bundle_file, raw_bytes_path)
            return LegacyPayload(
                key=key,
                data=bytes(values.astype("uint8").reshape(-1)),
                binary=True,
                source_kind="compatibility",
            )
        return None

    def _canonical_contract_for_key(self, key: str, fallback: IOContract) -> IOContract:
        for contract in self._contracts_by_key.get(key, ()):
            if contract.reverse_policy not in {"alias", "sidecar_only", "not_reversible"}:
                return contract
        return fallback

    def _materialize_dynamic_sidecars(self) -> None:
        root = self.case.root.resolve()
        for key, raw_path in self.case.commands.items():
            if not key.endswith("_in_file") or not raw_path.startswith("legacy_sidecars/"):
                continue
            if key in self._materialized_keys:
                continue
            source_path = (root / raw_path).resolve()
            try:
                source_path.relative_to(root)
            except ValueError as exc:
                raise BundlePathError(
                    f"dynamic legacy sidecar for {key} escapes bundle root: {raw_path}"
                ) from exc
            if not source_path.is_file():
                message = f"dynamic legacy sidecar for {key} does not exist: {source_path}"
                if self.strict:
                    raise BundleValidationError(message)
                self.manifest.warnings.append(message)
                continue
            data = source_path.read_bytes()
            try:
                rendered: str | bytes = data.decode("utf-8")
                binary = False
            except UnicodeDecodeError:
                rendered = data
                binary = True
            filename = f"{self.prefix}_{key.removesuffix('_in_file')}.txt"
            target = self.materializer.plan_payload(
                LegacyPayload(
                    key=key,
                    data=rendered,
                    binary=binary,
                    source_kind="sidecar",
                    source_path=str(source_path),
                ),
                filename,
            )
            self._materialized_keys.add(key)
            self.manifest.add(
                ManifestEntry(
                    contract_id=f"topology.dynamic_sidecar.{key}",
                    status="sidecar_restored",
                    source_key=key,
                    source_path=str(source_path),
                    target_path=str(target),
                    bundle_file="topology.spgt.h5",
                    direction="input",
                    component="topology",
                    payload_kind="file",
                    source_kind="sidecar",
                )
            )

    def _materialize_xponge_metadata(self, reader: BundleReader) -> None:
        for key, path in (
            ("resname", "/parameters/xponge/residues/name"),
            ("atom_name", "/parameters/xponge/atoms/name"),
            ("atom_type_name", "/parameters/xponge/atoms/type_name"),
        ):
            if not reader.contains("topology.spgt.h5", path):
                continue
            values = reader.read("topology.spgt.h5", path).reshape(-1)
            rendered_values = [
                value.decode("utf-8") if isinstance(value, bytes) else str(value)
                for value in values
            ]
            payload = LegacyPayload(
                key=key,
                data=f"{len(rendered_values)}\n" + "\n".join(rendered_values) + "\n",
                source_kind="typed",
                bind_in_mdin=False,
            )
            target = self.materializer.plan_payload(payload, f"{self.prefix}_{key}.txt")
            self.manifest.add(
                ManifestEntry(
                    contract_id=f"topology.xponge_metadata.{key}",
                    status="typed_exported",
                    source_key=key,
                    source_path=str(self.case.topology_path),
                    target_path=str(target),
                    bundle_file="topology.spgt.h5",
                    bundle_path=path,
                    direction="input",
                    component="topology_metadata",
                    payload_kind="file",
                    source_kind="typed",
                )
            )

    def _filename_for_contract(self, contract: IOContract) -> str:
        stem = contract.legacy_filename_stem or contract.legacy_keys[0].removesuffix("_in_file")
        return f"{self.prefix}_{stem}.txt"

    def _plan_legacy_mdin(self) -> Path:
        bindings = {
            payload.key: path.relative_to(self.output_dir).as_posix()
            for path, payload in self.materializer.planned.items()
            if not payload.key.startswith("__") and payload.bind_in_mdin
        }
        omit_keys = set(_H5_INPUT_KEYS)
        omit_keys.update(
            key
            for contract in CONTRACTS
            if contract.direction == "input" and contract.payload_kind == "file"
            for key in contract.legacy_keys
        )
        omit_keys.update(
            key
            for key, value in self.case.commands.items()
            if key.endswith("_in_file") and value.startswith("legacy_sidecars/")
        )

        root_lines = []
        section_lines: dict[str, list[str]] = {}
        for key, filename in sorted(bindings.items()):
            known_contracts = self._contracts_by_key.get(key)
            contract = (
                self._canonical_contract_for_key(key, known_contracts[0])
                if known_contracts
                else None
            )
            line_key = key
            if contract is not None and contract.legacy_section:
                prefix = contract.legacy_section + "_"
                line_key = key[len(prefix) :] if key.startswith(prefix) else key
                section_lines.setdefault(contract.legacy_section, []).append(
                    f'{line_key} = "{filename}"'
                )
            else:
                root_lines.append(f'{line_key} = "{filename}"')
        for section, lines in sorted(section_lines.items()):
            root_lines.append(f"[{section}]")
            root_lines.extend(lines)

        rendered = render_mdin_without_keys(self.case.mdin_text, omit_keys, root_lines)
        return self.materializer.plan_payload(
            LegacyPayload(
                key="__mdin__",
                data=rendered,
                source_kind="bundled_mdin",
            ),
            "mdin.legacy.spg.toml",
        )


def convert_bundle_to_legacy(
    bundle_root: str | Path,
    output_dir: str | Path,
    *,
    mdin: str | Path = "mdin.bundled.spg.toml",
    prefix: str | None = None,
    strict: bool = True,
    overwrite: bool = False,
    dry_run: bool = False,
) -> ReverseConversionManifest:
    """Convert a bundled SPONGE input case to direct/legacy files."""

    case = scan_bundle_case(bundle_root, mdin=mdin, strict=strict)
    return BundleToLegacyConverter(
        case,
        output_dir,
        prefix=prefix,
        strict=strict,
        overwrite=overwrite,
    ).convert(dry_run=dry_run)


def _safe_h5_name(name: str) -> str:
    safe = "".join(char if char.isalnum() or char == "_" else "_" for char in name)
    return safe or "unnamed"
