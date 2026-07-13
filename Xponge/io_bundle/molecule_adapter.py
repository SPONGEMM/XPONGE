"""Adapter from XPONGE molecule serializers to typed bundled datasets."""

from __future__ import annotations

import numpy as np

from Xponge.build import _iter_sponge_input_payloads

from .bundle_builder import BundleBuilder
from .contracts import classify_sponge_serializer_key, contracts_by_legacy_key
from .errors import BundleCapabilityError, BundleExportError
from .payloads import MemoryPayload
from .topology_parsers import (
    listed_force_schemas,
    parse_listed_force_data_file,
    parse_topology_file,
)
from .trajectory_parsers import box_to_edges


_METADATA_PATHS = {
    "resname": "/parameters/xponge/residues/name",
    "atom_name": "/parameters/xponge/atoms/name",
    "atom_type_name": "/parameters/xponge/atoms/type_name",
}


def add_molecule_to_bundle(mol, builder: BundleBuilder) -> None:
    """Materialize registered XPONGE serializer payloads into a bundle."""

    contracts = contracts_by_legacy_key()
    serializer_payloads = list(_iter_sponge_input_payloads(mol))
    listed_schemas = {}
    for serializer_key, text in serializer_payloads:
        if serializer_key != "listed_forces":
            continue
        payload = MemoryPayload("listed_forces.txt", text)
        listed_schemas.update(
            {
                name: (parameter_types, parameter_names)
                for name, parameter_types, parameter_names in listed_force_schemas(payload)  # type: ignore[arg-type]
            }
        )
    saw_coordinate = False
    saw_topology = False
    for serializer_key, text in serializer_payloads:
        if serializer_key in listed_schemas:
            parameter_types, parameter_names = listed_schemas[serializer_key]
            datasets = parse_listed_force_data_file(
                serializer_key,
                MemoryPayload(f"{serializer_key}.txt", text),  # type: ignore[arg-type]
                parameter_types=parameter_types,
                parameter_names=parameter_names,
            )
            builder.add_datasets("topology.spgt.h5", datasets)
            saw_topology = True
            continue
        classification, legacy_key = classify_sponge_serializer_key(serializer_key)
        if classification == "metadata":
            values = _parse_counted_strings(serializer_key, text)
            builder.add_dataset(
                "topology.spgt.h5",
                _METADATA_PATHS[serializer_key],
                np.asarray(values, dtype=object),
            )
            saw_topology = True
            continue
        if classification == "unsupported":
            raise BundleCapabilityError(
                f"save_sponge_input_bundle does not support serializer {serializer_key!r}"
            )
        if classification == "unknown" or legacy_key is None:
            raise BundleCapabilityError(
                f"save_sponge_input_bundle encountered unclassified serializer {serializer_key!r}"
            )
        if legacy_key == "coordinate_in_file":
            _add_coordinate_payload(builder, text)
            saw_coordinate = True
            continue

        contract = next(
            (
                item
                for item in contracts.get(legacy_key, ())
                if item.component == "topology" and item.reverse_policy != "alias"
            ),
            None,
        )
        if contract is None:
            raise BundleCapabilityError(
                f"serializer {serializer_key!r} has no bundled topology contract"
            )
        payload = MemoryPayload(f"{serializer_key}.txt", text)
        datasets = parse_topology_file(legacy_key, payload)  # type: ignore[arg-type]
        if datasets is None:
            raise BundleCapabilityError(
                f"serializer {serializer_key!r} has no typed bundled parser"
            )
        builder.add_datasets(contract.bundle_file, datasets)
        saw_topology = True

    if not saw_coordinate:
        raise BundleExportError("Molecule serializers did not produce a coordinate payload")
    if not saw_topology:
        raise BundleExportError("Molecule serializers did not produce any topology payload")


def _add_coordinate_payload(builder: BundleBuilder, text: str) -> None:
    lines = [line.split() for line in text.splitlines() if line.strip()]
    if not lines:
        raise BundleExportError("coordinate serializer returned an empty payload")
    atom_count = int(lines[0][0])
    if len(lines) != atom_count + 2:
        raise BundleExportError(
            f"coordinate serializer declared {atom_count} atoms but returned {len(lines) - 2}"
        )
    coordinate = np.asarray(
        [[float(value) for value in row[:3]] for row in lines[1 : atom_count + 1]],
        dtype=np.float32,
    )
    box = np.asarray([float(value) for value in lines[-1][:6]], dtype=np.float32)
    if coordinate.shape != (atom_count, 3) or box.size != 6:
        raise BundleExportError("coordinate serializer returned an invalid coordinate or box shape")
    builder.add_dataset(
        "restart.spgr.h5",
        "/particles/all/position/value",
        coordinate[np.newaxis, ...],
    )
    builder.add_dataset(
        "restart.spgr.h5",
        "/particles/all/box/edges/value",
        box_to_edges(box)[np.newaxis, ...],
    )
    builder.add_dataset(
        "restart.spgr.h5",
        "/particles/all/step",
        np.asarray([0], dtype=np.int64),
    )
    builder.add_dataset(
        "restart.spgr.h5",
        "/particles/all/time",
        np.asarray([0.0], dtype=np.float64),
    )


def _parse_counted_strings(key: str, text: str) -> list[str]:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise BundleExportError(f"{key} serializer returned an empty payload")
    count = int(lines[0])
    values = lines[1:]
    if len(values) != count:
        raise BundleExportError(
            f"{key} serializer declared {count} values but returned {len(values)}"
        )
    return values
