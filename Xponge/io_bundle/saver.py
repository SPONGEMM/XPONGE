"""Direct XPONGE object to SPONGE bundled input saver."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile
from uuid import NAMESPACE_URL, uuid5

import numpy as np

from Xponge.build import (
    _capture_sponge_source_atom_ids,
    _prepare_sponge_input_molecule,
    _sponge_save_result,
)

from .bundle_builder import BundleBuilder, BundleMetadata, BundlePaths
from .errors import BundlePathError
from .molecule_adapter import add_molecule_to_bundle
from .protocol import SpongeProtocol, add_protocol_to_bundle


def save_sponge_input_bundle(
    cls,
    prefix=None,
    dirname=".",
    *,
    simulation_mapping=None,
    topology_hash=None,
    atom_order_hash=None,
    mapping_hash=None,
    source_atom_ids=None,
    return_mapping=False,
    protocol: SpongeProtocol | None = None,
):
    """Save an XPONGE object as topology, protocol, and restart HDF5 inputs.

    The returned value matches :func:`Xponge.save_sponge_input`: the prepared
    ``Molecule`` instance used for serialization.
    """

    captured_source_atom_ids = _capture_sponge_source_atom_ids(cls, source_atom_ids)
    mol = _prepare_sponge_input_molecule(cls)
    prefix = prefix or mol.name
    output_root = Path(dirname).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    final_paths = _prefixed_bundle_paths(output_root, str(prefix))
    final_paths.topology.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".xponge-bundle-", dir=output_root) as tempdir:
        temporary_paths = BundlePaths.canonical(tempdir)
        builder = (
            BundleBuilder(
                temporary_paths,
                identity_uuid=str(uuid5(NAMESPACE_URL, f"xponge-simulation-mapping:{mapping_hash}")),
            )
            if mapping_hash else BundleBuilder(temporary_paths)
        )
        add_molecule_to_bundle(mol, builder)
        protocol_summary = add_protocol_to_bundle(
            protocol,
            builder,
            atom_count=len(mol.atoms),
        )
        if simulation_mapping is not None:
            if not topology_hash or not atom_order_hash or not mapping_hash:
                raise ValueError("simulation mapping requires topology, atom-order, and mapping hashes")
            records = tuple(simulation_mapping)
            if len(records) != len(mol.atoms):
                raise ValueError("simulation mapping must cover every serialized atom")
            builder.add_dataset(
                "topology.spgt.h5",
                "/topology/mapping/simulation_atom_index",
                np.asarray([record["simulation_index"] for record in records], dtype=np.int64),
            )
            builder.add_dataset(
                "topology.spgt.h5",
                "/topology/mapping/external_id",
                np.asarray([record["external_id"] for record in records], dtype=object),
            )
            builder.add_dataset(
                "topology.spgt.h5",
                "/topology/mapping/canonical_atom_id",
                np.asarray([record["canonical_atom_id"] for record in records], dtype=np.int64),
            )
            builder.add_dataset(
                "topology.spgt.h5",
                "/topology/mapping/simulation_residue_index",
                np.asarray([record["simulation_residue_index"] for record in records], dtype=np.int64),
            )
            builder.add_dataset(
                "topology.spgt.h5",
                "/topology/mapping/simulation_residue_id",
                np.asarray([record["simulation_residue_id"] for record in records], dtype=object),
            )
            builder.add_dataset(
                "topology.spgt.h5",
                "/topology/mapping/mapping_hash",
                np.asarray([mapping_hash], dtype=object),
            )
        metadata = BundleMetadata(
            atom_count=len(mol.atoms),
            topology_hash=topology_hash or builder.content_hash(bundle_file="topology.spgt.h5"),
            atom_order_hash=atom_order_hash or builder.content_hash(
                bundle_file="topology.spgt.h5", path_prefixes=("/atoms/", "/residues/"),
            ),
            forcefield_hash=builder.content_hash(
                bundle_file="topology.spgt.h5",
                path_prefixes=("/forcefield/", "/manybody/", "/qc/"),
            ),
            protocol_hash=builder.content_hash(bundle_file="protocol.spgp.h5"),
            cv_count=protocol_summary.cv_count,
            restraint_count=protocol_summary.restraint_count,
            enhanced_methods=protocol_summary.enhanced_methods,
            creator_version="xponge-bundled-saver",
        )
        touched = {"topology.spgt.h5", "protocol.spgp.h5", "restart.spgr.h5"}
        builder.finalize(touched, metadata)
        for source, target in (
            (temporary_paths.topology, final_paths.topology),
            (temporary_paths.protocol, final_paths.protocol),
            (temporary_paths.restart, final_paths.restart),
        ):
            target.parent.mkdir(parents=True, exist_ok=True)
            os.replace(source, target)
    return _sponge_save_result(mol, captured_source_atom_ids, return_mapping)


def _prefixed_bundle_paths(root: Path, prefix: str) -> BundlePaths:
    relative = Path(prefix)
    if relative.is_absolute():
        raise BundlePathError(f"bundled input prefix must be relative: {prefix}")
    base = (root / relative).resolve()
    try:
        base.relative_to(root)
    except ValueError as exc:
        raise BundlePathError(f"bundled input prefix escapes output directory: {prefix}") from exc
    return BundlePaths(
        topology=Path(str(base) + "_topology.spgt.h5"),
        protocol=Path(str(base) + "_protocol.spgp.h5"),
        restart=Path(str(base) + "_restart.spgr.h5"),
        trajectory=None,
    )
