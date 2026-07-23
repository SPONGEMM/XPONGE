"""Direct XPONGE object to SPONGE bundled input saver."""

from __future__ import annotations

import os
from pathlib import Path
import tempfile

from Xponge.build import _prepare_sponge_input_molecule

from .bundle_builder import BundleBuilder, BundleMetadata, BundlePaths
from .errors import BundlePathError
from .molecule_adapter import add_molecule_to_bundle
from .protocol import SpongeProtocol, add_protocol_to_bundle


def save_sponge_input_bundle(
    cls,
    prefix=None,
    dirname=".",
    *,
    protocol: SpongeProtocol | None = None,
):
    """Save an XPONGE object as topology, protocol, and restart HDF5 inputs.

    The returned value matches :func:`Xponge.save_sponge_input`: the prepared
    ``Molecule`` instance used for serialization.
    """

    mol = _prepare_sponge_input_molecule(cls)
    prefix = prefix or mol.name
    output_root = Path(dirname).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    final_paths = _prefixed_bundle_paths(output_root, str(prefix))
    final_paths.topology.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=".xponge-bundle-", dir=output_root) as tempdir:
        temporary_paths = BundlePaths.canonical(tempdir)
        builder = BundleBuilder(temporary_paths)
        add_molecule_to_bundle(mol, builder)
        protocol_summary = add_protocol_to_bundle(
            protocol,
            builder,
            atom_count=len(mol.atoms),
        )
        metadata = BundleMetadata(
            atom_count=len(mol.atoms),
            topology_hash=builder.content_hash(bundle_file="topology.spgt.h5"),
            atom_order_hash=builder.content_hash(
                bundle_file="topology.spgt.h5",
                path_prefixes=("/atoms/", "/residues/"),
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
    return mol


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
