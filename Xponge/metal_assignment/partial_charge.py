"""Hash-closed composition of explicit partial-charge artifacts."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from .contracts import (
    ChargeOverlay,
    MetalAssignmentInput,
    PartialChargeArtifact,
    PartialChargeArtifactBundle,
    _canonicalize,
    validate_input,
    validate_partial_charge_artifacts,
)


def compose_partial_charge_artifacts(
    request: MetalAssignmentInput,
    additional_artifacts: Sequence[PartialChargeArtifact] = (),
) -> tuple[ChargeOverlay | None, Mapping[str, Any]]:
    """Resolve per-atom charges by explicit precedence without implicit fallback."""

    validate_input(request)
    base_artifacts = (
        request.partial_charge_artifacts.artifacts
        if request.partial_charge_artifacts is not None
        else ()
    )
    artifacts = tuple(base_artifacts) + tuple(additional_artifacts)
    if not artifacts:
        return None, {"status": "not_provided", "artifact_count": 0}
    combined = PartialChargeArtifactBundle(
        schema_version=request.schema_version,
        topology_hash=request.topology.topology_hash,
        graph_revision=request.graph_revision,
        input_hash=request.input_hash,
        artifacts=artifacts,
    ).with_computed_hash()
    validate_partial_charge_artifacts(request, combined)

    selected: dict[str, tuple[PartialChargeArtifact, float]] = {}
    overridden: list[dict[str, Any]] = []
    for artifact in artifacts:
        for atom_id, charge in zip(artifact.atom_ids, artifact.charges):
            previous = selected.get(atom_id)
            if previous is None or artifact.precedence > previous[0].precedence:
                if previous is not None:
                    overridden.append({
                        "atom_id": atom_id,
                        "superseded_artifact_id": previous[0].artifact_id,
                        "selected_artifact_id": artifact.artifact_id,
                    })
                selected[atom_id] = (artifact, float(charge))
    stable_order = {atom.external_id: atom.stable_order for atom in request.topology.atoms}
    atom_ids = tuple(sorted(selected, key=lambda atom_id: (stable_order[atom_id], atom_id)))
    charges = {atom_id: selected[atom_id][1] for atom_id in atom_ids}
    atom_sources = {
        atom_id: {
            "artifact_id": selected[atom_id][0].artifact_id,
            "artifact_hash": selected[atom_id][0].artifact_hash,
            "precedence": selected[atom_id][0].precedence,
            "scope_kind": selected[atom_id][0].scope_kind,
            "scope_id": selected[atom_id][0].scope_id,
            "atomic_charge_role": selected[atom_id][0].atomic_charge_role,
            "provider": selected[atom_id][0].provider,
            "provider_version": selected[atom_id][0].provider_version,
            "method": selected[atom_id][0].method,
            "source": selected[atom_id][0].source,
        }
        for atom_id in atom_ids
    }
    overlay = ChargeOverlay(
        topology_hash=request.topology.topology_hash,
        charges=charges,
        source="partial-charge-artifact-composition:v1",
        atom_sources=atom_sources,
        artifact_hashes=tuple(sorted(artifact.artifact_hash for artifact in artifacts)),
    ).with_computed_hash()
    return overlay, {
        "status": "composed",
        "bundle_hash": combined.bundle_hash,
        "artifact_count": len(artifacts),
        "artifact_hashes": overlay.artifact_hashes,
        "covered_atom_ids": atom_ids,
        "atom_sources": _canonicalize(atom_sources),
        "overridden": tuple(overridden),
    }


__all__ = ["compose_partial_charge_artifacts"]
