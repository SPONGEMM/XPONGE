"""Isolated analytical-Hessian worker consuming one locked small model."""

from __future__ import annotations

from contextlib import redirect_stdout
import json
import math
import os
import sys
from typing import Any, Mapping

from ._ion_worker import _canonical_hash, _strict_object
from ._resp_worker import WorkerInputError, _finite_number, _package_version, _validate_model


WORKER_PROTOCOL_VERSION = 1


def _validate_fit_protocol(
    value: Any,
    elements: tuple[str, ...],
    metal_elements: tuple[str, ...],
) -> tuple[Mapping[str, Any], Any]:
    protocol = _strict_object(
        value,
        {
            "fit_input_hash", "backend", "method", "basis_family", "metal_basis_policy",
            "basis_source", "optimize_geometry", "scf_convergence_tolerance",
            "scf_max_cycles", "threads", "memory_limit_bytes", "source",
        },
        "fit_protocol",
    )
    if not isinstance(protocol["fit_input_hash"], str) or len(protocol["fit_input_hash"]) != 64:
        raise WorkerInputError("fit_protocol.fit_input_hash: expected sha256")
    if protocol["backend"] != "pyscf":
        raise WorkerInputError("fit_protocol.backend: analytical Hessian provider requires pyscf")
    if protocol["method"] != "scf":
        raise WorkerInputError("fit_protocol.method: expected scf")
    if not isinstance(protocol["basis_family"], str) or not protocol["basis_family"].strip():
        raise WorkerInputError("fit_protocol.basis_family: expected non-empty string")
    if protocol["metal_basis_policy"] != "require_ecp":
        raise WorkerInputError("fit_protocol.metal_basis_policy: expected require_ecp")
    if not isinstance(protocol["basis_source"], str) or not protocol["basis_source"].strip():
        raise WorkerInputError("fit_protocol.basis_source: expected non-empty string")
    if protocol["optimize_geometry"] is not False:
        raise WorkerInputError("fit_protocol.optimize_geometry: locked derived-model coordinates are required")
    if _finite_number(
        protocol["scf_convergence_tolerance"], "fit_protocol.scf_convergence_tolerance"
    ) <= 0:
        raise WorkerInputError("fit_protocol.scf_convergence_tolerance: expected positive value")
    for name in ("scf_max_cycles", "threads", "memory_limit_bytes"):
        number = protocol[name]
        if isinstance(number, bool) or not isinstance(number, int) or number <= 0:
            raise WorkerInputError(f"fit_protocol.{name}: expected positive integer")
    if not isinstance(protocol["source"], str) or not protocol["source"].strip():
        raise WorkerInputError("fit_protocol.source: expected non-empty string")

    from Xponge.qm.resp_basis import resolve_resp_basis

    try:
        resolved_basis = resolve_resp_basis(protocol["basis_family"], elements)
    except ValueError as exc:
        raise WorkerInputError(f"fit_protocol.basis_family: {exc}") from exc
    if metal_elements:
        missing = sorted(set(metal_elements) - set(resolved_basis.ecp or {}))
        if missing:
            raise WorkerInputError(
                "fit_protocol.metal_basis_policy: resolved ECP does not cover " + ",".join(missing)
            )
    return protocol, resolved_basis


def _configure_thread_environment(threads: int) -> None:
    value = str(threads)
    for name in (
        "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
    ):
        os.environ[name] = value


def _execute(value: Any) -> dict[str, Any]:
    request = _strict_object(value, {"protocol_version", "model", "fit_protocol"}, "request")
    if request["protocol_version"] != WORKER_PROTOCOL_VERSION:
        raise WorkerInputError("unsupported protocol_version")
    model, atom_ids, elements, metal_elements = _validate_model(
        request["model"], required_purpose="small"
    )
    protocol, resolved_basis = _validate_fit_protocol(
        request["fit_protocol"], elements, metal_elements
    )
    _configure_thread_environment(protocol["threads"])

    from Xponge.assign import Assign
    from Xponge.qm.scheduler import compute_hessian

    atom_index = {atom_id: index for index, atom_id in enumerate(atom_ids)}
    assign = Assign(name=model["external_id"])
    for atom in model["atoms"]:
        x, y, z = atom["coordinates"]
        assign.add_atom(
            atom["element"], x, y, z, name=atom["model_atom_id"],
            charge=0.0 if atom["initial_charge"] is None else atom["initial_charge"],
        )
    for edge in model["bonds"]:
        first, second = edge["model_atom_ids"]
        assign.add_bond(atom_index[first], atom_index[second], edge["order"])
    for edge in model["links"]:
        first, second = edge["model_atom_ids"]
        assign.add_bond(atom_index[first], atom_index[second], 1.0)

    state = model["electronic_state"]
    result = compute_hessian(
        assign,
        backend=protocol["backend"],
        basis=resolved_basis.basis,
        ecp=resolved_basis.ecp,
        cart=resolved_basis.cart,
        charge=state["net_charge"],
        spin=state["spin_multiplicity"] - 1,
        threads=protocol["threads"],
        memory_limit_bytes=protocol["memory_limit_bytes"],
        scf_convergence_tolerance=protocol["scf_convergence_tolerance"],
        scf_max_cycles=protocol["scf_max_cycles"],
        return_timings=True,
    )
    if result.scf_converged is not True:
        raise RuntimeError("QM backend returned an unconverged SCF Hessian")
    if list(result.atom_symbols) != list(elements):
        raise RuntimeError("QM backend changed the locked atom order")
    coordinates = [list(coordinate) for coordinate in result.coordinates_angstrom]
    expected_coordinates = [list(atom["coordinates"]) for atom in model["atoms"]]
    if len(coordinates) != len(expected_coordinates):
        raise RuntimeError("QM backend returned the wrong coordinate count")
    for index, (actual, expected) in enumerate(zip(coordinates, expected_coordinates)):
        if len(actual) != 3 or any(
            not math.isclose(float(left), float(right), rel_tol=0.0, abs_tol=1.0e-8)
            for left, right in zip(actual, expected)
        ):
            raise RuntimeError(f"QM backend changed locked coordinates at atom {index}")

    try:
        hessian = result.cartesian_hessian_au.tolist()
    except AttributeError:
        hessian = result.cartesian_hessian_au
    backend_version = _package_version(protocol["backend"])
    xponge_version = _package_version("Xponge")
    provider = f"xponge-hessian:{protocol['backend']}"
    provider_version = f"xponge={xponge_version};{protocol['backend']}={backend_version}"
    artifact = {
        "model_id": model["external_id"],
        "model_hash": model["model_hash"],
        "atom_order": list(atom_ids),
        "coordinates_angstrom": coordinates,
        "cartesian_hessian_au": hessian,
        "provider": provider,
        "provider_version": provider_version,
    }
    request_hash = _canonical_hash(request)
    fit_report = {
        "model_id": model["external_id"],
        "model_hash": model["model_hash"],
        "fit_input_hash": protocol["fit_input_hash"],
        "request_hash": request_hash,
        "backend": protocol["backend"],
        "backend_version": backend_version,
        "xponge_version": xponge_version,
        "method": protocol["method"],
        "basis_family": protocol["basis_family"],
        "metal_basis_policy": protocol["metal_basis_policy"],
        "basis_source": protocol["basis_source"],
        "metal_elements": list(metal_elements),
        "net_charge": state["net_charge"],
        "spin_multiplicity": state["spin_multiplicity"],
        "backend_spin_2s": state["spin_multiplicity"] - 1,
        "coordinate_unit": "angstrom",
        "geometry_locked": True,
        "atom_count": len(atom_ids),
        "setup_metadata": {
            "method": "SCF_HESSIAN",
            "scf_converged": True,
            "total_energy_hartree": result.total_energy,
            "basis_family": resolved_basis.label,
            "basis": resolved_basis.basis,
            "ecp": resolved_basis.ecp,
            "cart": resolved_basis.cart,
            "references": list(resolved_basis.references),
            "scf_convergence_tolerance": float(protocol["scf_convergence_tolerance"]),
            "scf_max_cycles": protocol["scf_max_cycles"],
            "threads": protocol["threads"],
            "memory_limit_bytes": protocol["memory_limit_bytes"],
            "timings": dict(result.timings),
        },
        "source": protocol["source"],
    }
    response = {
        "ok": True,
        "protocol_version": WORKER_PROTOCOL_VERSION,
        "model_id": model["external_id"],
        "model_hash": model["model_hash"],
        "fit_input_hash": protocol["fit_input_hash"],
        "request_hash": request_hash,
        "artifact": artifact,
        "fit_report": fit_report,
    }
    response["response_hash"] = _canonical_hash(response)
    return response


def main() -> int:
    try:
        value = json.load(sys.stdin)
        with redirect_stdout(sys.stderr):
            response = _execute(value)
        json.dump(response, sys.stdout, sort_keys=True, separators=(",", ":"), allow_nan=False)
        sys.stdout.write("\n")
        return 0
    except Exception as exc:
        json.dump(
            {
                "ok": False,
                "protocol_version": WORKER_PROTOCOL_VERSION,
                "error": {"type": type(exc).__name__, "message": str(exc)},
            },
            sys.stdout,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        sys.stdout.write("\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
