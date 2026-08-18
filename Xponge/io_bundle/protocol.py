"""Native object model for ``sponge.input.v2`` protocol bundles.

These objects describe protocol semantics directly.  They intentionally do not
compile legacy SPONGE configuration text before writing HDF5.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import re
from typing import Any, Mapping, Sequence

import numpy as np

from .bundle_builder import BundleBuilder
from .errors import BundleValidationError


_OBJECT_NAME = re.compile(r"^[A-Za-z_][A-Za-z0-9_.-]*$")
_CV_RESERVED_PARAMETERS = frozenset(
    {
        "CV_type",
        "atom",
        "coordinate",
        "period",
        "sigma",
        "rotate",
        "function",
        "min_padding",
        "max_padding",
    }
)
_METADYNAMICS_FLAGS = frozenset(
    {"subhill", "kde", "mask", "sink", "convmeta", "grw", "dip"}
)


@dataclass(frozen=True)
class ProtocolCollectiveVariable:
    """One named scalar collective variable."""

    name: str
    type: str
    atom_indices: tuple[int, ...] = ()
    atom_refs: tuple[int | str, ...] = ()
    parameters: Mapping[str, Any] = field(default_factory=dict)
    period: tuple[float, ...] = ()
    sigma: tuple[float, ...] = ()
    reference_coordinates: tuple[tuple[float, float, float], ...] = ()
    enabled: bool = True
    dimension: int = 1
    rotate: bool | None = None
    function: str | None = None
    min_padding: float | None = None
    max_padding: float | None = None


@dataclass(frozen=True)
class ProtocolVirtualAtom:
    """One named center used as an atom by collective variables."""

    name: str
    type: str
    atom_indices: tuple[int, ...]
    weight: tuple[float, ...] = ()
    enabled: bool = True


@dataclass(frozen=True)
class ProtocolDistanceConstraints:
    """Extra distance-constraint pairs used by the runtime constraint module."""

    atoms: tuple[tuple[int, int], ...]
    r0: tuple[float, ...]
    name: str = "default"


@dataclass(frozen=True)
class ProtocolPositionalRestraint:
    """One harmonic positional restraint with a full-system launch reference."""

    name: str
    atom_indices: tuple[int, ...]
    reference_coordinates: tuple[tuple[float, float, float], ...]
    weight: tuple[tuple[float, float, float], ...] = ()
    single_weight_default: float | None = None
    refcoord_scaling_default: str | None = None
    calc_virial_default: bool | None = None
    enabled: bool = True


@dataclass(frozen=True)
class ProtocolCVRestraint:
    """A harmonic restraint over named scalar collective variables."""

    name: str
    cv_refs: tuple[str, ...]
    weight: tuple[float, ...]
    reference: tuple[float, ...]
    period: tuple[float, ...] = ()
    start_step: tuple[int, ...] = ()
    max_step: tuple[int, ...] = ()
    reduce_step: tuple[int, ...] = ()
    stop_step: tuple[int, ...] = ()
    enabled: bool = True


@dataclass(frozen=True)
class ProtocolMetadynamics:
    """One named metadynamics bias over scalar collective variables."""

    name: str
    cv_refs: tuple[str, ...]
    grid_min: tuple[float, ...]
    grid_max: tuple[float, ...]
    grid_count: tuple[int, ...]
    sigma: tuple[float, ...] = ()
    period: tuple[float, ...] = ()
    cutoff: tuple[float, ...] = ()
    hill_height_default: float | None = None
    sumhill_freq_default: int | None = None
    potential_update_interval_default: int | None = None
    well_tempered_factor_default: float | None = None
    wall_height_default: float | None = None
    max_force_default: float | None = None
    method_flags: Mapping[str, int | float | bool] = field(default_factory=dict)
    enabled: bool = True


@dataclass(frozen=True)
class ProtocolSteering:
    """Static linear steering bias over named scalar collective variables."""

    cv_refs: tuple[str, ...]
    weight: tuple[float, ...]
    enabled: bool = True


@dataclass(frozen=True)
class ProtocolSoftWall:
    """One named external coordinate potential compiled by SPONGE."""

    name: str
    potential: str


@dataclass(frozen=True)
class ProtocolHardWall:
    """Axis-aligned reflection bounds owned by the simulation protocol."""

    bounds_low: tuple[float | None, float | None, float | None] = (
        None,
        None,
        None,
    )
    bounds_high: tuple[float | None, float | None, float | None] = (
        None,
        None,
        None,
    )
    allow_npt: bool = False


@dataclass(frozen=True)
class ProtocolSITS:
    """One native SITS enhanced-sampling protocol and its launch state."""

    mode: str
    atom_indices: tuple[int, ...] = ()
    atom_numbers_policy: str | int | None = None
    k_numbers: int | None = None
    temperature_ladder: tuple[float, ...] = ()
    temperature_low: float | None = None
    temperature_high: float | None = None
    pe_a: float | None = None
    pe_b: float | None = None
    fb_interval: int | None = None
    fb_bias: float | None = None
    record_interval: int | None = None
    update_interval: int | None = None
    nk_rest: bool | None = None
    nk_fix: bool | None = None
    cross_enhance_factor: float | None = None
    initial_nk: tuple[float, ...] = ()
    initial_log_norm: tuple[float, ...] = ()
    initial_log_nk: tuple[float, ...] = ()
    enabled: bool = True


@dataclass(frozen=True)
class SpongeProtocol:
    """Composable native protocol attached to a serialized XPONGE system."""

    collective_variables: tuple[ProtocolCollectiveVariable, ...] = ()
    virtual_atoms: tuple[ProtocolVirtualAtom, ...] = ()
    distance_constraints: tuple[ProtocolDistanceConstraints, ...] = ()
    positional_restraints: tuple[ProtocolPositionalRestraint, ...] = ()
    cv_restraints: tuple[ProtocolCVRestraint, ...] = ()
    metadynamics: tuple[ProtocolMetadynamics, ...] = ()
    steering: ProtocolSteering | None = None
    sits: ProtocolSITS | None = None
    hard_wall: ProtocolHardWall | None = None
    soft_walls: tuple[ProtocolSoftWall, ...] = ()


@dataclass(frozen=True)
class ProtocolSummary:
    """Metadata derived while serializing a native protocol."""

    cv_count: int = 0
    restraint_count: int = 0
    enhanced_methods: tuple[str, ...] = ()


def add_protocol_to_bundle(
    protocol: SpongeProtocol | None,
    builder: BundleBuilder,
    *,
    atom_count: int,
) -> ProtocolSummary:
    """Validate and serialize a protocol directly into canonical HDF5 groups."""

    if protocol is None:
        return ProtocolSummary()
    if not isinstance(protocol, SpongeProtocol):
        raise TypeError(
            "protocol must be a SpongeProtocol instance or None, not "
            f"{type(protocol).__name__}"
        )

    cvs = tuple(protocol.collective_variables)
    virtual_atoms = tuple(protocol.virtual_atoms)
    constraints = tuple(protocol.distance_constraints)
    positional = tuple(protocol.positional_restraints)
    cv_restraints = tuple(protocol.cv_restraints)
    metadynamics = tuple(protocol.metadynamics)
    steering = protocol.steering
    sits = protocol.sits
    hard_wall = protocol.hard_wall
    soft_walls = tuple(protocol.soft_walls)
    _validate_protocol(
        cvs,
        virtual_atoms,
        constraints,
        positional,
        cv_restraints,
        metadynamics,
        steering,
        sits,
        hard_wall,
        soft_walls,
        atom_count=atom_count,
    )

    for virtual_atom in virtual_atoms:
        _write_virtual_atom(builder, virtual_atom)
    for cv in cvs:
        _write_cv(builder, cv)
    for constraint in constraints:
        _write_constraint(builder, constraint)
    for restraint in positional:
        _write_positional_restraint(builder, restraint)
    for restraint in cv_restraints:
        _write_cv_restraint(builder, restraint)
    for bias in metadynamics:
        _write_metadynamics(builder, bias)
    if steering is not None:
        _write_steering(builder, steering)
    if sits is not None:
        _write_sits(builder, sits)
    if hard_wall is not None:
        _write_hard_wall(builder, hard_wall)
    if soft_walls:
        _write_soft_walls(builder, soft_walls)

    enhanced_methods = []
    if any(item.enabled for item in metadynamics):
        enhanced_methods.append("metadynamics")
    if sits is not None and sits.enabled:
        enhanced_methods.append("SITS")
    if soft_walls:
        enhanced_methods.append("soft_walls")
    return ProtocolSummary(
        cv_count=sum(cv.enabled for cv in cvs),
        restraint_count=sum(item.enabled for item in positional + cv_restraints),
        enhanced_methods=tuple(enhanced_methods),
    )


def _validate_protocol(
    cvs: tuple[ProtocolCollectiveVariable, ...],
    virtual_atoms: tuple[ProtocolVirtualAtom, ...],
    constraints: tuple[ProtocolDistanceConstraints, ...],
    positional: tuple[ProtocolPositionalRestraint, ...],
    cv_restraints: tuple[ProtocolCVRestraint, ...],
    metadynamics: tuple[ProtocolMetadynamics, ...],
    steering: ProtocolSteering | None,
    sits: ProtocolSITS | None,
    hard_wall: ProtocolHardWall | None,
    soft_walls: tuple[ProtocolSoftWall, ...],
    *,
    atom_count: int,
) -> None:
    if atom_count <= 0:
        raise BundleValidationError("a native protocol requires a positive atom count")
    _require_unique_names(cvs, "collective variable")
    _require_unique_names(virtual_atoms, "virtual atom")
    _require_unique_names(constraints, "distance constraint")
    _require_unique_names(positional + cv_restraints, "restraint")
    _require_unique_names(metadynamics, "metadynamics")

    cv_by_name = {cv.name: cv for cv in cvs if cv.enabled}
    virtual_atom_names = {item.name for item in virtual_atoms if item.enabled}
    overlapping_names = {cv.name for cv in cvs} & {item.name for item in virtual_atoms}
    if overlapping_names:
        raise BundleValidationError(
            f"collective variables and virtual atoms must use distinct names: {sorted(overlapping_names)}"
        )
    for item in virtual_atoms:
        _validate_name(item.name, "virtual atom")
        if item.type not in {"center", "center_of_mass"}:
            raise BundleValidationError(
                f"virtual atom {item.name!r} type must be 'center' or 'center_of_mass'"
            )
        if not item.atom_indices:
            raise BundleValidationError(f"virtual atom {item.name!r} requires atom_indices")
        _validate_atom_indices(item.atom_indices, atom_count, f"virtual atom {item.name!r}")
        if item.type == "center":
            _validate_vector_length(
                item.weight,
                len(item.atom_indices),
                f"virtual atom {item.name!r} weight",
                required=True,
            )
            _require_finite(item.weight, f"virtual atom {item.name!r} weight")
        elif item.weight:
            raise BundleValidationError(
                f"virtual atom {item.name!r} center_of_mass must not define weight"
            )
    for cv in cvs:
        _validate_name(cv.name, "collective variable")
        if cv.name == "virtual_atom":
            raise BundleValidationError(
                "collective variable name 'virtual_atom' is reserved by /cv/virtual_atom"
            )
        _validate_name(cv.type, f"collective variable {cv.name!r} type")
        if cv.dimension != 1:
            raise BundleValidationError(
                f"collective variable {cv.name!r} dimension must be 1 for the current SPONGE runtime"
            )
        if cv.atom_indices and cv.atom_refs:
            raise BundleValidationError(
                f"collective variable {cv.name!r} atom_indices and atom_refs are mutually exclusive"
            )
        _validate_atom_indices(cv.atom_indices, atom_count, f"collective variable {cv.name!r}")
        physical_atom_refs = [item for item in cv.atom_refs if isinstance(item, (int, np.integer))]
        invalid_atom_refs = [item for item in cv.atom_refs if not isinstance(item, (str, int, np.integer))]
        if invalid_atom_refs:
            raise BundleValidationError(
                f"collective variable {cv.name!r} atom_refs must contain atom indices or virtual atom names"
            )
        normalized_atom_refs = tuple(str(item) for item in cv.atom_refs)
        if len(normalized_atom_refs) != len(set(normalized_atom_refs)):
            raise BundleValidationError(
                f"collective variable {cv.name!r} atom_refs must be unique"
            )
        _validate_atom_indices(physical_atom_refs, atom_count, f"collective variable {cv.name!r}")
        missing_atom_refs = [
            name for name in cv.atom_refs if isinstance(name, str) and name not in virtual_atom_names
        ]
        if missing_atom_refs:
            raise BundleValidationError(
                f"collective variable {cv.name!r} references missing or disabled virtual atoms: "
                f"{missing_atom_refs}"
            )
        _validate_vector_length(cv.period, cv.dimension, f"collective variable {cv.name!r} period")
        _validate_vector_length(cv.sigma, cv.dimension, f"collective variable {cv.name!r} sigma")
        if cv.sigma and any(value <= 0 or not np.isfinite(value) for value in cv.sigma):
            raise BundleValidationError(
                f"collective variable {cv.name!r} sigma values must be finite and positive"
            )
        if cv.reference_coordinates:
            _validate_xyz(
                cv.reference_coordinates,
                len(cv.atom_indices or cv.atom_refs),
                f"collective variable {cv.name!r} reference coordinates",
            )
        for key, value in cv.parameters.items():
            _validate_name(key, f"collective variable {cv.name!r} parameter")
            if key in _CV_RESERVED_PARAMETERS:
                raise BundleValidationError(
                    f"collective variable {cv.name!r} parameter {key!r} is a reserved native field"
                )
            _normalize_parameter(value, f"collective variable {cv.name!r} parameter {key!r}")

    for constraint in constraints:
        _validate_name(constraint.name, "distance constraint")
        if constraint.name != "default":
            raise BundleValidationError(
                "sponge.input.v2 currently supports only the 'default' distance-constraint object"
            )
        if not constraint.atoms or len(constraint.atoms) != len(constraint.r0):
            raise BundleValidationError(
                "distance constraint atoms and r0 must have the same non-zero length"
            )
        for pair in constraint.atoms:
            if len(pair) != 2 or pair[0] == pair[1]:
                raise BundleValidationError("distance constraint atom pairs must contain two distinct atoms")
            _validate_atom_indices(pair, atom_count, "distance constraint")
        _require_finite(constraint.r0, "distance constraint r0", positive=True)

    for restraint in positional:
        _validate_name(restraint.name, "positional restraint")
        if not restraint.atom_indices:
            raise BundleValidationError(
                f"positional restraint {restraint.name!r} requires atom_indices"
            )
        _validate_atom_indices(
            restraint.atom_indices, atom_count, f"positional restraint {restraint.name!r}"
        )
        _validate_xyz(
            restraint.reference_coordinates,
            atom_count,
            f"positional restraint {restraint.name!r} reference coordinates",
        )
        if restraint.weight:
            _validate_xyz(
                restraint.weight,
                len(restraint.atom_indices),
                f"positional restraint {restraint.name!r} weight",
                nonnegative=True,
            )
        elif restraint.single_weight_default is None:
            raise BundleValidationError(
                f"positional restraint {restraint.name!r} requires weight or single_weight_default"
            )
        if restraint.single_weight_default is not None and (
            not np.isfinite(restraint.single_weight_default)
            or restraint.single_weight_default < 0
        ):
            raise BundleValidationError(
                f"positional restraint {restraint.name!r} single_weight_default must be finite and non-negative"
            )
        if restraint.refcoord_scaling_default not in {
            None,
            "no",
            "none",
            "all",
            "com_ug",
            "com_res",
            "com_mol",
        }:
            raise BundleValidationError(
                f"positional restraint {restraint.name!r} refcoord_scaling_default "
                "must be no/all/com_ug/com_res/com_mol"
            )

    for restraint in cv_restraints:
        _validate_name(restraint.name, "CV restraint")
        count = len(restraint.cv_refs)
        if count == 0 or len(set(restraint.cv_refs)) != count:
            raise BundleValidationError(
                f"CV restraint {restraint.name!r} requires unique non-empty cv_refs"
            )
        _validate_cv_refs(restraint.cv_refs, cv_by_name, f"CV restraint {restraint.name!r}")
        for label, values in (
            ("weight", restraint.weight),
            ("reference", restraint.reference),
        ):
            _validate_vector_length(values, count, f"CV restraint {restraint.name!r} {label}", required=True)
            _require_finite(values, f"CV restraint {restraint.name!r} {label}")
        _validate_vector_length(restraint.period, count, f"CV restraint {restraint.name!r} period")
        for label, values in (
            ("start_step", restraint.start_step),
            ("max_step", restraint.max_step),
            ("reduce_step", restraint.reduce_step),
            ("stop_step", restraint.stop_step),
        ):
            _validate_vector_length(values, count, f"CV restraint {restraint.name!r} {label}")
            if any(value < 0 for value in values):
                raise BundleValidationError(
                    f"CV restraint {restraint.name!r} {label} values must be non-negative"
                )

    if sum(item.enabled for item in positional) > 1:
        raise BundleValidationError(
            "the current SPONGE runtime supports at most one enabled positional restraint"
        )
    if sum(item.enabled for item in metadynamics) > 1:
        raise BundleValidationError(
            "the current SPONGE runtime supports at most one enabled metadynamics object"
        )
    for bias in metadynamics:
        _validate_name(bias.name, "metadynamics")
        count = len(bias.cv_refs)
        if count == 0 or len(set(bias.cv_refs)) != count:
            raise BundleValidationError(
                f"metadynamics {bias.name!r} requires unique non-empty cv_refs"
            )
        _validate_cv_refs(bias.cv_refs, cv_by_name, f"metadynamics {bias.name!r}")
        for label, values in (
            ("grid_min", bias.grid_min),
            ("grid_max", bias.grid_max),
            ("grid_count", bias.grid_count),
        ):
            _validate_vector_length(values, count, f"metadynamics {bias.name!r} {label}", required=True)
        _require_finite(bias.grid_min, f"metadynamics {bias.name!r} grid_min")
        _require_finite(bias.grid_max, f"metadynamics {bias.name!r} grid_max")
        for minimum, maximum in zip(bias.grid_min, bias.grid_max):
            if maximum <= minimum:
                raise BundleValidationError(
                    f"metadynamics {bias.name!r} grid_max must exceed grid_min"
                )
        if any(value <= 1 for value in bias.grid_count):
            raise BundleValidationError(
                f"metadynamics {bias.name!r} grid_count values must be greater than 1"
            )
        for label, values, positive in (
            ("sigma", bias.sigma, True),
            ("period", bias.period, False),
            ("cutoff", bias.cutoff, True),
        ):
            _validate_vector_length(values, count, f"metadynamics {bias.name!r} {label}")
            _require_finite(values, f"metadynamics {bias.name!r} {label}", positive=positive)
        if not bias.sigma:
            missing = [name for name in bias.cv_refs if not cv_by_name[name].sigma]
            if missing:
                raise BundleValidationError(
                    f"metadynamics {bias.name!r} has no sigma and referenced CVs lack sigma: {missing}"
                )
        for key, value in bias.method_flags.items():
            if key not in _METADYNAMICS_FLAGS:
                raise BundleValidationError(
                    f"metadynamics {bias.name!r} has unsupported method flag {key!r}"
                )
            if not np.isfinite(value):
                raise BundleValidationError(
                    f"metadynamics {bias.name!r} method flag {key!r} must be finite"
                )
        if (
            bias.well_tempered_factor_default is not None
            and bias.well_tempered_factor_default <= 1
        ):
            raise BundleValidationError(
                f"metadynamics {bias.name!r} well_tempered_factor_default must be greater than 1"
            )
    if steering is not None:
        count = len(steering.cv_refs)
        if count == 0 or len(set(steering.cv_refs)) != count:
            raise BundleValidationError("steering requires unique non-empty cv_refs")
        _validate_cv_refs(steering.cv_refs, cv_by_name, "steering")
        _validate_vector_length(steering.weight, count, "steering weight", required=True)
        _require_finite(steering.weight, "steering weight")
    if sits is not None:
        _validate_sits(sits, atom_count)
    if hard_wall is not None:
        _validate_hard_wall(hard_wall)
    _require_unique_names(soft_walls, "soft wall")
    for wall in soft_walls:
        _validate_name(wall.name, "soft wall")
        if not wall.potential or "\0" in wall.potential:
            raise BundleValidationError(
                f"soft wall {wall.name!r} potential must be non-empty and contain no null byte"
            )
        for line in wall.potential.splitlines():
            if line.lstrip().startswith("[["):
                raise BundleValidationError(
                    f"soft wall {wall.name!r} potential contains a configuration delimiter"
                )


def _validate_sits(sits: ProtocolSITS, atom_count: int) -> None:
    allowed_modes = {"observation", "iteration", "production", "empirical", "amd", "gamd"}
    if sits.mode not in allowed_modes:
        raise BundleValidationError(
            f"SITS mode must be one of {sorted(allowed_modes)}, got {sits.mode!r}"
        )
    has_indices = bool(sits.atom_indices)
    has_policy = sits.atom_numbers_policy is not None
    if has_indices == has_policy:
        raise BundleValidationError(
            "SITS requires exactly one of atom_indices or atom_numbers_policy"
        )
    if has_indices:
        _validate_atom_indices(sits.atom_indices, atom_count, "SITS")
    else:
        policy = sits.atom_numbers_policy
        if isinstance(policy, str):
            if policy.upper() not in {"ITS", "ALL"}:
                raise BundleValidationError("SITS atom_numbers_policy must be ITS, ALL, or an integer")
        elif (
            isinstance(policy, (bool, np.bool_))
            or not isinstance(policy, (int, np.integer))
            or policy <= 0
            or policy > atom_count
        ):
            raise BundleValidationError(
                f"SITS integer atom_numbers_policy must be in [1, {atom_count}]"
            )

    for label, value in (
        ("temperature_low", sits.temperature_low),
        ("temperature_high", sits.temperature_high),
        ("pe_a", sits.pe_a),
        ("pe_b", sits.pe_b),
        ("fb_bias", sits.fb_bias),
        ("cross_enhance_factor", sits.cross_enhance_factor),
    ):
        if value is not None and not np.isfinite(value):
            raise BundleValidationError(f"SITS {label} must be finite")
    for label, value in (
        ("fb_interval", sits.fb_interval),
        ("record_interval", sits.record_interval),
        ("update_interval", sits.update_interval),
    ):
        if value is not None and (not isinstance(value, (int, np.integer)) or value <= 0):
            raise BundleValidationError(f"SITS {label} must be a positive integer")
    for label, value in (("nk_rest", sits.nk_rest), ("nk_fix", sits.nk_fix)):
        if value is not None and not isinstance(value, (bool, np.bool_)):
            raise BundleValidationError(f"SITS {label} must be boolean")

    if sits.mode in {"amd", "gamd"} and (sits.pe_a is None or sits.pe_b is None):
        raise BundleValidationError(f"SITS {sits.mode} mode requires pe_a and pe_b")
    if sits.mode == "empirical":
        if sits.temperature_low is None or sits.temperature_high is None:
            raise BundleValidationError("SITS empirical mode requires temperature_low and temperature_high")
        if sits.temperature_low <= 0 or sits.temperature_high <= sits.temperature_low:
            raise BundleValidationError(
                "SITS empirical temperatures must be positive with temperature_high > temperature_low"
            )

    adaptive = sits.mode in {"iteration", "production"}
    if adaptive:
        has_ladder = bool(sits.temperature_ladder)
        has_bounds = sits.temperature_low is not None or sits.temperature_high is not None
        if has_ladder == has_bounds:
            raise BundleValidationError(
                "SITS iteration/production requires exactly one temperature_ladder or low/high pair"
            )
        if has_bounds and (sits.temperature_low is None or sits.temperature_high is None):
            raise BundleValidationError("SITS temperature_low and temperature_high must be provided together")
        if has_ladder:
            _require_finite(sits.temperature_ladder, "SITS temperature_ladder", positive=True)
        elif sits.temperature_low <= 0 or sits.temperature_high <= sits.temperature_low:
            raise BundleValidationError(
                "SITS temperatures must be positive with temperature_high > temperature_low"
            )

        effective_k = sits.k_numbers or (len(sits.temperature_ladder) if has_ladder else 40)
        if (
            isinstance(effective_k, (bool, np.bool_))
            or not isinstance(effective_k, (int, np.integer))
            or effective_k <= 1
        ):
            raise BundleValidationError("SITS k_numbers must be an integer greater than 1")
        if has_ladder and len(sits.temperature_ladder) != effective_k:
            raise BundleValidationError("SITS temperature_ladder length must match k_numbers")
        for label, values in (
            ("initial_nk", sits.initial_nk),
            ("initial_log_norm", sits.initial_log_norm),
            ("initial_log_nk", sits.initial_log_nk),
        ):
            if values and len(values) != effective_k:
                raise BundleValidationError(f"SITS {label} length must match k_numbers")
        _require_finite(sits.initial_nk, "SITS initial_nk", positive=True)
        _require_finite(sits.initial_log_norm, "SITS initial_log_norm")
        _require_finite(sits.initial_log_nk, "SITS initial_log_nk")
        needs_restart = sits.nk_rest if sits.nk_rest is not None else sits.mode == "production"
        if needs_restart and not sits.initial_nk:
            raise BundleValidationError(
                "SITS production/nk_rest protocol requires initial_nk in the restart bundle"
            )
    else:
        if (
            sits.k_numbers is not None
            or sits.temperature_ladder
            or sits.initial_nk
            or sits.initial_log_norm
            or sits.initial_log_nk
        ):
            raise BundleValidationError(
                f"SITS {sits.mode} mode does not accept adaptive ladder or Nk state"
            )


def _validate_hard_wall(hard_wall: ProtocolHardWall) -> None:
    if not isinstance(hard_wall.allow_npt, (bool, np.bool_)):
        raise BundleValidationError("hard wall allow_npt must be boolean")
    low = _hard_wall_bounds(hard_wall.bounds_low, low=True)
    high = _hard_wall_bounds(hard_wall.bounds_high, low=False)
    if not np.any(np.isfinite(low)) and not np.any(np.isfinite(high)):
        raise BundleValidationError("hard wall requires at least one finite bound")
    invalid = np.flatnonzero(low >= high)
    if invalid.size:
        axis = "xyz"[int(invalid[0])]
        raise BundleValidationError(
            f"hard wall {axis}_low must be smaller than {axis}_high"
        )


def _hard_wall_bounds(values: Sequence[float | None], *, low: bool) -> np.ndarray:
    if len(values) != 3:
        label = "bounds_low" if low else "bounds_high"
        raise BundleValidationError(f"hard wall {label} must contain exactly three values")
    default = -np.inf if low else np.inf
    normalized = []
    for value in values:
        if value is None:
            normalized.append(default)
            continue
        if isinstance(value, (bool, np.bool_)):
            raise BundleValidationError("hard wall bounds must be numeric or None")
        try:
            number = float(value)
        except (TypeError, ValueError) as exc:
            raise BundleValidationError(
                "hard wall bounds must be numeric or None"
            ) from exc
        if np.isnan(number):
            raise BundleValidationError("hard wall bounds must not contain NaN")
        normalized.append(number)
    return np.asarray(normalized, dtype=np.float32)


def _write_cv(builder: BundleBuilder, cv: ProtocolCollectiveVariable) -> None:
    root = f"/cv/{cv.name}"
    _add_scalar(builder, root + "/type", cv.type)
    _add_scalar(builder, root + "/dimension", cv.dimension, np.int64)
    _add_scalar(builder, root + "/enabled_default", int(cv.enabled), np.int32)
    if cv.atom_indices:
        _add_array(builder, root + "/atom_indices", cv.atom_indices, np.int32)
    if cv.atom_refs:
        _add_array(builder, root + "/atom_refs", tuple(str(item) for item in cv.atom_refs), object)
    if cv.period:
        _add_array(builder, root + "/period", cv.period, np.float32)
    if cv.sigma:
        _add_array(builder, root + "/sigma", cv.sigma, np.float32)
    for name, value in sorted(cv.parameters.items()):
        builder.add_dataset(
            "protocol.spgp.h5",
            root + "/parameter/" + name,
            _normalize_parameter(value, f"collective variable {cv.name!r} parameter {name!r}"),
        )
    for field_name, value, dtype in (
        ("rotate", cv.rotate, np.int32),
        ("function", cv.function, None),
        ("min_padding", cv.min_padding, np.float32),
        ("max_padding", cv.max_padding, np.float32),
    ):
        if value is not None:
            _add_scalar(builder, root + "/" + field_name, value, dtype)
    if cv.reference_coordinates:
        _add_array(
            builder,
            f"/parameters/restart/references/cv/{cv.name}/coordinate",
            cv.reference_coordinates,
            np.float32,
            bundle_file="restart.spgr.h5",
        )


def _write_virtual_atom(builder: BundleBuilder, virtual_atom: ProtocolVirtualAtom) -> None:
    root = f"/cv/virtual_atom/{virtual_atom.name}"
    _add_scalar(builder, root + "/type", virtual_atom.type)
    _add_scalar(builder, root + "/enabled_default", int(virtual_atom.enabled), np.int32)
    _add_array(builder, root + "/atom_indices", virtual_atom.atom_indices, np.int32)
    if virtual_atom.weight:
        _add_array(builder, root + "/weight", virtual_atom.weight, np.float32)


def _write_constraint(builder: BundleBuilder, constraint: ProtocolDistanceConstraints) -> None:
    root = f"/constraint/{constraint.name}"
    _add_scalar(builder, root + "/schema_version", 1, np.int64)
    _add_array(builder, root + "/pairs/atoms", constraint.atoms, np.int32)
    _add_array(builder, root + "/pairs/r0", constraint.r0, np.float32)


def _write_positional_restraint(
    builder: BundleBuilder, restraint: ProtocolPositionalRestraint
) -> None:
    root = f"/restraint/{restraint.name}"
    _add_scalar(builder, root + "/type", "harmonic_positional")
    _add_scalar(builder, root + "/schema_version", 1, np.int64)
    _add_scalar(builder, root + "/enabled_default", int(restraint.enabled), np.int32)
    _add_array(builder, root + "/atom_indices", restraint.atom_indices, np.int32)
    if restraint.weight:
        _add_array(builder, root + "/weight", restraint.weight, np.float32)
    if restraint.single_weight_default is not None:
        _add_scalar(
            builder,
            root + "/single_weight_default",
            restraint.single_weight_default,
            np.float32,
        )
    if restraint.refcoord_scaling_default is not None:
        _add_scalar(
            builder,
            root + "/refcoord_scaling_default",
            "no"
            if restraint.refcoord_scaling_default == "none"
            else restraint.refcoord_scaling_default,
        )
    if restraint.calc_virial_default is not None:
        _add_scalar(
            builder,
            root + "/calc_virial_default",
            int(restraint.calc_virial_default),
            np.int32,
        )
    _add_array(
        builder,
        f"/parameters/restart/references/restraint/{restraint.name}/coordinate",
        restraint.reference_coordinates,
        np.float32,
        bundle_file="restart.spgr.h5",
    )


def _write_cv_restraint(builder: BundleBuilder, restraint: ProtocolCVRestraint) -> None:
    root = f"/restraint/{restraint.name}"
    _add_scalar(builder, root + "/type", "cv_harmonic")
    _add_scalar(builder, root + "/schema_version", 1, np.int64)
    _add_scalar(builder, root + "/enabled_default", int(restraint.enabled), np.int32)
    _add_array(builder, root + "/cv_refs", restraint.cv_refs, object)
    _add_array(builder, root + "/weight", restraint.weight, np.float32)
    _add_array(builder, root + "/reference", restraint.reference, np.float32)
    for name, values in (
        ("period", restraint.period),
        ("schedule/start_step", restraint.start_step),
        ("schedule/max_step", restraint.max_step),
        ("schedule/reduce_step", restraint.reduce_step),
        ("schedule/stop_step", restraint.stop_step),
    ):
        if values:
            dtype = np.float32 if name == "period" else np.int64
            _add_array(builder, root + "/" + name, values, dtype)


def _write_metadynamics(builder: BundleBuilder, bias: ProtocolMetadynamics) -> None:
    root = f"/meta/{bias.name}"
    _add_scalar(builder, root + "/schema_version", 1, np.int64)
    _add_scalar(builder, root + "/enabled_default", int(bias.enabled), np.int32)
    _add_array(builder, root + "/cv_refs", bias.cv_refs, object)
    _add_scalar(builder, root + "/ndim", len(bias.cv_refs), np.int64)
    _add_array(builder, root + "/grid/min", bias.grid_min, np.float32)
    _add_array(builder, root + "/grid/max", bias.grid_max, np.float32)
    _add_array(builder, root + "/grid/count", bias.grid_count, np.int64)
    for name, values in (
        ("sigma", bias.sigma),
        ("period", bias.period),
        ("cutoff", bias.cutoff),
    ):
        if values:
            _add_array(builder, root + "/" + name, values, np.float32)
    for name, value, dtype in (
        ("hill_height_default", bias.hill_height_default, np.float32),
        ("sumhill_freq_default", bias.sumhill_freq_default, np.int64),
        (
            "potential_update_interval_default",
            bias.potential_update_interval_default,
            np.int64,
        ),
        (
            "well_tempered_factor_default",
            bias.well_tempered_factor_default,
            np.float32,
        ),
        ("wall_height_default", bias.wall_height_default, np.float32),
        ("max_force_default", bias.max_force_default, np.float32),
    ):
        if value is not None:
            _add_scalar(builder, root + "/" + name, value, dtype)
    for name, value in sorted(bias.method_flags.items()):
        dtype = np.float32 if name == "dip" else np.int64
        _add_scalar(builder, root + "/method_flags/" + name, value, dtype)


def _write_steering(builder: BundleBuilder, steering: ProtocolSteering) -> None:
    _add_scalar(builder, "/steer/enabled_default", int(steering.enabled), np.int32)
    _add_array(builder, "/steer/cv_refs", steering.cv_refs, object)
    _add_array(builder, "/steer/weight", steering.weight, np.float32)


def _write_sits(builder: BundleBuilder, sits: ProtocolSITS) -> None:
    root = "/sits"
    method = root + "/method"
    _add_scalar(builder, root + "/enabled_default", int(sits.enabled), np.int32)
    _add_scalar(builder, method + "/schema_version", 1, np.int64)
    _add_scalar(builder, method + "/mode", sits.mode)
    if sits.atom_indices:
        _add_array(builder, root + "/atom_indices", sits.atom_indices, np.int32)
    elif isinstance(sits.atom_numbers_policy, (int, np.integer)):
        _add_scalar(builder, root + "/atom_numbers_policy", sits.atom_numbers_policy, np.int64)
    else:
        _add_scalar(builder, root + "/atom_numbers_policy", sits.atom_numbers_policy.upper())

    if sits.k_numbers is not None:
        _add_scalar(builder, method + "/k_numbers", sits.k_numbers, np.int64)
    elif sits.temperature_ladder:
        _add_scalar(builder, method + "/k_numbers", len(sits.temperature_ladder), np.int64)
    if sits.temperature_ladder:
        _add_array(builder, method + "/temperature_ladder", sits.temperature_ladder, np.float32)
    for name, value, dtype in (
        ("temperature_low", sits.temperature_low, np.float32),
        ("temperature_high", sits.temperature_high, np.float32),
        ("pe_a", sits.pe_a, np.float32),
        ("pe_b", sits.pe_b, np.float32),
        ("fb_interval", sits.fb_interval, np.int64),
        ("fb_bias", sits.fb_bias, np.float32),
        ("record_interval", sits.record_interval, np.int64),
        ("update_interval", sits.update_interval, np.int64),
        ("nk_rest", None if sits.nk_rest is None else int(sits.nk_rest), np.int32),
        ("nk_fix", None if sits.nk_fix is None else int(sits.nk_fix), np.int32),
        ("cross_enhance_factor", sits.cross_enhance_factor, np.float32),
    ):
        if value is not None:
            _add_scalar(builder, method + "/" + name, value, dtype)

    restart_root = "/parameters/restart/bias/sits/SITS"
    if sits.initial_nk:
        _add_scalar(
            builder,
            restart_root + "/schema_version",
            1,
            np.int64,
            bundle_file="restart.spgr.h5",
        )
        _add_array(
            builder,
            restart_root + "/nk",
            sits.initial_nk,
            np.float32,
            bundle_file="restart.spgr.h5",
        )
        if sits.initial_log_norm:
            _add_array(
                builder,
                restart_root + "/log_norm",
                sits.initial_log_norm,
                np.float32,
                bundle_file="restart.spgr.h5",
            )
        if sits.initial_log_nk:
            _add_array(
                builder,
                restart_root + "/log_nk",
                sits.initial_log_nk,
                np.float32,
                bundle_file="restart.spgr.h5",
            )


def _write_soft_walls(
    builder: BundleBuilder, soft_walls: tuple[ProtocolSoftWall, ...]
) -> None:
    _add_scalar(builder, "/wall/soft/count", len(soft_walls), np.int64)
    _add_array(
        builder,
        "/wall/soft/name",
        tuple(wall.name for wall in soft_walls),
        object,
    )
    _add_array(
        builder,
        "/wall/soft/potential",
        tuple(wall.potential for wall in soft_walls),
        object,
    )


def _write_hard_wall(builder: BundleBuilder, hard_wall: ProtocolHardWall) -> None:
    _add_array(
        builder,
        "/wall/hard/bounds_low",
        _hard_wall_bounds(hard_wall.bounds_low, low=True),
        np.float32,
    )
    _add_array(
        builder,
        "/wall/hard/bounds_high",
        _hard_wall_bounds(hard_wall.bounds_high, low=False),
        np.float32,
    )
    allow_npt = np.int32(bool(hard_wall.allow_npt))
    _add_scalar(builder, "/wall/hard/allow_npt", allow_npt, np.int32)
    builder.io.set_attrs(
        builder.paths.protocol,
        "/wall/hard",
        {"allow_npt": allow_npt},
    )


def _add_scalar(
    builder: BundleBuilder,
    path: str,
    value: Any,
    dtype=None,
    *,
    bundle_file: str = "protocol.spgp.h5",
) -> None:
    data = value if dtype is None else np.asarray(value, dtype=dtype)
    builder.add_dataset(bundle_file, path, data)


def _add_array(
    builder: BundleBuilder,
    path: str,
    values: Sequence[Any],
    dtype,
    *,
    bundle_file: str = "protocol.spgp.h5",
) -> None:
    builder.add_dataset(bundle_file, path, np.asarray(values, dtype=dtype))


def _validate_name(value: str, label: str) -> None:
    if not isinstance(value, str) or not _OBJECT_NAME.fullmatch(value):
        raise BundleValidationError(
            f"{label} name must match {_OBJECT_NAME.pattern!r}, got {value!r}"
        )


def _require_unique_names(values: Sequence[Any], label: str) -> None:
    names = [value.name for value in values]
    if len(names) != len(set(names)):
        raise BundleValidationError(f"duplicate {label} object name")


def _validate_atom_indices(values: Sequence[int], atom_count: int, label: str) -> None:
    if len(values) != len(set(values)):
        raise BundleValidationError(f"{label} atom indices must be unique")
    if any(not isinstance(value, (int, np.integer)) or value < 0 or value >= atom_count for value in values):
        raise BundleValidationError(
            f"{label} atom indices must be integers in [0, {atom_count})"
        )


def _validate_cv_refs(
    references: Sequence[str],
    cv_by_name: Mapping[str, ProtocolCollectiveVariable],
    label: str,
) -> None:
    missing = [name for name in references if name not in cv_by_name]
    if missing:
        raise BundleValidationError(f"{label} references missing or disabled CVs: {missing}")


def _validate_vector_length(
    values: Sequence[Any],
    count: int,
    label: str,
    *,
    required: bool = False,
) -> None:
    if required and not values:
        raise BundleValidationError(f"{label} is required")
    if values and len(values) != count:
        raise BundleValidationError(f"{label} must contain {count} values")


def _validate_xyz(
    values: Sequence[Sequence[float]],
    count: int,
    label: str,
    *,
    nonnegative: bool = False,
) -> None:
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (count, 3):
        raise BundleValidationError(f"{label} must have shape ({count}, 3)")
    if not np.all(np.isfinite(array)):
        raise BundleValidationError(f"{label} must contain only finite values")
    if nonnegative and np.any(array < 0):
        raise BundleValidationError(f"{label} must contain non-negative values")


def _require_finite(values: Sequence[float], label: str, *, positive: bool = False) -> None:
    if any(not np.isfinite(value) or (positive and value <= 0) for value in values):
        suffix = " finite positive" if positive else " finite"
        raise BundleValidationError(f"{label} must contain only{suffix} values")


def _normalize_parameter(value: Any, label: str):
    if isinstance(value, str):
        if not value or any(token in value for token in ("\n", "\r", "{", "}")):
            raise BundleValidationError(f"{label} contains an invalid string")
        return value
    array = np.asarray(value)
    if array.ndim > 1 or array.dtype.kind not in {"b", "i", "u", "f", "U", "S", "O"}:
        raise BundleValidationError(
            f"{label} must be a scalar or one-dimensional numeric/string value"
        )
    if array.dtype.kind in {"f"} and not np.all(np.isfinite(array)):
        raise BundleValidationError(f"{label} contains a non-finite value")
    if array.dtype.kind in {"U", "S", "O"}:
        values = [str(item) for item in array.reshape(-1)]
        if any(not item or any(token in item for token in ("\n", "\r", "{", "}")) for item in values):
            raise BundleValidationError(f"{label} contains an invalid string")
        if array.ndim == 0:
            return values[0]
        return np.asarray(values, dtype=object)
    return array
