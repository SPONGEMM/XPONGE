"""Shared errors for SPONGE bundled input conversion."""

from __future__ import annotations


class BundleError(RuntimeError):
    """Base class for bundled input failures."""


class BundleSchemaError(BundleError):
    """Raised when an HDF5 artifact has an unsupported schema."""


class BundleValidationError(BundleError):
    """Raised when bundled artifacts are internally inconsistent."""


class BundleExportError(BundleError):
    """Raised when typed bundled data cannot be exported."""


class BundleCapabilityError(BundleError):
    """Raised when a known input has no implemented bundled codec."""


class BundlePathError(BundleError):
    """Raised when an input or sidecar path escapes its allowed root."""


class BundleConflictError(BundleError):
    """Raised when reverse conversion would overwrite an existing target."""


class BundleMDAnalysisError(BundleValidationError):
    """Base class for bundle-to-MDAnalysis integration failures."""


class BundleTopologyError(BundleMDAnalysisError):
    """Raised when a bundled topology cannot be exposed to MDAnalysis."""


class BundleTrajectoryError(BundleMDAnalysisError):
    """Raised when a SPONGE H5MD trajectory is invalid or unsupported."""


class BundleUnitError(BundleMDAnalysisError):
    """Raised when bundle units cannot be converted to MDAnalysis units."""


class UnverifiedBundlePairError(BundleMDAnalysisError):
    """Raised when topology/trajectory compatibility cannot be proven."""


class IncompleteBundleError(BundleMDAnalysisError):
    """Raised when a trajectory bundle is not finalized or is truncated."""


class AmbiguousH5MDLayoutError(BundleTrajectoryError):
    """Raised when a H5MD file matches conflicting SPONGE layouts."""
