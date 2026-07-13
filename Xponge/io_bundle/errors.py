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
