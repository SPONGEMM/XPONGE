"""QM subsystem exceptions."""

from __future__ import annotations


class QMError(RuntimeError):
    """Base error for the shared QM subsystem."""


class QMBackendImportError(ImportError, QMError):
    """Raised when an optional QM backend dependency is unavailable."""


class QMBackendSelectionError(ValueError, QMError):
    """Raised when a requested backend is unknown or invalid."""


class QMCapabilityError(NotImplementedError, QMError):
    """Raised when a backend lacks a requested capability."""
