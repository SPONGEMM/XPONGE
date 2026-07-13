"""
SPONGE legacy-to-bundle conversion helpers.
"""

from .converter import ConversionError, LegacyToBundleConverter, convert_legacy_to_bundle
from .errors import (
    BundleCapabilityError,
    BundleConflictError,
    BundleError,
    BundleExportError,
    BundlePathError,
    BundleSchemaError,
    BundleValidationError,
)
from .output_writer import LegacyOutputBundleWriter, LegacyOutputConversionError, convert_legacy_outputs_to_bundle
from .reverse_converter import BundleToLegacyConverter, convert_bundle_to_legacy
from .saver import save_sponge_input_bundle

__all__ = [
    "ConversionError",
    "BundleError",
    "BundleSchemaError",
    "BundleValidationError",
    "BundleExportError",
    "BundleCapabilityError",
    "BundlePathError",
    "BundleConflictError",
    "LegacyToBundleConverter",
    "BundleToLegacyConverter",
    "LegacyOutputBundleWriter",
    "LegacyOutputConversionError",
    "convert_legacy_to_bundle",
    "convert_bundle_to_legacy",
    "convert_legacy_outputs_to_bundle",
    "save_sponge_input_bundle",
]
