"""Unit tests for explicit base-component charge protocols."""

from dataclasses import replace
import unittest

from Xponge.metal_assignment import (
    BaseChargeInput,
    ValidationError,
    base_charge_input_dumps,
    base_charge_input_loads,
    validate_base_charge_input,
)


class BaseChargeInputTests(unittest.TestCase):
    def test_hash_closed_round_trip(self):
        value = BaseChargeInput(
            schema_version=1,
            method="tpacm4",
            source="unit-test",
        ).with_computed_hash()

        restored = base_charge_input_loads(base_charge_input_dumps(value))

        self.assertEqual(restored, value)
        validate_base_charge_input(restored)

    def test_stale_hash_and_implicit_method_are_rejected(self):
        value = BaseChargeInput(
            schema_version=1,
            method="tpacm4",
            source="unit-test",
        ).with_computed_hash()
        with self.assertRaisesRegex(ValidationError, "base charge input hash mismatch"):
            validate_base_charge_input(replace(value, source="changed"))
        with self.assertRaisesRegex(ValidationError, "unsupported_base_charge_method"):
            validate_base_charge_input(replace(value, method="auto", charge_input_hash=""))


if __name__ == "__main__":
    unittest.main()
