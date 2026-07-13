"""Static contract checks for bidirectional SPONGE bundled input support."""

from dataclasses import replace

import pytest

from Xponge.io_bundle.contracts import (
    CONTRACTS,
    LEGACY_KEY_ALIASES,
    REVERSE_POLICIES,
    SPONGE_SERIALIZER_METADATA_KEYS,
    SPONGE_SERIALIZER_TO_LEGACY_KEY,
    SPONGE_SERIALIZER_UNSUPPORTED_KEYS,
    classify_sponge_serializer_key,
    contracts_by_bundle_file,
    contracts_by_legacy_key,
    reversible_contracts,
    validate_contract_registry,
)
from Xponge.io_bundle.exporters import EXPORTERS


def test_contract_registry_is_valid():
    validate_contract_registry()


def test_contract_indexes_contain_every_contract():
    by_legacy_key = contracts_by_legacy_key()
    by_bundle_file = contracts_by_bundle_file()

    for contract in CONTRACTS:
        assert contract in by_bundle_file[contract.bundle_file]
        for key in contract.legacy_keys:
            assert contract in by_legacy_key[key]


def test_reversible_contracts_are_canonical_and_materializable():
    for mode in ("normal", "rerun"):
        selected = reversible_contracts(mode)
        assert selected
        for contract in selected:
            assert contract.direction == "input"
            assert contract.reverse_policy in REVERSE_POLICIES - {"alias", "not_reversible"}
            if contract.reverse_policy in {"typed_or_sidecar", "typed_required"}:
                assert contract.exporter_id
                assert contract.required_bundle_paths
                assert contract.legacy_filename_stem


def test_alias_contracts_point_to_canonical_legacy_keys():
    by_legacy_key = contracts_by_legacy_key()
    for alias, canonical in LEGACY_KEY_ALIASES.items():
        assert alias in by_legacy_key
        assert canonical in by_legacy_key
        assert all(contract.reverse_policy == "alias" for contract in by_legacy_key[alias])
        assert any(alias in contract.aliases for contract in by_legacy_key[canonical])


def test_save_sponge_serializer_registry_is_fully_classified():
    by_legacy_key = contracts_by_legacy_key()

    for serializer_key, legacy_key in SPONGE_SERIALIZER_TO_LEGACY_KEY.items():
        assert classify_sponge_serializer_key(serializer_key) == ("contract", legacy_key)
        assert legacy_key in by_legacy_key
    for serializer_key in SPONGE_SERIALIZER_METADATA_KEYS:
        assert classify_sponge_serializer_key(serializer_key) == ("metadata", None)
    for serializer_key in SPONGE_SERIALIZER_UNSUPPORTED_KEYS:
        assert classify_sponge_serializer_key(serializer_key) == ("unsupported", None)
    assert classify_sponge_serializer_key("future_force") == ("unknown", None)


def test_all_canonical_typed_contracts_have_exporters():
    for contract in CONTRACTS:
        if contract.reverse_policy not in {"typed_or_sidecar", "typed_required"}:
            continue
        assert contract.exporter_id in EXPORTERS, contract.contract_id


def test_improper_is_typed_required_with_alias_only_for_input():
    by_legacy_key = contracts_by_legacy_key()
    canonical = by_legacy_key["improper_dihedral_in_file"]
    alias = by_legacy_key["improper_in_file"]

    assert len(canonical) == 1
    assert canonical[0].reverse_policy == "typed_required"
    assert canonical[0].required_bundle_paths == ("/forcefield/improper",)
    assert len(alias) == 1
    assert alias[0].reverse_policy == "alias"


def test_contract_validation_rejects_duplicate_ids():
    duplicate = replace(CONTRACTS[0], bundle_path="/duplicate")
    with pytest.raises(ValueError, match="duplicate I/O contract id"):
        validate_contract_registry(CONTRACTS + (duplicate,))


def test_contract_validation_rejects_unknown_reverse_policy():
    invalid = replace(CONTRACTS[0], contract_id="invalid.policy", reverse_policy="guess")
    with pytest.raises(ValueError, match="unknown reverse policy"):
        validate_contract_registry(CONTRACTS + (invalid,))


def test_contract_validation_rejects_missing_exporter_metadata():
    invalid = replace(
        CONTRACTS[0],
        contract_id="invalid.exporter",
        exporter_id=None,
        reverse_policy="typed_required",
    )
    with pytest.raises(ValueError, match="no reverse exporter id"):
        validate_contract_registry(CONTRACTS + (invalid,))
