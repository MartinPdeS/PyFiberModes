"""Regression tests for the supported package-level API contract."""

import json
from pathlib import Path

import PyFiberModes


SNAPSHOT = Path(__file__).with_name("api_snapshot.json")
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_public_api_matches_reviewed_snapshot():
    """Require explicit review for every supported import-surface change."""
    expected = json.loads(SNAPSHOT.read_text(encoding="utf-8"))

    assert expected["policy_version"] == 1
    assert expected["public_import"] == "PyFiberModes"
    assert sorted(PyFiberModes.__all__) == expected["exports"]


def test_every_public_export_is_importable():
    """Ensure every promised root export resolves to a concrete object."""
    for name in PyFiberModes.__all__:
        assert getattr(PyFiberModes, name) is not None


def test_version_is_metadata_not_part_of_symbol_snapshot():
    """Expose package metadata without treating its value as a stable symbol."""
    assert isinstance(PyFiberModes.__version__, str)
    assert "__version__" not in PyFiberModes.__all__


def test_compatibility_documents_are_published_in_navigation():
    """Keep the policy and migration guide visible in built documentation."""
    documentation = PROJECT_ROOT / "docs" / "source"
    navigation = (documentation / "index.rst").read_text(encoding="utf-8")
    policy = (documentation / "api_policy.rst").read_text(encoding="utf-8")
    migration = (documentation / "migration_guide.rst").read_text(encoding="utf-8")

    assert "api_policy.rst" in navigation
    assert "migration_guide.rst" in navigation
    assert "at least two minor releases and six months" in policy
    assert "Unreleased" in migration
