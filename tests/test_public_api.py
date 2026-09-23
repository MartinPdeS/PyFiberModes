"""Regression tests for the supported package-level API contract."""

import json
import inspect
from pathlib import Path

import PyFiberModes


SNAPSHOT = Path(__file__).with_name("api_snapshot.json")
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_public_api_matches_reviewed_snapshot():
    """Require explicit review for every supported import-surface change."""
    expected = json.loads(SNAPSHOT.read_text(encoding="utf-8"))

    assert expected["policy_version"] == 2
    assert expected["public_import"] == "PyFiberModes"
    assert sorted(PyFiberModes.__all__) == expected["exports"]


def _resolve_public_name(name):
    """Resolve a dotted name from the supported package namespace."""
    value = PyFiberModes
    for component in name.split("."):
        value = getattr(value, component)
    return value


def _stable_default(value):
    """Represent reviewed defaults without implementation-specific reprs."""
    if value is inspect.Parameter.empty:
        raise AssertionError("required parameters do not have defaults")
    if isinstance(value, (str, int, float, bool, type(None))):
        return value
    if isinstance(value, tuple):
        return list(value)
    if type(value).__name__ == "_HAS_DEFAULT_FACTORY_CLASS":
        return "<factory>"
    return f"<{type(value).__name__}>"


def test_public_signatures_and_defaults_match_reviewed_snapshot():
    """Freeze parameter order and important defaults for the 1.0 contract."""
    expected = json.loads(SNAPSHOT.read_text(encoding="utf-8"))
    for name, contract in expected["signatures"].items():
        signature = inspect.signature(_resolve_public_name(name))
        assert list(signature.parameters) == contract["parameters"], name
        defaults = {
            parameter.name: _stable_default(parameter.default)
            for parameter in signature.parameters.values()
            if parameter.default is not inspect.Parameter.empty
        }
        assert defaults == contract["defaults"], name


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
    assert "public_api_contract.rst" in navigation
    assert "migration_guide.rst" in navigation
    assert "at least two minor releases and six months" in policy
    assert "Unreleased" in migration
