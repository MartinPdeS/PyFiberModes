import importlib.util
import json
from pathlib import Path

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]


def load_tool(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_release_metadata_is_synchronized_and_changelog_has_current_version():
    checker = load_tool("check_release")
    assert set(checker.versions().values()) == {"0.10.0"}
    assert "0.10.0" in checker.changelog_versions()


def test_citation_and_zenodo_metadata_agree():
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    assert citation["title"] == "PyFiberModes"
    assert citation["license"] == zenodo["license"] == "MIT"
    assert citation["repository-code"] == "https://github.com/MartinPdeS/PyFiberModes"
    assert citation["authors"][0]["family-names"] in zenodo["creators"][0]["name"]


@pytest.mark.parametrize(
    ("kind", "expected"),
    [("major", (2, 0, 0)), ("minor", (1, 3, 0)), ("patch", (1, 2, 4))],
)
def test_semantic_version_increment(kind, expected):
    release = load_tool("next_release_version")
    assert release.next_version((1, 2, 3), kind) == expected


def test_release_tag_validation_is_strict():
    release = load_tool("release_tag")
    assert release.validate_tag("v1.2.3") == "1.2.3"
    for invalid in ("1.2.3", "v1.2", "v01.2.3", "v1.2.3rc1"):
        with pytest.raises(ValueError):
            release.validate_tag(invalid)


def test_release_changelog_rolls_unreleased_entries_forward(tmp_path, monkeypatch):
    release = load_tool("release_tag")
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        "# Changelog\n\n## [Unreleased]\n\n### Added\n\n- Feature.\n\n"
        "[Unreleased]: https://example.test/compare/v1.0.0...HEAD\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(release, "ROOT", tmp_path)
    release.update_changelog("1.1.0")
    updated = changelog.read_text(encoding="utf-8")
    assert "## [Unreleased]\n\n## [1.1.0]" in updated
    assert "compare/v1.1.0...HEAD" in updated
    assert "[1.1.0]:" in updated
