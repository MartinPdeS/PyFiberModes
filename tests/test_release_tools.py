import json
from pathlib import Path

from martinpdes_maintenance import cli
import yaml


ROOT = Path(__file__).resolve().parents[1]


def test_release_metadata_is_synchronized_and_changelog_has_current_version():
    assert cli.main(["--root", str(ROOT), "check"]) == 0


def test_citation_and_zenodo_metadata_agree():
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    assert citation["title"] == "PyFiberModes"
    assert citation["license"] == zenodo["license"] == "MIT"
    assert citation["repository-code"] == "https://github.com/MartinPdeS/PyFiberModes"
    assert citation["authors"][0]["family-names"] in zenodo["creators"][0]["name"]


def test_release_changelog_rolls_unreleased_entries_forward():
    text = (
        "# Changelog\n\n## [Unreleased]\n\n### Added\n\n- Feature.\n\n"
        "[Unreleased]: https://example.test/compare/v1.0.0...HEAD\n"
    )
    updated = cli.update_changelog(text, "1.1.0", "2026-10-25", "MartinPdeS/PyFiberModes")
    assert "## [Unreleased]\n\n## [1.1.0]" in updated
    assert "compare/v1.1.0...HEAD" in updated
    assert "[1.1.0]:" in updated
