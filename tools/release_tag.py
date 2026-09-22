#!/usr/bin/env python3
"""Synchronize release metadata, commit it, and create an annotated Git tag."""

import argparse
from datetime import date
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
TAG_PATTERN = re.compile(r"v(?P<version>(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)\.(?:0|[1-9]\d*))$")


def run(*command: str, capture_output: bool = False, env=None) -> str:
    completed = subprocess.run(
        command, cwd=ROOT, check=True, text=True, capture_output=capture_output, env=env
    )
    return completed.stdout.strip() if capture_output else ""


def validate_tag(tag: str) -> str:
    match = TAG_PATTERN.fullmatch(tag)
    if match is None:
        raise ValueError("tag must use the form vMAJOR.MINOR.PATCH, for example v0.11.0")
    return match.group("version")


def require_clean_worktree() -> None:
    if run("git", "status", "--porcelain", capture_output=True):
        raise RuntimeError("working tree is not clean; commit or stash changes first")


def require_unused_tag(tag: str) -> None:
    if run("git", "tag", "--list", tag, capture_output=True):
        raise RuntimeError(f"tag {tag} already exists")


def replace_once(path: Path, pattern: str, replacement: str, description: str) -> None:
    text = path.read_text(encoding="utf-8")
    updated, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"could not update {description} in {path.relative_to(ROOT)}")
    path.write_text(updated, encoding="utf-8")


def update_metadata(version: str) -> None:
    zenodo_path = ROOT / ".zenodo.json"
    zenodo = json.loads(zenodo_path.read_text(encoding="utf-8"))
    zenodo.update(version=version, publication_date=date.today().isoformat())
    zenodo_path.write_text(f"{json.dumps(zenodo, indent=2)}\n", encoding="utf-8")
    replace_once(
        ROOT / "pyproject.toml", r'^(version\s*=\s*)["\'][^"\']+["\']',
        rf'\g<1>"{version}"', "project version",
    )
    replace_once(
        ROOT / "conda.recipe" / "meta.yaml", r"^(\s*version:\s*)\S+", rf"\g<1>{version}",
        "Conda version",
    )


def update_changelog(version: str) -> None:
    path = ROOT / "CHANGELOG.md"
    text = path.read_text(encoding="utf-8")
    marker = "## [Unreleased]\n"
    if marker not in text:
        raise RuntimeError("CHANGELOG.md has no [Unreleased] section")
    text = text.replace(marker, f"{marker}\n## [{version}] - {date.today().isoformat()}\n", 1)
    text, count = re.subn(
        r"^\[Unreleased\]:.*$",
        f"[Unreleased]: https://github.com/MartinPdeS/PyFiberModes/compare/v{version}...HEAD\n"
        f"[{version}]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v{version}",
        text, count=1, flags=re.MULTILINE,
    )
    if count != 1:
        raise RuntimeError("CHANGELOG.md has no [Unreleased] comparison link")
    path.write_text(text, encoding="utf-8")


def generate_version_file(version: str) -> None:
    environment = os.environ.copy()
    environment["SETUPTOOLS_SCM_PRETEND_VERSION"] = version
    run(sys.executable, "-m", "vcs_versioning", "--force-write-version-files", env=environment)


def create_release(tag: str, version: str) -> None:
    update_metadata(version)
    update_changelog(version)
    generate_version_file(version)
    run(
        "git", "add", ".zenodo.json", "CHANGELOG.md", "conda.recipe/meta.yaml",
        "pyproject.toml", "PyFiberModes/_version.py",
    )
    run("git", "commit", "-m", f"Release {tag}")
    run("git", "tag", "-a", tag, "-m", f"Release {tag}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag", help="annotated release tag, for example v0.11.0")
    arguments = parser.parse_args()
    try:
        version = validate_tag(arguments.tag)
        require_clean_worktree()
        require_unused_tag(arguments.tag)
        create_release(arguments.tag, version)
    except (RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"release aborted: {error}", file=sys.stderr)
        return 1
    print(f"created release commit and annotated tag {arguments.tag}")
    print("Push it when ready with: git push origin HEAD --tags")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
