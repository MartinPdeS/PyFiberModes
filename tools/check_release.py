#!/usr/bin/env python3
"""Check that PyFiberModes release metadata agrees across project files."""

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]


def normalized_version(version: str) -> str:
    return version.removeprefix("v")


def extract(pattern: str, path: Path, description: str) -> str:
    match = re.search(pattern, path.read_text(encoding="utf-8"), flags=re.MULTILINE)
    if match is None:
        raise RuntimeError(f"could not find {description} in {path.relative_to(ROOT)}")
    return match.group(1)


def versions() -> dict[str, str]:
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    zenodo = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    return {
        "pyproject.toml": str(project["project"]["version"]),
        ".zenodo.json": str(zenodo["version"]),
        "conda.recipe/meta.yaml": extract(
            r"^\s*version:\s*['\"]?([^'\"\s]+)",
            ROOT / "conda.recipe" / "meta.yaml",
            "Conda version",
        ),
        "PyFiberModes/_version.py": extract(
            r"^__version__\s*=\s*version\s*=\s*['\"]([^'\"]+)",
            ROOT / "PyFiberModes" / "_version.py", "generated package version",
        ),
    }


def changelog_versions() -> set[str]:
    text = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    return set(re.findall(r"^## \[(\d+\.\d+\.\d+)\]", text, flags=re.MULTILINE))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", help="expected version, optionally prefixed with v")
    arguments = parser.parse_args()
    try:
        declared = versions()
        releases = changelog_versions()
    except (KeyError, OSError, RuntimeError, json.JSONDecodeError) as error:
        print(f"release check failed: {error}", file=sys.stderr)
        return 1
    expected = (
        normalized_version(arguments.version)
        if arguments.version
        else declared["pyproject.toml"]
    )
    failures = []
    for source, version in declared.items():
        status = "OK" if version == expected else "MISMATCH"
        print(f"{status:8} {source}: {version}")
        if version != expected:
            failures.append(f"{source} declares {version}; expected {expected}")
    if expected not in releases:
        failures.append(f"CHANGELOG.md has no release section for {expected}")
    tag = f"v{expected}"
    tag_exists = subprocess.run(
        ["git", "rev-parse", "--verify", "--quiet", f"refs/tags/{tag}"], cwd=ROOT,
        check=False, stdout=subprocess.DEVNULL,
    ).returncode == 0
    tag_status = "exists" if tag_exists else "not created"
    print(f"{'OK' if tag_exists else 'INFO':8} Git tag: {tag} {tag_status}")
    if failures:
        for failure in failures:
            print(f"ERROR: {failure}", file=sys.stderr)
        return 1
    print(f"Release metadata consistently declares {expected}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
