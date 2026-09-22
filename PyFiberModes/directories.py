"""Canonical filesystem locations used by PyFiberModes and its documentation."""

from pathlib import Path


__all__ = [
    "PACKAGE_PATH",
    "PROJECT_PATH",
    "DOCS_PATH",
    "EXAMPLES_PATH",
    "DOC_CSS_PATH",
]

PACKAGE_PATH = Path(__file__).resolve().parent
PROJECT_PATH = PACKAGE_PATH.parent
DOCS_PATH = PROJECT_PATH / "docs"
EXAMPLES_PATH = DOCS_PATH / "examples"
DOC_CSS_PATH = DOCS_PATH / "source" / "_static" / "default.css"
