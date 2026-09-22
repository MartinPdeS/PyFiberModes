"""Canonical filesystem locations used by PyFiberModes and its documentation."""

from pathlib import Path


__all__ = [
    "PACKAGE_PATH",
    "PROJECT_PATH",
    "DOCS_PATH",
    "EXAMPLES_PATH",
    "DOC_CSS_PATH",
    "root_path",
    "project_path",
    "examples_path",
    "doc_path",
    "doc_css_path",
]

PACKAGE_PATH = Path(__file__).resolve().parent
PROJECT_PATH = PACKAGE_PATH.parent
DOCS_PATH = PROJECT_PATH / "docs"
EXAMPLES_PATH = DOCS_PATH / "examples"
DOC_CSS_PATH = DOCS_PATH / "source" / "_static" / "default.css"

# Backwards-compatible aliases for the original public names.
root_path = PACKAGE_PATH
project_path = PROJECT_PATH
doc_path = DOCS_PATH
examples_path = EXAMPLES_PATH
doc_css_path = DOC_CSS_PATH
