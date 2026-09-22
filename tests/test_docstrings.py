"""Structural documentation checks for the Python package."""

import ast
from pathlib import Path


PACKAGE_PATH = Path(__file__).resolve().parents[1] / "PyFiberModes"


def _documentable_nodes(tree):
    """Yield every class and function definition from an AST."""
    return (
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
    )


def test_all_modules_classes_and_functions_have_docstrings():
    """Require documentation for every maintained Python callable."""
    missing = []
    for path in sorted(PACKAGE_PATH.rglob("*.py")):
        if path.name == "_version.py":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        if ast.get_docstring(tree) is None:
            missing.append(f"{path.relative_to(PACKAGE_PATH)}: module")
        for node in _documentable_nodes(tree):
            if ast.get_docstring(node) is None:
                missing.append(
                    f"{path.relative_to(PACKAGE_PATH)}:{node.lineno}: {node.name}"
                )

    assert not missing, "Missing docstrings:\n" + "\n".join(missing)
