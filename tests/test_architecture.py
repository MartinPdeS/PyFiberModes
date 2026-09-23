from pathlib import Path
import subprocess
import sys

import PyFiberModes
from PyFiberModes import directories
from PyFiberModes import solver

import tomllib


def test_directory_locations_are_derived_from_module_location():
    package_path = Path(PyFiberModes.__file__).resolve().parent

    assert directories.PACKAGE_PATH == package_path
    assert directories.PROJECT_PATH == package_path.parent
    assert directories.DOCS_PATH == package_path.parent / "docs"
    assert directories.EXAMPLES_PATH == directories.DOCS_PATH / "examples"
    assert directories.DOC_CSS_PATH.is_file()


def test_package_has_no_generic_tools_namespace():
    tools_path = directories.PACKAGE_PATH / "tools"
    assert not tools_path.exists() or not any(tools_path.glob("*.py"))


def test_mpsplots_is_not_a_declared_dependency():
    pyproject = tomllib.loads((directories.PROJECT_PATH / "pyproject.toml").read_text())
    conda_recipe = (directories.PROJECT_PATH / "conda.recipe" / "meta.yaml").read_text().lower()

    dependencies = [dependency.lower() for dependency in pyproject["project"]["dependencies"]]

    assert not any("mpsplots" in dependency for dependency in dependencies)
    assert "mpsplots" not in conda_recipe
    assert not any(dependency.startswith("matplotlib") for dependency in dependencies)
    assert any(
        dependency.startswith("matplotlib")
        for dependency in pyproject["project"]["optional-dependencies"]["plotting"]
    )


def test_importing_core_does_not_import_plotting_module():
    """Keep the package plotting module outside the core import graph."""
    script = """
import sys
import types

pyoptik = types.ModuleType("PyOptik")
pyoptik.MaterialBank = object
sys.modules["PyOptik"] = pyoptik

import PyFiberModes

assert "PyFiberModes.plotting" not in sys.modules
"""
    subprocess.run(
        [sys.executable, "-c", script],
        cwd=directories.PROJECT_PATH,
        check=True,
    )


def test_solver_names_describe_supported_geometries():
    """Expose geometry names instead of historical initialisms."""
    assert hasattr(solver, "two_layer")
    assert hasattr(solver, "three_layer")
    assert hasattr(solver, "multilayer")
    assert not any(hasattr(solver, name) for name in ("ssif", "tlsif", "mlsif"))


def test_each_documentation_gallery_has_multiple_examples():
    """Prevent example galleries from regressing to a single page."""
    galleries = (
        "basic",
        "SMF28",
        "specialty_fibers",
        "fields",
        "analysis",
        "validation",
        "benchmarks",
    )
    for gallery in galleries:
        examples = tuple((directories.EXAMPLES_PATH / gallery).glob("*.py"))
        assert len(examples) >= 2, f"{gallery} needs at least two examples"
