from pathlib import Path
import sys

import PyFiberModes
from PyFiberModes import directories

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib


def test_directory_locations_are_derived_from_module_location():
    package_path = Path(PyFiberModes.__file__).resolve().parent

    assert directories.PACKAGE_PATH == package_path
    assert directories.PROJECT_PATH == package_path.parent
    assert directories.DOCS_PATH == package_path.parent / "docs"
    assert directories.EXAMPLES_PATH == directories.DOCS_PATH / "examples"
    assert directories.DOC_CSS_PATH.is_file()


def test_legacy_directory_names_are_aliases():
    assert directories.root_path is directories.PACKAGE_PATH
    assert directories.project_path is directories.PROJECT_PATH
    assert directories.doc_path is directories.DOCS_PATH
    assert directories.examples_path is directories.EXAMPLES_PATH
    assert directories.doc_css_path is directories.DOC_CSS_PATH


def test_package_has_no_generic_tools_namespace():
    tools_path = directories.PACKAGE_PATH / "tools"
    assert not tools_path.exists() or not any(tools_path.glob("*.py"))


def test_mpsplots_is_not_a_declared_dependency():
    pyproject = tomllib.loads((directories.PROJECT_PATH / "pyproject.toml").read_text())
    conda_recipe = (directories.PROJECT_PATH / "conda.recipe" / "meta.yaml").read_text().lower()

    dependencies = [dependency.lower() for dependency in pyproject["project"]["dependencies"]]

    assert not any("mpsplots" in dependency for dependency in dependencies)
    assert "mpsplots" not in conda_recipe
    assert any(dependency.startswith("matplotlib") for dependency in dependencies)
