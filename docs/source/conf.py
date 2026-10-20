"""PyFiberModes documentation, using the shared MartinPdeS-docs theme."""

import os
from pathlib import Path
import sys
import tomllib

from martinpdes_docs import configure


ROOT = Path(__file__).resolve().parents[2]
# Gallery examples and autodoc must describe this checkout.
sys.path.insert(0, str(ROOT))

examples_sections = [
    "basic",
    "SMF28",
    "specialty_fibers",
    "fields",
    "analysis",
    "validation",
    "benchmarks",
]

globals().update(
    configure(
        project="PyFiberModes",
        package="PyFiberModes",
        repository="MartinPdeS/PyFiberModes",
        logo="_static/logo.svg",
        favicon="_static/favicon.svg",
        version=os.getenv("tag", "latest"),
        html_theme_options={
            "logo": {"text": "", "image": "_static/logo.svg"},
        },
        sphinx_gallery_conf={
            "examples_dirs": [f"../examples/{name}" for name in examples_sections],
            "gallery_dirs": [f"gallery/{name}" for name in examples_sections],
            "default_thumb_file": str(
                Path(__file__).resolve().parent / "_static" / "favicon.png"
            ),
            "within_subsection_order": "FileNameSortKey",
            "abort_on_example_error": True,
        },
    )
)

extensions += ["sphinx.ext.mathjax", "sphinx.ext.autosectionlabel"]
autodoc_default_options["members"] = True
html_context["doc_path"] = "docs/source"
html_context["default_mode"] = "dark"

# Keep page substitutions independent of importing the numerical package.
metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
version = release = metadata["project"]["version"]
