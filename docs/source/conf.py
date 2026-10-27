"""PyFiberModes documentation, using the shared MartinPdeS-docs theme."""

from pathlib import Path
import sys

from martinpdes_docs import configure


ROOT = Path(__file__).resolve().parents[2]
# Gallery examples and autodoc must describe this checkout.
sys.path.insert(0, str(ROOT))

gallery_sections = [
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
        version="auto",
        repository="MartinPdeS/PyFiberModes",
        logo="_static/logo.png",
        favicon="_static/favicon.png",
        gallery_sections=gallery_sections,
        gallery_default_thumb="_static/thumbnail.png",
        gallery_options={"abort_on_example_error": True},
        logo_text="",
        autodoc_default_options={"members": True},
    )
)
