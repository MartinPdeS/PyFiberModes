# Contributing to PyFiberModes

Use Python 3.10 or newer. Install development dependencies with `python -m pip install -e ".[testing,documentation]"`, then run `python -m pytest`.

Releases use semantic version tags. Run `make patch`, `make minor`, or `make major` from a clean `master` checkout; the command creates the release commit, tags it, and pushes both. Publishing workflows run only for `v*` tags.

The Conda recipe is at `conda.recipe/meta.yaml`. Do not commit generated documentation, virtual environments, build artefacts, or caches.
