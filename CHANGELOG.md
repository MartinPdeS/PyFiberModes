# Changelog

All notable changes to PyFiberModes are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and releases
use [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.14.3] - 2026-09-23

### Fixed

- Build the noarch Conda package once instead of merging nine identically
  named matrix artifacts, preventing corrupted archives during Anaconda
  upload.

## [0.14.2] - 2026-09-23

### Changed

- Rebuilt the documentation around a connected user guide covering units and
  conventions, fiber and material definitions, mode selection, field
  convergence, parameter sweeps, inverse design, coupled-mode propagation,
  validation, and reproducibility.
- Expanded the getting-started and theory pages with installation variants,
  numerical checks, normalized propagation quantities, and guidance near
  cutoff.

## [0.14.1] - 2026-09-23

## [0.14.0] - 2026-09-23

### Added

- Expanded the specialty-fiber, field, analysis, validation, and performance
  galleries so every documentation series contains multiple examples.
- Defined a stable root-level public API, compatibility and deprecation policy,
  API snapshot regression test, and migration guide.

### Changed

- Replaced the opaque `ssif`, `tlsif`, `mlsif`, and `NeffSolver` names with
  descriptive two-layer, three-layer, multilayer, and effective-index solver
  namespaces and class names.

## [0.13.2] - 2026-09-23

### Changed

- Raised the minimum supported Python version to 3.11 and aligned CI,
  contributor guidance, Ruff, and Conda metadata with that policy.
- Moved Matplotlib rendering into the optional `PyFiberModes.plotting` module
  and added the `plotting` installation extra, keeping it out of core runtime
  dependencies while preserving the `Field.plot` convenience method.
- Replaced generic package descriptions and keywords with fiber-mode,
  photonics, waveguide, dispersion, and coupled-mode metadata.

## [0.13.1] - 2026-09-22

### Added

- Added validated immutable layer and solver settings, structural solver
  protocols, domain-specific exceptions, and structured solver diagnostics.
- Added deterministic physical-reference, convergence, caching, architecture,
  and opt-in `pytest-benchmark` performance suites.

### Changed

- Unified cylindrical coordinates under one canonical model, split modal and
  field orchestration into focused services, and added geometry-aware caching.
- Vectorized radial index lookup and routed field construction through a
  single batched radial-field interface.
- Standardized the group-velocity API on `get_group_velocity` and removed
  unreferenced experimental modules, unfinished solver branches, misspelled
  methods, lowercase directory aliases, and obsolete compatibility wrappers.

## [0.13.0] - 2026-09-22

### Changed

- Replaced MPSPlots styling and colormaps with native Matplotlib equivalents,
  removing MPSPlots from all runtime and documentation dependencies.
- Field maps now share one radial solver pass across all electric and magnetic
  components and cache interpolated, azimuthal, and batched component arrays.
- Added NumPyDoc documentation across modules, public APIs, data models, and
  numerical solver helpers, with structural coverage enforced by tests.
- Reorganized the documentation into learning, theory, workflow, performance,
  troubleshooting, gallery, and task-oriented API sections; added four new
  example series and a unified logo and favicon identity.

## [0.12.0] - 2026-09-22

### Changed

- Consolidated filesystem locations in the root `directories` module and removed
  the duplicate internal `tools` package.
- Moved solver-test helpers out of the distributed package and into the test suite.

## [0.11.0] - 2026-09-22

### Added

- Automatic guided-mode discovery and structured parameter sweeps.
- Constrained inverse design for layer radii, refractive indices, and custom parameters.
- Mode-overlap, coupling-matrix, and coupled-mode propagation utilities.
- Vectorized field mapping, Poynting vectors, guided power, and confinement factors.
- Citation metadata for Citation File Format and Zenodo.
- Release-metadata consistency checks and expanded behavioral unit tests.

### Changed

- Field maps reuse cached radial solutions rather than solving independently at every pixel.
- Effective-area integration now includes the physical grid-cell area consistently.

### Fixed

- Field intensity and normalization now use the current `Fiber.get_effective_index` API.
- Radial-field caches are cleared when the fiber wavelength changes.

## [0.10.0] - 2026-02-13

### Changed

- Established the current tagged baseline for subsequent documented releases.

[Unreleased]: https://github.com/MartinPdeS/PyFiberModes/compare/v0.14.3...HEAD
[0.14.3]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.14.3
[0.14.2]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.14.2
[0.14.1]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.14.1
[0.14.0]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.14.0
[0.13.2]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.13.2
[0.13.1]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.13.1
[0.13.0]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.13.0
[0.12.0]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.12.0
[0.11.0]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.11.0
[0.10.0]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.10.0
