# Changelog

All notable changes to PyFiberModes are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and releases
use [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/MartinPdeS/PyFiberModes/compare/v0.11.0...HEAD
[0.11.0]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.11.0
[0.10.0]: https://github.com/MartinPdeS/PyFiberModes/releases/tag/v0.10.0
