Architecture
============

PyFiberModes separates physical state, orchestration, and numerical backends
while retaining ``Fiber`` as the convenient public façade.

Physical models
---------------

``Fiber`` owns wavelength and ordered layers. ``LayerSpec`` and
``SolverSettings`` provide immutable validated inputs for applications that
need stricter construction boundaries. Cartesian and cylindrical coordinates
have one canonical implementation in ``PyFiberModes.coordinates``.

Services
--------

``fiber.analysis`` owns geometry-aware modal and cutoff caches and exposes
structured ``SolverResult`` diagnostics. ``fiber.fields`` constructs sampled
field grids. Existing methods such as ``fiber.get_effective_index`` delegate
to these services, keeping established code compatible.

Solver contracts
----------------

Runtime-checkable protocols describe effective-index, cutoff, and radial-field
backends. Third-party implementations can satisfy these interfaces without
inheriting package internals. Package exceptions distinguish invalid physical
input from solver and convergence failures.

Concrete solvers use descriptive namespaces: ``solver.two_layer`` provides
effective-index, cutoff, and radial-field calculations;
``solver.three_layer`` provides three-layer cutoff calculations; and
``solver.multilayer`` provides the general effective-index calculation.

Performance
-----------

Radial index profiles use vectorized boundary lookup. Field components share a
single batched radial request, interpolation maps, and azimuthal factors.
Install the ``benchmarking`` extra and run ``pytest benchmarks
--benchmark-only`` to compare performance between revisions.
