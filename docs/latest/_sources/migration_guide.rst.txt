.. _migration_guide:

Migration guide
===============

This page records user-visible name and behavior changes. Changes to internal
solver modules are included when they are likely to affect advanced users,
even though those modules are outside the stable public API.

Unreleased
----------

Explicit physical units
~~~~~~~~~~~~~~~~~~~~~~~

Physical inputs now require Pint quantities from ``PyFiberModes.units.ureg``,
which is the shared ``TypedUnit.ureg`` registry also used by PyMieSim. Bare
wavelengths, radii, field angles and propagation constants raise ``TypeError``;
incompatible dimensions raise ``pint.DimensionalityError``. Indices require
explicit dimensionless quantities. Wavelength is now required by ``load_fiber``.

For example, replace ``wavelength=1550e-9`` with
``wavelength=1550 * ureg.nanometer`` and ``radius=4.1e-6`` with
``radius=4.1 * ureg.micrometer``. Use ``index=1.45 * ureg.dimensionless``.
An unbounded outer layer uses ``radius=float("inf") * ureg.meter``.

Physical outputs are Pint quantities, including dimensionless modal indices.
Use ``result.to("desired_unit").magnitude`` when a plain number or NumPy array
is needed. Sweep arrays and records retain units. Optimizer objectives return
plain scalar costs; explicitly convert dimensionful observables in the objective.

Default modal output units are meters for cutoff wavelength, inverse meters
for propagation constants, meters per second for velocities, seconds squared
per meter for GVD, ps/(nm km) for dispersion, and ps/(nm² km) for the historical
S-parameter. Fields use V/m and A/m, coordinates use meters, effective areas
use square meters, and integrated powers use watts. Analytical fields retain
an arbitrary overall amplitude; absolute power requires choosing normalization.

Coupled-mode inputs ``betas`` and ``coupling`` require inverse-length units.
A position-dependent coupling callback receives a length quantity and returns
an inverse-length matrix. Normalized field overlaps are dimensionless;
``coupling_matrix`` inherits the units of its perturbation weighting.

The numerical solvers and YAML geometry remain in SI internally. Private
``_fiber``, ``_field`` and ``_coordinates`` modules are implementation details.
The public ``Fiber`` constructor now takes only wavelength; add geometry with
``add_layer`` rather than supplying the former dataclass storage lists.

Mode validation and failure semantics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``Mode`` now accepts either strings or ``Family`` members and normalizes the
stored family to its string value. Radial order ``m`` must be at least one;
invalid families and orders raise ``ValidationError`` instead of relying on
Python assertions.

Unsupported cutoff or normalized-parameter geometries raise
``UnsupportedGeometryError``. Non-finite wavelength derivatives raise
``ConvergenceError``. Structured solver methods continue to represent an
unguided or unconverged mode with ``SolverResult(converged=False)``; scalar
effective-index and cutoff conveniences continue to return ``NaN`` for that
expected absence.

Material loading
~~~~~~~~~~~~~~~~

Fiber loading no longer depends on ``PyOptik.MaterialBank``. Use the public
material protocol and registry for fixed or dispersive materials::

   from PyFiberModes import ConstantIndex, MaterialRegistry, load_fiber

   materials = MaterialRegistry({"custom_glass": ConstantIndex(value=1.47 * ureg.dimensionless)})
   fiber = load_fiber("custom_fiber", wavelength=1.55e-6 * ureg.meter, materials=materials)

The built-in ``fused_silica`` and ``silica`` names continue to work in YAML
fiber definitions.

Structured solver failures
~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``Fiber.solve_effective_index(mode)`` when failure details matter. It
returns ``SolverResult`` with ``converged``, ``value``, ``residual``,
``iterations``, ``bracket``, and ``message``. The historical
``Fiber.get_effective_index(mode)`` convenience method remains available and
returns ``NaN`` on failure.

Solver namespace names
~~~~~~~~~~~~~~~~~~~~~~

Internal solver namespaces were renamed for clarity:

.. list-table::
   :header-rows: 1

   * - Former internal name
     - Current internal name
   * - ``solver.ssif``
     - ``solver.two_layer``
   * - ``solver.tlsif``
     - ``solver.three_layer``
   * - ``solver.mlsif``
     - ``solver.multilayer``
   * - ``NeffSolver``
     - ``EffectiveIndexSolver``

These solver modules remain internal. Application code should normally use
``Fiber.get_effective_index`` or ``Fiber.solve_effective_index``.

Plotting dependency
~~~~~~~~~~~~~~~~~~~

Matplotlib is now optional. Install ``PyFiberModes[plotting]`` when using
``Field.plot``. The plotting module is loaded lazily. TypedUnit also depends on Matplotlib for
quantity-aware plotting support.

0.13
----

Python 3.10 support ended in 0.13.2. PyFiberModes now requires Python 3.11 or
newer. The misspelled group-velocity aliases and obsolete solver compatibility
wrappers were removed before the stable public API policy was introduced.
