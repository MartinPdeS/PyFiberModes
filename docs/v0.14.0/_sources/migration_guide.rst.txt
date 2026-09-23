.. _migration_guide:

Migration guide
===============

This page records user-visible name and behavior changes. Changes to internal
solver modules are included when they are likely to affect advanced users,
even though those modules are outside the stable public API.

Unreleased
----------

Material loading
~~~~~~~~~~~~~~~~

Fiber loading no longer depends on ``PyOptik.MaterialBank``. Use the public
material protocol and registry for fixed or dispersive materials::

   from PyFiberModes import ConstantIndex, MaterialRegistry, load_fiber

   materials = MaterialRegistry({"custom_glass": ConstantIndex(1.47)})
   fiber = load_fiber("custom_fiber", wavelength=1.55e-6, materials=materials)

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
``Field.plot``. Computational use no longer installs a plotting stack.

0.13
----

Python 3.10 support ended in 0.13.2. PyFiberModes now requires Python 3.11 or
newer. The misspelled group-velocity aliases and obsolete solver compatibility
wrappers were removed before the stable public API policy was introduced.
