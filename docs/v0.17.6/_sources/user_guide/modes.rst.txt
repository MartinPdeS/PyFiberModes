Modes and modal properties
==========================

Choosing a family
-----------------

LP modes are appropriate for weakly guiding fibers, where the refractive-index
contrast is small and a scalar description is sufficient. HE, EH, TE, and TM
modes retain the vector family structure and are the better choice for higher
contrast or when longitudinal field components matter.

Create arbitrary identifiers with :class:`~PyFiberModes.mode.Mode`, or import
common instances from the package root:

.. code-block:: python

   from PyFiberModes import HE11, LP01, Mode

   lp21 = Mode("LP", nu=2, m=1)
   assert LP01.nu == 0

Solving one mode
----------------

The effective index is the central solved quantity. Other longitudinal modal
properties derive from it and its wavelength dependence:

.. code-block:: python

   from PyFiberModes import LP01, load_fiber

   fiber = load_fiber("SMF28", wavelength=1.55e-6)
   values = {
       "effective index": fiber.get_effective_index(LP01),
       "propagation constant": fiber.get_propagation_constant(LP01),
       "group index": fiber.get_group_index(LP01),
       "dispersion": fiber.get_dispersion(LP01),
   }

A guided solution should have an effective index between the exterior index
and the largest layer index. Close to cutoff, repeat the calculation at nearby
wavelengths because derivative-based group and dispersion quantities become
sensitive to step size and a mode may disappear.

Discovering supported modes
---------------------------

Use :meth:`~PyFiberModes.fiber.Fiber.find_modes` when the supported set is not
known in advance:

.. code-block:: python

   modes = fiber.find_modes(families=("LP",), max_nu=4, max_m=3)
   for mode in modes:
       print(mode, fiber.get_effective_index(mode))

The bounds describe candidate orders, not the number of returned modes. Keep
them close to the physical range of interest: every candidate requires a root
search. Search scalar and vector families separately to avoid interpreting two
representations of the same weak-guidance mode as independent solutions.

Cutoff and normalized quantities
--------------------------------

For two-layer fibers, :meth:`~PyFiberModes.fiber.Fiber.V_number` describes the
normalized frequency and :meth:`~PyFiberModes.fiber.Fiber.get_normalized_beta`
places the propagation constant between the core and cladding limits.
``get_mode_cutoff_v0`` and ``get_mode_cutoff_wavelength`` expose the analytical
cutoff supported by the applicable solver. See :doc:`../theory` for equations
and assumptions.
