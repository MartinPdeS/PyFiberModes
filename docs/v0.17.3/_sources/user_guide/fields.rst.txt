Fields, power, and overlap
==========================

Sampling a field
----------------

A :class:`~PyFiberModes.field.Field` samples an analytical radial solution on a
square Cartesian grid:

.. code-block:: python

   from PyFiberModes import LP01, load_fiber

   fiber = load_fiber("SMF28", wavelength=1.55e-6)
   field = fiber.get_mode_field(LP01, limit=20e-6, n_point=251)
   components = field.get_components(("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"))

Use an odd ``n_point`` to include the fiber axis. ``limit`` is the positive
half-width, so the full domain is twice that value. Request related components
together with ``get_components`` and reuse the same ``Field`` instance; radial
solutions and interpolated arrays are cached.

Derived quantities
------------------

.. code-block:: python

   area_um2 = field.get_effective_area() * 1e12
   core_fraction = field.get_confinement_factor(radius=4.1e-6)
   sx, sy, sz = field.get_poynting_vector()
   guided_power = field.get_power()

Field amplitudes represent modal shapes rather than a launched absolute power.
Normalized overlap and confinement are therefore usually more portable than
raw amplitudes. Power-like integrals depend on the sampled window and grid.

Convergence procedure
---------------------

For any reported field integral:

#. Increase ``limit`` until the field is negligible at all four boundaries.
#. Hold that window fixed and increase ``n_point``.
#. Compare the final quantity, not only the image, between refinements.
#. Apply the same grid to both fields in an overlap calculation.

A smooth-looking plot can still have a biased integral if the domain clips the
evanescent tail. Near a material interface, resolution must also be sufficient
to represent the discontinuity.

Overlap and coupling
--------------------

The convenience method ``field_a.overlap(field_b)`` returns normalized complex
electric-field overlap. The standalone :func:`~PyFiberModes.propagation.overlap`
also supports an unnormalized result. Both fields must have identical grid
shapes and spacing. :func:`~PyFiberModes.propagation.coupling_matrix` constructs
a Hermitian matrix from several fields and accepts an optional scalar or
sampled spatial perturbation.

Plotting
--------

Install the ``plotting`` extra, then render one or more named components:

.. code-block:: python

   figure = field.plot(plot_type=["Ex", "Ey", "Emod"], show=False)
   figure.savefig("mode-field.png", dpi=200, bbox_inches="tight")

Use ``show=False`` in scripts, tests, and headless environments.
