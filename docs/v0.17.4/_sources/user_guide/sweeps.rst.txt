Sweeps, design, and propagation
===============================

Wavelength sweeps
-----------------

Structured sweeps preserve the source fiber and return dense arrays. Rows
correspond to parameter values and columns to modes:

.. code-block:: python

   import numpy as np
   from PyFiberModes import LP01, LP11, load_fiber

   fiber = load_fiber("SMF28", wavelength=1.55e-6)
   wavelengths = np.linspace(1.30e-6, 1.65e-6, 71)
   result = fiber.sweep(
       wavelengths,
       modes=(LP01, LP11),
       metrics=("effective_index", "group_index", "dispersion"),
   )

   lp01 = result.for_mode(LP01)
   table_rows = result.as_records()

Unguided or failed solutions are ``NaN``. Do not silently interpolate through
these values: they may mark a physical cutoff. If ``modes`` is omitted, mode
discovery runs across the parameter grid and can be bounded with ``families``,
``max_nu``, and ``max_m``.

Custom parameter sweeps
-----------------------

For a parameter other than wavelength, provide a name and a setter operating
on each copied fiber:

.. code-block:: python

   def set_core_radius(candidate, radius):
       candidate.layers[0].radius_out = radius
       candidate.layers[1].radius_in = radius
       candidate.clear_caches()

   radii = np.linspace(3.5e-6, 5.0e-6, 31)
   result = fiber.sweep(
       radii,
       parameter_name="core_radius",
       setter=set_core_radius,
       modes=(LP01,),
       metrics=("effective_index",),
   )

Inverse design
--------------

:func:`~PyFiberModes.optimization.optimize_fiber` minimizes a scalar objective
over bounded design parameters. Constraint callables are feasible when their
return value is non-negative:

.. code-block:: python

   from PyFiberModes import DesignParameter, optimize_fiber

   radius = DesignParameter.layer_radius(0, bounds=(3.5e-6, 5.0e-6))
   design = optimize_fiber(
       fiber,
       parameters=(radius,),
       objective=lambda candidate: abs(candidate.get_dispersion(LP01)),
       constraints=(
           lambda candidate: candidate.get_effective_index(LP01)
           - candidate.last_layer.refractive_index,
       ),
       options={"seed": 7},
   )

Inspect ``design.success``, ``design.message``, and the unpenalized physical
constraints before accepting a result. Repeat stochastic optimization with
several seeds and validate the returned ``design.fiber`` independently.

Coupled-mode propagation
------------------------

:class:`~PyFiberModes.propagation.CoupledModeSystem` integrates
:math:`dA/dz=-i(\mathrm{diag}(\beta)+\kappa)A`. Subtracting a common propagation
constant reduces rapid common phase rotation without changing modal powers:

.. code-block:: python

   from PyFiberModes import CoupledModeSystem

   betas = np.array([fiber.get_propagation_constant(LP01),
                     fiber.get_propagation_constant(LP11)])
   relative_betas = betas - betas[0]
   coupling = np.array([[0.0, 0.25], [0.25, 0.0]])
   system = CoupledModeSystem(relative_betas, coupling)
   z = np.linspace(0.0, 20.0, 501)
   result = system.propagate([1.0, 0.0], z)

``coupling`` may be a constant matrix or a callable ``coupling(z)``. Use a
Hermitian matrix for lossless reciprocal exchange; check that the sum of
``result.powers`` remains constant to the integration tolerance.
