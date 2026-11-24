Sweeps, design, and propagation
===============================

Wavelength sweeps
-----------------

Structured sweeps preserve the source fiber and return dense arrays. Rows
correspond to parameter values and columns to modes:

.. code-block:: python

   from PyFiberModes.units import ureg

   import numpy as np
   from PyFiberModes import LP01, LP11, load_fiber

   fiber = load_fiber("SMF28", wavelength=1.55e-6 * ureg.meter)
   wavelengths = np.linspace(1.30e-6, 1.65e-6, 71) * ureg.meter
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

   from PyFiberModes.units import ureg

   def set_core_radius(candidate, radius):
       candidate.layers[0].radius_out = radius
       candidate.layers[1].radius_in = radius
       candidate.clear_caches()

   radii = np.linspace(3.5e-6, 5.0e-6, 31) * ureg.meter
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

   from PyFiberModes.units import ureg

   from PyFiberModes import DesignParameter, optimize_fiber

   radius = DesignParameter.layer_radius(0, bounds=(3.5e-6 * ureg.meter, 5.0e-6 * ureg.meter))
   design = optimize_fiber(
       fiber,
       parameters=(radius,),
       objective=lambda candidate: abs(candidate.get_dispersion(LP01).to("ps / nm / km").magnitude),
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

For a phase-matched pair, choose zero relative propagation constants and a
real, symmetric coupling matrix. With light launched into mode 1, its power
fraction is :math:`\cos^2(\kappa z)` and the fraction in mode 2 is
:math:`\sin^2(\kappa z)`. This gives an independent check on both the launch
condition and the exchange distance.

.. code-block:: python

   import numpy as np
   from PyFiberModes import CoupledModeSystem
   from PyFiberModes.units import ureg

   kappa = 0.25 / ureg.meter
   betas = np.zeros(2) / ureg.meter
   coupling = np.array([[0, 1], [1, 0]]) * kappa
   system = CoupledModeSystem(betas=betas, coupling=coupling)
   transfer_length = np.pi / (2 * kappa)
   z = np.linspace(0, 1, 101) * transfer_length
   result = system.propagate(initial=[1, 0], z=z)
   np.testing.assert_allclose(result.powers[-1].magnitude, [0, 1], atol=1e-7)

Each row of ``result.amplitudes`` corresponds to a position and each column
to a mode, ordered as in ``betas``. For this normalized launch, the sum of
each row of ``result.powers`` is one within the integration tolerance.
The first complete transfer occurs about 6.28 m from launch.

For detuned modes, supply their individual propagation constants. Subtracting
a common constant removes a common phase and preserves the modal powers.
For varying coupling, pass a function that takes a length quantity and returns
a finite square matrix with inverse-length units at each position. Supply at
least two strictly ordered positions. Decreasing positions propagate backward
from the initial amplitudes at the first sample.
