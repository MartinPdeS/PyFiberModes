Validation and reproducibility
==============================

Numerical output is most useful when its assumptions can be reconstructed.
Record the following with published or shared results:

* the exact PyFiberModes version and Python version;
* the fiber YAML file or every layer radius, index, and material model;
* vacuum wavelength and selected mode family and orders;
* field window, grid resolution, and polarization angles;
* sweep spacing or derivative settings near cutoff;
* optimization method, bounds, constraints, seed, and termination message.

Minimum validation checklist
----------------------------

#. Confirm ``n_cladding < n_effective <= n_max`` for every guided mode.
#. Repeat field integrals with a larger window and a finer grid.
#. Repeat derivative-based quantities with a different wavelength sampling.
#. Check expected limiting behavior, symmetry, and conservation laws.
#. Compare a simple two-layer case with an analytical reference before using a
   complicated multilayer structure.

Failure handling
----------------

Automatic discovery and sweeps represent unavailable modes with omission or
``NaN`` so that a broad exploration can continue. A single requested mode may
instead raise :class:`~PyFiberModes.exceptions.SolverError` or
:class:`~PyFiberModes.exceptions.ConvergenceError`. Catch these specific
exceptions only when a missing point is an expected outcome:

.. code-block:: python

   import numpy as np
   from PyFiberModes import ConvergenceError, SolverError

   try:
       effective_index = fiber.get_effective_index(mode)
   except (SolverError, ConvergenceError):
       effective_index = np.nan

Do not replace all exceptions with ``NaN``: invalid geometry, units, or array
shapes should remain visible and be corrected.

Archiving results
-----------------

``ModeSweepResult.as_records`` produces plain row-oriented dictionaries useful
for CSV or dataframe export. Preserve the original parameter array as well as
the records, and save configuration and environment metadata beside numerical
data. Cite the software using :doc:`../references` and use a tagged release,
not a moving branch, when preparing long-lived results.
