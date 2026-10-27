Troubleshooting
===============

No effective-index solution
---------------------------

The requested mode may be below cutoff or outside the solver bracket. Confirm
the wavelength, layer ordering, radii, indices, and mode numbering. Search a
small candidate set with ``fiber.find_modes`` before requesting high orders.

Fields appear clipped
---------------------

Increase ``limit`` until the outer field decays well before the grid boundary,
then increase ``n_point`` until integrated quantities converge.

Dispersion is noisy near cutoff
-------------------------------

Dispersion uses numerical wavelength derivatives. Repeat the calculation with
several wavelength spacings and do not interpret points where a mode ceases to
be guided.

Documentation examples fail locally
-----------------------------------

Install the documentation extra, select a non-interactive Matplotlib backend,
and rebuild from the repository root:

.. code-block:: bash

   python -m pip install -e ".[documentation]"
   MPLBACKEND=Agg make -C docs html
