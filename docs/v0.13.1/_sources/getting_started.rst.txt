Getting started
===============

Install PyFiberModes from PyPI and verify the interpreter that will run your
analysis:

.. code-block:: bash

   python -m pip install PyFiberModes
   python -c "import PyFiberModes; print(PyFiberModes.__version__)"

Your first guided mode
----------------------

Fiber dimensions and wavelengths use SI units. Load a bundled definition,
solve the fundamental LP mode, and inspect its effective index:

.. code-block:: python

   from PyFiberModes import LP01
   from PyFiberModes.fiber import load_fiber

   fiber = load_fiber("SMF28", wavelength=1550e-9)
   neff = fiber.get_effective_index(LP01)
   field = fiber.get_mode_field(LP01, limit=15e-6, n_point=151)

   print(neff)
   print(field.get_effective_area())

Where to go next
----------------

Use :doc:`theory` for conventions, :doc:`workflows` for task-oriented guides,
:doc:`examples` for executable notebooks, and :doc:`code` for the complete API.
