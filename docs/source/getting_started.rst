Getting started
===============

PyFiberModes requires Python 3.11 or newer. Create an isolated environment,
install the package from PyPI, and verify the interpreter that will run your
analysis:

.. code-block:: bash

   python -m pip install PyFiberModes
   python -c "import PyFiberModes; print(PyFiberModes.__version__)"

Install the optional plotting dependency if you intend to render field maps:

.. code-block:: bash

   python -m pip install "PyFiberModes[plotting]"

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

``LP01`` identifies the scalar fundamental mode. ``limit`` is the positive
half-width of the square sampling window in metres, and ``n_point`` is the
number of samples along each Cartesian axis.

Check the result
----------------

A useful calculation includes a physical and a numerical check:

.. code-block:: python

   n_cladding = fiber.last_layer.refractive_index
   n_core = fiber.maximum_index
   assert n_cladding < neff < n_core

   coarse = fiber.get_mode_field(LP01, limit=15e-6, n_point=151)
   fine = fiber.get_mode_field(LP01, limit=20e-6, n_point=251)
   relative_change = abs(
       fine.get_effective_area() / coarse.get_effective_area() - 1
   )
   print(f"effective-area change: {relative_change:.2%}")

Increase the window until the field decays before its boundary, then increase
the resolution until the derived quantity changes less than the accuracy your
application requires.

Where to go next
----------------

Continue with :doc:`user_guide/fibers` to define a fiber,
:doc:`user_guide/modes` to calculate modal properties, or
:doc:`user_guide/fields` to analyze vector fields. Use :doc:`examples` for
downloadable examples and :doc:`code` for the complete API.
