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

Load the bundled SMF28 step-index model at a vacuum wavelength of 1550 nm.
The scalar fundamental mode, LP01, has an effective refractive index between
the cladding and core indices.

.. code-block:: python

   from PyFiberModes import LP01, load_fiber
   from PyFiberModes.units import ureg

   fiber = load_fiber("SMF28", wavelength=1550 * ureg.nanometer)
   neff = fiber.get_effective_index(LP01)
   assert fiber.last_layer.refractive_index < neff < fiber.maximum_index
   print(f"{neff.magnitude:.6f}")

For this fiber definition the result is approximately **1.447491**.
The effective index is dimensionless. To extract a number from a result
with dimensions, first convert it to the unit you intend to report.

Check the sampled mode area
---------------------------

The effective area depends on how much of the field fits inside the window
and how finely the field is sampled. Use the same fiber for these three grids.

.. code-block:: python

   for half_width, samples in [(10, 101), (20, 201), (20, 401)]:
       field = fiber.get_mode_field(
           mode=LP01, limit=half_width * ureg.micrometer, n_point=samples
       )
       area = field.get_effective_area().to("micrometer ** 2")
       print(half_width, samples, round(area.magnitude, 3))

The first two grids keep the spacing at 0.2 µm and double the window width.
The area increases from about 64.417 to 64.519 µm² as the window includes more
of the evanescent tail. The last grid keeps the window fixed and halves the
spacing, giving about 64.515 µm². Repeat these checks for the mode and
accuracy you need, especially near cutoff.

Where to go next
----------------

Continue with :doc:`user_guide/fibers` to define a fiber,
:doc:`user_guide/modes` to calculate modal properties, or
:doc:`user_guide/fields` to analyze vector fields. Use :doc:`examples` for
downloadable examples and :doc:`code` for the complete API.
