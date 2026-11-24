PyFiberModes documentation
==========================

.. image:: _static/logo.png
   :width: 620
   :align: center
   :alt: PyFiberModes logo

Calculate guided modes in circular step-index fibers, sample their electric
and magnetic fields, and propagate coupled modal amplitudes. The examples use
Pint quantities so you can enter a core radius in micrometers and a wavelength
in nanometers in the same calculation.

.. toctree::
    :maxdepth: 2
    :hidden:

    Start <getting_started.rst>
    Guide <guide.rst>
    Examples <examples.rst>
    API <code.rst>
    Resources <resources.rst>

Choose your path
----------------

.. grid:: 3
   :gutter: 2

   .. grid-item-card:: New to PyFiberModes?
      :link: getting_started
      :link-type: doc

      Install the package and complete a guided mode calculation in a few
      minutes.

   .. grid-item-card:: Solving a real problem?
      :link: user_guide/index
      :link-type: doc

      Follow task-oriented guides for fiber models, modal properties, fields,
      sweeps, optimization, and propagation.

   .. grid-item-card:: Looking up an object?
      :link: code
      :link-type: doc

      Browse the generated API reference, signatures, and source links.

See :doc:`user_guide/units` for input dimensions and output units, or
:doc:`troubleshooting` if a solve fails or a field integral changes with grid size.
