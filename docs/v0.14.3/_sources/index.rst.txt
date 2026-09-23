**Date**: |today|, **Version**: |version|

.. image:: _static/logo.png
   :width: 620
   :align: center
   :alt: PyFiberModes logo

.. include:: ../../README.rst
    :start-line: 5

.. toctree::
    :maxdepth: 2
    :hidden:

    Getting started <getting_started.rst>
    User guide <user_guide/index.rst>
    Theory <theory.rst>
    Architecture <architecture.rst>
    API policy <api_policy.rst>
    Migration <migration_guide.rst>
    Workflows <workflows.rst>
    Performance <performance.rst>
    Troubleshooting <troubleshooting.rst>
    Examples <examples.rst>
    API reference <code.rst>
    References <references.rst>

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

PyFiberModes uses SI units throughout. Wavelengths and radii are in metres,
propagation constants are in radians per metre, and the dispersion helper
returns ps/(nm km). See :doc:`user_guide/units` for a complete conventions
table.
