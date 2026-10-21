Workflows
=========

Choose a workflow based on the question you are asking.

.. grid:: 2
   :gutter: 2

   .. grid-item-card:: Characterize a known fiber

      Load a YAML definition, discover guided modes, then sweep wavelength for
      effective index, group index, and dispersion.

   .. grid-item-card:: Inspect modal fields

      Build a :class:`~PyFiberModes.field.Field`, evaluate components in one
      batch, and compute power, confinement, or overlap.

   .. grid-item-card:: Explore a design space

      Use structured sweeps for maps and ``optimize_fiber`` for bounded inverse
      design without mutating the source fiber.

   .. grid-item-card:: Propagate coupled modes

      Assemble a Hermitian coupling matrix and integrate constant or
      longitudinally varying coupled-mode equations.

Recommended numerical checks
----------------------------

* Repeat field calculations with a larger domain and finer sampling.
* Inspect solutions close to cutoff with more than one wavelength step.
* Compare LP and vector families only within their applicable guidance regime.
* Preserve the exact package version and fiber definition with published data.
