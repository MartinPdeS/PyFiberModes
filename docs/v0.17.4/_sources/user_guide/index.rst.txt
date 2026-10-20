User guide
==========

This guide connects the short examples to the complete API reference. It
follows the lifecycle of a fiber calculation: define the structure, choose and
solve modes, validate sampled fields, explore parameters, and preserve enough
context to reproduce the result.

.. toctree::
   :maxdepth: 2

   units
   fibers
   modes
   fields
   sweeps
   reproducibility

The package is designed for circularly symmetric, concentric step-index
structures. It is not a general-purpose solver for non-circular cores,
longitudinally varying cross-sections, stress birefringence, or arbitrary
anisotropic materials. Those cases require a numerical full-vector solver;
PyFiberModes can still provide reference solutions or input modes for a larger
workflow.
