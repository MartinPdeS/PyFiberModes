Performance
===========

Field maps exploit cylindrical symmetry: analytical radial fields are solved
once on a radial line, interpolated over the two-dimensional grid, and cached.
Request several components through ``Field.get_components`` so they share the
same radial pass and azimuthal factors.

Practical guidance
------------------

* Reuse a ``Field`` instance when evaluating several derived quantities.
* Request only the components required by an analysis.
* Use modest grids while exploring, then verify final results at higher
  resolution.
* Bound ``max_nu`` and ``max_m`` during automatic mode discovery.
* Prefer ``fiber.sweep`` to hand-written mutation loops; it restores state and
  returns dense, structured arrays.

See :doc:`examples` for reproducible scaling examples.
