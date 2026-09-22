Theory and conventions
======================

PyFiberModes solves guided modes of circularly symmetric, concentric
step-index fibers. In each homogeneous layer, Maxwell's equations reduce to
Bessel or modified-Bessel radial solutions. Interface continuity determines
the allowed propagation constants.

Mode families
-------------

``LP`` modes use the weak-guidance approximation. ``HE`` and ``EH`` are hybrid
vector modes, while ``TE`` and ``TM`` are the transverse families. A mode is
identified by its family, azimuthal order ``nu``, and radial order ``m``.

Normalized frequency
--------------------

For a two-layer fiber, the normalized frequency is

.. math::

   V = k_0 a \sqrt{n_\mathrm{core}^2 - n_\mathrm{clad}^2},

where :math:`a` is the core radius and :math:`k_0=2\pi/\lambda`. Guided
solutions satisfy :math:`n_\mathrm{clad} < n_\mathrm{eff} < n_\mathrm{core}`.

Field normalization
-------------------

Field amplitudes are modal shapes. Derived overlap, confinement, effective
area, Poynting vector, and guided-power methods apply consistent numerical
normalization on the sampled Cartesian grid. Increase both the window and the
grid resolution when checking convergence.
