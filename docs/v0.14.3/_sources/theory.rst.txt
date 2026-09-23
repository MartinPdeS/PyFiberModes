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

The effective index and propagation constant are related by

.. math::

   \beta = k_0 n_\mathrm{eff}, \qquad k_0 = \frac{2\pi}{\lambda}.

Under weak guidance, a vector-mode group may be represented by one LP mode.
LP labels are consequently an approximation and not an additional independent
set of physical solutions alongside HE, EH, TE, and TM labels.

Normalized frequency
--------------------

For a two-layer fiber, the normalized frequency is

.. math::

   V = k_0 a \sqrt{n_\mathrm{core}^2 - n_\mathrm{clad}^2},

where :math:`a` is the core radius and :math:`k_0=2\pi/\lambda`. Guided
solutions satisfy :math:`n_\mathrm{clad} < n_\mathrm{eff} < n_\mathrm{core}`.

The normalized propagation constant is commonly written

.. math::

   b = \frac{n_\mathrm{eff}^2-n_\mathrm{clad}^2}
            {n_\mathrm{core}^2-n_\mathrm{clad}^2},

so guided two-layer solutions lie between zero and one. Multilayer fibers use
the actual ordered interfaces and do not in general reduce to a single
core-cladding :math:`V` value.

Wavelength derivatives
----------------------

Group index, group-velocity dispersion, chromatic dispersion, and dispersion
slope depend on numerical derivatives of the solved propagation constant.
They are less robust than a single effective-index evaluation near cutoff or
an avoided crossing. Validate them by changing the wavelength sampling and by
checking that the same modal branch is followed.

Field normalization
-------------------

Field amplitudes are modal shapes. Derived overlap, confinement, effective
area, Poynting vector, and guided-power methods apply consistent numerical
normalization on the sampled Cartesian grid. Increase both the window and the
grid resolution when checking convergence.
