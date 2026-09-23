Public API contract
===================

This page defines the supported PyFiberModes 1.0 surface. Only names imported
from the package root and listed in ``PyFiberModes.__all__`` are stable. Public
signatures and important defaults are recorded in ``tests/api_snapshot.json``
and checked on every supported Python version.

Core models
-----------

.. list-table:: Stable models
   :header-rows: 1
   :widths: 24 28 48

   * - Object
     - Construction
     - Contract
   * - :class:`~PyFiberModes.mode.Mode`
     - ``Mode(family, nu, m)``
     - Immutable identifier. ``family`` accepts a string or
       :class:`~PyFiberModes.mode.Family`; ``nu >= 0`` and ``m >= 1``.
   * - :class:`~PyFiberModes.fiber.Fiber`
     - ``Fiber(wavelength, ...)``
     - Mutable physical model; wavelength and radii are metres. Direct layer
       mutation requires ``clear_caches()``.
   * - :class:`~PyFiberModes.models.LayerSpec`
     - ``LayerSpec(name, radius, refractive_index)``
     - Immutable validated layer input; radius is metres and index is
       dimensionless.
   * - :class:`~PyFiberModes.models.SolverSettings`
     - ``SolverSettings(...)``
     - Immutable numerical settings. Tolerances must be positive.
   * - ``CartesianCoordinates``, ``CylindricalCoordinates``
     - Coordinate arrays
     - Coordinate distances are metres and angular coordinates are radians.

Fiber construction and analysis
-------------------------------

.. list-table:: Stable operations
   :header-rows: 1
   :widths: 26 28 20 26

   * - Operation
     - Result
     - Mutates input?
     - Failure
   * - :func:`~PyFiberModes.fiber.load_fiber`
     - Initialized ``Fiber``
     - No
     - ``ValidationError``, ``KeyError``, or file-loading exception for invalid
       definitions
   * - ``fiber.solve_effective_index(mode)``
     - ``SolverResult[float]``; effective index is dimensionless
     - No
     - Failed numerical or unguided solve is represented by
       ``converged=False``
   * - ``fiber.get_effective_index(mode)``
     - Effective index or ``NaN``
     - No
     - ``NaN`` only for an unsuccessful numerical/unguided solve
   * - ``fiber.get_mode_cutoff_v0(mode)``
     - Dimensionless cutoff V-number or ``NaN``
     - No
     - ``UnsupportedGeometryError`` when the geometry has no cutoff backend
   * - ``fiber.get_mode_field(mode, limit=None, n_point=101)``
     - Sampled ``Field``
     - No
     - ``ValidationError`` for invalid input; ``ConvergenceError`` when the
       mode cannot be solved
   * - :func:`~PyFiberModes.analysis.find_modes`
     - Tuple of guided modes
     - No
     - Expected non-convergence is omitted; invalid input propagates
   * - :func:`~PyFiberModes.analysis.sweep_modes`
     - ``ModeSweepResult``
     - No; copied fibers are varied
     - Expected per-point convergence failure becomes ``NaN``; invalid input
       and unsupported operations propagate
   * - :func:`~PyFiberModes.optimization.optimize_fiber`
     - ``DesignResult``
     - No; candidate fibers are copies
     - Invalid configuration raises ``ValueError`` or ``ValidationError``

Modal-property names and units
------------------------------

.. list-table:: ``Fiber`` property façade
   :header-rows: 1
   :widths: 36 25 39

   * - Name
     - Unit
     - Meaning
   * - ``get_effective_index``
     - dimensionless
     - Modal effective refractive index
   * - ``get_normalized_beta``
     - dimensionless
     - Normalized propagation constant
   * - ``get_propagation_constant``
     - rad/m
     - Longitudinal propagation constant
   * - ``get_phase_velocity``, ``get_group_velocity``
     - m/s
     - Phase and group velocities
   * - ``get_group_index``
     - dimensionless
     - Group refractive index
   * - ``get_group_velocity_dispersion``
     - s\ :sup:`2`/m
     - Second angular-frequency derivative of propagation constant
   * - ``get_dispersion``
     - ps/(nm km)
     - Chromatic dispersion
   * - ``get_S_parameter``
     - ps/(nm\ :sup:`2` km)
     - Dispersion slope; historical spelling retained for compatibility
   * - ``V_number``
     - dimensionless
     - Normalized frequency; historical property spelling retained
   * - ``M_number``
     - dimensionless
     - Approximate mode count :math:`V^2/2`; historical property spelling
       retained

Failure model
-------------

The exception hierarchy has stable meanings:

* :class:`~PyFiberModes.exceptions.ValidationError` means the request or
  physical input is invalid. It is also a ``ValueError``.
* :class:`~PyFiberModes.exceptions.UnsupportedGeometryError` means the input is
  valid but that calculation has no backend for its geometry.
* :class:`~PyFiberModes.exceptions.ConvergenceError` means a supported numerical
  calculation failed to produce a finite result.
* :class:`~PyFiberModes.exceptions.SolverError` is the common base for supported
  solver failures.

``SolverResult.unwrap()`` returns a converged value and otherwise raises
``ConvergenceError``. Derivative-based scalar methods also raise
``ConvergenceError`` rather than returning an unexplained non-finite value.
Convenience effective-index and cutoff methods retain ``NaN`` for an expected
absence of a root so array workflows can represent mode cutoff naturally.

Fields, design, and propagation
-------------------------------

``Field`` and coordinate objects expose sampled arrays and derived integrals;
the sampling grid is immutable by convention, while internal arrays are normal
NumPy arrays. ``DesignResult``, ``PropagationResult``, ``ModeSweepResult``, and
``SolverResult`` are result containers. Coupled-mode propagation does not
modify initial amplitudes, and overlap/coupling operations do not mutate their
fields.

Compatibility names
-------------------

``V_number``, ``M_number``, and ``get_S_parameter`` predate the current naming
policy. They remain supported public ``Fiber`` members for the 1.x series.
New APIs use descriptive snake-case names, but these three spellings are not
scheduled for removal and do not emit deprecation warnings.
