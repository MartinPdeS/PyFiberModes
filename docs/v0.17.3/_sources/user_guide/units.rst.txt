Units and conventions
=====================

PyFiberModes accepts bare numbers and uses SI units internally. Unit-aware
objects are not converted automatically, so convert values before passing them
to the package.

.. list-table:: Quantity conventions
   :header-rows: 1
   :widths: 28 25 47

   * - Quantity
     - Unit
     - Example
   * - Wavelength, radius, distance
     - m
     - ``1550e-9`` for 1550 nm
   * - Refractive index, numerical aperture
     - dimensionless
     - ``1.444`` or ``0.14``
   * - Propagation constant :math:`\beta`
     - rad/m
     - ``fiber.get_propagation_constant(mode)``
   * - Phase or group velocity
     - m/s
     - ``fiber.get_group_velocity(mode)``
   * - Group-velocity dispersion
     - s\ :sup:`2`/m
     - ``fiber.get_group_velocity_dispersion(mode)``
   * - Dispersion
     - ps/(nm km)
     - ``fiber.get_dispersion(mode)``
   * - Effective area
     - m\ :sup:`2`
     - Multiply by ``1e12`` for :math:`\mu\mathrm{m}^2`

Geometry convention
-------------------

Layer radii are interface radii measured from the fiber axis. The final layer
is the exterior medium and is treated as extending outward; its stored outer
radius is not a second finite interface. Layers are indexed from the innermost
layer in the constructed :class:`~PyFiberModes.fiber.Fiber`, even though YAML
files are written from the outside toward the center for readability and
compatibility with the loader.

Mode convention
---------------

A :class:`~PyFiberModes.mode.Mode` combines a family, azimuthal order ``nu``,
and radial order ``m``. Radial numbering starts at one. TE and TM families use
``nu = 0``; HE and EH hybrid modes use positive azimuthal orders. LP modes are
the scalar weak-guidance representation and should not be mixed with vector
families when counting physically distinct modes without accounting for their
correspondence and degeneracy.

Angles passed as ``phi`` and ``theta`` to field component methods are in
radians. ``phi`` selects the azimuthal phase of the modal pattern; ``theta``
rotates the transverse observation axes.
