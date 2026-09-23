Defining fibers and materials
=============================

Bundled definitions
-------------------

The simplest route is :func:`~PyFiberModes.fiber.load_fiber`. Pass a bundled
YAML filename without its extension and always provide the wavelength when a
layer uses a dispersive material:

.. code-block:: python

   from PyFiberModes import load_fiber

   fiber = load_fiber("SMF28", wavelength=1.55e-6)
   for layer in fiber.layers:
       print(layer.name, layer.radius_in, layer.radius_out, layer.refractive_index)

The bundled catalog lives in ``PyFiberModes/fiber_files``. Definitions may use
an explicit ``index``, a registered ``material``, or a numerical aperture
relative to the adjacent outer layer. A compact two-layer definition is:

.. code-block:: yaml

   # Layers are listed from outside to inside.
   layers:
     cladding:
       name: cladding
       radius: 62.5e-6
       material: fused_silica
     core:
       name: core
       radius: 4.1e-6
       NA: 0.14

Custom material models
----------------------

Use :class:`~PyFiberModes.materials.MaterialRegistry` when a definition names a
material that is not built in. A model can be any callable receiving wavelength
in metres and returning refractive index:

.. code-block:: python

   from PyFiberModes import MaterialRegistry, load_fiber

   materials = MaterialRegistry({
       "my_glass": lambda wavelength: 1.46 - 0.003 * (wavelength / 1e-6 - 1.0),
   })
   fiber = load_fiber("my_fiber", wavelength=1.55e-6, materials=materials)

Programmatic construction
--------------------------

For a single structure, create a :class:`~PyFiberModes.fiber.Fiber`, add layers
from the center outward, and initialize it once all interfaces are defined:

.. code-block:: python

   from PyFiberModes import Fiber

   fiber = Fiber(wavelength=1.55e-6)
   fiber.add_layer(name="core", radius=4.1e-6, index=1.450)
   fiber.add_layer(name="cladding", radius=62.5e-6, index=1.444)
   fiber.initialize_layers()

Use :class:`~PyFiberModes.factory.FiberFactory` when radii or indices are
iterables and you want the Cartesian product of several designs.

Mutation and caches
-------------------

Call :meth:`~PyFiberModes.fiber.Fiber.update_wavelength` when changing
wavelength; it updates dependent state and clears cached modal and radial-field
solutions. After directly changing a layer radius or refractive index, call
:meth:`~PyFiberModes.fiber.Fiber.clear_caches`. Prefer
:meth:`~PyFiberModes.fiber.Fiber.sweep` or the optimization helpers for repeated
calculations because they work on copies and keep the input fiber unchanged.
