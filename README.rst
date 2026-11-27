PyFiberModes
============

.. image:: https://raw.githubusercontent.com/MartinPdeS/PyFiberModes/HEAD/docs/source/_static/logo.png
   :alt: PyFiberModes logo
   :width: 100%

.. list-table::
   :width: 100%
   :widths: 35 65
   :header-rows: 1

   * - Badge
     - Status
   * - Python versions
     - |mps-python|
   * - Documentation
     - |mps-docs|
   * - Tests
     - |mps-ci| |mps-coverage|
   * - PyPI
     - |mps-pypi| |mps-downloads|
   * - Conda
     - |mps-conda| |mps-conda-downloads| |mps-conda-release|

Guided-mode analysis and design for circular optical fibers.

Find the modes supported by a circular optical fiber at a given wavelength. Compare
effective indices, inspect field profiles, and follow modes as the wavelength or
fiber geometry changes.

See the `usage guide <https://github.com/MartinPdeS/PyFiberModes/blob/HEAD/USAGE.rst>`__ for further examples and implementation details.

Installation
------------

Requires Python ``>=3.11``.

.. code-block:: console

   $ python -m pip install PyFiberModes

Install from the ``martinpdes`` Conda channel:

.. code-block:: console

   $ conda install -c martinpdes pyfibermodes

For an editable source installation:

.. code-block:: console

   $ git clone https://github.com/MartinPdeS/PyFiberModes.git
   $ cd PyFiberModes
   $ python -m pip install -e .

Discover the guided modes of a fiber
------------------------------------

.. code-block:: python

   """
   Guided-mode discovery
   =====================

   Discover the supported linearly polarized modes and compare their effective
   indices at one wavelength.
   """
   from PyFiberModes.units import ureg

   import matplotlib.pyplot as plt

   from PyFiberModes.analysis import find_modes
   from PyFiberModes.fiber import load_fiber

   fiber = load_fiber("SMF28", wavelength=1310e-9 * ureg.meter)
   modes = find_modes(fiber, families=("LP",), max_nu=3, max_m=3)
   effective_indices = [fiber.get_effective_index(mode) for mode in modes]

   figure, axis = plt.subplots()
   axis.bar([str(mode) for mode in modes], effective_indices)
   axis.set(xlabel="Guided mode", ylabel="Effective index", title="SMF-28 at 1310 nm")
   axis.ticklabel_format(axis="y", useOffset=False)
   plt.show()

Documentation
-------------

Use the `PyFiberModes documentation <https://martinpdes.github.io/PyFiberModes/>`__ to work through
examples and look up the inputs and return values of individual functions.

Citation
--------

Cite the version used in your calculations with the `software citation metadata <https://github.com/MartinPdeS/PyFiberModes/blob/HEAD/CITATION.cff>`__.

Contributing
------------

Bug reports, worked examples, and documentation corrections are welcome.
To set up a development checkout and run the tests, follow the `contributing guide <https://github.com/MartinPdeS/PyFiberModes/blob/HEAD/CONTRIBUTING.md>`__.

Contact
-------

For questions, contact `Martin Poinsinet de Sivry-Houle <mailto:martin.poinsinet.de.sivry@gmail.com?subject=PyFiberModes>`__.

License
-------

See the `PyFiberModes license <https://github.com/MartinPdeS/PyFiberModes/blob/HEAD/LICENSE>`__ for terms of use.

.. |mps-python| image:: https://img.shields.io/pypi/pyversions/pyfibermodes.svg
   :alt: Supported Python versions
   :target: https://pypi.org/project/pyfibermodes/

.. |mps-pypi| image:: https://img.shields.io/pypi/v/pyfibermodes.svg
   :alt: PyPI version
   :target: https://pypi.org/project/pyfibermodes/

.. |mps-downloads| image:: https://img.shields.io/pypi/dm/pyfibermodes.svg
   :alt: Monthly PyPI downloads
   :target: https://pypistats.org/packages/pyfibermodes

.. |mps-docs| image:: https://github.com/MartinPdeS/PyFiberModes/actions/workflows/deploy_documentation.yml/badge.svg
   :alt: Documentation
   :target: https://martinpdes.github.io/PyFiberModes/

.. |mps-ci| image:: https://github.com/MartinPdeS/PyFiberModes/actions/workflows/tests.yml/badge.svg
   :alt: Tests
   :target: https://github.com/MartinPdeS/PyFiberModes/actions/workflows/tests.yml

.. |mps-coverage| image:: https://raw.githubusercontent.com/MartinPdeS/PyFiberModes/python-coverage-comment-action-data/badge.svg
   :alt: Test coverage
   :target: https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html

.. |mps-conda| image:: https://anaconda.org/martinpdes/pyfibermodes/badges/version.svg
   :alt: Anaconda version
   :target: https://anaconda.org/martinpdes/pyfibermodes

.. |mps-conda-downloads| image:: https://anaconda.org/martinpdes/pyfibermodes/badges/downloads.svg
   :alt: Anaconda downloads
   :target: https://anaconda.org/martinpdes/pyfibermodes

.. |mps-conda-release| image:: https://anaconda.org/martinpdes/pyfibermodes/badges/latest_release_relative_date.svg
   :alt: Latest Anaconda release
   :target: https://anaconda.org/martinpdes/pyfibermodes
