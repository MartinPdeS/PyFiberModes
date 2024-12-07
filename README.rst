.. list-table::
   :widths: 10 25 25
   :header-rows: 0

   * - Meta
     - |python|
     - |docs|
   * - Testing
     - |ci/cd|
     - |coverage|
   * - PyPi
     - |PyPi|
     - |PyPi_download|
   * - Anaconda
     - |anaconda|
     - |anaconda_download|


PyFiberModes
============

PyFiberModes is a Python package designed to simulate the propagation of modes in optical fibers. It supports all circularly symmetric geometries and provides tools for in-depth analysis of optical fiber properties and modal characteristics.

Key Features
------------
- **Comprehensive Simulations**: Supports step-index fibers, multilayer fibers, and custom refractive index profiles.
- **Extensive Mode Analysis**: Calculate mode profiles, propagation constants, effective indices, and more.
- **Open Source**: Designed for customization and extensibility.

----

Documentation
**************
The latest documentation is always available `here <https://martinpdes.github.io/PyFiberModes/>`_ or via the badge below:

|docs|

----


Installation
************

Using Pip
---------

Install the PyFiberModes package directly from PyPi. Ensure you have Python 3.10 or later:

.. code-block:: bash

   pip install PyFiberModes


Manual Installation
-------------------

To manually install the package:

.. code-block:: bash

   git clone https://github.com/MartinPdeS/PyFiberModes.git
   cd PyFiberModes
   pip install .


----

Testing
*******

PyFiberModes includes comprehensive tests. To run the tests locally:

.. code-block:: bash

   git clone https://github.com/MartinPdeS/PyFiberModes.git
   cd PyFiberModes
   pip install .
   coverage run --source=PyFiberModes --module pytest --verbose tests
   coverage report --show-missing


This will generate a coverage report detailing untested portions of the code.

----


Examples and Usage
******************
Learn how to use PyFiberModes with detailed examples in the `Examples Section <https://martinpdes.github.io/PyFiberModes>`_ of the documentation. Examples include:
- Calculating the effective index for a given mode.
- Visualizing mode profiles in optical fibers.
- Evaluating fiber parameters like V-number and dispersion.

----

Get Involved
************
PyFiberModes is an actively maintained project, and contributions are highly encouraged! Whether it's a bug fix, feature request, or enhancement, your input is valuable.

- Report issues or request features on the `GitHub Issue Tracker <https://github.com/MartinPdeS/PyFiberModes/issues>`_.
- Open pull requests to improve the codebase.

Contact
*******
If you would like to collaborate, please reach out:

Author: `Martin Poinsinet de Sivry-Houle <https://github.com/MartinPdS>`_

Email: `martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet-de-sivry@polymtl.ca?subject=PyFiberModes>`_

----

.. |python| image:: https://img.shields.io/pypi/pyversions/pyfibermodes.svg
   :target: https://www.python.org/

.. |docs| image:: https://github.com/martinpdes/pyfibermodes/actions/workflows/deploy_documentation.yml/badge.svg
   :target: https://martinpdes.github.io/PyFiberModes/
   :alt: Documentation Status

.. |coverage| image:: https://raw.githubusercontent.com/MartinPdeS/PyFiberModes/python-coverage-comment-action-data/badge.svg
   :alt: Unittest coverage
   :target: https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html

.. |PyPi| image:: https://badge.fury.io/py/PyFiberModes.svg
   :target: https://pypi.org/project/PyFiberModes/

.. |PyPi_download| image:: https://img.shields.io/pypi/dm/PyFiberModes.svg
   :target: https://pypistats.org/packages/pyfibermodes

.. |ci/cd| image:: https://github.com/martinpdes/pyfibermodes/actions/workflows/deploy_coverage.yml/badge.svg
   :target: https://martinpdes.github.io/PyFiberModes/actions
   :alt: Unittest Status

.. |anaconda_download| image:: https://anaconda.org/martinpdes/pyfibermodes/badges/downloads.svg
   :alt: Anaconda downloads
   :target: https://anaconda.org/martinpdes/pyfibermodes

.. |anaconda| image:: https://anaconda.org/martinpdes/pyfibermodes/badges/version.svg
   :alt: Anaconda version
   :target: https://anaconda.org/martinpdes/pyfibermodes
