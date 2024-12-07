
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

A numerical tool for simulating propagating modes in fiber optics, supporting all kinds of circular-symmetric geometries.


----

Documentation
**************
The latest documentation is available `here <https://martinpdes.github.io/PyFiberModes/>`_ or by clicking the badge below:

|docs|


----


Installation
************

Pip Installation
================
PyFiberModes is available as a Python wheel for Linux and macOS. It requires Python 3.10. To install, simply run:

.. code-block:: bash

   pip install PyFiberModes

Manual Installation
===================
For manual installation, clone the repository and install dependencies:

.. code-block:: bash

   git clone https://github.com/MartinPdeS/PyFiberModes.git
   cd PyFiberModes
   pip install -r requirements/requirements.txt
   pip install .

----

Testing
*******

To run tests locally after cloning the repository, follow these steps:

.. code-block:: bash

   git clone https://github.com/MartinPdeS/PyFiberModes.git
   cd PyFiberModes
   pip install .
   coverage run --source=PyFiberModes --module pytest --verbose tests
   coverage report --show-missing

This ensures the package is thoroughly tested and provides a coverage report.

----

Coding Examples
***************
Explore plenty of examples in the `examples section <https://martinpdes.github.io/PyFiberModes/docs/v0.5.9.3/gallery/index.html>`_ of the documentation.

----

Contact
*******

PyFiberModes is actively developed, and contributions are welcome! If you'd like to collaborate or provide feedback, please reach out.

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
