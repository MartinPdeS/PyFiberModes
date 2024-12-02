PyFiberModes
============

|python|
|docs|
|coverage|
|PyPi|
|PyPi_download|




This project aims to develop an useful tool to simulate propagating mode in fiber optics for all kind of circular-symmetric geometries.

----

Documentation
**************
All the latest available documentation is available `here <https://pyfibermodes.readthedocs.io/en/latest/>`_ or you can click the following badge:

|docs|


----


Installation
------------


Pip installation
****************

The package have been uploaded as wheel for a few OS (Linux, MacOS) and need Python 3.10.
As such, with the adequate configuration one can simply do

.. code-block:: python

   >>> pip3 install PyFiberModes



Manual installation
*******************
The following shell commands should do the trick.

.. code-block:: python

    >>> git clone https://github.com/MartinPdeS/PyFiberModes.git
    >>> cd PyFiberModes
    >>> pip install -r requirements/requirements.txt
    >>> pip install .

----

Testing
*******

To test localy (with cloning the GitHub repository) you'll need to install the dependencies and run the coverage command as

.. code:: python

   >>> git clone https://github.com/MartinPdeS/PyFiberModes.git
   >>> cd PyFiberModes
   >>> pip install -r requirements/requirements.txt
   >>> coverage run --source=PyFiberModes --module pytest --verbose tests
   >>> coverage report --show-missing

----



Coding examples
***************
Plenty of examples are available online, I invite you to check the `examples <https://pyfibermodes.readthedocs.io/en/master/gallery/index.html>`_
section of the documentation.


----


Contact Information
*******************

As of 2023 the project is still under development if you want to collaborate it would be a pleasure. I encourage you to contact me.

PyMieSim was written by `Martin Poinsinet de Sivry-Houle <https://github.com/MartinPdS>`_  .

Email:`martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet-de-sivry@polymtl.ca?subject=PyFiberModes>`_ .


.. |python| image:: https://img.shields.io/pypi/pyversions/pyfibermodes.svg
   :target: https://www.python.org/

.. |docs| image:: https://readthedocs.org/projects/pyfibermodes/badge/?version=latest
   :target: https://pyfibermodes.readthedocs.io/en/latest/
   :alt: Documentation Status

.. |coverage| image:: https://raw.githubusercontent.com/MartinPdeS/PyFiberModes/python-coverage-comment-action-data/badge.svg
   :alt: Unittest coverage
   :target: https://htmlpreview.github.io/?https://github.com/MartinPdeS/PyFiberModes/blob/python-coverage-comment-action-data/htmlcov/index.html

.. |PyPi| image:: https://badge.fury.io/py/PyFiberModes.svg
   :target: https://pypi.org/project/PyFiberModes/

.. |PyPi_download| image:: https://img.shields.io/pypi/dm/PyFiberModes.svg
   :target: https://pypistats.org/packages/pyfibermodes




