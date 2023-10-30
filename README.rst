PyFiberModes
============

|python|
|docs|
|Citation|
|Unittest|
|PyPi|
|PyPi_download|


..  figure:: https://github.com/MartinPdeS/PyFiberModes/blob/master/docs/images/mode_propagation.gif?raw=true
   :alt: some image
   :class: with-shadow float-left
   :width: 800px

   Propagation of mode in an adiabatic 2x1 modally-specific photonic lantern.




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

To install manually (os independent) you will need to install:

1. cmake (3.16+)

In order to use the PyFiberModes Simulator Library, one must have installed the python dependencies:

Then, download and install the PyFiberModes package:

.. code-block:: python

    >>> git clone https://github.com/MartinPdeS/PyFiberModes.git
    >>> cd PyFiberModes
    >>> pip install -r requirements/requirements.txt
    >>> pip install .

----


Contact Information
*******************

As of 2021 the project is still under development if you want to collaborate it would be a pleasure. I encourage you to contact me.

PyMieSim was written by `Martin Poinsinet de Sivry-Houle <https://github.com/MartinPdS>`_  .

Email:`martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet-de-sivry@polymtl.ca?subject=PyFiberModes>`_ .


.. |python| image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
   :target: https://www.python.org/

.. |docs| image:: https://readthedocs.org/projects/pyfibermodes/badge/?version=latest
   :target: https://pyfibermodes.readthedocs.io/en/latest/
   :alt: Documentation Status

.. |Citation| image:: https://zenodo.org/badge/366930899.svg
   :target: https://zenodo.org/badge/latestdoi/366930899

.. |Unittest| image:: https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/MartinPdeS/2cac8ebb51cd9ce07dc7a955648301d5/raw/d75a89f19acc11302d81dacefe5be207beee24a8/PyFiberModes_coverage_badge.json

.. |PyPi| image:: https://badge.fury.io/py/PyFiberModes.svg
   :target: https://pypi.org/project/PyFiberModes/

.. |PyPi_download| image:: https://img.shields.io/pypi/dm/PyFiberModes.svg
   :target: https://pypi.org/project/PyFiberModes/



