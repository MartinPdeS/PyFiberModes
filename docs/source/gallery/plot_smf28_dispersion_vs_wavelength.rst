
.. DO NOT EDIT.
.. THIS FILE WAS AUTOMATICALLY GENERATED BY SPHINX-GALLERY.
.. TO MAKE CHANGES, EDIT THE SOURCE PYTHON FILE:
.. "gallery/plot_smf28_dispersion_vs_wavelength.py"
.. LINE NUMBERS ARE GIVEN BELOW.

.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        :ref:`Go to the end <sphx_glr_download_gallery_plot_smf28_dispersion_vs_wavelength.py>`
        to download the full example code

.. rst-class:: sphx-glr-example-title

.. _sphx_glr_gallery_plot_smf28_dispersion_vs_wavelength.py:


Modal dispersion VS core index
==============================

.. GENERATED FROM PYTHON SOURCE LINES 8-10

Imports
~~~~~~~

.. GENERATED FROM PYTHON SOURCE LINES 10-15

.. code-block:: python3

    import numpy
    from PyFiberModes import LP01
    from PyFiberModes.fiber import load_fiber
    from MPSPlots.render2D import SceneList








.. GENERATED FROM PYTHON SOURCE LINES 16-19

Generating the fiber structures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Here we create the different fiber design that we want to explore

.. GENERATED FROM PYTHON SOURCE LINES 19-39

.. code-block:: python3

    wavelength_list = numpy.linspace(1000e-9, 1800e-9, 100)
    data = []
    for wavelegnth in wavelength_list:
        smf28 = load_fiber(fiber_name='SMF28', wavelength=wavelegnth)
        dispersion = smf28.get_effective_index(mode=LP01)
        data.append(dispersion)


    figure = SceneList()

    ax = figure.append_ax(
        x_label=r'Wavelength [$\mu m$]',
        y_label='Dispersion',
        x_scale_factor=1e6
    )

    ax.add_scatter(x=wavelength_list, y=data)

    _ = figure.show()
    # -



.. image-sg:: /gallery/images/sphx_glr_plot_smf28_dispersion_vs_wavelength_001.png
   :alt: plot smf28 dispersion vs wavelength
   :srcset: /gallery/images/sphx_glr_plot_smf28_dispersion_vs_wavelength_001.png
   :class: sphx-glr-single-img






.. rst-class:: sphx-glr-timing

   **Total running time of the script:** (0 minutes 0.268 seconds)


.. _sphx_glr_download_gallery_plot_smf28_dispersion_vs_wavelength.py:

.. only:: html

  .. container:: sphx-glr-footer sphx-glr-footer-example




    .. container:: sphx-glr-download sphx-glr-download-python

      :download:`Download Python source code: plot_smf28_dispersion_vs_wavelength.py <plot_smf28_dispersion_vs_wavelength.py>`

    .. container:: sphx-glr-download sphx-glr-download-jupyter

      :download:`Download Jupyter notebook: plot_smf28_dispersion_vs_wavelength.ipynb <plot_smf28_dispersion_vs_wavelength.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
