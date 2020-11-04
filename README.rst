MODAPE
=====

|CI| |version| |pyversions| |downloads| |license| |documentation|

The **MOD**\ IS **A**\ ssimilation and **P**\ rocessing **E**\ ngine combines a state-of-the art Whittaker smoother, implemented as fast C-extension through Cython and including a V-curve optimization of the smoothing parameter, with a HDF5 based processing chain optimized for MODIS data.

The sub-module ``modape.whittaker`` includes the following variations of the Whittaker smoother with 2nd order differences:

- **ws2d**: Whittaker with fixed smoothing parameter (``s``)
- **ws2dp**: Whittaker with fixed smoothing parameter (``s``) and expectile smoothing using asymmetric weights
- **ws2doptv**: Whittaker with V-curve optimization of the smoothing parameter
- **ws2doptvp**: Whittaker with V-curve optimization of the smoothing parameter and expectile smoothing using asymmetric weights

The MODIS processing chain consists of the following executables, which can be called through commandline:

- ``modis_download``: Query and download raw MODIS products (requires Earthdata credentials)
- ``modis_collect``: Collect raw MODIS data into daily datacubes stored in an HDF5 file
- ``modis_smooth``: Smooth, gapfill and interpolate raw MODIS data using the implemented Whittaker smoother
- ``modis_window``: Extract mosaic(s) of multiple MODIS tiles, or subset(s) of a global/tiled MODIS product and export it as GeoTIFF raster in WGS1984 coordinate system

Additional executables:

- ``csv_smooth``: Smooth timeseries stored within a CSV file
- ``modis_info``: Retrieve metadata from created HDF5 files

For a more information please check out the `documentation <https://wfp-vam.github.io/modape/>`_!

Installation
------------
**Dependencies:**

modape depends on these packages:

- click
- gdal
- h5py
- numpy
- pandas
- python-cmr
- requests

Some of these packages (eg. GDAL) can be difficult to build, especially on windows machines. In the latter case it's advisable to download an unofficial binary wheel from `Christoph Gohlke's Unofficial Windows Binaries for Python Extension Packages <https://www.lfd.uci.edu/~gohlke/pythonlibs/>`_ and install it locally with ``pip install`` before installing modape.

**Installation from github:**

.. code:: bash

    $ git clone https://github.com/WFP-VAM/modape
    $ cd modape
    $ pip3 install .

**Installation from PyPi:**

.. code:: bash

    $ pip3 install modape


Bugs, typos & feature requests
-----

If you find a bug, see a typo, have some kind of troubles running the module or just simply want to have a feature added, please `submit an issue! <https://github.com/WFP-VAM/modape/issues/new>`_


-----

References:

P. H. C. Eilers, V. Pesendorfer and R. Bonifacio, "Automatic smoothing of remote sensing data," 2017 9th International Workshop on the Analysis of Multitemporal Remote Sensing Images (MultiTemp), Brugge, 2017, pp. 1-3.
doi: 10.1109/Multi-Temp.2017.8076705
URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8076705&isnumber=8035194

Core Whittaker function adapted from ``whit2`` function from `R` package `ptw <https://cran.r-project.org/package=ptw>`_:

Bloemberg, T. G. et al. (2010) "Improved Parametric Time Warping for Proteomics", Chemometrics and Intelligent Laboratory Systems, 104 (1), 65-74

Wehrens, R. et al. (2015) "Fast parametric warping of peak lists", Bioinformatics, in press.


.. |CI| image:: https://github.com/WFP-VAM/modape/workflows/build/badge.svg
             :target: https://github.com/WFP-VAM/modape/actions/

.. |version| image:: https://img.shields.io/pypi/v/modape.svg
                  :target: https://pypi.org/project/modape/

.. |pyversions| image:: https://img.shields.io/pypi/pyversions/modape.svg
                     :target: https://pypi.org/project/modape/

.. |downloads| image:: https://img.shields.io/pypi/dm/modape.svg
                    :target: https://pypi.org/project/modape/

.. |license| image:: https://img.shields.io/github/license/WFP-VAM/modape.svg
                  :target: https://github.com/WFP-VAM/modape/blob/master/LICENSE

.. |documentation| image:: https://img.shields.io/badge/documentation-passing-brightgreen
   :target: https://wfp-vam.github.io/modape/
|
