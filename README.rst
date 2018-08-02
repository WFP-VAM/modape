wsmtk
=====

The **w**\ hittaker **sm**\ oothing **t**\ ool\ **k**\ it combines a state-of-the art whittaker smoother, implemented as fast C-extension through Cython and including a V-curve optimization of the smoothing parameter, with a HDF5 based processing chain optimized for MODIS data.

The sub-module ``wsmtk.whittaker`` includes the following variations of the whittaker smoother with 2nd order differences:

- Whittaker with fixed smoothing parameter (``s``)
- Whittaker with V-curve optimization of the smoothing parameter (``s``)
- Whittaker with V-curve optimization of the smoothing parameter (``s``) and expectile smoothing using asymmetric weights

The MODIS processing chain consists of the following executables, which can be called through commandline:

- ``downloadMODIS``: Query and download raw MODIS products (requires Earthdata credentials)
- ``processMODIS``: Collect raw MODIS data into daily datacubes stored in an HDF5 file
- ``smoothMODIS``: Smooth, gapfill and interpolate raw MODIS data using the implemented whittaker smoother
- ``windowMODIS``: Extract mosaic(s) of multiple MODIS tiles, or subset(s) of a global/tiled MODIS product and export it as GeoTIFF raster in WGS1984 coordinate system

Additional executables:

- ``smoothCSV``: Smooth timeseries stored within a CSV file
- ``smoothRTS``: Smooth a series of raster files stored in a local directory

Installation
------------
**Dependencies:**

wsmtk depends on these packages:

- numpy
- gdal
- h5py
- beautifulsoup4
- requests
- progress
- pandas

Some of these packages (eg. GDAL) can be difficult to build, especially on windows machines. In the latter case it's advisable to download an unofficial binary wheel from `Christoph Gohlke's Unofficial Windows Binaries for Python Extension Packages <https://www.lfd.uci.edu/~gohlke/pythonlibs/>`_ and install it locally with ``pip install`` before installing wsmtk.

**Installation from github:**

.. code:: bash

    $ git clone https://github.com/WFP-VAM/wsmtk
    $ cd wsmtk
    $ pip install .

**Installation from PyPi:**

.. code:: bash

    $ pip install wsmtk

Usage tutorial
-----

All executables can be called with a ``-h`` flag for detailed usage.

For a more detailed tutorial on how to use the executables, please visit `WFP-VAM.github.io/wsmtk <http://WFP-VAM.github.io/wsmtk>`_.


CHANGES
-----

TBD: Initial release

TODO
-----
-
