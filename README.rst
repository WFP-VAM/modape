wsmtk
=====

The **w**\ hittaker **sm**\ oothing **t**\ ool\ **k**\ it

Implemented functions are:

-  Download MODIS data (earthdata account required)
-  Process downloaded data

   -  Extract subdatasets
   -  Store data in HDF5 files

-  Extract windows over one or more tiles, reproject to WGS1984 and save
   as GeoTIFF

These functions can be called from command-line, outside of the python
interpreter.

Installation
------------

You can install wsmtk from github with:

.. code:: bash

    $ git clone https://github.com/WFP-VAM/wsmtk
    $ cd wsmtk
    $ python setup.py install

Usage
-----

The following example show exemplary usage of the wsmtk functionality:



Get help
^^^^^^^^

All command-line functions can be called with a ``-h`` flag to view the
usage:

::

    $ downloadMODIS -h

    usage: downloadMODIS [-h] --roi ROI [ROI ...] [-c] [-b] [-e] [--username] [--password] [-d] [-v] [--download]
                         product [product ...]

    Query and download MODIS products (earthdata accound required for download)

    positional arguments:
      product              MODIS product ID(s)

    optional arguments:
      -h, --help           show this help message and exit
      --roi ROI [ROI ...]  Region of interest. Can be LAT/LON point or bounding box in format llx,lly,urx,ury
      -c , --collection    MODIS collection
      -b , --begin-date    Start date (YYYY-MM-DD)
      -e , --end-date      End date (YYYY-MM-DD)
      --username           Earthdata username (required for download)
      --password           Earthdata password (required for download)
      -d , --dest          Destination directory
      -v, --verbose        Destination directory
      --download           Download data

Download Data
^^^^^^^^^^^^^

.. code:: bash

    # query number of products available

    $ downloadMODIS MOD13A2 --roi 31.6931 10.869 34.093 12.562 --begin-date 2010-01-01

    Checking for MODIS products ...... done.
    190 results found.

    # downloading to destination directory (earthdata credentials required)

    $ downloadMODIS MOD13A2 --roi 31.6931 10.869 34.093 12.562 --begin-date 2010-01-01 --dest /data/modis/ --username user --password pass --download



Process Data
^^^^^^^^^^^^

.. code:: bash

    # Raw data downloaded to /data/modis/

    $ processMODIS /data/modis/ --prcdir /data/modis/processed/



Extract (sub-)windows from processed data as GeoTIFF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash

    # bounding box for South Africa

    $ windowMODIS MOD13A2 --roi -2.746470 7.897562 7.610793 12.939542 --reg SAF --dataset Smooth --prcdir /data/modis/processed/ --targetdir /data/modis/tiff_extract/

    # output naming-convention is: regparYYYjDOY.tif (reg = region, par = VAM parameter)
