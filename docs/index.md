---
title: Document Center
---

# Whittaker smoothing toolkit

The Whittaker smoothing toolkit (wsmtk) is a python based smoothing engine for EO-data, combining a state-of-the art Whittaker smoother with a fast processing chain for handling MODIS data.

### Main features


Whittaker smoother:

- Whittaker with 2nd order differences
- fast implementation with Cython
- Pixel-wise optimization of the smoothing paramter using the V-curve
- V-curve optimization with expectile smoothing using asymmetric weights
- User defined temporal interpolation


MODIS processing chain:

- Easy query and download of MODIS products
- Collection of raw data into daily HDF5 files for fast processing
- Automatic interleaving of MOD13 and MYD13 products
- Fast smoothing, gapfilling and temporal interpolation
- Extraction of mosaics and subsets from smoothed data as GeoTIFFs in GCS WGS1984

### Installation

Dependencies:

wsmtk depends on these packages:

- numpy
- gdal
- h5py
- beautifulsoup4
- requests
- progress
- pandas

Some of these packages (eg. GDAL) can be difficult to build, especially on windows machines. In the latter case it's advisable to download an unofficial binary wheel from [Christoph Gohlke's Unofficial Windows Binaries for Python Extension Packages](https://www.lfd.uci.edu/~gohlke/pythonlibs/) and install it locally with ``pip install`` before installing wsmtk.

**Installation from github:**



    $ git clone https://github.com/WFP-VAM/wsmtk
    $ cd wsmtk
    $ pip install .

**Installation from PyPi:**



    $ pip install wsmtk


### Executables

The wsmtk package features some executables, which are mainly for handling and smoothing MODIS data, as well as smoothing arbitrary raster timeseries and pixel timeseries inside CSV files.

#### downloadMODIS

**Description:**

Query and download MODIS data using product ID(s).

For querying/downloading tile-based MODIS product, a region of interest as ```--roi``` is required. This can be either a point location, with latitude and longitude, or a bounding box, with lower left x/y and upper right x/y coordinates.

Querying and downloading global MODIS products (5km) does not require a region of interest.

If both MOD and MYD products are required, using M?D will query/download both MODIS Aqua and Terra.

Multiple product IDs can be queried and downloaded at the same time. They need to be separated by space (e.g. MOD13A2 MOD11A1)

Downloading data requires valid Earthdata user credentials (register at https://urs.earthdata.nasa.gov/users/new). The download can be performed with Python's request module (default), or with WGET by adding ```--wget``` (WGET must be available in PATH - this option can be better for global or a big number of products).

The products are downloaded by default to the current working directory, if no target directory has been supplied by specifying ```-d, --targetdir```.

**downloadMODIS help:**

```
$ downloadMODIS -h

usage: downloadMODIS [-h] [--roi ROI [ROI ...]] [-c] [-b] [-e] [--username]
                     [--password] [-d] [--download] [--wget]
                     product [product ...]

Query and download MODIS products (Earthdata account required for download)

positional arguments:
  product              MODIS product ID(s)

optional arguments:
  -h, --help           show this help message and exit
  --roi ROI [ROI ...]  Region of interest. Can be LAT/LON point or bounding
                       box in format llx,lly,urx,ury
  -c , --collection    MODIS collection
  -b , --begin-date    Start date (YYYY-MM-DD)
  -e , --end-date      End date (YYYY-MM-DD)
  --username           Earthdata username (required for download)
  --password           Earthdata password (required for download)
  -d , --targetdir     Destination directory
  --download           Download data
  --wget               Use WGET for downloading
```
**Usage example:**

![downloadMODIS]


[downloadMODIS]: img\download.gif
