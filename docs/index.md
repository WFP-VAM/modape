---
title: Document Center
---

# Whittaker smoothing toolkit

The Whittaker smoothing toolkit (wsmtk) is a python based smoothing engine for EO-data, combining a state-of-the art Whittaker smoother with a fast processing chain for handling MODIS data.

## Main features


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

## Installation

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


## Executables

The wsmtk package features some executables, which are mainly for handling and smoothing MODIS data, as well as smoothing arbitrary raster timeseries and pixel timeseries inside CSV files.

### downloadMODIS

**Description:**

Query and download MODIS data using product ID(s).

For querying/downloading tile-based MODIS product, a region of interest as `--roi` is required. This can be either a point location, with latitude and longitude, or a bounding box, with lower left x/y and upper right x/y coordinates.

Querying and downloading global MODIS products (5km) does not require a region of interest.

If both MOD and MYD products are required, using M?D will query/download both MODIS Aqua and Terra.

Multiple product IDs can be queried and downloaded at the same time. They need to be separated by space (e.g. MOD13A2 MOD11A1)

Downloading data requires valid Earthdata user credentials (register at https://urs.earthdata.nasa.gov/users/new). The download can be performed with Python's request module (default), or with WGET by adding `--wget` (WGET must be available in PATH - this option can be better for global or a big number of products).

The products are downloaded by default to the current working directory, if no target directory has been supplied by specifying `-d, --targetdir`.

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

### processMODIS

**Description:**

Collect raw MODIS hdf files into raw daily HDF5 files for further processing.

After downloading raw MODIS files, they need to be collected into an HDF5 file for smoothing.

During this step, the composite files are converted back to daily data. If information about the compositing day-of-the-year (DOY) is available (e.g. MOD/MYD13A2, the observation is inserted at the respective position. If such information is not available, the observation is set to the midpoint of the compositing period.

For MODIS vegetation products, 16-day TERRA (MOD) and AQUA (MYD) products are interleaved to a combined 8-day product (MXD).

For continuous processing chain, the target directory should remain constant. If the processed HDF5 file already exists, the new raw data is ingested at the correct place, otherwise a new file is created.

The HDF5 files are created in a subdirectory of the target directory (default current working directory), VIM and LTD for vegetation indices and land surface temperature, respectively.

The only mandatory input is `srcdir`, the directory containing the raw hdf files to be processed.

Currently processing is implemented for MODIS vegetation and MODIS LST products. By default, NDVI and LST_Day will be extracted and processed from the raw hdf files. Additionally all parameters (EVI, LST_Night) can be extracted too by adding `--all-parameters`.

The default blocksize for processing is set to 120 by 120. When specifying a custom blocksize, make sure the number of rows and columns divides evenly by the respective blocksize!


**processMODIS help:**

```
$ processMODIS -h

usage: processMODIS [-h] [-d] [-c] [--all-parameters] [-b ] srcdir

Process downloaded RAW MODIS hdf files

positional arguments:
  srcdir               directory with raw MODIS .hdf files

optional arguments:
  -h, --help           show this help message and exit
  -d , --targetdir     Target directory for PROCESSED MODIS files (default is
                       scrdir)
  -c , --compression   Compression for HDF5 files
  --all-parameters     Flag to process all possible VAM parameters
  -b  , --blocksize    Minimum values for row & columns per processing block
                       (default 120 120)
```
**Usage example:**

![processMODIS]

### smoothMODIS

**Description:**

Smooth, gapfill and interpolate raw MODIS HDF5 files.

As mandatory input, raw MODIS HDF5 files (as processed by `processMODIS`) are smoothed, gapfilled and interpolated, using a whittaker filter with 2nd order differences.

By default, the standard whittaker smoother will be applied, and if no `s` value is provided (to be provided as `log10(s)` - e.g. 0.1), the function will read a pixel based `s` value from the s-grid.

If the `--soptimze` flag is set, the `s` value will be per-pixel optimized using the V-curve, and the resoluting values stored as s-grid.

If `-p, --pvalue` is set, expectile smoothing with V-curve optimization is performed. The results of the optimized `s` are again stored as s-grid.

The smoothed dataset can be interpolated to any integer timestep by passing `-t, --tempint` to the function. If no value is set, no interpolation is performed and the smoothed data will have the same temporal resolution as the raw data.

If a smoothed HDF5 already exists in the targetdirectory, it will be updated, otherwise a new one will be created.

If not all raw timesteps should be used for smoothing, or if not all smooted timesteps should be written back to disk, this can be controlled by specifying `-n, --nsmooth` and `-u, --nupdate` respectively.

By default, the smoothing will be processed in serial. To leverage parallel processing, add the `--parallel` flag, and optionally specify `--nworkers`.

Naming convention for smoothed HDF5 files:

_ProductID . Tile ID . Version number . VAM parameter . Temporal flag_

e.g. MXD13A2.h18v04.006.VIM.txn.h5

The temporal flag can have the following values:
- txn: native temporal resolution
- txp: pentadal (5-daily) data
- txd: dekadal (10 daily) data
- txc: custom interpolation (user input)

**smoothMODIS help:**

```
$ smoothMODIS -h

usage: smoothMODIS [-h] [-s] [-S  [...]] [-t] [-n] [-u] [-p] [-d]
                   [--soptimize] [--parallel] [--nworkers]
                   RAW HDF5

Smooth, gapfill and interpolate processed raw MODIS HDF5 files

positional arguments:
  RAW HDF5              Raw MODIS HDF5 file

optional arguments:
  -h, --help            show this help message and exit
  -s , --svalue         S value for smoothing (has to be log10(s)
  -S  [ ...], --srange  [ ...]
                        S value range for V-curve (float log10(s) values as
                        smin smax sstep - default 0.0 4.0 0.1)
  -t , --tempint        Value for temporal interpolation (integer required -
                        default is native temporal resolution AKA no
                        interpolation)
  -n , --nsmooth        Number of raw timesteps used for smoothing
  -u , --nupdate        Number of smoothed timesteps to be updated in HDF5
                        file
  -p , --pvalue         Value for asymmetric smoothing (float required)
  -d , --targetdir      Target directory for smoothed output
  --soptimize           Use V-curve for s value optimization
  --parallel            Parallel processing
  --nworkers            Number of workers used for parallel processing
                        (default is number of cores available minus 1)
```
**Usage example:**

![smoothMODIS]

### windowMODIS

**Description:**

Create mosaics (or subsets) from smoothed MODIS HDF5 files.

Given a directory containing smoothed MODIS HDF5 files `path`, this tool creates mosaics or subsets of smoothed data, depending on the input.

The desired product can be specified with the `-p, --product flag`.

For tile-based products, specifying a region of interest is mandatory. If a point location is specified, the entire intersecting tile will be used as window. If a bounding box is specified, multiple tiles of the same product will be mosaicked. If the bounding box is smaller than one tile, a subset of the intersecting tile will be created.

For global MODIS products, no region of interest is required. If a bounding box is specified, a subset of the global product is returned as window.

The time range required can be controlled with `-b, --begin-date` and `-e, --end-date` where the input should be a date in the format YYYYMM (e.g. 201808).

Adding the `--sgrid` flag will extract or mosaic the s-grid rather than the smoothed data.

All mosaics/subsets will be saved to the target directory (default current working directory) as GeoTIFF files, in a geographic coordinate system (EPSG:4326 - WGS84) with the naming convention:

_REGION PARAMETER yyyyjddd_

e.g. regvim2018j169.tif

where REGION should be a three letter code for the region (can be defined with `--region`, default is "reg"), PARAMETER is the VAM parameter code and yyyyjddd is the julian date of the image.


**windowMODIS help:**

```
$ windowMODIS -h

usage: windowMODIS [-h] [-p] [--roi ROI [ROI ...]] [--region] [-b] [-e]
                   [--parameter] [-d] [--sgrid]
                   path

Extract a window from MODIS products

positional arguments:
  path                 Path to processed MODIS h5 files

optional arguments:
  -h, --help           show this help message and exit
  -p , --product       MODIS product ID (can be parial match with *)
  --roi ROI [ROI ...]  Region of interest. Can be LAT/LON point or bounding
                       box in format llx lly urx ury
  --region             region 3 letter region code (default is "reg")
  -b , --begin-date    Start date (YYYYMM)
  -e , --end-date      End date (YYYYMM)
  --parameter          VAM parameter code
  -d , --targetdir     Target directory for GeoTIFFs (default current
                       directory)
  --sgrid              Extract (mosaic of) s value grid(s))
```
**Usage example:**

![windowMODIS]


### infoMODIS

**Description:**

Retrieves and displays metadata for raw/smooth MODIS HDF5 files.

Information includes:

- Dimensions
- Type
- Startdate
- Enddate
- NoData value
- Temporal resolution
- Spatial resolution
- Last modified
- Last smoothing run and parameters

**infoMODIS help:**

```
$ infoMODIS -h

usage: infoMODIS [-h] file

Get MODIS raw/smooth file info

positional arguments:
  file        Full path to MODIS h5 file

optional arguments:
  -h, --help  show this help message and exit

```
**Usage example:**

![infoMODIS]

### producttableMODIS

**Description:**

Display the MODIS product table, optionally filter by product, resolution and parameter.

**producttableMODIS help:**

```
$ producttableMODIS -h

usage: producttableMODIS [-h] [--product] [--resolution] [--parameter]

MODIS product table

optional arguments:
  -h, --help     show this help message and exit
  --product      MODIS product ID
  --resolution   Filter for pixel size
  --parameter    Filter for VAM parameter

```
**Usage example:**

![producttableMODIS]

### smoothRTS

**Description:**

Smooth and gapfill local raster timeseries which is located in `path`.

If no `-P, --pattern` is specified, the function assumes that only rasters beloging to the timeseries are in the `path` directory.

The user has to make sure that the rasters, which define the timeseries, have a file naming that puts them in the right order when applying a simple sort to the filenames.

If `-S, --srange` or `--soptimize` is specified, the whittaker filter with V-curve will be used. If no s-range value is supplied, the default range 0 to 4 with 0.1 step will be applied.

Nodata values can be specified using the `--nodata` flag. If none is specified, the function will try to read it from input. If that fails, nodata is assumed to be 0.

The resulting smoothed raster files have the same naming as the input files, and will be saved to disk in a subfolder corresponding to the smoothing method applied, in the target directory (default is current working directory).

Additionally, a config.txt file is saved to disk with the smoothing parameters.

If s-optimization is selected, the resulting s-grid is also saved to disk.

The subdirectories are:

- _./filt0_ for fixed s
- _./filtvc_ for s-optimization
- _./filtvcp_ for s-optimization and expectile smoothing

The smoothing currently only runs sequential, so big timeseries or big raster datasets can take some time!

**smoothRTS help:**

```
$ smoothRTS -h

usage: smoothRTS [-h] [-P] [-d] [-s] [-S  [...]] [-p] [-b] [--nodata]
                 [--soptimize]
                 path

Extract a window from MODIS products

positional arguments:
  path                  Path containing raster files

optional arguments:
  -h, --help            show this help message and exit
  -P , --pattern        Pattern to filter file names
  -d , --targetdir      Target directory for GeoTIFFs (default current
                        directory)
  -s , --svalue         S value for smoothing (has to be log10(s)
  -S  [ ...], --srange  [ ...]
                        S range for V-curve (float log10(s) values as smin
                        smax sstep - default 0.0 4.0 0.1)
  -p , --pvalue         Value for asymmetric smoothing (float required)
  -b , --blocksize      Processing block side length (default 256)
  --nodata              NoData value
  --soptimize           Use V-curve for s value optimization

```

**Usage example:**

![smoothRTS]


### smoothCSV

**Description:**

Smooth pixel-timeseries from a CSV file.

The CSV file for `file` is expected to have a specific naming convention and structure:

The first three letters of the filename are expected to be the region code and will be included in the output name.

The first row is skipped by default.

The first column is reserved for metadata, the first 3 rows of the first columns are "ID", "Lon", "Lat".

All columns starting number two are considered separate pixel timeseries, where each timeseries should start with an ID, and with the 2 coordinates lon and lat.

Example:

| -      	| -          	| -          	| -          	|
|--------	|------------	|------------	|------------	|
| ID     	| Location 1 	| Location 2 	| Location n 	|
| Lon    	| 25.65      	| 26.75      	| ...        	|
| Lat    	| 7.92       	| 7.78       	| ...        	|
| Date 1 	| 0.62       	| 0.62       	| ...        	|
| Date 2 	| 0.62       	| 0.56       	| ...        	|
| Date 3 	| 0.65       	| 0.5        	| ...        	|
| ...    	| ...        	| ...        	| ...        	|

The output is saved to CSV in the same directory as the input.

The filename of the output is combined of region code (first three letters) and smoothing method:

- _filt0.csv_ : for fixed s
- _filtvc.csv_ : for s-optimized
- _filtvcp.csv_ : for s-optimized with expectile

The optimal s (sopt) and log10(sopt) are appended to each column, if a fixed s is chosen, it's assumed to be sopt.

If expectile smoothing is chosen, the p value is also be appended.

**smoothCSV help:**

```
$ smooothCSV -h

usage: smoothCSV [-h] [-s] [-S  [...]] [-p] file

Smooth CSV file

positional arguments:
  file                  CSV file

optional arguments:
  -h, --help            show this help message and exit
  -s , --svalue         S value for smoothing (has to be log10(s)
  -S  [ ...], --srange  [ ...]
                        S range for V-curve (float log10(s) values as smin
                        smax sstep - default 0.0 4.0 0.1)
  -p , --pvalue         Value for asymmetric smoothing (float required)
```
**Usage example:**


![smoothCSV]

---

Last updated on 2018-08-03

Terminal recordings created with [asciinema](https://github.com/asciinema/asciinema/blob/develop/LICENSE)

created by:

Valentin Pesendorfer, 2018  
valentin.pesendofer@wfp.org




[downloadMODIS]: img/download.gif
[processMODIS]: img/process.gif
[smoothMODIS]: img/smooth.gif
[windowMODIS]: img/window.gif
[infoMODIS]: img/info.gif
[producttableMODIS]: img/ptable.gif
[smoothRTS]: img/smoothRTS.gif
[smoothCSV]: img/smoothCSV.gif
