# MODAPE

## Overview

The **MOD**IS **A**ssimilation and **P**rocessing **E**ngine combines a state-of-the-art whittaker smoother, implemented as fast C-extension through Cython, with a HDF5 based processing chain optimized for MODIS data.

The implementation of the Whittaker filter includes a V-curve based optimization of the smoothing parameter, which allows a pixel by pixel variation in the degree of smoothing applied, directly derived from the pixel’s timeseries. In addition, MODAPE implements expectile smoothing to estimate a smoothly varying envelope of the input signal.

Most users will want to use the out-of-the-box processing chain, which features pre-tuned parameters for running the smoother and command line executables, that enable the user to run the entire processing chain with minimal input from the command line.

The chain breaks down into 4 separate steps, each executable with a dedicated command line script, that handle everything from downloading raw MODIS data, to filtering the data and exporting it as GeoTIFF for downstream analysis and visualization:

<br>

<figure>
  <img src="img/overview.png" alt="overview">
  <br>
  <figcaption>MODIS processing chain in MODAPE</figcaption>
</figure>

<ol>
<li style='color: #2688FB;'>
  <span style='color:black;'>
  <code>modis_download</code>
  <br> Query and download MODIS data to local storage (requires <a href="https://urs.earthdata.nasa.gov/" target="_blank">Earthdata credentials</a>)
  </span>
</li>
<li style='color: #6FDD82;'>
  <span style='color:black;'>
  <code>modis_collect</code>
  <br> Collect raw NASA HDF files into raw HDF5 file(s)</a>
  </span>
</li>
</li>
<li style='color: #FED138;'>
  <span style='color:black;'>
  <code>modis_smooth</code>
  <br> Run Whittaker smoother on raw HDF5 file(s)</a>
  </span>
</li>
<li style='color: #FC4B5A;'>
  <span style='color:black;'>
  <code>modis_window</code>
  <br> Export data from raw/smooth HDF5 file(s) as GeoTIFF</a>
  </span>
</li>
</ol>


## Installation

MODAPE builds on [GDAL](https://gdal.org/) and **expects working python bindings**. Both can be a bit tricky to set up correctly, but here's a couple of tips for installing it poperly depending on the system:

- **Ubuntu**: Check out [UbuntuGIS](https://wiki.ubuntu.com/UbuntuGIS)
- **Windows**: Chris Gohlke has an amazing collection of unofficial binary wheels, one of them is [GDAL for various version of python](https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal)

### Installation with python/pip

_Note: It's recommended to use a virtual environment_

#### from [PyPI](https://pypi.org/project/modape/)

```
pip3 install modape
```

#### from [GitHub](https://github.com/WFP-VAM/modape.git)

```
git clone https://github.com/WFP-VAM/modape.git
cd modape
pip3 install .
```

or with `setup.py`:

```
git clone https://github.com/WFP-VAM/modape.git
cd modape
python3 setup.py install
```

### Installation with Docker

```
git clone https://github.com/WFP-VAM/modape.git
cd modape
docker build -t modape .

# check if it's working
docker run --rm -it modape modape_version
```

## Quick guide to naming convention

Since MODAPE was developed for WFP VAM's operational needs, it _heavily_ uses an established naming convention for different variables etc.


### VAM parameter codes

- **VIM**: MODIS NDVI
- **VEM**: MODIS EVI
- **LTD**:
    + **TDA**
    + **TDT**
- **LTN**:
    + **TNA**
    + **TNT**
