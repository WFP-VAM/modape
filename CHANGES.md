# CHANGES

- ## v 1.0
  - #### v 1.0.0
    - New Dockerfile (fix #82)
    - Improved downloading based on python requests #93
      - Better catching of failed downloads #29
      - new executable script
      - add option for strict date handling #87
      - allows downloading all available data from NASA CMR #86
      - adds option to fail if target directory is not empty #90
    - Reworked collection into MODIS raw HDF5 (#98)
     - new executable script
     - new `io` module with HDF5Base class
     - reworked `ModisRawH5` class with inheritance from HDF5Base
     - enable cleanup of collected HDF files with tracefile (#91)
     - enable check on last_collected
    - Reworked smoothing of raw HDF5 files (#100)
      - new executable script
      - reworked `ModisSmoothH5` class with inheritance from HDF5Base and single smoothing method which covers all Whittaker options
      - enable check on last_collected
      - fail if smoothing from non-initialized sgrid is requested (#80)
    - Reworked mosaicing to GeoTiff from HDF5 (#104)
      - new executable script
      - reworked `ModisMosaic` leveraging `gdal` VRT and in-memory rasters for performing the mosaicing
      - Warping can now be performed to user defined target spatial reference instead of just EPSG:4326
      - Optional clipping is now performed after warp, generating coherent results independent of the medthod (#85)
      - Improved control over gdal's `creationOptions`, including ability to pass kwargs directly to `gdal.Translate` (#89)
      - Optional clipping to valid data range or MODIS NDVI and LST (#88)
      - Optional rounding of integers to exponents of 10 (#88)

- ## v 0.3
  - #### v 0.3.0:
    - Split `modis.py` into separate sub-modules
    - MODIS download only possible using `aria2`
    - `modis_collect` now handles duplicated `hdf` files by taking file with most recent processing timestamp (issue #68)
    - When using `--interleave` in `modis_collect`, acquisitions before `2002185` (start of AQUA) are ignored and not collected (issue #65)
    - Fix calculation of asymmetric weights in `ws2optvp`
    - New smoothing function `ws2dp` runs the Whittaker filter with asymmetric weights and fixed S
    - Filtering with S values from grid now runs with the asymmetric filter `ws2dp`
    - minor enhancements

  - #### v 0.3.1:
    - fix `numpy` requirement in `setup.py`
  - #### v 0.3.2:
      - fix issues with `np.linspace` and `srange` (issue #75)
      - fix deprecation warning in `modis_download` (issue #76)
      - fix bug when windowing sgrid without AOI (issue #77)

- ## v 0.2
  - #### v 0.2.1:
    - Re-designing update of datasets (fix issue #66)
    - Remove incomplete downloads with aria2 in case of fail
    - Return filenames in `ModisQuery` only when files are on disk
    - Make overwrite of tiffs in `modis_window` optional with a flag
  - #### v 0.2.0:
    - greater changes and fixes to updating smooth datasets (fix issue #58)
    - change from `os` to `pathlib` for most path operations
    - `modis_window` to process all available products if not filtered (issue #63)
    - enable pentadal & dekadal labels for file naming in `modis_window` (issue #59)
    - address issues #57 & #61
    - simplify imports


- ## v 0.1
  - #### v 0.1.7:
    - fix bug when updating datasets
  - #### v 0.1.6:
    - fix bug handling vector AOIs in `modis_download`
    - add docstring to `modis_download`

  - #### v0.1.1/2/3/4/5:
    - changes to MANIFEST.in
    - fix issues with pytest and dates in HDF5 for PYTHON 2.7
    - fix bug for updating smoothed datasets
    - patch for handling duplicate raw files
    - patch for srange in sequential processing mode, minor indentation fixes

  - #### v0.1.0:
    - initial release
