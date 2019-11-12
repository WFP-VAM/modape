# CHANGES

- ## v 0.4
  - #### v 0.4.0:
    - New Dockerfile based on Ubuntu bionic

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
