# CHANGES

- ## v 0.2
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