
Additional commandline scripts which are not strictly part of the MODIS processing chain.

### `csv_smooth`

This executable allows to smooth timeseries saved in CSV format.

!!! Note
    For the different ways the smoother can be applied, please refer to the [`modis_smooth`](../modis_executables/#modis_smooth) documentation.

`csv_smooth` expects the timeseries for smoothing to be stored in the columns, therefore each column will be smoothed independently. If a column with dates or another index is present, specify it in `--index-column`.

Header rows can be skipped using the `--skip-header-rows` flag. If there are column names or other metainfo (e.g. coordinates) that should not be filtered, the index of the first data cell can be specified using `--data-start`.

#### Usage

```
Usage: csv_smooth [OPTIONS] CSV_FILE

Options:
  -s, --svalue FLOAT          S value for smoothing (has to be log10(s))
  -S, --srange FLOAT...       S value range for V-curve (float log10(s) values
                              as smin smax sstep - default 0 4 0.1)

  -p, --pvalue FLOAT          P value for asymmetric smoothing
  -n, --nodata FLOAT          nodata value
  --skip-header-rows INTEGER  Number of header rows to skip
  --data-start INTEGER        Number of row with first datapoint
  --index-column INTEGER      Index column number
  --soptimize                 Use V-curve for s value optimization
  --help                      Show this message and exit.
```


### `modis_info`

This executable can be used to inspect the raw/smooth HDF5 files. When supplied with the path to a file, a selection of metadata is printed to console, such as dimensions, begin/end date, last smoothing run (for smooth files), ... etc.

#### Usage

```
Usage: modis_info [OPTIONS] FILE

  Info tool for processed MODIS HDF5 files.

  Returns metadata on processed MODIS HDF5 files, both for raw and smoothed
  files.

Options:
  --help  Show this message and exit.
```
