
The [MODIS processing chain](../#overview) can be executed at a high level using CLI scrips with minimal user input. All of the scripts are created using the Python module [click](https://click.palletsprojects.com/en/7.x/).

### `modis_download`

With `modis_download` you can query and download MODIS raw data, where the downloading to local disk requires valid [Earthdata credentials](https://urs.earthdata.nasa.gov/).

The products to be queried for need to be supplied as a list of MODIS product IDs (case insensitive), each item separated by a single space, e.g. `MOD13A2 MYD11A1 ...`.

If you want to query for both MODIS satellites for the same product, the letter indicating which satellite can be replaced by a `?`, e.g. `M?d13A2`.

The query can be limited by coordinates (as `--roi` point or bounding box), MODIS tile IDs (supplied as comma separated list to `--tile-filter`, e.g `h20v08,h20v09`) and a date bracket (supplied to `--begin-date` and `--end-date`). Due to the nature of composited products, results can be returned slightly outside the specified date bracket. To avoid this and strictly enforce the bracket, pass the `--strict-dates` flag.

The query is sent to the NASA Common Metadata Repository (CMR) using the [`python-cmr`](https://github.com/jddeal/python-cmr) module. To look at the details of the returned results, pass the `--return-results` flag.

The download is only performed when both credential flags are passed (`--username` & `--password`) as well as the `--download` flag. By default, the download is performed sequentially using Python's `requests`. To speed up the download, you may use multiple threads (`--multithread`) where the number of used threads can be adjusted with the `--nthreads` flag (default is 4 threads).

A processing chain might require the target directory to not contain any previous HDF files, which can be enforced using `--target-empty`. By default, existing files won't be overwritten. If this is required, pass the `--overwrite` flag.

#### Usage

```
Usage: modis_download [OPTIONS] [PRODUCTS]...

Options:
  --roi TEXT                   Region of interest. Either LAT,LON or
                               xmin,ymin,xmax,ymax

  -b, --begin-date [%Y-%m-%d]  Start date for query
  -e, --end-date [%Y-%m-%d]    End date for query
  -d, --targetdir PATH         Destination directory for downloaded files
  --target-empty               Fail if there are hdf files in the target
                               directory

  --tile-filter TEXT           Filter tiles - supplied as csv list
  --username TEXT              Earthdata username
  --password TEXT              Earthdata password
  --strict-dates               Don't allow files with timestamps outside of
                               provided date(s)

  --return-results             Print results to console
  --download                   Download data
  --overwrite                  Overwrite existing files
  --multithread                Use multiple threads for downloading
  --nthreads INTEGER           Number of threads to use
  -c, --collection TEXT        MODIS collection
  --help                       Show this message and exit.
```

### `modis_collect`

Before applying the Whittake smoother to the data, the raw HDF files need to be collected into a 3D spatio-temporal cube inside an HDF5 file, using `modis_collect`.

The only required input is `src_dir`, a directory with the raw MODIS HDF files to be collected. These can be multiple products and tiles.

If no target directory (`--target-dir`) is specified, the output will be stored in the `src_dir`. The HDF5 files themselves get stored in a sub-directory, which is attributed to the extracted parameter (indicated by the VAM parameter code). This needs to be considered carefully, as updating the existing files depends on the executable being able to locate them.

To extract a specific parameter, a [VAM parameter code](../#vam-parameter-codes) can be specified using `--vam-code`. If not specified, the defaults will be extracted, which are NDVI (VIM) and LST Daytime (TDA/TDT).

In the case of 16 day MODIS NDVI products, both satellites can be interleaved to form a combined 8-day product, indicated by a new _MXD_ product code.

Each combination of MODIS product and tile get collected into the same HDF5 file. Using `--parallel-tiles`, these can be collected in parallel (up to 4 collection processes). If an HDF5 file does not exist, it will be created, otherwise updated with the new data.

It might be required to perform a check on temporal continuity, which can be done using `--last-collected`. Here a MODIS julian date can be specified that needs to be the last timestep collected before the new collection starts. If that's not the case, the process will fail on exception.

_Note: It's the user's responsibility to ensure temporal continuity when creating and updating files. New timesteps are appended, and there's no internal checks if a timestep might be missing etc. This follows the "garbage in - garbage out" principle._

#### Usage

```
Usage: modis_collect [OPTIONS] SRC_DIR

Options:
  -d, --targetdir PATH      Destination for raw HDF5 files
  -x, --compression TEXT    Compression for HDF5 files
  --vam-code TEXT           VAM code for dataset to process
  --interleave              Interleave MOD13 & MYD13 products to MXD (only
                            works for VIM!)

  --parallel-tiles INTEGER  Number of tiles processed in parallel (default =
                            None)

  --cleanup                 Remove collected HDF files
  --last-collected [%Y%j]   Last collected date in julian format (YYYYDDD -
                            %Y%j)

  --help                    Show this message and exit.

```

### `modis_smooth`

This step applies the Whittaker smoother to the raw data. The input (`src`) can be either a single raw HDF5 file, or a directory containing multiple.

There's a one-to-many relationship between raw and smooth HDF5 files. One raw file can have multiple smooth files, depending on the settings used for smoothing. As with the raw HDF5 files, the smooth files get created if they don't exist or updated if they do. So it's again crucial to be aware of where the data gets written to, and adjust accordingly with `--targetdir`.

The smoothing can be performed in multiple ways, which is determined by the user's input:

- if the `--soptimize` flag is passed, the V-curve optimization of the smoothing parameter S is performed. This always takes precedence over other methods, independent of other inputs supplied.
- if `--soptimize` is passed alongside a P value (`--pvalue`), then an asymmetric V-curve optimization is performed using the supplied P value.
- if only the S Value (`--svalue`) is supplied, then the smoothing is performed using a fixed S (as supplied - should be in log10)
- if none of the above is supplied, the smoothing will be performed reading S values from a previously optimized S-grid. **Note:** This requires that `--soptimize` has been run at least once previously.

The V-curve optimization will look for an optimized S value for each pixel within a given range. This range can be defined using the `--srange` flag (expects 3 space separated log10 values). If the range is not supplied, the algorithm will use one of two for each pixel, depending on the lag-1 correlation of the pixel's timeseries:

- lag-1 correlation > 0.5: `smin=-2`, `smax=1`, `sstep=0.2`
- lag-1 correlation <= 0.5: `smin=0`, `smax=3`, `sstep=0.2`

The smoother can also perform temporal interpolation, to a desired temporal grid. The desired timestep should be defined with `--tempint`. Dekad (10-day) and pentad (5-day) interpolation follow a fixed date cycle. For other more custom interpolations, it's required to specify a start date (`--tempint-start`) from which the timestep will be calculated.

For subsequent runs after optimization, the number of timesteps used for smoothing (`--nsmooth`) and the number updated in the target dataset (`--nupdate`) can be limited to speed up the process.

Similar to raw HDF5 files, the `--last-collected` flag can be passed to enforce a check on the last raw date used for smoothing in the previous run.


#### Usage

```
Usage: modis_smooth [OPTIONS] SRC

Options:
  -d, --targetdir PATH            Target directory for smoothed output
  -s, --svalue FLOAT              S value for smoothing (in log10)
  -S, --srange FLOAT...           S value range for V-curve (float log10(s)
                                  values as smin smax sstep - default -1 1 0)

  -p, --pvalue FLOAT              P value for asymmetric smoothing
  -t, --tempint INTEGER           Value for temporal interpolation (integer
                                  required - default is native temporal
                                  resolution i.e. no interpolation)

  --tempint-start [%Y-%m-%d|%Y%j]
                                  Startdate for temporal interpolation
  -n, --nsmooth INTEGER           Number of raw timesteps used for smoothing
  -u, --nupdate INTEGER           Number of smoothed timesteps to be updated
                                  in HDF5 file

  --soptimize                     Use V-curve for s value optimization
  --parallel-tiles INTEGER        Number of tiles processed in parallel
  --last-collected [%Y%j]         Last collected date in julian format
                                  (YYYYDDD - %Y%j)

  --help                          Show this message and exit.
```

### `modis_window`

`modis_window` enables the export of data from HDF5 files as GeoTIFF. The data can be a mosaic of many tiles, or a subset of a single tile / global file.

The input (`src`) can be a single HDF5 file or a directory containing multiple files. In the second case, the input can be filtered by MODIS product (`--filter-product`) and VAM parameter code (`--filter-vampc`).

For subsets, a region of interest (`--roi`) needs to be specified as comma separated list. The exported dates can be limited by passing a date bracket. (`--begin-date` & `--end-date`).

In case of smooth HDF5 files, the grid containing the V-curve optimized S values for each pixel can be exported by passing the `--sgrid` flag. In this case any input dates won't have an effect (there's only one S-grid).

By default, all outputs are warped to `EPSG:4326`. For warping to a different spatial reference, specify a `GDAL` readable reference to `--target-srs`.

By default, all exported GeoTIFFs get compressed using `LZW` with predictor 2 (using the creation options `["COMPRESS=LZW", "PREDICTOR=2"]`). To specify different `GDAL` creation options, a space separated list in form of `KEY=VALUE` can be passed to `--co`. For further customization of the exported GeoTIFFs, additional GDAL options can be supplied to `--gdal-kwarg` in the same format.

To further improve the compression, the data can be rounded to a significant integer. E.g. a NDVI value of `3141` can be rounded to `3100` by passing `2` to `--round-int`.

The file naming for exported GeoTIFFs is detailed in the [main documentation](../#exported-geotiffs). A custom region code can be specified using `--region`. When exporting dekad / pentad data, `--force-doy` can be used to force a `YYYYjDDD` timestamp instead of the dekad / pentad format.
<br>_Note: please read the section in the main documentation carefully to avoid confusion with indentical named GeoTIFFs._

By default existing GeoTIFFs won't be overwritten. This can be forced with passwing the `--overwrite` flag.

#### Usage

```
Usage: modis_window [OPTIONS] SRC

Options:
  -d, --targetdir PATH         Target directory for Tiffs
  -b, --begin-date [%Y-%m-%d]  Begin date for Tiffs
  -e, --end-date [%Y-%m-%d]    End date for Tiffs
  --roi TEXT                   AOI for clipping (as ULX,ULY,LRX,LRY)
  --region TEXT                Region prefix for Tiffs (default is reg)
  --sgrid                      Extract sgrid instead of data
  --force-doy                  Force DOY filenaming
  --filter-product TEXT        Filter by product
  --filter-vampc TEXT          Filter by VAM parameter code
  --target-srs TEXT            Target spatial reference for warping
  --co TEXT                    GDAL creationOptions
  --clip-valid                 clip values to valid range for product
  --round-int INTEGER          Round to integer places (either decimals or
                               exponent of 10)

  --gdal-kwarg TEXT            Addition kwargs for GDAL in form KEY=VALUE
                               (multiple allowed)

  --overwrite                  Overwrite existsing Tiffs
  --help                       Show this message and exit.
```
