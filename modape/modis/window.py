"""
MODIS window functions & classes.

This file contains the classes and functions for extracting GeoTIFFs from HDF5 files.

Author: Valentin Pesendorfer, April 2019
"""
from __future__ import absolute_import, division, print_function

from contextlib import contextmanager
from datetime import datetime
from os.path import basename
try:
    from pathlib2 import Path
except ImportError:
    from pathlib import Path
import re


import numpy as np
try:
    import gdal
except ImportError:
    from osgeo import gdal
import h5py # pylint: disable=import-error

from modape.utils import date2label, dtype_GDNP, fromjulian, ldom

class ModisMosaic(object):
    """Class for mosaic of MODIS tiles.

    Moisaics tiles per Product, parameter and timestep. Enables extraction as GeoTiff.
    """

    def __init__(self, files, datemin, datemax, global_flag):
        """ Creates ModisMosaic object.

        Args:
            files: List of paths to files used for creating the mosaic
            datemin: Datestring for date of earliest mosaic (format YYYYMM)
            datemax: Datestring for date of latest mosaic (format YYYYMM)
            global_flag: Boolean flag if mosaic is global product
        """

        tile_re = re.compile(r'.+(h\d+v\d+).+') # Regular expression for tile ID
        self.global_flag = global_flag
        self.tiles = [re.sub(tile_re, '\\1', basename(x)) for x in files]
        self.tiles.sort()
        self.files = files

        # Extract tile IDs
        self.h_ix = list({re.sub(r'(h\d+)(v\d+)', '\\1', x) for x in self.tiles})
        self.h_ix.sort()
        self.v_ix = list({re.sub(r'(h\d+)(v\d+)', '\\2', x) for x in self.tiles})
        self.v_ix.sort()

        # get referece tile identifiers
        ref_tile_h = min([x for x in self.tiles if min(self.h_ix) in x])
        ref_tile_v = min([x for x in self.tiles if min(self.v_ix) in x])

        # vertical reference tile is top (left)
        ref = [x for x in self.files if ref_tile_v in x][0]

        # Read metadata from HDF5
        try:
            with h5py.File(ref, 'r') as h5f:
                dset = h5f.get('data')
                self.tile_rws = dset.attrs['RasterYSize'].item()
                self.tile_cls = dset.attrs['RasterXSize'].item()
                temporalresolution = dset.attrs['temporalresolution']
                self.datatype = dset.dtype
                gt_temp_v = dset.attrs['geotransform']
                self.prj = dset.attrs['projection']
                self.nodata = dset.attrs['nodata'].item()
                self.labels = []

                # If file is global, resolution is already in degrees, otherwhise it's resolution divided with 112000
                if self.global_flag:
                    self.resolution_degrees = dset.attrs['resolution']
                else:
                    self.resolution = dset.attrs['resolution']
                    self.resolution_degrees = self.resolution/112000
                dset = None
                self.dates = [x.decode() for x in h5f.get('dates')[...]]

                # if dates are either 5 or 10 days, retrive pentad or dekad labels
                if int(temporalresolution) == 5 or int(temporalresolution) == 10:
                    self.labels = date2label(self.dates, temporalresolution)

            # checking referece tile for h
            ref = [x for x in self.files if ref_tile_h in x][0]

            with h5py.File(ref, 'r') as h5f:
                dset = h5f.get('data')
                gt_temp_h = dset.attrs['geotransform']
                dset = None
            self.gt = [y for x in [gt_temp_h[0:3], gt_temp_v[3:6]] for y in x]
        except Exception as e:
            print('\nError reading referece file {} for mosaic! Error message: {}\n'.format(ref, e))
            raise

        # Create temporal index from dates available and min max input
        dates_dt = [fromjulian(x) for x in self.dates]

        if datemin:
            try:
                datemin_p = datetime.strptime(datemin, '%Y%m').date()
            except ValueError:
                try:
                    datemin_p = datetime.strptime(datemin, '%Y-%m-%d').date()
                except ValueError:
                    raise ValueError('Invalid begin time specified for mosaic. Accepted formats: {%Y%m, %Y-%m-%d}')

        else:
            datemin_p = dates_dt[0]

        if datemax:
            try:
                datemax_p = ldom(datetime.strptime(datemax, '%Y%m'))
            except ValueError:
                try:
                    datemax_p = datetime.strptime(datemax, '%Y-%m-%d').date()
                except ValueError:
                    raise ValueError('Invalid end time specified for mosaic. Accepted formats: {%Y%m, %Y-%m-%d}')
        else:
            datemax_p = dates_dt[-1]

        self.temp_index = np.flatnonzero(np.array([datemin_p <= x <= datemax_p for x in dates_dt]))

    def get_array(self, dataset, ix, dt):
        """Reads values for mosaic into array.

        Args:
            dataset: Defines dataset to be read from HDF5 file (default is 'data')
            ix: Temporal index
            dt: Datatype (default will be read from file)

        Returns
            Array for mosaic
        """

        # Initialize array
        tiles_array = np.full(((len(self.v_ix) * self.tile_rws), len(self.h_ix) * self.tile_cls), self.nodata, dtype=dt)

        # read data from intersecting HDF5 files
        for h5f in self.files:
            # Extract tile ID from filename
            t_h = re.sub(r'.+(h\d+)(v\d+).+', '\\1', basename(h5f))
            t_v = re.sub(r'.+(h\d+)(v\d+).+', '\\2', basename(h5f))

            # Caluclate row/column offset
            xoff = self.h_ix.index(t_h) * self.tile_cls
            yoff = self.v_ix.index(t_v) * self.tile_rws

            try:
                with h5py.File(h5f, 'r') as h5f_o:
                    # Dataset 'sgrid' is 2D, so no idex needed
                    if dataset == 'sgrid':
                        tiles_array[yoff:(yoff+self.tile_rws),
                                    xoff:(xoff+self.tile_cls)] = h5f_o.get(dataset)[...].reshape(self.tile_rws, self.tile_cls)
                    else:
                        tiles_array[yoff:(yoff+self.tile_rws),
                                    xoff:(xoff+self.tile_cls)] = h5f_o.get(dataset)[..., ix].reshape(self.tile_rws, self.tile_cls)
            except Exception as e:
                print('Error reading data from file {} to array! Error message {}:\n'.format(h5f, e))
                raise
        return tiles_array

    def get_array_global(self, dataset, ix, dt):
        """Reads values for global mosaic into array.

        Since files are global, the array will be a spatial and temporal subset rather than a mosaic.

        Args:
            dataset: Defines dataset to be read from HDF5 file (default is 'data')
            ix: Temporal index
            dt: Datatype (default will be read from file)

        Returns
            Array for mosaic
        """

        global_array = np.full((self.tile_rws, self.tile_cls), self.nodata, dtype=dt)

        for h5f in self.files:
            try:
                with h5py.File(h5f, 'r') as h5f_o:
                    if dataset == 'sgrid':
                        global_array[...] = h5f_o.get(dataset)[...].reshape(self.tile_rws, self.tile_cls)
                    else:
                        global_array[...] = h5f_o.get(dataset)[..., ix].reshape(self.tile_rws, self.tile_cls)
            except Exception as e:
                print('Error reading data from file {} to array! Error message {}:\n'.format(h5f, e))
                raise
        return global_array

    @contextmanager
    def get_raster(self, dataset, ix):
        """Generator for mosaic raster.

        This generator can be used within a context manager and will yield an in-memory raster.

        Args:
            dataset: Defines dataset to be read from HDF5 file (default is 'data')
            ix: Temporal index

        Yields:
            in-memory raster to be passed to GDAL warp
        """

        try:
            if dataset == 'sgrid':
                self.dt_gdal = dtype_GDNP('float32') # dtype for sgrid is set to float32
            else:
                self.dt_gdal = dtype_GDNP(self.datatype.name)
        except IndexError:
            print('\n\n Couldn\'t read data type from dataset. Using default Int16!\n')
            self.dt_gdal = (3, 'int16')

        # Use the corresponding getArray function if global_flag
        if self.global_flag:
            value_array = self.get_array_global(dataset, ix, self.dt_gdal[1])
        else:
            value_array = self.get_array(dataset, ix, self.dt_gdal[1])
        height, width = value_array.shape
        driver = gdal.GetDriverByName('GTiff')

        # Create in-memory dataset with virtual filename driver
        self.raster = driver.Create('/vsimem/inmem.tif', width, height, 1, self.dt_gdal[0])
        self.raster.SetGeoTransform(self.gt)
        self.raster.SetProjection(self.prj)
        rb = self.raster.GetRasterBand(1)
        rb.SetNoDataValue(self.nodata)

        # Write array
        rb.WriteArray(value_array)
        yield self

        # Cleanup to be exectuted when context manager closes after yield
        gdal.Unlink('/vsimem/inmem.tif')
        self.raster = None
        driver = None
        del value_array

def modis_tiles(aoi):
    """Function for querying MODIS tiles.

    Converts AOI coordinates to MODIS tile numbers by extracting values from MODIS_TILES.tif.

    Args:
        aoi: AOI coordinates, either LAT LON or XMIN, YMAX, XMAX, YMIN

    Returns:
        List of MODIS tile IDs intersecting AOI
    """

    # Load MODIS_TILES.tif from data directory
    this_dir = Path(__file__).parent
    ds = gdal.Open(this_dir.parent.joinpath('data', 'MODIS_TILES.tif').as_posix())

    # Try to catch TIFF issues
    try:
        gt = ds.GetGeoTransform()
    except AttributeError:
        raise SystemExit('Could not find \'MODIS_TILES.tif\' index raster. Try re-installing the package.')

    # Indices fpr point AOI
    if len(aoi) == 2:
        xo = int(round((aoi[1]-gt[0])/gt[1]))
        yo = int(round((gt[3]-aoi[0])/gt[1]))
        xd = 1
        yd = 1
    # Indices for bounding box AOI
    elif len(aoi) == 4:
        xo = int(round((aoi[0]-gt[0])/gt[1]))
        yo = int(round((gt[3]-aoi[1])/gt[1]))
        xd = int(round((aoi[2] - aoi[0])/gt[1]))
        yd = int(round((aoi[1] - aoi[3])/gt[1]))
    tile_extract = ds.ReadAsArray(xo, yo, xd, yd) # Read
    ds = None
    tile_tmp = np.unique(tile_extract/100) # Tile IDs are stored as H*100+V
    tiles = ['{:05.2f}'.format(x) for x in tile_tmp[tile_tmp != 0]]

    return ['h{}v{}'.format(*x.split('.')) for x in tiles]
