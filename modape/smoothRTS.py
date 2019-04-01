#!/usr/bin/env python
# pylint: disable=line-too-long, too-many-statements, wildcard-import, C0103, R0205, E0602

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import array
import glob
import os
import sys
import time

import numpy as np
try:
    import gdal
except ImportError:
    from osgeo import gdal

from modape.utils import dtype_GDNP
from modape.whittaker import *

def init_gdal(x, p, fn=None, dt=None):
    '''Initializes empty GeoTIFF based on template.

    Args:
        x (str): Path to template
        p (str): Target directory
        fn (str): Output filename (optional)
        dt (str): Output datatype (optional)
    '''

    if not fn:
        fn = os.path.basename(x) # same filename if none is supplied
    ds = gdal.Open(x)
    dr = ds.GetDriver()
    if not dt:
        dt_new = ds.GetRasterBand(1).DataType # same datatype if none is supplied
    else:
        dt_new = dtype_GDNP(dt)[0] # Parse datatype

    # Create empty copy
    ds_new = dr.Create(p + fn, ds.RasterXSize, ds.RasterYSize, ds.RasterCount, dt_new)
    ds_new.SetGeoTransform(ds.GetGeoTransform())
    ds_new.SetProjection(ds.GetProjection())
    ds_new.GetRasterBand(1).SetNoDataValue(ds.GetRasterBand(1).GetNoDataValue())
    ds = None
    dr = None
    ds_new = None

def iterateBlocks(rows, cols, n):
    '''Generator for blockwise iteration over array.

    Args:
        rows (int): Number of rows
        cols (int): Number of columns
        n (int): Side length of block

    Yields:
        Tuple with values
            - Start row
            - Number of rows
            - Start column
            - Number of columns
    '''

    # Iterate over rows then columns
    for ri in range(0, rows, n):
        for ci in range(0, cols, n):
            yield (ri, min(n, rows-ri), ci, min(n, cols-ci))


class RTS(object):
    '''Class for raster timeseries for smoothing.'''

    def __init__(self, files, targetdir, bsize=256, nodata=None):
        '''Creates instance of raster timeseries class.

        The metadata for the timeseries is extracted from the first file in
        the directoryself.

        Args:
            files ([str]): List or filepaths to process
            targetdir (str): Target directory for smoothed files
            bsize (int): Side length of processing blocks (default = 256)
            nodata (int): Nodata value (default is read from reference file)
        '''

        # Select reference file and sort
        self.files = files
        self.files.sort()
        self.ref_file = self.files[0]
        self.nfiles = len(self.files)
        self.bsize = bsize
        ds = gdal.Open(self.ref_file)
        self.nrows = ds.RasterYSize
        self.ncols = ds.RasterXSize

        if nodata:
            self.nodata = nodata
        else:
            self.nodata = ds.GetRasterBand(1).GetNoDataValue() # nodata from file
            if not self.nodata:
                self.nodata = 0 # Set to 0 if read fails
                print('Failed to read NoData value from files. NoData set to 0.')
        ds = None
        self.targetdir = targetdir

    def initRasters(self, tdir):
        '''Intitialize empty rasters for smoothed data.

        Args:
            tdir (str): Target directory
        '''

        # Iterate over files
        for f in self.files:
            try:
                init_gdal(f, tdir) # Initialize empty copy
            except AttributeError:
                print('Error initializing {}! Please check data'.format(f))
                raise

    def ws2d(self, s):
        '''Apply whittaker smoother with fixed s-value to data.

        Args:
            s (float): log10 value of s
        '''

        tdir = self.targetdir + '/filt0/'

        # Create full path filenames for smoothed rasters
        outfiles = [tdir + '/' + os.path.basename(x) for x in self.files]
        if not os.path.exists(tdir):
            try:
                os.makedirs(tdir)
            except:
                print('Issues creating subdirectory {}'.format(tdir))
                raise
        self.initRasters(tdir) # Initialize rasters

        # Iterate over blocks
        for yo, ys, xo, xs in iterateBlocks(self.nrows, self.ncols, self.bsize):
            arr = np.zeros((ys*xs, self.nfiles), dtype='double') # values array
            wts = arr.copy() # weights array
            arr_helper = arr.view() # helper
            arr_helper.shape = (ys, xs, self.nfiles)

            # Iterate files to read data
            for fix in range(self.nfiles):
                ds = gdal.Open(self.files[fix])
                arr_helper[..., fix] = ds.ReadAsArray(xoff=xo, xsize=xs, yoff=yo, ysize=ys)
                ds = None

            # Data which is not nodata gets weight 1, others 0
            wts[...] = (arr != self.nodata) * 1
            ndix = np.sum(arr != self.nodata, 1) > 0 #70
            mapIX = np.where(ndix)[0]

            if mapIX.size == 0:
                continue # skip bc no data in block

            arr[np.logical_not(ndix), :] = self.nodata

            for r in mapIX:
                arr[r, ...] = ws2d(arr[r, ...], 10**s, wts[r, ...])

            # Write smoothed data to disk
            for fix in range(self.nfiles):
                ds = gdal.Open(outfiles[fix], gdal.GA_Update)
                ds_b = ds.GetRasterBand(1)
                ds_b.WriteArray(arr_helper[..., fix].round(), xo, yo)
                ds_b.FlushCache()

                ds_b = None
                ds = None

        # Write config text file to disk with processing parameters and info
        with open(tdir + 'filt0_config.txt', 'w') as thefile:
            thefile.write('Running whittaker smoother with fixed s value\n')
            thefile.write('\n')
            thefile.write('Timestamp: {}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
            thefile.write('\n')
            thefile.write('Sopt: {}\n'.format(s))
            thefile.write('log10(Sopt): {}\n'.format(np.log10(s)))
            thefile.write('Nodata value: {}\n'.format(self.nodata))
            thefile.write('\n')

    def ws2dopt(self, srange, p=None):
        '''Apply whittaker smoother with V-curve optimization of s to data.

        If a p-value is supplied, the asymmetric whittaker smoother will be
        applied.

        Args:
            srange (arr): Float32 array of s-values to apply
            p (float): P-value for percentile
        '''

        srange_arr = array.array('d', srange)
        if p:
            tdir = self.targetdir + '/filtoptvp/'
        else:
            tdir = self.targetdir + '/filtoptv/'
        outfiles = [tdir + '/' + os.path.basename(x) for x in self.files]

        if not os.path.exists(tdir):
            try:
                os.makedirs(tdir)
            except:
                print('Issues creating subdirectory {}'.format(tdir))
                raise

        self.sgrid = tdir + 'sgrid.tif' # Path to s-grid
        self.initRasters(tdir)

        # S-grid needs to be initialized separately
        init_gdal(self.ref_file, tdir, 'sgrid.tif', dt='float32')

        for yo, ys, xo, xs in iterateBlocks(self.nrows, self.ncols, self.bsize):
            arr = np.zeros((ys*xs, self.nfiles), dtype='double')
            wts = arr.copy()
            sarr = np.zeros((ys*xs), dtype='double')
            arr_helper = arr.view()
            arr_helper.shape = (ys, xs, self.nfiles)

            for fix in range(self.nfiles):
                ds = gdal.Open(self.files[fix])
                arr_helper[..., fix] = ds.ReadAsArray(xoff=xo, xsize=xs, yoff=yo, ysize=ys)
                ds = None
            wts[...] = (arr != self.nodata)*1

            ndix = np.sum(arr != self.nodata, 1) > 0 #70
            mapIX = np.where(ndix)[0]

            if mapIX.size == 0:
                continue # skip bc no data in block

            arr[np.logical_not(ndix), :] = self.nodata

            for r in mapIX:
                if p:
                    arr[r, ...], sarr[r] = ws2doptvp(arr[r, ...], wts[r, ...], srange_arr, p)
                else:
                    arr[r, ...], sarr[r] = ws2doptv(arr[r, ...], wts[r, ...], srange_arr)

            for fix in range(self.nfiles):
                ds = gdal.Open(outfiles[fix], gdal.GA_Update)
                ds_b = ds.GetRasterBand(1)
                ds_b.WriteArray(arr_helper[..., fix].round(), xo, yo)
                ds_b.FlushCache()
                ds_b = None
                ds = None

            # Convert s values in grid to log10(s)
            sarr[sarr > 0] = np.log10(sarr[sarr > 0])

            # Write s-values to grid
            ds = gdal.Open(self.sgrid, gdal.GA_Update)
            ds_b = ds.GetRasterBand(1)
            ds_b.WriteArray(sarr.reshape(ys, xs), xo, yo)
            ds_b.FlushCache()
            ds_b = None
            ds = None

            if p:
                with open(tdir + '/filtoptvp_config.txt', 'w') as thefile:
                    thefile.write('Running asymmetric whittaker smoother with V-curve optimization\n')
                    thefile.write('\n')
                    thefile.write('Timestamp: {}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
                    thefile.write('\n')
                    thefile.write('Sgrid: {}\n'.format(self.sgrid))
                    thefile.write('P value: {}\n'.format(p))
                    thefile.write('Nodata value: {}\n'.format(self.nodata))
                    thefile.write('\n')
            else:
                with open(tdir + '/filtoptv_config.txt', 'w') as thefile:
                    thefile.write('Running whittaker smoother with V-curve optimization\n')
                    thefile.write('\n')
                    thefile.write('Timestamp: {}\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
                    thefile.write('\n')
                    thefile.write('Sgrid: {}\n'.format(self.sgrid))
                    thefile.write('Nodata value: {}\n'.format(self.nodata))
                    thefile.write('\n')

def main():
    '''Apply whittaker smoother to a timeseries of local raster files.

    Raster files in path are combined to a timeseries and smoothed using the whittaker smoother,
    optionally with V-curve optimization of s.

    The user needs to make sure the raster files to be smoothed are identical in dimensions and type.

    Parallel processing is currently not implemented, so big timeseries might take some time!
    '''

    parser = argparse.ArgumentParser(description='Extract a window from MODIS products')
    parser.add_argument('path', help='Path containing raster files')
    parser.add_argument('-P', '--pattern', help='Pattern to filter file names', default='*', metavar='')
    parser.add_argument('-d', '--targetdir', help='Target directory for GeoTIFFs (default current directory)', default=os.getcwd(), metavar='')
    parser.add_argument('-s', '--svalue', help='S value for smoothing (has to be log10(s)', metavar='', type=float)
    parser.add_argument('-S', '--srange', help='S range for V-curve (float log10(s) values as smin smax sstep - default 0.0 4.0 0.1)', nargs='+', metavar='', type=float)
    parser.add_argument('-p', '--pvalue', help='Value for asymmetric smoothing (float required)', metavar='', type=float)
    parser.add_argument('-b', '--blocksize', help='Processing block side length (default 256)', default=256, metavar='', type=int)
    parser.add_argument('--nodata', help='NoData value', metavar='', type=float)
    parser.add_argument('--soptimize', help='Use V-curve (with p if supplied) for s value optimization', action='store_true')

    # fail and print help if no arguments supplied
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    print('\n[{}]: Starting smoothRTS.py ... \n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
    if not os.path.exists(args.path):
        raise SystemExit('directory PATH does not exist!')

    # Find files in path
    fls = [x for x in glob.glob('{}/{}'.format(args.path, args.pattern)) if os.path.isfile(x)]

    if not fls:
        raise ValueError('No files found in {} with pattern {}, please check input.'.format(args.path, args.pattern))

    # Create raster timeseries object
    rts = RTS(files=fls,
              targetdir=args.targetdir,
              bsize=args.blocksize,
              nodata=args.nodata)

    # V-curve optimization is triggered by either supplying the soptimize flag or a s-range
    if args.soptimize:
        # Parse s-range or use default
        if args.srange:
            try:
                assert len(args.srange) == 3
                srange = np.linspace(float(args.srange[0]),
                                     float(args.srange[1]),
                                     abs((args.srange[0]-args.srange[1]))/float(args.srange[2]) + 1.0)
            except (IndexError, TypeError, AssertionError):
                raise SystemExit('Error with s value array values. Expected three values of float log10(s) -  smin smax sstep !')
        else:
            srange = np.linspace(0.0, 4.0, 41.0)
        if args.pvalue:
            print('\nRunning asymmetric whittaker smoother with V-curve optimization ... \n')
            rts.ws2dopt(srange=srange, p=args.pvalue)
        else:
            print('\nRunning whittaker smoother with V-curve optimization ... \n')
            rts.ws2dopt(srange=srange)
    else:
        ## insert if-clause for s value if grid option is needed (see smoothMODIS)
        try:
            s = 10**float(args.svalue) # Convert s value from log10(s)
        except:
            raise SystemExit('Error with s value. Expected float log10(s)!')

        print('\nRunning whittaker smoother with fixed s value ... \n')
        rts.ws2d(s=s)

    print('\n[{}]: smoothMODIS.py finished successfully.\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))

if __name__ == '__main__':
    main()
