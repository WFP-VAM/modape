#!/usr/bin/env python
from __future__ import print_function
from .whittaker import *
import glob
import re
import os
import sys
import gdal
import argparse
import time
import numpy as np
import array


def initGDAL(x,p,fn=None):

    if not fn:
        fn = os.path.basename(x)

    ds = gdal.Open(x)
    dr = ds.GetDriver()
    ds_new = dr.Create(p+fn,ds.RasterXSize,ds.RasterYSize,ds.RasterCount,ds.GetRasterBand(1).DataType)
    ds_new.SetGeoTransform(ds.GetGeoTransform())
    ds_new.SetProjection(ds.GetProjection())
    ds = None
    dr = None
    ds_new = None


def iterateBlocks(rows,cols,n):
    for ri in range(0,rows,n):
        for ci in range(0,cols,n):
            yield (ri,min(n,rows-ri),ci,min(n,cols-ci))


class RTS:

    def __init__(self,files,targetdir,bsize=256,nodata=0):

        self.files = files
        self.files.sort()
        self.ref_file = self.files[0]
        self.nfiles = len(self.files)
        self.bsize = bsize
        self.nodata = nodata

        ds = gdal.Open(self.ref_file)

        self.nrows = ds.RasterYSize
        self.ncols = ds.RasterXSize
        ds = None

        #self.parameters_ = {}

        self.targetdir = targetdir


        self.outfiles = [self.targtedir + '/' + os.path.basename(x) for x in self.files]

    def initRasters(self):

        for f in self.files:
            try:
                initGDAL(f,self.targetdir)
            except AttributeError:
                print("Error initializing {}! Please check data".format(f))
                raise


    def ws2d(self,s):

        for yo, ys, xo, xs in iterateBlocks(self.nrows,self.ncols,self.bsize):

            arr = np.zeros((ys*xs,self.nfiles),dtype='float32')
            wts = arr.copy()

            arr_helper = arr.view()
            arr_helper.shape = (ys,xs,self.nfiles)

            for fix in range(self.nfiles):

                ds = gdal.Open(self.files[fix])

                arr_helper[...,fix] = ds.ReadAsArray(xoff=xo,xsize=xs,yoff=yo,ysize=ys)

                ds = None

            wts[...] = (arr != self.nodata) * 1

            for r in range(arr_helper.shape[0]):
                if wts[r,...].sum().item() != 0.0:
                    arr[r,...] = ws2d(arr[r,...],10**s,wts[r,...])


            for fix in range(self.nfiles):

                ds = gdal.Open(self.outfiles[fix],gdal.GA_Update)

                ds_b = ds.GetRasterBand(1)

                ds_b.WriteArray(arr_helper[...,fix].round(),xo,yo)
                ds_b.FlushCache()

                ds_b = None
                ds = None

            with open(self.targetdir + '/filt0_config.txt') as thefile:

                thefile.write('Running whittaker smoother with fixed s value\n')
                thefile.write('\n')
                thefile.write('Start: {}\n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
                thefile.write('\n')
                thefile.write('Sopt: {}\n'.format(10**s))
                thefile.write('log10(Sopt): {}\n'.format(s))








        # create array

        # read and broadcast

        # write with write Band


        #loop key in dict

        # write smoothing method, timesteps and params to txt file











def main():

    parser = argparse.ArgumentParser(description="Extract a window from MODIS products")
    parser.add_argument("path", help='Path to processed MODIS h5 files')
    parser.add_argument("-P","--pattern", help='Pattern to filter file names',default = '*' ,metavar='')
    parser.add_argument("-d","--targetdir", help='Target directory for GeoTIFFs (default current directory)',default=os.getcwd(),metavar='')
    parser.add_argument("-s","--svalue", help='S value for smoothing (has to be log10(s)', metavar='', type = float)
    parser.add_argument("-S","--srange", help='S range for V-curve (float log10(s) values as smin smax sstep - default 0.0 4.0 0.1)',nargs='+',metavar='', type=float)
    parser.add_argument("-p","--pvalue", help='Value for asymmetric smoothing (float required)', metavar='', type = float)
    parser.add_argument("--soptimize", help='Use V-curve for s value optimization',action='store_true')
    #parser.add_argument("--sgrid", help='Extract (mosaic of) s value grid(s))',action='store_true')

    # fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    print('\n[{}]: Starting smoothRTS.py ... \n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    if not os.path.exists(args.path):
        raise SystemExit('directory PATH does not exist!')


    fls = glob.glob('{}/{}'.format(args.path,args.pattern))

    if not len(fls) > 0:
        raise SystemExit('No files found in {} with pattern {}, please check input.'.format(args.path,args.pattern))

    rts = RTS(files=fls,targetdir = args.targetdir)

    if args.soptimize:

        if args.srange:
            try:
                assert len(args.srange) == 3
                srange = array.array('f',np.linspace(float(args.srange[0]),float(args.srange[1]),float(args.srange[1])/float(args.srange[2]) + 1.0))
            except (IndexError,TypeError,AssertionError):
                raise SystemExit('Error with s value array values. Expected three values of float log10(s) -  smin smax sstep !')
        else:
            srange = array.array('f',np.linspace(0.0,4.0,41.0))

        if args.pvalue:

            print('\nRunning asymmetric whittaker smoother with v-curve optimization ... \n')

            try:
                os.makedirs(args.targetdir + '/filtvcp/')
            except:
                print('Issues creating subdirectory in {}'.fromat(args.path))
                raise

            #ws2d_vc_asy(srange=srange,p=p)

        else:

            print('\nRunning whittaker smoother with v-curve optimization ... \n')

            try:
                os.makedirs(args.targetdir + '/filtvc/')
            except:
                print('Issues creating subdirectory in {}'.fromat(args.path))
                raise


            #ws2d_vc(srange=srange)

    else:

        ## insert if-clause for s value if grid option is needed (see smoothMODIS)

        try:
            s = 10**float(args.svalue)
        except:
            raise SystemExit('Error with s value. Expected float log10(s)!')

        print('\nRunning whittaker smoother with fixed s value ... \n')


        try:
            os.makedirs(args.targetdir + '/filt0/')
        except:
            print('Issues creating subdirectory in {}'.fromat(args.path))
            raise

        #ws2d(s=s)


    print('\n[{}]: smoothMODIS.py finished successfully.\n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))





if __name__=='__main__':
    main()
