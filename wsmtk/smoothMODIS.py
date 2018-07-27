#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import argparse
import multiprocessing as mp
import array
import numpy as np
from .modis import MODISsmth5
import time

def main():

    parser = argparse.ArgumentParser(description="Smooth, gapfill and interpolate processed raw MODIS HDF5 files")
    parser.add_argument("rawfile", help='Raw MODIS HDF5 file',metavar='RAW HDF5')
    parser.add_argument("-s","--svalue", help='S value for smoothing (has to be log10(s)', metavar='', type = float)
    parser.add_argument("-S","--srange", help='S value range for V-curve (float log10(s) values as smin smax sstep - default 0.0 4.0 0.1)',nargs='+',metavar='')
    parser.add_argument("-t","--tempint", help='Value for temporal interpolation (integer required - default is native temporal resolution AKA no interpolation)', metavar='',type = int)
    parser.add_argument("-n","--nsmooth", help='Number of raw timesteps used for smoothing', metavar='',type = int)
    parser.add_argument("-u","--nupdate", help='Number of smoothed timesteps to be updated in HDF5 file', metavar='',type = int)
    parser.add_argument("-p","--pvalue", help='Value for asymmetric smoothing (float required)', metavar='', type = float)
    parser.add_argument("-d","--targetdir", help='Target directory for smoothed output',default=os.getcwd(),metavar='')
    parser.add_argument("--soptimize", help='Use V-curve for s value optimization',action='store_true')
    parser.add_argument("--parallel", help='Parallel processing',action='store_true')
    parser.add_argument("--nworkers", help='Number of cores used for parallel processing (default is number available minus 1)',default=mp.cpu_count()-1, metavar='', type = int)


    # fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        print('\n')

    args = parser.parse_args()

    print('\n[{}]: Starting smoothMODIS.py ... \n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # check if raw file exisis

    if not os.path.isfile(args.rawfile):
        raise SystemExit('Raw HDF5 {} not found! Please check path.'.format(args.rawfile))

    print('\nInput file: {}\n'.format(args.rawfile))

    smt_h5 = MODISsmth5(rawfile = args.rawfile, tempint = args.tempint, nsmooth = args.nsmooth, nupdate = args.nupdate, targetdir = args.targetdir, parallel = args.parallel, nworkers = args.nworkers)

    if not smt_h5.exists:
        smt_h5.create()

    # check if v-curve optimization is true

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

            try:
                p = float(args.pvalue)
            except ValueError:
                raise SystemExit('Error using p-value ... expected float, got {}.'.format(type(args.pvalue)))

            print('\nRunning asymmetric whittaker smoother with v-curve optimization ... \n')

            smt_h5.ws2d_vc_asy(srange=srange,p=p)

        else:

            print('\nRunning whittaker smoother with v-curve optimization ... \n')

            smt_h5.ws2d_vc(srange=srange)

    else:

        if args.svalue:

            try:
                s = 10**float(args.svalue)
            except:
                raise SystemExit('Error with s value. Expected float log10(s)!')

            print('\nRunning whittaker smoother with fixed s value ... \n')

            smt_h5.ws2d(s=s)

        else:

            print('\nRunning whittaker smoother with s value from grid ... \n')

            smt_h5.ws2d_sgrid()

    print('\n[{}]: smoothMODIS.py finished successfully.\n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

if __name__ == '__main__':
    mp.freeze_support()
    main()
