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
    parser.add_argument("-l","--Lambda", help='LAMBDA value for smoothing (has to be log10(l)', metavar='')
    parser.add_argument("-L","--llas", help='LAMBDA range for V-curve (tuple of float log10(l) values where (lmin, lmax, lstep) - default (0,4,0.1))', metavar='')
    parser.add_argument("-t","--tempint", help='Value for temporal interpolation (integer required - default is native temporal resolution AKA no interpolation)', metavar='')
    parser.add_argument("-p","--pvalue", help='Value for asymmetric smoothing (float required)', metavar='')
    parser.add_argument("-d","--targetdir", help='Target directory for smoothed output',default=os.getcwd(),metavar='')
    parser.add_argument("--loptimize", help='Use V-curve for lambda optimization',action='store_true')
    parser.add_argument("--parallel", help='Parallel processing',action='store_true')
    parser.add_argument("--ncores", help='Number of cores used for parallel processing (defaul is number available minus 1)',default=mp.cpu_count()-1, metavar='')


    # fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)

    print('\n[{}]: Starting smoothMODIS.py ... \n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    args = parser.parse_args()

    # check if raw file exisis

    if not os.path.isfile(args.rawfile):
        raise SystemExit('Raw HDF5 {} not found! Please check path.'.format(args.rawfile))

    print('\nInput file: {}\n'.format(args.rawfile))

    smt_h5 = MODISsmth5(rawfile = args.rawfile, tempint = args.tempint, targetdir = args.targetdir, parallel = args.parallel, ncores = args.ncores)

    if not smt_h5.exists:
        smt_h5.create()

    # check if v-curve optimization is true

    if args.loptimize:

        if args.llas:

            try:
                llas = array.array('f',np.linspace(float(args.llas[0]),float(args.llas[1]),float(args.llas[1])/float(args.llas[2]) + 1.0))
            except (IndexError,TypeError):
                raise SystemExit('Error with lambda array values. Expected tuple of float log10(lambda) - (lmin,lmax,lstep)!')
        else:
            llas = array.array('f',np.linspace(0.0,4.0,41.0))

        if args.pvalue:

            try:
                p = float(args.pvalue)
            except ValueError:
                raise SystemExit('Error using p-value ... expected float, got {}.'.format(type(args.pvalue)))

            print('\nRunning asymmetric whittaker smoother with v-cuve optimization ... \n')

            smt_h5.ws2d_vc_asy(llas=llas,p=p)

        else:

            print('\nRunning whittaker smoother with v-cuve optimization ... \n')

            smt_h5.ws2d_vc(llas=llas)

    else:

        if args.Lambda:

            try:
                l = 10**float(args.Lambda)
            except:
                raise SystemExit('Error with lambda value. Expected float log10(lambda)!')

            print('\nRunning whittaker smoother with fixed lambda ... \n')

            smt_h5.ws2d(l=l)

        else:

            print('\nRunning whittaker smoother with lambda from grid ... \n')

            smt_h5.ws2d_lgrid()

    print('\n[{}]: smoothMODIS.py finished successfully.\n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

if __name__ == '__main__':
    main()
