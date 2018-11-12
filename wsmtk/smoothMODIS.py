#!/usr/bin/env python
from __future__ import print_function
from wsmtk.modis import MODISsmth5
from wsmtk.utils import init_parameters
import shutil
import os
import sys
import glob
import argparse
import multiprocessing as mp
import numpy as np
import time
from progress.bar import Bar

def initfun(pdict_):
    '''Initfun for worker'''

    global pdict
    pdict = pdict_

def run_ws2d(h5):
    '''Run smoother with fixed s

    Args:
        pdict: dictionary with processing parameters for tile
    '''

    if not os.path.isfile(h5):

        print('Raw HDF5 {} not found! Please check path.'.format(h5))

    else:

        smt_h5 = MODISsmth5(rawfile = h5, tempint = pdict['tempint'], nsmooth = pdict['nsmooth'], nupdate = pdict['nupdate'], targetdir = pdict['targetdir'], nworkers = pdict['nworkers'])

        if not smt_h5.exists:
            smt_h5.create()

        smt_h5.ws2d(pdict['s'])


def run_ws2d_sgrid(h5):
    '''Run smoother with fixed s from grid

    Args:
        pdict: dictionary with processing parameters for tile
    '''
    if not os.path.isfile(h5):

        print('Raw HDF5 {} not found! Please check path.'.format(h5))

    else:

        smt_h5 = MODISsmth5(rawfile = h5, tempint = pdict['tempint'], nsmooth = pdict['nsmooth'], nupdate = pdict['nupdate'], targetdir = pdict['targetdir'], nworkers = pdict['nworkers'])

        if not smt_h5.exists:
            smt_h5.create()

        smt_h5.ws2d_sgrid()


def run_ws2d_vOpt(h5):
    '''Run smoother with V-curve optimization of s

    Args:
        pdict: dictionary with processing parameters for tile
    '''
    if not os.path.isfile(h5):

        print('Raw HDF5 {} not found! Please check path.'.format(h5))

    else:

        smt_h5 = MODISsmth5(rawfile = h5, tempint = pdict['tempint'], nsmooth = pdict['nsmooth'], nupdate = pdict['nupdate'], targetdir = pdict['targetdir'], nworkers = pdict['nworkers'])

        if not smt_h5.exists:
            smt_h5.create()

        smt_h5.ws2d_vOpt(pdict['srange'],pdict['pvalue'])


def main():
    '''Smooth, gapfill and interpolate processed raw MODIS HDF5 files.

    The smoothing function takes a previously created raw MODIS HDF file (as created by processMODIS) as input.
    The raw data can be smoothed with eiter a fixed s value, a pixel-by-pixel s value read from a previously computed grid or
    V-curve optimization of s (creates or updates the s-grid)

    The desired temporal resolution of the ouput file can be defined with tempint, allowing for seamless interpolation.

    If a smooth MODIS HDF5 file for a given product, tile (if not global) and temporal interpolation is already in
    the targetdir, it will be updated.

    By default, the entire temporal range of the raw data is used for smoothing, and the entire smoothed data is updated.
    The parameters nsmooth and nupdate can modify this behaviour.

    To speed up processing time, parallel processing can be leveraged, processing a user defined number of tiles in parallel, and for each tile splitting up the task into n worker processes (by default,
    no concurrency is enabled)
    '''

    parser = argparse.ArgumentParser(description="Smooth, gapfill and interpolate processed raw MODIS HDF5 files")
    parser.add_argument("input", help='Smoothing input - either one or more raw MODIS HDF5 file(s) or path containing raw MODIS HDF5 file(s)',metavar='input')
    parser.add_argument("-s","--svalue", help='S value for smoothing (has to be log10(s)', metavar='', type = float)
    parser.add_argument("-S","--srange", help='S value range for V-curve (float log10(s) values as smin smax sstep - default -1.0 2.0 0.2)',nargs='+',metavar='')
    parser.add_argument("-t","--tempint", help='Value for temporal interpolation (integer required - default is native temporal resolution AKA no interpolation)', metavar='',type = int)
    parser.add_argument("-n","--nsmooth", help='Number of raw timesteps used for smoothing', metavar='',type = int)
    parser.add_argument("-u","--nupdate", help='Number of smoothed timesteps to be updated in HDF5 file', metavar='',type = int)
    parser.add_argument("-p","--pvalue", help='Value for asymmetric smoothing (float required)', metavar='', type = float)
    parser.add_argument("-d","--targetdir", help='Target directory for smoothed output',default=os.getcwd(),metavar='')
    parser.add_argument("--soptimize", help='Use V-curve for s value optimization',action='store_true')
    parser.add_argument("--parallel-tiles", help='Number of tiles processed in parallel (default = None)',default=1,type=int,metavar='')
    parser.add_argument("--nworkers", help='Number of worker processes used per tile (default is number is 1 - no concurrency)',default=1, metavar='', type = int)
    parser.add_argument("--quiet", help='Be quiet',action='store_true')

    # Fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        print('\n')

    args = parser.parse_args()

    if len(args.input) == 1 & os.path.isdir(args.input[0]):
        files = glob.glob(args.input[0] + '/*.h5')
    else:
        files = args.input

    assert len(files) > 0; "No files found to process"

    if not os.path.isdir(args.targetdir):
        print("Targetdir {} doesn't exist. Creating ... ",end='')
        shutil.os.mkdir(args.targetdir)
        print('done.')

    if args.srange:
        try:
            assert len(args.srange) == 3
            args.srange = np.linspace(float(args.srange[0]),float(args.srange[1]),float(args.srange[1])/float(args.srange[2]) + 1.0)
        except (IndexError,TypeError,AssertionError):
            raise SystemExit('Error with s value array values. Expected three values of float log10(s) -  smin smax sstep !')

    if args.svalue:
        try:
            processing_dict['s'] = 10**float(args.svalue)
        except:
            raise SystemExit('Error with s value. Expected float log10(s)!')

    # prepare processing dict

    processing_dict = init_parameters(tempint=args.tempint,nsmooth=args.nsmooth,nupdate=nupdate,targetdir=agrs.targetdir,nworkers=args.nworkers)

    if not args.quiet:
        print('\n[{}]: Starting smoothMODIS.py ... \n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    if args.parallel_tiles > 1:

        # Check if V-curve optimization is true
        if args.soptimize:

            if not args.srange:

                processing_dict['srange'] = np.linspance(-1.0,2.0,16.0)

            else:

                processing_dict['srange'] = args.srange

            processing_dict['pvalue'] = args.pvalue

            if not args.quiet:
                print('\nRunning whittaker smoother V-curve optimization ... \n')

            with closing(mp.Pool(processes=args.parallel_tiles,initializer = initfun, initargs = (processing_dict,))) as pool:

                res = pool.map(run_ws2d_vOpt,files)

            pool.close()
            pool.join()

            if not args.quiet:
                print('[{}]: Done.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        elif args.svalue:

            processing_dict['s'] = args.svalue

            if not args.quiet:
                print('\nRunning whittaker smoother with fixed s value ... \n')

            with closing(mp.Pool(processes=args.parallel_tiles,initializer = initfun, initargs = (processing_dict,))) as pool:

                res = pool.map(run_ws2d,files)

            pool.close()
            pool.join()

            if not args.quiet:
                print('[{}]: Done.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        else:

            if not args.quiet:
                print('\nRunning whittaker smoother with s value from grid ... \n')

            with closing(mp.Pool(processes=args.parallel_tiles,initializer = initfun, initargs = (processing_dict,))) as pool:

                res = pool.map(run_ws2d_sgrid,files)

            pool.close()
            pool.join()

            if not args.quiet:
                print('[{}]: Done.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    else:

        if args.soptimize:

            if not args.quiet:
                print('\nRunning whittaker smoother V-curve optimization ... \n')
                bar = Bar('Processing',fill='=',max=len(files),suffix='%(percent)d%%  ')
                bar.goto(0)

            for h5 in files:

                if not os.path.isfile(h5):
                    print('Raw HDF5 {} not found! Please check path.'.format(h5))
                    continue

                smt_h5 = MODISsmth5(rawfile = h5, tempint = args.tempint, nsmooth = args.nsmooth, nupdate = args.nupdate, targetdir = args.targetdir, nworkers = args.nworkers)

                if not smt_h5.exists:
                    smt_h5.create()

                smt_h5.ws2d_vOpt(args.sgrid,args.pvalue)

                if not args.quiet:
                    bar.next()

            if not args.quiet:
                bar.finish()
                print('[{}]: Done.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        elif args.svalue:

            if not args.quiet:
                print('\nRunning whittaker smoother with fixed s value ... \n')
                bar = Bar('Processing',fill='=',max=len(files),suffix='%(percent)d%%  ')
                bar.goto(0)

            for h5 in files:

                if not os.path.isfile(h5):
                    print('Raw HDF5 {} not found! Please check path.'.format(h5))
                    continue

                smt_h5 = MODISsmth5(rawfile = h5, tempint = args.tempint, nsmooth = args.nsmooth, nupdate = args.nupdate, targetdir = args.targetdir, nworkers = args.nworkers)

                if not smt_h5.exists:
                    smt_h5.create()

                smt_h5.ws2d(args.svalue)

                if not args.quiet:
                    bar.next()

            if not args.quiet:
                bar.finish()
                print('[{}]: Done.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        else:

            if not args.quiet:
                print('\nRunning whittaker smoother with s value from grid ... \n')
                bar = Bar('Processing',fill='=',max=len(files),suffix='%(percent)d%%  ')
                bar.goto(0)

            for h5 in files:

                if not os.path.isfile(h5):
                    print('Raw HDF5 {} not found! Please check path.'.format(h5))
                    continue

                smt_h5 = MODISsmth5(rawfile = h5, tempint = args.tempint, nsmooth = args.nsmooth, nupdate = args.nupdate, targetdir = args.targetdir, nworkers = args.nworkers)

                if not smt_h5.exists:
                    smt_h5.create()

                smt_h5.ws2d_sgrid()

                if not args.quiet:
                    bar.next()

            if not args.quiet:
                bar.finish()
                print('[{}]: Done.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    print('\n[{}]: smoothMODIS.py finished successfully.\n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

if __name__ == '__main__':
    mp.freeze_support()
    main()
