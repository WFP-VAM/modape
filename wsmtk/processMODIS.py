#!/usr/bin/env python
from __future__ import print_function
from wsmtk.modis import MODISrawh5
import os
import sys
import glob
import argparse
import datetime
import re
import traceback
from contextlib import contextmanager, closing
import multiprocessing as mp


def execute_process(pdict):
    '''Execute processing into raw HDF5 file

    Little wrapper to execute the processing into a raw HDF5 file. This enables multi-process concurrency.

    Args:
        pdict: dictionary with processing parameters for tile
    '''

    for par in pdict['parameters']:

        try:

            rh5 = MODISrawh5(pdict['files'],param=par,targetdir=pdict['targetdir'])

            # Creatre if file doesn't exist yet
            if not rh5.exists:
                rh5.create(compression=pdict['compression'], chunk=pdict['chunksize'])
            rh5.update()

        except Exception as e:

            print('\nError processing product {}, parameter {}. \n\n Traceback:\n'.format(rh5.product,par))

        traceback.print_exc()
        print('\n')


def main():
    '''Collect raw MODIS hdf files into a raw MODIS HDF5 file.

    All MODIS hdf files within srcdir will be collected into a raw MODIS HDF5 file, corresponding to product type and tile.
    If the respective HDF5 file does not exists in the target directory, it will be created. Otherwhise, the file will be
    updated and the data inserted at the proper temporal location within the HDF5 file.

    By default, 16-day MOD13* and MYD13* products will be interleaved into an 8-day product with the new product ID MXD*.
    '''

    parser = argparse.ArgumentParser(description="Process downloaded RAW MODIS hdf files")
    parser.add_argument("srcdir", help='directory with raw MODIS .hdf files',default=os.getcwd(),metavar='srcdir')
    #parser.add_argument("-p","--parameter", help='VAM parameter code',metavar='') ## paramter selection not implemented
    parser.add_argument("-d","--targetdir", help='Target directory for PROCESSED MODIS files (default is scrdir)',metavar='')
    parser.add_argument("-x","--compression", help='Compression for HDF5 files',default='gzip',metavar='')
    parser.add_argument("--all-parameters", help='Flag to process all possible VAM parameters',action='store_true')
    parser.add_argument("-c","--chunksize", help='Number of pixels per block (value needs to result in integer number of blocks)',type=int,metavar='')
    parser.add_argument("--parallel-tiles", help='Number of tiles processed in parallel (default = None)',default=1,type=int,metavar='')

    # Fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Check if srcdir exists
    if not os.path.isdir(args.srcdir):
        raise SystemExit('Source directory not a valid path! Please check inputs.targetdir')

    # Default for targetdir is srcdir
    if not args.targetdir:
        args.targetdir = args.srcdir

    # Find files
    files = glob.glob('{}/*.hdf'.format(args.srcdir))

    # Regex patterns
    ppatt = re.compile(r'M\w{6}')
    vpatt = re.compile('.+\.(\d{3})\..+')
    tpatt = re.compile(r'h\d+v\d+')
    vimvem = re.compile('M.D13')
    lst = re.compile('M.D11')


    # processin dictionary

    processing_dict = {}

    # Seperate input files into group
    groups = ['.*'.join(re.findall(ppatt,os.path.basename(x)) + re.findall(tpatt,os.path.basename(x)) + [re.sub(vpatt,'\\1',os.path.basename(x))])  for x in files]

    # Join MOD13/MYD13
    groups = list(set([re.sub('(M.{1})(D.+)','M.'+'\\2',x) if re.match(vimvem,x) else x for x in groups]))


    for g in groups:

        gpatt = re.compile(g + '.*hdf')

        processing_dict[g] = {}

        processing_dict[g]['targetdir'] = args.targetdir

        processing_dict[g]['files'] = [x for x in files if re.search(gpatt,x)]

        if args.chunksize:

            processing_dict[g]['chunksize'] = args.chunksize

        else:

            processing_dict[g]['chunksize'] = None

        processing_dict[g]['compression'] = args.compression


        if args.all_parameters:

            if re.match(vimvem,g.split('.*')[0]):

                processing_dict[g]['parameters'] = ['VIM','VEM']

            elif re.match(lst,g.split('.*')[0]):

                processing_dict[g]['parameters'] = ['LTD','LTN']

            else:
                raise SystemExit('No parameters implemented for {}'.format(g.split('.*')[0]))
        else:
            processing_dict[g]['parameters'] = [None]


    if args.parallel_tiles > 1:

        with closing(mp.Pool(processes=args.parallel_tiles)) as pool:

            res = pool.map(execute_process,[processing_dict[g] for g in groups])

        pool.close()
        pool.join()

    else:

        for g in groups:

            execute_process(processing_dict[g])



if __name__ == '__main__':
    main()
