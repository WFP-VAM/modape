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

def main():
    '''Collect raw MODIS hdf files into a daily raw MODIS HDF5 file.

    All MODIS hdf files within srcdir will be collected into a raw MODIS HDF5 file, corresponding to product type and tile.
    If the respective HDF5 file does not exists in the target directory, it will be created. Otherwhise, the file will be
    updated and the data inserted at the proper temporal location within the HDF5 file.

    The 2D raw data will be converted into 3D daily data, where the 3rd dimension of the array corresponds to the temporal
    resolution of the input product. If available, the compositing DOYs will be used to assing the values to the corresponding
    pixel within the 3D array. If such information is not availabe, the data will be set to the mid-point of the compositing period.

    By default, 16-day MOD13* and MYD13* products will be interleaved into an 8-day product with the new product ID MXD*.
    '''

    parser = argparse.ArgumentParser(description="Process downloaded RAW MODIS hdf files")
    parser.add_argument("srcdir", help='directory with raw MODIS .hdf files',default=os.getcwd(),metavar='srcdir')
    #parser.add_argument("-p","--parameter", help='VAM parameter code',metavar='') ## paramter selection not implemented
    parser.add_argument("-d","--targetdir", help='Target directory for PROCESSED MODIS files (default is scrdir)',metavar='')
    parser.add_argument("-c","--compression", help='Compression for HDF5 files',default='gzip',metavar='')
    parser.add_argument("--all-parameters", help='Flag to process all possible VAM parameters',action='store_true')
    parser.add_argument("-b","--blocksize", help='Minimum values for row & columns per processing block (default 120 120)',nargs=2,default=[120,120],type=int,metavar='')

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


    # Seperate input files into group
    groups = ['.*'.join(re.findall(ppatt,os.path.basename(x)) + re.findall(tpatt,os.path.basename(x)) + [re.sub(vpatt,'\\1',os.path.basename(x))])  for x in files]

    # Join MOD13/MYD13
    groups = list(set([re.sub('(M.{1})(D.+)','M.'+'\\2',x) if re.match(vimvem,x) else x for x in groups]))

    # If all parameters are requested
    if args.all_parameters:

        # Iterate over groups
        for g in groups:

            gpatt = re.compile(g + '.*hdf')

            # Subset file list
            files_sub = [x for x in files if re.search(gpatt,x)]

            if re.match(vimvem,g.split('.*')[0]):
                allps = ['VIM','VEM']
            elif re.match(lst,g.split('.*')[0]):
                allps = ['LTD','LTN']
            else:
                ## maybe change to warning?
                raise SystemExit('No parameters implemented for {}'.format(g.split('.*')[0]))

            # Iterate over parameters
            for p in allps:

                # Create raw MODIS HDF5 file object
                try:
                    h5 = MODISrawh5(files_sub,param=p,targetdir=args.targetdir,compression=args.compression,crow=args.blocksize[0],ccol=args.blocksize[1])

                    # Creatre if file doesn't exist yet
                    if not h5.exists:
                        h5.create()
                    h5.update()

                except Exception as e:
                    print('\nError processing file group {}, parameter {}. Skipping to next group/parameter!. \n\n Traceback:\n'.format(g,p))
                    traceback.print_exc()
                    print('\n')
                    continue
    else:

        for g in groups:

            gpatt = re.compile(g + '.*hdf')

            files_sub = [x for x in files if re.search(gpatt,x)]

            try:

                h5 = MODISrawh5(files_sub,targetdir=args.targetdir,compression=args.compression,crow=args.blocksize[0],ccol=args.blocksize[1])

                if not h5.exists:
                    h5.create()

                h5.update()

            except Exception as e:
                print('\nError processing file group {}. Skipping to next group!. \n\n Traceback:\n'.format(g))
                traceback.print_exc()
                print('\n')
                continue

if __name__ == '__main__':
    main()
