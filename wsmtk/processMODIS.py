#!/usr/bin/env python
from __future__ import print_function
from wsmtk.modis import MODIShdf5
import os
import sys
import glob
import argparse
import datetime
import re

def main():

    parser = argparse.ArgumentParser(description="Process downloaded RAW MODIS hdf files")
    parser.add_argument("srcdir", help='directory with raw MODIS .hdf files',default=os.getcwd(),metavar='DIR')
    #parser.add_argument("-p","--parameter", help='VAM parameter code',metavar='') ## paramter selection not implemented
    parser.add_argument("--prcdir", help='Storage directory for PROCESSED MODIS files',default=os.getcwd(),metavar='')
    #parser.add_argument("--rawdir", help='Storage directory for RAW MODIS files',metavar='')
    parser.add_argument("-c","--compression", help='Compression for HDF5 files',default=32001,metavar='')
    parser.add_argument("--all-parameters", help='Flag to process all possible VAM parameters',action='store_true')

    # fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if not os.path.isdir(args.srcdir):
        raise SystemExit('Source directory not a valid path! Please check inputs.targetdir')


    files = glob.glob('{}/*.hdf'.format(args.srcdir))

    # seperate input files into group

    ppatt = re.compile(r'M\w{6}')
    vpatt = re.compile('.+\.(\d{3})\..+')
    tpatt = re.compile(r'h\d+v\d+')

    groups = list(set(['.*'.join(re.findall(ppatt,os.path.basename(x)) + re.findall(tpatt,os.path.basename(x)) + [re.sub(vpatt,'\\1',os.path.basename(x))])  for x in files]))

    # if all parameters are requested
    if args.all_parameters:

        vimvem = re.compile('M.D13')
        lst = re.compile('M.D11')

        for g in groups:

            gpatt = re.compile(g + '.*hdf')

            files_sub = [x for x in files if re.search(gpatt,x)]

            if re.match(vimvem,g.split('.*')[0]):
                allps = ['VIM','VEM']
            elif re.match(lst,g.split('.*')[0]):
                allps = ['LTD','LTN']
            else:
                ## maybe change to warning?
                raise SystemExit('No parameters implemented for {}'.format(g.split('.*')[0]))

            for p in allps:

                try:
                    h5 = MODIShdf5(files_sub,param=p,targetdir=args.prcdir,compression=args.compression)
                    if not h5.exists:
                        h5.create()
                    h5.update()

                except Exception as e:
                    print('\nError processing file group {}, parameter {}. Skipping to next group/parameter!. \n\n (Exception raised: {})\n)'.format(g,p,e))
                    continue
    else:

        for g in groups:

            gpatt = re.compile(g + '.*hdf')

            files_sub = [x for x in files if re.search(gpatt,x)]

            try:

                h5 = MODIShdf5(files_sub,targetdir=args.prcdir,compression=args.compression)

                if not h5.exists:
                    h5.create()

                h5.update()

            except Exception as e:
                print('\nError processing file group {}. Skipping to next group!. \n\n (Exception raised: {})\n)'.format(g,e))
                continue

if __name__ == '__main__':
    main()
