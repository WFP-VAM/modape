#!/usr/bin/env python
from __future__ import print_function
import h5py
import sys,os
import time
import argparse



def main():
    '''Info tool for processed MODIS HDF5 files.

    Returns metadata on processed MODIS HDF5 files, both for raw and smoothed files.
    '''

    parser = argparse.ArgumentParser(description="Get MODIS raw/smooth file info")
    parser.add_argument("file", help='Full path to MODIS h5 file')

    # Fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Check if file exists
    if not os.path.isfile(args.file):
        raise SystemExit('File not found!')

    # Message head
    msg_head = "WSMTK info tool - {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    # Read metadata
    try:

        with h5py.File(args.file,'r') as h5f:

            dset = h5f.get('data')
            dts  = h5f.get('dates')

            dim = dset.shape
            startdate = dts[0].decode()
            enddate = dts[-1].decode()
            temporalresolution = dset.attrs['temporalresolution']
            nodata_value = dset.attrs['nodata']
            ptstmp = dset.attrs['processingtimestamp']

            # If reading lastrun fails, it's assumed the product is a raw HDF5 file
            try:
                lr = dset.attrs['lastrun']
            except KeyError:
                lr = None
    except:
        raise SystemExit('Error reading file information.')

    # If lastrun attribute is not none (aka product is smooth HDF5)
    if lr:

        msg = '''
File: {}

Type: MODIS smoothed HDF5

Dimensions:

    - {} rows

    - {} columns

    - {} temporal units

Start date: {}

End date: {}

Temporal resolution: {} daily

NoData value: {}

Last modified: {}

Last smoothing run: Whittaker smoother with {}\n'''.format(args.file,dim[0],dim[1],dim[2],startdate,enddate,temporalresolution,nodata_value,ptstmp,lr)


    else:

        msg = '''
File: {}

Type: MODIS raw daily HDF5

Dimensions:

    - {} rows

    - {} columns

    - {} temporal units

Start date: {}

End date: {}

Temporal resolution: {} daily

NoData value: {}

Last modified: {}\n'''.format(args.file,dim[0],dim[1],dim[2],startdate,enddate,temporalresolution,nodata_value,ptstmp)


    # Print message - header is centered
    print(' ',msg_head.center(os.get_terminal_size().columns),' ',sep='\n')
    print(msg)




if __name__ == '__main__':
    main()
