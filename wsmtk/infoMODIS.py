#!/usr/bin/env python
from __future__ import print_function
import h5py
import sys,os
import time
import argparse



def main():

    parser = argparse.ArgumentParser(description="Get MODIS raw/smooth file info")
    parser.add_argument("file", help='Full path to MODIS h5 file')

    # fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if not os.path.isfile(args.file):
        raise SystemExit('File not found!')

    msg_head = "WSMTK info tool - {}".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

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

            try:
                lr = dset.attrs['lastrun']
            except KeyError:
                lr = None
    except:
        raise SystemExit('Error reading file information.')

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



    print(' ',msg_head.center(os.get_terminal_size().columns),' ',sep='\n')
    print(msg)




if __name__ == '__main__':
    main()
