#!/usr/bin/env python
"""modis_info.py: Return metadata stored in MODIS HDF5 files."""

from __future__ import absolute_import, division, print_function
import argparse
import os
import sys
import time

import h5py ## pylint: disable=import-error

def main():
    """Info tool for processed MODIS HDF5 files.

    Returns metadata on processed MODIS HDF5 files, both for raw and smoothed files.
    """

    parser = argparse.ArgumentParser(description='Get MODIS raw/smooth file info')
    parser.add_argument('file', help='Full path to MODIS h5 file')

    # Fail and print help if no arguments supplied
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    # Check if file exists
    if not os.path.isfile(args.file):
        raise SystemExit('File not found!')

    # Message head
    message_head = 'MODAPE info tool - {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()))

    # Read metadata
    try:
        with h5py.File(args.file, 'r') as h5f:
            dset = h5f.get('data')
            dates = h5f.get('dates')
            dim = dset.shape
            startdate = dates[0].decode()
            enddate = dates[-1].decode()
            temporalresolution = dset.attrs['temporalresolution']
            resolution = dset.attrs['resolution']
            nodata_value = dset.attrs['nodata']
            processing_timestamp = dset.attrs['processingtimestamp']
            ncols = dset.attrs['RasterXSize'].item()
            nrows = dset.attrs['RasterYSize'].item()

            # If reading lastrun fails, it's assumed the product is a raw HDF5 file
            try:
                last_run = dset.attrs['lastrun']
            except KeyError:
                last_run = None
    except:
        raise SystemExit('Error reading file information.')

    # If lastrun attribute is not none (aka product is smooth HDF5)
    if last_run:
        message = '''
File: {}

Type: MODIS smoothed HDF5

Dimensions:

    - {} rows

    - {} columns

    - {} temporal units

Start date: {}

End date: {}

Temporal resolution: {} daily

Spatial resolution: {} m

NoData value: {}

Last modified: {}

Last smoothing run: Whittaker smoother with {}\n'''.format(args.file, nrows, ncols, dim[1],
                                                           startdate, enddate, temporalresolution, resolution,
                                                           nodata_value, processing_timestamp, last_run)

    else:
        message = '''
File: {}

Type: MODIS raw daily HDF5

Dimensions:

    - {} rows

    - {} columns

    - {} temporal units

Start date: {}

End date: {}

Temporal resolution: {} daily

Spatial resolution: {} m

NoData value: {}

Last modified: {}\n'''.format(args.file, nrows, ncols,
                              dim[1], startdate, enddate, temporalresolution,
                              resolution, nodata_value, processing_timestamp)

    # Print message - header is centered
    print(' ', message_head.center(os.get_terminal_size().columns), ' ', sep='\n')
    print(message)

if __name__ == '__main__':
    main()
