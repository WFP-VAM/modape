#!/usr/bin/env python
"""csv_smooth.py: Smooth timeseries in a CSV file."""

from __future__ import absolute_import, division, print_function

import argparse
import array
import os
import sys
import time

import numpy as np
import pandas as pd # pylint: disable=E0401

from modape.whittaker import ws2d, ws2doptv, ws2doptvp # pylint: disable=E0611

def main():
    """Smooth timeseries in a CSV file.

    By default, the first row and the first column in the CSV file is skipped. The following three rows should be NAME/ID, LON and LAT.
    Starting row 4, the rows are interpreted as raw data values. Each column starting from number 2, is interpreted as separate timeseries.

    The defined or determined sopt and log10(sopt) values are appended to the end of the smoothed timeseries.

    The first three letters of the CSV filename are interpreted as region code and are included in the output filename.

    In addition to the region code, the output filename contains a suffix which indicates the smoothing method used:

        - fixed s: filt0.csv
        - V-curve: filtoptv.csv
        - asymmetric V-curve: filtoptvp.csv

    The resulting CSV is created in the directory the input file is located.
    """

    parser = argparse.ArgumentParser(description='Smooth CSV file')
    parser.add_argument('file', help='CSV file')
    parser.add_argument('-s', '--svalue', help='S value for smoothing (has to be log10(s)', metavar='', type=float)
    parser.add_argument('-S', '--srange', help='S range for V-curve (float log10(s) values as smin smax sstep - default 0.0 4.0 0.1)', nargs='+', metavar='', type=float)
    parser.add_argument('-p', '--pvalue', help='Value for asymmetric smoothing (float required)', metavar='', type=float)

    # Fail and print help if no arguments supplied
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.isfile(args.file):
        raise SystemExit('Input CSV file {} not found! Please check path.'.format(args.file))

    print('\n[{}]: Starting smoothCSV.py ... \n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))

    # Create filename of output CSV
    outname = '{}/{}'.format(os.path.dirname(os.path.abspath(args.file)), os.path.basename(args.file)[0:3])

    df = pd.read_csv(args.file, header=1) # Read input
    resdf = pd.DataFrame(index=range(len(df)+2)) #result dataframe, +2 for Sopt

    # Add ID column
    resdf['ID'] = pd.concat([pd.Series(['Lon', 'Lat']),
                             pd.Series(np.linspace(1, len(df)-2, len(df)-2)),
                             pd.Series(['Sopt', 'logSopt'])],
                            ignore_index=True)
    tmparr = np.zeros(len(df)-2, dtype='double') # Initialize array


    if args.svalue:
        s = 10**args.svalue
        outname = outname + 'filt0.csv'
        print('\nSmoothing using fixed S value {}. Writing to file: {}\n'.format(s, outname))

        # Iterate columns (skip 1st)
        for col in df.columns[1:]:
            val = df[col].values
            tmparr[...] = val[2:]
            val[2:] = ws2d(tmparr, s, np.array((tmparr > 0)*1, dtype='double'))
            resdf[col] = pd.concat([pd.Series(val),
                                    pd.Series([s, np.log10(s)])],
                                   ignore_index=True)
    else:
        if args.srange:
            if len(args.srange) != 3:
                raise ValueError('Expected 3 inputs for S range: smin smax step!')
            try:
                srange = array.array('d',
                                     np.linspace(args.srange[0],
                                                 args.srange[1],
                                                 args.srange[1]/args.srange[2]+1))
            except:
                print('Error parsing S range values')
                raise
        else:
            srange = array.array('d', np.linspace(0.0, 4.0, 41))
            args.srange = [0.0, 4.0, 0.1]

        if args.pvalue:
            outname = outname + 'filtoptvp.csv'
            resdf = pd.DataFrame(resdf['ID'].append(pd.Series('pvalue'),
                                                    ignore_index=True),
                                 columns=['ID'])

            print('\nSmoothing using asymmetric V-curve optimization with smin:{}, smax:{}, sstep:{} and pvalue:{}.\n\nWriting to file: {}\n'
                  .format(args.srange[0], args.srange[1], args.srange[2], args.pvalue, outname))

            for col in df.columns[1:]:
                val = df[col].values
                tmparr[...] = val[2:]

                val[2:], sopt = ws2doptvp(tmparr,
                                          np.array((tmparr > 0)*1, dtype='double'),
                                          srange,
                                          args.pvalue)

                resdf[col] = pd.concat([pd.Series(val),
                                        pd.Series([sopt, np.log10(sopt)]),
                                        pd.Series(args.pvalue)],
                                       ignore_index=True)
        else:
            outname = outname + 'filtoptv.csv'
            print('\nSmoothing using V-curve optimization with smin:{}, smax:{}, sstep:{}.\n\nWriting to file: {}\n'
                  .format(args.srange[0], args.srange[1], args.srange[2], outname))

            for col in df.columns[1:]:
                val = df[col].values
                tmparr[...] = val[2:]
                val[2:], sopt = ws2doptv(tmparr, np.array((tmparr > 0)*1, dtype='double'), srange)
                resdf[col] = pd.concat([pd.Series(val),
                                        pd.Series([sopt, np.log10(sopt)])],
                                       ignore_index=True)

    # Write to disk
    resdf.to_csv(outname, index=False)
    print('\n[{}]:smoothCSV.py finished.\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))

if __name__ == '__main__':
    main()
