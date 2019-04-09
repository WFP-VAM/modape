#!/usr/bin/env python
# pylint: disable=broad-except
"""modis_collect.py: Collect raw MODIS data into HDF5 file."""

from __future__ import absolute_import, division, print_function

import argparse
import glob
import multiprocessing as mp
import os
import re
import sys
import time
import traceback

from modape.modis import ModisRawH5

def run_process(pdict):
    """Execute processing into raw HDF5 file

    Little wrapper to execute the processing into a raw HDF5 file. This enables multi-process concurrency.

    Args:
        pdict: dictionary with processing parameters for tile
    """

    for vam_product_code in pdict['vam_product_code']:
        try:
            rh5 = ModisRawH5(pdict['files'],
                             vam_product_code=vam_product_code,
                             targetdir=pdict['targetdir'],
                             interleave=pdict['interleave'])
            if not rh5.exists:
                rh5.create(compression=pdict['compression'],
                           chunk=pdict['chunksize'])
            rh5.update()
        except Exception as e: # pylint: disable=unused-variable
            print('\nError processing product {}, product code {}. \n\n Traceback:\n'.format(rh5.product, vam_product_code))
            traceback.print_exc()
        print('\n')

def main():
    """Collect raw MODIS hdf files into a raw MODIS HDF5 file.

    All MODIS hdf files within srcdir will be collected into a raw MODIS HDF5 file, corresponding to product type and tile.
    If the respective HDF5 file does not exists in the target directory, it will be created. Otherwhise, the file will be
    updated and the data inserted at the proper temporal location within the HDF5 file.

    16-day MOD13* and MYD13* products can be interleaved into an 8-day product with the new product ID MXD* by adding the `--interleave` flag.
    """

    parser = argparse.ArgumentParser(description='Process downloaded RAW MODIS hdf files')
    parser.add_argument('srcdir', help='directory with raw MODIS .hdf files', default=os.getcwd(), metavar='srcdir')
    parser.add_argument('-d', '--targetdir', help='Target directory for PROCESSED MODIS files (default is scrdir)', metavar='')
    parser.add_argument('-x', '--compression', help='Compression for HDF5 files', default='gzip', metavar='')
    parser.add_argument('-c', '--chunksize', help='Number of pixels per block (value needs to result in integer number of blocks)', type=int, metavar='')
    parser.add_argument('--all-vampc', help='Flag to process all possible VAM product codes', action='store_true')
    parser.add_argument('--interleave', help='Interleave MOD13 & MYD13 products to MXD (only works for VIM!)', action='store_true')
    parser.add_argument('--parallel-tiles', help='Number of tiles processed in parallel (default = None)', default=1, type=int, metavar='')
    parser.add_argument('--quiet', help='Be quiet', action='store_true')

    # Fail and print help if no arguments supplied
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()
    if not os.path.isdir(args.srcdir):
        raise SystemExit('Source directory not a valid path! Please check inputs.targetdir')
    if not args.targetdir:
        args.targetdir = args.srcdir # default is srcdir

    files = glob.glob('{}/*.hdf'.format(args.srcdir))

    # Regex patterns
    ppatt = re.compile(r'M\w{6}')
    vpatt = re.compile(r'.+\.(\d{3})\..+')
    tpatt = re.compile(r'h\d+v\d+')
    vimvem = re.compile('M.D13')
    lst = re.compile('M.D11')
    processing_dict = {} # processing dictionary

    # Seperate input files into group
    groups = ['.*'.join(re.findall(ppatt, os.path.basename(x)) + re.findall(tpatt, os.path.basename(x)) + [re.sub(vpatt, '\\1', os.path.basename(x))])  for x in files]

    if args.interleave:
        groups = list({re.sub('(M.{1})(D.+)', 'M.'+'\\2', x) if re.match(vimvem, x) else x for x in groups}) # Join MOD13/MYD13
    else:
        groups = list(set(groups))
    for group in groups:
        gpatt = re.compile(group + '.*hdf')
        processing_dict[group] = {}
        processing_dict[group]['targetdir'] = args.targetdir
        processing_dict[group]['files'] = [x for x in files if re.search(gpatt, x)]
        processing_dict[group]['interleave'] = args.interleave

        if args.chunksize:
            processing_dict[group]['chunksize'] = args.chunksize
        else:
            processing_dict[group]['chunksize'] = None
        processing_dict[group]['compression'] = args.compression

        if args.all_vampc:
            if re.match(vimvem, group.split('.*')[0]):
                processing_dict[group]['vam_product_code'] = ['VIM', 'VEM']
            elif re.match(lst, group.split('.*')[0]):
                processing_dict[group]['vam_product_code'] = ['LTD', 'LTN']
            else:
                raise ValueError('No VAM product code implemented for {}'.format(group.split('.*')[0]))
        else:
            processing_dict[group]['vam_product_code'] = [None]

    if args.parallel_tiles > 1:
        if not args.quiet:
            print('\n\n[{}]: Start processing - {} tiles in parallel ...'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()), args.parallel_tiles))

        pool = mp.Pool(processes=args.parallel_tiles)
        _ = pool.map(run_process, [processing_dict[group] for group in groups])
        pool.close()
        pool.join()

        if not args.quiet:
            print('[{}]: Done.'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
    else:
        if not args.quiet:
            print('\n\n[{}]: Start processing ... \n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))

        for group in groups:
            run_process(processing_dict[group])
        if not args.quiet:
            print('\n[{}]: Done.'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))

if __name__ == '__main__':
    main()
