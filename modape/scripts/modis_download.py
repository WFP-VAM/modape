#!/usr/bin/env python
"""modis_download.py: Query and download MODIS HDF files."""

import argparse
from datetime import datetime
import os
from pathlib import Path
import sys
import warnings

from modape.modis import ModisQuery

warnings.filterwarnings("default", category=DeprecationWarning)

def main():
    '''Query and download MODIS products.

    This function allows for querying and downloading MODIS products in bulk.
    Multiple products can be queried and downloaded with one
    function call. For downloading data, valid earthdata credentials are required
    (to register, visit https://urs.earthdata.nasa.gov/users/new).
    Data download can be performed with python's request module
    or with external ARIA2 (needs to be available in PATH) if --aria2 flag is added.

    To query for both MODIS AQUA and TERRA, replace MOD/MYD with M?D.
    Product IDs also accepted in lowercase.
    '''

    parser = argparse.ArgumentParser(description='Query and download MODIS products (Earthdata account required for download)')
    parser.add_argument('product', help='MODIS product ID(s)', nargs='+')
    parser.add_argument('--roi', help='Region of interest. Can be LAT/LON point, bounding box in format (xmin, ymin, xmax, ymax)', nargs='+', required=False)
    parser.add_argument('--tile-filter', help='MODIS tile filter (download only specified tiles)', nargs='+', required=False, metavar='')
    parser.add_argument('-c', '--collection', help='MODIS collection', default='006', metavar='')
    parser.add_argument('-b', '--begin-date', help='Start date (YYYY-MM-DD)', metavar='')
    parser.add_argument('-e', '--end-date', help='End date (YYYY-MM-DD)', metavar='')
    parser.add_argument('--username', help='Earthdata username (required for download)', metavar='')
    parser.add_argument('--password', help='Earthdata password (required for download)', metavar='')
    parser.add_argument('-d', '--targetdir', help='Destination directory', default=os.getcwd(), metavar='')
    parser.add_argument('--download', help='Download data', action='store_true')
    parser.add_argument('--strict-dates', help='Enforce dates strictly', action='store_true')
    parser.add_argument('--multithread', help='Use multiple threads for downloading', action='store_true')

    # fail and print help if no arguments supplied
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    if args.download:

        if not args.username or not args.password:
            raise ValueError("Download was requested, but credentials are missing.")

    # handle targetdir
    target_directory = Path(args.targetdir)
    target_directory.mkdir(exist_ok=True)

    assert target_directory.exists()
    assert target_directory.is_dir()

    products = []

    for product_code in args.product:

        if "M?D" in product_code:
            products.append(product_code.replace("?", "O").upper())
            products.append(product_code.replace("?", "Y").upper())
        else:
            products.append(product_code.upper())

    if args.begin_date:
        try:
            start = datetime.strptime(args.begin_date, "%Y-%m-%d")
        except ValueError:
            print('Error parsing begin-date!')
            raise
    else:
        start = None

    if args.end_date:
        try:
            stop = datetime.strptime(args.end_date, "%Y-%m-%d")
        except ValueError:
            print('Error parsing begin-date!')
            raise
    else:
        stop = None

    print('Running query!')

    query = ModisQuery(
        products=products,
        aoi=args.roi,
        begindate=start,
        enddate=stop,
        tile_filter=args.tile_filter,
        version=args.collection,
    )

    query.search(strict_dates=args.strict_dates)

    print(f'Found {query.nresults} results!')

    if args.download:

        print('Downloading!')

        if query.nresults > 0:

            query.download(
                targetdir=target_directory,
                username=args.username,
                password=args.password,
                multithread=args.multithread
            )

        else:
            print('Nothing to download! Exiting ... ')


    print('Done!')

if __name__ == '__main__':
    main()
