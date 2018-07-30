#!/usr/bin/env python
from __future__ import print_function
from wsmtk.modis import MODISquery
import os
import sys
import argparse
import datetime
import pickle
import re

def main():
    '''Query and download MODIS products.

    This function allows for querying and downloading MODIS products in bulk. Multiple products can be queried and downloaded with one
    function call. For downloading data, valid earthdata credentials are required (to register, visit https://urs.earthdata.nasa.gov/users/new).
    Data download can be performed with python's request module or with external WGET (needs to be available in PATH) if --wget flag is added.

    To query for both MODIS AQUA and TERRA, replace MOD/MYD with M?D. Product IDs also accepted in lowercase.
    '''

    parser = argparse.ArgumentParser(description="Query and download MODIS products (earthdata accound required for download)")
    parser.add_argument("product", help='MODIS product ID(s)',nargs='+')
    parser.add_argument("--roi", help='Region of interest. Can be LAT/LON point or bounding box in format llx,lly,urx,ury',nargs='+',required=False)
    parser.add_argument("-c","--collection", help='MODIS collection',default=6,metavar='')
    parser.add_argument("-b","--begin-date", help='Start date (YYYY-MM-DD)',default='2000-01-01',metavar='')
    parser.add_argument("-e","--end-date", help='End date (YYYY-MM-DD)',default=datetime.date.today().strftime("%Y-%m-%d"),metavar='')
    parser.add_argument("--username", help='Earthdata username (required for download)',metavar='')
    parser.add_argument("--password", help='Earthdata password (required for download)',metavar='')
    parser.add_argument("-d","--dest", help='Destination directory',default=os.getcwd(),metavar='')
    parser.add_argument("-v","--verbose", help='Verbosity',action='store_true')
    parser.add_argument("--download", help='Download data',action='store_true')
    parser.add_argument("--wget", help='Use WGET for downloading',action='store_true')

    # fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Check for credentials if download is True
    if args.download & (not args.username or not args.password):
        raise SystemExit('Downloading requires username and password!')

    args.product = [x.upper() for x in args.product]

    # Load product table
    this_dir, this_filename = os.path.split(__file__)
    with open(os.path.join(this_dir, "data", "MODIS_V6_PT.pkl"),'rb') as table_raw:
        product_table = pickle.load(table_raw)

    for p in args.product:

        # Handle ? wildcard
        if '?' in p:
            patt = re.compile(p.replace('?','.{1,2}') + '$')
            p = [x for x in product_table if re.match(patt,x)]
        else:
            p = [p]

        for p2 in p:

            # Load product info from table
            try:
                product_table_sub = product_table[p2]

            except KeyError:

                if len(args.product) > 1:
                    print('Product {} not recognized. Skipping ...'.format(p2))
                    continue
                else:
                    raise SystemExit('Product {} not recognized!'.format(p2))

            # If resolution is bigger than 1km, the product is global
            global_flag = int(product_table_sub['pixel_size']) > 1000

            # If global, select corresponding base URL
            if global_flag:

                if 'MOD' in p2:
                    queryURL = 'https://e4ftl01.cr.usgs.gov/MOLT/{}.006/'.format(p2)
                elif 'MYD' in p2:
                    queryURL = 'https://e4ftl01.cr.usgs.gov/MOLA/{}.006/'.format(p2)
                elif 'MCD' in p2:
                    queryURL = 'https://e4ftl01.cr.usgs.gov/MOTA/{}.006/'.format(p2)
                else:
                    raise SystemExit('Product {} not recognized. Available: MOD*, MYD*, MCD*')

            else:

                query = []

                # Construct query URL
                try:

                    if len(args.roi) is 2:
                        query.append('latitude={}&longitude={}'.format(*args.roi))
                    elif len(args.roi) is 4:
                        query.append('bbox={},{},{},{}'.format(*args.roi))
                    else:
                        raise SystemExit('ROI is expected to be point or bounding box coordinates!')

                except TypeError:
                    raise SystemExit('Download of tiled MODIS products requires ROI!')


                query.append('version={}'.format(args.collection))
                query.append('date={}/{}'.format(args.begin_date,args.end_date))

                queryURL = 'https://lpdaacsvc.cr.usgs.gov/services/inventory?product={}&{}'.format(p2,'&'.join(query))


            # Run query
            print('\nPRODUCT: {}\n'.format(p2))

            res = MODISquery(queryURL,rawdir=args.dest,begindate=args.begin_date,enddate=args.end_date,global_flag=global_flag,wget=args.wget)

            # If download is True and at least one result, download data
            if args.download and res.results > 0:
                res.setCredentials(args.username,args.password)
                res.download()


if __name__ == '__main__':
    main()
