#!/usr/bin/env python
from __future__ import print_function
from wsmtk.modis import MODISquery, MODIStiles
import os
import argparse
import datetime
import pickle


def main():

    parser = argparse.ArgumentParser(description="Query and download MODIS products (earthdata accound required for download)")
    parser.add_argument("product", help='MODIS product ID(s)',nargs='+')
    parser.add_argument("--roi", help='Region of interest. Can be LAT/LON point or bounding box in format llx,lly,urx,ury',nargs='+',required=False)
    parser.add_argument("-c","--collection", help='MODIS collection',default=6,metavar='')
    parser.add_argument("-b","--begin-date", help='Start date (YYYY-MM-DD)',default='2000-01-01',metavar='')
    parser.add_argument("-e","--end-date", help='End date (YYYY-MM-DD)',default=datetime.date.today().strftime("%Y-%m-%d"),metavar='')
    parser.add_argument("--username", help='Earthdata username (required for download)',metavar='')
    parser.add_argument("--password", help='Earthdata password (required for download)',metavar='')
    parser.add_argument("-d","--dest", help='Destination directory',default=os.getcwd(),metavar='')
    parser.add_argument("-v","--verbose", help='Destination directory',action='store_true')
    parser.add_argument("--download", help='Download data',action='store_true')
    args = parser.parse_args()

    '''
    if args.roi:
        if len(args.roi) is 2:
            tiles = MODIStiles(args.roi).tiles
        elif len(args.roi) is 4:
            tiles = MODIStiles([args.roi[i] for i in [0,3,2,1]]).tiles
        else:
            print('Error parsing ROI - querying globally.')
            tiles = None
    else:
        tiles = None

    '''

    if args.download & (not args.username or not args.password):
        raise SystemExit('Downloading requires username and password!')

    this_dir, this_filename = os.path.split(__file__)

    with open(os.path.join(this_dir, "data", "MODIS_V6_PT.pkl")) as table_raw:
        product_table = pickle.load(table_raw)

    for p in args.product:

        try:
            product_table_sub = product_table[p]

        except KeyError:

            if len(args.product) > 1:
                print('Product {} not recognized. Skipping ...'.format(p))
                continue
            else:
                raise SystemExit('Product {} not recognized!'.format(p))

        global_flag = int(product_table_sub['pixel_size']) > 1000

        if global_flag:

            if 'MOD' in p:
                queryURL = 'https://e4ftl01.cr.usgs.gov/MOLT/{}.006/'.format(p)
            elif 'MYD' in p:
                queryURL = 'https://e4ftl01.cr.usgs.gov/MOLA/{}.006/'.format(p)
            elif 'MCD' in p:
                queryURL = 'https://e4ftl01.cr.usgs.gov/MOTA/{}.006/'.format(p)
            else:
                raise SystemExit('Product {} not recognized. Available: MOD*, MYD*, MCD*')


        else:

            query = []

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

            queryURL = 'https://lpdaacsvc.cr.usgs.gov/services/inventory?product={}&{}'.format(p,'&'.join(query))


        print('\nPRODUCT: {}\n'.format(p))

        res = MODISquery(queryURL,rawdir=args.dest,begindate=args.begin_date,enddate=args.end_date,global_flag=global_flag)

        if args.download and res.results > 0:
            res.setCredentials(args.username,args.password)
            res.download()


if __name__ == '__main__':
    main()
