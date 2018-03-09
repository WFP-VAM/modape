#!/usr/bin/env python
from __future__ import print_function
from wsmtk.modis import MODISquery
import os
import argparse
import datetime


def main():

    parser = argparse.ArgumentParser(description="Query and download MODIS products (earthdata accound required for download)")
    parser.add_argument("product", help='MODIS product ID(s)',nargs='+')
    parser.add_argument("--roi", help='Region of interest. Can be LAT/LON point or bounding box in format llx,lly,urx,ury',nargs='+',required=True)
    parser.add_argument("-c","--collection", help='MODIS collection',default=6,metavar='')
    parser.add_argument("-b","--begin-date", help='Start date (YYY-MM-DD)',default='2000-01-01',metavar='')
    parser.add_argument("-e","--end-date", help='End date (YYY-MM-DD)',default=datetime.date.today().strftime("%Y-%m-%d"),metavar='')
    parser.add_argument("--username", help='Earthdata username (required for download)',metavar='')
    parser.add_argument("--password", help='Earthdata password (required for download)',metavar='')
    parser.add_argument("-d","--dest", help='Destindation directory',default=os.getcwd(),metavar='')
    parser.add_argument("-v","--verbose", help='Destindation directory',action='store_true')
    parser.add_argument("--download", help='Download data',action='store_true')
    args = parser.parse_args()

    query = []


    if len(args.roi) is 2:
        query.append('latitude={}&longitude={}'.format(*args.roi))
    elif len(args.roi) is 4:
        query.append('bbox={},{},{},{}'.format(*args.roi))
    else:
        raise SystemExit('ROI is expected to be point or bounding box coordinates!')

    if args.download & (not args.username or not args.password):
        raise SystemExit('Downloading requires username and password!')

    query.append('version={}'.format(args.collection))
    query.append('date={}/{}'.format(args.begin_date,args.end_date))

    for p in args.product:

        queryURL = 'https://lpdaacsvc.cr.usgs.gov/services/inventory?product={}&{}'.format(p,'&'.join(query))

        res = MODISquery(queryURL,rawdir=args.dest)

        if args.download:
            res.setCredentials(args.username,args.password)
            res.download()


if __name__ == '__main__':
    main()
