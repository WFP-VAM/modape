#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import datetime
import os
import pickle
import re
import sys

from modape.modis import MODISquery
from modape.utils import Credentials, pload
import ogr

try:
    range = xrange
except NameError:
    pass

def main():
    '''Query and download MODIS products.

    This function allows for querying and downloading MODIS products in bulk. Multiple products can be queried and downloaded with one
    function call. For downloading data, valid earthdata credentials are required (to register, visit https://urs.earthdata.nasa.gov/users/new).
    Data download can be performed with python's request module or with external ARIA2 (needs to be available in PATH) if --aria2 flag is added.

    To query for both MODIS AQUA and TERRA, replace MOD/MYD with M?D. Product IDs also accepted in lowercase.
    '''

    parser = argparse.ArgumentParser(description="Query and download MODIS products (Earthdata account required for download)")
    parser.add_argument("product", help='MODIS product ID(s)',nargs='+')
    parser.add_argument("--roi", help='Region of interest. Can be LAT/LON point, bounding box in format llx,lly,urx,ury or OGR file (shp, geojson - convex hull will be used)',nargs='+',required=False)
    parser.add_argument("--tile-filter", help='MODIS tile filter (download only specified tiles)',nargs='+',required=False,metavar='')
    parser.add_argument("-c","--collection", help='MODIS collection',default=6,metavar='')
    parser.add_argument("-b","--begin-date", help='Start date (YYYY-MM-DD)',default='2000-01-01',metavar='')
    parser.add_argument("-e","--end-date", help='End date (YYYY-MM-DD)',default=datetime.date.today().strftime("%Y-%m-%d"),metavar='')
    parser.add_argument("--username", help='Earthdata username (required for download)',metavar='')
    parser.add_argument("--password", help='Earthdata password (required for download)',metavar='')
    parser.add_argument("-d","--targetdir", help='Destination directory',default=os.getcwd(),metavar='')
    parser.add_argument("--store-credentials", help='Store Earthdata credentials on disk to be used for future downloads (unsecure!)',action='store_true')
    #parser.add_argument("-v","--verbose", help='Verbosity',action='store_true')
    parser.add_argument("--download", help='Download data',action='store_true')
    parser.add_argument("--aria2", help='Use ARIA2 for downloading',action='store_true')

    # fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    credentials = Credentials(args.username, args.password)

    # Check for credentials if download is True
    if args.download & (not credentials.complete):

        try:

            credentials.retrieve()

        except:

            raise SystemExit('\nError: Earthdata credentials not found!\n')

    elif args.store_credentials:

        credentials.store()

    else:

        pass

    args.product = [x.upper() for x in args.product]

    # Load product table
    this_dir, this_filename = os.path.split(__file__)
    product_table = pload(os.path.join(this_dir, "data", "MODIS_V6_PT.pkl"))

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

                # if tile_filter and no ROI, limit spatial query

                if args.tile_filter and not args.roi:

                    # cast to lowercase
                    tiles = [x.lower() for x in args.tile_filter]

                    with open(os.path.join(this_dir, "data", "ModlandTiles_bbx.pkl"),'rb') as bbx_raw:
                        bbx = pickle.load(bbx_raw)


                    if len(tiles) == 1:

                        h,v = re.findall('\d+',tiles[0])

                        bbx_sel = bbx[(bbx.ih == int(h)) & (bbx.iv == int(v))]

                        # roi is approx. center point of tile
                        args.roi = [bbx_sel.lat_max.values[0] - 5, bbx_sel.lon_max.values[0] - (bbx_sel.lon_max.values[0] - bbx_sel.lon_min.values[0])/2]

                    elif len(tiles) > 1:

                        h = list(set([re.findall('\d+',x.split('v')[0])[0] for x in tiles]))
                        v = list(set([re.findall('\d+',x.split('v')[1])[0] for x in tiles]))

                        # aoi is bbox including all tiles from tile_filter, plus 1 degree buffer
                        bbx_sel = bbx.query('|'.join(['ih == {}'.format(int(x)) for x in h])).query('|'.join(['iv == {}'.format(int(x)) for x in v]))

                        args.roi = [min(bbx_sel.lon_min.values)-1, min(bbx_sel.lat_min.values)-1, max(bbx_sel.lon_max.values)+1,max(bbx_sel.lat_max.values)+1]

                    else:

                        # not happening
                        pass

                # Construct query URL

                try:

                    if len(args.roi) == 1 and os.path.isfile(args.roi[0]):

                        try:

                            ds = ogr.Open(args.roi[0])
                            lyr = ds.GetLayer()

                            geomcol = ogr.Geometry(ogr.wkbGeometryCollection)

                            for feature in lyr:
                                geomcol.AddGeometry(feature.GetGeometryRef())

                            hull = geomcol.ConvexHull()
                            geom = hull.GetGeometryRef(0)
                            pointcount = geom.GetPointCount()

                            crds = []

                            for pt in range(pointcount):
                                lat, lon, z = geom.GetPoint(pt)
                                crds.append('{},{}'.format(lon,lat))

                            query.append('polygon=' + ','.join(crds))

                            ds = None
                            lyr = None

                        except:
                            print('\nError reading polygon file. Traceback:\n\n')
                            raise


                    elif len(args.roi) is 2:
                        query.append('latitude={}&longitude={}'.format(*args.roi))
                    elif len(args.roi) is 4:
                        query.append('bbox={},{},{},{}'.format(*args.roi))
                    else:
                        raise SystemExit('ROI is expected to be point or bounding box coordinates!')

                except TypeError:
                    #query.append('bbox=-180,-90,180,90')
                    raise SystemExit('\nDownload of tiled MODIS products requires ROI or tile-filter!\n')


                query.append('version={}'.format(args.collection))
                query.append('date={}/{}'.format(args.begin_date,args.end_date))

                queryURL = 'https://lpdaacsvc.cr.usgs.gov/services/inventory?product={}&{}'.format(p2,'&'.join(query))


            # Run query

            print('\nPRODUCT: {}\n'.format(p2))

            res = MODISquery(queryURL,targetdir=args.targetdir,begindate=args.begin_date,enddate=args.end_date,global_flag=global_flag,aria2=args.aria2, tile_filter = args.tile_filter)

            # If download is True and at least one result, download data
            if args.download and res.results > 0:
                res.setCredentials(credentials.username,credentials.password)
                res.download()


if __name__ == '__main__':
    main()
