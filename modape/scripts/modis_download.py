#!/usr/bin/env python
"""modis_download.py: Query and download MODIS HDF files."""

from __future__ import absolute_import, division, print_function

import argparse
import datetime
import os
try:
    from pathlib2 import Path
except ImportError:
    from pathlib import Path
import pickle
import re
from subprocess import check_output
import sys
import warnings

import ogr
from modape.modis import ModisQuery
from modape.utils import Credentials

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
    parser.add_argument('--roi', help='Region of interest. Can be LAT/LON point, bounding box in format llx,lly,urx,ury or OGR file (expects epsg:4326! - shp, geojson - convex hull will be used)', nargs='+', required=False)
    parser.add_argument('--tile-filter', help='MODIS tile filter (download only specified tiles)', nargs='+', required=False, metavar='')
    parser.add_argument('-c', '--collection', help='MODIS collection', default=6, metavar='')
    parser.add_argument('-b', '--begin-date', help='Start date (YYYY-MM-DD)', default='2000-01-01', metavar='')
    parser.add_argument('-e', '--end-date', help='End date (YYYY-MM-DD)', default=datetime.date.today().strftime("%Y-%m-%d"), metavar='')
    parser.add_argument('--username', help='Earthdata username (required for download)', metavar='')
    parser.add_argument('--password', help='Earthdata password (required for download)', metavar='')
    parser.add_argument('-d', '--targetdir', help='Destination directory', default=os.getcwd(), metavar='')
    parser.add_argument('--store-credentials', help='Store Earthdata credentials on disk to be used for future downloads (unsecure!)', action='store_true')
    parser.add_argument('--download', help='Download data', action='store_true')
    parser.add_argument('--aria2', help='DEPRACATED! Use ARIA2 for downloading', action='store_true')

    # fail and print help if no arguments supplied
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    credentials = Credentials(args.username, args.password)

    if args.aria2:
        warnings.warn("\nARIA2 is now standard download tool. The flag is therefore depracated and will raise an error in a future release!\n Please consider removing it,", DeprecationWarning)

    if args.download:

        #fail if ARIA2 not installed
        try:
            _ = check_output(['aria2c', '--version'])
        except:
            raise SystemExit('For downloading, ARIA2 to be available in PATH! Please make sure it\'s installed and available in PATH!')


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
    this_dir = Path(__file__).parent
    with open(this_dir.parent.joinpath('data', 'MODIS_V6_PT.pkl').as_posix(), 'rb') as table_raw:
        product_table = pickle.load(table_raw)

    for product in args.product:

        # Handle ? wildcard
        if '?' in product:
            pattern = re.compile(product.replace('?', '.{1,2}') + '$')
            product = [x for x in product_table if re.match(pattern, x)]
        else:
            product = [product] # enable iteration

        for product_subset in product:

            # Load product info from table
            try:
                product_table_subset = product_table[product_subset]
            except KeyError:
                if len(args.product) > 1:
                    print('Product {} not recognized. Skipping ...'.format(product_subset))
                    continue
                else:
                    raise SystemExit('Product {} not recognized!'.format(product_subset))

            # If resolution is bigger than 1km, the product is global
            global_flag = int(product_table_subset['pixel_size']) > 1000

            # If global, select corresponding base URL
            if global_flag:
                if 'MOD' in product_subset:
                    query_url = 'https://e4ftl01.cr.usgs.gov/MOLT/{}.006/'.format(product_subset)
                elif 'MYD' in product_subset:
                    query_url = 'https://e4ftl01.cr.usgs.gov/MOLA/{}.006/'.format(product_subset)
                elif 'MCD' in product_subset:
                    query_url = 'https://e4ftl01.cr.usgs.gov/MOTA/{}.006/'.format(product_subset)
                else:
                    raise SystemExit('Product {} not recognized.'
                                     ' Available: MOD*, MYD*, MCD*')
            else:
                query = []
                # if tile_filter and no ROI, limit spatial query
                if args.tile_filter and not args.roi:
                    # cast to lowercase
                    tiles = [x.lower() for x in args.tile_filter]

                    with open(this_dir.parent.joinpath('data', 'ModlandTiles_bbx.pkl').as_posix(), 'rb') as bbox_raw:
                        bbox = pickle.load(bbox_raw)

                    if len(tiles) == 1:
                        h_indicator, v_indicator = re.findall(r'\d+', tiles[0])
                        bbox_selection = bbox[(bbox.ih == int(h_indicator)) &
                                              (bbox.iv == int(v_indicator))]

                        if bbox_selection.empty:
                            raise SystemExit('No tile with ID {} found. Please check!'.format(tiles[0]))

                        # roi is approx. center point of tile
                        args.roi = [bbox_selection.lat_max.values[0] - 5,
                                    bbox_selection.lon_max.values[0] - (bbox_selection.lon_max.values[0] - bbox_selection.lon_min.values[0])/2]

                    elif len(tiles) > 1:
                        h_indicator = list({re.findall(r'\d+', x.split('v')[0])[0] for x in tiles})
                        v_indicator = list({re.findall(r'\d+', x.split('v')[1])[0] for x in tiles})

                        # aoi is bbox including all tiles from tile_filter, plus 1 degree buffer
                        bbox_selection = bbox.query('|'.join(['ih == {}'.format(int(x)) for x in h_indicator])).query('|'.join(['iv == {}'.format(int(x)) for x in v_indicator]))

                        if bbox_selection.empty:
                            raise SystemExit('No tiles with IDs {} found. Please check!'.format(' '.join(tiles)))

                        args.roi = [min(bbox_selection.lon_min.values)-1,
                                    min(bbox_selection.lat_min.values)-1,
                                    max(bbox_selection.lon_max.values)+1,
                                    max(bbox_selection.lat_max.values)+1]
                    else:
                        pass # not happening

                # Construct query URL
                try:
                    if len(args.roi) == 1 and os.path.isfile(args.roi[0]):
                        try:
                            ds = ogr.Open(args.roi[0])
                            lyr = ds.GetLayer()
                            geometry_collection = ogr.Geometry(ogr.wkbGeometryCollection)

                            for feature in lyr:
                                geometry_collection.AddGeometry(feature.GetGeometryRef())

                            hull = geometry_collection.ConvexHull()
                            geometry = hull.GetGeometryRef(0)
                            point_count = geometry.GetPointCount()
                            coordinates = []

                            for point in range(point_count):
                                lat, lon, _ = geometry.GetPoint(point)
                                coordinates.append('{},{}'.format(lat, lon))

                            query.append('polygon=' + ','.join(coordinates))
                            ds = None
                            lyr = None
                        except:
                            print('\nError reading polygon file. Traceback:\n\n')
                            raise

                    elif len(args.roi) == 2:
                        query.append('latitude={}&longitude={}'.format(*args.roi))
                    elif len(args.roi) == 4:
                        query.append('bbox={},{},{},{}'.format(*args.roi))
                    else:
                        raise ValueError('ROI is expected to be point or bounding box coordinates!')

                except TypeError:
                    raise ValueError('\nDownload of tiled MODIS products requires ROI or tile-filter!\n')

                query.append('version={}'.format(args.collection))
                query.append('date={}/{}'.format(args.begin_date, args.end_date))

                query_url = 'https://lpdaacsvc.cr.usgs.gov/services/inventory?product={}&{}'.format(product_subset,
                                                                                                    '&'.join(query))

            # Run query
            print('\nPRODUCT: {}\n'.format(product_subset))
            query_result = ModisQuery(query_url,
                                      targetdir=args.targetdir,
                                      begindate=args.begin_date,
                                      enddate=args.end_date,
                                      global_flag=global_flag,
                                      tile_filter=args.tile_filter)

            # If download is True and at least one result, download data
            if args.download and query_result.results > 0:
                query_result.set_credentials(credentials.username, credentials.password)
                query_result.download()

if __name__ == '__main__':
    main()
