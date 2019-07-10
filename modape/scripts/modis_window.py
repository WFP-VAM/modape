#!/usr/bin/env python
# pylint: disable=broad-except
"""modis_window.py: Create mosaics from smooth MODIS HDF5 files and
   save them as GeoTIFFs.
"""

from __future__ import absolute_import, division, print_function

import argparse
import datetime
import os
from pathlib import Path
import pickle
import re
import sys

from modape.modis import modis_tiles, ModisMosaic

try:
    import gdal
except ImportError:
    from osgeo import gdal

def main():
    """Create mosaics (or subsets) from smoothed MODIS HDF5 files.

    Given an ROI and smoothed MODIS files in path, either a mosaic or a subset of smoothed MODIS file(s)
    is created. Depending on the begin-date and end-date parameters, timesteps within the range are written
    to disk as GeoTIFF files.

    If ROI is only a point location, the entire intersecting tile or the entire global file will be returned.
    """

    parser = argparse.ArgumentParser(description='Extract a window from MODIS products')
    parser.add_argument('path', help='Path to processed MODIS h5 files')
    parser.add_argument('-p', '--product', help='MODIS product ID (can be parial match with *)', metavar='')
    parser.add_argument('--roi', help='Region of interest. Can be LAT/LON point or bounding box in format llx lly urx ury', nargs='+', type=float)
    parser.add_argument('--region', help='region 3 letter region code (default is \'reg\')', default='reg', metavar='')
    parser.add_argument('-b', '--begin-date', help='Start date (YYYYMM)', default=datetime.date(2000, 1, 1).strftime('%Y%m'), metavar='')
    parser.add_argument('-e', '--end-date', help='End date (YYYYMM)', default=datetime.date.today().strftime('%Y%m'), metavar='')
    parser.add_argument('--vampc', help='VAM product code', metavar='')
    parser.add_argument('-d', '--targetdir', help='Target directory for GeoTIFFs (default current directory)', default=os.getcwd(), metavar='')
    parser.add_argument('--sgrid', help='Extract (mosaic of) s value grid(s)', action='store_true')
    parser.add_argument('--force-doy', help='Force filenaming with DOY for 5 & 10 day data', action='store_true')

    # fail and print help if no arguments supplied
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)

    args = parser.parse_args()

    input_dir = Path(args.path)
    output_dir = Path(args.targetdir)

    if not input_dir.exists():
        raise SystemExit('directory PATH does not exist!')

    # Check if targetdir exists, create if not
    if not output_dir.exists():
        print('\nTarget directory {} does not exist! Creating ... '.format(output_dir.as_posix()), end='')
        output_dir.mkdir(parents=True)
        print('done.\n')

    if args.product:
        args.product = args.product.upper() # force product code upper
    else:
        args.product = ''

    # Select dataset
    if args.sgrid:
        dset = 'sgrid'
    else:
        dset = 'data'

    # If ROI is a bounding box, change order or corner coordinates for modis_tiles
    if args.roi and len(args.roi) == 4:
        args.roi = [args.roi[i] for i in [0, 3, 2, 1]]

    # Load product table
    this_dir, _ = os.path.split(__file__)
    package_dir = os.path.abspath(os.path.join(this_dir, os.pardir))

    with open(os.path.join(package_dir, 'data', 'MODIS_V6_PT.pkl'), 'rb') as table_raw:
        product_table = pickle.load(table_raw)

    # List HDF5 files in path

    h5files = list(input_dir.glob(args.product + '*h5'))

    # Make sure only compatible files are used for the mosaic
    if not h5files:
        raise ValueError('\nNo products found in specified path (for specified product) - please check input!\n')

    #if len({re.sub(r'h\d{2}v\d{2}.', '', x.name) for x in h5files}) > 1:
    #raise ValueError('\nMultiple product types found! Please specify and/or check product parameter!\n')

    groups = {re.sub(r'h\d{2}v\d{2}.', '*', x.name) for x in h5files}

    for grp in groups:

        h5files_select = [x for x in h5files if re.match(grp, x.name)]

        # Extract product ID from referece
        product_ = re.sub(r'(M\w{6}).+', '\\1', h5files_select[0].name)

        print('\nProcessing product {}'.format(product_))

        # Check if product is global
        if 'MXD' in product_:
            global_flag = int(product_table[re.sub('MXD', 'MOD', product_)]['pixel_size']) == 5600
        else:
            global_flag = int(product_table[product_]['pixel_size']) == 5600

        # If the product is not global and there's an ROI, we need to query the intersecting tiles
        if not global_flag and args.roi:
            tiles = modis_tiles(args.roi)
            if not tiles:
                raise ValueError('\nNo MODIS tile(s) found for location. Please check coordinates!')

            # Regexp for tiles
            tile_regexp = re.compile('|'.join(tiles))

            # Files for tile result
            h5files_select = [x for x in h5files_select if re.search(tile_regexp, x.name)]

        # Filter for product code
        if args.vampc:
            h5files_select = [x for x in h5files_select if args.vampc in x.name]

        # Assert that there are results
        assert h5files_select, '\nNo processed MODIS HDF5 files found for combination of product/tile (and VAM product code)'

        # Iterate over VPCs (could be multiple if unspecified)
        for vam_code in {re.sub(r'.+([^\W\d_]{3}).h5', '\\1', x.name) for x in h5files_select}:
            print('\n')

            h5files_select_vpc = [x.as_posix() for x in h5files_select if vam_code in x.name] # Subset files for vam_code

            # Get mosaic
            mosaic = ModisMosaic(files=h5files_select_vpc,
                                 datemin=args.begin_date,
                                 datemax=args.end_date,
                                 global_flag=global_flag)

            # Extract s-grid if True
            if args.sgrid:
                filename = output_dir.joinpath(args.region.lower() + vam_code.lower() + '_sgrid.tif').as_posix()

                print('Processing file {}'.format(filename))
                with mosaic.get_raster(dset, None) as mosaic_ropen:
                    # Subset if bbox was supplied
                    try:
                        if args.roi and len(args.roi) > 2:
                            wopt = gdal.WarpOptions(
                                dstSRS='EPSG:4326',
                                outputType=mosaic_ropen.dt_gdal[0],
                                xRes=mosaic_ropen.resolution_degrees,
                                yRes=mosaic_ropen.resolution_degrees,
                                srcNodata=mosaic.nodata,
                                dstNodata=mosaic.nodata,
                                outputBounds=(args.roi[0], args.roi[3], args.roi[2], args.roi[1]),
                                resampleAlg='near',
                                multithread=True,
                                creationOptions=['COMPRESS=LZW', 'PREDICTOR=2'],
                            )

                            _ = gdal.Warp(
                                filename,
                                mosaic_ropen.raster,
                                options=wopt,
                            )

                            _ = None
                        else:
                            wopt = gdal.WarpOptions(
                                dstSRS='EPSG:4326',
                                outputType=mosaic_ropen.dt_gdal[0],
                                xRes=mosaic_ropen.resolution_degrees,
                                yRes=mosaic_ropen.resolution_degrees,
                                srcNodata=mosaic.nodata,
                                dstNodata=mosaic.nodata,
                                resampleAlg='near',
                                multithread=True,
                                creationOptions=['COMPRESS=LZW', 'PREDICTOR=2'],
                            )

                            _ = gdal.Warp(
                                filename,
                                mosaic_ropen.raster,
                                options=wopt,
                            )

                            _ = None
                    except Exception as e:
                        print('Error while reading {} data for {}! Please check if dataset exits within file. \n\n Error message:\n\n {}'.format(args.dataset, filename, e))
                del mosaic_ropen
            else:
                # If the dataset is not s-grid, we need to iterate over dates
                for ix in mosaic.temp_index:

                    if mosaic.labels and not args.force_doy:

                        filename = output_dir.joinpath(args.region.lower() + vam_code.lower() + mosaic.labels[ix] + '.tif').as_posix()

                    else:

                        filename = output_dir.joinpath(args.region.lower() + vam_code.lower() + mosaic.dates[ix][0:4] + 'j' + mosaic.dates[ix][4:7] + '.tif').as_posix()


                    with mosaic.get_raster(dset, ix) as mosaic_ropen:
                        try:
                            if args.roi and len(args.roi) > 2:
                                wopt = gdal.WarpOptions(
                                    dstSRS='EPSG:4326',
                                    outputType=mosaic_ropen.dt_gdal[0],
                                    xRes=mosaic_ropen.resolution_degrees,
                                    yRes=mosaic_ropen.resolution_degrees,
                                    srcNodata=mosaic.nodata,
                                    dstNodata=mosaic.nodata,
                                    outputBounds=(args.roi[0], args.roi[3], args.roi[2], args.roi[1]),
                                    resampleAlg='near',
                                    multithread=True,
                                    creationOptions=['COMPRESS=LZW', 'PREDICTOR=2']
                                )

                                _ = gdal.Warp(
                                    filename,
                                    mosaic_ropen.raster,
                                    options=wopt,
                                )

                                _ = None
                            else:
                                wopt = gdal.WarpOptions(
                                    dstSRS='EPSG:4326',
                                    outputType=mosaic_ropen.dt_gdal[0],
                                    xRes=mosaic_ropen.resolution_degrees,
                                    yRes=mosaic_ropen.resolution_degrees,
                                    srcNodata=mosaic.nodata,
                                    dstNodata=mosaic.nodata,
                                    resampleAlg='near',
                                    multithread=True,
                                    creationOptions=['COMPRESS=LZW', 'PREDICTOR=2'],
                                )

                                _ = gdal.Warp(
                                    filename,
                                    mosaic_ropen.raster,
                                    options=wopt,
                                )

                                _ = None
                        except Exception as e:
                            print('Error while reading {} data for {}! Please check if dataset exits within file. \n\n Error message:\n\n {}'.format(args.dataset, filename, e))
                    del mosaic_ropen

if __name__ == '__main__':
    main()
