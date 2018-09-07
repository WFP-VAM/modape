#!/usr/bin/env python
from __future__ import print_function
from wsmtk.modis import MODIStiles, MODISmosaic
import glob
import re
import os
import sys
import gdal
import argparse
import datetime
import numpy as np
import pickle


def main():
    '''Create mosaics (or subsets) from smoothed MODIS HDF5 files.

    Given an ROI and smoothed MODIS files in path, either a mosaic or a subset of smoothed MODIS file(s)
    is created. Depending on the begin-date and end-date parameters, timesteps within the range are written
    to disk as GeoTIFF files.

    If ROI is only a point location, the entire intersecting tile or the entire global file will be returned.
    '''

    parser = argparse.ArgumentParser(description="Extract a window from MODIS products")
    parser.add_argument("path", help='Path to processed MODIS h5 files')
    parser.add_argument("-p","--product", help='MODIS product ID (can be parial match with *)',metavar='')
    parser.add_argument("--roi", help='Region of interest. Can be LAT/LON point or bounding box in format llx lly urx ury',nargs='+',type=float)
    parser.add_argument("--region", help='region 3 letter region code (default is "reg")',default='reg',metavar='')
    parser.add_argument("-b","--begin-date", help='Start date (YYYYMM)',default=datetime.date(2000,1,1).strftime("%Y%m"),metavar='')
    parser.add_argument("-e","--end-date", help='End date (YYYYMM)',default=datetime.date.today().strftime("%Y%m"),metavar='')
    parser.add_argument("--parameter", help='VAM parameter code',metavar='')
    parser.add_argument("-d","--targetdir", help='Target directory for GeoTIFFs (default current directory)',default=os.getcwd(),metavar='')
    parser.add_argument("--sgrid", help='Extract (mosaic of) s value grid(s))',action='store_true')

    # fail and print help if no arguments supplied
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    # Check if path exists
    if not os.path.exists(args.path):
        raise SystemExit('directory PATH does not exist!')


    # If supplied, force product code upper
    if args.product:
        args.product = args.product.upper()
    else:
        args.product = ''

    # Select dataset
    if args.sgrid:
        dset = 'sgrid'
    else:
        dset = 'data'


    # If ROI is a bounding box, change order or corner coordinates for MODIStiles
    if args.roi and len(args.roi) == 4:
        args.roi = [args.roi[i] for i in [0,3,2,1]]

    # Load product table
    this_dir, this_filename = os.path.split(__file__)

    with open(os.path.join(this_dir, "data", "MODIS_V6_PT.pkl"),'rb') as table_raw:
        product_table = pickle.load(table_raw)

    # List HDF5 files in path
    h5files = glob.glob('{}/{}*h5'.format(args.path,args.product))

    # Make sure only compatible files are used for the mosaic
    if len(h5files) > 0:
        if len(set([re.sub('h\d{2}v\d{2}.','',os.path.basename(x)) for x in h5files])) > 1:
            raise SystemExit("\nMultiple product types found! Please specify and/or check product parameter!\n")
    else:
        raise SystemExit('\nNo products found in specified path (for specified product) - please check input!\n')

    # Extract product ID from referece
    product_ = re.sub('(M\w{6}).+','\\1', os.path.basename(h5files[0]))

    # Check if product is global
    if 'MXD' in product_:
        global_flag = int(product_table[re.sub('MXD','MOD',product_)]['pixel_size']) == 5600
    else:
        global_flag = int(product_table[product_]['pixel_size']) == 5600

    # Check if targetdir exists, create if not
    if not os.path.exists(args.targetdir):
        print('\nTarget directory {} does not exist! Creating ... '.format(args.targetdir),end='')
        try:
            os.makedirs(args.targetdir)
            print('done.\n')
        except:
            raise

    # If the product is not global and there's an ROI, we need to query the intersecting tiles
    if not global_flag and args.roi:

        tiles = MODIStiles(args.roi)

        if len(tiles.tiles) is 0:
            raise SystemExit("\nNo MODIS tile(s) found for location. Please check coordinates!")

        # Regexp for tiles
        tileRX = re.compile('|'.join(tiles.tiles))

        # Files for tile result
        h5files = [x for x in h5files if re.search(tileRX,x)]

    # Filter for parameter
    if args.parameter:
        h5files = [x for x in h5files if args.parameter in x]

    # Assert that there are results
    assert len(h5files) > 0, "\nNo processed MODIS HDF5 files found for combination of product/tile (and parameter)"

    # Iterate over parameters (could be multiple if unspecified)
    for par in set([re.sub('.+([^\W\d_]{3}).h5','\\1',x) for x in h5files]):

        print('\n')

        # Subset files for parameter
        h5files_par = [x for x in h5files if par in x]

        # Get mosaic
        mosaic = MODISmosaic(files=h5files_par,datemin=args.begin_date,datemax=args.end_date,global_flag=global_flag)

        # Extract s-grid if True
        if args.sgrid:

            filename = '{}/{}{}_sgrid.tif'.format(args.targetdir,args.region.lower(),par.lower())

            print('Processing file {}'.format(filename))

            with mosaic.getRaster(dset,None) as mosaic_ropen:

                # Subset if bbox was supplied
                try:
                    if args.roi and len(args.roi) > 2:

                        wopt = gdal.WarpOptions(
                        dstSRS='EPSG:4326',
                        outputType=mosaic_ropen.dt_gdal[0],
                        xRes = mosaic_ropen.resolution_degrees,
                        yRes = mosaic_ropen.resolution_degrees,
                        srcNodata=mosaic.nodata,
                        dstNodata=mosaic.nodata,
                        outputBounds=(args.roi[0],args.roi[3],args.roi[2],args.roi[1]),
                        resampleAlg='near',
                        multithread=True,
                        creationOptions=['COMPRESS=LZW','PREDICTOR=2'])

                        ds = gdal.Warp(filename,mosaic_ropen.raster,
                        options=wopt)

                        ds = None

                    else:

                        wopt = gdal.WarpOptions(
                        dstSRS='EPSG:4326',
                        outputType=mosaic_ropen.dt_gdal[0],
                        xRes = mosaic_ropen.resolution_degrees,
                        yRes = mosaic_ropen.resolution_degrees,
                        srcNodata=mosaic.nodata,
                        dstNodata=mosaic.nodata,
                        resampleAlg='near',
                        multithread=True,
                        creationOptions=['COMPRESS=LZW','PREDICTOR=2'])

                        ds = gdal.Warp(filename,mosaic_ropen.raster,
                        options=wopt)

                        ds = None

                except Exception as e:
                    print('Error while reading {} data for {}! Please check if dataset exits within file. \n\n Error message:\n\n {}'.format(args.dataset,filename,e))


            del mosaic_ropen

        else:


            # If the dataset is not s-grid, we need to iterate over dates
            for ix in mosaic.tempIX:

                filename = '{}/{}{}{}j{}.tif'.format(args.targetdir,args.region.lower(),par.lower(),mosaic.dates[ix][0:4],mosaic.dates[ix][4:7])

                print('Processing file {}'.format(filename))

                with mosaic.getRaster(dset,ix) as mosaic_ropen:

                    try:

                        if args.roi and len(args.roi) > 2:

                            wopt = gdal.WarpOptions(
                            dstSRS='EPSG:4326',
                            outputType=mosaic_ropen.dt_gdal[0],
                            xRes = mosaic_ropen.resolution_degrees,
                            yRes = mosaic_ropen.resolution_degrees,
                            srcNodata=mosaic.nodata,
                            dstNodata=mosaic.nodata,
                            outputBounds=(args.roi[0],args.roi[3],args.roi[2],args.roi[1]),
                            resampleAlg='near',
                            multithread=True,
                            creationOptions=['COMPRESS=LZW','PREDICTOR=2'])

                            ds = gdal.Warp(filename,mosaic_ropen.raster,
                            options=wopt)

                            ds = None

                        else:

                            wopt = gdal.WarpOptions(
                            dstSRS='EPSG:4326',
                            outputType=mosaic_ropen.dt_gdal[0],
                            xRes = mosaic_ropen.resolution_degrees,
                            yRes = mosaic_ropen.resolution_degrees,
                            srcNodata=mosaic.nodata,
                            dstNodata=mosaic.nodata,
                            resampleAlg='near',
                            multithread=True,
                            creationOptions=['COMPRESS=LZW','PREDICTOR=2'])

                            ds = gdal.Warp(filename,mosaic_ropen.raster,
                            options=wopt)

                            ds = None

                    except Exception as e:
                        print('Error while reading {} data for {}! Please check if dataset exits within file. \n\n Error message:\n\n {}'.format(args.dataset,filename,e))


                del mosaic_ropen


if __name__=='__main__':
    main()
