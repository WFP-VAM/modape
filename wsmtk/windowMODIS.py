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

    if not os.path.exists(args.path):
        raise SystemExit('directory PATH does not exist!')

    if args.product:
        # force product code upper
        args.product = args.product.upper()
    else:
        args.product = ''

    if args.sgrid:
        dset = 'sgrid'
    else:
        dset = 'data'

    if args.roi and len(args.roi) == 4:
        # change order or corner coordinates for MODIStiles
        args.roi = [args.roi[i] for i in [0,3,2,1]]

    this_dir, this_filename = os.path.split(__file__)

    with open(os.path.join(this_dir, "data", "MODIS_V6_PT.pkl"),'rb') as table_raw:
        product_table = pickle.load(table_raw)

    h5files = glob.glob('{}/{}*h5'.format(args.path,args.product))

    if len(h5files) > 0:
        if len(set([re.sub('h\d{2}v\d{2}.','',os.path.basename(x)) for x in h5files])) > 1:
            raise SystemExit("\nMultiple product types found! Please specify and/or check product parameter!\n")
    else:
        raise SystemExit('\nNo products found in specified path (for specified product) - please check input!\n')

    product_ = re.sub('(M\w{6}).+','\\1', os.path.basename(h5files[0]))

    if 'MXD' in product_:
        global_flag = int(product_table[re.sub('MXD','MOD',product_)]['pixel_size']) == 5600
    else:
        global_flag = int(product_table[product_]['pixel_size']) == 5600

    # check if targetdir exists
    if not os.path.exists(args.targetdir):
        print('\nTarget directory {} does not exist! Creating ... '.format(args.targetdir),end='')
        try:
            os.makedirs(args.targetdir)
            print('done.\n')
        except:
            raise

    if not global_flag and args.roi:

        tiles = MODIStiles(args.roi)

        if len(tiles.tiles) is 0:
            raise SystemExit("\nNo MODIS tile(s) found for location. Please check coordinates!")

        tileRX = re.compile('|'.join(tiles.tiles))

        # files for tile result
        h5files = [x for x in h5files if re.search(tileRX,x)]

    # filter for parameter
    if args.parameter:
        h5files = [x for x in h5files if args.parameter in x]

    assert len(h5files) > 0, "\nNo processed MODIS HDF5 files found for combination of product/tile (and parameter)"

    # loop over parameters (could be multiple if unspecified)

    for par in set([re.sub('.+([^\W\d_]{3}).h5','\\1',x) for x in h5files]):

        print('\n')

        h5files_par = [x for x in h5files if par in x]

        # get mosaic
        mosaic = MODISmosaic(files=h5files_par,datemin=args.begin_date,datemax=args.end_date,global_flag=global_flag)

        if args.sgrid:

            filename = '{}/{}{}_sgrid.tif'.format(args.targetdir,args.region.lower(),par.lower())

            print('Processing file {}'.format(filename))

            with mosaic.getRaster(dset,None) as mosaic_ropen:

                try:

                    if args.roi and len(args.roi) > 2:
                        ds = gdal.Warp(filename,mosaic_ropen.raster,
                        dstSRS='EPSG:4326',
                        outputType=mosaic_ropen.dt_gdal[0],
                        xRes=mosaic_ropen.resolution_degrees,
                        yRes=mosaic_ropen.resolution_degrees,
                        srcNodata = mosaic.nodata,
                        dstNodata = mosaic.nodata,
                        outputBounds=(args.roi[0],args.roi[3],args.roi[2],args.roi[1]),
                        resampleAlg='near')

                        ds = None

                    else:

                        ds = gdal.Warp(filename,mosaic_ropen.raster,
                        dstSRS='EPSG:4326',
                        outputType=mosaic_ropen.dt_gdal[0],
                        xRes=mosaic_ropen.resolution_degrees,
                        yRes=mosaic_ropen.resolution_degrees,
                        srcNodata = mosaic.nodata,
                        dstNodata = mosaic.nodata,
                        resampleAlg='near')

                        ds = None

                except Exception as e:
                    print('Error while reading {} data for {}! Please check if dataset exits within file. \n\n Error message:\n\n {}'.format(args.dataset,filename,e))


            del mosaic_ropen

        else:

            # loop over dates

            for ix in mosaic.tempIX:

                filename = '{}/{}{}{}j{}.tif'.format(args.targetdir,args.region.lower(),par.lower(),mosaic.dates[ix][0:4],mosaic.dates[ix][4:7])

                print('Processing file {}'.format(filename))

                with mosaic.getRaster(dset,ix) as mosaic_ropen:

                    try:

                        if args.roi and len(args.roi) > 2:

                            ds = gdal.Warp(filename,mosaic_ropen.raster,
                            dstSRS='EPSG:4326',
                            outputType=mosaic_ropen.dt_gdal[0],
                            xRes=mosaic_ropen.resolution_degrees,
                            yRes=mosaic_ropen.resolution_degrees,
                            srcNodata = mosaic.nodata,
                            dstNodata = mosaic.nodata,
                            outputBounds=(args.roi[0],args.roi[3],args.roi[2],args.roi[1]),
                            resampleAlg='near')

                            ds = None

                        else:

                            ds = gdal.Warp(filename,mosaic_ropen.raster,
                            dstSRS='EPSG:4326',
                            outputType=mosaic_ropen.dt_gdal[0],
                            xRes=mosaic_ropen.resolution_degrees,
                            yRes=mosaic_ropen.resolution_degrees,
                            srcNodata = mosaic.nodata,
                            dstNodata = mosaic.nodata,
                            resampleAlg='near')

                            ds = None

                    except Exception as e:
                        print('Error while reading {} data for {}! Please check if dataset exits within file. \n\n Error message:\n\n {}'.format(args.dataset,filename,e))


                del mosaic_ropen


if __name__=='__main__':
    main()
