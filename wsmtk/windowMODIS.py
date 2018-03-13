#!/usr/bin/env python
from __future__ import print_function
from wsmtk.modis import MODISquery
import os
import gdal
import argparse
import datetime


def main():

    # parser = argparse.ArgumentParser(description="Extract a window from MODIS products")
    # parser.add_argument("product", help='MODIS product ID(s)',nargs='+')
    # parser.add_argument("--roi", help='Region of interest. Can be LAT/LON point or bounding box in format llx,lly,urx,ury',nargs='+',required=True)
    # parser.add_argument("--prcdir", help='Storage directory for PROCESSED MODIS files',default=os.getcwd(),metavar='')
    # parser.add_argument("--region", help='region 3 letter region code (default is "reg")',default='reg',metavar='')
    # parser.add_argument("-b","--begin-date", help='Start date (YYY-MM-DD)',metavar='')
    # parser.add_argument("-e","--end-date", help='End date (YYY-MM-DD)',metavar='')
    # parser.add_argument("--parameter", help='VAM parameter code',metavar='')
    # parser.add_argument("--dataset", help='Dataset to extract (either Raw or Smoothed [DEFAULT = Smoothed])',metavar='')
    # parser.add_argument("--targetdir", help='Target directory for GeoTIFFs (default current directory)',default=os.getcwd(),metavar='')
    # args = parser.parse_args()

    pass


if __name__=='__main__':
    main()
