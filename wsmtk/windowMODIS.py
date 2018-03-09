#!/usr/bin/env python
from __future__ import print_function
from wsmtk.modis import MODISquery
import os
import gdal
import argparse
import datetime


def main():

    ##test!

    this_dir, this_filename = os.path.split(__file__)

    DATA_PATH = os.path.join(this_dir, "data", "MODIS_TILES.tif")
    print(DATA_PATH)

    ds = gdal.Open(DATA_PATH)

    print(ds)

    ds = None


if __name__=='__main__':
    main()
