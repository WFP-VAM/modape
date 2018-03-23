#!/usr/bin/env python
from __future__ import print_function
from wsmtk.modis import MODIStiles, MODISwindow
import os
import gdal
import argparse
import datetime
import numpy as np
from .utils import aoi2ix
from .hdf5 import h5_file, h5_readArr

def main():

    parser = argparse.ArgumentParser(description="Extract a window from MODIS products")
    parser.add_argument("product", help='MODIS product ID')
    parser.add_argument("--roi", help='Region of interest. Can be LAT/LON point or bounding box in format llx,lly,urx,ury',nargs='+',required=True)
    parser.add_argument("--prcdir", help='Storage directory for PROCESSED MODIS files',default=os.getcwd(),metavar='')
    parser.add_argument("--region", help='region 3 letter region code (default is "reg")',default='reg',metavar='')
    parser.add_argument("-b","--begin-date", help='Start date (YYYYMM)',default=datetime.date(2000,1,1).strftime("%Y%m"),metavar='')
    parser.add_argument("-e","--end-date", help='End date (YYYYMM)',default=datetime.date.today().strftime("%Y%m"),metavar='')
    parser.add_argument("--parameter", help='VAM parameter code',metavar='')
    parser.add_argument("--dataset", help='Dataset to extract (either Raw or Smoothed [DEFAULT = Smoothed])',default='Smoothed',metavar='')
    parser.add_argument("--targetdir", help='Target directory for GeoTIFFs (default current directory)',default=os.getcwd(),metavar='')
    args = parser.parse_args()

    # change order or corner coordinates for MODIStiles
    if len(aoi) is 4:
        aoi = [aoi[i] for i in [0,3,2,1]]

    tiles = MODIStiles(aoi)
    tileRX = re.compile('|'.join(tiles.tiles))

    # files for tile result
    h5files = [y for x in os.walk(args.prcdir) for y in glob.glob(os.path.join(x[0], '*.h5')) if re.search(tileRX,y)]


    # filter for product (and parameter)
    if args.parameter:
        h5files_fil = [x for x in h5files if args.product in x and agrs.parameter in x]
    else:
        h5files_fil = [x for x in h5files if args.product in x]

    # get window
    win = MODISwindow(aoi=aoi,datemin=args.datemin,datemax=args.datemax,files=h5files)

    try:
        arr = np.zeros((win.height,win.width,sum(win.ix)),dtype='float32')

    except MemoryError:
        raise SystemExit('Error allocating memory. Please specify smaller AOI and/or shorter date range.')

    for f in h5files:

        h5f = h5_file(f)

        xr,xrd,yr,yrd = aoi2ix(h5f.bbox(),aoi,h5f.res)

        xw,xwd,yw,ywd = aoi2ix(h5f.bbox(),aoi,h5f.res)

        arr[yw:(yw+ywd),xw:(xw+xwd),...] = h5_readArr(f,xr,xrd,yr,yrd,np.flatnonzero(win.ix),args.dataset)


if __name__=='__main__':
    main()
