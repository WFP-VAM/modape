#!/usr/bin/env python
import pickle
import os
import pandas as pd
import argparse

def main():
    '''Print MODIS product table.

    Prints either the entire MODIS product table, or a subset filtered by either
    product, resolution or parameter.
    '''

    parser = argparse.ArgumentParser(description="MODIS product table")
    parser.add_argument("--product", help='MODIS product ID',metavar='')
    parser.add_argument("--resolution", help='Filter for pixel size',type=str,metavar='')
    parser.add_argument("--parameter", help='Filter for VAM parameter',metavar='')

    args = parser.parse_args()

    # Load product table

    this_dir, this_filename = os.path.split(__file__)

    with open(os.path.join(this_dir, "data", "MODIS_V6_PT.pkl"),'rb') as table_raw:
        product_table = pickle.load(table_raw)

    tbl = pd.DataFrame(product_table).T

    print('\n')

    # Subset table by product
    if args.product:
        try:
            tbl = tbl.loc[args.product.upper()]
        except:
            raise SystemExit('Product ID not recognized!')

    # Subset table by resolution
    if args.resolution and not args.product:
        try:
            tbl = tbl[tbl['pixel_size'] == args.resolution]
        except:
            raise SystemExit('Resolution not valid! (possible: 250,500,1000,5600)')

    # Subset table by parameter
    if args.parameter and not args.product:
        d = dict(zip(['VIM','VEM','LTD','LTN'],['Vegetation Indices','Vegetation Indices','Temperature, Emissivity','Temperature, Emissivity']))
        try:
            tbl = tbl[tbl['product'] == d[args.parameter.upper()]]
        except:
            raise SystemExit('Parameter not valid! (possible: VIM, VEM, LTD, LTN)')


    # print results
    with pd.option_context('display.max_rows', None, 'display.max_columns', 4):
        print(tbl)
        print('\n')



if __name__=='__main__':
    main()
