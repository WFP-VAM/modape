#!/usr/bin/env python

import pickle
import os
import pandas as pd

def main():

    this_dir, this_filename = os.path.split(__file__)

    with open(os.path.join(this_dir, "data", "MODIS_V6_PT.pkl"),'rb') as table_raw:
        product_table = pickle.load(table_raw)

    tbl = pd.DataFrame(product_table).T

    with pd.option_context('display.max_rows', None, 'display.max_columns', 4):
        print(tbl)



if __name__=='__main__':
    main()
