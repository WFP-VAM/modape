from __future__ import print_function, division
import numpy as np
import re
import sys, os
import time
import datetime
import tables
import h5py
from .utils import aoi2ix

class h5_file:

    '''MODIS h5 file class'''

    def __init__(self,file):
        with h5py.File(file,'r') as h5file:
            ds = h5file.get('Raw')
            gt = ds.attrs['Geotransform']
            res = ds.attrs['Resolution']
            r,c,t = ds.shape

        #self.gt = gt
        self.res = res
        self.cols = c
        self.rows = r
        self.xmin = gt[0]
        self.xmax = gt[0] + (c*res)
        self.ymin = gt[3] - (r*res)
        self.ymax = gt[3]

    def bbox(self):
        return((self.xmin,self.ymax,self.xmax,self.ymin))


def h5_readArr(file,xoff,xd,yoff,yd,ix,dataset):

    with h5py.File(file,'r') as h5open:
        ds = h5open.get(dataset)
        arr = ds[yoff:(yoff+yd),xoff:(xoff+xd),ix]

    return(arr)
