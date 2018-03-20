from __future__ import print_function, division
import numpy as np
import re
import sys, os
import time
import datetime
import tables
import h5py

def h5_extractAOI(file,aoi,ix,dataset):

'''Extract AOI intersection from h5 file'''

    with h5py.File(file,'r') as h5file:
        ds = h5file.get('Raw')
        gt = ds.attrs['Geotransform']
        res = ds.attrs['Resolution']
        r,c,t = ds.shape


    x2 = gt[0] + (c*res)
    y2 = gt[3] - (r*res)

    isect = [max(gt[0],aoi[0]),min(gt[3],aoi[1]),min(x2,aoi[2]),max(y2,aoi[3])]

    xoff = int(round((isect[0] - gt[0])/res))
    yoff = int(round((gt[3] - isect[1])/res))

    xd = int(round((isect[2] - isect[0])/res))+1
    yd = int(round((isect[1] - isect[3])/res))+1


    with h5py.File(file,'r') as h5file:
        ds = h5file.get(dataset)
        arr = ds[yoff:(yoff+yd),xoff:(xoff+xd),ix]

    return(arr)
