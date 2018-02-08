import numpy as np
import gdal
import h5py
import glob
import time
from .utils import block_view
#create hdf5 file

def createH5(files, convfun, name = None, minrows = 10, compression = 'gzip'):
    print('Initializing createH5 - start processing at %s ...' % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    fopen = gdal.Open(files[0])
    arr = fopen.ReadAsArray()
    [rows,cols] = arr.shape
    trans = fopen.GetGeoTransform()
    proj = fopen.GetProjection()
    nblocks = len(range(0,rows,minrows))
    print('Template loaded! Creating HDF5 file ...')
    if name is None:
        name = ('dataset-%s.h5' % time.strftime("%Y_%m_%d-%H%M%S", time.localtime()))
    with h5py.File(name,'w',libver='latest') as h5f:
        chunks = (minrows,cols,len(files))
        dset = h5f.create_dataset('Raw',shape=(rows,cols,len(files)),dtype='f',maxshape=(rows,cols,None),chunks=chunks,compression=compression)#32001
        h5f.create_dataset('Smooth',shape=(rows,cols,len(files)),dtype='f',maxshape=(rows,cols,None),chunks=chunks,compression=compression)#32001
        h5f.create_dataset('lgrd',shape=(rows,cols),dtype='f',maxshape=(rows,cols),chunks=chunks[0:2],compression=compression)
        dset.attrs['Extent'] = trans
        dset.attrs['Projection'] = proj
        print('Writing data in blocks ...')
        for blk in range(nblocks):
            print("Block %i of %i" %(blk+1,nblocks))
            block = np.empty([minrows,cols,len(files)])
            for ix,f in enumerate(files):
                currfl = gdal.Open(f)
                imarr = block_view(currfl.ReadAsArray(),block=(minrows,cols))
                im = convfun(im)
                block[:,:,ix] = imarr[blk,0,...]
            dset[blk*minrows:blk*minrows+minrows,...] = block

#update hdf5 file

def updateH5(files,convfun ,h5file):
    with h5py.File(h5file,'r+') as h5f:
        dset = h5f.get('Raw')
        refsize = dset.shape[2]
        dset.resize((dset.shape[0],dset.shape[1],dset.shape[2]+len(files)))
        for ix,f in enumerate(files):
            currfl = gdal.Open(f)
            im = currfl.ReadAsArray()
            im = convfun(im)
            dset[:,:,refsize+ix:refsize+(ix+1)] = im
