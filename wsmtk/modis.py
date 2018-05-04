from __future__ import print_function, division
import numpy as np
import requests
from bs4 import BeautifulSoup
import re
import sys, os
import time
import datetime
from subprocess import Popen, check_output
import tables
import h5py
import osr
from progress.bar import Bar
from .utils import LDOM
from contextlib import contextmanager
import pickle
import warnings
import itertools
try:
    import gdal
except ImportError:
    from osgeo import gdal

# turn off BeautifulSoup warnings
warnings.filterwarnings("ignore", category=UserWarning, module='bs4')

class MODISquery:

    def __init__(self,url,begindate,enddate,username=None,password=None,rawdir=os.getcwd(),global_flag=None):

        self.queryURL = url
        self.username = username
        self.password = password
        self.rawdir = rawdir
        self.files = []
        self.modisURLs = []
        self.begin = datetime.datetime.strptime(begindate,'%Y-%m-%d').date()
        self.end = datetime.datetime.strptime(enddate,'%Y-%m-%d').date()

        print('Checking for MODIS products ...',end='')
        try:
            response = requests.get(self.queryURL)
            self.statuscode = response.status_code
            response.raise_for_status()

        except requests.exceptions.RequestException as e:
            print(e)
            sys.exit(1)

        soup = BeautifulSoup(response.content)#,"html5lib")

        if global_flag:

            r = re.compile('.*.hdf$')

            dates = np.array([x.getText() for x in soup.findAll('a',href=True) if re.match('\d{4}\.\d{2}\.\d{2}',x.getText())])
            dates_parsed = [datetime.datetime.strptime(x,'%Y.%m.%d/').date() for x in dates]

            dates_ix = np.flatnonzero(np.array([x >= self.begin and x < self.end for x in dates_parsed]))

            self.modisURLs = [self.queryURL + x for x in dates[dates_ix]]

        else:

            r = re.compile(".+(h\d+v\d+).+")

            self.modisURLs = [x.getText() for x in soup.find_all('url')]
            self.tiles = list(set([r.search(x).group(1) for x in self.modisURLs]))


        self.results = len(self.modisURLs)
        print('... done.\n')

        if self.results > 0:
            print('{} results found.'.format(self.results))
        else:
            print('0 results found. Please check query!')


    def setCredentials(self,username,password):
        self.username=username
        self.password=password

    def download(self):
        try:
            temp = check_output(['wget', '--version'])
        except:
            raise SystemExit("Download needs wget to be available in PATH! Please make sure it's installed and available in PATH!")


        if self.username is None or self.password is None:
            raise SystemExit('No credentials found. Please run .setCredentials(username,password)!')

        args = ['wget','-q','-nd','-nc','-np','-r','-l1','-A','hdf','--show-progress','--progress=bar:force','--no-check-certificate','--user',self.username,'--password',self.password,'-P',self.rawdir]

        print('[%s]: Downloading products to %s ...\n' % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),self.rawdir))
        for ix,u in enumerate(self.modisURLs):
            print('%s of %s' %(ix+1,self.results))
            p = Popen(args + [u])
            p.wait()
            if p.returncode is not 0:
                print("Couldn't download {} - continuing.".format(u))
                continue
            self.files = self.files + [self.rawdir + os.path.basename(u)]

        print('\n[{}]: Downloading finished.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))



class MODIShdf5:

    def __init__(self,files,param=None,targetdir=os.getcwd(),compression=32001):

        self.targetdir = targetdir
        #self.resdict = dict(zip(['250m','500m','1km','0.05_Deg'],[x/112000 for x in [250,500,1000,5600]])) ## commented for original resolution
        self.paramdict = dict(zip(['VIM','VEM','LTD','LTN'],['NDVI','EVI','LST_Day','LST_Night']))
        self.minrows = 120
        self.compression = compression
        self.dts_regexp = re.compile(r'.+A(\d{7}).+')
        self.dates = [re.findall(self.dts_regexp,x)[0] for x in files]

        self.files = [x for (y,x) in sorted(zip(self.dates,files))]
        self.dates.sort()
        self.nfiles = len(self.files)
        self.ref_file = self.files[0]
        self.ref_file_basename = os.path.basename(self.ref_file)


        ppatt = re.compile(r'M\w{6}')
        vpatt = re.compile('.+\.(\d{3})\..+')
        tpatt = re.compile(r'h\d+v\d+')

        ref = gdal.Open(self.ref_file)

        if not param:
            ref_sds = ref_sds = [x[0] for x in ref.GetSubDatasets() if self.paramdict['VIM'] in x[0] or self.paramdict['LTD'] in x[0]][0]
            self.param = [key for key, value in self.paramdict.items() if value in ref_sds][0]
            ref_sds = None
        elif param in self.paramdict.keys():
            self.param = param
        else:
            raise ValueError('Parameter string not recognized. Available parameters are %s.' % [x for x in self.paramdict.keys()])

        ref = None

        self.outname = '{}/{}/{}_{}.h5'.format(
                                    self.targetdir,
                                    self.param,
                                    '.'.join(re.findall(ppatt,self.ref_file_basename) + re.findall(tpatt,self.ref_file_basename) + [re.sub(vpatt,'\\1',self.ref_file_basename)]),
                                    self.param)

        self.exists = os.path.isfile(self.outname)

    def create(self):
        print('\nCreating file: %s ... ' % self.outname, end='')

        ref = gdal.Open(self.ref_file)
        ref_sds = [x[0] for x in ref.GetSubDatasets() if self.paramdict[self.param] in x[0]][0]
        #res = [value for key, value in self.resdict.items() if key in ref_sds][0] ## commented for original resolution

        ## commented for original resolution
        #rst = gdal.Warp('', ref_sds, dstSRS='EPSG:4326', format='VRT',
        #                              outputType=gdal.GDT_Float32, xRes=res, yRes=res)

        ## TODO real raster resolution??

        rst = gdal.Open(ref_sds)

        ref = None
        ref_sds = None

        self.rows = rst.RasterYSize
        self.cols = rst.RasterXSize
        self.chunks = (self.minrows,self.cols,self.nfiles)

        trans = rst.GetGeoTransform()
        prj = rst.GetProjection()

        rst = None

        if not os.path.exists(os.path.dirname(self.outname)):
            os.makedirs(os.path.dirname(self.outname))

        try:

            with h5py.File(self.outname,'x',libver='latest') as h5f:
                dset = h5f.create_dataset('Raw',shape=(self.rows,self.cols,self.nfiles),dtype='int16',maxshape=(self.rows,self.cols,None),chunks=self.chunks,compression=self.compression)
                h5f.create_dataset('Smooth',shape=(self.rows,self.cols,self.nfiles),dtype='int16',maxshape=(self.rows,self.cols,None),chunks=self.chunks,compression=self.compression)
                h5f.create_dataset('lgrd',shape=(self.rows,self.cols),dtype='float32',maxshape=(self.rows,self.cols),chunks=self.chunks[0:2],compression=self.compression)
                h5f.create_dataset('Dates',shape=(self.nfiles,),maxshape=(None,),dtype='S8',compression=self.compression)
                dset.attrs['Geotransform'] = trans
                dset.attrs['Projection'] = prj
                dset.attrs['Resolution'] = trans[1] # res ## commented for original resolution
                dset.attrs['flag'] = False

            self.exists = True
            print('done.\n')

        except:
            print('\n\nError creating {}! Check input parameters (especially if compression/filter is available) and try again. Corrupt file will be removed now. \n\nError message: \n'.format(self.outname))
            os.remove(self.outname)
            raise

    def update(self):
        print('Processing MODIS files ...\n')

        try:


            with h5py.File(self.outname,'r+',libver='latest') as h5f:
                dset = h5f.get('Raw')
                dts  = h5f.get('Dates')
                self.chunks = dset.chunks
                self.rows = dset.shape[0]
                self.cols = dset.shape[1]
                #res  = dset.attrs['Resolution'] ## comment for original resolution

                if dset.attrs['flag']:
                    uix = dset.shape[2]
                    dset.resize((dset.shape[0],dset.shape[1],dset.shape[2]+self.nfiles))
                else:
                    uix = 0
                    dset.attrs['flag'] = True

                dts[uix:uix+self.nfiles] = [n.encode("ascii", "ignore") for n in self.dates]

                blks = itertools.product(range(0,self.rows,self.chunks[0]),range(0,self.cols,self.chunks[1]),range(0,self.nfiles,self.chunks[2]))

                nblcks = len(list(range(0,self.rows,self.chunks[0]))) * len(list(range(0,self.nfiles,self.chunks[2])))

                bar = Bar('Processing',fill='=',max=nblcks,suffix='%(percent)d%%')
                bar.goto(0)

                for blk in blks:

                    arr = np.zeros((self.chunks[0],self.chunks[1],min(self.nfiles,self.chunks[2])),dtype='int16')

                    for fix,fl in enumerate(self.files[blk[2]:(blk[2]+arr.shape[2])]):

                        try:

                            fl_o = gdal.Open(fl)

                            ref_sds = [x[0] for x in fl_o.GetSubDatasets() if self.paramdict[self.param] in x[0]][0]

                            ## comment for original resolution
                            #rst = gdal.Warp('', ref_sds, dstSRS='EPSG:4326', format='VRT', outputType=gdal.GDT_Float32, xRes=res, yRes=res)

                            rst = gdal.Open(ref_sds)

                            arr[...,fix] = rst.ReadAsArray(yoff=blk[0],ysize=self.minrows)

                            fl_o = None
                            ref_sds = None
                            rst = None

                        except AttributeError:

                            print('Error reading {} ... using empty array.'.format(fl))

                            arr[...,fix] = np.zeros((self.chunks[0],self.chunks[1]),dtype='int16')

                    dset[blk[0]:(blk[0]+self.chunks[0]),:,uix:(uix+arr.shape[2])] = arr[...]
                    bar.next()
            bar.finish()

            print('\ndone.\n')

        except:
            print('Error updating {}! File may be corrupt, consider creating the file from scratch, or closer investigation. \n\nError message: \n'.format(self.outname))
            raise

    def __str__(self):
        return("MODIShdf5 object: %s - %s files - exists on disk: %s" % (self.outname, self.nfiles, self.exists))



class MODIStiles:

    def __init__(self,aoi):

        self.aoi = aoi

        this_dir, this_filename = os.path.split(__file__)
        ds = gdal.Open(os.path.join(this_dir, "data", "MODIS_TILES.tif"))

        try:
            gt = ds.GetGeoTransform()

        except AttributeError:
            raise SystemExit("Could not find 'MODIS_TILES.tif' index raster. Try reinstalling the package.")


        if len(self.aoi) is 2:

            xo = int(round((self.aoi[1]-gt[0])/gt[1]))
            yo = int(round((gt[3]-self.aoi[0])/gt[1]))

            xd = 1
            yd = 1

        elif len(self.aoi) is 4:

            xo = int(round((self.aoi[0]-gt[0])/gt[1]))
            yo = int(round((gt[3]-self.aoi[1])/gt[1]))

            xd = int(round((self.aoi[2] - self.aoi[0])/gt[1]))
            yd = int(round((self.aoi[1] - self.aoi[3])/gt[1]))

        tile_extract = ds.ReadAsArray(xo,yo,xd,yd)
        ds = None
        tile_tmp = np.unique(tile_extract/100)
        tiles = ["{:05.2f}".format(x) for x in tile_tmp[tile_tmp != 0]]

        self.tiles = ["h{}v{}".format(*x.split('.')) for x in tiles]


class MODISmosaic:

    def __init__(self,files,datemin,datemax):
        tile_re = re.compile('.+(h\d+v\d+).+')

        self.tiles = [re.sub(tile_re,'\\1',os.path.basename(x)) for x in files]
        self.tiles.sort()
        self.files = files
        self.h_ix = list(set([re.sub('(h\d+)(v\d+)','\\1',x) for x in self.tiles]))
        self.h_ix.sort()
        self.v_ix = list(set([re.sub('(h\d+)(v\d+)','\\2',x) for x in self.tiles]))
        self.v_ix.sort()

        # reference tile is top left
        ref = [x for x in self.files if self.tiles[0] in x][0]

        try:

            with h5py.File(ref,'r') as h5f:
                dset = h5f.get('Raw')
                r,c,t = dset.shape
                self.tile_rws = r
                self.tile_cls = c
                self.resolution = dset.attrs['Resolution']
                self.resolution_degrees = self.resolution/112000
                self.gt = dset.attrs['Geotransform']
                self.pj = dset.attrs['Projection']
                dset = None
                self.dates = [x.decode() for x in h5f.get('Dates')[...]]
        except Exception as e:
            print('\nError reading refreferece file {} for mosaic! Error message: {}\n'.format(ref,e))
            raise


        dts_dt = [datetime.datetime.strptime(x,'%Y%j').date() for x in self.dates]
        datemin_p = datetime.datetime.strptime(datemin,'%Y%m').date()
        datemax_p = datetime.datetime.strptime(datemax,'%Y%m').date()

        self.tempIX = np.flatnonzero(np.array([x >= datemin_p and x <= datemax_p for x in dts_dt]))

    def getArray(self,dataset,ix):

        array = np.zeros(((len(self.v_ix) * self.tile_rws),len(self.h_ix) * self.tile_cls),dtype='int16')

        for h5f in self.files:

            t_h = re.sub('.+(h\d+)(v\d+).+','\\1',os.path.basename(h5f))
            t_v = re.sub('.+(h\d+)(v\d+).+','\\2',os.path.basename(h5f))

            xoff = self.h_ix.index(t_h) * self.tile_cls
            yoff = self.v_ix.index(t_v) * self.tile_rws

            try:

                with h5py.File(h5f,'r') as h5f_o:
                    arr = h5f_o.get(dataset)[...,ix]
                array[yoff:(yoff+self.tile_rws),xoff:(xoff+self.tile_cls)] = arr[...]

            except Exception as e:
                print('Error reading data from file {} to array! Error message {}:\n'.format(h5f,e))
                raise

        return(array)

    @contextmanager
    def getRaster(self,dataset,ix):

        array = self.getArray(dataset,ix)

        height, width = array.shape

        driver = gdal.GetDriverByName('GTiff')

        self.raster = driver.Create('/vsimem/inmem.tif', width, height, 1, gdal.GDT_Int16)

        self.raster.SetGeoTransform(self.gt)
        self.raster.SetProjection(self.pj)

        self.raster.GetRasterBand(1).WriteArray(array)

        yield self

        gdal.Unlink('/vsimem/inmem.tif')
        self.raster = None
        driver = None
        del array
