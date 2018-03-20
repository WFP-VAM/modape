from __future__ import print_function, division
import numpy as np
import requests
from bs4 import BeautifulSoup
import re
import sys, os
import time
from subprocess import Popen, check_output
import tables
import h5py
from progress.bar import Bar
import osr
try:
    import gdal
except ImportError:
    from osgeo import gdal


class MODISquery:

    def __init__(self,url,username=None,password=None,rawdir=os.getcwd()):

        self.queryURL = url
        self.username = username
        self.password = password
        self.rawdir = rawdir
        self.files = []
        self.minrows = 112

        r = re.compile(".+(h\d+v\d+).+")

        print('Checking for MODIS products ...',end='')
        try:
            response = requests.get(url)
            self.statuscode = response.status_code
            response.raise_for_status()

        except requests.exceptions.RequestException as e:
            print(e)
            sys.exit(1)

        soup = BeautifulSoup(response.content,"html5lib")

        self.modisURLs = [x.getText() for x in soup.find_all('url')]
        self.results = len(self.modisURLs)
        self.tiles = list(set([r.search(x).group(1) for x in self.modisURLs]))

        print('... done.\n')

        if self.results > 0:
            print('%s results found.' % self.results)
        else:
            print('0 results found. Please check query!')


    def setCredentials(self,username,password):
        self.username=username
        self.password=password

    def download(self):

        try:
            temp = check_output(['wget', '--version'])
        except:
            print("Download needs wget to be available in PATH! Please make sure it's installed and available in PATH!")
            sys.exit(1)

        if self.username is None or self.password is None:
            print('No credentials found. Please run .setCredentials(username,password)!')
            sys.exit(1)

        args = ['wget','-q','--show-progress','--progress=bar:force','--no-check-certificate','--user',self.username,'--password',self.password,'-P',self.rawdir]

        print('[%s]: Downloading products to %s ...\n' % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),self.rawdir))
        for ix,u in enumerate(self.modisURLs):
            print('%s of %s' %(ix+1,self.results))
            p = Popen(args + [u])
            p.wait()
            if p.returncode is not 0:
                print("Couldn't download %s - continuing." % u)
                continue
            self.files = self.files + [self.rawdir + os.path.basename(u)]

        print('\n[%s]: Downloading finished.' % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))



class MODIShdf5:

    def __init__(self,files,param=None,targetdir=os.getcwd(),compression=32001):

        self.targetdir = targetdir
        self.resdict = dict(zip(['250m','500m','1km','0.05_Deg'],[x/112000 for x in [250,500,1000,5600]]))
        self.paramdict = dict(zip(['VIM','VEM','LTD','LTN'],['NDVI','EVI','LST_Day','LST_Night']))
        self.minrows = 112
        self.compression = compression
        self.dts_regexp = re.compile(r'.+A(\d{7}).+')
        self.dates = [re.findall(self.dts_regexp,x)[0] for x in files]
        self.dates.sort()

        self.files = [x for (y,x) in sorted(zip(self.dates,files))]
        self.nfiles = len(self.files)
        self.ref_file = self.files[0]

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
                                    '.'.join([os.path.basename(self.ref_file).split('.')[i] for i in [0,2,3]]),
                                    self.param)

        self.exists = os.path.isfile(self.outname)

    def create(self):
        print('\nCreating file: %s ... ' % self.outname, end='')

        ref = gdal.Open(self.ref_file)
        ref_sds = [x[0] for x in ref.GetSubDatasets() if self.paramdict[self.param] in x[0]][0]
        res = [value for key, value in self.resdict.items() if key in ref_sds][0]

        rst = gdal.Warp('', ref_sds, dstSRS='EPSG:4326', format='VRT',
                                      outputType=gdal.GDT_Float32, xRes=res, yRes=res)

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
                dset = h5f.create_dataset('Raw',shape=(self.rows,self.cols,self.nfiles),dtype='float32',maxshape=(self.rows,self.cols,None),chunks=self.chunks,compression=self.compression)
                h5f.create_dataset('Smooth',shape=(self.rows,self.cols,self.nfiles),dtype='float32',maxshape=(self.rows,self.cols,None),chunks=self.chunks,compression=self.compression)
                h5f.create_dataset('lgrd',shape=(self.rows,self.cols),dtype='float32',maxshape=(self.rows,self.cols),chunks=self.chunks[0:2],compression=self.compression)
                h5f.create_dataset('Dates',shape=(self.nfiles,),maxshape=(None,),dtype='S8',compression=self.compression)
                dset.attrs['Extent'] = trans
                dset.attrs['Projection'] = prj
                dset.attrs['Resolution'] = res
                dset.attrs['flag'] = False

            self.exists = True
            print('done.\n')

        except OSError as err:
            #raise RuntimeError ("{0}".format(err))
            raise

    def update(self):
        print('Processing MODIS files ...\n')
        bar = Bar('Processing',fill='=',max=self.nfiles,suffix='%(percent)d%%')
        with h5py.File(self.outname,'r+',libver='latest') as h5f:
            dset = h5f.get('Raw')
            dts  = h5f.get('Dates')
            res  = dset.attrs['Resolution']

            if dset.attrs['flag']:
                uix = dset.shape[2]
                dset.resize((dset.shape[0],dset.shape[1],dset.shape[2]+self.nfiles))
            else:
                uix = 0
                dset.attrs['flag'] = True

            dts[uix:uix+self.nfiles] = [n.encode("ascii", "ignore") for n in self.dates]

            for fix,fl in enumerate(self.files):
                fl_o = gdal.Open(fl)

                ref_sds = [x[0] for x in fl_o.GetSubDatasets() if self.paramdict[self.param] in x[0]][0]

                rst = gdal.Warp('', ref_sds, dstSRS='EPSG:4326', format='VRT', outputType=gdal.GDT_Float32, xRes=res, yRes=res)

                arr = rst.ReadAsArray()

                fl_o = None
                ref_sds = None
                rst = None

                if self.param in ['LTD','LTN']:
                    arr = (arr * 0.02) + (-273)
                elif self.param in ['VIM']:
                    arr[arr < 0] = 0
                    arr = arr * 0.0001
                else:
                    pass

                dset[...,uix+fix] = arr[...]
                bar.next()
        bar.finish()

        print('\ndone.\n')


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

            yo = round((self.aoi[0]-gt[3])/gt[5])
            yd = 1
            xo = round((self.aoi[0]-gt[3])/gt[5])
            xd = 1


        elif len(self.aoi) is 4:

            y = [round((x-gt[3])/gt[5]) for x in [aoi[1],aoi[3]]]
            yo = y[0]
            yd = abs(y[0] - y[1])
            x = [round((x-gt[0])/gt[1]) for x in [aoi[0],aoi[2]]]
            xo = x[0]
            xd = abs(x[0] - x[1])

        tile_extract = ds.ReadAsArray(xo,yo,xd,yd)
        ds = None
        tile_tmp = np.unique(tile_extract/100)
        tiles = ["{:05.2f}".format(x) for x in tile_tmp[tile_tmp != 0]]

        self.tiles = ["h{}v{}".format(*x.split('.')) for x in tiles]


class MODISwindow:

    def __init__(self,aoi,files):

        ## TODO docstring aoi order

        with h5py.File(files[0],'r') as h5f:
            dat = h5f.get('Raw')
            self.resolution = dat.attrs['Resolution']
            dat = None

        self.width = int(round(abs(aoi[2] - aoi[0]) / self.resolution))
        self.height = int(round(abs(aoi[3] - aoi[1]) / self.resolution))
        self.geotransform = [aoi[0],self.resolution,0.0,aoi[1],0.0,-self.resolution]
        self.projection = osr.SRS_WKT_WGS84
