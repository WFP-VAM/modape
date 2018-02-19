from __future__ import print_function, division
import requests
from bs4 import BeautifulSoup
import re
import sys, os
import glob
import time
from subprocess import Popen, check_output
from .utils import block_view
import h5py
try:
    import gdal
except ImportError:
    from osgeo import gdal


class MODISquery:

    def __init__(self,url,username=None,password=None,rawdir=os.getcwd(),targetdir=os.getcwd()):

        self.queryURL = url
        self.username = username
        self.password = password
        self.rawdir = rawdir
        self.targetdir = targetdir
        self.files = []
        self.res_m = ['250m','500m','1km','0.05_Deg']
        self.res_dg = [x/112000 for x in [250,500,1000,5600]]
        self.minrows = 112


        r = re.compile(".+(h\d+v\d+).+")

        print('Checking for MODIS products ...')
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


        args = ['wget','--user',self.username,'--password',self.password,'-P',self.rawdir]

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


    def process(self):
        pass

        ## CHANGE DESIGN:

        ## BLOCKWISE PROCESSING -> LOOP OVER TEMPORAL DIMENSION PER BLOCK
        ## PROCESS FUNCTION AS CLASS FUNCTION OR GLOBAL?
        ## GDALWARP IN MEMORY BENEFICIAL?

'''

        print('[%s]: Starting processing ...\n' % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        for f in self.files:

            fopen = gdal.Open(f)

            fopen_sds = [x[0] for x in fopen.GetSubDatasets() if "NDVI" in x[0] or 'LST_' in x[0]]

            for fs in fopen_sds:

                param = [['VIM','LTD','LTN'][ix] for ix,x in enumerate(['NDVI','LST_Day','LST_Night']) if x in fs]

                outname = '{}/{}/{}._{}.h5'.format(self.targetdir,
                                            param,
                                            [os.path.basename(f).split('.')[i] for i in [0,2,3]],
                                            param)

                res = [self.res_dg[ix] for ix,x in enumerate(self.res_m) if x in fs]

                ds = gdal.Warp('', fs, dstSRS='EPSG:4326', format='VRT',
                              outputType=gdal.GDT_Int16, xRes=res, yRes=res)


                rows = ds.RasterYSize
                cols = ds.RasterXSize

                nblocks = len(range(0,rows,self.minrows))
                block = np.empty([minrows,cols,len(files)])

                if not os.path.isfile(outname):
                    if not os.path.exists(os.dirname(outname)):
                        os.mkdir(os.dirname(outname))

                    trans = ds.GetGeoTransform()
                    proj = ds.GetProjection()

'''


def MODISprocess(modisHDFs,dstorage=None,rdstorage=None):
    pass
