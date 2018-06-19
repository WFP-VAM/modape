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
from progress.spinner import Spinner
from .utils import LDOM, dtype_GDNP, SessionWithHeaderRedirection
from contextlib import contextmanager
import pickle
import warnings
import itertools
import bisect
import gc
try:
    import gdal
except ImportError:
    from osgeo import gdal

# turn off BeautifulSoup warnings
warnings.filterwarnings("ignore", category=UserWarning, module='bs4')

class MODISquery:

    def __init__(self,url,begindate,enddate,username=None,password=None,rawdir=os.getcwd(),global_flag=None,wget=False):

        self.queryURL = url
        self.username = username
        self.password = password
        self.rawdir = rawdir
        self.files = []
        self.modisURLs = []
        self.begin = datetime.datetime.strptime(begindate,'%Y-%m-%d').date()
        self.end = datetime.datetime.strptime(enddate,'%Y-%m-%d').date()
        self.global_flag = global_flag
        self.wget = wget

        with requests.Session() as sess:

            print('Checking for MODIS products ...',end='')
            try:
                response = sess.get(self.queryURL)
                self.statuscode = response.status_code
                response.raise_for_status()

            except requests.exceptions.RequestException as e:
                print(e)
                sys.exit(1)

            soup = BeautifulSoup(response.content)

            if self.global_flag:

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
            print('{} results found.\n'.format(self.results))
        else:
            print('0 results found. Please check query!')


    def setCredentials(self,username,password):
        self.username=username
        self.password=password

    def download(self):

        if self.username is None or self.password is None:
            raise SystemExit('No credentials found. Please run .setCredentials(username,password)!')

        print('[%s]: Downloading products to %s ...\n' % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),self.rawdir))

        if self.wget:

            try:
                temp = check_output(['wget', '--version'])
            except:
                raise SystemExit("WGET download needs WGET to be available in PATH! Please make sure it's installed and available in PATH!")

            with open(self.rawdir + '/MODIS_filelist.txt','w') as flist:
                for item in self.modisURLs:
                    flist.write("%s\n" % item)

            args = ['wget','-q','-nd','-nc','-np','-r','-l1','-A','hdf','--show-progress','--progress=bar:force','--no-check-certificate','--user',self.username,'--password',self.password,'-P',self.rawdir]

            p = Popen(args + ['-i','{}/MODIS_filelist.txt'.format(self.rawdir)])
            p.wait()
            if p.returncode is not 0:
                print("Error occured during download, please check files against MODIS_filelist.txt!")
            else:
                os.remove(self.rawdir + '/MODIS_filelist.txt')


            self.files = [self.rawdir + os.path.basename(x) for x in self.modisURLs]


        else:

            r = re.compile('.*.hdf$')

            session = SessionWithHeaderRedirection(self.username, self.password)

            for ix,u in enumerate(self.modisURLs):
                print('%s of %s' %(ix+1,self.results))

                if self.global_flag:

                    try:
                        resp_temp = session.get(u)

                    except requests.exceptions.RequestException as e:
                        print(e)
                        print('Error accessing {} - skipping.'.format(u))
                        continue

                    soup_temp = BeautifulSoup(resp_temp.content)

                    hrefs = soup_temp.find_all('a',href=True)

                    hdf_file = [x.getText() for x in hrefs if re.match(r,x.getText())]

                    try:
                        u = u + hdf_file[0]

                    except IndexError:
                        print('No HDF file found in {} - skipping.'.format(d_url))
                        continue

                fname = u[u.rfind('/')+1:]

                if os.path.exists('{}/{}'.format(self.rawdir,fname)):
                    print('\nSkipping {} - {} already exists in {}!\n'.format(u,fname,self.rawdir))
                    continue


                try:
                    response = session.get(u, stream=True)
                    response.raise_for_status()

                    spinner = Spinner('Downloading {} ... '.format(fname))

                    with open(fname, 'wb') as fopen:
                        for chunk in response.iter_content(chunk_size=1024*1024):
                            fopen.write(chunk)
                            spinner.next()

                    self.files = self.files + [self.rawdir + fname]
                    print(' done.\n')
                except requests.exceptions.HTTPError as e:
                    print('Error downloading {} - skipping. Error message: {}'.format(u,e))
                    continue

        print('\n[{}]: Downloading finished.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))



class MODIShdf5:

    def __init__(self,files,param=None,targetdir=os.getcwd(),compression='gzip',chunks=(120,120,16)):

        self.targetdir = targetdir
        #self.resdict = dict(zip(['250m','500m','1km','0.05_Deg'],[x/112000 for x in [250,500,1000,5600]])) ## commented for original resolution
        self.paramdict = dict(zip(['VIM','VEM','LTD','LTN'],['NDVI','EVI','LST_Day','LST_Night']))
        self.compression = compression
        self.chunks = chunks
        self.dts_regexp = re.compile(r'.+A(\d{7}).+')
        self.dates = [re.findall(self.dts_regexp,x)[0] for x in files]

        self.files = [x for (y,x) in sorted(zip(self.dates,files))]
        self.dates.sort()
        self.nfiles = len(self.files)
        self.ref_file = self.files[0]
        self.ref_file_basename = os.path.basename(self.ref_file)


        if not re.match(r'M.D13\w\d',self.ref_file_basename) and not re.match(r'M.D11\w\d',self.ref_file_basename):
            raise SystemExit("Processing only implemented for M*D11 or M*13 products!")

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

        if self.param is 'VIM' and any(['MOD' in os.path.basename(x) for x in files]) and any(['MYD' in os.path.basename(x) for x in files]):
            self.product = ['MXD']
        else:
            self.product = re.findall(ppatt,self.ref_file_basename)

        self.outname = '{}/{}/{}_{}.h5'.format(
                                    self.targetdir,
                                    self.param,
                                    '.'.join(self.product + re.findall(tpatt,self.ref_file_basename) + [re.sub(vpatt,'\\1',self.ref_file_basename)]),
                                    self.param)

        self.exists = os.path.isfile(self.outname)
        ref = None


    def create(self):
        print('\nCreating file: %s ... ' % self.outname, end='')

        ref = gdal.Open(self.ref_file)
        ref_sds = [x[0] for x in ref.GetSubDatasets() if self.paramdict[self.param] in x[0]][0]
        ref_doy = [x[0] for x in ref.GetSubDatasets() if 'day of the year' in x[0]][0]

        #res = [value for key, value in self.resdict.items() if key in ref_sds][0] ## commented for original resolution

        ## commented for original resolution
        #rst = gdal.Warp('', ref_sds, dstSRS='EPSG:4326', format='VRT',
        #                              outputType=gdal.GDT_Float32, xRes=res, yRes=res)

        ## TODO real raster resolution??

        rst = gdal.Open(ref_sds)

        ref_sds = None

        self.rows = rst.RasterYSize
        self.cols = rst.RasterXSize
        self.nodata_value = int(rst.GetMetadataItem('_FillValue'))

        if re.match(r'M.{1}D13\w\d',self.ref_file_basename):
            self.numberofdays = 16
        elif re.match(r'M.{1}D11\w\d',self.ref_file_basename):
            self.numberofdays = 8
        dt = rst.GetRasterBand(1).DataType

        try:
            self.datatype = dtype_GDNP(dt)
        except IndexError:
            print("\n\n Couldn't read data type from dataset. Using default Int16!\n")
            self.datatype = (3,'int16')

        trans = rst.GetGeoTransform()
        prj = rst.GetProjection()

        rst = None

        if not os.path.exists(os.path.dirname(self.outname)):
            os.makedirs(os.path.dirname(self.outname))

        if ref_doy:
            doyflag = True
        else:
            doyflag = False

        try:

            with h5py.File(self.outname,'x',libver='latest') as h5f:
                dset = h5f.create_dataset('Raw',shape=(self.rows,self.cols,self.nfiles * self.numberofdays),dtype=self.datatype[1],maxshape=(self.rows,self.cols,None),chunks=self.chunks,compression=self.compression)
                #h5f.create_dataset('Smooth',shape=(self.rows,self.cols,self.nfiles),dtype=self.datatype[1],maxshape=(self.rows,self.cols,None),chunks=self.chunks,compression=self.compression)
                #h5f.create_dataset('lgrd',shape=(self.rows,self.cols),dtype='float32',maxshape=(self.rows,self.cols),chunks=self.chunks[0:2],compression=self.compression)
                h5f.create_dataset('Dates',shape=(self.nfiles,),maxshape=(None,),dtype='S8',compression=self.compression)
                dset.attrs['Geotransform'] = trans
                dset.attrs['Projection'] = prj
                dset.attrs['Resolution'] = trans[1] # res ## commented for original resolution
                dset.attrs['doyflag'] = doyflag
                dset.attrs['nodata'] = self.nodata_value
                dset.attrs['numberofdays'] = self.numberofdays

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
                self.nodata_value = dset.attrs['nodata']
                self.numberofdays = dset.attrs['numberofdays']
                self.doyindex = int(self.numberofdays / 2)
                self.datatype = dtype_GDNP(dset.dtype.name)
                #res  = dset.attrs['Resolution'] ## comment for original resolution

                dates = [x.decode() for x in dts[...] if len(x) > 0]
                [bisect.insort_left(dates,x) for x in self.dates if x not in dates]

                if len(dates) > dts.shape[0]:
                    dts.resize((len(dates),))
                    dset.resize((dset.shape[0],dset.shape[1],dset.shape[2] + len(dates) * self.numberofdays))

                dts[...] = [n.encode("ascii", "ignore") for n in dates]

                [gc.collect() for x in range(3)]

                try:
                    arr = np.zeros(self.chunks,dtype=self.datatype[1])

                except MemoryError:
                        print("\n\n Can't allocate arrays for block due to memory restrictions! Make sure enough RAM is availabe, consider using a 64bit PYTHON version or reduce block size.\n\n Traceback:")
                        raise

                if dset.attrs['doyflag']:

                    valarr = np.zeros(self.chunks[0:2],dtype=self.datatype[1])
                    doyarr = np.zeros(self.chunks[0:2],dtype=self.datatype[1])


                    I,J = np.ogrid[:self.chunks[0],:self.chunks[1]]

                    bar = Bar('Processing',fill='=',max=self.nfiles,suffix='%(percent)d%%  ')
                    bar.goto(0)

                    for fl in self.files:

                        try:

                            flix = dates.index(re.sub(self.dts_regexp,'\\1',fl))

                            ix = [int((datetime.datetime.strptime(dates[flix],'%Y%j').date() + datetime.timedelta(x)).strftime('%j')) for x in range(16)]

                            fl_o = gdal.Open(fl)

                            val_sds = [x[0] for x in fl_o.GetSubDatasets() if self.paramdict[self.param] in x[0]][0]

                            doy_sds = [x[0] for x in fl_o.GetSubDatasets() if 'day of the year' in x[0]][0]

                            blks = itertools.product(range(0,self.rows,self.chunks[0]),range(0,self.cols,self.chunks[1]))

                            val_rst = gdal.Open(val_sds)

                            doy_rst = gdal.Open(doy_sds)

                            for blk in blks:

                                valarr[...] = val_rst.ReadAsArray(xoff=blk[1],yoff=blk[1],xsize=self.chunks[1],ysize=self.chunks[1])

                                doyarr[...] = doy_rst.ReadAsArray(xoff=blk[1],yoff=blk[1],xsize=self.chunks[1],ysize=self.chunks[1])

                                for doy_ix, doy in enumerate(ix):

                                    doyarr[doyarr == doy] = doy_ix
                                doyarr[doyarr < 0] = self.doyindex # set -1 to mid value
                                doyarr[doyarr > doy_ix] = doy_ix # clip to max doy

                                arr[...] =  dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix*self.numberofdays:flix*self.numberofdays+self.numberofdays]

                                arr[arr == 0] = self.nodata_value

                                arr[I,J,doyarr] = np.maximum.reduce([arr[I,J,doyarr],valarr[...]])

                                dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix*self.numberofdays:flix*self.numberofdays+self.numberofdays] = arr[...]


                        except AttributeError:

                            print('Error reading {} ... using NoData value {}.'.format(fl,self.nodata_value))

                            blks = itertools.product(range(0,self.rows,self.chunks[0]),range(0,self.cols,self.chunks[1]))

                            ndarr = np.full((self.chunks[0],self.chunks[1]),self.nodata_value,dtype=self.datatype[1])

                            for blk in blks:


                                arr[...] =  dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix*self.numberofdays:flix*self.numberofdays+self.numberofdays]

                                arr[arr == 0] = self.nodata_value

                                arr[...,self.doyindex] = np.maximum.reduce([arr[...,self.doyindex],ndarr])

                                dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix*self.numberofdays:flix*self.numberofdays+self.numberofdays] = arr[...]

                            del ndarr


                        fl_o = None
                        val_sds = None
                        doy_sds = None
                        val_rst = None
                        doy_rst = None

                        del fl_o, val_sds, val_rst, doy_sds, doy_rst

                        bar.next()
                    bar.finish()

                else:

                    valarr = np.zeros(self.chunks[0:2],dtype=self.datatype[1])

                    bar = Bar('Processing',fill='=',max=self.nfiles,suffix='%(percent)d%%  ')
                    bar.goto(0)

                    for fl in self.files:

                        try:

                            flix = dates.index(re.sub(self.dts_regexp,'\\1',fl))

                            ix = [int((datetime.datetime.strptime(dates[flix],'%Y%j').date() + datetime.timedelta(x)).strftime('%j')) for x in range(16)]

                            fl_o = gdal.Open(fl)

                            val_sds = [x[0] for x in fl_o.GetSubDatasets() if self.paramdict[self.param] in x[0]][0]

                            blks = itertools.product(range(0,self.rows,self.chunks[0]),range(0,self.cols,self.chunks[1]))

                            val_rst = gdal.Open(val_sds)

                            for blk in blks:

                                valarr[...] = val_rst.ReadAsArray(xoff=blk[1],yoff=blk[1],xsize=self.chunks[1],ysize=self.chunks[1])

                                arr[...] =  dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix*self.numberofdays:flix*self.numberofdays+self.numberofdays]

                                arr[arr == 0] = self.nodata_value

                                arr[...,self.doyindex] = np.maximum.reduce([arr[...,self.doyindex],valarr[...]])

                                dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix*self.numberofdays:flix*self.numberofdays+self.numberofdays] = arr[...]

                        except AttributeError:

                            print('Error reading {} ... using NoData value {}.'.format(fl,self.nodata_value))

                            blks = itertools.product(range(0,self.rows,self.chunks[0]),range(0,self.cols,self.chunks[1]))

                            ndarr = np.full((self.chunks[0],self.chunks[1]),self.nodata_value,dtype=self.datatype[1])

                            for blk in blks:

                                arr[...] =  dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix*self.numberofdays:flix*self.numberofdays+self.numberofdays]

                                arr[arr == 0] = self.nodata_value

                                arr[...,self.doyindex] = np.maximum.reduce([arr[...,self.doyindex],ndarr])

                                dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix*self.numberofdays:flix*self.numberofdays+self.numberofdays] = arr[...]

                            del ndarr

                        fl_o = None
                        val_sds = None
                        val_rst = None

                        del fl_o, val_sds, val_rst

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
                self.datatype = dset.dtype
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

        array = np.zeros(((len(self.v_ix) * self.tile_rws),len(self.h_ix) * self.tile_cls),dtype=self.datatype)

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

        try:
            dt_gdal = dtype_GDNP(self.datatype)
        except IndexError:
            print("\n\n Couldn't read data type from dataset. Using default Int16!\n")
            dt_gdal = (3,'int16')

        driver = gdal.GetDriverByName('GTiff')

        self.raster = driver.Create('/vsimem/inmem.tif', width, height, 1, dt_gdal[0])

        self.raster.SetGeoTransform(self.gt)
        self.raster.SetProjection(self.pj)

        self.raster.GetRasterBand(1).WriteArray(array)

        yield self

        gdal.Unlink('/vsimem/inmem.tif')
        self.raster = None
        driver = None
        del array
