from __future__ import print_function, division
import numpy as np
import requests
from bs4 import BeautifulSoup
import re
import sys, os
import time
import datetime
from subprocess import Popen, check_output
import h5py
from progress.bar import Bar
from progress.spinner import Spinner
from .utils import *
from .whittaker import ws2d, ws2d_vc, ws2d_vc_asy
from contextlib import contextmanager, closing
import warnings
import itertools
import bisect
import gc
import array
import multiprocessing as mp
import traceback
try:
    import gdal
except ImportError:
    from osgeo import gdal

# turn off BeautifulSoup warnings
warnings.filterwarnings("ignore", category=UserWarning, module='bs4')

class MODISquery:
    '''Class for querying and downloading MODIS data.'''

    def __init__(self,url,begindate,enddate,username=None,password=None,targetdir=os.getcwd(),global_flag=None,aria2=False):
        '''Creates a MODISquery object.

        Args:
            url (str): Query URL as created by downloadMODIS.py
            begindate (str): Begin of period for query as ISO 8601 date string (YYYY-MM-DD)
            enddate (str): End of period for query as ISO 8601 date string (YYYY-MM-DD)
            username (str): Earthdata username (only required for download)
            password (str): Earthdata password (only required for download)
            targetdir (str): Path to target directory for downloaded files (default cwd)
            global_flag (bool):Flag indictaing queried product si global file instead of tiled product
            aria2 (bool): Use aria2 for downloading rather then python's requests
        '''

        self.queryURL = url
        self.username = username
        self.password = password
        self.targetdir = targetdir
        self.files = []
        self.modisURLs = []
        self.begin = datetime.datetime.strptime(begindate,'%Y-%m-%d').date()
        self.end = datetime.datetime.strptime(enddate,'%Y-%m-%d').date()
        self.global_flag = global_flag
        self.aria2 = aria2

        # query for products using session object
        with requests.Session() as sess:

            print('Checking for MODIS products ...',end='')
            try:
                response = sess.get(self.queryURL)
                self.statuscode = response.status_code
                response.raise_for_status()

            except requests.exceptions.RequestException as e:
                print(e)
                sys.exit(1)

            soup = BeautifulSoup(response.content,features="html.parser")

            # results for global products are date directories on server, so need to query separately
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

        # check for results
        if self.results > 0:
            print('{} results found.\n'.format(self.results))
        else:
            print('0 results found. Please check query!')


    def setCredentials(self,username,password):
        '''Set Earthdata credentials.

        Sets Earthdata username and password in created MODISquery object.

        Args:
            username (str): Earthdata username
            password (str): Earthdata password
        '''

        self.username=username
        self.password=password

    def download(self):
        '''Downloads MODIS products.

        Download of files found through query, Earthdata username and password required!
        '''

        if self.username is None or self.password is None:
            raise SystemExit('No credentials found. Please run .setCredentials(username,password)!')

        print('[%s]: Downloading products to %s ...\n' % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),self.targetdir))

        # download using WGET if True
        if self.aria2:

            # fail if WGET not installed
            try:
                temp = check_output(['aria2c', '--version'])
            except:
                raise SystemExit("ARIA2 download needs ARIA2 to be available in PATH! Please make sure it's installed and available in PATH!")

            # write URLs of query resuls to disk for WGET
            with open(self.targetdir + '/MODIS_filelist.txt','w') as flist:
                for item in self.modisURLs:
                    flist.write("%s\n" % item)

            args = ['aria2c','--file-allocation=none','-c','-x','10','-s','10','--http-user',self.username,'--http-passwd',self.password,'-d',self.targetdir]

            # execute subprocess
            p = Popen(args + ['-i','{}/MODIS_filelist.txt'.format(self.targetdir)])
            p.wait()

            # remove filelist.txt if all downloads are successful
            if p.returncode is not 0:
                print("Error occured during download, please check files against MODIS_filelist.txt!")
            else:
                os.remove(self.targetdir + '/MODIS_filelist.txt')


            self.files = [self.targetdir + os.path.basename(x) for x in self.modisURLs]

        # download with requests
        else:

            r = re.compile('.*.hdf$')

            session = SessionWithHeaderRedirection(self.username, self.password)

            for ix,u in enumerate(self.modisURLs):
                print('%s of %s' %(ix+1,self.results))

                # for global files, URLs of HDF files have to be parsed for each date before download
                if self.global_flag:

                    try:
                        resp_temp = session.get(u)

                    except requests.exceptions.RequestException as e:
                        print(e)
                        print('Error accessing {} - skipping.'.format(u))
                        continue

                    soup_temp = BeautifulSoup(resp_temp.content,features="html.parser")

                    hrefs = soup_temp.find_all('a',href=True)

                    hdf_file = [x.getText() for x in hrefs if re.match(r,x.getText())]

                    try:
                        u = u + hdf_file[0]

                    except IndexError:
                        print('No HDF file found in {} - skipping.'.format(d_url))
                        continue

                fname = u[u.rfind('/')+1:]

                if os.path.exists('{}/{}'.format(self.targetdir,fname)):
                    print('\nSkipping {} - {} already exists in {}!\n'.format(u,fname,self.targetdir))
                    continue


                try:
                    response = session.get(u, stream=True)
                    response.raise_for_status()

                    spinner = Spinner('Downloading {} ... '.format(fname))

                    with open(fname, 'wb') as fopen:
                        for chunk in response.iter_content(chunk_size=1024*1024):
                            fopen.write(chunk)
                            spinner.next()

                    self.files = self.files + [self.targetdir + fname]
                    print(' done.\n')
                except requests.exceptions.HTTPError as e:
                    print('Error downloading {} - skipping. Error message: {}'.format(u,e))
                    continue

        print('\n[{}]: Downloading finished.'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))



class MODISrawh5:
    '''Class for raw MODIS data collected into HDF5 file, ready for smoothing.

    For MOD/MYD 13 products, MOD and MYD are interleaved into a combined MXD.
    '''

    def __init__(self,files,param=None,targetdir=os.getcwd(),compression='gzip',chunk=120*120):
        '''Create a MODISrawh5 class

        Args:
            files ([str]): A list of absolute paths to MODIS raw hdf files to be processed
            param (str): VAM parameter to be processed (default VIM/LTD)
            targetdir (str): Target directory for raw MODIS HDF5 file
            compression (str): Compression method to be used (default = gzip)
        '''

        self.targetdir = targetdir
        #self.resdict = dict(zip(['250m','500m','1km','0.05_Deg'],[x/112000 for x in [250,500,1000,5600]])) ## commented for original resolution
        self.paramdict = dict(zip(['VIM', 'VEM', 'LTD', 'LTN'], ['NDVI', 'EVI', 'LST_Day', 'LST_Night']))
        self.compression = compression
        self.dts_regexp = re.compile(r'.+A(\d{7}).+')
        self.dates = [re.findall(self.dts_regexp,x)[0] for x in files]
        self.chunk = chunk
        self.files = [x for (y,x) in sorted(zip(self.dates,files))]
        self.dates.sort()
        self.nfiles = len(self.files)
        self.ref_file = self.files[0]
        self.ref_file_basename = os.path.basename(self.ref_file)

        # class works only for M.D11 and M.D13 products
        if not re.match(r'M.D13\w\d',self.ref_file_basename) and not re.match(r'M.D11\w\d',self.ref_file_basename):
            raise SystemExit("Processing only implemented for M*D11 or M*13 products!")

        # Patterns for string extraction
        ppatt = re.compile(r'M\w{6}')
        vpatt = re.compile('.+\.(\d{3})\..+')
        tpatt = re.compile(r'h\d+v\d+')

        # Open reference file
        ref = gdal.Open(self.ref_file)

        # When no parameter is selected, the default is VIM and LTD
        if not param:
            ref_sds = ref_sds = [x[0] for x in ref.GetSubDatasets() if self.paramdict['VIM'] in x[0] or self.paramdict['LTD'] in x[0]][0]
            self.param = [key for key, value in self.paramdict.items() if value in ref_sds][0]
            ref_sds = None
        elif param in self.paramdict.keys():
            self.param = param
        else:
            raise ValueError('Parameter string not recognized. Available parameters are %s.' % [x for x in self.paramdict.keys()])

        ref = None

        # check for MOD/MYD interleaving
        if self.param is 'VIM' and any(['MOD' in os.path.basename(x) for x in files]) and any(['MYD' in os.path.basename(x) for x in files]):
            self.product = [re.sub(r'M[O|Y]D','MXD',re.findall(ppatt,self.ref_file_basename)[0])]
            self.temporalresolution = 8
        else:
            self.product = re.findall(ppatt,self.ref_file_basename)
            self.temporalresolution = None

        # Name of file to be created/updated
        self.outname = '{}/{}/{}.{}.h5'.format(
                                    self.targetdir,
                                    self.param,
                                    '.'.join(self.product + re.findall(tpatt,self.ref_file_basename) + [re.sub(vpatt,'\\1',self.ref_file_basename)]),
                                    self.param)

        self.exists = os.path.isfile(self.outname)
        ref = None


    def create(self):
        '''Creates the HDF5 file.'''

        print('\nCreating file: {} ... '.format(self.outname), end='')

        ref = gdal.Open(self.ref_file)
        ref_sds = [x[0] for x in ref.GetSubDatasets() if self.paramdict[self.param] in x[0]][0]

        ## to be removed
        #res = [value for key, value in self.resdict.items() if key in ref_sds][0] ## commented for original resolution

        ## commented for original resolution
        #rst = gdal.Warp('', ref_sds, dstSRS='EPSG:4326', format='VRT',
        #                              outputType=gdal.GDT_Float32, xRes=res, yRes=res)

        # reference raster
        rst = gdal.Open(ref_sds)

        ref_sds = None

        self.rows = rst.RasterYSize
        self.cols = rst.RasterXSize
        self.nodata_value = int(rst.GetMetadataItem('_FillValue'))

        if re.match(r'M[O|Y]D13\w\d',self.ref_file_basename):
            if not self.temporalresolution:
                self.numberofdays = 16
                self.temporalresolution = self.numberofdays

            else:
                self.numberofdays = 16

        elif re.match(r'M[O|Y]D11\w\d',self.ref_file_basename):
            self.numberofdays = 8
            self.temporalresolution = self.numberofdays

        # Read datatype
        dt = rst.GetRasterBand(1).DataType

        # Parse datatype - on error use default Int16
        try:
            self.datatype = dtype_GDNP(dt)
        except IndexError:
            print("\n\n Couldn't read data type from dataset. Using default Int16!\n")
            self.datatype = (3,'int16')

        self.chunks = (self.chunk,1)

        trans = rst.GetGeoTransform()
        prj = rst.GetProjection()

        rst = None

        # Create directory if necessary
        if not os.path.exists(os.path.dirname(self.outname)):
            os.makedirs(os.path.dirname(self.outname))

        # Create HDF5 file
        try:

            with h5py.File(self.outname,'x',libver='latest') as h5f:
                dset = h5f.create_dataset('data',shape=(self.rows*self.cols,self.nfiles),dtype=self.datatype[1],maxshape=(self.rows*self.cols,None),chunks=self.chunks,compression=self.compression,fillvalue=self.nodata_value)
                h5f.create_dataset('dates',shape=(self.nfiles,),maxshape=(None,),dtype='S8',compression=self.compression)
                dset.attrs['geotransform'] = trans
                dset.attrs['projection'] = prj
                dset.attrs['resolution'] = trans[1] # res ## commented for original resolution
                dset.attrs['nodata'] = self.nodata_value
                dset.attrs['numberofdays'] = self.numberofdays
                dset.attrs['temporalresolution'] = self.temporalresolution

            self.exists = True
            print('done.\n')

        except:
            print('\n\nError creating {}! Check input parameters (especially if compression/filter is available) and try again. Corrupt file will be removed now. \n\nError message: \n'.format(self.outname))
            os.remove(self.outname)
            raise

    def update(self):
        '''Ingest raw data into MODIS raw HDF5 file.

        When a new HDF5 file is created, uodate will also handle the first data ingest.
        '''

        print('Processing MODIS files ...\n')

        try:

            with h5py.File(self.outname,'r+',libver='latest') as h5f:
                dset = h5f.get('data')
                dts  = h5f.get('dates')
                self.chunks = dset.chunks
                self.rows = dset.shape[0]
                self.cols = dset.shape[1]
                self.nodata_value = dset.attrs['nodata'].item()
                self.numberofdays = dset.attrs['numberofdays'].item()
                self.temporalresolution = dset.attrs['temporalresolution'].item()
                self.datatype = dtype_GDNP(dset.dtype.name)
                #res  = dset.attrs['Resolution'] ## comment for original resolution
                dset.attrs['processingtimestamp'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

                # Insert new dates into existing date list
                dates = [x.decode() for x in dts[...] if len(x) > 0]
                [bisect.insort_left(dates,x) for x in self.dates if x not in dates]

                # if date list is bigger after insert, datasets need to be resized for additional data
                if len(dates) > dts.shape[0]:
                    diff = len(dts) - dts.shape[0]
                    dts.resize((len(dates),))
                    dset.resize((dset.shape[0],len(dates))

                # Write back date list
                dts[...] = [n.encode("ascii", "ignore") for n in dates]

                # Manual garbage collect to prevent out of memory
                [gc.collect() for x in range(3)]

                ### DIFF -> insert after read and then discard the DIFF!

                # Create array and catch MemoryError
                try:
                    arr = np.zeros((self.chunks[0],self.chunks[1]),dtype=self.datatype[1])

                except MemoryError:
                        print("\n\n Can't allocate arrays for block due to memory restrictions! Make sure enough RAM is availabe or consider using a 64bit PYTHON version.\n\n Traceback:")
                        raise


                valarr = np.zeros(self.chunks[0:2],dtype=self.datatype[1])

                bar = Bar('Processing',fill='=',max=self.nfiles,suffix='%(percent)d%%  ')
                bar.goto(0)

                for fl in self.files:

                    try:

                        flix = dates_daily.index(re.sub(self.dts_regexp,'\\1',fl))

                        fl_o = gdal.Open(fl)

                        val_sds = [x[0] for x in fl_o.GetSubDatasets() if self.paramdict[self.param] in x[0]][0]

                        blks = itertools.product(range(0,self.rows,self.chunks[0]),range(0,self.cols,self.chunks[1]))

                        val_rst = gdal.Open(val_sds)

                        for blk in blks:

                            valarr[...] = val_rst.ReadAsArray(xoff=blk[1],yoff=blk[0],xsize=self.chunks[1],ysize=self.chunks[0])

                            arr[...] =  dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix:flix+self.numberofdays]

                            # Doyindex is a previosly defiend point, by default the temporal midpoint
                            arr[...,self.doyindex] = np.maximum.reduce([arr[...,self.doyindex],valarr[...]])

                            dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix:flix+self.numberofdays] = arr[...]

                    except AttributeError:

                        print('Error reading {} ... using NoData value {}.'.format(fl,self.nodata_value))

                        blks = itertools.product(range(0,self.rows,self.chunks[0]),range(0,self.cols,self.chunks[1]))

                        ndarr = np.full((self.chunks[0],self.chunks[1]),self.nodata_value,dtype=self.datatype[1])

                        for blk in blks:

                            arr[...] =  dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix:flix+self.numberofdays]

                            arr[...,self.doyindex] = np.maximum.reduce([arr[...,self.doyindex],ndarr])

                            dset[blk[0]:(blk[0]+self.chunks[0]),blk[1]:(blk[1]+self.chunks[1]),flix:flix+self.numberofdays] = arr[...]

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
            traceback.print_exc()
            raise

    def __str__(self):
        '''String to be displayed wen printing an instance of the class object'''
        return("MODISrawh5 object: %s - %s files - exists on disk: %s" % (self.outname, self.nfiles, self.exists))


class MODISsmth5:
    '''Class for smoothed MODIS data collected into HDF5 file.'''

    def __init__(self,rawfile,tempint=None,nsmooth=None,nupdate=None,targetdir=os.getcwd(),parallel=False,nworkers=mp.cpu_count()-1):
        '''Create MODISsmth5 object.

        Args:
            rawfile (str): Full path to a MODISrawh5 file
            tempint (int): Integer specifying temporal interpolation (default is None, so native temporal resolution)
            nsmooth (int): Number of raw timesteps used for smoothing (default is all)
            nupdate (int): Number of smoothed timesteps to be updated (default is all)
            targetdir (str): Path to target directory for smoothed HDF5 file
            parallel (bool): Flag for use of parallel processing (default is False)
            nworkers (int): Number of worker processes used in parallel (requres parallel to be True)
        '''
        self.targetdir = targetdir
        self.rawfile = rawfile
        self.parallel = parallel
        self.nworkers = nworkers
        self.nupdate = nupdate

        # Get info from raw HDF5
        with h5py.File(self.rawfile,'r') as h5f:
            dset = h5f.get('data')
            dts = h5f.get('dates')
            self.numberofdays = int(dset.attrs['numberofdays'])
            rawshape = dset.shape

            self.rawdaily = [(fromjulian(dts[0].decode()) + datetime.timedelta(x)).strftime('%Y%j') for x in range(rawshape[2]+1)]

            # Number of timesteps for smoothing must be bigger than for updating
            if nsmooth and nupdate:
                assert nsmooth >= nupdate, "nsmooth >= nupdate!!!!"

            # If no nsmooth set, take all available timesteps
            if not nsmooth:
                firstday = dts[0].decode()
            else:
                firstday = dts[-nsmooth].decode()

            self.ndays = ((fromjulian(dts[-1].decode()) + datetime.timedelta(self.numberofdays)) - fromjulian(firstday)).days
            self.startix = self.rawdaily.index(firstday)
            self.daily = [(fromjulian(firstday) + datetime.timedelta(x)).strftime('%Y%j') for x in range(self.ndays+1)]

        # Parse tempint to get flag for filename
        try:
            txflag = txx(tempint)
        except ValueError:
            raise SystemExit('Value for temporal interpolation not valid (interger for number of days required)!')

        # Filename for smoothed HDF5
        self.outname = '{}/{}.tx{}.{}.h5'.format(
                                    self.targetdir,
                                    '.'.join(os.path.basename(rawfile).split('.')[:-2]),
                                    txflag,
                                    os.path.basename(rawfile).split('.')[-2:-1][0])

        self.exists = os.path.isfile(self.outname)

        # If txtflag == 'n', then native resolution will be used
        if txflag == 'n':
            self.temporalresolution = None
        else:
            self.temporalresolution = tempint


    def create(self):
        '''Creates smoothed HDF5 file on disk.'''

        # Try reading info from raw HDF5
        try:
            with h5py.File(self.rawfile,'r') as h5f:
                dset = h5f.get('data')
                dts = h5f.get('dates')
                dt = dset.dtype.name
                cmpr = dset.compression
                rows,cols,rawdays = dset.shape
                rawchunks = dset.chunks
                rgt = dset.attrs['geotransform']
                rpj = dset.attrs['projection']
                rres = dset.attrs['resolution']
                rnd = dset.attrs['nodata']
                rtres = dset.attrs['temporalresolution']
                firstday = fromjulian(dts[0].decode())

        except Exception as e:
            raise SystemExit('Error reading rawfile {}. File may be corrupt. \n\n Error message: \n\n {}'.format(self.rawfile,e))

        # Read native temporal resolution if no user input was supplied
        if not self.temporalresolution:
            self.temporalresolution = rtres

        # Create date list for smoothed timesteps
        dates = [self.rawdaily[ix] for ix in range(0,len(self.rawdaily),self.temporalresolution)]
        days = len(dates)

        self.chunks = (rawchunks[0],rawchunks[1],1)

        print('\nCreating file: {} ... '.format(self.outname), end='')

        try:
            with h5py.File(self.outname,'x',libver='latest') as h5f:
                dset = h5f.create_dataset('data',shape=(rows,cols,days),dtype=dt,maxshape=(rows,cols,None),chunks=self.chunks,compression=cmpr,fillvalue=rnd)
                h5f.create_dataset('sgrid',shape=(rows,cols),dtype='float32',maxshape=(rows,cols),chunks=self.chunks[0:2],compression=cmpr)
                h5f.create_dataset('dates',shape=(days,),maxshape=(None,),dtype='S8',compression=cmpr,data = [x.encode('ascii','ignore') for x in dates])
                dset.attrs['geotransform'] = rgt
                dset.attrs['projection'] = rpj
                dset.attrs['resolution'] = rres
                dset.attrs['nodata'] = rnd
                dset.attrs['temporalresolution'] = self.temporalresolution

        except:
            print('\n\nError creating {}! Check input parameters (especially if compression/filter is available) and try again. Corrupt file will be removed now. \n\nError message: \n'.format(self.outname))
            os.remove(self.outname)
            raise


        self.exists = True
        print('done.\n')

    def ws2d(self,s):
        '''Apply whittaker smoother with fixed s-value to data.

        Args:
            s (float): log10 value of s
        '''

        with h5py.File(self.rawfile,'r') as rawh5, h5py.File(self.outname,'r+') as smth5:

            raw_ds = rawh5.get('data')
            raw_dts = rawh5.get('dates')
            smt_ds = smth5.get('data')
            smt_dts = smth5.get('dates')

            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks

            nodata = raw_ds.attrs['nodata']
            t_interval = smt_ds.attrs['temporalresolution']

            # Store run parameters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            smt_ds.attrs['lastrun'] = "fixed s: log10(sopt) = {}".format(s)
            smt_ds.attrs['log10sopt'] = s


            # Check if file needs to be resized
            dates_check = [self.rawdaily[ix] for ix in range(0,len(self.rawdaily),t_interval)]

            # Resize if date list is bigger than shape of smoothed data
            if len(dates_check) > smt_dts.shape[0]:
                smt_dts.resize((len(dates_check),))
                smt_dts[...] = [x.encode("ascii", "ignore") for x in dates_check]
                smt_ds.resize((smt_ds.shape[0],smt_ds.shape[1],len(dates_check)))

            # Calculate update index and offsets
            if not self.nupdate:

                for i,d in enumerate(self.daily):
                    if d in dates_check:
                        self.smtoffset = dates_check.index(d)
                        self.rawoffset = i
                        break
            else:
                self.rawoffset = self.daily.index(dates_check[-self.nupdate])
                self.smtoffset = len(dates_check[:-self.nupdate])

            barmax = (rawshape[0]/rawchunks[0]) * (rawshape[1]/rawchunks[1])
            bar = Bar('Processing',fill='=',max=barmax,suffix='%(percent)d%%  ')
            bar.goto(0)

            # Process in parallel
            if self.parallel:

                # Initialize parameters for workers
                params = init_parameters(s=s,nd=nodata,dim=(rawchunks[0]*rawchunks[1],self.ndays))

                # Create shared memory array
                shared_array = init_shared(rawchunks[0] * rawchunks[1] * self.ndays)

                # Create numpy array from shared memory
                arr = tonumpyarray(shared_array)

                arr.shape = (rawchunks[0] * rawchunks[1],self.ndays)

                # 3d helper for read/write
                arr_helper = arr.view()

                arr_helper.shape = (rawchunks[0],rawchunks[1],self.ndays)

                blks = itertools.product(range(0,rawshape[0],rawchunks[0]),range(0,rawshape[1],rawchunks[1]))

                # Create worker pool
                with closing(mp.Pool(processes=self.nworkers,initializer = init_worker, initargs = (shared_array,params))) as pool:

                    # Iterate blocks
                    for b in blks:

                        # Read data in raw chunks
                        for ii in range(0,arr_helper.shape[2],rawchunks[2]):

                            arr_helper[...,ii:ii+rawchunks[2]] = raw_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.startix+ii:self.startix+ii+rawchunks[2]]

                        del ii

                        # Execute smoother in parallel, split indices by number of workers
                        res = pool.map(execute_ws2d,np.array_split(range(arr.shape[0]),self.nworkers))

                        # Write back data, float is rounded to INT
                        for ii,jj in enumerate(range(self.rawoffset,arr_helper.shape[2],t_interval)):

                            smt_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.smtoffset+ii] = arr_helper[...,jj].round()

                        del ii,jj

                        bar.next()
                    bar.finish()

            # Serial processing
            else:

                arr = np.zeros((rawchunks[0]*rawchunks[1],self.ndays),dtype='float32')

                # Create weights array
                wts = arr.copy()

                arr_helper = arr.view()
                arr_helper.shape = (rawchunks[0],rawchunks[1],self.ndays)

                blks = itertools.product(range(0,rawshape[0],rawchunks[0]),range(0,rawshape[1],rawchunks[1]))

                for b in blks:

                    for ii in range(0,arr_helper.shape[2],rawchunks[2]):

                        arr_helper[...,ii:ii+rawchunks[2]] = raw_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.startix+ii:self.startix+ii+rawchunks[2]]

                    del ii

                    # All ovbservations which are not nodata get weight 1, others 0
                    wts[...] = (arr != nodata) * 1

                    # Iterate pixels, skip pixels with only nodata
                    for r in range(arr.shape[0]):
                        if wts[r,...].sum().item() != 0.0:
                            arr[r,...] = ws2d(y = arr[r,...],lmda = s, w = wts[r,...])

                    # Write back data
                    for ii,jj in enumerate(range(self.rawoffset,arr_helper.shape[2],t_interval)):

                        smt_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.smtoffset+ii] = arr_helper[...,jj].round()

                    del ii,jj

                    bar.next()
                bar.finish()


    def ws2d_sgrid(self):
        '''Apply whittaker smootehr with fixed s to data.

        This fixed s version reads a pixel based s value from file, so it needs
        a previous run of V-curve s-optimization.
        '''

        with h5py.File(self.rawfile,'r') as rawh5, h5py.File(self.outname,'r+') as smth5:

            raw_ds = rawh5.get('data')
            raw_dts = rawh5.get('dates')
            smt_ds = smth5.get('data')
            smt_dts = smth5.get('dates')
            sgrid_ds = smth5.get('sgrid')

            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks

            nodata = raw_ds.attrs['nodata']
            t_interval = smt_ds.attrs['temporalresolution']

            # Store run parameters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            smt_ds.attrs['lastrun'] = "fixed s from grid"

            # Check if file needs to be resized
            dates_check = [self.rawdaily[ix] for ix in range(0,len(self.rawdaily),t_interval)]

            if len(dates_check) > smt_dts.shape[0]:
                smt_dts.resize((len(dates_check),))
                smt_dts[...] = [x.encode("ascii", "ignore") for x in dates_check]
                smt_ds.resize((smt_ds.shape[0],smt_ds.shape[1],len(dates_check)))

            # Calculate update index and offsets
            if not self.nupdate:

                for i,d in enumerate(self.daily):
                    if d in dates_check:
                        self.smtoffset = dates_check.index(d)
                        self.rawoffset = i
                        break
            else:
                self.rawoffset = self.daily.index(dates_check[-self.nupdate])
                self.smtoffset = len(dates_check[:-self.nupdate])


            barmax = (rawshape[0]/rawchunks[0]) * (rawshape[1]/rawchunks[1])
            bar = Bar('Processing',fill='=',max=barmax,suffix='%(percent)d%%  ')
            bar.goto(0)

            if self.parallel:

                params = init_parameters(nd=nodata,dim=(rawchunks[0]*rawchunks[1],self.ndays))

                shared_array = init_shared(rawchunks[0] * rawchunks[1] * self.ndays)

                # Add shared memory array for sgrid
                params['shared_sarr'] = init_shared(rawchunks[0] * rawchunks[1])

                arr = tonumpyarray(shared_array)

                sarr = tonumpyarray(params['shared_sarr'])

                arr.shape = (rawchunks[0] * rawchunks[1],self.ndays)
                sarr.shape = (rawchunks[0], rawchunks[1])

                arr_helper = arr.view()

                arr_helper.shape = (rawchunks[0],rawchunks[1],self.ndays)

                blks = itertools.product(range(0,rawshape[0],rawchunks[0]),range(0,rawshape[1],rawchunks[1]))

                with closing(mp.Pool(processes=self.nworkers,initializer = init_worker, initargs = (shared_array,params))) as pool:

                    for b in blks:

                        for ii in range(0,arr_helper.shape[2],rawchunks[2]):

                            arr_helper[...,ii:ii+rawchunks[2]] = raw_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.startix+ii:self.startix+ii+rawchunks[2]]

                        # Read s values
                        sarr[...] = sgrid_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1]]

                        del ii

                        res = pool.map(execute_ws2d_sgrid,np.array_split(range(arr.shape[0]),self.nworkers))

                        for ii,jj in enumerate(range(self.rawoffset,arr_helper.shape[2],t_interval)):

                            smt_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.smtoffset+ii] = arr_helper[...,jj].round()

                        del ii,jj

                        bar.next()
                    bar.finish()

            else:

                arr = np.zeros((rawchunks[0]*rawchunks[1],self.ndays),dtype='float32')
                wts = arr.copy()
                sarr = np.zeros((rawchunks[0]*rawchunks[1]),dtype='float32')

                arr_helper = arr.view()
                arr_helper.shape = (rawchunks[0],rawchunks[1],self.ndays)

                blks = itertools.product(range(0,rawshape[0],rawchunks[0]),range(0,rawshape[1],rawchunks[1]))

                for b in blks:

                    for ii in range(0,arr_helper.shape[2],rawchunks[2]):

                        arr_helper[...,ii:ii+rawchunks[2]] = raw_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.startix+ii:self.startix+ii+rawchunks[2]]

                    wts[...] = (arr != nodata) * 1

                    sarr[...] = sgrid_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1]].reshape(rawchunks[0] * rawchunks[1])

                    for r in range(arr.shape[0]):
                        if wts[r,...].sum().item() != 0.0:
                            arr[r,...] = ws2d(arr[r,...],lmda = 10**sarr[r],w = wts[r,...])

                    for ii,jj in enumerate(range(self.rawoffset,arr_helper.shape[2],t_interval)):

                        smt_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.smtoffset+ii] = arr_helper[...,jj].round()

                    del ii,jj

                    bar.next()
                bar.finish()

    def ws2d_vc(self,srange):
        '''Apply whittaker smoother with V-curve optimization of s to data.

        Args:
            srange (arr): Float32 array of s-values to apply
        '''

        with h5py.File(self.rawfile,'r') as rawh5, h5py.File(self.outname,'r+') as smth5:

            raw_ds = rawh5.get('data')
            raw_dts = rawh5.get('dates')
            smt_ds = smth5.get('data')
            smt_dts = smth5.get('dates')
            sgrid_ds = smth5.get('sgrid')

            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks

            nodata = raw_ds.attrs['nodata']
            t_interval = smt_ds.attrs['temporalresolution']

            # Store run parameters for infotool

            smt_ds.attrs['processingtimestamp'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            smt_ds.attrs['lastrun'] = "V-curve optimization of s"

            # Check if file needs to be resized

            dates_check = [self.rawdaily[ix] for ix in range(0,len(self.rawdaily),t_interval)]

            if len(dates_check) > smt_dts.shape[0]:
                smt_dts.resize((len(dates_check),))
                smt_dts[...] = [x.encode("ascii", "ignore") for x in dates_check]
                smt_ds.resize((smt_ds.shape[0],smt_ds.shape[1],len(dates_check)))

            # Calculate update index and offsets

            if not self.nupdate:

                for i,d in enumerate(self.daily):
                    if d in dates_check:
                        self.smtoffset = dates_check.index(d)
                        self.rawoffset = i
                        break
            else:
                self.rawoffset = self.daily.index(dates_check[-self.nupdate])
                self.smtoffset = len(dates_check[:-self.nupdate])

            barmax = (rawshape[0]/rawchunks[0]) * (rawshape[1]/rawchunks[1])
            bar = Bar('Processing',fill='=',max=barmax,suffix='%(percent)d%%  ')
            bar.goto(0)

            if self.parallel:

                params = init_parameters(nd=nodata,dim=(rawchunks[0]*rawchunks[1],self.ndays),srange=srange)

                shared_array = init_shared(rawchunks[0] * rawchunks[1] * self.ndays)

                params['shared_sarr'] = init_shared(rawchunks[0] * rawchunks[1])

                arr = tonumpyarray(shared_array)

                sarr = tonumpyarray(params['shared_sarr'])

                arr.shape = (rawchunks[0] * rawchunks[1],self.ndays)
                sarr.shape = (rawchunks[0], rawchunks[1])

                arr_helper = arr.view()

                arr_helper.shape = (rawchunks[0],rawchunks[1],self.ndays)

                blks = itertools.product(range(0,rawshape[0],rawchunks[0]),range(0,rawshape[1],rawchunks[1]))

                with closing(mp.Pool(processes=self.nworkers,initializer = init_worker, initargs = (shared_array,params))) as pool:

                    for b in blks:

                        for ii in range(0,arr_helper.shape[2],rawchunks[2]):

                            arr_helper[...,ii:ii+rawchunks[2]] = raw_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.startix+ii:self.startix+ii+rawchunks[2]]

                        del ii

                        # Set s values to zero
                        sarr[...] = 0

                        pool.map(execute_ws2d_vc,np.array_split(range(arr.shape[0]),self.nworkers))

                        sarr[sarr>0] = np.log10(sarr[sarr>0])

                        sgrid_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1]] = sarr[...]

                        for ii,jj in enumerate(range(self.rawoffset,arr_helper.shape[2],t_interval)):

                            smt_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.smtoffset+ii] = arr_helper[...,jj].round()

                        del ii,jj

                        bar.next()
                    bar.finish()

            else:

                arr = np.zeros((rawchunks[0]*rawchunks[1],self.ndays),dtype='float32')
                wts = arr.copy()
                sarr = np.zeros((rawchunks[0]*rawchunks[1]),dtype='float32')

                arr_helper = arr.view()
                arr_helper.shape = (rawchunks[0],rawchunks[1],self.ndays)

                blks = itertools.product(range(0,rawshape[0],rawchunks[0]),range(0,rawshape[1],rawchunks[1]))

                for b in blks:

                    for ii in range(0,arr_helper.shape[2],rawchunks[2]):

                        arr_helper[...,ii:ii+rawchunks[2]] = raw_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.startix+ii:self.startix+ii+rawchunks[2]]

                    wts[...] = (arr != nodata) * 1

                    sarr[...] = 0

                    for r in range(arr.shape[0]):
                        if wts[r,...].sum().item() != 0.0:
                            arr[r,...], sarr[r] = ws2d_vc(arr[r,...],w = wts[r,...],llas = srange)

                    sarr[sarr>0] = np.log10(sarr[sarr>0])

                    sgrid_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1]] = sarr.reshape(rawchunks[0],rawchunks[1])

                    for ii,jj in enumerate(range(self.rawoffset,arr_helper.shape[2],t_interval)):

                        smt_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.smtoffset+ii] = arr_helper[...,jj].round()

                    del ii,jj

                    bar.next()
                bar.finish()

    def ws2d_vc_asy(self,srange,p):
        '''Apply asymmetric whittaker smoother with V-curve optimization of s to data.

        Args:
            srange (arr): Float32 array of s-values to apply
            p (float): Percentile value
        '''

        with h5py.File(self.rawfile,'r') as rawh5, h5py.File(self.outname,'r+') as smth5:

            raw_ds = rawh5.get('data')
            raw_dts = rawh5.get('dates')
            smt_ds = smth5.get('data')
            smt_dts = smth5.get('dates')
            sgrid_ds = smth5.get('sgrid')

            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks

            nodata = raw_ds.attrs['nodata']
            t_interval = smt_ds.attrs['temporalresolution']

            # Store run parameters for infotool

            smt_ds.attrs['processingtimestamp'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            smt_ds.attrs['lastrun'] = "asymmetric V-curve optimization of s - p = {}".format(p)
            smt_ds.attrs['pvalue'] = p

            # Check if file needs to be resized

            dates_check = [self.rawdaily[ix] for ix in range(0,len(self.rawdaily),t_interval)]

            if len(dates_check) > smt_dts.shape[0]:
                smt_dts.resize((len(dates_check),))
                smt_dts[...] = [x.encode("ascii", "ignore") for x in dates_check]
                smt_ds.resize((smt_ds.shape[0],smt_ds.shape[1],len(dates_check)))

            # Calculate update index and offsets

            if not self.nupdate:

                for i,d in enumerate(self.daily):
                    if d in dates_check:
                        self.smtoffset = dates_check.index(d)
                        self.rawoffset = i
                        break
            else:
                self.rawoffset = self.daily.index(dates_check[-self.nupdate])
                self.smtoffset = len(dates_check[:-self.nupdate])


            barmax = (rawshape[0]/rawchunks[0]) * (rawshape[1]/rawchunks[1])
            bar = Bar('Processing',fill='=',max=barmax,suffix='%(percent)d%%  ')
            bar.goto(0)

            if self.parallel:

                params = init_parameters(nd=nodata,dim=(rawchunks[0]*rawchunks[1],self.ndays),srange=srange,p=p)

                shared_array = init_shared(rawchunks[0] * rawchunks[1] * self.ndays)

                params['shared_sarr'] = init_shared(rawchunks[0] * rawchunks[1])

                arr = tonumpyarray(shared_array)

                sarr = tonumpyarray(params['shared_sarr'])

                arr.shape = (rawchunks[0] * rawchunks[1],self.ndays)
                sarr.shape = (rawchunks[0], rawchunks[1])

                arr_helper = arr.view()

                arr_helper.shape = (rawchunks[0],rawchunks[1],self.ndays)

                blks = itertools.product(range(0,rawshape[0],rawchunks[0]),range(0,rawshape[1],rawchunks[1]))

                with closing(mp.Pool(processes=self.nworkers,initializer = init_worker, initargs = (shared_array,params))) as pool:

                    for b in blks:

                        for ii in range(0,arr_helper.shape[2],rawchunks[2]):

                            arr_helper[...,ii:ii+rawchunks[2]] = raw_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.startix+ii:self.startix+ii+rawchunks[2]]

                        del ii

                        # Set s values to zero
                        sarr[...] = 0

                        pool.map(execute_ws2d_vc_asy,np.array_split(range(arr.shape[0]),self.nworkers))

                        sarr[sarr>0] = np.log10(sarr[sarr>0])

                        sgrid_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1]] = sarr[...]

                        for ii,jj in enumerate(range(self.rawoffset,arr_helper.shape[2],t_interval)):

                            smt_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.smtoffset+ii] = arr_helper[...,jj].round()

                        del ii,jj

                        bar.next()
                    bar.finish()

            else:

                arr = np.zeros((rawchunks[0]*rawchunks[1],self.ndays),dtype='float32')
                wts = arr.copy()
                sarr = np.zeros((rawchunks[0]*rawchunks[1]),dtype='float32')

                arr_helper = arr.view()
                arr_helper.shape = (rawchunks[0],rawchunks[1],self.ndays)

                blks = itertools.product(range(0,rawshape[0],rawchunks[0]),range(0,rawshape[1],rawchunks[1]))

                for b in blks:

                    for ii in range(0,arr_helper.shape[2],rawchunks[2]):

                        arr_helper[...,ii:ii+rawchunks[2]] = raw_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.startix+ii:self.startix+ii+rawchunks[2]]

                    del ii

                    wts[...] = (arr != nodata) * 1

                    sarr[...] = 0

                    for r in range(arr.shape[0]):
                        if wts[r,...].sum().item() != 0.0:
                            arr[r,...], sarr[r] = ws2d_vc_asy(arr[r,...],w = wts[r,...],llas = srange,p = p)

                    sarr[sarr>0] = np.log10(sarr[sarr>0])

                    sgrid_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1]] = sarr.reshape(rawchunks[0],rawchunks[1])

                    for ii,jj in enumerate(range(self.rawoffset,arr_helper.shape[2],t_interval)):

                        smt_ds[b[0]:b[0]+rawchunks[0],b[1]:b[1]+rawchunks[1],self.smtoffset+ii] = arr_helper[...,jj].round()

                    del ii,jj

                    bar.next()
                bar.finish()

class MODIStiles:
    '''Class for MODIS tiles.

    Converts AOI coordinates to MODIS tile numbers by extracting values from MODIS_TILES.tif.
    '''

    def __init__(self,aoi):
        '''Creates MODIStiles object.

        Args:
            aoi (str): AOI coordinates, eiher LAT LON or XMIN, YMAX, XMAX, YMIN
        '''


        self.aoi = aoi

        # Load MODIS_TILES.tif from data directory
        this_dir, this_filename = os.path.split(__file__)
        ds = gdal.Open(os.path.join(this_dir, "data", "MODIS_TILES.tif"))

        # Try to except TIFF issues
        try:
            gt = ds.GetGeoTransform()

        except AttributeError:
            raise SystemExit("Could not find 'MODIS_TILES.tif' index raster. Try reinstalling the package.")

        # Indices fpr point AOI
        if len(self.aoi) is 2:

            xo = int(round((self.aoi[1]-gt[0])/gt[1]))
            yo = int(round((gt[3]-self.aoi[0])/gt[1]))

            xd = 1
            yd = 1

        # Indices for bounding box AOI
        elif len(self.aoi) is 4:

            xo = int(round((self.aoi[0]-gt[0])/gt[1]))
            yo = int(round((gt[3]-self.aoi[1])/gt[1]))

            xd = int(round((self.aoi[2] - self.aoi[0])/gt[1]))
            yd = int(round((self.aoi[1] - self.aoi[3])/gt[1]))

        # Read
        tile_extract = ds.ReadAsArray(xo,yo,xd,yd)
        ds = None

        # Tile IDs are stored as H*100+V
        tile_tmp = np.unique(tile_extract/100)
        tiles = ["{:05.2f}".format(x) for x in tile_tmp[tile_tmp != 0]]

        self.tiles = ["h{}v{}".format(*x.split('.')) for x in tiles]


class MODISmosaic:
    '''Class for mosaic of MODIS tiles.

    Moisaics tiles per Product, parameter and timestep. Enables extraction as GeoTiff.
    '''

    def __init__(self,files,datemin,datemax,global_flag):
        ''' Creates MODISmosaic object.

        Args:
            files ([str]): List of paths to files used for creating the mosaic
            datemin (str): Datestring for date of earliest mosaic (format YYYYMM)
            datemax (str): Datestring for date of latest mosaic (format YYYYMM)
            global_flag (bool): Flag if mosaic is global product
        '''

        # Regular expression for tile ID
        tile_re = re.compile('.+(h\d+v\d+).+')

        self.global_flag = global_flag

        self.tiles = [re.sub(tile_re,'\\1',os.path.basename(x)) for x in files]
        self.tiles.sort()
        self.files = files

        # Extract tile IDs
        self.h_ix = list(set([re.sub('(h\d+)(v\d+)','\\1',x) for x in self.tiles]))
        self.h_ix.sort()
        self.v_ix = list(set([re.sub('(h\d+)(v\d+)','\\2',x) for x in self.tiles]))
        self.v_ix.sort()

        # get referece tile identifiers

        ref_tile_h = min([x for x in self.tiles if min(self.h_ix) in x])
        ref_tile_v = min([x for x in self.tiles if min(self.v_ix) in x])


        # vertical reference tile is top (left)
        ref = [x for x in self.files if ref_tile_v in x][0]

        # Read metadata from HDF5
        try:

            with h5py.File(ref,'r') as h5f:
                dset = h5f.get('data')
                r,c,t = dset.shape
                self.tile_rws = r
                self.tile_cls = c
                self.datatype = dset.dtype
                gt_temp_v = dset.attrs['geotransform']
                self.pj = dset.attrs['projection']
                self.nodata = dset.attrs['nodata'].item()

                # If file is global, resolution is already in degrees, otherwhise it's resolution divided with 112000
                if self.global_flag:
                    self.resolution_degrees = dset.attrs['resolution']
                else:
                    self.resolution = dset.attrs['resolution']
                    self.resolution_degrees = self.resolution/112000


                dset = None
                self.dates = [x.decode() for x in h5f.get('dates')[...]]


            # checking referece tile for h

            ref = [x for x in self.files if ref_tile_h in x][0]

            with h5py.File(ref,'r') as h5f:
                dset = h5f.get('data')
                gt_temp_h = dset.attrs['geotransform']
                dset = None


            self.gt = [y for x in [gt_temp_h[0:3],gt_temp_v[3:6]] for y in x]


        except Exception as e:
            print('\nError reading referece file {} for mosaic! Error message: {}\n'.format(ref,e))
            raise

        # Create temporal index from dates available and min max input
        dts_dt = [fromjulian(x) for x in self.dates]
        datemin_p = datetime.datetime.strptime(datemin,'%Y%m').date()
        datemax_p = datetime.datetime.strptime(datemax,'%Y%m').date()

        self.tempIX = np.flatnonzero(np.array([x >= datemin_p and x <= datemax_p for x in dts_dt]))

    def getArray(self,dataset,ix,dt):
        '''Reads values for mosaic into array.

        Args:
            dataset (str): Defines dataset to be read from HDF5 file (default is 'data')
            ix (int): Temporal index
            dt (str): Datatype (default will be read from file)

        Returns
            Array for mosaic
        '''

        # Initialize array
        array = np.zeros(((len(self.v_ix) * self.tile_rws),len(self.h_ix) * self.tile_cls),dtype=dt)

        # read data from intersecting HDF5 files
        for h5f in self.files:

            # Extract tile ID from filename
            t_h = re.sub('.+(h\d+)(v\d+).+','\\1',os.path.basename(h5f))
            t_v = re.sub('.+(h\d+)(v\d+).+','\\2',os.path.basename(h5f))

            # Caluclate row/column offset
            xoff = self.h_ix.index(t_h) * self.tile_cls
            yoff = self.v_ix.index(t_v) * self.tile_rws

            try:

                with h5py.File(h5f,'r') as h5f_o:

                    # Dataset 'sgrid' is 2D, so no idex needed
                    if dataset == 'sgrid':
                        array[yoff:(yoff+self.tile_rws),xoff:(xoff+self.tile_cls)] = h5f_o.get(dataset)[...]
                    else:
                        array[yoff:(yoff+self.tile_rws),xoff:(xoff+self.tile_cls)] = h5f_o.get(dataset)[...,ix]

            except Exception as e:
                print('Error reading data from file {} to array! Error message {}:\n'.format(h5f,e))
                raise

        return(array)

    def getArrayGlobal(self,dataset,ix,dt):
        '''Reads values for global mosaic into array.

        Since files are global, the array will be a spatial and temporal subset rather than a mosaic.

        Args:
            dataset (str): Defines dataset to be read from HDF5 file (default is 'data')
            ix (int): Temporal index
            dt (str): Datatype (default will be read from file)

        Returns
            Array for mosaic
        '''

        array = np.zeros((self.tile_rws,self.tile_cls),dtype=dt)

        for h5f in self.files:

            try:

                with h5py.File(h5f,'r') as h5f_o:

                    if dataset == 'sgrid':
                        array[...] = h5f_o.get(dataset)[...]
                    else:
                        array[...] = h5f_o.get(dataset)[...,ix]

            except Exception as e:
                print('Error reading data from file {} to array! Error message {}:\n'.format(h5f,e))
                raise

        return(array)

    @contextmanager
    def getRaster(self,dataset,ix):
        '''Generator for mosaic raster.

        This generator can be used within a context manager and will yield an in-memory raster.

        Args:
            dataset (str): Defines dataset to be read from HDF5 file (default is 'data')
            ix (int): Temporal index
        '''

        # The datatype for sgrid is set to float32
        try:
            if dataset == 'sgrid':
                self.dt_gdal = dtype_GDNP('float32')
            else:
                self.dt_gdal = dtype_GDNP(self.datatype.name)
        except IndexError:
            print("\n\n Couldn't read data type from dataset. Using default Int16!\n")
            dt_gdal = (3,'int16')

        # Use the corresponding getArray function if global_flag
        if self.global_flag:
            array = self.getArrayGlobal(dataset,ix,self.dt_gdal[1])
        else:
            array = self.getArray(dataset,ix,self.dt_gdal[1])

        height, width = array.shape

        driver = gdal.GetDriverByName('GTiff')

        # Create in-memory dataset with virtual filename driver
        self.raster = driver.Create('/vsimem/inmem.tif', width, height, 1, self.dt_gdal[0])

        # Set metadata
        self.raster.SetGeoTransform(self.gt)
        self.raster.SetProjection(self.pj)

        rb = self.raster.GetRasterBand(1)

        rb.SetNoDataValue(self.nodata)

        # Write array
        rb.WriteArray(array)

        yield self

        # Cleanup to be exectuted when context manager closes after yield
        gdal.Unlink('/vsimem/inmem.tif')
        self.raster = None
        driver = None
        del array
