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

    def __init__(self,files,param=None,targetdir=os.getcwd()):
        '''Create a MODISrawh5 class

        Args:
            files ([str]): A list of absolute paths to MODIS raw hdf files to be processed
            param (str): VAM parameter to be processed (default VIM/LTD)
            targetdir (str): Target directory for raw MODIS HDF5 file
        '''

        self.targetdir = targetdir
        #self.resdict = dict(zip(['250m','500m','1km','0.05_Deg'],[x/112000 for x in [250,500,1000,5600]])) ## commented for original resolution
        self.paramdict = dict(zip(['VIM', 'VEM', 'LTD', 'LTN'], ['NDVI', 'EVI', 'LST_Day', 'LST_Night']))
        self.dts_regexp = re.compile(r'.+A(\d{7}).+')
        self.rawdates = [re.findall(self.dts_regexp,x)[0] for x in files]
        self.files = [x for (y,x) in sorted(zip(self.rawdates,files))]
        self.rawdates.sort()
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
            self.tshift = 8
        else:
            self.product = re.findall(ppatt,self.ref_file_basename)
            if re.match(r'M[O|Y]D13\w\d',self.product[0]):
                self.temporalresolution = 16
                self.tshift = 8

            elif re.match(r'M[O|Y]D11\w\d',self.product[0]):
                self.temporalresolution = 8
                self.tshift = 4

        # Name of file to be created/updated
        self.outname = '{}/{}/{}.{}.h5'.format(
                                    self.targetdir,
                                    self.param,
                                    '.'.join(self.product + re.findall(tpatt,self.ref_file_basename) + [re.sub(vpatt,'\\1',self.ref_file_basename)]),
                                    self.param)

        self.exists = os.path.isfile(self.outname)
        ref = None


    def create(self,compression='gzip',chunk=None):
        '''Creates the HDF5 file.

        Args:
            compression (str): Compression method to be used (default = gzip)
            chunk (int): Number of pixels per chunk (needs to define equal sized chunks!)
        '''

        #print('\nCreating file: {} ... '.format(self.outname), end='')

        ref = gdal.Open(self.ref_file)
        ref_sds = [x[0] for x in ref.GetSubDatasets() if self.paramdict[self.param] in x[0]][0]

        # reference raster
        rst = gdal.Open(ref_sds)
        ref_sds = None
        ref = None

        nrows = rst.RasterYSize
        ncols = rst.RasterXSize

        if not chunk:
            # default chunksize (nrows, cols)
            self.chunks = ((nrows*ncols)//25,10)
        else:
            self.chunks = (chunk,10)

        # Check if chunksize is OK
        try:
            assert ((nrows*ncols)/self.chunks[0]).is_integer(), "Number of chunks not equal!"
        except AssertionError:
            print('\n\nChunksize must result in equal number of chunks. Please adjust chunksize!')
            rst = None
            raise


        self.nodata_value = int(rst.GetMetadataItem('_FillValue'))

        # Read datatype
        dt = rst.GetRasterBand(1).DataType

        # Parse datatype - on error use default Int16
        try:
            self.datatype = dtype_GDNP(dt)
        except IndexError:
            print("\n\n Couldn't read data type from dataset. Using default Int16!\n")
            self.datatype = (3,'int16')

        trans = rst.GetGeoTransform()
        prj = rst.GetProjection()

        rst = None

        # Create directory if necessary
        if not os.path.exists(os.path.dirname(self.outname)):
            os.makedirs(os.path.dirname(self.outname))

        # Create HDF5 file
        try:

            with h5py.File(self.outname,'x',libver='latest') as h5f:
                dset = h5f.create_dataset('data',shape=(nrows*ncols,self.nfiles),dtype=self.datatype[1],maxshape=(nrows*ncols,None),chunks=self.chunks,compression=compression,fillvalue=self.nodata_value)
                h5f.create_dataset('dates',shape=(self.nfiles,),maxshape=(None,),dtype='S8',compression=compression)
                dset.attrs['geotransform'] = trans
                dset.attrs['projection'] = prj
                dset.attrs['resolution'] = trans[1] # res ## commented for original resolution
                dset.attrs['nodata'] = self.nodata_value
                dset.attrs['temporalresolution'] = self.temporalresolution
                dset.attrs['tshift'] = self.tshift
                dset.attrs['RasterXSize'] = ncols
                dset.attrs['RasterYSize'] = nrows

            self.exists = True
            #print('done.\n')

        except:
            print('\n\nError creating {}! Check input parameters (especially if compression/filter is available) and try again. Corrupt file will be removed now. \n\nError message: \n'.format(self.outname))
            os.remove(self.outname)
            raise

    def update(self):
        '''Ingest raw data into MODIS raw HDF5 file.

        When a new HDF5 file is created, uodate will also handle the first data ingest.
        '''

        #print('Processing MODIS files ...\n')

        try:

            with h5py.File(self.outname,'r+',libver='latest') as h5f:
                dset = h5f.get('data')
                dts  = h5f.get('dates')
                self.chunks = dset.chunks
                self.nodata_value = dset.attrs['nodata'].item()
                self.ncols = dset.attrs['RasterXSize'].item()
                self.nrows = dset.attrs['RasterYSize'].item()
                self.datatype = dtype_GDNP(dset.dtype.name)
                dset.attrs['processingtimestamp'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

                # Load any existing dates and combine with new dates
                dates_combined = [x.decode() for x in dts[...] if len(x) > 0 and x.decode() not in self.rawdates]

                [dates_combined.append(x) for x in self.rawdates]

                # New total temporal length
                n = len(dates_combined)

                # if new total temporal length is bigger than dataset, datasets need to be resized for additional data
                if n > dset.shape[1]:
                    dts.resize((n,))
                    dset.resize((dset.shape[0],n))

                # Sorting index to ensure temporal continuity
                sort_ix = np.argsort(dates_combined)

                # Manual garbage collect to prevent out of memory
                [gc.collect() for x in range(3)]


                # preallocate array
                arr = np.zeros((self.chunks[0],n),dtype=self.datatype[1])

                #bar = Bar('Processing',fill='=',max=self.nfiles,suffix='%(percent)d%%  ')

                #bar.goto(0)

                # Open all files and keep reference in handler

                handler = FileHandler(self.files,self.paramdict[self.param])

                handler.open()

                ysize = self.chunks[0]//self.ncols

                for b in range(0, dset.shape[0], self.chunks[0]):

                    yoff = b//self.ncols

                    for b1 in range(0, n, self.chunks[1]):

                        arr[..., b1:b1+self.chunks[1]] = dset[b:b+self.chunks[0], b1:b1+self.chunks[1]]

                    del b1

                    for fix,f in enumerate(self.files):

                        try:

                            arr[...,dates_combined.index(self.rawdates[fix])] = handler.handles[fix].ReadAsArray(xoff=0,yoff=yoff,xsize=self.ncols,ysize=ysize).flatten()

                        except AttributeError:

                            print('Error reading from {}. Using nodata ({}) value.'.format(f,self.nodata_value))

                            arr[...,dates_combined.index(self.rawdates[fix])] = self.nodata_value

                    arr = arr[...,sort_ix]

                    for b1 in range(0, n, self.chunks[1]):

                        dset[b:b+self.chunks[0], b1:b1+self.chunks[1]] = arr[..., b1:b1+self.chunks[1]]

                handler.close()

                # Write back date list
                dates_combined.sort()
                dts[...] = [n.encode("ascii", "ignore") for n in dates_combined]

            #print('\ndone.\n')

        except:
            print('Error updating {}! File may be corrupt, consider creating the file from scratch, or closer investigation. \n\nError message: \n'.format(self.outname))
            traceback.print_exc()
            raise

    def __str__(self):
        '''String to be displayed wen printing an instance of the class object'''
        return("MODISrawh5 object: %s - %s files - exists on disk: %s" % (self.outname, self.nfiles, self.exists))


class MODISsmth5:
    '''Class for smoothed MODIS data collected into HDF5 file.'''

    def __init__(self,rawfile,startdate=None,tempint=None,nsmooth=0,nupdate=0,targetdir=os.getcwd(),nworkers=1):
        '''Create MODISsmth5 object.

        Args:
            rawfile (str): Full path to a MODISrawh5 file
            tempint (int): Integer specifying temporal interpolation (default is None, so native temporal resolution)
            nsmooth (int): Number of raw timesteps used for smoothing (default is all)
            nupdate (int): Number of smoothed timesteps to be updated (default is all)
            targetdir (str): Path to target directory for smoothed HDF5 file
            nworkers (int): Number of worker processes used in parallel
        '''
        if nsmooth and nupdate:
            assert nsmooth >= nupdate, "nsmooth >= nupdate!!!!"

        self.targetdir = targetdir
        self.rawfile = rawfile
        self.nworkers = nworkers
        self.nupdate = nupdate
        self.nsmooth = nsmooth
        self.startdate = startdate

        # Get info from raw HDF5
        with h5py.File(self.rawfile,'r') as h5f:
            dset = h5f.get('data')
            dts = h5f.get('dates')

            rawshape = dset.shape

            self.rawdates = [x.decode() for x in dts[-self.nsmooth:]]

            # Number of timesteps for smoothing must be bigger than for updating

        # Parse tempint to get flag for filename
        try:
            txflag = txx(tempint)
        except ValueError:
            raise SystemExit('Value for temporal interpolation not valid (interger for number of days required)!')

        if txflag != 'n':
            self.tinterpolate = True
            self.temporalresolution = tempint
            if self.startdate:
                txflag = 'c'
        else:
            self.tinterpolate = False
            self.temporalresolution = None


        # Filename for smoothed HDF5
        self.outname = '{}/{}.tx{}.{}.h5'.format(
                                    self.targetdir,
                                    '.'.join(os.path.basename(rawfile).split('.')[:-2]),
                                    txflag,
                                    os.path.basename(rawfile).split('.')[-2:-1][0])

        self.exists = os.path.isfile(self.outname)

    def create(self):
        '''Creates smoothed HDF5 file on disk.'''

        # Try reading info from raw HDF5
        try:
            with h5py.File(self.rawfile,'r') as h5f:
                dset = h5f.get('data')
                dts = h5f.get('dates')
                dt = dset.dtype.name
                cmpr = dset.compression
                rawshape = dset.shape
                ncols = dset.attrs['RasterXSize']
                nrows = dset.attrs['RasterYSize']
                rawchunks = dset.chunks
                rgt = dset.attrs['geotransform']
                rpj = dset.attrs['projection']
                rres = dset.attrs['resolution']
                rnd = dset.attrs['nodata']
                rtres = dset.attrs['temporalresolution'].item()
                tshift = dset.attrs['tshift'].item()
                firstday = fromjulian(dts[0].decode())

        except Exception as e:
            raise SystemExit('Error reading rawfile {}. File may be corrupt. \n\n Error message: \n\n {}'.format(self.rawfile,e))

        # Read native temporal resolution if no user input was supplied
        if not self.temporalresolution:
            self.temporalresolution = rtres

        dates = DateHelper(rawdates=self.rawdates, rtres=rtres,stres=self.temporalresolution, start = self.startdate, tshift=tshift,nupdate=self.nupdate)

        if not self.tinterpolate:
            dates.target = self.rawdates

        n = len(dates.target)

        #print('\nCreating file: {} ... '.format(self.outname), end='')

        try:
            with h5py.File(self.outname,'x',libver='latest') as h5f:
                dset = h5f.create_dataset('data',shape=(rawshape[0],n),dtype=dt,maxshape=(rawshape[0],None),chunks=rawchunks,compression=cmpr,fillvalue=rnd)
                h5f.create_dataset('sgrid',shape=(nrows*ncols,),dtype='float32',maxshape=(nrows*ncols,),chunks=(rawchunks[0],),compression=cmpr)
                h5f.create_dataset('dates',shape=(n,),maxshape=(None,),dtype='S8',compression=cmpr,data = [x.encode('ascii','ignore') for x in dates.target])
                dset.attrs['geotransform'] = rgt
                dset.attrs['projection'] = rpj
                dset.attrs['resolution'] = rres
                dset.attrs['nodata'] = rnd
                dset.attrs['temporalresolution'] = self.temporalresolution
                dset.attrs['RasterYSize'] = nrows
                dset.attrs['RasterXSize'] = ncols

        except:
            print('\n\nError creating {}! Check input parameters (especially if compression/filter is available) and try again. Corrupt file will be removed now. \n\nError message: \n'.format(self.outname))
            os.remove(self.outname)
            raise


        self.exists = True
        #print('done.\n')

    def ws2d(self,s):

        '''Apply whittaker smoother with fixed s-value to data.

        Args:
            s (float): log10 value of s
        '''

        with h5py.File(self.rawfile,'r') as rawh5, h5py.File(self.outname,'r+') as smth5:
            raw_ds = rawh5.get('data')
            raw_dts = rawh5.get('dates')
            rtres = raw_ds.attrs['temporalresolution'].item()

            smt_ds = smth5.get('data')
            smt_dts = smth5.get('dates')

            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks

            smoothshape = smt_ds.shape
            smoothchunks = smt_ds.chunks

            nodata = raw_ds.attrs['nodata'].item()

            self.temporalresolution = smt_ds.attrs['temporalresolution'].item()
            tshift = raw_ds.attrs['tshift'].item()

            # Store run parameters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            smt_ds.attrs['lastrun'] = "fixed s: log10(sopt) = {}".format(s)
            smt_ds.attrs['log10sopt'] = s

            dates = DateHelper(rawdates=self.rawdates, rtres=rtres,stres=self.temporalresolution, start = self.startdate, tshift=tshift,nupdate=self.nupdate)
            if not self.tinterpolate:
                dates.target = self.rawdates

            dix = dates.getDIX()

            # Resize if date list is bigger than shape of smoothed data
            if len(dates.target) > smoothshape[0]:
                smt_dts.resize((len(dates.target),))
                smt_ds.resize((smoothshape[0],len(dates.target)))
                smt_dts[...] = [x.encode("ascii", "ignore") for x in dates.target]

            # calculate offsets

            rawoffset = [x.decode() for x in raw_dts[...]].index(self.rawdates[0])
            smoothoffset = [x.decode() for x in smt_dts[...]].index(dates.target[0])

            if self.nworkers > 1:

                if self.tinterpolate:

                    shared_array_smooth = init_shared(smoothchunks[0] * len(dates.target))
                    arr_smooth = tonumpyarray(shared_array_smooth)
                    arr_smooth.shape = (smoothchunks[0],len(dates.target))
                    arr_smooth[...] = nodata


                    vec_dly = dates.getDV(nodata)

                    # Shift for interpolation
                    for d in self.rawdates:
                        vec_dly[dates.daily.index((fromjulian(d) + datetime.timedelta(tshift)).strftime('%Y%j'))] = 0

                else:
                    vec_dly = None
                    shared_array_smooth = None
                    arr_smooth = None

                shared_array_raw = init_shared(rawchunks[0] * len(self.rawdates))
                parameters = init_parameters(rdim=(rawchunks[0], len(self.rawdates)),sdim=(smoothchunks[0], len(dates.target)), nd=nodata, s=s, shared_array_smooth=shared_array_smooth, vec_dly=vec_dly, dix=dix)

                arr_raw = tonumpyarray(shared_array_raw)

                arr_raw.shape = (rawchunks[0], len(self.rawdates))

                with closing(mp.Pool(processes=self.nworkers,initializer = init_worker, initargs = (shared_array_raw,parameters))) as pool:

                    # load raw data
                    for br in range(0,rawshape[0],rawchunks[0]):

                        for bc in range(0,len(self.rawdates),rawchunks[1]):
                            bco = bc + rawoffset

                            arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]

                        ndix = np.sum(arr_raw!=-3000,1)>0 #70
                        mapIX = np.where(ndix)[0]

                        if len(mapIX) == 0:
                            #no data points, skipping to next block
                            continue

                        res = pool.map(execute_ws2d,mapIX)

                        # write back data
                        if self.tinterpolate:

                            for bc in range(0,len(dates.target),smoothchunks[1]):
                                bco = bc + smoothoffset

                                smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]

                            arr_smooth[...] = nodata

                        else:

                            for bc in range(0,len(dates.target),smoothchunks[1]):
                                bco = bc + smoothoffset

                                smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]

            else:

                arr_raw = np.zeros((rawchunks[0], len(self.rawdates)),dtype='float32')

                # Create weights array
                wts = arr_raw.copy()

                if self.tinterpolate:

                    arr_smooth = np.zeros((smoothchunks[0],len(dates.target)),dtype='float32')

                    vec_dly = dates.getDV(nodata)

                    # Shift for interpolation
                    for d in self.rawdates:
                        vec_dly[dates.daily.index((fromjulian(d) + datetime.timedelta(tshift)).strftime('%Y%j'))] = 0

                else:
                    arr_smooth = None

                for br in range(0,rawshape[0],rawchunks[0]):

                    arr_smooth[...] = nodata
                    wts[...] = 0

                    for bc in range(0,len(self.rawdates),rawchunks[1]):
                        bco = bc + rawoffset

                        arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]

                    wts[...] = (arr_raw != nodata) * 1

                    ndix = np.sum(wts,1)>0 #70
                    mapIX = np.where(ndix)[0]

                    if len(mapIX) == 0:
                        #no data points, skipping to next block
                        continue

                    for r in mapIX:

                        arr_raw[r,:] = ws2d(y = arr_raw[r,:],lmda = 10**s, w = wts[r,:])

                        if self.tinterpolate:

                            z2 = vec_dly.copy()
                            z2[z2 != nodata] = arr_raw[r,:]
                            z2[...] = ws2d(y = z2, lmda = 0.0001, w = np.array((z2 != nodata) * 1,dtype='float32'))
                            arr_smooth[r,:] = z2[dix]

                        else:
                            pass


                    # write back data
                    if self.tinterpolate:

                        for bc in range(0,len(dates.target),smoothchunks[1]):
                            bco = bc + smoothoffset

                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]

                        arr_smooth[...] = nodata

                    else:

                        for bc in range(0,len(dates.target),smoothchunks[1]):
                            bco = bc + smoothoffset

                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]

    def ws2d_sgrid(self):

        '''Apply whittaker smootehr with fixed s to data.

        This fixed s version reads a pixel based s value from file, so it needs
        a previous run of V-curve s-optimization.
        '''

        with h5py.File(self.rawfile,'r') as rawh5, h5py.File(self.outname,'r+') as smth5:
            raw_ds = rawh5.get('data')
            raw_dts = rawh5.get('dates')
            rtres = raw_ds.attrs['temporalresolution'].item()

            smt_ds = smth5.get('data')
            smt_dts = smth5.get('dates')
            smt_sgrid = smth5.get('sgrid')

            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks

            smoothshape = smt_ds.shape
            smoothchunks = smt_ds.chunks

            nodata = raw_ds.attrs['nodata'].item()

            self.temporalresolution = smt_ds.attrs['temporalresolution'].item()
            tshift = raw_ds.attrs['tshift'].item()

            # Store run parameters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            smt_ds.attrs['lastrun'] = "fixed s from grid"

            dates = DateHelper(rawdates=self.rawdates, rtres=rtres,stres=self.temporalresolution, start = self.startdate, tshift=tshift,nupdate=self.nupdate)

            if not self.tinterpolate:
                dates.target = self.rawdates

            dix = dates.getDIX()

            # Resize if date list is bigger than shape of smoothed data
            if len(dates.target) > smoothshape[0]:
                smt_dts.resize((len(dates.target),))
                smt_ds.resize((smoothshape[0],len(dates.target)))
                smt_dts[...] = [x.encode("ascii", "ignore") for x in dates.target]

            # calculate offsets

            rawoffset = [x.decode() for x in raw_dts[...]].index(self.rawdates[0])
            smoothoffset = [x.decode() for x in smt_dts[...]].index(dates.target[0])

            if self.nworkers > 1:

                if self.tinterpolate:

                    shared_array_smooth = init_shared(smoothchunks[0] * len(dates.target))
                    arr_smooth = tonumpyarray(shared_array_smooth)
                    arr_smooth.shape = (smoothchunks[0],len(dates.target))
                    arr_smooth[...] = nodata


                    vec_dly = dates.getDV(nodata)

                    # Shift for interpolation
                    for d in self.rawdates:
                        vec_dly[dates.daily.index((fromjulian(d) + datetime.timedelta(tshift)).strftime('%Y%j'))] = 0

                else:
                    vec_dly = None
                    shared_array_smooth = None
                    arr_smooth = None

                shared_array_raw = init_shared(rawchunks[0] * len(self.rawdates))
                parameters = init_parameters(rdim=(rawchunks[0], len(self.rawdates)),sdim=(smoothchunks[0], len(dates.target)), nd=nodata, shared_array_smooth=shared_array_smooth, vec_dly=vec_dly, dix=dix)

                parameters['shared_array_sgrid'] = init_shared(rawchunks[0])

                arr_raw = tonumpyarray(shared_array_raw)
                arr_raw.shape = (rawchunks[0], len(self.rawdates))

                arr_sgrid = tonumpyarray(parameters['shared_array_sgrid'])

                with closing(mp.Pool(processes=self.nworkers,initializer = init_worker, initargs = (shared_array_raw,parameters))) as pool:

                    # load raw data
                    for br in range(0,rawshape[0],rawchunks[0]):

                        for bc in range(0,len(self.rawdates),rawchunks[1]):
                            bco = bc + rawoffset

                            arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]

                        ndix = np.sum(arr_raw!=-3000,1)>0 #70
                        mapIX = np.where(ndix)[0]

                        if len(mapIX) == 0:
                            #no data points, skipping to next block
                            continue

                        arr_sgrid[...] = smt_sgrid[br:br+rawchunks[0]]

                        res = pool.map(execute_ws2d_sgrid,mapIX)

                        # write back data
                        if self.tinterpolate:

                            for bc in range(0,len(dates.target),smoothchunks[1]):
                                bco = bc + smoothoffset

                                smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]

                            arr_smooth[...] = nodata

                        else:

                            for bc in range(0,len(dates.target),smoothchunks[1]):
                                bco = bc + smoothoffset

                                smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]

            else:

                arr_raw = np.zeros((rawchunks[0], len(self.rawdates)),dtype='float32')
                arr_sgrid = np.zeros((rawchunks[0],),dtype='float32')

                # Create weights array
                wts = arr_raw.copy()

                if self.tinterpolate:

                    arr_smooth = np.zeros((smoothchunks[0],len(dates.target)),dtype='float32')

                    vec_dly = dates.getDV(nodata)

                    # Shift for interpolation
                    for d in self.rawdates:
                        vec_dly[dates.daily.index((fromjulian(d) + datetime.timedelta(tshift)).strftime('%Y%j'))] = 0

                else:
                    arr_smooth = None

                for br in range(0,rawshape[0],rawchunks[0]):

                    arr_smooth[...] = nodata
                    wts[...] = 0

                    for bc in range(0,len(self.rawdates),rawchunks[1]):
                        bco = bc + rawoffset

                        arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]

                    wts[...] = (arr_raw != nodata) * 1

                    ndix = np.sum(wts,1)>0 #70
                    mapIX = np.where(ndix)[0]

                    if len(mapIX) == 0:
                        #no data points, skipping to next block
                        continue

                    arr_sgrid[...] = smt_sgrid[br:br+rawchunks[0]]

                    for r in mapIX:

                        arr_raw[r,:] = ws2d(y = arr_raw[r,:],lmda = 10**arr_sgrid[r], w = wts[r,:])

                        if self.tinterpolate:

                            z2 = vec_dly.copy()
                            z2[z2 != nodata] = arr_raw[r,:]
                            z2[...] = ws2d(y = z2, lmda = 0.0001, w = np.array((z2 != nodata) * 1,dtype='float32'))
                            arr_smooth[r,:] = z2[dix]

                        else:
                            pass


                    # write back data
                    if self.tinterpolate:

                        for bc in range(0,len(dates.target),smoothchunks[1]):
                            bco = bc + smoothoffset

                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]
                        arr_smooth[...] = nodata

                    else:

                        for bc in range(0,len(dates.target),smoothchunks[1]):
                            bco = bc + smoothoffset

                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]


    def ws2d_vc(self,srange,p=None):

        '''Apply whittaker smoother V-curve optimization of s.

        Optionally, a p value can be specified to use asymmetric smoothing.

        Args:
            srange (arr): Float32 array of s-values to apply
            p (float): Percentile value
        '''

        with h5py.File(self.rawfile,'r') as rawh5, h5py.File(self.outname,'r+') as smth5:
            raw_ds = rawh5.get('data')
            raw_dts = rawh5.get('dates')
            rtres = raw_ds.attrs['temporalresolution'].item()

            smt_ds = smth5.get('data')
            smt_dts = smth5.get('dates')
            smt_sgrid = smth5.get('sgrid')

            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks

            smoothshape = smt_ds.shape
            smoothchunks = smt_ds.chunks

            nodata = raw_ds.attrs['nodata'].item()

            self.temporalresolution = smt_ds.attrs['temporalresolution'].item()
            tshift = raw_ds.attrs['tshift'].item()

            # Store run parameters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

            if p:
                smt_ds.attrs['lastrun'] = "V-curve optimization of s with p = {}".format(p)
                smt_ds.attrs['pvalue'] = p
            else:
                smt_ds.attrs['lastrun'] = "V-curve optimization of s"

            dates = DateHelper(rawdates=self.rawdates, rtres=rtres,stres=self.temporalresolution, start = self.startdate,tshift=tshift,nupdate=self.nupdate)

            if not self.tinterpolate:
                dates.target = self.rawdates

            dix = dates.getDIX()

            # Resize if date list is bigger than shape of smoothed data
            if len(dates.target) > smoothshape[0]:
                smt_dts.resize((len(dates.target),))
                smt_ds.resize((smoothshape[0],len(dates.target)))
                smt_dts[...] = [x.encode("ascii", "ignore") for x in dates.target]

            # calculate offsets

            rawoffset = [x.decode() for x in raw_dts[...]].index(self.rawdates[0])
            smoothoffset = [x.decode() for x in smt_dts[...]].index(dates.target[0])

            if self.nworkers > 1:

                if self.tinterpolate:

                    shared_array_smooth = init_shared(smoothchunks[0] * len(dates.target))
                    arr_smooth = tonumpyarray(shared_array_smooth)
                    arr_smooth.shape = (smoothchunks[0],len(dates.target))
                    arr_smooth[...] = nodata


                    vec_dly = dates.getDV(nodata)

                    # Shift for interpolation
                    for d in self.rawdates:
                        vec_dly[dates.daily.index((fromjulian(d) + datetime.timedelta(tshift)).strftime('%Y%j'))] = 0

                else:
                    vec_dly = None
                    shared_array_smooth = None
                    arr_smooth = None

                shared_array_raw = init_shared(rawchunks[0] * len(self.rawdates))
                parameters = init_parameters(rdim=(rawchunks[0], len(self.rawdates)),sdim=(smoothchunks[0], len(dates.target)), nd=nodata, p=p, shared_array_smooth=shared_array_smooth, vec_dly=vec_dly, dix=dix, srange=srange)

                parameters['shared_array_sgrid'] = init_shared(rawchunks[0])

                arr_raw = tonumpyarray(shared_array_raw)
                arr_raw.shape = (rawchunks[0], len(self.rawdates))

                arr_sgrid = tonumpyarray(parameters['shared_array_sgrid'])

                with closing(mp.Pool(processes=self.nworkers,initializer = init_worker, initargs = (shared_array_raw,parameters))) as pool:

                    # load raw data
                    for br in range(0,rawshape[0],rawchunks[0]):

                        for bc in range(0,len(self.rawdates),rawchunks[1]):
                            bco = bc + rawoffset

                            arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]

                        ndix = np.sum(arr_raw!=-3000,1)>0 #70
                        mapIX = np.where(ndix)[0]

                        if len(mapIX) == 0:
                            #no data points, skipping to next block
                            continue

                        res = pool.map(execute_ws2d_vc,mapIX)

                        # write back data

                        arr_sgrid[arr_sgrid>0] = np.log10(arr_sgrid[arr_sgrid>0])

                        smt_sgrid[br:br+rawchunks[0]] = arr_sgrid[...]
                        arr_sgrid[...] = 0

                        if self.tinterpolate:

                            for bc in range(0,len(dates.target),smoothchunks[1]):
                                bco = bc + smoothoffset

                                smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]

                            arr_smooth[...] = nodata

                        else:

                            for bc in range(0,len(dates.target),smoothchunks[1]):
                                bco = bc + smoothoffset

                                smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]

            else:

                arr_raw = np.zeros((rawchunks[0], len(self.rawdates)),dtype='float32')
                arr_sgrid = np.zeros((rawchunks[0],),dtype='float32')

                # Create weights array
                wts = arr_raw.copy()

                if self.tinterpolate:

                    arr_smooth = np.zeros((smoothchunks[0],len(dates.target)),dtype='float32')

                    vec_dly = dates.getDV(nodata)

                    # Shift for interpolation
                    for d in self.rawdates:
                        vec_dly[dates.daily.index((fromjulian(d) + datetime.timedelta(tshift)).strftime('%Y%j'))] = 0

                else:
                    arr_smooth = None

                for br in range(0,rawshape[0],rawchunks[0]):

                    arr_smooth[...] = nodata
                    wts[...] = 0

                    for bc in range(0,len(self.rawdates),rawchunks[1]):
                        bco = bc + rawoffset

                        arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]

                    wts[...] = (arr_raw != nodata) * 1

                    ndix = np.sum(wts,1)>0 #70
                    mapIX = np.where(ndix)[0]

                    if len(mapIX) == 0:
                        #no data points, skipping to next block
                        continue

                    for r in mapIX:

                        if p:
                            arr_raw[r,:] , arr_sgrid[r] = ws2d_vc_asy(y = arr_raw[r,:], w = np.array((arr_raw[r,:] != nodata) * 1,dtype='float32'), llas = array.array('f',srange),p=p)
                        else:
                            arr_raw[r,:] , arr_sgrid[r] = ws2d_vc(y = arr_raw[r,:], w = np.array((arr_raw[r,:] != nodata) * 1,dtype='float32'), llas = array.array('f',srange))

                        if self.tinterpolate:

                            z2 = vec_dly.copy()
                            z2[z2 != nodata] = arr_raw[r,:]
                            z2[...] = ws2d(y = z2, lmda = 0.0001, w = np.array((z2 != nodata) * 1,dtype='float32'))
                            arr_smooth[r,:] = z2[dix]

                        else:
                            pass

                    # write back data

                    arr_sgrid[arr_sgrid>0] = np.log10(arr_sgrid[arr_sgrid>0])

                    smt_sgrid[br:br+rawchunks[0]] = arr_sgrid[...]
                    arr_sgrid[...] = 0

                    if self.tinterpolate:

                        for bc in range(0,len(dates.target),smoothchunks[1]):
                            bco = bc + smoothoffset

                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]

                        arr_smooth[...] = nodata


                    else:

                        for bc in range(0,len(dates.target),smoothchunks[1]):
                            bco = bc + smoothoffset

                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]



    def ws2d_vcOpt(self,srange,p=None):

        '''Apply whittaker smoother V-curve optimization of s.

        This current implementation runs the optimization two times, using the first sOpt to further constrain the
        srange.

        Optionally, a p value can be specified to use asymmetric smoothing.

        Args:
            srange (arr): Float32 array of s-values to apply
            p (float): Percentile value
        '''

        with h5py.File(self.rawfile,'r') as rawh5, h5py.File(self.outname,'r+') as smth5:
            raw_ds = rawh5.get('data')
            raw_dts = rawh5.get('dates')
            rtres = raw_ds.attrs['temporalresolution'].item()

            smt_ds = smth5.get('data')
            smt_dts = smth5.get('dates')
            smt_sgrid = smth5.get('sgrid')

            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks

            smoothshape = smt_ds.shape
            smoothchunks = smt_ds.chunks

            nodata = raw_ds.attrs['nodata'].item()

            self.temporalresolution = smt_ds.attrs['temporalresolution'].item()
            tshift = raw_ds.attrs['tshift'].item()

            # Store run parameters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())

            if p:
                smt_ds.attrs['lastrun'] = "V-curve 2-step optimization of s with p = {}".format(p)
                smt_ds.attrs['pvalue'] = p
            else:
                smt_ds.attrs['lastrun'] = "V-curve 2-step optimization of s"

            dates = DateHelper(rawdates=self.rawdates, rtres=rtres,stres=self.temporalresolution, start = self.startdate,tshift=tshift,nupdate=self.nupdate)

            if not self.tinterpolate:
                dates.target = self.rawdates

            dix = dates.getDIX()

            # Resize if date list is bigger than shape of smoothed data
            if len(dates.target) > smoothshape[0]:
                smt_dts.resize((len(dates.target),))
                smt_ds.resize((smoothshape[0],len(dates.target)))
                smt_dts[...] = [x.encode("ascii", "ignore") for x in dates.target]

            # calculate offsets

            rawoffset = [x.decode() for x in raw_dts[...]].index(self.rawdates[0])
            smoothoffset = [x.decode() for x in smt_dts[...]].index(dates.target[0])

            if self.nworkers > 1:

                if self.tinterpolate:

                    shared_array_smooth = init_shared(smoothchunks[0] * len(dates.target))
                    arr_smooth = tonumpyarray(shared_array_smooth)
                    arr_smooth.shape = (smoothchunks[0],len(dates.target))
                    arr_smooth[...] = nodata


                    vec_dly = dates.getDV(nodata)

                    # Shift for interpolation
                    for d in self.rawdates:
                        vec_dly[dates.daily.index((fromjulian(d) + datetime.timedelta(tshift)).strftime('%Y%j'))] = 0

                else:
                    vec_dly = None
                    shared_array_smooth = None
                    arr_smooth = None

                shared_array_raw = init_shared(rawchunks[0] * len(self.rawdates))
                parameters = init_parameters(rdim=(rawchunks[0], len(self.rawdates)),sdim=(smoothchunks[0], len(dates.target)), nd=nodata, p=p, shared_array_smooth=shared_array_smooth, vec_dly=vec_dly, dix=dix, srange=srange)

                parameters['shared_array_sgrid'] = init_shared(rawchunks[0])

                arr_raw = tonumpyarray(shared_array_raw)
                arr_raw.shape = (rawchunks[0], len(self.rawdates))

                arr_sgrid = tonumpyarray(parameters['shared_array_sgrid'])

                with closing(mp.Pool(processes=self.nworkers,initializer = init_worker, initargs = (shared_array_raw,parameters))) as pool:

                    # load raw data
                    for br in range(0,rawshape[0],rawchunks[0]):

                        for bc in range(0,len(self.rawdates),rawchunks[1]):
                            bco = bc + rawoffset

                            arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]

                        ndix = np.sum(arr_raw!=-3000,1)>0 #70
                        mapIX = np.where(ndix)[0]

                        if len(mapIX) == 0:
                            #no data points, skipping to next block
                            continue

                        res = pool.map(execute_ws2d_vcOpt,mapIX)

                        # write back data

                        arr_sgrid[arr_sgrid>0] = np.log10(arr_sgrid[arr_sgrid>0])

                        smt_sgrid[br:br+rawchunks[0]] = arr_sgrid[...]
                        arr_sgrid[...] = 0

                        if self.tinterpolate:

                            for bc in range(0,len(dates.target),smoothchunks[1]):
                                bco = bc + smoothoffset

                                smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]

                            arr_smooth[...] = nodata

                        else:

                            for bc in range(0,len(dates.target),smoothchunks[1]):
                                bco = bc + smoothoffset

                                smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]

            else:

                arr_raw = np.zeros((rawchunks[0], len(self.rawdates)),dtype='float32')
                arr_sgrid = np.zeros((rawchunks[0],),dtype='float32')

                # Create weights array
                wts = arr_raw.copy()

                if self.tinterpolate:

                    arr_smooth = np.zeros((smoothchunks[0],len(dates.target)),dtype='float32')

                    vec_dly = dates.getDV(nodata)

                    # Shift for interpolation
                    for d in self.rawdates:
                        vec_dly[dates.daily.index((fromjulian(d) + datetime.timedelta(tshift)).strftime('%Y%j'))] = 0

                else:
                    arr_smooth = None

                for br in range(0,rawshape[0],rawchunks[0]):

                    arr_smooth[...] = nodata
                    wts[...] = 0

                    for bc in range(0,len(self.rawdates),rawchunks[1]):
                        bco = bc + rawoffset

                        arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]

                    wts[...] = (arr_raw != nodata) * 1

                    ndix = np.sum(wts,1)>0 #70
                    mapIX = np.where(ndix)[0]

                    if len(mapIX) == 0:
                        #no data points, skipping to next block
                        continue

                    for r in mapIX:

                        z , lopt = ws2d_vc(y = arr_raw[r,:], w = np.array((arr_raw[r,:] != nodata) * 1,dtype='float32'), llas = array.array('f',srange))

                        srange_lim = srange[srange <= np.log10(lopt)]

                        if len(srange_lim)==1:
                            srange_lim = np.concatenate([srange_lim-0.2,srange_lim])

                        if p:
                            arr_raw[r,:] , arr_sgrid[r] = ws2d_vc_asy(y = arr_raw[r,:], w = np.array((arr_raw[r,:] != nodata) * 1,dtype='float32'), llas = array.array('f',srange_lim),p=p)
                        else:
                            arr_raw[r,:] , arr_sgrid[r] = ws2d_vc(y = arr_raw[r,:], w = np.array((arr_raw[r,:] != nodata) * 1,dtype='float32'), llas = array.array('f',srange_lim))

                        if self.tinterpolate:

                            z2 = vec_dly.copy()
                            z2[z2 != nodata] = arr_raw[r,:]
                            z2[...] = ws2d(y = z2, lmda = 0.0001, w = np.array((z2 != nodata) * 1,dtype='float32'))
                            arr_smooth[r,:] = z2[dix]

                        else:
                            pass

                    # write back data

                    arr_sgrid[arr_sgrid>0] = np.log10(arr_sgrid[arr_sgrid>0])

                    smt_sgrid[br:br+rawchunks[0]] = arr_sgrid[...]
                    arr_sgrid[...] = 0

                    if self.tinterpolate:

                        for bc in range(0,len(dates.target),smoothchunks[1]):
                            bco = bc + smoothoffset

                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]

                        arr_smooth[...] = nodata


                    else:

                        for bc in range(0,len(dates.target),smoothchunks[1]):
                            bco = bc + smoothoffset

                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]


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
                self.tile_rws = dset.attrs['RasterYSize'].item()
                self.tile_cls = dset.attrs['RasterXSize'].item()
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
                        array[yoff:(yoff+self.tile_rws),xoff:(xoff+self.tile_cls)] = h5f_o.get(dataset)[...].reshape(self.tile_rws,self.tile_cls)
                    else:
                        array[yoff:(yoff+self.tile_rws),xoff:(xoff+self.tile_cls)] = h5f_o.get(dataset)[...,ix].reshape(self.tile_rws,self.tile_cls)

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
                        array[...] = h5f_o.get(dataset)[...].reshape(self.tile_rws,self.tile_cls)
                    else:
                        array[...] = h5f_o.get(dataset)[...,ix].reshape(self.tile_rws,self.tile_cls)

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
