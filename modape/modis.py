"""
MODIS processing chain classes.

This file contains all the classes used in the MODIS processing chain:

 - ModisQuery: query and download raw MODIS HDF files
 - ModisRawH5: HDF5 file object containing raw MODIS data
 - ModisSmoothH5: HDF5 file object containing smoothed MODIS data
 - modis_tiles: MODIS tiles over AOI
 - ModisMosaic: mosaic of multiple smooth HDF5 files

Author: Valentin Pesendorfer, April 2019
"""
from __future__ import absolute_import, division, print_function

import array
from contextlib import contextmanager
import datetime
import gc
import multiprocessing as mp
import os
import re
from subprocess import Popen, check_output
import sys
import time
import traceback
import uuid
import warnings

import numpy as np
from progress.spinner import Spinner
import requests
try:
    import gdal
except ImportError:
    from osgeo import gdal

from bs4 import BeautifulSoup # pylint: disable=import-error
import h5py # pylint: disable=import-error
from modape.utils import (SessionWithHeaderRedirection, FileHandler, DateHelper,
                          dtype_GDNP, txx, fromjulian, init_shared, tonumpyarray,
                          init_parameters, init_worker, execute_ws2d, execute_ws2d_sgrid, execute_ws2d_vc)
from modape.whittaker import lag1corr, ws2d, ws2doptv, ws2doptvp # pylint: disable=no-name-in-module

__all__ = ['ModisQuery', 'ModisRawH5', 'ModisSmoothH5', 'modis_tiles', 'ModisMosaic']

# turn off BeautifulSoup warnings
warnings.filterwarnings('ignore', category=UserWarning, module='bs4')

class ModisQuery(object):
    """Class for querying and downloading MODIS data."""

    def __init__(self, url, begindate,
                 enddate, username=None, password=None,
                 targetdir=os.getcwd(), global_flag=None,
                 aria2=False, tile_filter=None):
        """Creates a ModisQuery object.

        Args:
            url: Query URL as created by modis_download.py
            begindate: Begin of period for query as ISO 8601 date string (YYYY-MM-DD)
            enddate: End of period for query as ISO 8601 date string (YYYY-MM-DD)
            username: Earthdata username (only required for download)
            password: Earthdata password (only required for download)
            targetdir: Path to target directory for downloaded files (default cwd)
            global_flag: Boolean flag indictaing queried product is global file instead of tiled product
            aria2: Boolean flag to use aria2 for downloading instead of python's requests
            tile_filter: List of MODIS files to query and optionally download
        """

        self.query_url = url
        self.username = username
        self.password = password
        self.targetdir = targetdir
        self.files = []
        self.modis_urls = []
        self.begin = datetime.datetime.strptime(begindate, '%Y-%m-%d').date()
        self.end = datetime.datetime.strptime(enddate, '%Y-%m-%d').date()
        self.global_flag = global_flag
        self.aria2 = aria2

        # query for products using session object
        with requests.Session() as sess:

            print('Checking for MODIS products ...', end='')
            try:
                response = sess.get(self.query_url)
                self.statuscode = response.status_code
                response.raise_for_status()
            except requests.exceptions.RequestException as e:
                print(e)
                sys.exit(1)

            soup = BeautifulSoup(response.content, features='html.parser')

            # results for global products are date directories on server, so need to query separately
            if self.global_flag:
                regex = re.compile('.*.hdf$')
                dates = np.array([x.getText() for x in soup.findAll('a', href=True) if re.match(r'\d{4}\.\d{2}\.\d{2}', x.getText())])
                dates_parsed = [datetime.datetime.strptime(x, '%Y.%m.%d/').date() for x in dates]
                dates_ix = np.flatnonzero(np.array([self.begin <= x < self.end for x in dates_parsed]))

                for date_sel in dates[dates_ix]:
                    try:
                        response = sess.get(self.query_url + date_sel)
                    except requests.exceptions.RequestException as e:
                        print(e)
                        print('Error accessing {} - skipping.'.format(self.query_url + date_sel))

                    soup = BeautifulSoup(response.content, features='html.parser')
                    hrefs = soup.find_all('a', href=True)
                    hdf_file = [x.getText() for x in hrefs if re.match(regex, x.getText())]

                    try:
                        self.modis_urls.append(self.query_url + date_sel + hdf_file[0])
                    except IndexError:
                        print('No HDF file found in {} - skipping.'.format(self.query_url + date_sel))
                        continue
            else:
                regex = re.compile(r'.+(h\d+v\d+).+')
                urls = [x.getText() for x in soup.find_all('url')]

                if tile_filter:
                    tiles = [x.lower() for x in tile_filter]
                    self.modis_urls = [x for x in urls if any(t in x for t in tiles)]
                else:
                    self.modis_urls = urls

                self.tiles = list({regex.search(x).group(1) for x in self.modis_urls})

        self.results = len(self.modis_urls)
        print('... done.\n')

        # check for results
        if self.results > 0:
            print('{} results found.\n'.format(self.results))
        else:
            print('0 results found. Please check query!')

    def set_credentials(self, username, password):
        """Set Earthdata credentials.

        Sets Earthdata username and password in created ModisQuery object.

        Args:
            username (str): Earthdata username
            password (str): Earthdata password
        """

        self.username = username
        self.password = password

    def download(self):
        """Downloads MODIS products.

        Download of files found through query, Earthdata username and password required!
        """

        if self.username is None or self.password is None:
            raise ValueError('No credentials found. Please run .setCredentials(username,password)!')

        print('[{}]: Downloading products to {} ...\n'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime()), self.targetdir))

        # download using ARIA2 if True
        if self.aria2:
            # fail if ARIA2 not installed
            try:
                _ = check_output(['aria2c', '--version'])
            except:
                raise SystemExit('ARIA2 download needs ARIA2 to be available in PATH! Please make sure it\'s installed and available in PATH!')

            # if targetdir doesn't exist, create
            if not os.path.exists(self.targetdir):
                try:
                    os.mkdir(self.targetdir)
                except FileNotFoundError:
                    print('\nCould not create target directory {} (Trying to create directories recursively?)\n'.format(self.targetdir))
                    raise

            flist = self.targetdir + '/{}'.format(str(uuid.uuid4()))

            # write URLs of query resuls to disk for WGET
            with open(flist, 'w') as thefile:
                for item in self.modis_urls:
                    thefile.write('%s\n' % item)

            args = [
                'aria2c',
                '--file-allocation=none',
                '-m', '50',
                '--retry-wait', '2',
                '-c',
                '-x', '10',
                '-s', '10',
                '--http-user', self.username,
                '--http-passwd', self.password,
                '-d', self.targetdir,
            ]

            # execute subprocess
            process_output = Popen(args + ['-i', flist])
            process_output.wait()

            # remove filelist.txt if all downloads are successful
            if process_output.returncode != 0:
                print('\nError (error code {}) occured during download, please check files against MODIS URL list ({})!\n'.format(process_output.returncode, flist))
            else:
                os.remove(flist)

            self.files = [self.targetdir + os.path.basename(x) for x in self.modis_urls]

        # download with requests
        else:
            session = SessionWithHeaderRedirection(self.username, self.password)

            for ix, url in enumerate(self.modis_urls):
                print('{} of {}'.format(ix+1, self.results))
                fname = url[url.rfind('/')+1:]

                if os.path.exists('{}/{}'.format(self.targetdir, fname)):
                    print('\nSkipping {} - {} already exists in {}!\n'.format(url, fname, self.targetdir))
                    continue

                try:
                    response = session.get(url, stream=True)
                    response.raise_for_status()
                    spinner = Spinner('Downloading {} ... '.format(fname))

                    with open(fname, 'wb') as fopen:
                        for chunk in response.iter_content(chunk_size=1024*1024):
                            fopen.write(chunk)
                            spinner.next()

                    self.files = self.files + [self.targetdir + fname]
                    print(' done.\n')
                except requests.exceptions.HTTPError as e:
                    print('Error downloading {} - skipping. Error message: {}'.format(url, e))
                    continue

        print('\n[{}]: Downloading finished.'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))


class ModisRawH5(object):
    """Class for raw MODIS data collected into HDF5 file, ready for smoothing.

    MOD/MYD 13 products can be interleaved into a combined MXD.
    """

    def __init__(self, files, vam_product_code=None,
                 targetdir=os.getcwd(), interleave=False):
        """Create a ModisRawH5 class

        Args:
            files: A list of absolute paths to MODIS raw hdf files to be processed
            vam_product_code: VAM product code to be processed (default VIM/LTD)
            targetdir: Target directory for raw MODIS HDF5 file
            interleave: Boolean flag if MOD/MYD 13  products should be interleaved
        """

        self.targetdir = targetdir
        #self.resdict = dict(zip(['250m','500m','1km','0.05_Deg'],[x/112000 for x in [250,500,1000,5600]])) ## commented for original resolution
        self.vam_product_code_dict = dict(zip(['VIM', 'VEM', 'LTD', 'LTN'],
                                              ['NDVI', 'EVI', 'LST_Day', 'LST_Night']))
        self.dates_regexp = re.compile(r'.+A(\d{7}).+')
        self.rawdates = [re.findall(self.dates_regexp, x)[0] for x in files]
        self.files = [x for (y, x) in sorted(zip(self.rawdates, files))]
        self.rawdates.sort()
        self.nfiles = len(self.files)
        self.reference_file = self.files[0]
        self.reference_file_basename = os.path.basename(self.reference_file)

        # class works only for M.D11 and M.D13 products
        if not re.match(r'M.D13\w\d', self.reference_file_basename) and not re.match(r'M.D11\w\d', self.reference_file_basename):
            raise SystemExit("Processing only implemented for M*D11 or M*13 products!")

        # make sure number of dates is equal to number of files, so no duplicates!
        assert len(set(self.rawdates)) == self.nfiles, "Number of files not equal to number of derived dates - are there duplicate HDF files?"

        # Patterns for string extraction
        ppatt = re.compile(r'M\w{6}')
        vpatt = re.compile(r'.+\.(\d{3})\..+')
        tpatt = re.compile(r'h\d+v\d+')

        # Open reference file
        ref = gdal.Open(self.reference_file)

        # When no product is selected, the default is VIM and LTD
        if not vam_product_code:
            ref_sds = [x[0] for x in ref.GetSubDatasets() if self.vam_product_code_dict['VIM'] in x[0] or self.vam_product_code_dict['LTD'] in x[0]][0]
            self.vam_product_code = [key for key, value in self.vam_product_code_dict.items() if value in ref_sds][0]
            ref_sds = None
        elif vam_product_code in self.vam_product_code_dict.keys():
            self.vam_product_code = vam_product_code
        else:
            raise ValueError('VAM product code string not recognized. Available products are %s.' % [x for x in self.vam_product_code_dict.keys()])
        ref = None

        # VAM product code to be included in filename
        fname_vpc = self.vam_product_code

        # check for MOD/MYD interleaving
        if interleave and self.vam_product_code == 'VIM':
            self.product = [re.sub(r'M[O|Y]D', 'MXD', re.findall(ppatt, self.reference_file_basename)[0])]
            self.temporalresolution = 8
            self.tshift = 8
        else:
            self.product = re.findall(ppatt, self.reference_file_basename)
            if re.match(r'M[O|Y]D13\w\d', self.product[0]):
                self.temporalresolution = 16
                self.tshift = 8
            elif re.match(r'M[O|Y]D11\w\d', self.product[0]):
                self.temporalresolution = 8
                self.tshift = 4

                # LST has specific VAM product codes for the file naming
                if self.vam_product_code == 'LTD':
                    if re.match(r'MOD11\w\d', self.product[0]):
                        fname_vpc = 'TDT'
                    elif re.match(r'MYD11\w\d', self.product[0]):
                        fname_vpc = 'TDA'
                    else:
                        pass
                elif self.vam_product_code == 'LTN':
                    if re.match(r'MOD11\w\d', self.product[0]):
                        fname_vpc = 'TNT'
                    elif re.match(r'MYD11\w\d', self.product[0]):
                        fname_vpc = 'TNA'
                    else:
                        pass

        # Name of file to be created/updated
        self.outname = '{}/{}/{}.{}.h5'.format(
            self.targetdir,
            self.vam_product_code,
            '.'.join(self.product + re.findall(tpatt, self.reference_file_basename) + [re.sub(vpatt, '\\1', self.reference_file_basename)]),
            fname_vpc)

        self.exists = os.path.isfile(self.outname)
        ref = None

    def create(self, compression='gzip', chunk=None):
        """Creates the HDF5 file.

        Args:
            compression: Compression method to be used (default = gzip)
            chunk: Number of pixels per chunk (needs to define equal sized chunks!)
        """

        ref = gdal.Open(self.reference_file)
        ref_sds = [x[0] for x in ref.GetSubDatasets() if self.vam_product_code_dict[self.vam_product_code] in x[0]][0]

        # reference raster
        rst = gdal.Open(ref_sds)
        ref_sds = None
        ref = None
        nrows = rst.RasterYSize
        ncols = rst.RasterXSize

        # check chunksize
        if not chunk:
            self.chunks = ((nrows*ncols)//25, 10) # default
        else:
            self.chunks = (chunk, 10)

        # Check if chunksize is OK
        if not ((nrows*ncols)/self.chunks[0]).is_integer():
            rst = None # close ref file
            raise ValueError('\n\nChunksize must result in equal number of chunks. Please adjust chunksize!')

        self.nodata_value = int(rst.GetMetadataItem('_FillValue'))

        # Read datatype
        dt = rst.GetRasterBand(1).DataType

        # Parse datatype - on error use default Int16
        try:
            self.datatype = dtype_GDNP(dt)
        except IndexError:
            print('\n\n Couldn\'t read data type from dataset. Using default Int16!\n')
            self.datatype = (3, 'int16')

        trans = rst.GetGeoTransform()
        prj = rst.GetProjection()
        rst = None

        # Create directory if necessary
        if not os.path.exists(os.path.dirname(self.outname)):
            # Try statements caches possible error when multiple tiles are processed in parallel
            try:
                os.makedirs(os.path.dirname(self.outname))
            except FileExistsError:
                pass

        # Create HDF5 file
        try:
            with h5py.File(self.outname, 'x', libver='latest') as h5f:
                dset = h5f.create_dataset('data',
                                          shape=(nrows*ncols, self.nfiles),
                                          dtype=self.datatype[1],
                                          maxshape=(nrows*ncols, None),
                                          chunks=self.chunks,
                                          compression=compression,
                                          fillvalue=self.nodata_value)

                h5f.create_dataset('dates',
                                   shape=(self.nfiles,),
                                   maxshape=(None,),
                                   dtype='S8',
                                   compression=compression)

                dset.attrs['geotransform'] = trans
                dset.attrs['projection'] = prj
                dset.attrs['resolution'] = trans[1] # res ## commented for original resolution
                dset.attrs['nodata'] = self.nodata_value
                dset.attrs['temporalresolution'] = self.temporalresolution
                dset.attrs['tshift'] = self.tshift
                dset.attrs['RasterXSize'] = ncols
                dset.attrs['RasterYSize'] = nrows
            self.exists = True
        except:
            print('\n\nError creating {}! Check input parameters (especially if compression/filter is available) and try again. Corrupt file will be removed now. \n\nError message: \n'.format(self.outname))
            os.remove(self.outname)
            raise

    def update(self):
        """Update MODIS raw HDF5 file with raw data.

        When a new HDF5 file is created, update will also handle the first data ingest.
        """

        try:
            with h5py.File(self.outname, 'r+', libver='latest') as h5f:
                dset = h5f.get('data')
                dates = h5f.get('dates')
                self.chunks = dset.chunks
                self.nodata_value = dset.attrs['nodata'].item()
                self.ncols = dset.attrs['RasterXSize'].item()
                self.nrows = dset.attrs['RasterYSize'].item()
                self.datatype = dtype_GDNP(dset.dtype.name)
                dset.attrs['processingtimestamp'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

                # Load any existing dates and combine with new dates
                dates_combined = [x.decode() for x in dates[...] if len(x) > 0 and x.decode() not in self.rawdates]
                _ = [dates_combined.append(x) for x in self.rawdates]

                # New total temporal length
                dates_length = len(dates_combined)

                # if new total temporal length is bigger than dataset, datasets need to be resized for additional data
                if dates_length > dset.shape[1]:
                    dates.resize((dates_length,))
                    dset.resize((dset.shape[1], dates_length))

                # Sorting index to ensure temporal continuity
                sort_ix = np.argsort(dates_combined)

                # Manual garbage collect to prevent out of memory
                _ = [gc.collect() for x in range(3)]

                # preallocate array
                arr = np.zeros((self.chunks[0], dates_length), dtype=self.datatype[1])

                # Open all files and keep reference in handler
                handler = FileHandler(self.files, self.vam_product_code_dict[self.vam_product_code])
                handler.open()
                ysize = self.chunks[0]//self.ncols

                for b in range(0, dset.shape[0], self.chunks[0]):
                    yoff = b//self.ncols

                    for b1 in range(0, dates_length, self.chunks[1]):
                        arr[..., b1:b1+self.chunks[1]] = dset[b:b+self.chunks[0], b1:b1+self.chunks[1]]
                    del b1

                    for fix, f in enumerate(self.files):
                        try:
                            arr[..., dates_combined.index(self.rawdates[fix])] = handler.handles[fix].ReadAsArray(xoff=0,
                                                                                                                  yoff=yoff,
                                                                                                                  xsize=self.ncols,
                                                                                                                  ysize=ysize).flatten()
                        except AttributeError:
                            print('Error reading from {}. Using nodata ({}) value.'.format(f, self.nodata_value))
                            arr[..., dates_combined.index(self.rawdates[fix])] = self.nodata_value
                    arr = arr[..., sort_ix]

                    for b1 in range(0, dates_length, self.chunks[1]):
                        dset[b:b+self.chunks[0], b1:b1+self.chunks[1]] = arr[..., b1:b1+self.chunks[1]]
                handler.close()

                # Write back date list
                dates_combined.sort()
                dates[...] = np.array(dates_combined, dtype='S8')
        except:
            print('Error updating {}! File may be corrupt, consider creating the file from scratch, or closer investigation. \n\nError message: \n'.format(self.outname))
            traceback.print_exc()
            raise

    def __str__(self):
        """String to be displayed when printing an instance of the class object"""
        return 'ModisRawH5 object: {} - {} files - exists on disk: {}'.format(self.outname, self.nfiles, self.exists)


class ModisSmoothH5(object):
    """Class for smoothed MODIS data collected into HDF5 file."""

    def __init__(self, rawfile, startdate=None,
                 tempint=None, nsmooth=0, nupdate=0,
                 targetdir=os.getcwd(), nworkers=1):
        """Create ModisSmoothH5 object.

        Args:
            rawfile: Full path to a ModisRawH5 file
            tempint: Integer specifying temporal interpolation (default is None, so native temporal resolution)
            nsmooth: Number of raw timesteps used for smoothing (default is all)
            nupdate: Number of smoothed timesteps to be updated (default is all)
            targetdir: Path to target directory for smoothed HDF5 file
            nworkers: Number of worker processes used in parallel as integer
        """
        if nsmooth and nupdate:
            if nsmooth < nupdate:
                raise ValueError('nsmooth must be bigger or equal (>=) to nupdate!')

        self.targetdir = targetdir
        self.rawfile = rawfile
        self.nworkers = nworkers
        self.nupdate = nupdate
        self.nsmooth = nsmooth
        self.startdate = startdate

        # Get info from raw HDF5
        with h5py.File(self.rawfile, 'r') as h5f:
            dates = h5f.get('dates')
            self.rawdates = [x.decode() for x in dates[-self.nsmooth:]]

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
        """Creates smoothed HDF5 file on disk."""

        # Try reading info from raw HDF5
        try:
            with h5py.File(self.rawfile, 'r') as h5f:
                dset = h5f.get('data')
                dates = h5f.get('dates')
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
        except Exception as e:
            raise SystemExit('Error reading rawfile {}. File may be corrupt. \n\n Error message: \n\n {}'.format(self.rawfile, e))

        # Read native temporal resolution if no user input was supplied
        if not self.temporalresolution:
            self.temporalresolution = rtres

        dates = DateHelper(rawdates=self.rawdates,
                           rtres=rtres,
                           stres=self.temporalresolution,
                           start=self.startdate,
                           nupdate=self.nupdate)

        if not self.tinterpolate:
            dates.target = self.rawdates

        dates_length = len(dates.target)

        try:
            with h5py.File(self.outname, 'x', libver='latest') as h5f:
                dset = h5f.create_dataset('data',
                                          shape=(rawshape[0], dates_length),
                                          dtype=dt, maxshape=(rawshape[0], None),
                                          chunks=rawchunks,
                                          compression=cmpr,
                                          fillvalue=rnd)

                h5f.create_dataset('sgrid',
                                   shape=(nrows*ncols,),
                                   dtype='float32',
                                   maxshape=(nrows*ncols,),
                                   chunks=(rawchunks[0],),
                                   compression=cmpr)

                h5f.create_dataset('dates',
                                   shape=(dates_length,),
                                   maxshape=(None,),
                                   dtype='S8',
                                   compression=cmpr,
                                   data=np.array(dates.target, dtype='S8'))

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

    def ws2d(self, s):
        """Apply whittaker smoother with fixed s-value to data.

        Args:
            s: log10 value of s (float)
        """

        with h5py.File(self.rawfile, 'r') as rawh5, h5py.File(self.outname, 'r+') as smth5:
            raw_ds = rawh5.get('data')
            raw_dates = rawh5.get('dates')
            rtres = raw_ds.attrs['temporalresolution'].item()
            smt_ds = smth5.get('data')
            smt_dates = smth5.get('dates')
            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks
            smoothshape = smt_ds.shape
            smoothchunks = smt_ds.chunks
            nodata = raw_ds.attrs['nodata'].item()
            self.temporalresolution = smt_ds.attrs['temporalresolution'].item()
            tshift = raw_ds.attrs['tshift'].item()

            # Store run parameters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            smt_ds.attrs['lastrun'] = "fixed s: log10(sopt) = {}".format(s)
            smt_ds.attrs['log10sopt'] = s

            dates = DateHelper(rawdates=self.rawdates,
                               rtres=rtres,
                               stres=self.temporalresolution,
                               start=self.startdate,
                               nupdate=self.nupdate)

            if not self.tinterpolate:
                dates.target = self.rawdates
            dix = dates.getDIX()

            # Resize if date list is bigger than shape of smoothed data
            if len(dates.target) > smoothshape[1]:
                smt_dates.resize((len(dates.target),))
                smt_ds.resize((smoothshape[1], len(dates.target)))
                smt_dates[...] = np.array(dates.target, dtype='S8')

            # calculate offsets
            rawoffset = [x.decode() for x in raw_dates[...]].index(self.rawdates[0])
            smoothoffset = [x.decode() for x in smt_dates[...]].index(dates.target[0])

            if self.nworkers > 1:
                if self.tinterpolate:
                    shared_array_smooth = init_shared(smoothchunks[0] * len(dates.target))
                    arr_smooth = tonumpyarray(shared_array_smooth)
                    arr_smooth.shape = (smoothchunks[0], len(dates.target))
                    arr_smooth[...] = nodata
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates:
                        vector_daily[dates.daily.index((fromjulian(rdate) + datetime.timedelta(tshift)).strftime('%Y%j'))] = -1
                else:
                    vector_daily = None
                    shared_array_smooth = None
                    arr_smooth = None
                shared_array_raw = init_shared(rawchunks[0] * len(self.rawdates))

                parameters = init_parameters(rdim=(rawchunks[0], len(self.rawdates)),
                                             sdim=(smoothchunks[0], len(dates.target)),
                                             nd=nodata,
                                             s=s,
                                             shared_array_smooth=shared_array_smooth,
                                             vec_dly=vector_daily,
                                             dix=dix)

                arr_raw = tonumpyarray(shared_array_raw)
                arr_raw.shape = (rawchunks[0], len(self.rawdates))

                pool = mp.Pool(processes=self.nworkers, initializer=init_worker, initargs=(shared_array_raw, parameters))
                # load raw data
                for br in range(0, rawshape[0], rawchunks[0]):
                    for bc in range(0, arr_raw.shape[1], rawchunks[1]):
                        bco = bc + rawoffset
                        arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]
                    ndix = np.sum(arr_raw != nodata, 1) >= (arr_raw.shape[1] * 0.2) # 20%+ data
                    map_index = np.where(ndix)[0]
                    if map_index.size == 0:
                        continue #no data points, skipping to next block
                    _ = pool.map(execute_ws2d, map_index)

                    # write back data
                    if self.tinterpolate:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]
                # close pool
                pool.close()
                pool.join()

            else:
                arr_raw = np.zeros((rawchunks[0], len(self.rawdates)), dtype='double')

                # Create weights array
                wts = arr_raw.copy()

                if self.tinterpolate:
                    arr_smooth = np.zeros((smoothchunks[0], len(dates.target)), dtype='double')
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates:
                        vector_daily[dates.daily.index((fromjulian(rdate) + datetime.timedelta(tshift)).strftime('%Y%j'))] = -1
                else:
                    arr_smooth = None

                for br in range(0, rawshape[0], rawchunks[0]):
                    try:
                        arr_smooth[...] = nodata
                    except TypeError:
                        pass
                    wts[...] = 0

                    for bc in range(0, arr_raw.shape[1], rawchunks[1]):
                        bco = bc + rawoffset
                        arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]
                    wts[...] = (arr_raw != nodata)*1

                    ndix = np.sum(wts, 1) >= (arr_raw.shape[1] * 0.2) # 20%+ data
                    map_index = np.where(ndix)[0]

                    if map_index.size == 0:
                        continue #no data points, skipping to next block

                    for ix in map_index:
                        arr_raw[ix, :] = ws2d(y=arr_raw[ix, :], lmda=10**s, w=wts[ix, :])
                        if self.tinterpolate:
                            z2 = vector_daily.copy()
                            z2[z2 != nodata] = arr_raw[ix, :]
                            z2[...] = ws2d(y=z2, lmda=0.0001, w=np.array((z2 != nodata) * 1, dtype='double'))
                            arr_smooth[ix, :] = z2[dix]
                        else:
                            pass

                    # write back data
                    if self.tinterpolate:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]

    def ws2d_sgrid(self):
        """Apply whittaker smootehr with fixed s to data.

        This fixed s version reads a pixel based s value from file, so it needs
        a previous run of V-curve s-optimization.
        """

        with h5py.File(self.rawfile, 'r') as rawh5, h5py.File(self.outname, 'r+') as smth5:
            raw_ds = rawh5.get('data')
            raw_dates = rawh5.get('dates')
            rtres = raw_ds.attrs['temporalresolution'].item()
            smt_ds = smth5.get('data')
            smt_dates = smth5.get('dates')
            smt_sgrid = smth5.get('sgrid')
            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks
            smoothshape = smt_ds.shape
            smoothchunks = smt_ds.chunks
            nodata = raw_ds.attrs['nodata'].item()
            self.temporalresolution = smt_ds.attrs['temporalresolution'].item()
            tshift = raw_ds.attrs['tshift'].item()

            # Store run parameters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
            smt_ds.attrs['lastrun'] = 'fixed s from grid'

            dates = DateHelper(rawdates=self.rawdates,
                               rtres=rtres,
                               stres=self.temporalresolution,
                               start=self.startdate,
                               nupdate=self.nupdate)

            if not self.tinterpolate:
                dates.target = self.rawdates
            dix = dates.getDIX()

            # Resize if date list is bigger than shape of smoothed data
            if len(dates.target) > smoothshape[1]:
                smt_dates.resize((len(dates.target),))
                smt_ds.resize((smoothshape[1], len(dates.target)))
                smt_dates[...] = np.array(dates.target, dtype='S8')

            # calculate offsets
            rawoffset = [x.decode() for x in raw_dates[...]].index(self.rawdates[0])
            smoothoffset = [x.decode() for x in smt_dates[...]].index(dates.target[0])

            if self.nworkers > 1:
                if self.tinterpolate:
                    shared_array_smooth = init_shared(smoothchunks[0] * len(dates.target))
                    arr_smooth = tonumpyarray(shared_array_smooth)
                    arr_smooth.shape = (smoothchunks[0], len(dates.target))
                    arr_smooth[...] = nodata
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates:
                        vector_daily[dates.daily.index((fromjulian(rdate) + datetime.timedelta(tshift)).strftime('%Y%j'))] = -1
                else:
                    vector_daily = None
                    shared_array_smooth = None
                    arr_smooth = None
                shared_array_raw = init_shared(rawchunks[0] * len(self.rawdates))

                parameters = init_parameters(rdim=(rawchunks[0], len(self.rawdates)),
                                             sdim=(smoothchunks[0], len(dates.target)),
                                             nd=nodata,
                                             shared_array_smooth=shared_array_smooth,
                                             vec_dly=vector_daily,
                                             dix=dix)

                parameters['shared_array_sgrid'] = init_shared(rawchunks[0])
                arr_raw = tonumpyarray(shared_array_raw)
                arr_raw.shape = (rawchunks[0], len(self.rawdates))
                arr_sgrid = tonumpyarray(parameters['shared_array_sgrid'])

                pool = mp.Pool(processes=self.nworkers, initializer=init_worker, initargs=(shared_array_raw, parameters))
                # load raw data
                for br in range(0, rawshape[0], rawchunks[0]):
                    for bc in range(0, arr_raw.shape[1], rawchunks[1]):
                        bco = bc + rawoffset
                        arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]

                    ndix = np.sum(arr_raw != nodata, 1) >= (arr_raw.shape[1] * 0.2) # 20%+ data
                    map_index = np.where(ndix)[0]
                    if map_index.size == 0:
                        continue #no data points, skipping to next block

                    arr_sgrid[...] = smt_sgrid[br:br+rawchunks[0]]
                    _ = pool.map(execute_ws2d_sgrid, map_index)

                    # write back data
                    if self.tinterpolate:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]
                # close pool
                pool.close()
                pool.join()

            else:
                arr_raw = np.zeros((rawchunks[0], len(self.rawdates)), dtype='double')
                arr_sgrid = np.zeros((rawchunks[0],), dtype='double')

                # Create weights array
                wts = arr_raw.copy()
                if self.tinterpolate:
                    arr_smooth = np.zeros((smoothchunks[0], len(dates.target)), dtype='double')
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates:
                        vector_daily[dates.daily.index((fromjulian(rdate) + datetime.timedelta(tshift)).strftime('%Y%j'))] = -1
                else:
                    arr_smooth = None

                for br in range(0, rawshape[0], rawchunks[0]):
                    try:
                        arr_smooth[...] = nodata
                    except TypeError:
                        pass
                    wts[...] = 0

                    for bc in range(0, arr_raw.shape[1], rawchunks[1]):
                        bco = bc + rawoffset
                        arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]
                    wts[...] = (arr_raw != nodata)*1
                    ndix = np.sum(wts, 1) >= (arr_raw.shape[1] * 0.2) # 20%+ data
                    map_index = np.where(ndix)[0]

                    if map_index.size == 0:
                        continue #no data points, skipping to next block
                    arr_sgrid[...] = smt_sgrid[br:br+rawchunks[0]]

                    for ix in map_index:
                        arr_raw[ix, :] = ws2d(y=arr_raw[ix, :], lmda=10**arr_sgrid[ix], w=wts[ix, :])
                        if self.tinterpolate:
                            z2 = vector_daily.copy()
                            z2[z2 != nodata] = arr_raw[ix, :]
                            z2[...] = ws2d(y=z2, lmda=0.0001, w=np.array((z2 != nodata)*1, dtype='double'))
                            arr_smooth[ix, :] = z2[dix]
                        else:
                            pass

                    # write back data
                    if self.tinterpolate:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]

    def ws2d_vc(self, srange, p=None):
        """Apply whittaker smoother V-curve optimization of s.

        Optionally, p value can be specified to use asymmetric smoothing.

        Args:
            srange: array of s-values to apply
            p: Percentile value (float)
        """

        with h5py.File(self.rawfile, 'r') as rawh5, h5py.File(self.outname, 'r+') as smth5:
            raw_ds = rawh5.get('data')
            raw_dates = rawh5.get('dates')
            rtres = raw_ds.attrs['temporalresolution'].item()
            smt_ds = smth5.get('data')
            smt_dates = smth5.get('dates')
            smt_sgrid = smth5.get('sgrid')
            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks
            smoothshape = smt_ds.shape
            smoothchunks = smt_ds.chunks
            nodata = raw_ds.attrs['nodata'].item()
            self.temporalresolution = smt_ds.attrs['temporalresolution'].item()
            tshift = raw_ds.attrs['tshift'].item()

            # Store run vampceters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

            if p:
                smt_ds.attrs['lastrun'] = 'V-curve optimization of s with p = {}'.format(p)
                smt_ds.attrs['pvalue'] = p
            else:
                smt_ds.attrs['lastrun'] = 'V-curve optimization of s'

            dates = DateHelper(rawdates=self.rawdates,
                               rtres=rtres,
                               stres=self.temporalresolution,
                               start=self.startdate,
                               nupdate=self.nupdate)

            if not self.tinterpolate:
                dates.target = self.rawdates
            dix = dates.getDIX()

            # Resize if date list is bigger than shape of smoothed data
            if len(dates.target) > smoothshape[1]:
                smt_dates.resize((len(dates.target),))
                smt_ds.resize((smoothshape[1], len(dates.target)))
                smt_dates[...] = np.array(dates.target, dtype='S8')

            # calculate offsets
            rawoffset = [x.decode() for x in raw_dates[...]].index(self.rawdates[0])
            smoothoffset = [x.decode() for x in smt_dates[...]].index(dates.target[0])

            if self.nworkers > 1:
                if self.tinterpolate:
                    shared_array_smooth = init_shared(smoothchunks[0] * len(dates.target))
                    arr_smooth = tonumpyarray(shared_array_smooth)
                    arr_smooth.shape = (smoothchunks[0], len(dates.target))
                    arr_smooth[...] = nodata
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates:
                        vector_daily[dates.daily.index((fromjulian(rdate) + datetime.timedelta(tshift)).strftime('%Y%j'))] = -1
                else:
                    vector_daily = None
                    shared_array_smooth = None
                    arr_smooth = None
                shared_array_raw = init_shared(rawchunks[0] * len(self.rawdates))

                parameters = init_parameters(rdim=(rawchunks[0], len(self.rawdates)),
                                             sdim=(smoothchunks[0], len(dates.target)),
                                             nd=nodata,
                                             p=p,
                                             shared_array_smooth=shared_array_smooth,
                                             vec_dly=vector_daily,
                                             dix=dix,
                                             srange=srange)

                parameters['shared_array_sgrid'] = init_shared(rawchunks[0])
                arr_raw = tonumpyarray(shared_array_raw)
                arr_raw.shape = (rawchunks[0], len(self.rawdates))
                arr_sgrid = tonumpyarray(parameters['shared_array_sgrid'])

                pool = mp.Pool(processes=self.nworkers, initializer=init_worker, initargs=(shared_array_raw, parameters))
                # load raw data
                for br in range(0, rawshape[0], rawchunks[0]):
                    for bc in range(0, arr_raw.shape[1], rawchunks[1]):
                        bco = bc + rawoffset
                        arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]
                    ndix = np.sum(arr_raw != nodata, 1) >= (arr_raw.shape[1] * 0.2) # 20%+ data
                    map_index = np.where(ndix)[0]
                    if map_index.size == 0:
                        continue #no data points, skipping to next block

                    _ = pool.map(execute_ws2d_vc, map_index)

                    # write back data
                    arr_sgrid[arr_sgrid > 0] = np.log10(arr_sgrid[arr_sgrid > 0])
                    smt_sgrid[br:br+rawchunks[0]] = arr_sgrid[...]
                    arr_sgrid[...] = 0

                    if self.tinterpolate:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]
                # close pool
                pool.close()
                pool.join()

            else:
                arr_raw = np.zeros((rawchunks[0], len(self.rawdates)), dtype='double')
                arr_sgrid = np.zeros((rawchunks[0],), dtype='double')
                wts = arr_raw.copy() # Create weights array

                if self.tinterpolate:
                    arr_smooth = np.zeros((smoothchunks[0], len(dates.target)), dtype='double')
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates:
                        vector_daily[dates.daily.index((fromjulian(rdate) + datetime.timedelta(tshift)).strftime('%Y%j'))] = -1
                else:
                    arr_smooth = None
                for br in range(0, rawshape[0], rawchunks[0]):
                    try:
                        arr_smooth[...] = nodata
                    except TypeError:
                        pass
                    wts[...] = 0
                    for bc in range(0, arr_raw.shape[1], rawchunks[1]):
                        bco = bc + rawoffset
                        arr_raw[:, bc:bc+rawchunks[1]] = raw_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]]
                    wts[...] = (arr_raw != nodata)*1
                    ndix = np.sum(wts, 1) >= (arr_raw.shape[1] * 0.2) # 20%+ data
                    map_index = np.where(ndix)[0]
                    if map_index.size == 0:
                        continue #no data points, skipping to next block

                    for ix in map_index:
                        if not isinstance(srange, np.ndarray):
                            lag_correlation = lag1corr(arr_raw[ix, :-1], arr_raw[ix, 1:], nodata)
                            if lag_correlation > 0.5:
                                sr = np.linspace(-2, 1.0, 16)
                            elif lag_correlation <= 0.5:
                                sr = np.linspace(0, 3.0, 16)
                            else:
                                sr = np.linspace(-1, 1, 11)
                        else:
                            sr = srange
                        if p:
                            arr_raw[ix, :], arr_sgrid[ix] = ws2doptvp(y=arr_raw[ix, :],
                                                                      w=np.array((arr_raw[ix, :] != nodata)*1, dtype='double'),
                                                                      llas=array.array('d', sr),
                                                                      p=p)
                        else:
                            arr_raw[ix, :], arr_sgrid[ix] = ws2doptv(y=arr_raw[ix, :],
                                                                     w=np.array((arr_raw[ix, :] != nodata)*1, dtype='double'),
                                                                     llas=array.array('d', sr))

                        if self.tinterpolate:
                            z2 = vector_daily.copy()
                            z2[z2 != nodata] = arr_raw[ix, :]
                            z2[...] = ws2d(y=z2,
                                           lmda=0.0001,
                                           w=np.array((z2 != nodata)*1, dtype='double'))
                            arr_smooth[ix, :] = z2[dix]
                        else:
                            pass

                    # write back data
                    arr_sgrid[arr_sgrid > 0] = np.log10(arr_sgrid[arr_sgrid > 0])
                    smt_sgrid[br:br+rawchunks[0]] = arr_sgrid[...]
                    arr_sgrid[...] = 0

                    if self.tinterpolate:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_smooth[:, bc:bc+rawchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bc in range(0, len(dates.target), smoothchunks[1]):
                            bco = bc + smoothoffset
                            smt_ds[br:br+rawchunks[0], bco:bco+rawchunks[1]] = arr_raw[:, bc:bc+rawchunks[1]]


class ModisMosaic(object):
    """Class for mosaic of MODIS tiles.

    Moisaics tiles per Product, parameter and timestep. Enables extraction as GeoTiff.
    """

    def __init__(self, files, datemin, datemax, global_flag):
        """ Creates ModisMosaic object.

        Args:
            files: List of paths to files used for creating the mosaic
            datemin: Datestring for date of earliest mosaic (format YYYYMM)
            datemax: Datestring for date of latest mosaic (format YYYYMM)
            global_flag: Boolean flag if mosaic is global product
        """

        tile_re = re.compile(r'.+(h\d+v\d+).+') # Regular expression for tile ID
        self.global_flag = global_flag
        self.tiles = [re.sub(tile_re, '\\1', os.path.basename(x)) for x in files]
        self.tiles.sort()
        self.files = files

        # Extract tile IDs
        self.h_ix = list({re.sub(r'(h\d+)(v\d+)', '\\1', x) for x in self.tiles})
        self.h_ix.sort()
        self.v_ix = list({re.sub(r'(h\d+)(v\d+)', '\\2', x) for x in self.tiles})
        self.v_ix.sort()

        # get referece tile identifiers
        ref_tile_h = min([x for x in self.tiles if min(self.h_ix) in x])
        ref_tile_v = min([x for x in self.tiles if min(self.v_ix) in x])

        # vertical reference tile is top (left)
        ref = [x for x in self.files if ref_tile_v in x][0]

        # Read metadata from HDF5
        try:
            with h5py.File(ref, 'r') as h5f:
                dset = h5f.get('data')
                self.tile_rws = dset.attrs['RasterYSize'].item()
                self.tile_cls = dset.attrs['RasterXSize'].item()
                self.datatype = dset.dtype
                gt_temp_v = dset.attrs['geotransform']
                self.prj = dset.attrs['projection']
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

            with h5py.File(ref, 'r') as h5f:
                dset = h5f.get('data')
                gt_temp_h = dset.attrs['geotransform']
                dset = None
            self.gt = [y for x in [gt_temp_h[0:3], gt_temp_v[3:6]] for y in x]
        except Exception as e:
            print('\nError reading referece file {} for mosaic! Error message: {}\n'.format(ref, e))
            raise

        # Create temporal index from dates available and min max input
        dates_dt = [fromjulian(x) for x in self.dates]
        datemin_p = datetime.datetime.strptime(datemin, '%Y%m').date()
        datemax_p = datetime.datetime.strptime(datemax, '%Y%m').date()
        self.temp_index = np.flatnonzero(np.array([datemin_p <= x <= datemax_p for x in dates_dt]))

    def get_array(self, dataset, ix, dt):
        """Reads values for mosaic into array.

        Args:
            dataset: Defines dataset to be read from HDF5 file (default is 'data')
            ix: Temporal index
            dt: Datatype (default will be read from file)

        Returns
            Array for mosaic
        """

        # Initialize array
        tiles_array = np.zeros(((len(self.v_ix) * self.tile_rws), len(self.h_ix) * self.tile_cls), dtype=dt)

        # read data from intersecting HDF5 files
        for h5f in self.files:
            # Extract tile ID from filename
            t_h = re.sub(r'.+(h\d+)(v\d+).+', '\\1', os.path.basename(h5f))
            t_v = re.sub(r'.+(h\d+)(v\d+).+', '\\2', os.path.basename(h5f))

            # Caluclate row/column offset
            xoff = self.h_ix.index(t_h) * self.tile_cls
            yoff = self.v_ix.index(t_v) * self.tile_rws

            try:
                with h5py.File(h5f, 'r') as h5f_o:
                    # Dataset 'sgrid' is 2D, so no idex needed
                    if dataset == 'sgrid':
                        tiles_array[yoff:(yoff+self.tile_rws),
                                    xoff:(xoff+self.tile_cls)] = h5f_o.get(dataset)[...].reshape(self.tile_rws, self.tile_cls)
                    else:
                        tiles_array[yoff:(yoff+self.tile_rws),
                                    xoff:(xoff+self.tile_cls)] = h5f_o.get(dataset)[..., ix].reshape(self.tile_rws, self.tile_cls)
            except Exception as e:
                print('Error reading data from file {} to array! Error message {}:\n'.format(h5f, e))
                raise
        return tiles_array

    def get_array_global(self, dataset, ix, dt):
        """Reads values for global mosaic into array.

        Since files are global, the array will be a spatial and temporal subset rather than a mosaic.

        Args:
            dataset: Defines dataset to be read from HDF5 file (default is 'data')
            ix: Temporal index
            dt: Datatype (default will be read from file)

        Returns
            Array for mosaic
        """

        global_array = np.zeros((self.tile_rws, self.tile_cls), dtype=dt)
        for h5f in self.files:
            try:
                with h5py.File(h5f, 'r') as h5f_o:
                    if dataset == 'sgrid':
                        global_array[...] = h5f_o.get(dataset)[...].reshape(self.tile_rws, self.tile_cls)
                    else:
                        global_array[...] = h5f_o.get(dataset)[..., ix].reshape(self.tile_rws, self.tile_cls)
            except Exception as e:
                print('Error reading data from file {} to array! Error message {}:\n'.format(h5f, e))
                raise
        return global_array

    @contextmanager
    def get_raster(self, dataset, ix):
        """Generator for mosaic raster.

        This generator can be used within a context manager and will yield an in-memory raster.

        Args:
            dataset: Defines dataset to be read from HDF5 file (default is 'data')
            ix: Temporal index

        Yields:
            in-memory raster to be passed to GDAL warp
        """

        try:
            if dataset == 'sgrid':
                self.dt_gdal = dtype_GDNP('float32') # dtype for sgrid is set to float32
            else:
                self.dt_gdal = dtype_GDNP(self.datatype.name)
        except IndexError:
            print('\n\n Couldn\'t read data type from dataset. Using default Int16!\n')
            self.dt_gdal = (3, 'int16')

        # Use the corresponding getArray function if global_flag
        if self.global_flag:
            value_array = self.get_array_global(dataset, ix, self.dt_gdal[1])
        else:
            value_array = self.get_array(dataset, ix, self.dt_gdal[1])
        height, width = value_array.shape
        driver = gdal.GetDriverByName('GTiff')

        # Create in-memory dataset with virtual filename driver
        self.raster = driver.Create('/vsimem/inmem.tif', width, height, 1, self.dt_gdal[0])
        self.raster.SetGeoTransform(self.gt)
        self.raster.SetProjection(self.prj)
        rb = self.raster.GetRasterBand(1)
        rb.SetNoDataValue(self.nodata)

        # Write array
        rb.WriteArray(value_array)
        yield self

        # Cleanup to be exectuted when context manager closes after yield
        gdal.Unlink('/vsimem/inmem.tif')
        self.raster = None
        driver = None
        del value_array

def modis_tiles(aoi):
    """Function for querying MODIS tiles.

    Converts AOI coordinates to MODIS tile numbers by extracting values from MODIS_TILES.tif.

    Args:
        aoi: AOI coordinates, either LAT LON or XMIN, YMAX, XMAX, YMIN

    Returns:
        List of MODIS tile IDs intersecting AOI
    """

    # Load MODIS_TILES.tif from data directory
    this_dir, _ = os.path.split(__file__)
    ds = gdal.Open(os.path.join(this_dir, 'data', 'MODIS_TILES.tif'))

    # Try to catch TIFF issues
    try:
        gt = ds.GetGeoTransform()
    except AttributeError:
        raise SystemExit('Could not find \'MODIS_TILES.tif\' index raster. Try re-installing the package.')

    # Indices fpr point AOI
    if len(aoi) == 2:
        xo = int(round((aoi[1]-gt[0])/gt[1]))
        yo = int(round((gt[3]-aoi[0])/gt[1]))
        xd = 1
        yd = 1
    # Indices for bounding box AOI
    elif len(aoi) == 4:
        xo = int(round((aoi[0]-gt[0])/gt[1]))
        yo = int(round((gt[3]-aoi[1])/gt[1]))
        xd = int(round((aoi[2] - aoi[0])/gt[1]))
        yd = int(round((aoi[1] - aoi[3])/gt[1]))
    tile_extract = ds.ReadAsArray(xo, yo, xd, yd) # Read
    ds = None
    tile_tmp = np.unique(tile_extract/100) # Tile IDs are stored as H*100+V
    tiles = ['{:05.2f}'.format(x) for x in tile_tmp[tile_tmp != 0]]

    return ['h{}v{}'.format(*x.split('.')) for x in tiles]
