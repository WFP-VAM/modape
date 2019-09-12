"""
MODIS raw HDF5 class.

This file contains the class representing a raw MODIS HDF5 file.

Author: Valentin Pesendorfer, April 2019
"""
from __future__ import absolute_import, division, print_function

import gc
import os
from pathlib import Path
import re
import time
import traceback
import warnings

import numpy as np
try:
    import gdal
except ImportError:
    from osgeo import gdal
import h5py # pylint: disable=import-error
from modape.utils import dtype_GDNP, FileHandler

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

        # Patterns for string extraction
        ppatt = re.compile(r'M\w{6}')
        vpatt = re.compile(r'.+\.(\d{3})\..+')
        tpatt = re.compile(r'h\d+v\d+')
        dtspatt = re.compile(r'.+A(\d{7}).+')
        ptspatt = re.compile(r'.+(\d{13}).+')

        self.targetdir = Path(targetdir)
        #self.resdict = dict(zip(['250m','500m','1km','0.05_Deg'],[x/112000 for x in [250,500,1000,5600]])) ## commented for original resolution
        self.vam_product_code_dict = dict(zip(['VIM', 'VEM', 'LTD', 'LTN'],
                                              ['NDVI', 'EVI', 'LST_Day', 'LST_Night']))

        self.rawdates = [re.findall(dtspatt, x)[0] for x in files]
        self.files = [x for (y, x) in sorted(zip(self.rawdates, files))]
        self.rawdates.sort()
        self.nfiles = len(self.files)
        self.reference_file = Path(self.files[0])

        # class works only for M.D11 and M.D13 products
        if not re.match(r'M.D13\w\d', self.reference_file.name) and not re.match(r'M.D11\w\d', self.reference_file.name):
            raise SystemExit("Processing only implemented for M*D11 or M*13 products!")

        # make sure number of dates is equal to number of files, so no duplicates!
        processing_timestamps = [int(re.sub(ptspatt, '\\1', x)) for x in self.files]

        # check for duplicates
        if len(set(self.rawdates)) != self.nfiles:
            warnings.warn("Possibly duplicate files in {}! Using files with most recent processing date.".format(self.targetdir.as_posix()), Warning)

            dups = []
            dt_prev = None
            for dt in self.rawdates:
                if dt == dt_prev and dt not in dups:
                    dups.append(dt)
                dt_prev = dt

            # Loop over duplicates and exclude files based on processing timestamp
            to_pop = []
            for dup in dups:
                # max TS for duplicate
                max_ts = max([processing_timestamps[ix] for ix, x in enumerate(self.rawdates) if x == dup])

                for ix, x in enumerate(self.rawdates):
                    if x == dup and processing_timestamps[ix] != max_ts:
                        to_pop.append(ix)

            # re-set file-based instance variables after removal of duplicates
            for ix in sorted(to_pop, reverse=True):
                del self.files[ix]
                del self.rawdates[ix]
                del processing_timestamps[ix]

            self.nfiles = len(self.files)
            self.reference_file = Path(self.files[0])

            # assert that there are no duplicates now
            assert len(set(self.rawdates)) == self.nfiles

        # Open reference file
        ref = gdal.Open(self.reference_file.as_posix())

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
            self.product = [re.sub(r'M[O|Y]D', 'MXD', re.findall(ppatt, self.reference_file.name)[0])]
            self.temporalresolution = 8
            self.tshift = 8
        else:
            self.product = re.findall(ppatt, self.reference_file.name)
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
        self.outname = Path('{}/{}/{}.{}.h5'.format(
            self.targetdir.as_posix(),
            self.vam_product_code,
            '.'.join(self.product + re.findall(tpatt, self.reference_file.name) + [re.sub(vpatt, '\\1', self.reference_file.name)]),
            fname_vpc))

        self.exists = self.outname.exists()
        ref = None

    def create(self, compression='gzip', chunk=None):
        """Creates the HDF5 file.

        Args:
            compression: Compression method to be used (default = gzip)
            chunk: Number of pixels per chunk (needs to define equal sized chunks!)
        """

        ref = gdal.Open(self.reference_file.as_posix())
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
        self.outname.parent.mkdir(parents=True, exist_ok=True)

        # Create HDF5 file
        try:
            with h5py.File(self.outname.as_posix(), 'x', libver='latest') as h5f:
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
            print('\n\nError creating {}! Check input parameters (especially if compression/filter is available) and try again. Corrupt file will be removed now. \n\nError message: \n'.format(self.outname.as_posix()))
            self.outname.unlink()
            raise

    def update(self):
        """Update MODIS raw HDF5 file with raw data.

        When a new HDF5 file is created, update will also handle the first data ingest.
        """

        try:
            with h5py.File(self.outname.as_posix(), 'r+', libver='latest') as h5f:
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
                    dset.resize((dset.shape[0], dates_length))

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
            print('Error updating {}! File may be corrupt, consider creating the file from scratch, or closer investigation. \n\nError message: \n'.format(self.outname.as_posix()))
            traceback.print_exc()
            raise

    def __str__(self):
        """String to be displayed when printing an instance of the class object"""
        return 'ModisRawH5 object: {} - {} files - exists on disk: {}'.format(self.outname.as_posix(), self.nfiles, self.exists)
