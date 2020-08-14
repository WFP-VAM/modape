"""
MODIS raw HDF5 class.

This file contains the class representing a raw MODIS HDF5 file.

Author: Valentin Pesendorfer, July 2020
"""
#pylint: disable=E0401
from datetime import datetime
import gc
import logging
from pathlib import Path
import re
from typing import List, Tuple
import warnings

import h5py
import numpy as np
from modape.constants import (
    LST_NAME_LUD,
    REGEX_PATTERNS,
    TEMPORAL_DICT,
    VAM_PRODUCT_CODES,
)
from modape.exceptions import HDF5CreationError, HDF5WriteError
from modape.modis.io import HDF5Base, HDFHandler
from modape.utils import fromjulian

log = logging.getLogger(__name__)

class ModisRawH5(HDF5Base):
    """Class for raw MODIS data collected into HDF5 file, ready for smoothing.

    MOD/MYD 13 products can be interleaved into a combined MXD.
    """

    def __init__(self,
                 files: List[str],
                 targetdir: str,
                 vam_product_code: str = None,
                 interleave: bool = False) -> None:
        """Create a ModisRawH5 class

        Args:
            files: A list of absolute paths to MODIS raw hdf files to be processed
            vam_product_code: VAM product code to be processed (default VIM/LTD)
            targetdir: Target directory for raw MODIS HDF5 file
            interleave: Boolean flag if MOD/MYD 13  products should be interleaved
        """

        self.targetdir = Path(targetdir)

        # check if all files provided are either VIM or TDA
        for file in files:
            if not re.match(REGEX_PATTERNS['VIMLST'], file.split('/')[-1]):
                log.error("File %s not NDVI or LST", file)
                raise ValueError("MODIS collect processing only implemented for M{O|Y}D {11,13} products!")

        self.rawdates = [re.findall(REGEX_PATTERNS["date"], x)[0] for x in files]
        self.files = [x for (y, x) in sorted(zip(self.rawdates, files))]
        self.rawdates.sort()

        # extract product patterns
        products = []
        for file_tmp in self.files:
            products.append(
                re.findall(REGEX_PATTERNS["product"], file_tmp.split('/')[-1])[0]
            )

        # Make sure it's the same product
        assert len({x[3:] for x in products}) == 1, "Found different products in input files!"

        # remove duplicates
        if len(set(self.rawdates)) != len(self.files):
            duplicate_warning_msg = f"Possibly duplicate files in {self.targetdir}! Using files with most recent processing date."
            log.warning("Found duplicate files in %s", str(self.targetdir))
            warnings.warn(duplicate_warning_msg, Warning)

            # make sure number of dates is equal to number of files, so no duplicates!
            processing_timestamps = [int(re.sub(REGEX_PATTERNS["processing_timestamp"], '\\1', x)) for x in self.files]

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

            # assert that there are no duplicates now
            assert len(set(self.rawdates)) == len(self.files)

        if vam_product_code is None:

            refname = self.files[0].split('/')[-1]

            if re.match(REGEX_PATTERNS["VIM"], refname):
                self.vam_product_code = "VIM"

            elif re.match(REGEX_PATTERNS["LST"], refname):
                self.vam_product_code = "LTD"
            else:
                pass
        else:
            assert vam_product_code in VAM_PRODUCT_CODES, "Supplied product code not supported."
            self.vam_product_code = vam_product_code

        satset = {x[:3] for x in products}

        if not self.vam_product_code == "VIM":
            log.debug("Interleaving only possible for M{O|Y}D13 (VIM) products! Proceeding without interleave.")
            interleave = False

        if interleave:

            # assert we have both satellites
            assert len(satset) == 2, "Interleaving needs MOD & MYD products!"
            self.satellite = "MXD"
            self.product = f"MXD{products[0][3:]}"

            # for interleaving, exclude all dates before 1st aqua date
            ix = 0
            min_date = fromjulian("2002185")
            for x in self.rawdates:
                start = ix
                if fromjulian(x) >= min_date:
                    break
                ix += 1

            # update instance variables
            self.files = self.files[start:]
            self.rawdates = self.rawdates[start:]

        else:
            # if no interleave, only one satellite allowed
            assert len(satset) == 1, "Without interleave, only one satellite allowed!"
            self.satellite = list(satset)[0]
            self.product = products[0]

        tempinfo = TEMPORAL_DICT[self.product[:5]]

        self.globalproduct = False
        self.temporalresolution = tempinfo["temporalresolution"]
        self.tshift = tempinfo["tshift"]
        self.nfiles = len(self.files)
        self.reference_file = Path(self.files[0])

        ## Build filename

        # LST has specific VAM product codes for the file naming
        if self.vam_product_code in ["LTD", "LTN"]:
            fn_code = LST_NAME_LUD[self.vam_product_code][self.satellite]
            test_vampc = REGEX_PATTERNS["LST"]
        else:
            fn_code = self.vam_product_code
            test_vampc = REGEX_PATTERNS["VIM"]

        assert test_vampc.match(self.product), f"Incomaptible VAM code {self.vam_product_code} for product {self.product}"

        refname = self.reference_file.name

        tile = REGEX_PATTERNS["tile"].findall(refname)
        version = REGEX_PATTERNS["version"].findall(refname)
        tile_version = '.'.join(tile + version)

        if not tile:
            self.globalproduct = True

        filename = f"{self.targetdir}/{fn_code}/{self.product}.{tile_version}.{fn_code}.h5"

        super().__init__(filename=filename)

    def create(self,
               compression: str = 'gzip',
               chunks: Tuple[int] = None) -> None:
        """Creates HDF5 file for raw data.

        Args:
            compression (str): Compression for data (default = gzip).
            chunks (Tuple[int]): Chunksize for data (tuple of 2 int; default = (rows//25, 1)).
        """

        sds_indicator = VAM_PRODUCT_CODES[self.vam_product_code]

        ref_metadata = self._get_reference_metadata(
            reference_file=str(self.reference_file),
            sds_filter=sds_indicator
        )

        row_number = ref_metadata["RasterYSize"]
        col_number = ref_metadata["RasterXSize"]

        # check chunksize
        if chunks is None:
            chunks = ((row_number*col_number)//25, 1) # default
        else:
            assert isinstance(chunks, tuple), "Need chunks as tuple of length = 2"
            assert len(chunks) == 2, "Need chunks as tuple of length = 2"

        # Create directory if necessary
        self.filename.parent.mkdir(parents=True, exist_ok=True)

        log.info("Creating %s", str(self.filename))

        # Create HDF5 file
        try:

            with h5py.File(self.filename, 'x', libver='latest') as h5f:

                # create data array
                dset = h5f.create_dataset(
                    name='data',
                    shape=(row_number*col_number, self.nfiles),
                    dtype="int16",
                    maxshape=(row_number*col_number, None),
                    chunks=chunks,
                    compression=compression,
                    fillvalue=ref_metadata["nodata"])

                # create dates
                _ = h5f.create_dataset(
                    name='dates',
                    shape=(self.nfiles,),
                    maxshape=(None,),
                    dtype='S8',
                    compression=compression
                )

                # set attributes of data array
                # reference metadata
                dset.attrs.update(ref_metadata)

                if self.vam_product_code == "VIM":
                    valid_range = (-2000, 10000)
                elif self.vam_product_code in ["LTD", "LTN"]:
                    valid_range = (7500, 65535)
                else:
                    pass

                # teporal information
                dset.attrs.update({
                    "temporalresolution": self.temporalresolution,
                    "tshift": self.tshift,
                    "globalproduct": self.globalproduct,
                    "vamcode": self.filename.name.split('.')[-2],
                    "valid_range": valid_range,
                })

            self.exists = True
        except Exception as _:
            log.error("Error creating %s", str(self.filename))
            raise HDF5CreationError(f"Error creating {str(self.filename)}! Check if file exists, or if compression / chunksize is OK.")

    def update(self):
        """Update MODIS raw HDF5 file with raw data.

        When a new HDF5 file is created, update will also handle the first data ingest.
        """

        log.info("Updating %s", str(self.filename))

        with h5py.File(self.filename, 'r+', libver='latest') as h5f:
            dset = h5f.get('data')
            dates = h5f.get('dates')
            dataset_shape = dset.shape

            # set processing timestamp
            timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            dset.attrs['processingtimestamp'] = timestamp

            # Load any existing dates and combine with new dates
            dates_stored = [x.decode() for x in dates[...] if x and x.decode() not in self.rawdates]
            dates_combined = dates_stored + self.rawdates
            dates_combined.sort()

            # assert all dates are after last update
            assert dates_combined[:-len(self.rawdates)] == dates_stored, \
                "Files to be collected need to be sequential and AFTER dates of previous updates!"

            # New total temporal length
            dates_length = len(dates_combined)

            # if new total temporal length is bigger than dataset, datasets need to be resized for additional data
            if dates_length > dataset_shape[1]:
                dates.resize((dates_length,))
                dset.resize((dataset_shape[0], dates_length))

            # get required metadata
            attrs = {key: dset.attrs[key] for key in ["nodata", "RasterXSize", "RasterYSize"]}
            chunks = dset.chunks

        # start index for write
        start_index = dates_combined.index(self.rawdates[0])

        # Manual garbage collect to prevent out of memory
        _ = [gc.collect() for x in range(3)]

        # preallocate array
        arr = np.full((chunks[0], self.nfiles), attrs["nodata"], dtype="int16")

        sds_indicator = VAM_PRODUCT_CODES[self.vam_product_code]
        ysize = chunks[0]//attrs["RasterXSize"]

        hdf_datasets = HDFHandler(files=self.files, sds=sds_indicator)

        block_gen = ((x, x//attrs["RasterXSize"]) for x in range(0, dataset_shape[0], chunks[0]))

        log.debug("Opening HDF files ... ")
        with hdf_datasets.open_datasets():
            log.debug("Iterating chunks ... ")
            for yoff_ds, yoff_rst in block_gen:
                log.debug("Processing chunk (%s, %s)", yoff_ds, yoff_rst)

                for ii, hdf in hdf_datasets.iter_handles():
                    try:
                        arr[:, ii] = hdf_datasets.read_chunk(hdf, yoff=int(yoff_rst), ysize=int(ysize)).flatten()
                    except AttributeError:
                        log.warning("Error reading from %s, using nodata.", hdf_datasets.files[ii])
                        arr[:, ii] = attrs["nodata"]

                # write to HDF5
                write_check = self.write_chunk(
                    dataset="data",
                    arr_in=arr,
                    xchunk=10,
                    xoff=start_index,
                    yoff=yoff_ds
                )

                if not write_check:
                    msg = "Error writing to %s"
                    log.error(msg, self.filename)
                    raise HDF5WriteError(msg % self.filename)

        # Write back date list
        with h5py.File(self.filename, 'r+', libver='latest') as h5f:
            dates = h5f.get('dates')
            dates[...] = np.array(dates_combined, dtype='S8')

    @property
    def last_collected(self):
        """Last collected date in file"""
        assert self.exists, "File doesn't exist!"

        with h5py.File(self.filename, 'r') as h5_open:
            dates = h5_open.get("dates")
            last_date = dates[-1].decode()

        return last_date
