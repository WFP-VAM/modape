"""
MODIS smooth HDF5 class.

This file contains the class representing a smoothed MODIS HDF5 file.
"""
# pylint: disable=import-error
from array import array
from datetime import datetime, timedelta
import logging
from pathlib import Path

import h5py
from modape.constants import TEMPINT_LABELS
from modape.exceptions import HDF5CreationError, HDF5WriteError
from modape.modis.io import HDF5Base
from modape.utils import DateHelper, fromjulian
from modape.whittaker import lag1corr, ws2d, ws2dp, ws2doptv, ws2doptvp # pylint: disable=no-name-in-module
import numpy as np

log = logging.getLogger(__name__)

class ModisSmoothH5(HDF5Base):
    """Class representing HDF5 file containing smoothed MODIS data.

    A smooth HDF5 file is directly linked to it's raw HDF5 counterpart. When running
    the Whittaker smoother, data is read from the raw HDF5 file, smoothed and gapfilled
    (and if requested a temporal interpolation is performed), the resulting data is then written
    to the smooth HDF5 file.
    """

    def __init__(self,
                 rawfile: str,
                 targetdir: str,
                 startdate: str = None,
                 tempint: int = None) -> None:
        """Initialize instance of `ModisSmoothH5` class.

        To create an instance of `ModisSmoothH5`, a full path to
        a `ModisRawH5` HDF5 file needs to be specified.
        If the HDF5 file for the `ModisSmoothH5` already exists in the `targetdir`,
        the file will be updated, otherwise created.

        To perform temporal interpolation, a desired timestep
        for a temporal grid has to be specified with `tempint`.
        If `tempint` is 5 or 10, the grid is set to a specific default.
        Otherwise, `startdate` can be used to define a default grid in
        conjunction with `tempint`.

        Args:
            rawfile (str): Full path to raw HDF5 file.
            targetdir (str): Target directory for smooth HDF5 file.
            startdate (str): Start date for temporal interpolation (as julian date YYYYDDD).
            tempint (int): timesteps for temporal interpolation.

        Raises:
            AssertionError: If specified raw HDF5 file does not exist.
        """

        self.rawfile = Path(rawfile)
        self.startdate = startdate
        assert self.rawfile.exists(), f"Raw HDF5 file {self.rawfile} does not exist."

        # Parse tempint to get flag for filename
        if tempint is not None:
            try:
                txflag = TEMPINT_LABELS[int(tempint)]
            except KeyError:
                txflag = "c"

            self.tinterpolate = True
            self.temporalresolution = tempint

        else:
            txflag = "n"
            self.tinterpolate = False
            self.temporalresolution = None

        # Filename for smoothed HDF5
        rawfile_trunk = self.rawfile.name.split(".")
        smoothfile_trunk = ".".join(
            rawfile_trunk[:-2] + \
            ["tx"+txflag] + \
            rawfile_trunk[-2:-1]
        )

        filename = f"{targetdir}/{smoothfile_trunk}.h5"

        super().__init__(filename=filename)

    def create(self):
        """Creates HDF5 file.

        If the corresponding HDF5 is not found in the target directory,
        it's created.

        Raises:
            HDF5CreationError: If creation of HDF5 file fails.
        """

        # Try reading info from raw HDF5
        #pylint: disable=R1721
        with h5py.File(self.rawfile, "r") as h5f_raw:
            dset = h5f_raw.get("data")
            rawshape = dset.shape
            rawchunks = dset.chunks
            datatype = dset.dtype
            compression = dset.compression
            raw_dates_all = [x.decode() for x in h5f_raw.get("dates")[...]]
            raw_attrs = {key:value for key, value in dset.attrs.items()}

        if self.temporalresolution is None:
            tempres = raw_attrs["temporalresolution"]
        else:
            tempres = self.temporalresolution

        dates = DateHelper(rawdates=raw_dates_all,
                           rtres=int(raw_attrs["temporalresolution"]),
                           stres=int(tempres),
                           start=self.startdate)

        dates_length = len(dates.target)

        nrows = raw_attrs["RasterYSize"]
        ncols = raw_attrs["RasterXSize"]

        try:
            with h5py.File(self.filename, "x", libver="latest") as h5f:

                dset = h5f.create_dataset("data",
                                          shape=(rawshape[0], dates_length),
                                          dtype=datatype, maxshape=(rawshape[0], None),
                                          chunks=rawchunks,
                                          compression=compression,
                                          fillvalue=raw_attrs["nodata"])

                h5f.create_dataset("sgrid",
                                   shape=(nrows*ncols,),
                                   dtype="float32",
                                   maxshape=(nrows*ncols,),
                                   chunks=(rawchunks[0],),
                                   compression=compression)

                h5f.create_dataset("dates",
                                   shape=(dates_length,),
                                   maxshape=(None,),
                                   dtype="S8",
                                   compression=compression,
                                   data=np.array(dates.target, dtype="S8"))

                h5f.create_dataset("rawdates",
                                   shape=(len(raw_dates_all),),
                                   maxshape=(None,),
                                   dtype="S8",
                                   compression=compression,
                                   data=np.array(raw_dates_all, dtype="S8"))

                raw_attrs["temporalresolution"] = tempres
                dset.attrs.update(raw_attrs)

                self.exists = True
        except Exception as _:
            log.error("Error creating %s", str(self.filename))
            raise HDF5CreationError(f"Error creating {str(self.filename)}!")

    def smooth(self,
               svalue: float = None,
               p: float = None,
               soptimize: bool = None,
               srange: np.ndarray = None,
               nsmooth: int = 0,
               nupdate: int = 0,
               ) -> None:
        """Applies Whittaker smoother to the data.

        This method reads raw data from the raw HDF5 and applies the Whittaker
        filter according to the parameters supplied. The resulting filtered data
        is then stored in the smooth HDF5 file.

        The parameters relevant for derermining which Whittaker
        variant is supplied, are `svalue`, `soptimize` and `p`:

        If the `soptimize` flag is `True`, the V-curve optimization is performed.
        If in addition a `p` value is supplied, the assymmetric V-curve optimization is performed instead.
        In both cases, the range of S values to optimize in can be supplied with `srange`.

        If `soptimize` is `False` and a `svalue` is passed, smoothing using that fixed value (needs to be `log10` of S)
        for each pixel is performed. If a `p` value is passed, asymmetric smoothing with the fixed value is performed instead.

        If none of `svalue` and `soptimize` is specified, a previosuly initialized S value for each pixel
        will be read from grid and used for the smoothing.

        Using `nsmooth` and `nupdate`, the number or raw timesteps used for smoothing, and the number
        of filtered timesteps updated in the smooth HDF5 target file, can be adjusted.

        Args:
            svalue (float): Log10 value of smoothing parameter S (for fixed smoothing).
            p (float): P value for asymmetric smoothing.
            soptimize (bool): Flag for V-curve optimization.
            srange (np.ndarray): S-range for V-curve optimization.
            nsmooth (int): Number of raw timesteps for smoothing.
            nupdate (int): Number of smooth timesteps updated in file.

        Raises:
            AssertionError: If smooth HDF5 file does not exist (`create` needs to be executed before)
            ValueError: If `nsmooth` is smaller than `nupdate`.
            ValueError: If `srange` is not specified as a numpy array with expected dimensions.
            HDF5WriteError: If write to HDF5 file fails.
        """

        assert self.filename.exists(), "File doesn't exist! Can't run smoother."

        if (nsmooth != 0) and (nupdate != 0):
            if nsmooth < nupdate:
                raise ValueError("nsmooth must be bigger or equal (>=) to nupdate!")

        if soptimize and srange is not None:
            if not isinstance(srange, np.ndarray):
                raise ValueError("srange needs to be supplied as numpy array")

        log.info("Runnig smoother on %s", str(self.filename))

        processing_starttime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        log.debug("Reading metadata from HDF5s")
        with h5py.File(self.rawfile, "r") as h5f_open:
            raw_ds = h5f_open.get("data")
            raw_shape = raw_ds.shape
            raw_chunks = raw_ds.chunks
            raw_attrs = dict(raw_ds.attrs.items())
            raw_dates_all = [x.decode() for x in h5f_open.get("dates")[...]]
            raw_dates_nsmooth = raw_dates_all[-nsmooth:]

        with h5py.File(self.filename, "r+") as h5f_open:
            smt_ds = h5f_open.get("data")
            smt_attrs = smt_ds.attrs
            smt_shape = smt_ds.shape
            smt_chunks = smt_ds.chunks
            temporalresolution = smt_attrs["temporalresolution"]
            tshift = int(smt_attrs["tshift"])

            dates = DateHelper(rawdates=raw_dates_all,
                               rtres=int(raw_attrs["temporalresolution"]),
                               stres=int(temporalresolution),
                               start=self.startdate)

            # Resize if date list is bigger than shape of smoothed data
            if dates.target_length > smt_shape[1]:
                log.debug("Resizing dataset! Current %s, required %s", smt_shape[1], dates.target_length)
                smt_dates = h5f_open.get("dates")
                smt_dates.resize((dates.target_length,))
                smt_ds.resize((smt_shape[0], dates.target_length))
                smt_dates[...] = np.array(dates.target, dtype="S8")
                smt_shape = smt_ds.shape

        nodata = raw_attrs["nodata"]
        dix = dates.getDIX(nupdate)

        arr_raw = np.zeros((raw_chunks[0], len(raw_dates_nsmooth)), dtype="double")

        if nsmooth > 0:
            read_offset = raw_shape[1] - nsmooth
        else:
            read_offset = 0

        if nupdate > 0:
            write_offset = smt_shape[1] - nupdate
        else:
            write_offset = 0

        if self.tinterpolate:
            log.debug("Temporal interpolation triggered!")
            arr_smt = np.full((smt_chunks[0], len(dix)), fill_value=nodata, dtype="double")
            vector_daily = dates.getDV(nodata)

            # Shift for interpolation
            for rdate in raw_dates_nsmooth:
                rdate_shift = (fromjulian(rdate) + timedelta(tshift)).strftime("%Y%j")
                dd_index = dates.daily.index(rdate_shift)
                vector_daily[dd_index] = -1
        else:
            arr_smt = arr_raw

        # use HDF5 base for reading rawdata
        raw_h5 = HDF5Base(str(self.rawfile))

        chunk_generator = raw_h5.read_chunked(
            dataset="data",
            xoffset=read_offset,
            xchunk=10,
            arr_out=arr_raw,
            )

        if soptimize or svalue is None:
            sgrid_generator = self.read_chunked(dataset="sgrid")
        else:
            sgrid_generator = None
            arr_sgrid = None

        # counter
        chunk_counter = 0

        # iterate over chunks
        log.debug("Iterating over chunks")
        for arr_raw_chunk in chunk_generator:
            log.debug("Chunk %s", chunk_counter)

            # create weights
            wts = (arr_raw_chunk != nodata)
            wts = wts.astype("uint8")

            ndix = np.sum(wts, 1) >= (arr_raw.shape[1] * 0.2) # 20%+ data
            map_index = np.where(ndix)[0]

            if sgrid_generator:
                arr_sgrid = next(sgrid_generator)

            for ix in map_index:

                w = wts[ix, :].astype("double")
                if soptimize:
                    if srange is None:
                        lag_correlation = lag1corr(arr_raw[ix, :-1], arr_raw[ix, 1:], nodata)
                        if lag_correlation > 0.5:
                            sr = np.arange(-2, 1.2, 0.2).round(2)
                        elif lag_correlation <= 0.5:
                            sr = np.arange(0, 3.2, 0.2).round(2)
                        else:
                            sr = np.arange(-1, 1.2, 0.2).round(2)

                    else:
                        sr = srange

                    if p is not None:
                        arr_raw[ix, :], arr_sgrid[ix] = ws2doptvp(y=arr_raw[ix, :],
                                                                  w=w,
                                                                  llas=array("d", sr),
                                                                  p=p)
                    else:
                        arr_raw[ix, :], arr_sgrid[ix] = ws2doptv(y=arr_raw[ix, :],
                                                                 w=w,
                                                                 llas=array("d", sr))

                else:
                    if svalue is None:
                        s = 10 ** arr_sgrid[ix]
                    else:
                        s = 10 ** svalue

                    if p is None:
                        arr_raw[ix, :] = ws2d(y=arr_raw[ix, :], lmda=s, w=w)
                    else:
                        arr_raw[ix, :] = ws2dp(y=arr_raw[ix, :], lmda=s, w=w, p=p)


                if self.tinterpolate:

                    arr_smt[ix, :] = self._apply_tinterpolate(
                        z1=arr_raw[ix, :],
                        nodata=nodata,
                        vector_daily=vector_daily,
                        dix=dix,
                    )

            arr_smt = np.rint(arr_smt, out=arr_smt)

            write_check = self.write_chunk(
                dataset="data",
                arr_in=arr_smt,
                xoffset=write_offset,
                xchunk=10,
                yoffset=chunk_counter*raw_chunks[0],
            )

            if not write_check:
                msg = "Error writing to %s"
                log.error(msg, self.filename)
                raise HDF5WriteError(msg % self.filename)

            if arr_sgrid is not None:

                arr_sgrid[arr_sgrid > 0] = np.log10(arr_sgrid[arr_sgrid > 0])

                write_check = self.write_chunk(
                    dataset="sgrid",
                    arr_in=arr_sgrid,
                    yoffset=chunk_counter*raw_chunks[0],
                )

                if not write_check:
                    msg = "Error writing sgrid to %s"
                    log.error(msg, self.filename)
                    raise HDF5WriteError(msg % self.filename)

            # increment counter
            chunk_counter += 1

            if self.tinterpolate:
                # flush smooth array
                arr_smt[...] = nodata

        # processing information for modis_info
        if soptimize:
            processing_info = {"lastrun": "V-curve optimization of s"}

        else:
            if svalue is None:
                processing_info = {"lastrun": "fixed s from grid"}
            else:
                processing_info = {"lastrun": f"fixed s {svalue} (log10)"}

        if p is not None:
            processing_info.update({"pvalue": p})
            processing_info["lastrun"] = processing_info["lastrun"] + f" and with P-value of {p}"

        log.debug("Last run: %s", processing_info["lastrun"])

        with h5py.File(self.filename, "r+") as h5f_open:
            smt_ds = h5f_open.get("data")

            if p is None:
                try:
                    del smt_ds.attrs["pvalue"]
                except KeyError:
                    pass

            processing_info.update({
                "processingtimestamp": processing_starttime
                })

            smt_ds.attrs.update(processing_info)

            dates_ds = h5f_open.get("rawdates")

            if len(raw_dates_all) > dates_ds.size:
                dates_ds.resize((len(raw_dates_all),))

            dates_ds[...] = np.array(raw_dates_all, dtype="S8")


    @property
    def last_collected(self):
        """Last collected date in raw file"""
        with h5py.File(self.rawfile, "r") as h5f_raw:
            dates = h5f_raw.get("dates")
            last_date = dates[-1].decode()

        return last_date

    @property
    def last_smoothed(self):
        """Last smoothed date in file"""
        assert self.exists, "File doesn't exist!"

        with h5py.File(self.filename, "r") as h5_open:
            dates = h5_open.get("rawdates")
            last_date = dates[-1].decode()

        return last_date

    #pylint: disable=C0103
    @staticmethod
    def _apply_tinterpolate(z1, nodata, vector_daily, dix):
        z2 = vector_daily.copy()
        z2[z2 != nodata] = z1
        z2[...] = ws2d(
            y=z2,
            lmda=0.0001,
            w=np.array((z2 != nodata) * 1, dtype="double")
        )

        return z2[dix]
