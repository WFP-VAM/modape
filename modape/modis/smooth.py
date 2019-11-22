"""
MODIS smooth HDF5 class.

This file contains the class representing a smoothed MODIS HDF5 file.

Author: Valentin Pesendorfer, April 2019
"""
from __future__ import absolute_import, division, print_function

from array import array
import ctypes
from datetime import timedelta
import multiprocessing as mp
import os
from os.path import basename
try:
    from pathlib2 import Path
except ImportError:
    from pathlib import Path
import time
import warnings

import numpy as np
import h5py # pylint: disable=import-error
import xarray as xr # pylint: disable=import-error

from modape.utils import (DateHelper, execute_ws2d, execute_ws2d_sgrid, execute_ws2d_vc, execute_w_constraint, execute_tempint,
                          fromjulian, init_parameters, init_shared, init_worker, tonumpyarray, txx)
from modape.whittaker import lag1corr, w_constrain, ws2d, ws2dp, ws2doptv, ws2doptvp # pylint: disable=no-name-in-module

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

        self.targetdir = Path(targetdir)
        self.rawfile = rawfile
        self.nworkers = nworkers
        self.nupdate = nupdate
        self.nsmooth = nsmooth
        self.array_offset = nsmooth - nupdate
        self.startdate = startdate

        # Get info from raw HDF5
        with h5py.File(self.rawfile, 'r') as h5f:
            dates = h5f.get('dates')
            self.rawdates_nsmooth = [x.decode() for x in dates[-self.nsmooth:]]

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
        self.outname = Path('{}/{}.tx{}.{}.h5'.format(
            self.targetdir.as_posix(),
            '.'.join(basename(rawfile).split('.')[:-2]),
            txflag,
            basename(rawfile).split('.')[-2:-1][0]))

        self.exists = self.outname.exists()

    def create(self):
        """Creates smoothed HDF5 file on disk."""

        # Try reading info from raw HDF5
        try:
            with h5py.File(self.rawfile, 'r') as h5f:
                dset = h5f.get('data')
                raw_dates_all = [x.decode() for x in h5f.get('dates')[...]]
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

        dates = DateHelper(rawdates=raw_dates_all,
                           rtres=rtres,
                           stres=self.temporalresolution,
                           start=self.startdate)

        dates_length = len(dates.target)

        # create datasets

        try:
            with h5py.File(self.outname.as_posix(), 'x', libver='latest') as h5f:

                # data array
                dset = h5f.create_dataset('data',
                                          shape=(rawshape[0], dates_length),
                                          dtype=dt, maxshape=(rawshape[0], None),
                                          chunks=rawchunks,
                                          compression=cmpr,
                                          fillvalue=rnd)

                # grid with smoothing parameter
                h5f.create_dataset('sgrid',
                                   shape=(nrows*ncols,),
                                   dtype='float32',
                                   maxshape=(nrows*ncols,),
                                   chunks=(rawchunks[0],),
                                   compression=cmpr)

                # lower constraint
                h5f.create_dataset('clower',
                                   shape=(rawshape[0], dates.nmdays),
                                   dtype='Int16', maxshape=(rawshape[0], None),
                                   chunks=(rawchunks[0], 10),
                                   compression=cmpr,
                                   fillvalue=rnd)

                # upper constraint
                h5f.create_dataset('cupper',
                                   shape=(rawshape[0], dates.nmdays),
                                   dtype='Int16', maxshape=(rawshape[0], None),
                                   chunks=(rawchunks[0], 10),
                                   compression=cmpr,
                                   fillvalue=rnd)

                # dates
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
                dset.attrs['mdays_set'] = dates.mdays_set
                dset.attrs['RasterYSize'] = nrows
                dset.attrs['RasterXSize'] = ncols
        except:
            print('\n\nError creating {}! Check input parameters (especially if compression/filter is available) and try again. Corrupt file will be removed now. \n\nError message: \n'.format(self.outname.as_posix()))
            os.remove(self.outname.unlink())
            raise

        self.exists = True

    def ws2d(self, s):
        """Apply whittaker smoother with fixed s-value to data.

        Args:
            s: log10 value of s (float)
        """

        with h5py.File(self.rawfile, 'r') as rawh5, h5py.File(self.outname.as_posix(), 'r+') as smth5:
            raw_ds = rawh5.get('data')
            raw_dates_all = [x.decode() for x in rawh5.get('dates')[...]]
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
            try:
                del smt_ds.attrs['pvalue']
            except KeyError:
                pass

            dates = DateHelper(rawdates=raw_dates_all,
                               rtres=rtres,
                               stres=self.temporalresolution,
                               start=self.startdate)

            dix = dates.getDIX()[-self.nupdate:]

            # Resize if date list is bigger than shape of smoothed data
            if dates.target_length > smoothshape[1]:
                smt_dates.resize((dates.target_length,))
                smt_ds.resize((smoothshape[0], dates.target_length))
                smt_dates[...] = np.array(dates.target, dtype='S8')
                smoothshape = smt_ds.shape

            # calculate offsets
            rawoffset = raw_dates_all.index(self.rawdates_nsmooth[0])

            # if dataset is smaller or equal then nupdate, take index 0
            try:
                smoothoffset = [x.decode() for x in smt_dates[...]].index(dates.target[-self.nupdate])
            except IndexError:
                smoothoffset = [x.decode() for x in smt_dates[...]].index(dates.target[0])

            new_dim = smoothshape[1] - smoothoffset

            if self.nworkers > 1:
                if self.tinterpolate:
                    shared_array_smooth = init_shared(smoothchunks[0] * new_dim, dtype=ctypes.c_double)
                    arr_smooth = tonumpyarray(shared_array_smooth, dtype=ctypes.c_double)
                    arr_smooth.shape = (smoothchunks[0], new_dim)
                    arr_smooth[...] = nodata
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates_nsmooth:
                        vector_daily[dates.daily.index((fromjulian(rdate) + timedelta(tshift)).strftime('%Y%j'))] = -1
                else:
                    vector_daily = None
                    shared_array_smooth = None
                    arr_smooth = None
                shared_array_raw = init_shared(rawchunks[0] * len(self.rawdates_nsmooth), dtype=ctypes.c_double)

                parameters = init_parameters(rdim=(rawchunks[0], len(self.rawdates_nsmooth)),
                                             sdim=(smoothchunks[0], new_dim),
                                             nd=nodata,
                                             s=s,
                                             shared_array_smooth=shared_array_smooth,
                                             vec_dly=vector_daily,
                                             dix=dix)

                arr_raw = tonumpyarray(shared_array_raw, dtype=ctypes.c_double)
                arr_raw.shape = (rawchunks[0], len(self.rawdates_nsmooth))

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

                    if self.tinterpolate:
                        _ = pool.map(execute_tempint, map_index)

                    # write back data
                    if self.tinterpolate:
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(0, arr_smooth.shape[1], smoothchunks[1])):
                            smt_ds[br:br+smoothchunks[0], bcs:bcs+smoothchunks[1]] = arr_smooth[:, bcr:bcr+smoothchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(self.array_offset, arr_raw.shape[1], smoothchunks[1])):
                            smt_ds[br:br+smoothchunks[0], bcs:bcs+smoothchunks[1]] = arr_raw[:, bcr:bcr+smoothchunks[1]]

                # close pool
                pool.close()
                pool.join()

            else:
                arr_raw = np.zeros((rawchunks[0], len(self.rawdates_nsmooth)), dtype='double')

                # Create weights array
                wts = arr_raw.copy()

                if self.tinterpolate:
                    arr_smooth = np.zeros((smoothchunks[0], new_dim), dtype='double')
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates_nsmooth:
                        vector_daily[dates.daily.index((fromjulian(rdate) + timedelta(tshift)).strftime('%Y%j'))] = -1
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
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(0, arr_smooth.shape[1], smoothchunks[1])):
                            smt_ds[br:br+smoothchunks[0], bcs:bcs+smoothchunks[1]] = arr_smooth[:, bcr:bcr+smoothchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(self.array_offset, arr_raw.shape[1], smoothchunks[1])):
                            smt_ds[br:br+rawchunks[0], bcs:bcs+smoothchunks[1]] = arr_raw[:, bcr:bcr+smoothchunks[1]]

    def ws2d_sgrid(self, p=None, constrain=False):
        """Apply whittaker smootehr with fixed s to data.

        This fixed s version reads a pixel based s value from file, so it needs
        a previous run of V-curve s-optimization.
        """

        with h5py.File(self.rawfile, 'r') as rawh5, h5py.File(self.outname.as_posix(), 'r+') as smth5:
            raw_ds = rawh5.get('data')
            raw_dates_all = [x.decode() for x in rawh5.get('dates')[...]]
            rtres = raw_ds.attrs['temporalresolution'].item()
            rtres = raw_ds.attrs['temporalresolution'].item()
            smt_ds = smth5.get('data')
            smt_dates = smth5.get('dates')
            smt_sgrid = smth5.get('sgrid')
            smt_clower = smth5.get('clower')
            smt_cupper = smth5.get('cupper')
            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks
            smoothshape = smt_ds.shape
            smoothchunks = smt_ds.chunks
            nodata = raw_ds.attrs['nodata'].item()
            self.temporalresolution = smt_ds.attrs['temporalresolution'].item()
            tshift = raw_ds.attrs['tshift'].item()

            # if constrain datasets are not available, don't perform constrain
            if not smt_clower and smt_cupper:
                constrain = False

            if constrain:

                # constrain weights
                cweights = array('f', np.arange(0.9, -0.1, -0.1))

                if len(self.rawdates_nsmooth) < len(cweights):
                    raise ValueError('NSMOOTH can\'t be smaller than 10 when CONSTRAIN should be performed!')

            else:
                cweights = None


            dates = DateHelper(rawdates=raw_dates_all,
                               rtres=rtres,
                               stres=self.temporalresolution,
                               start=self.startdate)

            dix = dates.getDIX()[-self.nupdate:]

            if constrain:

                try:
                    mdays_set = smt_ds.attrs['mdays_set']
                    mdays = [int(fromjulian(x).strftime('%m%d')) for x in dates.target[-len(cweights):]]
                    constraint_ix = [list(mdays_set).index(x) for x in mdays]
                    constraint_ix_sorted = constraint_ix.copy()
                    constraint_ix_sorted.sort()

                except KeyError:
                    constrain = False
                    warnings.warn('Error reading constrain information from dataset. Continuing without constrain!')

                except ValueError:
                    raise ValueError('Constraints not for all dates in timeseries available! Please run smoother without --constrain\
                                     or re-initialize smooth HDF5 file with --optvp!')

            # Resize if date list is bigger than shape of smoothed data
            if len(dates.target) > smoothshape[1]:
                smt_dates.resize((dates.target_length,))
                smt_ds.resize((smoothshape[0], dates.target_length))
                smt_dates[...] = np.array(dates.target, dtype='S8')
                smoothshape = smt_ds.shape

            # calculate offsets
            rawoffset = raw_dates_all.index(self.rawdates_nsmooth[0])

            # if dataset is smaller or equal then nupdate, take index 0
            try:
                smoothoffset = [x.decode() for x in smt_dates[...]].index(dates.target[-self.nupdate])
            except IndexError:
                smoothoffset = [x.decode() for x in smt_dates[...]].index(dates.target[0])

            new_dim = smoothshape[1] - smoothoffset

            # Store run parameters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

            if p:
                smt_ds.attrs['lastrun'] = 'fixed s from grid with p = {} - constraining: {}'.format(p, constrain)
                smt_ds.attrs['pvalue'] = p
            else:
                smt_ds.attrs['lastrun'] = 'fixed s from grid - constraining: {}'.format(constrain)
                try:
                    del smt_ds.attrs['pvalue']
                except KeyError:
                    pass

            if self.nworkers > 1:
                if self.tinterpolate:
                    shared_array_smooth = init_shared(smoothchunks[0] * new_dim, dtype=ctypes.c_double)
                    arr_smooth = tonumpyarray(shared_array_smooth, dtype=ctypes.c_double)
                    arr_smooth.shape = (smoothchunks[0], new_dim)
                    arr_smooth[...] = nodata
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates_nsmooth:
                        vector_daily[dates.daily.index((fromjulian(rdate) + timedelta(tshift)).strftime('%Y%j'))] = -1
                else:
                    vector_daily = None
                    shared_array_smooth = None
                    arr_smooth = None
                shared_array_raw = init_shared(rawchunks[0] * len(self.rawdates_nsmooth), dtype=ctypes.c_double)

                parameters = init_parameters(rdim=(rawchunks[0], len(self.rawdates_nsmooth)),
                                             sdim=(smoothchunks[0], new_dim),
                                             nd=nodata,
                                             shared_array_smooth=shared_array_smooth,
                                             vec_dly=vector_daily,
                                             dix=dix,
                                             p=p)

                parameters['shared_array_sgrid'] = init_shared(rawchunks[0], dtype=ctypes.c_double)
                arr_raw = tonumpyarray(shared_array_raw, dtype=ctypes.c_double)
                arr_raw.shape = (rawchunks[0], len(self.rawdates_nsmooth))
                arr_sgrid = tonumpyarray(parameters['shared_array_sgrid'], dtype=ctypes.c_double)

                if constrain:
                    parameters['shared_array_clower'] = init_shared(rawchunks[0]*len(cweights), dtype=ctypes.c_int16)
                    arr_clower = tonumpyarray(parameters['shared_array_clower'], dtype=ctypes.c_int16)
                    arr_clower.shape = (rawchunks[0], len(cweights))

                    parameters['shared_array_cupper'] = init_shared(rawchunks[0]*len(cweights), dtype=ctypes.c_int16)
                    arr_cupper = tonumpyarray(parameters['shared_array_cupper'], dtype=ctypes.c_int16)
                    arr_cupper.shape = (rawchunks[0], len(cweights))

                    shared_cweights = init_shared(len(cweights), dtype=ctypes.c_float)
                    shared_cweights[0:len(cweights)] = cweights[0:len(cweights)]

                    parameters['cweights'] = shared_cweights

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

                    if self.tinterpolate:
                        _ = pool.map(execute_tempint, map_index)


                    if constrain:
                        arr_clower[...] = smt_clower[br:br+rawchunks[0], constraint_ix_sorted]
                        arr_clower[...] = arr_clower[:, np.argsort(constraint_ix)]

                        arr_cupper[...] = smt_cupper[br:br+rawchunks[0], constraint_ix_sorted]
                        arr_cupper[...] = arr_cupper[:, np.argsort(constraint_ix)]

                        _ = pool.map(execute_w_constraint, map_index)

                    # write back data
                    if self.tinterpolate:
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(0, arr_smooth.shape[1], smoothchunks[1])):
                            smt_ds[br:br+smoothchunks[0], bcs:bcs+smoothchunks[1]] = arr_smooth[:, bcr:bcr+smoothchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(self.array_offset, arr_raw.shape[1], smoothchunks[1])):
                            smt_ds[br:br+smoothchunks[0], bcs:bcs+smoothchunks[1]] = arr_raw[:, bcr:bcr+smoothchunks[1]]

                # close pool
                pool.close()
                pool.join()

            else:
                arr_raw = np.zeros((rawchunks[0], len(self.rawdates_nsmooth)), dtype='double')
                arr_sgrid = np.zeros((rawchunks[0],), dtype='double')

                if constrain:
                    arr_clower = np.zeros((rawchunks[0], len(cweights)), dtype='int16')
                    arr_cupper = arr_clower.copy()

                # Create weights array
                wts = arr_raw.copy()
                if self.tinterpolate:
                    arr_smooth = np.zeros((smoothchunks[0], new_dim), dtype='double')
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates_nsmooth:
                        vector_daily[dates.daily.index((fromjulian(rdate) + timedelta(tshift)).strftime('%Y%j'))] = -1
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

                    if constrain:
                        arr_clower[...] = smt_clower[br:br+rawchunks[0], constraint_ix]
                        arr_cupper[...] = smt_cupper[br:br+rawchunks[0], constraint_ix]

                    for ix in map_index:
                        if not p:
                            arr_raw[ix, :] = ws2d(y=arr_raw[ix, :], lmda=10**arr_sgrid[ix], w=wts[ix, :])
                        else:
                            arr_raw[ix, :] = ws2dp(y=arr_raw[ix, :], lmda=10**arr_sgrid[ix], w=wts[ix, :], p=p)
                        if self.tinterpolate:
                            z2 = vector_daily.copy()
                            z2[z2 != nodata] = arr_raw[ix, :]
                            z2[...] = ws2d(y=z2, lmda=0.0001, w=np.array((z2 != nodata)*1, dtype='double'))
                            arr_smooth[ix, :] = z2[dix]
                            if constrain:
                                arr_smooth[ix, :] = w_constrain(arr_smooth[ix, :], arr_clower[ix, :], arr_cupper[ix, :], cweights)
                        else:
                            if constrain:
                                arr_raw[ix, :] = w_constrain(arr_raw[ix, :], arr_clower[ix, :], arr_cupper[ix, :], cweights)
                            else:
                                pass

                    # write back data
                    if self.tinterpolate:
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(0, arr_smooth.shape[1], smoothchunks[1])):
                            smt_ds[br:br+smoothchunks[0], bcs:bcs+smoothchunks[1]] = arr_smooth[:, bcr:bcr+smoothchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(self.array_offset, arr_raw.shape[1], smoothchunks[1])):
                            smt_ds[br:br+rawchunks[0], bcs:bcs+smoothchunks[1]] = arr_raw[:, bcr:bcr+smoothchunks[1]]

    def ws2d_vc(self, srange, p=None):
        """Apply whittaker smoother V-curve optimization of s.

        Optionally, p value can be specified to use asymmetric smoothing.

        Args:
            srange: array of s-values to apply
            p: Percentile value (float)
        """

        with h5py.File(self.rawfile, 'r') as rawh5, h5py.File(self.outname.as_posix(), 'r+') as smth5:
            raw_ds = rawh5.get('data')
            raw_dates_all = [x.decode() for x in rawh5.get('dates')[...]]
            rtres = raw_ds.attrs['temporalresolution'].item()
            rtres = raw_ds.attrs['temporalresolution'].item()
            smt_ds = smth5.get('data')
            smt_dates = smth5.get('dates')
            smt_sgrid = smth5.get('sgrid')
            smt_cupper = smth5.get('cupper')
            smt_clower = smth5.get('clower')
            rawshape = raw_ds.shape
            rawchunks = raw_ds.chunks
            smoothshape = smt_ds.shape
            smoothchunks = smt_ds.chunks
            nodata = raw_ds.attrs['nodata'].item()
            self.temporalresolution = smt_ds.attrs['temporalresolution'].item()
            tshift = raw_ds.attrs['tshift'].item()
            ncols = smt_ds.attrs['RasterXSize'].item()

            # Store run vampceters for infotool
            smt_ds.attrs['processingtimestamp'] = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())

            if p:
                smt_ds.attrs['lastrun'] = 'V-curve optimization of s with p = {}'.format(p)
                smt_ds.attrs['pvalue'] = p
            else:
                smt_ds.attrs['lastrun'] = 'V-curve optimization of s'
                try:
                    del smt_ds.attrs['pvalue']
                except KeyError:
                    pass

            dates = DateHelper(rawdates=raw_dates_all,
                               rtres=rtres,
                               stres=self.temporalresolution,
                               start=self.startdate)

            dix = dates.getDIX()[-self.nupdate:]

            # Resize if date list is bigger than shape of smoothed data
            if dates.target_length > smoothshape[1]:
                smt_dates.resize((dates.target_length,))
                smt_ds.resize((smoothshape[0], dates.target_length))
                smt_dates[...] = np.array(dates.target, dtype='S8')
                smoothshape = smt_ds.shape

            # calculate offsets
            rawoffset = raw_dates_all.index(self.rawdates_nsmooth[0])

            # if dataset is smaller or equal then nupdate, take index 0
            try:
                smoothoffset = [x.decode() for x in smt_dates[...]].index(dates.target[-self.nupdate])
            except IndexError:
                smoothoffset = [x.decode() for x in smt_dates[...]].index(dates.target[0])

            new_dim = smoothshape[1] - smoothoffset

            if self.nworkers > 1:
                if self.tinterpolate:
                    shared_array_smooth = init_shared(smoothchunks[0] * new_dim, dtype=ctypes.c_double)
                    arr_smooth = tonumpyarray(shared_array_smooth, dtype=ctypes.c_double)
                    arr_smooth.shape = (smoothchunks[0], new_dim)
                    arr_smooth[...] = nodata
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates_nsmooth:
                        vector_daily[dates.daily.index((fromjulian(rdate) + timedelta(tshift)).strftime('%Y%j'))] = -1
                else:
                    vector_daily = None
                    shared_array_smooth = None
                    arr_smooth = None
                shared_array_raw = init_shared(rawchunks[0] * len(self.rawdates_nsmooth), dtype=ctypes.c_double)

                parameters = init_parameters(rdim=(rawchunks[0], len(self.rawdates_nsmooth)),
                                             sdim=(smoothchunks[0], new_dim),
                                             nd=nodata,
                                             p=p,
                                             shared_array_smooth=shared_array_smooth,
                                             vec_dly=vector_daily,
                                             dix=dix,
                                             srange=srange)

                parameters['shared_array_sgrid'] = init_shared(rawchunks[0], dtype=ctypes.c_double)
                arr_raw = tonumpyarray(shared_array_raw, dtype=ctypes.c_double)
                arr_raw.shape = (rawchunks[0], len(self.rawdates_nsmooth))
                arr_sgrid = tonumpyarray(parameters['shared_array_sgrid'], dtype=ctypes.c_double)

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

                    if self.tinterpolate:
                        _ = pool.map(execute_tempint, map_index)

                    # create constraints (only for full timeseries and if datasets are available)
                    if self.nsmooth == 0:

                        if smt_clower and smt_cupper:

                            if dates.nmdays > smt_cupper.shape[1]:
                                smt_cupper.resize((smt_cupper.shape[0], dates.nmdays))

                            if dates.nmdays > smt_clower.shape[1]:
                                smt_clower.resize((smt_clower.shape[0], dates.nmdays))

                            if self.tinterpolate:
                                smooth_view = arr_smooth.view()
                            else:
                                smooth_view = arr_raw.view()

                            smooth_view.shape = (-1, ncols, dates.target_length)


                            col_ix = np.arange(smooth_view.shape[1])
                            row_ix = np.arange(smooth_view.shape[0])

                            smooth_xarr = xr.Dataset({'data':(['x', 'y', 'time'], smooth_view)},
                                                     coords={'cols': (['x', 'y'], np.repeat([col_ix], smooth_view.shape[0], axis=0)),
                                                             'rows': (['x', 'y'], np.repeat([row_ix], smooth_view.shape[1], axis=1).reshape(smooth_view.shape[0], -1)),
                                                             'time': [int(fromjulian(x).strftime('%m%d')) for x in  dates.target]})

                            # no copy
                            if self.tinterpolate:
                                assert smooth_view.base is arr_smooth
                                assert smooth_xarr.data.values.base is arr_smooth
                            else:
                                assert smooth_view.base is arr_raw
                                assert smooth_xarr.data.values.base is arr_raw

                            smooth_diff_mean = smooth_xarr.diff('time').groupby('time').mean()
                            smooth_diff_sdev = smooth_xarr.diff('time').groupby('time').std()
                            # enforce min sdev of one NDVI value
                            smooth_diff_sdev = smooth_diff_sdev.where(smooth_diff_sdev.data > 100, 100)

                            # check dimensions are correct
                            assert smooth_diff_mean.data.shape[0] == smooth_diff_sdev.data.shape[0] == dates.nmdays

                            constraint = (smooth_diff_mean + smooth_diff_sdev).data.values

                            # flip for writing
                            constraint_view = np.moveaxis(constraint, 0, 2).reshape(constraint.shape[1]*constraint.shape[2], constraint.shape[0])

                            smt_cupper[br:br+smoothchunks[0], :] = constraint_view[...]

                            del constraint, constraint_view

                            constraint = (smooth_diff_mean - smooth_diff_sdev).data.values

                            # flip for writing
                            constraint_view = np.moveaxis(constraint, 0, 2).reshape(constraint.shape[1]*constraint.shape[2], constraint.shape[0])

                            smt_clower[br:br+smoothchunks[0], :] = constraint_view[...]

                            del constraint, constraint_view, smooth_xarr

                        else:
                            warnings.warn('Error reading constrain information from dataset. Won\'t create/update constrains!')

                    else:
                        # constraints won't be created/updated if
                        pass

                    # write back data
                    arr_sgrid[arr_sgrid > 0] = np.log10(arr_sgrid[arr_sgrid > 0])
                    smt_sgrid[br:br+rawchunks[0]] = arr_sgrid[...]
                    arr_sgrid[...] = 0

                    if self.tinterpolate:
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(0, arr_smooth.shape[1], smoothchunks[1])):
                            smt_ds[br:br+smoothchunks[0], bcs:bcs+smoothchunks[1]] = arr_smooth[:, bcr:bcr+smoothchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(self.array_offset, arr_raw.shape[1], smoothchunks[1])):
                            smt_ds[br:br+smoothchunks[0], bcs:bcs+smoothchunks[1]] = arr_raw[:, bcr:bcr+smoothchunks[1]]

                # close pool
                pool.close()
                pool.join()

            else:
                arr_raw = np.zeros((rawchunks[0], len(self.rawdates_nsmooth)), dtype='double')
                arr_sgrid = np.zeros((rawchunks[0],), dtype='double')
                wts = arr_raw.copy() # Create weights array

                if self.tinterpolate:
                    arr_smooth = np.zeros((smoothchunks[0], new_dim), dtype='double')
                    vector_daily = dates.getDV(nodata)

                    # Shift for interpolation
                    for rdate in self.rawdates_nsmooth:
                        vector_daily[dates.daily.index((fromjulian(rdate) + timedelta(tshift)).strftime('%Y%j'))] = -1
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

                    #pylint: disable=C0103
                    for ix in map_index:
                        if not isinstance(srange, np.ndarray):
                            lag_correlation = lag1corr(arr_raw[ix, :-1], arr_raw[ix, 1:], nodata)
                            if lag_correlation > 0.5:
                                sr = np.linspace(-2, 1.0, 16)
                                vp = 1
                            elif lag_correlation <= 0.5:
                                sr = np.linspace(0, 3.0, 16)
                                vp = 2
                            else:
                                sr = np.linspace(-1, 1, 11)
                                vp = 5
                        else:
                            sr = srange
                            vp = 5
                        if p:
                            arr_raw[ix, :], arr_sgrid[ix] = ws2doptvp(y=arr_raw[ix, :],
                                                                      w=np.array((arr_raw[ix, :] != nodata)*1, dtype='double'),
                                                                      llas=array('d', sr),
                                                                      p=p,
                                                                      vprec=vp)
                        else:
                            arr_raw[ix, :], arr_sgrid[ix] = ws2doptv(y=arr_raw[ix, :],
                                                                     w=np.array((arr_raw[ix, :] != nodata)*1, dtype='double'),
                                                                     llas=array('d', sr))

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
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(0, arr_smooth.shape[1], smoothchunks[1])):
                            smt_ds[br:br+smoothchunks[0], bcs:bcs+smoothchunks[1]] = arr_smooth[:, bcr:bcr+smoothchunks[1]]
                        arr_smooth[...] = nodata
                    else:
                        for bcs, bcr in zip(range(smoothoffset, smoothshape[1], smoothchunks[1]), range(self.array_offset, arr_raw.shape[1], smoothchunks[1])):
                            smt_ds[br:br+rawchunks[0], bcs:bcs+smoothchunks[1]] = arr_raw[:, bcr:bcr+smoothchunks[1]]
