"""
MODIS window functions & classes.

This file contains the classes and functions for extracting GeoTIFFs from HDF5 files.

Author: Valentin Pesendorfer, April 2019
"""
# pylint: disable=import-error,R0903,W0706
from contextlib import contextmanager
import datetime
import logging
from os.path import exists
from pathlib import Path
from typing import Union, List
from uuid import uuid4

from osgeo import gdal
import h5py
from modape.constants import DATE_LABELS
from modape.exceptions import HDF5MosaicError
from modape.utils import fromjulian
import numpy as np

log = logging.getLogger(__name__)

class ModisMosaic(object):
    """Class for subsetting or mosaicing HDF5 files.

    This class can be used to mosaic multiple raw or smooth HDF5 files or
    to subset a single HDF5 file (both tiled and global products).
    The generated data will exported as GeoTIFF using the Python GDAL bindings.
    """

    def __init__(self,
                 input_files: List[str]) -> None:
        """Initialize instance of `ModisMosaic` class.

        If multiple HDF5 files are specified as `input_files`, a mosaic of the files
        will be created (need to be of same product and tile/global). It's required that all of the files share the same temporal axis.
        If only one file is specified, either the full file is exported or a subset is created.

        Args:
            input_files (List[str]): List of paths to HDF5 files for mosaicing.

        Raises:
            HDF5MosaicError: If temporal axis is not conistent with multiple files (i.e. contained dates don't line up)

        """
        # check if files are temporal coherernt
        if isinstance(input_files, list):

            reference = None
            for file in input_files:
                with h5py.File(file, "r") as h5f_open:
                    dates = h5f_open.get("dates")[...]
                    if reference is None:
                        reference = dates
                    try:
                        np.testing.assert_array_equal(dates, reference)
                    except AssertionError:
                        raise HDF5MosaicError("Teporal axis of HDF5 input files incompatible for Mosaic")

            self.files = input_files
        else:
            self.files = [input_files]
            with h5py.File(input_files, "r") as h5f_open:
                dates = h5f_open.get("dates")[...]
        self.dates = [fromjulian(x.decode()) for x in dates]

    def generate_mosaics(self,
                         dataset,
                         targetdir,
                         target_srs: str = None,
                         aoi: List[float] = None,
                         overwrite: bool = False,
                         force_doy: bool = False,
                         prefix: str = None,
                         start: datetime.date = None,
                         stop: datetime.date = None,
                         clip_valid: bool = False,
                         round_int: int = None,
                         last_smoothed: str = None,
                         **kwargs,
                        ) -> List:
        """Generate GeoTIFF mosaics/subsets.

        This method is creating a GeoTiff mosaic/subsets from the HDF5 files
        passed to the class instance.
        Internally, a virtual raster for each timestep is created,
        warped to the desired SRS (`target_srs`) and then optionally clipped to the
        area of interest (`aoi` as list of corner coordinates). If the dataset is already in the desired SRS,
        the warping is skipped. To transform the resolution of meters to degrees in `EPSG:4326`,
        the original resolution of the dataset in meters is divided by `112000`.

        The exported dates can be limited with `start` and `stop.` By default, existing GeoTIFFs are
        not overwritten. This can be forced with setting `overwrite` to `True`.

        The data to be exported can be clipped to the valid data range of the MODIS product being
        processed (if applicable) to the setting `clip_valid` to `True`. The data can also be rounded
        by passing an integer to `round_int`. This integer will be passed on to `np.round`.

        To modify the output naming, a `prefix` can be specified and if `force_doy` is set to `True`,
        the timestamp will be forced to the format `YYYYjDDD` independent of temporal interpolation.

        More fine scaled adjustment of the output datasets can be achieved with passing keyword arguments
        to `gdal.WarpOptions` and `gdal.TranslateOptions`.

        Args:
            dataset (type): Dataset to mosaic (data or sgrid).
            targetdir (type): Target directory for output Tiffs.
            target_srs (str): Target spatial reference (as expected by gdalwarp).
            aoi (List[float]): AOI bouning box as ULX,ULY,LRX,LRY.
            overwrite (bool): Flag for overwriting existing Tiffs.
            force_doy (bool): Force DOY filenaming instead of VAM labels.
            prefix (str): Perfix for Tiff filenames.
            start (datetime.date): Start date for mosaics.
            stop (datetime.date): Stop date for mosaics.
            clip_valid (bool): Clip values to valid range for MODIS product.
            round_int (int): Round the output.
            last_smoothed (str): Rawdate (MODIS time step) that is checked to be the last in series at time of smoothing.
            **kwargs (type): **kwags passed on to `gdal.WarpOptions` and `gdal.TranslateOptions`.

        Raises:
            ValueError: If dataset supplied does not exists in files.
            AssertionError: If write fails.

        Returns:
            mosaics: List of mosaiced raster datasets

        """

        try:
            attrs = self._get_metadata(self.files[0], dataset)
        except AssertionError:
            raise ValueError(f"Supplied dataset {dataset} doesn't exist.")

        if prefix is None:
            prefix = ""

        try:
            labels = DATE_LABELS[attrs["temporalresolution"]]
        except KeyError:
            labels = None
            force_doy = True

        if "xRes" in kwargs and "yRes" in kwargs:
            output_res = [float(kwargs["xRes"]), float(kwargs["yRes"])]
            del kwargs["xRes"]
            del kwargs["yRes"]

        elif target_srs == "EPSG:4326":

            try:
                if not attrs["globalproduct"]:
                    output_res = attrs["resolution"] / 112000
                else:
                    output_res = attrs["resolution"]
            except KeyError:
                log.warning("Could not determine target resolution from file!")
                output_res = [None, None]
        else:
            output_res = [None, None]

        try:
            nodata = kwargs["noData"]
        except KeyError:
            nodata = attrs["nodata"]

        try:
            dtype = kwargs["outputType"]
        except KeyError:
            dtype = attrs["dtype"]

        try:
            resample = kwargs["resampleAlg"]
        except KeyError:
            resample = "near"

        if "multithread" in kwargs:
            gdal_multithread = bool(kwargs["multithread"])
            del kwargs["multithread"]
        else:
            gdal_multithread = False

        filename_root = f"{targetdir}/{prefix}{attrs['vamcode'].lower()}"

        if len(attrs["dataset_shape"]) == 1:
            date_index = [None]

        else:
            if start is None:
                start = self.dates[0]

            if stop is None:
                stop = self.dates[-1]

            date_index = [x for x, y in enumerate(self.dates) if start <= y <= stop]

        mosaics = []
        for ii in date_index:

            if ii is None:
                assert dataset == "sgrid"
                filename = filename_root + "_sgrid.tif"

            else:

                date = self.dates[ii]

                if force_doy:
                    filename = filename_root + date.strftime("%Yj%j") + ".tif"
                else:
                    filename = f"{filename_root}{date.year}{date.month:02}{labels[date.day]}.tif"

                if exists(filename):
                    if not overwrite:
                        log.info("%s already exists. Please set overwrite to True.")
                        continue

            log.info("Processing %s", filename)

            rasters = []
            for file in self.files:
                rasters.append(
                    self._get_raster(
                        file,
                        dataset,
                        clip_valid,
                        round_int=round_int,
                        ix=ii,
                        last_smoothed=last_smoothed
                    )
                )

            translate_options = {
                "outputType": attrs["dtype"],
                "noData": attrs["nodata"],
                "outputSRS": target_srs,
                "projWin": aoi
            }

            translate_options.update(kwargs)

            if aoi is not None and all(output_res):
                translate_options.update({
                    "outputBounds": aoi,
                    "width": abs(int(round((aoi[2] - aoi[0]) / output_res[0]))),
                    "height": abs(int(round((aoi[3] - aoi[1]) / output_res[1])))
                })

            with self._mosaic(rasters,
                              target_srs=target_srs,
                              resample=resample,
                              dtype=dtype,
                              nodata=nodata,
                              resolution=output_res,
                              gdal_multithread=gdal_multithread,
                             ) as warped_mosaic:

                log.debug("Writing to disk")

                write_check = self._translate(
                    src=warped_mosaic,
                    dst=filename,
                    **translate_options
                )

                try:
                    assert write_check, f"Error writing {filename}"
                    mosaics.append(filename)
                except:
                    raise
                finally:
                    _ = [gdal.Unlink(x) for x in rasters]

        return mosaics

    @staticmethod
    def _get_raster(file: Union[Path, str],
                    dataset: str,
                    clip_valid: bool,
                    round_int: int,
                    ix: int = None,
                    last_smoothed: str = None) -> str:

        if dataset not in ["data", "sgrid"]:
            raise NotImplementedError("_get_raster only implemented for datasetds 'data' and 'sgrid'")

        with h5py.File(file, "r") as h5f_open:

            if last_smoothed is not None:
                dates = h5f_open.get("rawdates")
                last_date = dates[-1].decode()
                assert last_smoothed == last_date, \
                    f"Last smoothed date in {file} is {last_date} not {last_smoothed}"
                
            ds = h5f_open.get(dataset)
            assert ds, "Dataset doesn't exist!"
            dataset_shape = ds.shape

            if len(dataset_shape) < 2:
                sgrid = True
                attrs = dict(h5f_open.get("data").attrs)

            else:
                if ix is None:
                    raise ValueError("Need index for 2-d dataset!")
                sgrid = False
                attrs = dict(ds.attrs)

            chunks = ds.chunks
            dtype = ds.dtype.num

            driver = gdal.GetDriverByName("GTiff")
            fn = f"/vsimem/{uuid4()}.tif"

            raster = driver.Create(
                fn,
                int(attrs["RasterXSize"]),
                int(attrs["RasterYSize"]),
                1,
                int(dtype),
            )
            raster.SetGeoTransform(attrs["geotransform"])
            raster.SetProjection(attrs["projection"])

            raster_band = raster.GetRasterBand(1)

            if not sgrid:
                raster_band.SetNoDataValue(int(attrs["nodata"]))

            block_gen = ((x, x//attrs["RasterXSize"]) for x in range(0, dataset_shape[0], chunks[0]))

            for yblock_ds, yblock in block_gen:

                if sgrid:
                    arr = ds[yblock_ds:(yblock_ds+chunks[0])]
                else:
                    arr = ds[yblock_ds:(yblock_ds+chunks[0]), ix]

                if clip_valid:
                    vmin, vmax = attrs["valid_range"]
                    arr = np.clip(arr, vmin, vmax, out=arr, where=arr != attrs["nodata"])

                if round_int is not None:
                    arr = np.round(arr, round_int)

                raster_band.WriteArray(
                    arr.reshape(-1, int(attrs["RasterXSize"])),
                    xoff=0,
                    yoff=int(yblock)
                )

            raster_band.FlushCache()
            raster.FlushCache()
            raster_band, raster, driver = (None, None, None)

            return fn

    @staticmethod
    @contextmanager
    def _mosaic(input_rasters,
                target_srs,
                resample,
                resolution,
                dtype,
                nodata,
                gdal_multithread):

        vrt_tempname = f"/vsimem/{uuid4()}.vrt"
        vrt = gdal.BuildVRT(vrt_tempname, input_rasters)
        assert vrt
        vrt.FlushCache()
        vrt = None
        log.debug("Created VRT")

        wrp_tempname = f"/vsimem/{uuid4()}.tif"
        wopt = gdal.WarpOptions(
            dstSRS=target_srs,
            outputType=dtype,
            resampleAlg=resample,
            xRes=resolution[0],
            yRes=resolution[1],
            srcNodata=nodata,
            dstNodata=nodata,
            multithread=gdal_multithread,
        )

        wrp = gdal.Warp(wrp_tempname, vrt_tempname, options=wopt)
        assert wrp

        log.debug("Created WRP")

        yield wrp

        log.debug("Performing Cleanup")

        wrp = None
        vrt = None

        rc1 = gdal.Unlink(vrt_tempname)
        rc2 = gdal.Unlink(wrp_tempname)
        if rc1 != 0 or rc2 != 0:
            log.warning("Received return codes [%s, %s] while removing MemRasters", rc1, rc2)

    @staticmethod
    def _translate(src, dst, **kwargs):

        topt = gdal.TranslateOptions(**kwargs)

        ds = gdal.Translate(dst, src, options=topt)
        assert ds
        ds = None
        return True

    @staticmethod
    def _get_metadata(reference, dataset):
        with h5py.File(reference, "r") as h5f_open:
            ds = h5f_open.get(dataset)
            assert ds

            attrs = dict(h5f_open.get("data").attrs)
            attrs["dataset_shape"] = ds.shape
            attrs["nodata"] = ds.fillvalue
            attrs["dtype"] = gdal.GetDataTypeByName(ds.dtype.name)

        return attrs
