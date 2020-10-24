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
    """Class for creating mosaics from multiple HDF5 files"""

    def __init__(self,
                 input_files: List[str]) -> None:
        """Create instance

        Args:
            input_files (List[str]): List of paths to HDF5 files for mosaicing.

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
                         target_srs: str,
                         aoi: List[float] = None,
                         overwrite: bool = False,
                         force_doy: bool = False,
                         prefix: str = None,
                         start: datetime.date = None,
                         stop: datetime.date = None,
                         clip_valid: bool = False,
                         round_int: int = None,
                         **kwargs,
                        ) -> list:
        """Generate TIFF mosaics.

        This method is creating a GeoTiff mosaic from the MDF5 files
        passed to the class instance.
        Internally, a virtual raster for each timestep is created,
        warped to the desired SRS and then optionally clipped to the
        area of interest.

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
            **kwargs (type): **kwags passed on to `gdal.WarpOptions` and `gdal.TranslateOptions`.

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

        if "xRes" in kwargs.keys() and "yRes" in kwargs.keys():
            output_res = [kwargs["xRes"], kwargs["yRes"]]

        elif target_srs == "EPSG:4326":

            if not attrs["globalproduct"]:
                output_res = attrs["resolution"] / 112000
            else:
                output_res = attrs["resolution"]

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
                    )
                )

            translate_options = {
                "outputType": attrs["dtype"],
                "noData": attrs["nodata"],
                "outputSRS": target_srs,
            }

            if aoi is not None:
                translate_options.update(
                    {"projWin": aoi}
                )

            translate_options.update(kwargs)

            if not attrs["globalproduct"] or target_srs != "EPSG:4326":

                with self._mosaic(rasters,
                                  target_srs=target_srs,
                                  resample=resample,
                                  dtype=dtype,
                                  nodata=nodata,
                                  resolution=output_res,
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


            else:

                log.debug("Processing global file with EPSG:4326! Skipping warp")
                assert len(rasters) == 1, "Expected only one raster!"

                log.debug("Writing to disk")

                write_check = self._translate(
                    src=rasters[0],
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
                    ix: int = None) -> str:

        if dataset not in ["data", "sgrid"]:
            raise NotImplementedError("_get_raster only implemented for datasetds 'data' and 'sgrid'")

        with h5py.File(file, "r") as h5f_open:
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
                nodata):

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
            multithread=True,
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
