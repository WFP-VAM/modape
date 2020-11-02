"""IO module for modape"""
# pylint: disable=E0401, C0103
from contextlib import contextmanager
import logging
from pathlib import Path
from typing import List, Tuple

from osgeo import gdal
import h5py
import numpy as np

log = logging.getLogger(__name__)

class HDF5Base(object):
    """Parent class for interaction with HDF5 files

    This class serves as a parent class for ModisRawH5 and ModisSmoothH5,
    enabling uniform chunked read and write of datasets and attributes to and from HDF5 files."""
    def __init__(self, filename: str) -> None:
        """Initialize HDF5Base instance.

        Creates an instance of the HDF5Base class. This is not stricly intended to be called
        by itself (although it can be), but rather through a `super` call of the child class.

        Args:
            filename (str): Full path to the HDF5 file.
        """

        self.filename = Path(filename)
        self.exists = self.filename.exists()

    def read_chunked(self,
                     dataset: str,
                     xoffset: int = 0,
                     xchunk: int = None,
                     arr_out: np.ndarray = None) -> np.ndarray:
        """Read data from dataset in a chunked manner.

        The chunks are iterated in a row by column pattern, where
        each chunk along row axis is yielded once the full
        column size is read into the array. The chunking of
        the colums (`xchunk`) can me modified, while the row chunking
        (`ychunk`) is strictly defined by the dataset.
        To enable the `nsmooth` functionality, an `xoffset` can be provided
        to skip datapoints along the time dimension.
        If no `arr_out` is provided, a new array will be created with the
        necessary dimensions.


        Args:
            dataset (str): Name of the dataset (expect 2d array).
            xoffset (int): Offset for time-dimension (xaxis) in file.
            xchunk (int): Chunking along time-dimension.
                          If not specified, it'll be read from the dataset.
            arr_out (np.ndarray): Output array.

        Returns:
            Yields the output chunk as np.ndarray

        Raises:
            AssertionError: Raised when `dataset` not found within HDF5 file and when provided `arr_out` is not correct object or shape

        """

        with h5py.File(self.filename, "r") as h5f_open:

            ds = h5f_open.get(dataset)
            assert ds, f"Dataset '{dataset}' not found!"
            ds_shape = ds.shape

            ychunk = ds.chunks[0]

            if arr_out is None:
                if len(ds_shape) == 1:
                    arr_out = np.full((ychunk,), fill_value=ds.fillvalue, dtype=ds.dtype.name)
                else:
                    arr_out = np.full((ychunk, ds_shape[1]-xoffset), fill_value=ds.fillvalue, dtype=ds.dtype.name)

            else:
                assert isinstance(arr_out, np.ndarray)
                assert arr_out.ndim <= 2, "Expected 1 or 2-d array as output!"

            if arr_out.ndim == 1:
                xsize = None
                xchunk = None
            else:
                xsize = arr_out.shape[1]
                if xchunk is None:
                    xchunk = ds.chunks[1]

        log.debug("arr_out shape: %s", arr_out.shape)

        for yb in range(0, ds_shape[0], ychunk):

            with h5py.File(self.filename, "r") as h5f_open:
                ds = h5f_open.get(dataset)
                log.debug("Reading chunk %s - %s", yb, yb+ychunk)
                if xsize is not None:
                    for xb in range(0, xsize, xchunk):
                        xb_data = xb + xoffset

                        log.debug("Reading arr_out[%s : %s, %s : %s] from dataset[:, %s : %s]", yb, yb+ychunk, xb, xb+xchunk, xb_data, xb_data+xchunk)
                        arr_out[:, xb:(xb+xchunk)] = ds[yb:(yb+ychunk), xb_data:(xb_data+xchunk)]
                else:
                    arr_out[...] = ds[yb:(yb+ychunk)]

            yield arr_out

    def write_chunk(self,
                    dataset: str,
                    arr_in: np.ndarray,
                    xoffset: int = 0,
                    xchunk: int = None,
                    yoffset: int = 0) -> bool:
        """Write chunk back to HDF5 file.

        Writes complete chunk back to HDF5 file, iterating
        over the time-dimension (x). The chunksize for x
        can be adjusted manually using `xchunk`.
        To implement `nupdate` behaviour, `xoffset` can be used to skip
        prior datapoints in the time-dimension.
        To write successive spatial chunks, `yoffset` has to be provided (the default is 0,
        as it starts at the top left).

        Args:
            dataset (str): Name of the dataset (expect 2d array).
            arr_in (np.ndarray): Input data to be written to file.
            xchunk (int): Chunking along row (x) axis.
                          If not specified, it'll be read from the dataset.
            xoffset (int): Offset for colums (xaxis) in file.
            yoffset (int): Offset for rows (yaxis) in file.

        Returns:
            Returns `True` if write was successful.

        Raises:
            AssertionError: Raised when `dataset` not found within HDF5 file and when provided `arr_in` is not correct object or shape


        """

        assert isinstance(arr_in, np.ndarray), "arr_in must be 2-d numpy array!"
        assert arr_in.ndim <= 2, "Expected 1-d or 2-d array as input!"

        with h5py.File(self.filename, "r+") as h5f_open:

            ds = h5f_open.get(dataset)
            assert ds, f"Dataset '{dataset}' not found!"
            assert ds.ndim <= 2, "Expected 1-d or 2-d array as dataset!"

            ychunk = ds.chunks[0]

            try:
                ysize, xsize = arr_in.shape
                if xchunk is None:
                    xchunk = xsize
            except ValueError:
                ysize, = arr_in.shape
                xsize = None

            ysize = max(ysize, ychunk)
            assert ysize <= ychunk

            log.debug("arr_in shape: %s", arr_in.shape)

            if xsize is not None:

                for xb in range(0, xsize, xchunk):
                    xb_data = xb + xoffset
                    log.debug("Writing dataset[%s : %s, %s : %s] from arr_in[:, %s : %s]", yoffset, yoffset+ysize, xb_data, xb_data+xchunk, xb, xb+xchunk)
                    ds[yoffset:(yoffset+ysize), xb_data:(xb_data+xchunk)] = arr_in[:, xb:(xb+xchunk)]
            else:
                log.debug("Writing dataset[%s : %s] from arr_in", yoffset, yoffset+ysize)
                ds[yoffset:(yoffset+ysize)] = arr_in[...]

        return True

    @staticmethod
    def _get_reference_metadata(reference_file: str,
                                sds_filter: str = None)-> dict:
        """Helper function to get metadata from reference file.

        Extracts metadata from subdataset, eitjer filtered
        by sds_filter or if None, using the first one.

        GDAL's geotransform is extracted following the
        Affine scheme:
            c = x-coordinate of the upper-left corner of the upper-left pixel
            a = width of a pixel
            b = row rotation (typically zero)
            f = y-coordinate of the of the upper-left corner of the upper-left pixel
            d = column rotation (typically zero)
            e = height of a pixel (typically negative)

        Args:
            reference_file (str): Full path to reference file.
            sds_filter (str): Indicator of Subdataset 9as returned by VAM_PRODUCT_CODES.

        Returns:
            dict: Dictionary with metadata

        """

        ds = gdal.Open(reference_file)
        sds_all = ds.GetSubDatasets()

        if sds_filter is None:
            sds = sds_all[0][0]
        else:
            sds = [x[0] for x in sds_all if sds_filter in x[0]][0]

        sds_open = gdal.Open(sds)
        c, a, b, f, d, e = sds_open.GetGeoTransform()

        metadata = dict(
            RasterXSize=sds_open.RasterXSize,
            RasterYSize=sds_open.RasterYSize,
            geotransform=(c, a, b, f, d, e),
            projection=sds_open.GetProjection(),
            resolution=(a, e),
            nodata=int(sds_open.GetMetadataItem("_FillValue")),
        )

        sds_open = None
        ds = None

        return metadata

class HDFHandler(object):
    """Class to handle reading from MODIS HDF files.

    This class enables reading specific subdatasets and attributes
    from the raw MODIS HDF files."""
    def __init__(self,
                 files: List[str],
                 sds: str) -> None:
        """Initialize HDFHandler instance.

        Reads the datasets, extracts the subdatasets and keeps
        a reference to the file handlers.

        Args:
            files (List[str]): List of file paths to hdf files.
            sds (str): subdataset as returned by `modape.constants.VAM_PRODUCT_CODES`.

        """

        self.files = files
        self.sds = sds
        self.handles = []

    @contextmanager
    def open_datasets(self) -> None:
        """Opens the selected subdataset from all files
        within a context manager and stores them in a class variable.
        When the context manager closes, the refereces are removed, closing
        all datasets.

        """
        for ds_handle in self._gen_sds_handle(self.files, self.sds):
            self.handles.append(ds_handle)

        yield

        for ii in range(len(self.handles)):
            self.handles[ii] = None
        self.handles = []

    def iter_handles(self) -> Tuple[int, "gdal.Dataset"]:
        """Iterates over all open dataset handles
        coming from `open_datasets` and returns a Tuple with index and a
        `gdal.Dataset` for each.

        Returns:
            Tuple with index and corresponding `gdal.Dataset`

        """
        ix = 0
        for handle in self.handles:
            yield (ix, handle)
            ix += 1

    @staticmethod
    def read_chunk(x: "gdal.Dataset", **kwargs: dict) -> np.ndarray:
        """Reads a chunk of an opened subdataset.

        The size of the chunk being read is defined by the
        values passed to `gdal.Dataset.ReadAsArray` with `kwargs` and can
        be as big as the entire dataset.

        Args:
            x (gdal.Dataset): GDAL Dataset.
            **kwargs (dict): kwargs passed on to `gdal.Dataset.ReadAsArray`

        Returns:
            Requested chunk as `np.ndarray`

        """
        return x.ReadAsArray(**kwargs)

    @staticmethod
    def _gen_sds_handle(x: str, sds: str):
        for xx in x:
            try:
                ds = gdal.Open(xx)
                ds_sds = [x[0] for x in ds.GetSubDatasets() if sds in x[0]][0]
                yield gdal.Open(ds_sds)
                ds = None

            except AttributeError:
                yield None
