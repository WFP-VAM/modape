"""IO module for modape"""
# pylint: disable=E0401, C0103
from pathlib import Path

import h5py
import numpy as np

class HDF5Base(object):
    """Base class for interaction with HDF5 files"""
    def __init__(self, filename: str) -> None:
        """initialize class instance.

        Args:
            filename (str): Filename of HDF5 file.
        """

        self.filename = Path(filename)
        self.exists = self.filename.exists()

    def read_chunked(self,
                     dataset: str,
                     xchunk: int = None,
                     arr_out: np.ndarray = None) -> np.ndarray:
        """Read data from dataset in a chunked manner.

        The chunks are iterated in a row by column pattern, where
        each chunk along row axis is yielded once the full
        column size is read into the array. The chunking of
        the colums (xchunk) can me modified, while the row chunking
        (ychunk) is strictly defined by the dataset.


        Args:
            dataset (str): Name of the dataset (expect 2d array).
            xchunk (int): Chunking along row (x) axis.
                          If not specified, it'll be read from the dataset.
            arr_out (np.ndarray): Output array.

        Yields:
            np.ndarray: Yields the output chunk as np.ndarray

        """


        with h5py.File(self.filename, 'r') as h5f_open:

            ds = h5f_open.get(dataset)
            assert ds, f"Dataset '{dataset}' not found!"
            assert ds.ndim == 2, "Expected 2-d array as dataset!"

            ysize, xsize = ds.shape

            if xchunk is None:
                ychunk, xchunk = ds.chunks
            else:
                ychunk, _ = ds.chunks

            if arr_out is None:
                arr_out = np.full((ychunk, xsize), fill_value=ds.fillvalue, dtype=ds.dtype.name)
            else:
                assert isinstance(arr_out, np.ndarray), "arr_out must be 2-d numpy array!"
                assert arr_out.ndim == 2, "Expected 2-d array as output!"

        for yb in range(0, ysize, ychunk):

            with h5py.File(self.filename, 'r') as h5f_open:
                ds = h5f_open.get(dataset)

                for xb in range(0, xsize, xchunk):
                    arr_out[:, xb:(xb+xchunk)] = ds[yb:(yb+ychunk), xb:(xb+xchunk)]

            yield arr_out

    def write_chunk(self,
                    dataset: str,
                    arr_in: np.ndarray,
                    xchunk: int = None,
                    xoff: int = 0,
                    yoff: int = 0) -> bool:
        """Write chunk back to HDF5 file.

        Writes complete chunk back to HDF5 file, iterating
        over the column (x) dimension. The chunksize for x
        can be adjusted manually.

        Args:
            dataset (str): Name of the dataset (expect 2d array).
            arr_in (np.ndarray): Input data to be written to file.
            xchunk (int): Chunking along row (x) axis.
                          If not specified, it'll be read from the dataset.
            xoff (int): Offset for colums (xaxis) in file.
            yoff (int): Offset for rows (yaxis) in file.

        Returns:
            bool: True if write successful.

        """

        assert isinstance(arr_in, np.ndarray), "arr_in must be 2-d numpy array!"
        assert arr_in.ndim == 2, "Expected 2-d array as input!"

        with h5py.File(self.filename, 'r+') as h5f_open:

            ds = h5f_open.get(dataset)
            assert ds, f"Dataset '{dataset}' not found!"
            assert ds.ndim == 2, "Expected 2-d array as dataset!"

            if xchunk is None:
                ychunk, xchunk = ds.chunks
            else:
                ychunk, _ = ds.chunks

            ysize, xsize = arr_in.shape

            ysize = max(ysize, ychunk)
            assert ysize <= ychunk

            for xb in range(0, xsize, xchunk):

                xb_data = xb + xoff
                ds[yoff:(yoff+ysize), xb_data:(xb_data+xchunk)] = arr_in[:, xb:(xb+xchunk)]

        return True
