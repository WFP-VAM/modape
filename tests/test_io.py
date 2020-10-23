"""test_modis.py: Test modape I/O and functions."""
# pylint: disable=E0401, W0702
import shutil
from tempfile import NamedTemporaryFile
import unittest

import h5py
import numpy as np

from modape.modis.io import HDF5Base

class TestHDF5io(unittest.TestCase):
    """Test class for HDF5 io."""

    @classmethod
    def setUpClass(cls):
        """Set up testing class"""

        with NamedTemporaryFile(suffix=".h5", delete=False) as temp:

            with h5py.File(temp.name, "a") as h5f:
                _ = h5f.create_dataset(
                    "data",
                    shape=(100, 100),
                    dtype="int16",
                    maxshape=(100, None),
                    chunks=(10, 1),
                    compression="gzip",
                    fillvalue=-1
                )

                _ = h5f.create_dataset(
                    "sgrid",
                    shape=(100,),
                    dtype="int16",
                    maxshape=(100,),
                    chunks=(10,),
                    compression="gzip",
                    fillvalue=-1
                )

                _ = h5f.create_dataset(
                    "dates",
                    shape=(100,),
                    maxshape=(None,),
                    dtype="int16",
                    compression="gzip")

            cls.testfile = temp.name

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree(cls.testfile)
        except:
            pass

    def test_rw(self):
        """Test read and write to HDF5 file"""

        test_array = np.arange(0, 10000, dtype="int16").reshape(100, 100)
        hdf5_file = HDF5Base(self.testfile)

        for ii in range(0, 100, 10):
            arr_sub = test_array[ii:ii+10, :]
            hdf5_file.write_chunk("data", arr_in=arr_sub, yoffset=ii)

        ii = 0
        for read_array in hdf5_file.read_chunked("data", xchunk=10):
            jj = ii * 10
            np.testing.assert_array_equal(read_array, test_array[jj:(jj+10)])
            ii += 1

        del test_array, read_array #pylint: disable=W0631

        test_array = np.full((10, 100), -1, dtype="int16")

        for read_array in hdf5_file.read_chunked("data", xchunk=10, arr_out=test_array):
            assert read_array is test_array
            assert np.alltrue(read_array != -1)

        # test atomatic array size for nsmooth
        # nsmooth = 36
        xoffset = 100 - 36
        for read_array in hdf5_file.read_chunked("data", xoffset=xoffset, xchunk=10):
            assert read_array.shape == (10, 36)
            break

        del read_array

        # test 1d read/write
        test_array = np.arange(100, dtype="int16")
        for ii in range(0, 100, 10):
            arr_sub = test_array[ii:ii+10]
            hdf5_file.write_chunk("sgrid", arr_in=arr_sub, yoffset=ii)

        ii = 0
        for read_array in hdf5_file.read_chunked("sgrid"):
            jj = ii*10
            np.testing.assert_array_equal(
                read_array,
                test_array[jj:(jj+10)]
            )
            ii += 1

    def test_assertions(self):
        """Test assertions in class"""
        hdf5_file = HDF5Base(self.testfile)

        with self.assertRaises(AssertionError):
            next(hdf5_file.read_chunked("not_a_dataset", xchunk=10))

        with self.assertRaises(AssertionError):
            next(hdf5_file.read_chunked("dates", xchunk=10, arr_out="not_an_array"))
