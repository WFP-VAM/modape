"""test_modis.py: Test MODIS classes and functions."""
# pylint: disable=E0401,E0611,W0702,W0613,C0103
from datetime import datetime
from pathlib import Path, PosixPath
import pickle
import shutil
from types import SimpleNamespace
import unittest
from unittest.mock import patch, MagicMock
from uuid import uuid4

import numpy as np
import h5py #pylint: disable=import-error
from osgeo import gdal

from modape.exceptions import DownloadError, HDF5WriteError
from modape.modis import ModisQuery, ModisRawH5, ModisSmoothH5, ModisMosaic
from modape.utils import SessionWithHeaderRedirection

class MockResponse:
    '''Mock response for testing'''
    def __init__(self, content, status_code):
        '''Create instance, setting content and status_code'''
        self._content = content
        self.status_code = status_code

    @property
    def content(self):
        '''Return content'''
        return self._content

    def raise_for_status(self):
        '''don't raise for status'''
        pass #pylint: disable=W0107

class MockedPath(PosixPath):
    '''Mocked version of PosixPath'''

    # file size for testing
    _filesize = 7526571
    def __init__(self, filename):
        super().__init__()

    @property
    def filesize(self):
        '''Get file size'''
        return self._filesize

    def is_dir(self):
        return True

    def exists(self):
        return True

    def stat(self):
        return SimpleNamespace(st_size=self.filesize, st_mode=33188)

def create_gdal(x, y):
    """Create in-memory gdal dataset for testing.

    Returns:
        GDAL dataset object
    """
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create("/vsimem/{}.tif".format(str(uuid4())), x, y, 1, 3)
    return ds

def create_h5temp(rows: int,
                  cols: int,
                  tr: int,
                  ts: int,
                  ) -> Path:
    """Create temporary HDF5 rawfile.

    Args:
        rows (int): Number of rows.
        cols (int): Number of columns.
        tr (int): Temporal resolution of rawfile.
        ts (int): Temporal shift.

    Returns:
        Path: Path object of filename.

    """

    fn = Path("/tmp/data/MXD13A2.h21v10.006.VIM.h5")


    with h5py.File(fn, "a", driver="core", backing_store=True) as h5f:

        dset = h5f.create_dataset("data",
                                  shape=(rows*cols, 4),
                                  dtype="int16",
                                  maxshape=(rows*cols, None),
                                  chunks=((rows*cols)//25, 10),
                                  compression="gzip",
                                  fillvalue=-3000)

        dset.attrs.update(
            dict(
                nodata=-3000,
                temporalresolution=tr,
                tshift=ts,
                RasterXSize=rows,
                RasterYSize=cols,
                geotransform=(0, 0, 0, 0, 0),
                projection="EPSG:4326",
                resolution=(1000, -1000),
                globalproduct=False,
                vamcode="VIM",
                )
        )

        h5f.create_dataset("dates",
                           shape=(4,),
                           data=np.array(["2002185", "2002193", "2002201", "2002209"], dtype="S8"),
                           maxshape=(None,),
                           dtype="S8",
                           compression="gzip")

    return fn

def create_h5temp_global() -> Path:
    """Create temporary global HDF5 rawfile.

    Returns:
        Path: Path object of filename.

    """

    fn = Path("/tmp/data/MXD13A2.006.VIM.h5")

    rows = 2000
    cols = 7200

    with h5py.File(fn, "a", driver="core", backing_store=True) as h5f:

        dset = h5f.create_dataset("data",
                                  shape=(rows*cols, 4),
                                  dtype="int16",
                                  maxshape=(rows*cols, None),
                                  chunks=((rows*cols)//25, 10),
                                  compression="gzip",
                                  fillvalue=-3000)

        dset.attrs.update(
            dict(
                nodata=-3000,
                temporalresolution=8,
                tshift=8,
                RasterXSize=rows,
                RasterYSize=cols,
                geotransform=(0, 0, 0, 0, 0),
                projection="EPSG:4326",
                resolution=(0.05, 0.05),
                globalproduct=True,
                vamcode="VIM",
                )
        )

        h5f.create_dataset("dates",
                           shape=(4,),
                           data=np.array(["2002185", "2002193", "2002201", "2002209"], dtype="S8"),
                           maxshape=(None,),
                           dtype="S8",
                           compression="gzip")

    return fn

class TestModisQuery(unittest.TestCase):
    """Test class for ModisQuery tests."""

    @classmethod
    def setUpClass(cls):
        """Set up testing class"""

        data_dir = str(Path(__file__).parent).replace("tests", "modape")

        with open(f"{data_dir}/data/cmr_api_response.pkl", "rb") as pkl:
            cls.api_response = pickle.load(pkl)

        with open(f"{data_dir}/data/hdf.xml", "rb") as fl:
            cls.hdfxml = fl.read()

        cls.testpath = Path(__name__).parent
        cls.query = ModisQuery(
            products=["MOD13A2", "MYD13A2"],
            aoi=(10, 10, 20, 20),
            begindate=datetime(2020, 1, 1),
            enddate=datetime(2020, 7, 24),
        )

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree("__pycache__")
        except:
            pass

    def test_query(self):
        """Test API query"""

        self.assertEqual(self.query.api.params["short_name"], ["MOD13A2", "MYD13A2"])
        self.assertEqual(self.query.api.params["bounding_box"], "10.0,10.0,20.0,20.0")
        self.assertEqual(self.query.api.params["temporal"], ["2020-01-01T00:00:00Z,2020-07-24T00:00:00Z"])

    def test_response_parse(self):
        """Test parsing of response"""

        with patch("modape.modis.download.GranuleQuery.get_all",
                   return_value=self.api_response):

            self.query.search(match_begin=False)

            self.assertTrue(self.query.results)
            self.assertEqual(len(self.query.results), len(self.api_response))

            tiles = list({values["tile"] for key, values in self.query.results.items()})
            tiles_select = tiles[:2]
            self.assertEqual(len(tiles), 6)

            self.setUp()
            self.query.tile_filter = tiles_select
            self.query.search(match_begin=True)

            self.assertEqual(len({values["tile"] for key, values in self.query.results.items()}), 2)
            self.assertTrue(all([values["time_start"] >= self.query.begin.date() for key, values in self.query.results.items()]))
            self.assertTrue(all([values["time_end"] <= self.query.end.date() for key, values in self.query.results.items()]))

    @patch("modape.modis.download.cksum", return_value=1534015008)
    @patch("modape.modis.download.shutil")
    @patch("modape.modis.download.open")
    @patch("modape.modis.download.ThreadPoolExecutor")
    @patch("modape.modis.download.GranuleQuery.get_all")
    def test_download(self, mock_response, mock_submit, mock_open, mock_shutil, mock_cksum):
        """Test download"""

        mock_response.return_value = self.api_response
        self.query.search()

        for key, values in self.query.results.items():
            test_results = {key: values}
            break
        self.query.results = test_results

        # Successful downloads
        fid = next(iter(test_results))

        future_result = (fid, None)
        mock_submit.return_value.__enter__.return_value.submit.return_value.result.return_value = future_result

        with patch("modape.modis.download.ModisQuery._fetch", return_value=future_result) as mocked_fetch:
            self.query.download(
                targetdir=(self.testpath),
                username="test",
                password="test",
                multithread=False,
            )

            mocked_fetch.assert_called_once()
            fetch_args = mocked_fetch.call_args[0]

            self.assertEqual(type(fetch_args[0]), SessionWithHeaderRedirection)
            self.assertEqual(fetch_args[1], self.query.results[fid]["link"])
            self.assertEqual(fetch_args[2], self.testpath)

        # test retry and error
        mock_submit.reset_mock()
        with patch("modape.modis.download.SessionWithHeaderRedirection"):
            future_result = (f"http://datalocation.com/{fid}", "Error")
            mock_submit.return_value.__enter__.return_value.submit.return_value.result.return_value = future_result

            with self.assertRaises(DownloadError):
                self.query.download(
                    targetdir=self.testpath,
                    username="test",
                    password="test",
                    multithread=True,
                    max_retries=5
                )

        # 1 + 5 retriess
        self.assertEqual(mock_submit.call_count, 6)

        # check robust download
        with patch("modape.modis.download.SessionWithHeaderRedirection") as session_mock:

            xml_mock = MockResponse(self.hdfxml, 200)
            session_mock.return_value.__enter__.return_value.get.return_value.__enter__.side_effect = [
                MagicMock(),
                xml_mock,
            ]

            mocked_path = MockedPath(self.testpath)

            dl = self.query.download(
                targetdir=mocked_path,
                username="test",
                password="test",
                multithread=False,
                max_retries=5,
                robust=True,
            )

            self.assertEqual(dl, ["MOD13A2.A2020001.h18v07.006.2020018001022.hdf"])

class TestModisCollect(unittest.TestCase):
    """Test class for ModisQuery tests."""

    @classmethod
    def setUpClass(cls):
        """Set up testing class"""


        cls.vim_files = ["MYD13A2.A2002201.h18v06.006.2015149071105.hdf",
                         "MYD13A2.A2002185.h18v06.006.2015149071113.hdf",
                         "MOD13A2.A2002177.h18v06.006.2015149001129.hdf",
                         "MOD13A2.A2002209.h18v06.006.2015149180726.hdf",
                         "MOD13A2.A2002193.h18v06.006.2015149022847.hdf"]

        cls.vim_files_terra = [x for x in cls.vim_files if "MOD13A2" in x]
        cls.vim_files_aqua = [x for x in cls.vim_files if "MYD13A2" in x]

        cls.lst_files = ["MYD11A2.A2002193.h18v06.006.2015146152945.hdf",
                         "MOD11A2.A2002209.h18v06.006.2015145152020.hdf",
                         "MYD11A2.A2002201.h18v06.006.2015146153241.hdf",
                         "MYD11A2.A2002185.h18v06.006.2015146152642.hdf",
                         "MOD11A2.A2002177.h18v06.006.2015144183717.hdf",
                         "MYD11A2.A2002209.h18v06.006.2015152152813.hdf",
                         "MOD11A2.A2002185.h18v06.006.2015145002847.hdf",
                         "MOD11A2.A2002193.h18v06.006.2015145055806.hdf",
                         "MOD11A2.A2002201.h18v06.006.2015145105749.hdf"]

        cls.lst_files_terra = [x for x in cls.lst_files if "MOD11A2" in x]
        cls.lst_files_aqua = [x for x in cls.lst_files if "MYD11A2" in x]

        cls.lst_duplicate = ["MOD11A2.A2002201.h18v06.006.2014145105749.hdf"]

        cls.referece_metadata = dict(
            RasterXSize=1200,
            RasterYSize=1200,
            geotransform=(0, 1000, 0, 0, 0, -1000),
            projection="EPSG:4326",
            resolution=(1000, -1000),
            nodata=0,
        )

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree("__pycache__")
        except:
            pass

        for file in Path("/tmp").glob("*h5"):
            file.unlink()

    def tearDown(self):
        try:
            shutil.rmtree("/tmp/VIM")
        except:
            pass

        try:
            shutil.rmtree("/tmp/TDA")
        except:
            pass

        try:
            shutil.rmtree("/tmp/TNA")
        except:
            pass

        try:
            shutil.rmtree("/tmp/TDT")
        except:
            pass

        try:
            shutil.rmtree("/tmp/TNT")
        except:
            pass

    def test_raw_instance(self):
        """Test creation of class instance"""

        raw_h5 = ModisRawH5(files=self.vim_files_aqua, targetdir="/tmp")
        self.assertEqual(raw_h5.vam_product_code, "VIM")
        self.assertEqual(raw_h5.product, "MYD13A2")
        self.assertEqual(raw_h5.temporalresolution, 16)
        self.assertEqual(raw_h5.tshift, 8)
        self.assertEqual(str(raw_h5.filename), "/tmp/VIM/MYD13A2.h18v06.006.VIM.h5")

        raw_h5 = ModisRawH5(files=self.vim_files_aqua, targetdir="/tmp", vam_product_code="VEM")
        self.assertEqual(raw_h5.vam_product_code, "VEM")
        self.assertEqual(raw_h5.product, "MYD13A2")
        self.assertEqual(str(raw_h5.filename), "/tmp/VEM/MYD13A2.h18v06.006.VEM.h5")

        with self.assertRaises(AssertionError):
            raw_h5 = ModisRawH5(files=self.vim_files_aqua + self.vim_files_terra, targetdir="/tmp")

        with self.assertRaises(AssertionError):
            raw_h5 = ModisRawH5(files=self.vim_files_aqua, targetdir="/tmp", vam_product_code="LTD")

        raw_h5 = ModisRawH5(files=self.vim_files, targetdir="/tmp", interleave=True)
        self.assertEqual(raw_h5.vam_product_code, "VIM")
        self.assertEqual(raw_h5.product, "MXD13A2")
        self.assertEqual(raw_h5.temporalresolution, 8)
        self.assertEqual(raw_h5.tshift, 8)
        self.assertEqual(str(raw_h5.filename), "/tmp/VIM/MXD13A2.h18v06.006.VIM.h5")
        self.assertNotEqual(raw_h5.nfiles, len(self.vim_files))
        self.assertEqual(raw_h5.nfiles, (len(self.vim_files) - 1))

        raw_h5 = ModisRawH5(files=self.lst_files_aqua, targetdir="/tmp")
        self.assertEqual(raw_h5.vam_product_code, "LTD")
        self.assertEqual(raw_h5.product, "MYD11A2")
        self.assertEqual(raw_h5.temporalresolution, 8)
        self.assertEqual(raw_h5.tshift, 4)
        self.assertEqual(str(raw_h5.filename), "/tmp/TDA/MYD11A2.h18v06.006.TDA.h5")

        raw_h5 = ModisRawH5(files=self.lst_files_aqua, targetdir="/tmp", vam_product_code="LTN")
        self.assertEqual(raw_h5.vam_product_code, "LTN")
        self.assertEqual(str(raw_h5.filename), "/tmp/TNA/MYD11A2.h18v06.006.TNA.h5")

        raw_h5 = ModisRawH5(files=self.lst_files_terra, targetdir="/tmp")
        self.assertEqual(raw_h5.vam_product_code, "LTD")
        self.assertEqual(raw_h5.product, "MOD11A2")
        self.assertEqual(str(raw_h5.filename), "/tmp/TDT/MOD11A2.h18v06.006.TDT.h5")

        with self.assertRaises(AssertionError):
            raw_h5 = ModisRawH5(files=self.lst_files_terra, targetdir="/tmp", vam_product_code="VIM")

        raw_h5 = ModisRawH5(files=self.lst_files_terra + self.lst_duplicate, targetdir="/tmp")
        self.assertEqual(raw_h5.files, sorted(self.lst_files_terra))

    def test_create(self):
        """Test creation of file"""
        h5f = ModisRawH5(
            files=self.vim_files,
            targetdir="/tmp",
            interleave=True,
        )

        self.assertFalse(h5f.exists)

        with patch.object(ModisRawH5, "_get_reference_metadata") as mocked_md:
            mocked_md.return_value = self.referece_metadata

            h5f.create()
            self.assertTrue(h5f.exists)

            with h5py.File(h5f.filename, "r") as hdf5_file:
                ds = hdf5_file.get("data")
                self.assertTrue(ds)

                attrs = ds.attrs
                ysize = attrs["RasterYSize"] * attrs["RasterXSize"]

                self.assertEqual(ds.shape, (ysize, len(self.vim_files) - 1))

                dates = hdf5_file.get("dates")
                self.assertTrue(dates)
                self.assertTrue(dates.shape, (len(h5f.rawdates,)))

            h5f.filename.unlink()

    @patch("modape.modis.collect.HDFHandler.open_datasets")
    @patch("modape.modis.collect.HDFHandler.read_chunk")
    def test_update(self, mocked_chunk, mocked_handles):
        """Test updating dataset"""

        ones = np.ones((48, 1200), dtype="int16")
        ones_view = ones.view()
        ones_view.shape = (48*1200,)

        h5f = ModisRawH5(
            files=self.vim_files,
            targetdir="/tmp",
            interleave=True,
        )

        with patch.object(ModisRawH5, "_get_reference_metadata") as mocked_md:
            mocked_md.return_value = self.referece_metadata
            h5f.create()

        mocked_chunk.return_value = ones

        with patch("modape.modis.collect.HDFHandler.iter_handles") as mocked_iter:
            mocked_iter.return_value = ((ii, None) for ii, x in enumerate(h5f.files))
            h5f.update()

        for test_arr in h5f.read_chunked("data", xchunk=10):
            for dim in range(test_arr.shape[1]):
                np.testing.assert_array_equal(test_arr[:, dim], ones_view)

        with h5py.File(h5f.filename, "r") as hdf5_file:
            dates = [x.decode() for x in hdf5_file.get("dates")]
            self.assertEqual(dates, h5f.rawdates)

        h5f.filename.unlink()
        del h5f

        files = self.lst_files_aqua
        files.sort()

        h5f = ModisRawH5(
            files=files[:2],
            targetdir="/tmp",
        )

        with patch.object(ModisRawH5, "_get_reference_metadata") as mocked_md:
            mocked_md.return_value = self.referece_metadata
            h5f.create()

        with patch("modape.modis.collect.HDFHandler.iter_handles") as mocked_iter:
            mocked_iter.return_value = ((ii, None) for ii, x in enumerate(h5f.files))
            h5f.update()

        dates_init = h5f.rawdates

        del h5f

        h5f = ModisRawH5(
            files=files[2:],
            targetdir="/tmp",
        )

        twos = ones.copy()
        twos[...] = 2
        twos_view = twos.view()
        twos_view.shape = (48*1200,)

        mocked_chunk.return_value = twos

        with patch("modape.modis.collect.HDFHandler.iter_handles") as mocked_iter:
            mocked_iter.return_value = ((ii, None) for ii, x in enumerate(h5f.files))
            h5f.update()

        for test_arr in h5f.read_chunked("data", xchunk=10):
            ii = 0
            for dim in range(test_arr.shape[1]):
                if ii > 1:
                    np.testing.assert_array_equal(test_arr[:, dim], twos_view)
                else:
                    np.testing.assert_array_equal(test_arr[:, dim], ones_view)
                ii += 1

        with h5py.File(h5f.filename, "r") as hdf5_file:
            dates = [x.decode() for x in hdf5_file.get("dates")]
            self.assertEqual(dates, dates_init + h5f.rawdates)


        h5f.filename.unlink()
        del h5f

        h5f = ModisRawH5(
            files=files[2:],
            targetdir="/tmp",
        )

        with patch.object(ModisRawH5, "_get_reference_metadata") as mocked_md:
            mocked_md.return_value = self.referece_metadata
            h5f.create()

        with patch("modape.modis.collect.HDFHandler.iter_handles") as mocked_iter:
            mocked_iter.return_value = ((ii, None) for ii, x in enumerate(h5f.files))
            h5f.update()

        del h5f

        h5f = ModisRawH5(
            files=files[:2],
            targetdir="/tmp",
        )

        with patch("modape.modis.collect.HDFHandler.iter_handles") as mocked_iter:
            mocked_iter.return_value = ((ii, None) for ii, x in enumerate(h5f.files))

            with self.assertRaises(AssertionError):
                h5f.update()

class TestModisSmooth(unittest.TestCase):
    """Test class for ModisSmooth tests"""

    @classmethod
    def setUpClass(cls):
        cls.testpath = Path("/tmp/data")
        cls.testpath.mkdir(exist_ok=True)
        cls.testfile = create_h5temp(12, 12, 8, 8)
        cls.y_chunksize = 12*12//25

    @classmethod
    def tearDownClass(cls):
        cls.testfile.unlink()
        try:
            shutil.rmtree(str(cls.testpath))
        except:
            pass

    def test_smooth_instance(self):
        """Test creation of class instance"""

        smtH5 = ModisSmoothH5(
            rawfile=self.testfile,
            targetdir="/tmp",
        )

        self.assertFalse(smtH5.exists)
        self.assertFalse(smtH5.tinterpolate)
        self.assertEqual(smtH5.temporalresolution, None)
        self.assertEqual(str(smtH5.filename), "/tmp/MXD13A2.h21v10.006.txn.VIM.h5")

        smtH5 = ModisSmoothH5(
            rawfile=self.testfile,
            targetdir="/tmp",
            tempint=10,
        )

        self.assertFalse(smtH5.exists)
        self.assertTrue(smtH5.tinterpolate)
        self.assertEqual(smtH5.temporalresolution, 10)
        self.assertEqual(str(smtH5.filename), "/tmp/MXD13A2.h21v10.006.txd.VIM.h5")


        smtH5 = ModisSmoothH5(
            rawfile=self.testfile,
            targetdir="/tmp",
            tempint=5,
        )

        self.assertFalse(smtH5.exists)
        self.assertTrue(smtH5.tinterpolate)
        self.assertEqual(smtH5.temporalresolution, 5)
        self.assertEqual(str(smtH5.filename), "/tmp/MXD13A2.h21v10.006.txp.VIM.h5")

    def test_create(self):
        """Test creation of smooth HDF5"""

        smtH5 = ModisSmoothH5(
            rawfile=self.testfile,
            targetdir="/tmp/data",
        )

        self.assertFalse(smtH5.exists)
        smtH5.create()
        self.assertTrue(smtH5.exists)

        with h5py.File(smtH5.filename, "r") as hdf5_file:
            ds = hdf5_file.get("data")
            self.assertTrue(ds)

            attrs = ds.attrs
            ysize = attrs["RasterYSize"] * attrs["RasterXSize"]
            self.assertEqual(ds.shape, (ysize, 4))
            self.assertEqual(attrs["temporalresolution"], 8)

            dates = hdf5_file.get("dates")
            self.assertTrue(dates)

        self.assertNotEqual(smtH5.last_smoothed, "2002201")
        self.assertEqual(smtH5.last_smoothed, "2002209")

        smtH5.filename.unlink()

    @patch.object(ModisSmoothH5, "read_chunked")
    @patch.object(ModisSmoothH5, "write_chunk", return_value=False)
    @patch("modape.modis.smooth.HDF5Base.read_chunked")
    def test_smooth(self, mock_hdf5base_read, mocked_write, mocked_read):
        """Test smoothing method"""

        ones = np.ones((self.y_chunksize, 4))

        mock_hdf5base_read.return_value = [ones]

        smtH5 = ModisSmoothH5(
            rawfile=self.testfile,
            targetdir="/tmp/data",
        )

        smtH5.create()
        self.assertTrue(smtH5.exists)

        with self.assertRaises(ValueError):
            smtH5.smooth(nsmooth=1, nupdate=2)

        with self.assertRaises(ValueError):
            smtH5.smooth(soptimize=True, srange=[1, 2, 3])

        ts_test = ones[0, :].copy()

        with patch("modape.modis.smooth.ws2d") as mocked_whit:
            with self.assertRaises(HDF5WriteError):
                mocked_whit.return_value = ts_test
                smtH5.smooth(svalue=1)

        mocked_write.return_value = True

        with patch("modape.modis.smooth.ws2d") as mocked_whit:
            mocked_whit.return_value = ts_test
            smtH5.smooth(svalue=1)

        mocked_whit.assert_called()
        self.assertEqual(mocked_whit.call_count, self.y_chunksize)
        _, mkwargs = mocked_whit.call_args
        np.testing.assert_array_equal(mkwargs["y"], ts_test)
        np.testing.assert_array_equal(mkwargs["w"], ts_test)
        self.assertEqual(mkwargs["lmda"], 10)

        with patch("modape.modis.smooth.ws2dp") as mocked_whit:
            mocked_whit.return_value = ts_test
            smtH5.smooth(svalue=1, p=0.90)

        mocked_whit.assert_called()
        self.assertEqual(mocked_whit.call_count, self.y_chunksize)
        _, mkwargs = mocked_whit.call_args
        np.testing.assert_array_equal(mkwargs["y"], ts_test)
        np.testing.assert_array_equal(mkwargs["w"], ts_test)
        self.assertEqual(mkwargs["lmda"], 10)
        self.assertEqual(mkwargs["p"], 0.90)

        with patch("modape.modis.smooth.ws2d") as mocked_whit:
            mocked_read.return_value = iter([ones[:, 0]])
            mocked_whit.return_value = ts_test
            smtH5.smooth()

        mocked_whit.assert_called()
        self.assertEqual(mocked_whit.call_count, self.y_chunksize)
        _, mkwargs = mocked_whit.call_args
        np.testing.assert_array_equal(mkwargs["y"], ts_test)
        np.testing.assert_array_equal(mkwargs["w"], ts_test)
        self.assertEqual(mkwargs["lmda"], 10)

        with patch("modape.modis.smooth.ws2doptv") as mocked_whit:
            mocked_read.return_value = iter([ones[:, 0]])
            mocked_whit.return_value = (ts_test, 10)
            smtH5.smooth(soptimize=True)

        mocked_whit.assert_called()
        self.assertEqual(mocked_whit.call_count, self.y_chunksize)
        _, mkwargs = mocked_whit.call_args
        np.testing.assert_array_equal(mkwargs["y"], ts_test)
        np.testing.assert_array_equal(mkwargs["w"], ts_test)
        np.testing.assert_array_equal(mkwargs["llas"], np.arange(-1, 1.2, 0.2).round(2))

        with patch("modape.modis.smooth.ws2doptv") as mocked_whit:
            mocked_read.return_value = iter([ones[:, 0]])
            mocked_whit.return_value = (ts_test, 10)
            with patch("modape.modis.smooth.lag1corr", return_value=0.1):
                smtH5.smooth(soptimize=True)

        mocked_whit.assert_called()
        self.assertEqual(mocked_whit.call_count, self.y_chunksize)
        _, mkwargs = mocked_whit.call_args
        np.testing.assert_array_equal(mkwargs["llas"], np.arange(0, 3.2, 0.2).round(2))

        with patch("modape.modis.smooth.ws2doptv") as mocked_whit1, patch("modape.modis.smooth.ws2doptvp") as mocked_whit2:
            mocked_read.return_value = iter([ones[:, 0]])
            mocked_whit2.return_value = (ts_test, 10)
            with patch("modape.modis.smooth.lag1corr", return_value=0.8):
                smtH5.smooth(soptimize=True, p=0.9)

        mocked_whit1.assert_not_called()
        mocked_whit2.assert_called()
        self.assertEqual(mocked_whit2.call_count, self.y_chunksize)
        _, mkwargs = mocked_whit2.call_args
        np.testing.assert_array_equal(mkwargs["llas"], np.arange(-2, 1.2, 0.2).round(2))
        np.testing.assert_array_equal(mkwargs["p"], 0.9)

        smtH5.filename.unlink()

class TestModisMosaic(unittest.TestCase):
    """Test class for ModisMosaic"""

    @classmethod
    def setUpClass(cls):
        cls.testpath = Path("/tmp/data")
        cls.testpath.mkdir(exist_ok=True)
        cls.testfile = create_h5temp(12, 12, 8, 8)
        cls.testfile_global = create_h5temp_global()

        smt_file = ModisSmoothH5(cls.testfile, targetdir="/tmp/data", tempint=5)
        smt_file.create()
        cls.testfile_smt = smt_file.filename


    @classmethod
    def tearDownClass(cls):
        cls.testfile.unlink()
        try:
            shutil.rmtree(str(cls.testpath))
        except:
            pass

    def test_instance(self):
        """Test instance creation"""

        mosaic = ModisMosaic([self.testfile]*5)

        with h5py.File(self.testfile, "r") as h5f_open:
            dts = [x.decode() for x in h5f_open.get("dates")[...]]

        self.assertEqual(
            [x.strftime("%Y%j") for x in mosaic.dates],
            dts
        )
        self.assertTrue(len(mosaic.files), 5)

        mosaic = ModisMosaic(self.testfile)
        self.assertEqual(
            [x.strftime("%Y%j") for x in mosaic.dates],
            dts
        )
        self.assertTrue(len(mosaic.files), 1)

    @patch("modape.modis.window.ModisMosaic._get_raster", return_value="/vsimem/inmem.tif")
    @patch("modape.modis.window.ModisMosaic._mosaic")
    @patch("modape.modis.window.ModisMosaic._translate")
    def test_generate_mosaics_basic(self, mock_translate, mock_mosaic, mock_raster):
        """Test mosaic creation"""
        mosaic = ModisMosaic(self.testfile)

        mock_mosaic.return_value.__enter__.return_value = "/inmem/warped.tif"

        with self.assertRaises(ValueError):
            mosaic.generate_mosaics("not_a_dataset", "/tmp/data", "EPSG:4326")

        mosaic.generate_mosaics("data", "/tmp/data", "EPSG:4326")

        mock_raster.assert_called()
        margs, mkwargs = mock_raster.call_args
        self.assertEqual(mock_raster.call_count, len(mosaic.dates))
        self.assertEqual(margs[0], mosaic.files[0])
        self.assertEqual(margs[1], "data")
        self.assertEqual(margs[2], False)
        self.assertEqual(mkwargs["round_int"], None)
        self.assertEqual(mkwargs["ix"], len(mosaic.dates)-1)

        mock_mosaic.assert_called()
        margs, mkwargs = mock_mosaic.call_args
        self.assertEqual(mock_mosaic.call_count, len(mosaic.dates))
        self.assertEqual(margs[0], ["/vsimem/inmem.tif"])
        self.assertEqual(mkwargs["target_srs"], "EPSG:4326")
        self.assertEqual(mkwargs["dtype"], 3)
        self.assertEqual(mkwargs["nodata"], -3000)
        np.testing.assert_almost_equal(mkwargs["resolution"], [1000/112000, (1000/112000)*-1])

        mock_translate.assert_called()
        _, mkwargs = mock_translate.call_args
        self.assertEqual(mock_translate.call_count, len(mosaic.dates))
        self.assertEqual(mkwargs["src"], "/inmem/warped.tif")
        self.assertEqual(mkwargs["dst"], "/tmp/data/vim2002j209.tif")
        self.assertEqual(mkwargs["outputSRS"], "EPSG:4326")
        self.assertEqual(mkwargs["noData"], -3000)
        self.assertEqual(mkwargs["outputType"], 3)


    @patch("modape.modis.window.ModisMosaic._get_raster", return_value="/vsimem/inmem.tif")
    @patch("modape.modis.window.ModisMosaic._mosaic")
    @patch("modape.modis.window.ModisMosaic._translate")
    def test_generate_mosaics_sgrid(self, mock_translate, mock_mosaic, mock_raster):
        """Test mosaic creation"""
        mosaic = ModisMosaic(self.testfile_smt)

        mock_mosaic.return_value.__enter__.return_value = "/inmem/warped.tif"

        mosaic.generate_mosaics("sgrid", "/tmp/data", "EPSG:4326")

        mock_raster.assert_called()
        margs, mkwargs = mock_raster.call_args
        self.assertEqual(mock_raster.call_count, 1)
        self.assertEqual(margs[0], mosaic.files[0])
        self.assertEqual(margs[1], "sgrid")
        self.assertEqual(margs[2], False)
        self.assertEqual(mkwargs["round_int"], None)
        self.assertEqual(mkwargs["ix"], None)

        mock_mosaic.assert_called()
        margs, mkwargs = mock_mosaic.call_args
        self.assertEqual(mock_mosaic.call_count, 1)
        self.assertEqual(margs[0], ["/vsimem/inmem.tif"])
        self.assertEqual(mkwargs["target_srs"], "EPSG:4326")
        self.assertEqual(mkwargs["dtype"], 6)
        self.assertEqual(mkwargs["nodata"], 0)
        np.testing.assert_almost_equal(mkwargs["resolution"], [1000/112000, (1000/112000)*-1])

        mock_translate.assert_called()
        _, mkwargs = mock_translate.call_args
        self.assertEqual(mock_translate.call_count, 1)
        self.assertEqual(mkwargs["src"], "/inmem/warped.tif")
        self.assertEqual(mkwargs["dst"], "/tmp/data/vim_sgrid.tif")
        self.assertEqual(mkwargs["outputSRS"], "EPSG:4326")
        self.assertEqual(mkwargs["noData"], 0)
        self.assertEqual(mkwargs["outputType"], 6)

    @patch("modape.modis.window.ModisMosaic._get_raster", return_value="/vsimem/inmem.tif")
    @patch("modape.modis.window.ModisMosaic._mosaic")
    @patch("modape.modis.window.ModisMosaic._translate")
    def test_generate_mosaics_kwarg_propagation(self, mock_translate, mock_mosaic, mock_raster):
        """Test mosaic creation"""
        mosaic = ModisMosaic(self.testfile)

        mock_mosaic.return_value.__enter__.return_value = "/inmem/warped.tif"

        cos = ["COMPRESS=LZW", "PREDICTOR=2"]
        aoi = [0, 0, 10, 10]

        mosaic.generate_mosaics("data",
                                "/tmp/data",
                                "EPSG:3857",
                                aoi=aoi,
                                overwrite=True,
                                force_doy=True,
                                prefix="test",
                                clip_valid=True,
                                round_int=-2,
                                xRes=10,
                                yRes=10,
                                noData=-1,
                                outputType=0,
                                creationOptions=cos,
                                resampleAlg="bilinear",
                                multithread=True,
                                )

        mock_raster.assert_called()
        margs, mkwargs = mock_raster.call_args
        self.assertEqual(margs[2], True)
        self.assertEqual(mkwargs["round_int"], -2)

        mock_mosaic.assert_called()
        margs, mkwargs = mock_mosaic.call_args

        self.assertEqual(margs[0], ["/vsimem/inmem.tif"])
        self.assertEqual(mkwargs["target_srs"], "EPSG:3857")
        self.assertEqual(mkwargs["dtype"], 0)
        self.assertEqual(mkwargs["nodata"], -1)
        self.assertEqual(mkwargs["resample"], "bilinear")
        self.assertEqual(mkwargs["gdal_multithread"], True)
        np.testing.assert_almost_equal(mkwargs["resolution"], [10, 10])

        mock_translate.assert_called()
        _, mkwargs = mock_translate.call_args
        self.assertEqual(mkwargs["src"], "/inmem/warped.tif")
        self.assertEqual(mkwargs["dst"], "/tmp/data/testvim2002j209.tif")
        self.assertEqual(mkwargs["outputSRS"], "EPSG:3857")
        self.assertEqual(mkwargs["noData"], -1)
        self.assertEqual(mkwargs["outputType"], 0)
        self.assertEqual(mkwargs["creationOptions"], cos)
        self.assertEqual(mkwargs["resampleAlg"], "bilinear")
        self.assertEqual(mkwargs["projWin"], aoi)
        self.assertEqual(mkwargs["outputBounds"], aoi)

    @patch("modape.modis.window.ModisMosaic._get_raster", return_value="/vsimem/inmem.tif")
    @patch("modape.modis.window.ModisMosaic._mosaic")
    @patch("modape.modis.window.ModisMosaic._translate")
    def test_generate_mosaics_global(self, mock_translate, mock_mosaic, mock_raster):
        """Test mosaic creation"""

        mosaic = ModisMosaic(self.testfile_global)
        mock_mosaic.return_value.__enter__.return_value = "/vsimem/warped.tif"
        mosaic.generate_mosaics("data", "/tmp/data", None)

        mock_raster.assert_called()
        self.assertEqual(mock_raster.call_count, len(mosaic.dates))
        mock_mosaic.assert_called()

        mock_translate.assert_called()
        _, mkwargs = mock_translate.call_args
        self.assertEqual(mock_translate.call_count, len(mosaic.dates))
        self.assertEqual(mkwargs["src"], "/vsimem/warped.tif")
        self.assertEqual(mkwargs["dst"], "/tmp/data/vim2002j209.tif")
        self.assertEqual(mkwargs["outputSRS"], None)

        mock_translate.reset_mock()
        mock_mosaic.reset_mock()

        mosaic.generate_mosaics("data", "/tmp/data", "EPSG:3857")
        mock_mosaic.assert_called()
        self.assertEqual(mock_mosaic.call_count, len(mosaic.dates))
        margs, mkwargs = mock_mosaic.call_args
        self.assertEqual(margs[0], ["/vsimem/inmem.tif"])
        self.assertEqual(mkwargs["target_srs"], "EPSG:3857")
        self.assertEqual(mkwargs["resolution"], [None, None])
        self.assertEqual(mkwargs["gdal_multithread"], False)

        _, mkwargs = mock_translate.call_args
        self.assertEqual(mock_translate.call_count, len(mosaic.dates))
        self.assertEqual(mkwargs["src"], "/vsimem/warped.tif")
        self.assertEqual(mkwargs["dst"], "/tmp/data/vim2002j209.tif")
        self.assertEqual(mkwargs["outputSRS"], "EPSG:3857")

    @patch("modape.modis.window.ModisMosaic._get_raster", return_value="/vsimem/inmem.tif")
    @patch("modape.modis.window.ModisMosaic._mosaic")
    @patch("modape.modis.window.ModisMosaic._translate")
    def test_generate_mosaics_dates(self, mock_translate, mock_mosaic, mock_raster):
        """Test mosaic creation"""
        mosaic = ModisMosaic(self.testfile)

        mock_mosaic.return_value.__enter__.return_value = "/vsimem/warped.tif"

        mosaic.generate_mosaics("data",
                                "/tmp/data",
                                "EPSG:4326",
                                start=datetime(2002, 7, 10).date(),
                                stop=datetime(2002, 7, 21).date())

        n = len(mosaic.dates[1:3])

        mock_raster.assert_called()
        self.assertEqual(mock_raster.call_count, n)
        mock_mosaic.assert_called()
        self.assertEqual(mock_mosaic.call_count, n)
        mock_translate.assert_called()
        self.assertEqual(mock_translate.call_count, n)


if __name__ == "__main__":
    unittest.main()
