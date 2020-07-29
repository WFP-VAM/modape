"""test_modis.py: Test MODIS classes and functions."""
# pylint: disable=E0401,E0611,W0702,W0613
from datetime import datetime
import os
from pathlib import Path
import pickle
import shutil
import unittest
from unittest.mock import patch
import uuid

import numpy as np
import h5py #pylint: disable=import-error
try:
    import gdal
except ImportError:
    from osgeo import gdal

from modape.exceptions import DownloadError
from modape.modis import ModisQuery, ModisRawH5, ModisSmoothH5, modis_tiles
from modape.utils import SessionWithHeaderRedirection

def create_gdal(x, y):
    """Create in-memory gdal dataset for testing.

    Returns:
        GDAL dataset object
    """
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create('/vsimem/{}.tif'.format(str(uuid.uuid4())), x, y, 1, 3)
    return ds

def create_h5(fn, x, y, tr, ts, r):
    '''Create HDF5 file for testing.'''

    with h5py.File(fn, 'a', driver='core', backing_store=True) as h5f:
        dset = h5f.create_dataset('data', shape=(x*y, 4), dtype='Int16', maxshape=(x*y, None), chunks=((x*y)//25, 10), compression='gzip', fillvalue=-3000)
        h5f.create_dataset('dates', shape=(4,), data=np.array(['2002185', '2002193', '2002201', '2002209'], dtype='S8'), maxshape=(None,), dtype='S8', compression='gzip')
        dset.attrs['nodata'] = -3000
        dset.attrs['temporalresolution'] = tr
        dset.attrs['tshift'] = ts
        dset.attrs['RasterXSize'] = x
        dset.attrs['RasterYSize'] = y
        dset.attrs['geotransform'] = (0, 0, 0, 0, 0)
        dset.attrs['projection'] = 'PRJ'
        dset.attrs['resolution'] = r


class TestModisQuery(unittest.TestCase):
    """Test class for ModisQuery tests."""

    @classmethod
    def setUpClass(cls):
        '''Set up testing class'''

        data_dir = str(Path(__file__).parent).replace('tests', 'modape')

        with open(f"{data_dir}/data/cmr_api_response.pkl", 'rb') as pkl:
            cls.api_response = pickle.load(pkl)

        cls.testpath = Path(__name__).parent

    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree('__pycache__')
        except:
            pass

    def setUp(self):
        '''Set up test'''

        # create query

        self.query = ModisQuery(
            products=['MOD13A2', 'MYD13A2'],
            aoi=(10, 10, 20, 20),
            begindate=datetime(2020, 1, 1),
            enddate=datetime(2020, 7, 24),
        )

    def test_query(self):
        '''Test API query'''

        self.assertEqual(self.query.api.params['short_name'], ['MOD13A2', 'MYD13A2'])
        self.assertEqual(self.query.api.params['bounding_box'], '10.0,10.0,20.0,20.0')
        self.assertEqual(self.query.api.params['temporal'], ['2020-01-01T00:00:00Z,2020-07-24T00:00:00Z'])

    def test_response_parse(self):
        '''Test parsing of response'''

        with patch('modape.modis.download.GranuleQuery.get_all',
                   return_value=self.api_response):

            self.query.search(strict_dates=False)

            self.assertTrue(self.query.results)
            self.assertEqual(len(self.query.results), len(self.api_response))

            tiles = list({x['tile'] for x in self.query.results})
            tiles_select = tiles[:2]
            self.assertEqual(len(tiles), 6)

            self.setUp()
            self.query.tile_filter = tiles_select
            self.query.search(strict_dates=True)

            self.assertEqual(len({x['tile'] for x in self.query.results}), 2)
            self.assertTrue(all([x['time_start'] >= self.query.begin.date() for x in self.query.results]))
            self.assertTrue(all([x['time_end'] <= self.query.end.date() for x in self.query.results]))

    @patch("modape.modis.download.ThreadPoolExecutor")#, new_callable=fake_fun)
    @patch("modape.modis.download.GranuleQuery.get_all")
    def test_download(self, mock_response, mock_submit):
        '''Test download'''

        mock_response.return_value = self.api_response
        self.query.search()

        self.query.results = [self.query.results[0]]

        # Successful downloads

        future_result = (self.testpath.joinpath(self.query.results[0]['file_id']), None)
        mock_submit.return_value.__enter__.return_value.submit.return_value.result.return_value = future_result

        try:
            future_result[0].unlink()
        except FileNotFoundError:
            pass

        # Raise AssertionError when file is missing
        with self.assertRaises(AssertionError):
            self.query.download(
                targetdir=self.testpath,
                username='test',
                password='test',
                multithread=True,
            )

        future_result[0].touch()

        self.query.download(
            targetdir=self.testpath,
            username='test',
            password='test',
            multithread=True,
        )

        with patch("modape.modis.download.ModisQuery._fetch_hdf", return_value=future_result) as mocked_fetch:
            self.query.download(
                targetdir=(self.testpath),
                username='test',
                password='test',
                multithread=False,
            )

            mocked_fetch.assert_called_once()
            fetch_args = mocked_fetch.call_args[0]

            self.assertEqual(type(fetch_args[0]), SessionWithHeaderRedirection)
            self.assertEqual(fetch_args[1], self.query.results[0]['link'])
            self.assertEqual(fetch_args[2], self.testpath)

        future_result[0].unlink()

        # Download Error
        future_result = (f"http://datalocation.com/{self.query.results[0]['file_id']}", "Error")
        mock_submit.return_value.__enter__.return_value.submit.return_value.result.return_value = future_result

        with self.assertRaises(DownloadError):
            self.query.download(
                targetdir=self.testpath,
                username='test',
                password='test',
                multithread=True
            )

class TestModisCollect(unittest.TestCase):
    """Test class for ModisQuery tests."""

    @classmethod
    def setUpClass(cls):
        '''Set up testing class'''


        cls.vim_files = ['MYD13A2.A2002201.h18v06.006.2015149071105.hdf',
                         'MYD13A2.A2002185.h18v06.006.2015149071113.hdf',
                         'MOD13A2.A2002177.h18v06.006.2015149001129.hdf',
                         'MOD13A2.A2002209.h18v06.006.2015149180726.hdf',
                         'MOD13A2.A2002193.h18v06.006.2015149022847.hdf']

        cls.vim_files_terra = [x for x in cls.vim_files if 'MOD13A2' in x]
        cls.vim_files_aqua = [x for x in cls.vim_files if 'MYD13A2' in x]

        cls.lst_files = ['MYD11A2.A2002193.h18v06.006.2015146152945.hdf',
                         'MOD11A2.A2002209.h18v06.006.2015145152020.hdf',
                         'MYD11A2.A2002201.h18v06.006.2015146153241.hdf',
                         'MYD11A2.A2002185.h18v06.006.2015146152642.hdf',
                         'MOD11A2.A2002177.h18v06.006.2015144183717.hdf',
                         'MYD11A2.A2002209.h18v06.006.2015152152813.hdf',
                         'MOD11A2.A2002185.h18v06.006.2015145002847.hdf',
                         'MOD11A2.A2002193.h18v06.006.2015145055806.hdf',
                         'MOD11A2.A2002201.h18v06.006.2015145105749.hdf']

        cls.lst_files_terra = [x for x in cls.lst_files if 'MOD11A2' in x]
        cls.lst_files_aqua = [x for x in cls.lst_files if 'MYD11A2' in x]

        cls.lst_duplicate = ['MOD11A2.A2002201.h18v06.006.2014145105749.hdf']


    @classmethod
    def tearDownClass(cls):
        try:
            shutil.rmtree('__pycache__')
        except:
            pass

    def test_raw_instance(self):
        """Test creation of class instance"""

        raw_h5 = ModisRawH5(files=self.vim_files_aqua, targetdir='/temp')
        self.assertEqual(raw_h5.vam_product_code, "VIM")
        self.assertEqual(raw_h5.product, "MYD13A2")
        self.assertEqual(raw_h5.temporalresolution, 16)
        self.assertEqual(raw_h5.tshift, 8)
        self.assertEqual(str(raw_h5.filename), "/tmp/VIM/MYD13A2.h18v06.006.VIM.h5")

        raw_h5 = ModisRawH5(files=self.vim_files_aqua, targetdir='/temp', vam_product_code="VEM")
        self.assertEqual(raw_h5.vam_product_code, "VEM")
        self.assertEqual(raw_h5.product, "MYD13A2")
        self.assertEqual(str(raw_h5.filename), "/tmp/VEM/MYD13A2.h18v06.006.VEM.h5")

        with self.assertRaises(AssertionError):
            raw_h5 = ModisRawH5(files=self.vim_files_aqua, targetdir='/temp', interleave=True)

        with self.assertRaises(AssertionError):
            raw_h5 = ModisRawH5(files=self.vim_files_aqua + self.vim_files_terra, targetdir='/temp')

        with self.assertRaises(AssertionError):
            raw_h5 = ModisRawH5(files=self.vim_files_aqua, targetdir='/temp', vam_product_code="LTD")

        raw_h5 = ModisRawH5(files=self.vim_files, targetdir='/temp', interleave=True)
        self.assertEqual(raw_h5.vam_product_code, "VIM")
        self.assertEqual(raw_h5.product, "MXD13A2")
        self.assertEqual(raw_h5.temporalresolution, 8)
        self.assertEqual(raw_h5.tshift, 8)
        self.assertEqual(str(raw_h5.filename), "/tmp/VIM/MXD13A2.h18v06.006.VIM.h5")
        self.assertNotEqual(raw_h5.nfiles, len(self.vim_files))
        self.assertEqual(raw_h5.nfiles, (len(self.vim_files) - 1))

        raw_h5 = ModisRawH5(files=self.lst_files_aqua, targetdir='/temp')
        self.assertEqual(raw_h5.vam_product_code, "LTD")
        self.assertEqual(raw_h5.product, "MYD11A2")
        self.assertEqual(raw_h5.temporalresolution, 8)
        self.assertEqual(raw_h5.tshift, 4)
        self.assertEqual(str(raw_h5.filename), "/tmp/TDA/MYD11A2.h18v06.006.TDA.h5")

        raw_h5 = ModisRawH5(files=self.lst_files_aqua, targetdir='/temp', vam_product_code="LTN")
        self.assertEqual(raw_h5.vam_product_code, "LTN")
        self.assertEqual(str(raw_h5.filename), "/tmp/TDN/MYD11A2.h18v06.006.TNA.h5")

        raw_h5 = ModisRawH5(files=self.lst_files_terra, targetdir='/temp')
        self.assertEqual(raw_h5.vam_product_code, "LTD")
        self.assertEqual(raw_h5.product, "MOD11A2")
        self.assertEqual(str(raw_h5.filename), "/tmp/TDN/MYD11A2.h18v06.006.TDT.h5")

        with self.assertRaises(AssertionError):
            raw_h5 = ModisRawH5(files=self.lst_files_terra, targetdir='/temp', interleave=True)

        with self.assertRaises(AssertionError):
            raw_h5 = ModisRawH5(files=self.lst_files_terra, targetdir='/temp', vam_product_code="VIM")

        raw_h5 = ModisRawH5(files=self.lst_files + self.lst_duplicate, targetdir='/temp')
        self.assertEqual(raw_h5.files, sorted(self.lst_files))

class TestMODIS(unittest.TestCase):
    """Test class for MODIS tests."""

    @classmethod
    def setUpClass(cls):
        pass

    @classmethod
    def tearDownClass(cls):
        pass

    def setUp(self):
        pass

    def tearDown(self):
        pass


    @patch('modape.modis.collect.gdal.Dataset.GetMetadataItem', return_value=-3000)
    @patch('modape.modis.collect.gdal.Dataset.GetSubDatasets', return_value=[['NDVI']])
    @patch('modape.modis.collect.gdal.Open', return_value=create_gdal(1200, 1200))
    def test_raw_hdf5(self, mock_ds, mock_sds, mock_nodata):
        """Test raw tiled NDVI with 8-day interleaving of MOD/MYD and raw global LST DAY."""
        rawfiles = [
            'MOD13A2.A2002193.h18v06.006.2019256103823.hdf',
            'MOD13A2.A2002209.h18v06.006.2019256103823.hdf',
            'MYD13A2.A2002185.h18v06.006.2019256103823.hdf',
            'MYD13A2.A2002201.h18v06.006.2019256103823.hdf',
        ]
        rawh5 = ModisRawH5(files=rawfiles, interleave=True)
        mock_ds.assert_called_with('MYD13A2.A2002185.h18v06.006.2019256103823.hdf')

        self.assertEqual(rawh5.nfiles, 4)
        self.assertFalse(rawh5.exists)
        self.assertEqual(rawh5.outname.name, 'MXD13A2.h18v06.006.VIM.h5')
        self.assertEqual(rawh5.temporalresolution, 8)
        self.assertEqual(rawh5.tshift, 8)
        self.assertEqual(rawh5.rawdates, [
            '2002185',
            '2002193',
            '2002201',
            '2002209',
        ])

        rawh5.create()
        self.assertTrue(rawh5.exists)
        self.assertEqual(rawh5.nodata_value, -3000)
        self.assertEqual(rawh5.chunks, ((1200*1200)//25, 10))

        shutil.rmtree(rawh5.outname.parent.name)

        # Test handling of duplicate files
        rawfiles = [
            'MOD13A2.A2002193.h18v06.006.2019256103823.hdf',
            'MOD13A2.A2002209.h18v06.006.2019256103823.hdf',
            'MOD13A2.A2002209.h18v06.006.2018256103823.hdf',
            'MYD13A2.A2002185.h18v06.006.2019256103823.hdf',
            'MYD13A2.A2002185.h18v06.006.2018256103823.hdf',
            'MYD13A2.A2002201.h18v06.006.2019256103823.hdf',
        ]
        rawh5 = ModisRawH5(files=rawfiles, interleave=True)
        mock_ds.assert_called_with('MYD13A2.A2002185.h18v06.006.2019256103823.hdf')

        self.assertEqual(rawh5.nfiles, 4)
        self.assertEqual(rawh5.temporalresolution, 8)
        self.assertEqual(rawh5.tshift, 8)
        self.assertEqual(rawh5.rawdates, [
            '2002185',
            '2002193',
            '2002201',
            '2002209',
        ])

        # Test raw global LST DAY
        rawfiles = [
            'MYD11C2.A2002193.*.006.2019256103823.hdf',
            'MYD11C2.A2002209.*.006.2019256103823.hdf',
            'MYD11C2.A2002185.*.006.2019256103823.hdf',
            'MYD11C2.A2002201.*.006.2019256103823.hdf',
        ]

        mock_ds.return_value = create_gdal(7200, 3600)
        mock_sds.return_value = [['LST_Day']]

        rawh5 = ModisRawH5(files=rawfiles)
        mock_ds.assert_called_with('MYD11C2.A2002185.*.006.2019256103823.hdf')
        self.assertEqual(rawh5.nfiles, 4)
        self.assertFalse(rawh5.exists)
        self.assertEqual(rawh5.outname.name, 'MYD11C2.006.TDA.h5')
        self.assertEqual(rawh5.temporalresolution, 8)
        self.assertEqual(rawh5.tshift, 4)
        self.assertEqual(rawh5.rawdates, [
            '2002185',
            '2002193',
            '2002201',
            '2002209',
        ])

        rawh5.create()
        self.assertTrue(rawh5.exists)
        self.assertEqual(rawh5.nodata_value, -3000)
        self.assertEqual(rawh5.chunks, ((3600*7200)//25, 10))

        shutil.rmtree(rawh5.outname.parent.name)

        # Test handling of duplicate files
        rawfiles = [
            'MYD11C2.A2002193.*.006.2019256103823.hdf',
            'MYD11C2.A2002209.*.006.2019256103823.hdf',
            'MYD11C2.A2002209.*.006.2018256103823.hdf',
            'MYD11C2.A2002185.*.006.2019256103823.hdf',
            'MYD11C2.A2002201.*.006.2019256103823.hdf',
            'MYD11C2.A2002201.*.006.2018256103823.hdf',
        ]

        rawh5 = ModisRawH5(files=rawfiles)
        mock_ds.assert_called_with('MYD11C2.A2002185.*.006.2019256103823.hdf')
        self.assertEqual(rawh5.nfiles, 4)
        self.assertEqual(rawh5.outname.name, 'MYD11C2.006.TDA.h5')
        self.assertEqual(rawh5.temporalresolution, 8)
        self.assertEqual(rawh5.tshift, 4)
        self.assertEqual(rawh5.rawdates, [
            '2002185',
            '2002193',
            '2002201',
            '2002209',
        ])

    def test_smoothHDF5(self):
        """Test smooth tiled 10-day NDVI and global 5-day LST Day."""
        try:
            create_h5(fn='MXD13A2.h18v06.006.VIM.h5', x=1200, y=1200, tr=8, ts=8, r=0.009)
            smth5 = ModisSmoothH5('MXD13A2.h18v06.006.VIM.h5', tempint=10)

            self.assertEqual(smth5.outname.name, 'MXD13A2.h18v06.006.txd.VIM.h5')
            self.assertEqual(smth5.rawdates_nsmooth, [
                '2002185',
                '2002193',
                '2002201',
                '2002209',
            ])
            self.assertTrue(smth5.tinterpolate)
            self.assertEqual(smth5.temporalresolution, 10)
            self.assertFalse(smth5.exists)

            smth5.create()
            self.assertTrue(smth5.exists)

            with h5py.File('MXD13A2.h18v06.006.txd.VIM.h5', 'r+') as h5f:
                self.assertEqual([x.decode() for x in h5f.get('dates')[...]],
                                 ['2002186', '2002196', '2002206', '2002217'])
        except:
            try:
                os.remove('MXD13A2.h18v06.006.VIM.h5')
                os.remove('MXD13A2.h18v06.006.txd.VIM.h5')
            except:
                pass
            raise
        else:
            os.remove('MXD13A2.h18v06.006.VIM.h5')
            os.remove('MXD13A2.h18v06.006.txd.VIM.h5')

        # Test smooth global 5-day LST Day
        try:
            create_h5(fn='MOD11C2.006.LTD.h5', x=3600, y=7200, tr=8, ts=4, r=0.05)
            smth5 = ModisSmoothH5('MOD11C2.006.LTD.h5', tempint=5)

            self.assertEqual(smth5.outname.name, 'MOD11C2.006.txp.LTD.h5')
            self.assertEqual(smth5.rawdates_nsmooth, [
                '2002185',
                '2002193',
                '2002201',
                '2002209',
            ])

            self.assertTrue(smth5.tinterpolate)
            self.assertEqual(smth5.temporalresolution, 5)
            self.assertFalse(smth5.exists)

            smth5.create()
            self.assertTrue(smth5.exists)

            with h5py.File('MOD11C2.006.txp.LTD.h5', 'r+') as h5f:
                self.assertEqual([x.decode() for x in h5f.get('dates')[...]],
                                 ['2002189', '2002194', '2002199', '2002204', '2002209', '2002215'])
        except:
            try:
                os.remove('MOD11C2.006.LTD.h5')
                os.remove('MOD11C2.006.txp.LTD.h5')
            except:
                pass
            raise
        else:
            os.remove('MOD11C2.006.LTD.h5')
            os.remove('MOD11C2.006.txp.LTD.h5')

    def test_modis_tiles(self):
        """Test modis_tiles."""
        tiles = modis_tiles([12, 19, 29, 1]) # xmin, ymax, xmax, ymin
        self.assertEqual(tiles, ['h19v07', 'h19v08', 'h20v07', 'h20v08'])

if __name__ == '__main__':
    unittest.main()
