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


from exceptions import DownloadError
from modape.modis import ModisQuery, ModisRawH5, ModisSmoothH5, modis_tiles
from utils import SessionWithHeaderRedirection

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

        # load response data
        with open('cmr_api_response.pkl', 'rb') as f:
            cls.api_response = pickle.load(f)

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

        with patch('modis.download.GranuleQuery.get_all',
                   return_value=self.api_response):

            self.query.query(strict_dates=False)

            self.assertTrue(self.query.results)
            self.assertEqual(len(self.query.results), len(self.api_response))

            tiles = list({x['tile'] for x in self.query.results})
            tiles_select = tiles[:2]
            self.assertEqual(len(tiles), 6)

            self.setUp()
            self.query.tile_filter = tiles_select
            self.query.query(strict_dates=True)

            self.assertEqual(len({x['tile'] for x in self.query.results}), 2)
            self.assertTrue(all([x['time_start'] >= self.query.begin.date() for x in self.query.results]))
            self.assertTrue(all([x['time_end'] <= self.query.end.date() for x in self.query.results]))

    @patch("modis.download.ThreadPoolExecutor")#, new_callable=fake_fun)
    @patch("modis.download.GranuleQuery.get_all")
    def test_download(self, mock_response, mock_submit):
        '''Test download'''

        mock_response.return_value = self.api_response
        self.query.query()

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

        with patch("modis.download.ModisQuery._fetch_hdf", return_value=future_result) as mocked_fetch:
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
