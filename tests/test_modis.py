from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os
import re
import shutil
import unittest
from unittest.mock import patch, Mock
import uuid

import h5py
try:
    import gdal
except ImportError:
    from osgeo import gdal

import fake
from modape.modis import MODISquery, MODISrawh5, MODISsmth5

def create_gdal(x, y):
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create('/vsimem/{}.tif'.format(str(uuid.uuid4())), x, y, 1, 3)
    return(ds)

def create_h5(fn, x, y, tr, ts, r):

    with h5py.File(fn, 'a', driver='core', backing_store=True) as h5f:
        dset = h5f.create_dataset('data', shape=(x*y, 4), dtype='Int16', maxshape=(x*y,None), chunks=((x*y)//25, 10), compression='gzip', fillvalue=-3000)
        h5f.create_dataset('dates', data=[x.encode('ascii') for x in ['2002185', '2002193', '2002201', '2002209']] ,shape=(4,), maxshape=(None,), dtype='S8', compression='gzip')
        dset.attrs['nodata'] = -3000
        dset.attrs['temporalresolution'] = tr
        dset.attrs['tshift'] = ts
        dset.attrs['RasterXSize'] = x
        dset.attrs['RasterYSize'] = y
        dset.attrs['geotransform'] = (0, 0, 0, 0, 0)
        dset.attrs['projection'] = 'PRJ'
        dset.attrs['resolution'] = r

class TestMODIS(unittest.TestCase):

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

    @patch('modape.modis.requests.Session')
    def test_query(self, mocked_get):

        class MockRSP:

            def raise_for_status():
                pass

            def get(self,*args):
                rsp = Mock()
                rsp.status_code.return_value = 200
                rsp.raise_for_status.return_value = None

                if 'tiled' in args[0]:
                    rsp.content = fake.tiled

                elif re.fullmatch('http://global-test.query/', args[0]):
                    rsp.content = fake.glob

                elif re.fullmatch('http://global-test.query/\\d.+\\.\\d.+\\.\\d.+/', args[0]):
                    rsp.content = fake.mola[args[0]]
                return(rsp)

        # Test query of MOD tiled NDVI
        mock_rsp = MockRSP()
        mocked_get.return_value.__enter__.return_value = mock_rsp

        urls = ["https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.02.18/MOD13A2.A2000049.h18v06.006.2015136104646.hdf",
                "https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.03.05/MOD13A2.A2000065.h18v06.006.2015136022922.hdf",
                "https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.03.21/MOD13A2.A2000081.h18v06.006.2015136035955.hdf",
                "https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.04.06/MOD13A2.A2000097.h18v06.006.2015136035959.hdf",
                "https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.04.22/MOD13A2.A2000113.h18v06.006.2015137034359.hdf"]

        query = MODISquery(url='http://tiled-test.query', begindate='2000-01-01', enddate='2000-04-30')

        self.assertFalse(query.global_flag)
        self.assertEqual(query.modis_urls, urls)
        self.assertEqual(query.tiles, ['h18v06'])
        self.assertEqual(query.results, 5)
        del query, urls

        urls = ["http://global-test.query/2002.07.04/MYD11C2.A2002185.006.2015168205556.hdf",
                "http://global-test.query/2002.07.12/MYD11C2.A2002193.006.2015149021321.hdf",
                "http://global-test.query/2002.07.20/MYD11C2.A2002201.006.2015149021446.hdf",
                "http://global-test.query/2002.07.28/MYD11C2.A2002209.006.2015149021240.hdf",
                "http://global-test.query/2002.08.05/MYD11C2.A2002217.006.2015149021044.hdf"]

        query = MODISquery(url='http://global-test.query/', begindate='2002-07-01',enddate='2002-08-15', global_flag=True)

        self.assertTrue(query.global_flag)
        self.assertEqual(query.modis_urls, urls)
        self.assertEqual(query.results, 5)

        try:
            shutil.rmtree('__pycache__')
        except:
            pass

    @patch('modape.modis.gdal.Dataset.GetMetadataItem', return_value=-3000)
    @patch('modape.modis.gdal.Dataset.GetSubDatasets', return_value=[['NDVI']])
    @patch('modape.modis.gdal.Open',return_value=create_gdal(1200, 1200))
    def test_rawHDF5(self, mock_ds, mock_sds, mock_nodata):
        # Test raw tiled NDVI with 8-day interleaving of MOD/MYD
        rawfiles = [
        'MOD13A2.A2002193.h18v06.006.*.hdf',
        'MOD13A2.A2002209.h18v06.006.*.hdf',
        'MYD13A2.A2002185.h18v06.006.*.hdf',
        'MYD13A2.A2002201.h18v06.006.*.hdf',
        ]
        rawh5 = MODISrawh5(files=rawfiles, interleave=True)
        mock_ds.assert_called_with('MYD13A2.A2002185.h18v06.006.*.hdf')

        self.assertEqual(rawh5.nfiles, 4)
        self.assertFalse(rawh5.exists)
        self.assertEqual(os.path.basename(rawh5.outname), 'MXD13A2.h18v06.006.VIM.h5')
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
        self.assertEqual(rawh5.chunks,((1200*1200)//25, 10))

        shutil.rmtree(os.path.dirname(rawh5.outname))

        # Test raw global LST DAY
        rawfiles = [
        'MYD11C2.A2002193.*.006.*.hdf',
        'MYD11C2.A2002209.*.006.*.hdf',
        'MYD11C2.A2002185.*.006.*.hdf',
        'MYD11C2.A2002201.*.006.*.hdf',
        ]

        mock_ds.return_value = create_gdal(7200, 3600)
        mock_sds.return_value = [['LST_Day']]

        rawh5 = MODISrawh5(files=rawfiles)
        mock_ds.assert_called_with('MYD11C2.A2002185.*.006.*.hdf')
        self.assertEqual(rawh5.nfiles, 4)
        self.assertFalse(rawh5.exists)
        self.assertEqual(os.path.basename(rawh5.outname), 'MYD11C2.006.TDA.h5')
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
        self.assertEqual(rawh5.chunks,((3600*7200)//25, 10))

        shutil.rmtree(os.path.dirname(rawh5.outname))

    def test_smoothHDF5(self):
        # Test smooth tiled 10-day NDVI
        try:
            create_h5(fn = 'MXD13A2.h18v06.006.VIM.h5', x=1200, y=1200, tr=8, ts=8, r=0.009)
            smth5 = MODISsmth5('MXD13A2.h18v06.006.VIM.h5', tempint=10)

            self.assertEqual(os.path.basename(smth5.outname), 'MXD13A2.h18v06.006.txd.VIM.h5')
            self.assertEqual(smth5.rawdates, [
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
            create_h5(fn = 'MOD11C2.006.LTD.h5', x=3600, y=7200, tr=8, ts=4, r=0.05)
            smth5 = MODISsmth5('MOD11C2.006.LTD.h5', tempint=5)

            self.assertEqual(os.path.basename(smth5.outname), 'MOD11C2.006.txp.LTD.h5')
            self.assertEqual(smth5.rawdates, [
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

if __name__ == '__main__':
    unittest.main()
