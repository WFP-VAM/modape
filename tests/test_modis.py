import unittest
from unittest.mock import patch
import uuid
import os
import shutil
import gdal
import h5py

from wsmtk.modis import MODISrawh5, MODISsmth5

def create_gdal(x,y):

    driver = gdal.GetDriverByName('GTiff')

    ds = driver.Create('/vsimem/{}.tif'.format(str(uuid.uuid4())),x,y,1,3)

    return(ds)

def create_h5():

    with h5py.File('MXD13A2.h18v06.006.VIM.h5','a',driver='core',backing_store=True) as h5f:

        dset = h5f.create_dataset('data',shape=(1200*1200,4),dtype='Int16',maxshape=(1200*1200,None),chunks=((1200*1200)//25,10),compression='gzip',fillvalue=-3000)
        h5f.create_dataset('dates',data = [x.encode('ascii') for x in ['2002185','2002193','2002201','2002209']],shape=(4,),maxshape=(None,),dtype='S8',compression='gzip')
        dset.attrs['nodata'] = -3000
        dset.attrs['temporalresolution'] = 8
        dset.attrs['tshift'] = 8
        dset.attrs['RasterXSize'] = 1200
        dset.attrs['RasterYSize'] = 1200
        dset.attrs['geotransform'] = (0,0,0,0,0)
        dset.attrs['projection'] = 'PRJ'
        dset.attrs['resolution'] = 0.009


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


    @patch('wsmtk.modis.gdal.Dataset.GetMetadataItem',return_value = -3000)
    @patch('wsmtk.modis.gdal.Dataset.GetSubDatasets',return_value = [['NDVI']])
    @patch('wsmtk.modis.gdal.Open',return_value = create_gdal(1200,1200))
    def test_rawHDF5(self,mock_ds,mock_sds,mock_nodata):

        # Test raw tiled NDVI with 8-day interleaving of MOD/MYD

        rawfiles = [
        'MOD13A2.A2002193.h18v06.006.*.hdf',
        'MOD13A2.A2002209.h18v06.006.*.hdf',
        'MYD13A2.A2002185.h18v06.006.*.hdf',
        'MYD13A2.A2002201.h18v06.006.*.hdf',
        ]


        rawh5 = MODISrawh5(files = rawfiles)
        mock_ds.assert_called_with('MYD13A2.A2002185.h18v06.006.*.hdf')

        self.assertEqual(rawh5.nfiles,4)
        self.assertFalse(rawh5.exists)
        self.assertEqual(os.path.basename(rawh5.outname),'MXD13A2.h18v06.006.VIM.h5')
        self.assertEqual(rawh5.temporalresolution,8)
        self.assertEqual(rawh5.tshift,8)
        self.assertEqual(rawh5.rawdates,['2002185','2002193','2002201','2002209'])

        rawh5.create()

        self.assertTrue(rawh5.exists)
        self.assertEqual(rawh5.nodata_value,-3000)
        self.assertEqual(rawh5.chunks,((1200*1200)//25,10))

        shutil.rmtree(os.path.dirname(rawh5.outname))


        # Test raw global LST DAY

        rawfiles = [
        'MYD11C2.A2002193.*.006.*.hdf',
        'MYD11C2.A2002209.*.006.*.hdf',
        'MYD11C2.A2002185.*.006.*.hdf',
        'MYD11C2.A2002201.*.006.*.hdf',
        ]


        mock_ds.return_value = create_gdal(7200,3600)
        mock_sds.return_value = [['LST_Day']]

        rawh5 = MODISrawh5(files = rawfiles)

        mock_ds.assert_called_with('MYD11C2.A2002185.*.006.*.hdf')

        self.assertEqual(rawh5.nfiles,4)
        self.assertFalse(rawh5.exists)
        self.assertEqual(os.path.basename(rawh5.outname),'MYD11C2.006.LTD.h5')
        self.assertEqual(rawh5.temporalresolution,8)
        self.assertEqual(rawh5.tshift,4)
        self.assertEqual(rawh5.rawdates,['2002185','2002193','2002201','2002209'])

        rawh5.create()

        self.assertTrue(rawh5.exists)
        self.assertEqual(rawh5.nodata_value,-3000)
        self.assertEqual(rawh5.chunks,((3600*7200)//25,10))

        shutil.rmtree(os.path.dirname(rawh5.outname))


    #@patch('wsmtk.modis.h5py.File')
    def test_smoothHDF5(self):

        create_h5()

        smth5 = MODISsmth5('MXD13A2.h18v06.006.VIM.h5',tempint = 10)

        self.assertEqual(os.path.basename(smth5.outname),'MXD13A2.h18v06.006.txd.VIM.h5')
        self.assertEqual(smth5.rawdates,['2002185','2002193','2002201','2002209'])
        self.assertTrue(smth5.tinterpolate)
        self.assertEqual(smth5.temporalresolution,10)
        self.assertFalse(smth5.exists)

        smth5.create()

        self.assertTrue(smth5.exists)

        with h5py.File('MXD13A2.h18v06.006.txd.VIM.h5','r+') as h5f:

            self.assertEqual([x.decode() for x in h5f.get('dates')[...]],['2002186', '2002196', '2002206', '2002217'])


        os.remove('MXD13A2.h18v06.006.VIM.h5')
        os.remove('MXD13A2.h18v06.006.txd.VIM.h5')


    ## TO TEST: lst, parameter,  global



    # def test_query(self):
    #
    #     with patch('wsmtk.modis.requests.session') as mocked_get:
    #
    #         mocked_get.get = 'test'
    #         mocked_get.statuscode = 200
    #         mocked_get.return_value.content = b'''<?xml version="1.0"?>
    #         <inventory>
    #             <url>https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.02.18/MOD13A2.A2000049.h18v06.006.2015136104646.hdf</url>
    #             <url>https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.03.05/MOD13A2.A2000065.h18v06.006.2015136022922.hdf</url>
    #             <url>https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.03.21/MOD13A2.A2000081.h18v06.006.2015136035955.hdf</url>
    #             <url>https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.04.06/MOD13A2.A2000097.h18v06.006.2015136035959.hdf</url>
    #             <url>https://e4ftl01.cr.usgs.gov//MODV6_Cmp_B/MOLT/MOD13A2.006/2000.04.22/MOD13A2.A2000113.h18v06.006.2015137034359.hdf</url>
    #         </inventory>
    #         '''
    #
    #         query = MODISquery(url='mock://test.query',begindate = '2002-07-04',enddate = '2003-12-31')
    #         print(query)





if __name__ == '__main__':
    unittest.main()
