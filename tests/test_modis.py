import unittest
from unittest.mock import patch
import uuid
import os
import shutil
import gdal
import h5py

from wsmtk.modis import MODISrawh5, MODISsmth5

def create_gdal():

    driver = gdal.GetDriverByName('GTiff')

    ds = driver.Create('/vsimem/{}.tif'.format(str(uuid.uuid4())),1200,1200,1,3)

    return(ds)

def create_h5():

    h5f = h5py.File('MXD13A2.h18v06.006.VIM.h5','x',driver='core',backing_store=False)


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

    return(h5f.fid)


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
    @patch('wsmtk.modis.gdal.Open',return_value = create_gdal())
    def test_rawHDF5(self,mock_ds,mock_sds,mock_nodata):

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
        self.assertEqual(rawh5.tshift,8)
        self.assertEqual(rawh5.rawdates,['2002185','2002193','2002201','2002209'])

        rawh5.create()

        self.assertTrue(rawh5.exists)
        self.assertEqual(rawh5.nodata_value,-3000)
        self.assertEqual(rawh5.chunks,((1200*1200)//25,10))

        shutil.rmtree(os.path.dirname(rawh5.outname))


if __name__ == '__main__':
    unittest.main()
