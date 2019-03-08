import unittest
import os
import pickle
import numpy as np
import array
from wsmtk.whittaker import *

class TestWhittaker(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        with open('{}/data/MXD_testdata.pkl'.format(os.path.dirname(__file__).replace('tests','wsmtk')),'rb') as pkl:
            cls.data = pickle.load(pkl)

    @classmethod
    def tearDownClass(cls):
        cls.data = None

    def setUp(self):

        self.y = self.data['y']
        self.w = self.data['w']

    def tearDown(self):

        self.y = None
        self.w = None

    def test_lag1corr(self):

        self.assertAlmostEqual(lag1corr(self.y[:-1],self.y[1:],-3000.0),self.data['lag1corr'])

    def test_ws2d(self):

        z = np.array(ws2d(self.y,10,self.w),dtype='double')

        self.assertTrue(np.all(z == self.data['z_ws2d']))


    def test_ws2dvc(self):

        z,sopt = ws2d_vc(self.y,self.w,array.array('d',np.linspace(-2,1,16)))

        self.assertTrue(np.all(np.array(z,dtype='double') == self.data['z_ws2dvc']))
        self.assertEqual(sopt,self.data['sopt_ws2dvc'])

    def test_ws2dvcp(self):

        z,sopt = ws2d_vc_asy(self.y,self.w,array.array('d',np.linspace(-2,1,16)),p = 0.90)

        self.assertTrue(np.all(np.array(z,dtype='double') == self.data['z_ws2dvcp']))
        self.assertEqual(sopt,self.data['sopt_ws2dvcp'])


if __name__ == '__main__':
    unittest.main()
