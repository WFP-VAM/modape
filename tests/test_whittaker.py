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

    def test_lag1corr(self):
        y = self.data['y']
        self.assertEqual(lag1corr(y[:-1],y[1:],-3000.0),self.data['lag1corr'])

    def test_ws2d(self):
        y = self.data['y']
        w = self.data['w']

        z = np.array(ws2d(y,10,w),dtype='double')

        self.assertTrue(np.all(z == self.data['z_ws2d']))


    def test_ws2dvc(self):
        y = self.data['y']
        w = self.data['w']

        z,sopt = ws2d_vc(y,w,array.array('d',np.linspace(-2.0,1.0,16.0)))

        self.assertTrue(np.all(np.array(z,dtype='double') == self.data['z_ws2dvc']))
        self.assertEqual(sopt,self.data['sopt_ws2dvc'])

    def test_ws2dvcp(self):
        y = self.data['y']
        w = self.data['w']

        z,sopt = ws2d_vc_asy(y,w,array.array('d',np.linspace(-2.0,1.0,16.0)),p = 0.90)

        self.assertTrue(np.all(np.array(z,dtype='double') == self.data['z_ws2dvcp']))
        self.assertEqual(sopt,self.data['sopt_ws2dvcp'])


if __name__ == '__main__':
    unittest.main()
