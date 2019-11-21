"""test_whittaker.py: Test core whittaker functions."""
# pylint: disable=invalid-name
from __future__ import absolute_import, division, print_function

from array import array
from datetime import date
import os
import pickle
import unittest

import numpy as np
import pandas as pd #pylint: disable=E0401

from modape.utils import fromjulian
from modape.whittaker import lag1corr, w_constrain, ws2d, ws2dp, ws2doptv, ws2doptvp # pylint: disable=E0611

class TestWhittaker(unittest.TestCase):
    """Test class for core whittaker functions."""

    @classmethod
    def setUpClass(cls):
        with open('{}/data/MXD_testdata.pkl'.format(os.path.dirname(__file__).replace('tests', 'modape')), 'rb') as pkl:
            cls.data = pickle.load(pkl)
        with open('{}/data/MXD_dates.pkl'.format(os.path.dirname(__file__).replace('tests', 'modape')), 'rb') as pkl:
            cls.dates = [fromjulian(x) for x in pickle.load(pkl)]

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
        """Test lag-1 correlation function."""
        self.assertAlmostEqual(lag1corr(self.y[:-1], self.y[1:], -3000.0), self.data['lag1corr'])

    def test_ws2d(self):
        """Test ws2d smoothing."""
        z = np.array(ws2d(self.y, 10, self.w), dtype='double')
        np.testing.assert_almost_equal(z, self.data['z_ws2d'], 5)

    def test_ws2dp(self):
        """Test ws2dp smoothing."""
        z = ws2dp(self.y, self.data['sopt_ws2dvcp'], self.w, p=0.90)
        np.testing.assert_almost_equal(z, self.data['z_ws2dp'], 5)

    def test_ws2dp_constrained_update(self):
        """Test constrained update with ws2dp."""

        z_rock = self.data['z_ws2dvcp'][[x < date(2014, 1, 1) for x in self.dates]]*10000
        dates_rock = [x for x in self.dates if x < date(2014, 1, 1)]

        df = pd.DataFrame({'z': np.diff((z_rock)), 'md': [int(x.strftime('%m%d')) for x in dates_rock[1:]]})
        mn = df.groupby('md').mean()
        sd = df.groupby('md').std()

        mdays = mn.index.values

        constraints = {'clower': np.array(mn.to_numpy().flatten() - sd.to_numpy().flatten(), dtype='int16'),
                       'cupper': np.array(mn.to_numpy().flatten() + sd.to_numpy().flatten(), dtype='int16'),
                       'mday_set': list(mdays)}

        wc = array('f', np.arange(0.9, -0.1, -0.1))

        np.testing.assert_equal(self.data['constraints']['clower'], constraints['clower'], 0)
        np.testing.assert_equal(self.data['constraints']['cupper'], constraints['cupper'], 0)
        np.testing.assert_equal(self.data['constraints']['mday_set'], constraints['mday_set'], 0)

        zall = []
        ii = 0
        for jj in range(len(self.y)):

            if self.dates[jj] < date(2014, 1, 1):
                ii += 1
                continue

            kk = jj - 10

            y = np.array(self.data['y'][0:jj], dtype='double') * 10000
            w = np.array((y != -3000)*1, dtype='double')

            z = np.array(ws2dp(y, self.data['sopt_ws2dvcp'], w, 0.90))

            constraint_ix = [constraints['mday_set'].index(x) for x in [int(x.strftime('%m%d')) for x in self.dates[kk:jj]]]

            z[...] = w_constrain(z, constraints['clower'][constraint_ix], constraints['cupper'][constraint_ix], wc)

            zall.append(z)

        self.assertEqual(len(zall), len(self.data['z_ws2dp_constrain_update']))

        for ix in range(len(zall)): #pylint: disable=C0200
            self.assertEqual(len(zall[ix]), len(self.data['z_ws2dp_constrain_update'][ix]))
            np.testing.assert_almost_equal(zall[ix], self.data['z_ws2dp_constrain_update'][ix], 5)

    def test_ws2dvc(self):
        """Test ws2doptv (V-CURVE) smoothing."""
        z, sopt = ws2doptv(self.y, self.w, array('d', np.linspace(-2, 1, 16)))
        np.testing.assert_almost_equal(z, self.data['z_ws2dvc'], 5)
        self.assertEqual(sopt, self.data['sopt_ws2dvc'])

    def test_ws2dvcp(self):
        """Test ws2doptvp (V-CURVE with p) smoothing."""
        z, sopt = ws2doptvp(self.y, self.w, array('d', np.linspace(-2, 1, 16)), p=0.90)
        np.testing.assert_almost_equal(z, self.data['z_ws2dvcp'], 5)
        self.assertEqual(sopt, self.data['sopt_ws2dvcp'])

if __name__ == '__main__':
    unittest.main()
