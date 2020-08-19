"""test_whittaker.py: Test core whittaker functions."""
# pylint: disable=invalid-name,E0401,E0611
from array import array
import os
import pickle
import unittest

import numpy as np
from modape.whittaker import lag1corr, ws2d, ws2dp, ws2doptv, ws2doptvp

class TestWhittaker(unittest.TestCase):
    """Test class for core whittaker functions."""

    @classmethod
    def setUpClass(cls):
        with open("{}/data/MXD_testdata.pkl".format(os.path.dirname(__file__).replace("tests", "modape")), "rb") as pkl:
            cls.data = pickle.load(pkl)

    @classmethod
    def tearDownClass(cls):
        cls.data = None

    def setUp(self):
        self.y = self.data["y"]
        self.w = self.data["w"]

    def tearDown(self):
        self.y = None
        self.w = None

    def test_lag1corr(self):
        """Test lag-1 correlation function."""
        self.assertAlmostEqual(lag1corr(self.y[:-1], self.y[1:], -3000.0), self.data["lag1corr"])

    def test_ws2d(self):
        """Test ws2d smoothing."""
        z = np.array(ws2d(self.y, 10, self.w), dtype="double")
        np.testing.assert_almost_equal(z, self.data["z_ws2d"], 5)

    def test_ws2dp(self):
        """Test ws2dp smoothing."""
        z = ws2dp(self.y, self.data["sopt_ws2dvcp"], self.w, p=0.90)
        np.testing.assert_almost_equal(z, self.data["z_ws2dp"], 5)

    def test_ws2dvc(self):
        """Test ws2doptv (V-CURVE) smoothing."""
        z, sopt = ws2doptv(self.y, self.w, array("d", np.arange(-2, 1.2, 0.2).round(2)))
        np.testing.assert_almost_equal(z, self.data["z_ws2dvc"], 5)
        np.testing.assert_almost_equal(sopt, self.data["sopt_ws2dvc"], 5)

    def test_ws2dvcp(self):
        """Test ws2doptvp (V-CURVE with p) smoothing."""
        z, sopt = ws2doptvp(self.y, self.w, array("d", np.arange(-2, 1.2, 0.2).round(2)), p=0.90)
        np.testing.assert_almost_equal(z, self.data["z_ws2dvcp"], 5)
        np.testing.assert_almost_equal(sopt, self.data["sopt_ws2dvcp"], 5)

if __name__ == "__main__":
    unittest.main()
