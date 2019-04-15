"""test_utils.py: Test uilility classes and functions."""
# pylint: disable=global-variable-undefined, invalid-name
from __future__ import absolute_import, division, print_function

import os
import pickle
import unittest

import numpy as np

from modape.utils import DateHelper, Credentials, ldom, tvec, fromjulian

class TestUtils(unittest.TestCase):
    """Test class for testing utils."""

    @classmethod
    def setUpClass(cls):
        with open('{}/data/MXD_dates.pkl'.format(os.path.dirname(__file__).replace('tests', 'modape')), 'rb') as pkl:
            cls.dates = pickle.load(pkl)

    @classmethod
    def tearDownClass(cls):
        cls.dates = None

    def test_datehelper(self):
        """Testing DateHelper class."""
        dh = DateHelper(self.dates, 8, 10)
        dix = dh.getDIX()
        dv = dh.getDV(-3000)

        self.assertEqual(len(dh.daily), 5893)
        self.assertEqual(len(dv), 5893)
        self.assertEqual(len(dh.target), 580)
        self.assertEqual(len(dix), 580)
        self.assertEqual([dh.daily[x] for x in dix], dh.target)
        self.assertTrue(np.all(dv == -3000))

    def test_credentials(self):
        """Testing Credentials class."""
        cred = Credentials(username='testuser', password='testpass')

        self.assertEqual(cred.username, 'testuser')
        self.assertEqual(cred.password, 'testpass')
        cred.store()
        self.assertTrue(os.path.exists('modape.cred.pkl'))
        self.assertTrue(os.path.exists('modape.key.pkl'))
        cred.username = None
        cred.password = None
        cred.retrieve()
        self.assertEqual(cred.username, 'testuser')
        self.assertEqual(cred.password, 'testpass')
        cred.destroy()
        self.assertFalse(os.path.exists('modape.cred.pkl'))
        self.assertFalse(os.path.exists('modape.key.pkl'))

    def test_fj_ldom(self):
        """Testing ldom and fromjulian."""
        test_day = ldom(fromjulian('2016032'))

        self.assertEqual(test_day.year, 2016)
        self.assertEqual(test_day.month, 2)
        self.assertEqual(test_day.day, 29)

    def test_tvec(self):
        """Testing tvec."""
        self.assertEqual(tvec(2003, 8), [x for x in self.dates if '2003' in x])

if __name__ == '__main__':
    unittest.main()
