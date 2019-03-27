from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import multiprocessing as mp
import os
import pickle
import unittest

import numpy as np

from modape.utils import *

class TestUtils(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        with open('{}/data/MXD_dates.pkl'.format(os.path.dirname(__file__).replace('tests','modape')),'rb') as pkl:
            cls.dates = pickle.load(pkl)

    @classmethod
    def tearDownClass(cls):
        cls.dates = None


    def test_datehelper(self):

        dh = DateHelper(self.dates,8,10)
        dix = dh.getDIX()
        dv = dh.getDV(-3000)

        self.assertEqual(len(dh.daily),5893)
        self.assertEqual(len(dv),5893)
        self.assertEqual(len(dh.target),580)
        self.assertEqual(len(dix),580)
        self.assertEqual([dh.daily[x] for x in dix],dh.target)
        self.assertTrue(np.all(dv == -3000))

    def test_credentials(self):

        cred = Credentials(username = 'testuser', password = 'testpass')

        self.assertEqual(cred.username,'testuser')
        self.assertEqual(cred.password,'testpass')

        cred.store()

        self.assertTrue(os.path.exists('modape.cred.pkl'))
        self.assertTrue(os.path.exists('modape.key.pkl'))

        cred.username = None
        cred.password = None

        cred.retrieve()

        self.assertEqual(cred.username,'testuser')
        self.assertEqual(cred.password,'testpass')

        cred.destroy()

        self.assertFalse(os.path.exists('modape.cred.pkl'))
        self.assertFalse(os.path.exists('modape.key.pkl'))


    def test_fj_ldom(self):

        test_day = LDOM(fromjulian('2016032'))

        self.assertEqual(test_day.year,2016)
        self.assertEqual(test_day.month,2)
        self.assertEqual(test_day.day,29)


    def test_tvec(self):

        self.assertEqual(tvec(2003,8),[x for x in self.dates if '2003' in x])

    def test_mp(self):

        def _init(arr_):

            global arr
            arr = tonumpyarray(arr_)

        def _work(ix):
            arr[ix] = ix
        return(mp.current_process().pid)

        arr_shared = init_shared(10)

        pool = mp.Pool(2,initializer=_init,initargs=(arr_shared,))

        res = pool.map(_work,range(10))

        pool.close()

        pool = None

        self.assertEqual(arr_shared[0:10],[x for x in range(10)])
        self.assertEqual(len(set(res)),2)

        del arr_shared


if __name__ == '__main__':
    unittest.main()
