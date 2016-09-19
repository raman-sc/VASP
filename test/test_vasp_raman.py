# -*- coding: utf-8 -*-
import os
import time
import unittest
import vasp_raman

class VaspRamanTester(unittest.TestCase):
    def testT(self):
        m = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        mref = [[1, 4, 7], [2, 5, 8], [3, 6, 9]]

        mres = vasp_raman.T(m)
        for i in range(len(m)):
            self.assertSequenceEqual(mref[i], mres[i])
