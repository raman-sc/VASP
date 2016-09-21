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

    def testMATmVEC(self):
        m = [[1.1, 2.1, -3.1], [-5.1, 2.3, 4.1], [0.0, 2.1, -10.1]]
        v = [1.1, -2.1, -3.2]
        mref = [6.72, -23.56, 27.91]

        mres = vasp_raman.MAT_m_VEC(m, v)
        for a, b in zip(mref, mres):
            self.assertAlmostEqual(a, b, 3)

    def testParsePoscar(self):
        with open(os.path.join('test', 'POSCAR_1')) as poscar_fh:
            nat, vol, b, positions, poscar_header = vasp_raman.parse_poscar(poscar_fh)

            natref = 92
            self.assertEqual(natref, nat)
            self.assertAlmostEqual(1188.44, vol, 2)

            bref = [[5.808, 0.0, 0.0], [0.0, 8.198, 0.0], [-3.360, 0.0, 24.962]]
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual(bref[i][j], b[i][j], 2)

            self.assertEqual(natref, len(positions))

    def testParseEnvParams(self):
        params = '1_2_3_-0.35'
        ref = [1, 2, 3, -0.35]
        self.assertSequenceEqual(ref, vasp_raman.parse_env_params(params))

        params = '-1_5_30_21.35'
        ref = [-1, 5, 30, 21.35]
        self.assertSequenceEqual(ref, vasp_raman.parse_env_params(params))
