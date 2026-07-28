# -*- coding: utf-8 -*-
import json
import os
import tempfile
import unittest
from io import StringIO
from pathlib import Path

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
            nat, vol, b, positions, frozen, poscar_header = \
                vasp_raman.parse_poscar(poscar_fh)

            natref = 92
            self.assertEqual(natref, nat)
            self.assertAlmostEqual(1188.44, vol, 2)

            bref = [[5.808, 0.0, 0.0], [0.0, 8.198, 0.0],
                    [-3.360, 0.0, 24.962]]
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual(bref[i][j], b[i][j], 2)

            self.assertEqual(natref, len(positions))
            self.assertEqual(natref, len(frozen))
            self.assertListEqual([False] * natref, frozen)

    def testParsePoscarFrozen(self):
        with open(os.path.join('test', 'POSCAR_2')) as poscar_fh:
            nat, vol, b, positions, frozen, poscar_header = \
                vasp_raman.parse_poscar(poscar_fh)

            natref = 92
            self.assertEqual(natref, nat)
            self.assertAlmostEqual(1188.44, vol, 2)

            bref = [[5.808, 0.0, 0.0], [0.0, 8.198, 0.0],
                    [-3.360, 0.0, 24.962]]
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual(bref[i][j], b[i][j], 2)

            self.assertEqual(natref, len(positions))
            self.assertEqual(natref, len(frozen))

            ref_frozen = [False] * natref
            ref_frozen[1] = ref_frozen[5] = True
            self.assertListEqual(ref_frozen, frozen)

    def testPposcarHeaderSelectiveDynamics(self):

        for poscar_fn, is_selective in zip(
            ('POSCAR_1', 'POSCAR_2', 'POSCAR_3'), (False, False, True),
                ):

            with open(os.path.join('test', poscar_fn)) as poscar_fh:
                nat, vol, b, positions, frozen, poscar_header = \
                    vasp_raman.parse_poscar(poscar_fh)

            self.assertTrue(
                vasp_raman.is_selective_dynamics(poscar_header) is is_selective
            )

    def testParseEnvParams(self):
        params = '1_2_3_-0.35'
        ref = [1, 2, 3, -0.35]
        self.assertSequenceEqual(ref, vasp_raman.parse_env_params(params))

        params = '-1_5_30_21.35'
        ref = [-1, 5, 30, 21.35]
        self.assertSequenceEqual(ref, vasp_raman.parse_env_params(params))
    def testParsePoscarSymbols(self):
        poscar = StringIO(
            """Si O test
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Si O
1 2
Direct
0.0 0.0 0.0
0.5 0.5 0.5
0.25 0.25 0.25
"""
        )

        symbols = vasp_raman.parse_poscar_symbols(poscar)

        self.assertEqual(
            symbols,
            ['Si', 'O', 'O'],
        )


    def testRealPhonopyEigenvector(self):
        eigenvector = [
            [
                [0.0, 1.0],
                [0.0, 0.0],
                [0.0, 0.0],
            ],
            [
                [0.0, -1.0],
                [0.0, 0.0],
                [0.0, 0.0],
            ],
        ]

        result = vasp_raman._real_phonopy_eigenvector(
            eigenvector
        )

        self.assertAlmostEqual(
            result[0][0],
            1.0,
        )
        self.assertAlmostEqual(
            result[1][0],
            -1.0,
        )


    def testParsePhonopyData(self):
        data = {
            'points': [
                {
                    'symbol': 'He',
                    'mass': 4.0,
                },
            ],
            'phonon': [
                {
                    'q-position': [0.0, 0.0, 0.0],
                    'band': [
                        {
                            'frequency': 3.0,
                            'eigenvector': [
                                [
                                    [2.0, 0.0],
                                    [0.0, 0.0],
                                    [0.0, 0.0],
                                ],
                            ],
                        },
                        {
                            'frequency': 2.0,
                            'eigenvector': [
                                [
                                    [0.0, 0.0],
                                    [2.0, 0.0],
                                    [0.0, 0.0],
                                ],
                            ],
                        },
                        {
                            'frequency': 0.0,
                            'eigenvector': [
                                [
                                    [0.0, 0.0],
                                    [0.0, 0.0],
                                    [2.0, 0.0],
                                ],
                            ],
                        },
                    ],
                },
            ],
        }

        modes = vasp_raman.parse_phonopy_data(
            data,
            expected_symbols=['He'],
        )

        self.assertEqual(
            len(modes),
            3,
        )
        self.assertEqual(
            modes[0]['source_mode'],
            1,
        )
        self.assertAlmostEqual(
            modes[0]['frequency_cm-1'],
            3.0 * vasp_raman.THZ_TO_CM,
        )

        # 2 / sqrt(4) = 1
        self.assertAlmostEqual(
            modes[0]['eigenvector'][0][0],
            1.0,
        )

        non_gamma = {
            'points': [
                {
                    'symbol': 'He',
                    'mass': 4.0,
                },
            ],
            'phonon': [
                {
                    'q-position': [0.5, 0.0, 0.0],
                    'band': [],
                },
            ],
        }

        with self.assertRaisesRegex(
            ValueError,
            'not Gamma',
        ):
            vasp_raman.parse_phonopy_data(
                non_gamma,
                expected_symbols=['He'],
            )


    def testParseAseData(self):
        data = [
            {
                'mode': 1,
                'frequency_cm-1': 503.8,
                'eigenvector': [[0.1, 0.0, 0.0]],
            },
            {
                'mode': 2,
                'frequency_cm-1': 200.0,
                'eigenvector': [[0.0, 0.1, 0.0]],
            },
            {
                'mode': 3,
                'frequency_cm-1': 0.0,
                'eigenvector': [[0.0, 0.0, 0.1]],
            },
        ]

        modes = vasp_raman.parse_ase_data(
            data,
            nat=1,
        )

        self.assertEqual(
            len(modes),
            3,
        )
        self.assertEqual(
            modes[0]['source_mode'],
            1,
        )
        self.assertAlmostEqual(
            modes[0]['frequency_cm-1'],
            503.8,
        )
        self.assertEqual(
            modes[0]['eigenvector'],
            [[0.1, 0.0, 0.0]],
        )

        ase_fh = StringIO(
            json.dumps({'modes': data})
        )

        file_modes = vasp_raman.read_ase_modes(
            ase_fh,
            nat=1,
        )

        self.assertEqual(
            file_modes,
            modes,
        )

        with self.assertRaisesRegex(
            ValueError,
            'expected 3',
        ):
            vasp_raman.parse_ase_data(
                data[:2],
                nat=1,
            )


    def testValidateExternalModes(self):
        modes = [
            {
                'frequency_cm-1': 100.0,
                'eigenvector': [[1.0, 0.0, 0.0]],
            },
            {
                'frequency_cm-1': 50.0,
                'eigenvector': [[0.0, 2.0, 0.0]],
            },
            {
                'frequency_cm-1': 0.0,
                'eigenvector': [[0.0, 0.0, 0.5]],
            },
        ]

        eigvals, eigvecs, norms = (
            vasp_raman.validate_external_modes(
                modes,
                nat=1,
            )
        )

        self.assertEqual(
            eigvals,
            [100.0, 50.0, 0.0],
        )
        self.assertEqual(
            eigvecs,
            [
                [[1.0, 0.0, 0.0]],
                [[0.0, 2.0, 0.0]],
                [[0.0, 0.0, 0.5]],
            ],
        )
        self.assertEqual(
            norms,
            [1.0, 2.0, 0.5],
        )

        zero_norm = [
            modes[0],
            modes[1],
            {
                'frequency_cm-1': 0.0,
                'eigenvector': [[0.0, 0.0, 0.0]],
            },
        ]

        with self.assertRaisesRegex(
            ValueError,
            'zero norm',
        ):
            vasp_raman.validate_external_modes(
                zero_norm,
                nat=1,
            )


    def testWriteExternalModeFiles(self):
        modes = [
            {
                'source_mode': 1,
                'frequency_cm-1': 0.0,
                'eigenvector': [[0.0, 0.0, 0.5]],
            },
            {
                'source_mode': 2,
                'frequency_cm-1': 100.0,
                'eigenvector': [[1.0, 0.0, 0.0]],
            },
            {
                'source_mode': 3,
                'frequency_cm-1': 50.0,
                'eigenvector': [[0.0, 2.0, 0.0]],
            },
        ]

        ordered = vasp_raman.order_external_modes(
            modes
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            freq_path, modes_path, order_path = (
                vasp_raman.write_external_mode_files(
                    ordered,
                    output_dir=tmpdir,
                )
            )

            self.assertTrue(
                Path(freq_path).is_file()
            )
            self.assertTrue(
                Path(modes_path).is_file()
            )
            self.assertTrue(
                Path(order_path).is_file()
            )

            with open(freq_path) as freq_fh:
                frequencies = vasp_raman.parse_freqdat(
                    freq_fh,
                    nat=1,
                )

            with open(modes_path) as modes_fh:
                eigenvectors, norms = (
                    vasp_raman.parse_modesdat(
                        modes_fh,
                        nat=1,
                    )
                )

            self.assertEqual(
                frequencies,
                [100.0, 50.0, 0.0],
            )
            self.assertEqual(
                eigenvectors,
                [
                    [[1.0, 0.0, 0.0]],
                    [[0.0, 2.0, 0.0]],
                    [[0.0, 0.0, 0.5]],
                ],
            )
            self.assertEqual(
                norms,
                [1.0, 2.0, 0.5],
            )

            order_lines = Path(
                order_path
            ).read_text().splitlines()

            self.assertEqual(
                order_lines[0],
                'new_mode,source_mode,frequency_cm-1',
            )
            self.assertTrue(
                order_lines[1].startswith(
                    '1,2,100.'
                )
            )
    def testReadPhonopyQpointsWithCompanionMetadata(self):
        qpoints_text = """nqpoint: 1
natom: 1
phonon:
- q-position: [0.0, 0.0, 0.0]
  band:
  - frequency: 3.0
    eigenvector:
    - [[2.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
  - frequency: 2.0
    eigenvector:
    - [[0.0, 0.0], [2.0, 0.0], [0.0, 0.0]]
  - frequency: 0.0
    eigenvector:
    - [[0.0, 0.0], [0.0, 0.0], [2.0, 0.0]]
"""

        metadata_text = """primitive_cell:
  points:
  - symbol: He
    coordinates: [0.0, 0.0, 0.0]
    mass: 4.0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            qpoints_path = Path(
                tmpdir,
                'qpoints.yaml',
            )
            metadata_path = Path(
                tmpdir,
                'phonopy.yaml',
            )

            qpoints_path.write_text(
                qpoints_text
            )
            metadata_path.write_text(
                metadata_text
            )

            with qpoints_path.open() as qpoints_fh:
                modes = (
                    vasp_raman.read_phonopy_qpoints(
                        qpoints_fh,
                        expected_symbols=['He'],
                    )
                )

        self.assertEqual(
            len(modes),
            3,
        )
        self.assertAlmostEqual(
            modes[0]['frequency_cm-1'],
            3.0 * vasp_raman.THZ_TO_CM,
        )

        # 2 / sqrt(4) = 1
        self.assertAlmostEqual(
            modes[0]['eigenvector'][0][0],
            1.0,
        )