# -*- coding: utf-8 -*-
import json
import math
import os
import sys
import tempfile
import types
import unittest

from unittest import mock
from contextlib import redirect_stdout
from io import StringIO

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
            (
                nat,
                vol,
                b,
                positions,
                frozen,
                poscar_header,
                symbols,
            ) = vasp_raman.parse_poscar(poscar_fh)

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
            self.assertIsNotNone(symbols)
            self.assertEqual(len(symbols), nat)

    def testParsePoscarFrozen(self):
        with open(os.path.join('test', 'POSCAR_2')) as poscar_fh:
            (
                nat,
                vol,
                b,
                positions,
                frozen,
                poscar_header,
                _symbols,
            ) = vasp_raman.parse_poscar(poscar_fh)

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
                (
                    nat,
                    vol,
                    b,
                    positions,
                    frozen,
                    poscar_header,
                    _symbols,
                ) = vasp_raman.parse_poscar(poscar_fh)

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

    def testParsePoscarReturnsSymbols(self):
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

        (
            nat,
            volume,
            lattice,
            positions,
            frozen,
            header,
            symbols,
        ) = vasp_raman.parse_poscar(poscar)

        self.assertEqual(nat, 3)
        self.assertEqual(symbols, ['Si', 'O', 'O'])
        self.assertEqual(len(positions), 3)
        self.assertEqual(len(frozen), 3)

    def testParsePoscarRejectsMissingCoordinates(self):
        poscar = StringIO(
            """Incomplete structure
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Si
2
Direct
0.0 0.0 0.0
"""
        )

        with self.assertRaisesRegex(
            ValueError,
            'coordinate lines',
        ):
            vasp_raman.parse_poscar(poscar)

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

        result = vasp_raman.real_phonopy_eigenvector(
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
            data
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
            ase_fh
        )

        self.assertEqual(
            file_modes,
            modes,
        )

        incomplete_modes = (
            vasp_raman.parse_ase_data(
                data[:2]
            )
        )

        self.assertEqual(
            len(incomplete_modes),
            2,
        )

        with self.assertRaisesRegex(
            ValueError,
            'expected 3',
        ):
            vasp_raman.validate_external_modes(
                incomplete_modes,
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
        wrong_atom_count = [
            modes[0],
            {
                'frequency_cm-1': 50.0,
                'eigenvector': [
                    [0.0, 1.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            },
            modes[2],
        ]

        with self.assertRaisesRegex(
            ValueError,
            'contains 2 atoms',
        ):
            vasp_raman.validate_external_modes(
                wrong_atom_count,
                nat=1,
            )

        wrong_vector_shape = [
            modes[0],
            {
                'frequency_cm-1': 50.0,
                'eigenvector': [[0.0, 1.0]],
            },
            modes[2],
        ]

        with self.assertRaisesRegex(
            ValueError,
            'non-3D vector',
        ):
            vasp_raman.validate_external_modes(
                wrong_vector_shape,
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
                os.path.isfile(freq_path)
            )
            self.assertTrue(
                os.path.isfile(modes_path)
            )
            self.assertTrue(
                os.path.isfile(order_path)
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

            with open(order_path) as order_fh:
                order_lines = (
                    order_fh.read().splitlines()
                )

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
            qpoints_path = os.path.join(
                tmpdir,
                'qpoints.yaml',
            )
            metadata_path = os.path.join(
                tmpdir,
                'phonopy.yaml',
            )

            with open(
                qpoints_path,
                'w',
            ) as qpoints_output:
                qpoints_output.write(
                    qpoints_text
                )

            with open(
                metadata_path,
                'w',
            ) as metadata_output:
                metadata_output.write(
                    metadata_text
                )

            with open(qpoints_path) as qpoints_fh:
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

    def testBuildParser(self):
        parser = vasp_raman.build_parser()

        args = parser.parse_args([
            '--ase-modes',
            'modes.json',
            '--phonon-poscar',
            'structure.vasp',
            '--keep-mode-order',
        ])

        self.assertEqual(
            args.ase_modes,
            'modes.json',
        )
        self.assertEqual(
            args.phonon_poscar,
            'structure.vasp',
        )
        self.assertTrue(
            args.keep_mode_order
        )


    def testValidateOptionsConversionOnly(self):
        parser = vasp_raman.build_parser()

        args = parser.parse_args([
            '--ase-modes',
            'modes.json',
            '--write-external-modes',
            'raman_inputs',
        ])

        configuration = (
            vasp_raman.validate_options(
                parser,
                args,
                environ={},
            )
        )

        self.assertTrue(
            configuration['conversion_only']
        )
        self.assertIsNone(
            configuration['vasp_raman_run']
        )
        self.assertIsNone(
            configuration['first']
        )


    def testValidateOptionsRejectsConflictingSources(self):
        parser = vasp_raman.build_parser()

        args = parser.parse_args([
            '--phonopy-qpoints',
            'qpoints.yaml',
            '--ase-modes',
            'modes.json',
        ])

        with self.assertRaises(SystemExit) as context:
            vasp_raman.validate_options(
                parser,
                args,
                environ={},
            )

        self.assertEqual(
            context.exception.code,
            2,
        )


    def testParseEnvParamsRejectsBadInput(self):
        with self.assertRaisesRegex(
            ValueError,
            'exactly four',
        ):
            vasp_raman.parse_env_params(
                '1_2_0.01'
            )

    def testMainAcceptsArgv(self):
        output = StringIO()

        with self.assertRaises(SystemExit) as context:
            with redirect_stdout(output):
                vasp_raman.main([
                    '--help',
                ])

        self.assertEqual(
            context.exception.code,
            0,
        )

        help_text = output.getvalue()

        self.assertIn(
            '--phonopy-qpoints',
            help_text,
        )
        self.assertIn(
            '--run-ase-vibrations',
            help_text,
        )

    def testLoadStructure(self):
        poscar_text = """Single Si atom
1.0
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
Si
1
Direct
0.0 0.0 0.0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(
                tmpdir,
                'POSCAR.phon',
            )

            with open(
                filename,
                'w',
            ) as poscar_fh:
                poscar_fh.write(
                    poscar_text
                )

            (
                nat,
                volume,
                lattice,
                positions,
                frozen,
                header,
                symbols,
            ) = vasp_raman.load_structure(
                filename
            )

        self.assertEqual(
            nat,
            1,
        )
        self.assertEqual(
            symbols,
            ['Si'],
        )
        self.assertEqual(
            len(positions),
            1,
        )
        self.assertEqual(
            frozen,
            [False],
        )

    def testLoadStructureRejectsMissingFile(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(
                tmpdir,
                'missing_POSCAR',
            )

            with self.assertRaisesRegex(
                ValueError,
                'Could not open input structure',
            ):
                vasp_raman.load_structure(
                    filename
                )

    def testLoadModesFromAseJson(self):
        ase_modes = [
            {
                'mode': 1,
                'frequency_cm-1': 10.0,
                'eigenvector': [[1.0, 0.0, 0.0]],
            },
            {
                'mode': 2,
                'frequency_cm-1': 30.0,
                'eigenvector': [[0.0, 1.0, 0.0]],
            },
            {
                'mode': 3,
                'frequency_cm-1': 20.0,
                'eigenvector': [[0.0, 0.0, 1.0]],
            },
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(
                tmpdir,
                'ase_modes.json',
            )

            with open(
                filename,
                'w',
            ) as ase_fh:
                json.dump(
                    ase_modes,
                    ase_fh,
                )

            parser = vasp_raman.build_parser()

            options = parser.parse_args([
                '--ase-modes',
                filename,
            ])

            (
                eigvals,
                eigvecs,
                norms,
                canonical_modes,
            ) = vasp_raman.load_modes(
                options,
                nat=1,
                frozen=[False],
                poscar_symbols=['Si'],
                direct_ase_settings=None,
            )

        self.assertEqual(
            eigvals,
            [30.0, 20.0, 10.0],
        )
        self.assertEqual(
            norms,
            [1.0, 1.0, 1.0],
        )
        self.assertEqual(
            len(canonical_modes),
            3,
        )

    def testLoadModesFromVtstFiles(self):
        old_cwd = os.getcwd()

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                os.chdir(tmpdir)

                with open(
                    'freq.dat',
                    'w',
                ) as freq_fh:
                    freq_fh.write(
                        "100.0\n"
                        "50.0\n"
                        "0.0\n"
                    )

                with open(
                    'modes_sqrt_amu.dat',
                    'w',
                ) as modes_fh:
                    modes_fh.write(
                        "1.0 0.0 0.0\n\n"
                        "0.0 2.0 0.0\n\n"
                        "0.0 0.0 0.5\n\n"
                    )

                parser = vasp_raman.build_parser()
                options = parser.parse_args([])

                (
                    eigvals,
                    eigvecs,
                    norms,
                    canonical_modes,
                ) = vasp_raman.load_modes(
                    options,
                    nat=1,
                    frozen=[False],
                    poscar_symbols=['Si'],
                    direct_ase_settings=None,
                )

        finally:
            os.chdir(old_cwd)

        self.assertEqual(
            eigvals,
            [100.0, 50.0, 0.0],
        )
        self.assertEqual(
            norms,
            [1.0, 2.0, 0.5],
        )
        self.assertIsNone(
            canonical_modes
        )

    def testLoadModesRejectsMissingSource(self):
        old_cwd = os.getcwd()

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                os.chdir(tmpdir)

                parser = vasp_raman.build_parser()
                options = parser.parse_args([])

                with self.assertRaisesRegex(
                    ValueError,
                    'Neither an external mode source',
                ):
                    vasp_raman.load_modes(
                        options,
                        nat=1,
                        frozen=[False],
                        poscar_symbols=['Si'],
                        direct_ase_settings=None,
                    )

        finally:
            os.chdir(old_cwd)

    def testBuildParserLegacyFlags(self):
        parser = vasp_raman.build_parser()

        options = parser.parse_args([
            '--gen',
            '--use-poscar',
        ])

        self.assertTrue(
            options.gen
        )
        self.assertTrue(
            options.use_poscar
        )

    def testWriteDisplacedPoscar(self):
        poscar_header = (
            "1.0\n"
            "1.0 0.0 0.0\n"
            "0.0 1.0 0.0\n"
            "0.0 0.0 1.0\n"
            "Si\n"
            "2\n"
            "Selective dynamics\n"
        )

        positions = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
        ]

        frozen = [
            False,
            True,
        ]

        eigenvector = [
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(
                tmpdir,
                'POSCAR',
            )

            vasp_raman.write_displaced_poscar(
                output_filename=filename,
                calculation_label='mode-test',
                step_size=0.01,
                poscar_header=poscar_header,
                positions=positions,
                frozen=frozen,
                eigenvector=eigenvector,
                norm=1.0,
                displacement=1,
            )

            (
                nat,
                volume,
                lattice,
                displaced_positions,
                displaced_frozen,
                header,
                symbols,
            ) = vasp_raman.load_structure(
                filename
            )

        self.assertEqual(
            nat,
            2,
        )
        self.assertEqual(
            symbols,
            ['Si', 'Si'],
        )

        self.assertAlmostEqual(
            displaced_positions[0][0],
            0.01,
        )
        self.assertAlmostEqual(
            displaced_positions[0][1],
            0.0,
        )

        # The frozen atom must not move.
        self.assertEqual(
            displaced_positions[1],
            [0.5, 0.5, 0.5],
        )
        self.assertEqual(
            displaced_frozen,
            [False, True],
        )

    def testGetDisplacedEpsilonReusesOutcar(self):
        old_cwd = os.getcwd()

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                os.chdir(tmpdir)

                outcar_filename = (
                    'OUTCAR.0001.+1.out'
                )

                with open(
                    outcar_filename,
                    'w',
                ) as outcar_fh:
                    outcar_fh.write(
                        "MACROSCOPIC STATIC DIELECTRIC TENSOR\n"
                        "------------------------------------\n"
                        "1.0 0.0 0.0\n"
                        "0.0 2.0 0.0\n"
                        "0.0 0.0 3.0\n"
                    )

                epsilon = (
                    vasp_raman.get_displaced_epsilon(
                        mode_number=1,
                        displacement=1,
                        step_size=0.01,
                        poscar_header='',
                        positions=[[0.0, 0.0, 0.0]],
                        frozen=[False],
                        eigenvector=[[1.0, 0.0, 0.0]],
                        norm=1.0,
                        # It must not run this command because
                        # the displaced OUTCAR already exists.
                        vasp_command='false',
                        generate_only=False,
                        use_poscar=False,
                    )
                )

        finally:
            os.chdir(old_cwd)

        self.assertEqual(
            epsilon,
            [
                [1.0, 0.0, 0.0],
                [0.0, 2.0, 0.0],
                [0.0, 0.0, 3.0],
            ],
        )


    def testGetDisplacedEpsilonGenerateOnly(self):
        old_cwd = os.getcwd()

        poscar_header = (
            "1.0\n"
            "1.0 0.0 0.0\n"
            "0.0 1.0 0.0\n"
            "0.0 0.0 1.0\n"
            "Si\n"
            "1\n"
        )

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                os.chdir(tmpdir)

                result = (
                    vasp_raman.get_displaced_epsilon(
                        mode_number=1,
                        displacement=1,
                        step_size=0.01,
                        poscar_header=poscar_header,
                        positions=[[0.0, 0.0, 0.0]],
                        frozen=[False],
                        eigenvector=[[1.0, 0.0, 0.0]],
                        norm=1.0,
                        vasp_command=None,
                        generate_only=True,
                        use_poscar=False,
                    )
                )

                generated_filename = (
                    'POSCAR.+1.out'
                )

                self.assertIsNone(
                    result
                )
                self.assertTrue(
                    os.path.isfile(
                        generated_filename
                    )
                )
                self.assertFalse(
                    os.path.exists('POSCAR')
                )

                (
                    nat,
                    volume,
                    lattice,
                    positions,
                    frozen,
                    header,
                    symbols,
                ) = vasp_raman.load_structure(
                    generated_filename
                )

        finally:
            os.chdir(old_cwd)

        self.assertEqual(
            nat,
            1,
        )
        self.assertEqual(
            symbols,
            ['Si'],
        )
        self.assertAlmostEqual(
            positions[0][0],
            0.01,
        )

    def testCalculateRamanProperties(self):
        raman_tensor = [
            [1.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 3.0],
        ]

        (
            alpha,
            beta2,
            activity,
        ) = vasp_raman.calculate_raman_properties(
            raman_tensor
        )

        self.assertAlmostEqual(
            alpha,
            2.0,
        )
        self.assertAlmostEqual(
            beta2,
            3.0,
        )
        self.assertAlmostEqual(
            activity,
            201.0,
        )

    def testRunRamanCalculationReusesOutcars(self):
        old_cwd = os.getcwd()

        def write_outcar(
            filename,
            diagonal,
        ):
            with open(
                filename,
                'w',
            ) as outcar_fh:
                outcar_fh.write(
                    "MACROSCOPIC STATIC DIELECTRIC TENSOR\n"
                    "------------------------------------\n"
                    f"{diagonal[0]} 0.0 0.0\n"
                    f"0.0 {diagonal[1]} 0.0\n"
                    f"0.0 0.0 {diagonal[2]}\n"
                )

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                os.chdir(tmpdir)

                write_outcar(
                    'OUTCAR.0001.-1.out',
                    [0.0, 0.0, 0.0],
                )
                write_outcar(
                    'OUTCAR.0001.+1.out',
                    [2.0, 4.0, 6.0],
                )

                parser = vasp_raman.build_parser()
                options = parser.parse_args([])

                step_size = 0.01

                configuration = {
                    'first': 1,
                    'last': 1,
                    'step_size': step_size,
                    'disps': [-1, 1],
                    'coeffs': [-0.5, 0.5],
                    # This command must not run because both
                    # displaced OUTCAR files already exist.
                    'vasp_raman_run': 'false',
                }

                results = (
                    vasp_raman.run_raman_calculation(
                        options=options,
                        configuration=configuration,
                        # This choice makes the Raman tensor
                        # equal to diag(1, 2, 3).
                        volume=(
                            4.0
                            * math.pi
                            * step_size
                        ),
                        positions=[
                            [0.0, 0.0, 0.0],
                        ],
                        frozen=[
                            False,
                        ],
                        poscar_header='',
                        eigvals=[
                            123.45,
                        ],
                        eigvecs=[
                            [
                                [1.0, 0.0, 0.0],
                            ],
                        ],
                        norms=[
                            1.0,
                        ],
                    )
                )

                with open(
                    'vasp_raman.dat',
                    'r',
                ) as output_fh:
                    output_lines = (
                        output_fh.read().splitlines()
                    )

        finally:
            os.chdir(old_cwd)

        self.assertEqual(
            len(results),
            1,
        )
        self.assertEqual(
            results[0]['mode'],
            1,
        )
        self.assertAlmostEqual(
            results[0]['frequency_cm-1'],
            123.45,
        )
        self.assertAlmostEqual(
            results[0]['alpha'],
            2.0,
        )
        self.assertAlmostEqual(
            results[0]['beta2'],
            3.0,
        )
        self.assertAlmostEqual(
            results[0]['activity'],
            201.0,
        )

        self.assertEqual(
            len(output_lines),
            2,
        )
        self.assertTrue(
            output_lines[1].startswith(
                '001'
            )
        )

    def testValidateAseCacheReuseRejectsChangedDelta(self):
        recorded_manifest = {
            'schema_version': 1,
            'structure_sha256': 'same-structure',
            'calculator_settings': {
                'encut': 500,
            },
            'delta_angstrom': 0.01,
            'nfree': 2,
            'vibrated_atom_indices': [
                0,
                1,
                2,
            ],
        }

        requested_manifest = dict(
            recorded_manifest
        )
        requested_manifest[
            'delta_angstrom'
        ] = 0.02

        with tempfile.TemporaryDirectory() as tmpdir:
            vibration_directory = os.path.join(
                tmpdir,
                'vib',
            )

            os.makedirs(
                vibration_directory
            )

            # A non-empty cache file represents a completed
            # displacement from an earlier calculation.
            cache_filename = os.path.join(
                vibration_directory,
                'cache.eq.json',
            )

            with open(
                cache_filename,
                'w',
            ) as cache_fh:
                json.dump(
                    {
                        'forces': [
                            [0.0, 0.0, 0.0],
                        ],
                    },
                    cache_fh,
                )

            manifest_filename = os.path.join(
                tmpdir,
                'ase_vibration_cache_manifest.json',
            )

            with open(
                manifest_filename,
                'w',
            ) as manifest_fh:
                json.dump(
                    recorded_manifest,
                    manifest_fh,
                    indent=2,
                )

            with self.assertRaisesRegex(
                ValueError,
                'delta_angstrom',
            ):
                vasp_raman.validate_ase_cache_reuse(
                    workdir=tmpdir,
                    expected_manifest=requested_manifest,
                )

    def testValidateAseCacheReuseWithoutCache(self):
        expected_manifest = {
            'schema_version': 1,
            'structure_sha256': 'structure-a',
            'calculator_settings': {
                'encut': 500,
            },
            'delta_angstrom': 0.01,
            'nfree': 2,
            'vibrated_atom_indices': [
                0,
                1,
            ],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            reused = (
                vasp_raman.validate_ase_cache_reuse(
                    workdir=tmpdir,
                    expected_manifest=expected_manifest,
                )
            )

        self.assertFalse(
            reused
        )


    def testValidateAseCacheReuseAcceptsMatchingManifest(self):
        expected_manifest = {
            'schema_version': 1,
            'structure_sha256': 'structure-a',
            'calculator_settings': {
                'encut': 500,
            },
            'delta_angstrom': 0.01,
            'nfree': 2,
            'vibrated_atom_indices': [
                0,
                1,
            ],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            vibration_directory = os.path.join(
                tmpdir,
                'vib',
            )

            os.makedirs(
                vibration_directory
            )

            cache_filename = os.path.join(
                vibration_directory,
                'cache.eq.json',
            )

            with open(
                cache_filename,
                'w',
            ) as cache_fh:
                json.dump(
                    {
                        'forces': [
                            [0.0, 0.0, 0.0],
                        ],
                    },
                    cache_fh,
                )

            manifest_filename = os.path.join(
                tmpdir,
                'ase_vibration_cache_manifest.json',
            )

            with open(
                manifest_filename,
                'w',
            ) as manifest_fh:
                json.dump(
                    expected_manifest,
                    manifest_fh,
                    indent=2,
                )

            reused = (
                vasp_raman.validate_ase_cache_reuse(
                    workdir=tmpdir,
                    expected_manifest=expected_manifest,
                )
            )

        self.assertTrue(
            reused
        )


    def testValidateAseCacheReuseRejectsMissingManifest(self):
        expected_manifest = {
            'schema_version': 1,
            'structure_sha256': 'structure-a',
            'calculator_settings': {
                'encut': 500,
            },
            'delta_angstrom': 0.01,
            'nfree': 2,
            'vibrated_atom_indices': [
                0,
            ],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            vibration_directory = os.path.join(
                tmpdir,
                'vib',
            )

            os.makedirs(
                vibration_directory
            )

            cache_filename = os.path.join(
                vibration_directory,
                'cache.eq.json',
            )

            with open(
                cache_filename,
                'w',
            ) as cache_fh:
                json.dump(
                    {
                        'forces': [
                            [0.0, 0.0, 0.0],
                        ],
                    },
                    cache_fh,
                )

            with self.assertRaisesRegex(
                ValueError,
                'no compatibility manifest',
            ):
                vasp_raman.validate_ase_cache_reuse(
                    workdir=tmpdir,
                    expected_manifest=expected_manifest,
                )


    def testValidateAseCacheReuseRejectsEmptyCacheFile(self):
        expected_manifest = {
            'schema_version': 1,
            'structure_sha256': 'structure-a',
            'calculator_settings': {
                'encut': 500,
            },
            'delta_angstrom': 0.01,
            'nfree': 2,
            'vibrated_atom_indices': [
                0,
            ],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            vibration_directory = os.path.join(
                tmpdir,
                'vib',
            )

            os.makedirs(
                vibration_directory
            )

            cache_filename = os.path.join(
                vibration_directory,
                'cache.eq.json',
            )

            # An interrupted ASE calculation may leave an
            # empty JSON cache file.
            with open(
                cache_filename,
                'w',
            ):
                pass

            manifest_filename = os.path.join(
                tmpdir,
                'ase_vibration_cache_manifest.json',
            )

            with open(
                manifest_filename,
                'w',
            ) as manifest_fh:
                json.dump(
                    expected_manifest,
                    manifest_fh,
                    indent=2,
                )

            with self.assertRaisesRegex(
                ValueError,
                'empty files',
            ):
                vasp_raman.validate_ase_cache_reuse(
                    workdir=tmpdir,
                    expected_manifest=expected_manifest,
                )

    def testBuildAseCacheManifest(self):
        class FakeConstraint:
            def todict(self):
                return {
                    'name': 'FixAtoms',
                    'kwargs': {
                        'indices': [
                            1,
                        ],
                    },
                }

        class FakeAtoms:
            def __init__(self):
                self.positions = [
                    [0.0, 0.0, 0.0],
                    [1.0, 2.0, 3.0],
                ]

                self.constraints = [
                    FakeConstraint(),
                ]

            def __len__(self):
                return 2

            def get_atomic_numbers(self):
                return [
                    14,
                    8,
                ]

            def get_positions(self):
                return self.positions

            def get_cell(self):
                return [
                    [5.0, 0.0, 0.0],
                    [0.0, 5.0, 0.0],
                    [0.0, 0.0, 5.0],
                ]

            def get_pbc(self):
                return [
                    True,
                    True,
                    True,
                ]

            def get_masses(self):
                return [
                    28.085,
                    15.999,
                ]

            def get_initial_magnetic_moments(self):
                return [
                    0.0,
                    0.0,
                ]

            def get_initial_charges(self):
                return [
                    0.0,
                    0.0,
                ]

        atoms = FakeAtoms()

        settings_a = {
            'encut': 500,
            'kpts': [
                1,
                1,
                1,
            ],
        }

        # The insertion order differs, but the manifest
        # should remain deterministic.
        settings_b = {
            'kpts': [
                1,
                1,
                1,
            ],
            'encut': 500,
        }

        manifest_a = (
            vasp_raman.build_ase_cache_manifest(
                atoms=atoms,
                settings=settings_a,
                delta=0.01,
                nfree=2,
                vibrated_atom_indices=[
                    0,
                ],
            )
        )

        manifest_b = (
            vasp_raman.build_ase_cache_manifest(
                atoms=FakeAtoms(),
                settings=settings_b,
                delta=0.01,
                nfree=2,
                vibrated_atom_indices=[
                    0,
                ],
            )
        )

        self.assertEqual(
            manifest_a,
            manifest_b,
        )

        settings_with_command = dict(
            settings_a
        )
        settings_with_command[
            'command'
        ] = 'false'

        manifest_with_command = (
            vasp_raman.build_ase_cache_manifest(
                atoms=FakeAtoms(),
                settings=settings_with_command,
                delta=0.01,
                nfree=2,
                vibrated_atom_indices=[
                    0,
                ],
            )
        )

        self.assertEqual(
            manifest_a,
            manifest_with_command,
        )

        self.assertRegex(
            manifest_a['structure_sha256'],
            r'^[0-9a-f]{64}$',
        )

        self.assertEqual(
            manifest_a[
                'vibrated_atom_indices'
            ],
            [
                0,
            ],
        )

        self.assertEqual(
            manifest_a[
                'structure'
            ][
                'constraints'
            ][0][
                'name'
            ],
            'FixAtoms',
        )

        # Confirm that the full manifest can be written as JSON.
        json.dumps(
            manifest_a,
            sort_keys=True,
            allow_nan=False,
        )

        changed_atoms = FakeAtoms()
        changed_atoms.positions[1][0] = 1.25

        changed_structure_manifest = (
            vasp_raman.build_ase_cache_manifest(
                atoms=changed_atoms,
                settings=settings_a,
                delta=0.01,
                nfree=2,
                vibrated_atom_indices=[
                    0,
                ],
            )
        )

        self.assertNotEqual(
            manifest_a['structure_sha256'],
            changed_structure_manifest[
                'structure_sha256'
            ],
        )

        changed_settings_manifest = (
            vasp_raman.build_ase_cache_manifest(
                atoms=FakeAtoms(),
                settings={
                    'encut': 600,
                    'kpts': [
                        1,
                        1,
                        1,
                    ],
                },
                delta=0.01,
                nfree=2,
                vibrated_atom_indices=[
                    0,
                ],
            )
        )

        self.assertEqual(
            manifest_a['structure_sha256'],
            changed_settings_manifest[
                'structure_sha256'
            ],
        )

        self.assertNotEqual(
            manifest_a[
                'calculator_settings'
            ],
            changed_settings_manifest[
                'calculator_settings'
            ],
        )

        with self.assertRaisesRegex(
            ValueError,
            'outside the structure',
        ):
            vasp_raman.build_ase_cache_manifest(
                atoms=FakeAtoms(),
                settings=settings_a,
                delta=0.01,
                nfree=2,
                vibrated_atom_indices=[
                    2,
                ],
            )

    def testPrepareAseCacheFreshAndReuse(self):
        class FakeAtoms:
            def __init__(self):
                self.constraints = []

            def __len__(self):
                return 1

            def get_atomic_numbers(self):
                return [
                    14,
                ]

            def get_positions(self):
                return [
                    [0.0, 0.0, 0.0],
                ]

            def get_cell(self):
                return [
                    [5.0, 0.0, 0.0],
                    [0.0, 5.0, 0.0],
                    [0.0, 0.0, 5.0],
                ]

            def get_pbc(self):
                return [
                    True,
                    True,
                    True,
                ]

            def get_masses(self):
                return [
                    28.085,
                ]

            def get_initial_magnetic_moments(self):
                return [
                    0.0,
                ]

            def get_initial_charges(self):
                return [
                    0.0,
                ]

        settings = {
            'encut': 500,
            'kpts': [
                1,
                1,
                1,
            ],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            (
                first_manifest,
                first_reused,
                manifest_filename,
            ) = vasp_raman.prepare_ase_cache(
                atoms=FakeAtoms(),
                settings=settings,
                delta=0.01,
                nfree=2,
                vibrated_atom_indices=[
                    0,
                ],
                workdir=tmpdir,
            )

            self.assertFalse(
                first_reused
            )

            self.assertTrue(
                os.path.isfile(
                    manifest_filename
                )
            )

            with open(
                manifest_filename,
                'r',
            ) as manifest_fh:
                written_manifest = json.load(
                    manifest_fh
                )

            self.assertEqual(
                written_manifest,
                first_manifest,
            )

            vibration_directory = os.path.join(
                tmpdir,
                'vib',
            )

            os.makedirs(
                vibration_directory
            )

            cache_filename = os.path.join(
                vibration_directory,
                'cache.eq.json',
            )

            with open(
                cache_filename,
                'w',
            ) as cache_fh:
                json.dump(
                    {
                        'forces': [
                            [0.0, 0.0, 0.0],
                        ],
                    },
                    cache_fh,
                )

            (
                second_manifest,
                second_reused,
                second_manifest_filename,
            ) = vasp_raman.prepare_ase_cache(
                atoms=FakeAtoms(),
                settings=settings,
                delta=0.01,
                nfree=2,
                vibrated_atom_indices=[
                    0,
                ],
                workdir=tmpdir,
            )

            self.assertTrue(
                second_reused
            )

            self.assertEqual(
                second_manifest,
                first_manifest,
            )

            self.assertEqual(
                second_manifest_filename,
                manifest_filename,
            )

    def testRunAseVibrationsRejectsIncompatibleCacheBeforeRun(self):
        class FakeAtoms:
            def __init__(self):
                self.constraints = []
                self.calc = None

            def __len__(self):
                return 1

            def get_atomic_numbers(self):
                return [
                    14,
                ]

            def get_positions(self):
                return [
                    [0.0, 0.0, 0.0],
                ]

            def get_cell(self):
                return [
                    [5.0, 0.0, 0.0],
                    [0.0, 5.0, 0.0],
                    [0.0, 0.0, 5.0],
                ]

            def get_pbc(self):
                return [
                    True,
                    True,
                    True,
                ]

            def get_masses(self):
                return [
                    28.085,
                ]

            def get_initial_magnetic_moments(self):
                return [
                    0.0,
                ]

            def get_initial_charges(self):
                return [
                    0.0,
                ]

        class FakeVasp:
            def __init__(
                self,
                *args,
                **kwargs,
            ):
                pass

        class FakeVibrations:
            run_called = False

            def __init__(
                self,
                atoms,
                name,
                delta,
                nfree,
                indices=None,
            ):
                self.atoms = atoms
                self.name = name
                self.delta = delta
                self.nfree = nfree
                self.indices = (
                    [0]
                    if indices is None
                    else list(indices)
                )

            def run(self):
                type(self).run_called = True

        fake_ase = types.ModuleType(
            'ase'
        )
        fake_ase.__version__ = 'test'
        fake_ase.__path__ = []

        fake_calculators = types.ModuleType(
            'ase.calculators'
        )
        fake_calculators.__path__ = []

        fake_vasp_module = types.ModuleType(
            'ase.calculators.vasp'
        )
        fake_vasp_module.Vasp = FakeVasp

        fake_io_module = types.ModuleType(
            'ase.io'
        )
        fake_io_module.read = (
            lambda *args, **kwargs: FakeAtoms()
        )

        fake_vibrations_module = types.ModuleType(
            'ase.vibrations'
        )
        fake_vibrations_module.Vibrations = (
            FakeVibrations
        )

        fake_constraints_module = types.ModuleType(
            'ase.constraints'
        )
        fake_constraints_module.FixAtoms = type(
            'FakeFixAtoms',
            (),
            {},
        )

        fake_ase.calculators = (
            fake_calculators
        )
        fake_ase.io = fake_io_module
        fake_ase.vibrations = (
            fake_vibrations_module
        )
        fake_ase.constraints = (
            fake_constraints_module
        )
        fake_calculators.vasp = (
            fake_vasp_module
        )

        fake_modules = {
            'ase': fake_ase,
            'ase.calculators': (
                fake_calculators
            ),
            'ase.calculators.vasp': (
                fake_vasp_module
            ),
            'ase.io': fake_io_module,
            'ase.vibrations': (
                fake_vibrations_module
            ),
            'ase.constraints': (
                fake_constraints_module
            ),
        }

        settings = {
            'encut': 500,
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            structure_filename = os.path.join(
                tmpdir,
                'POSCAR.phon',
            )

            with open(
                structure_filename,
                'w',
            ) as structure_fh:
                structure_fh.write(
                    'placeholder structure\n'
                )

            recorded_manifest = (
                vasp_raman.build_ase_cache_manifest(
                    atoms=FakeAtoms(),
                    settings=settings,
                    delta=0.01,
                    nfree=2,
                    vibrated_atom_indices=[
                        0,
                    ],
                )
            )

            vasp_raman.write_ase_cache_manifest(
                workdir=tmpdir,
                manifest=recorded_manifest,
            )

            vibration_directory = os.path.join(
                tmpdir,
                'vib',
            )

            os.makedirs(
                vibration_directory
            )

            cache_filename = os.path.join(
                vibration_directory,
                'cache.eq.json',
            )

            with open(
                cache_filename,
                'w',
            ) as cache_fh:
                json.dump(
                    {
                        'forces': [
                            [0.0, 0.0, 0.0],
                        ],
                    },
                    cache_fh,
                )

            with mock.patch.dict(
                sys.modules,
                fake_modules,
            ):
                with self.assertRaisesRegex(
                    ValueError,
                    'delta_angstrom',
                ):
                    vasp_raman.run_ase_vibrations(
                        structure_filename=(
                            structure_filename
                        ),
                        settings=settings,
                        workdir=tmpdir,
                        delta=0.02,
                        nfree=2,
                    )

        self.assertFalse(
            FakeVibrations.run_called
        )


    def testValidateExternalModesAllowsReducedModeCountForFixedAtoms(self):
        # Two atoms are present, but atom 0 is fully fixed.
        # ASE therefore returns three modes for the one
        # vibrated atom rather than six modes for both atoms.
        modes = [
            {
                'frequency_cm-1': 100.0,
                'eigenvector': [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                ],
            },
            {
                'frequency_cm-1': 50.0,
                'eigenvector': [
                    [0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ],
            },
            {
                'frequency_cm-1': 25.0,
                'eigenvector': [
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            },
        ]

        (
            eigvals,
            eigvecs,
            norms,
        ) = vasp_raman.validate_external_modes(
            modes=modes,
            nat=2,
            expected_mode_count=3,
        )

        self.assertEqual(
            eigvals,
            [
                100.0,
                50.0,
                25.0,
            ],
        )

        self.assertEqual(
            len(eigvecs),
            3,
        )

        self.assertEqual(
            eigvecs[0][0],
            [0.0, 0.0, 0.0],
        )

        self.assertEqual(
            norms,
            [
                1.0,
                1.0,
                1.0,
            ],
        )


    def testLoadModesDirectAseUsesReducedModeCount(self):
        direct_ase_modes = [
            {
                'source_mode': 1,
                'frequency_cm-1': 100.0,
                'eigenvector': [
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                ],
            },
            {
                'source_mode': 2,
                'frequency_cm-1': 50.0,
                'eigenvector': [
                    [0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ],
            },
            {
                'source_mode': 3,
                'frequency_cm-1': 25.0,
                'eigenvector': [
                    [0.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                ],
            },
        ]

        parser = vasp_raman.build_parser()

        options = parser.parse_args([
            '--run-ase-vibrations',
            '--ase-vasp-settings',
            'ase_vasp_settings.json',
            '--phonon-poscar',
            'POSCAR.phon',
        ])

        with mock.patch.object(
            vasp_raman,
            'run_ase_vibrations',
            return_value=(
                direct_ase_modes,
                'ase_gamma_modes.json',
                'ase_vibrations_summary.txt',
                [1],
            ),
        ) as mocked_run:
            (
                eigvals,
                eigvecs,
                norms,
                canonical_modes,
            ) = vasp_raman.load_modes(
                options=options,
                nat=2,
                frozen=[
                    True,
                    False,
                ],
                poscar_symbols=[
                    'Si',
                    'O',
                ],
                direct_ase_settings={
                    'encut': 500,
                },
            )

        mocked_run.assert_called_once()

        self.assertEqual(
            eigvals,
            [
                100.0,
                50.0,
                25.0,
            ],
        )

        self.assertEqual(
            len(canonical_modes),
            3,
        )

        self.assertEqual(
            len(eigvecs[0]),
            2,
        )

        self.assertEqual(
            eigvecs[0][0],
            [0.0, 0.0, 0.0],
        )

        self.assertEqual(
            norms,
            [
                1.0,
                1.0,
                1.0,
            ],
        )


    def testLoadModesFromAseJsonUsesFixedAtomMetadata(self):
        ase_data = {
            'metadata': {
                'schema_version': 1,
                'number_of_atoms': 2,
                'vibrated_atom_indices': [
                    1,
                ],
                'expected_mode_count': 3,
            },
            'modes': [
                {
                    'mode': 1,
                    'frequency_cm-1': 100.0,
                    'eigenvector': [
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                    ],
                },
                {
                    'mode': 2,
                    'frequency_cm-1': 50.0,
                    'eigenvector': [
                        [0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ],
                },
                {
                    'mode': 3,
                    'frequency_cm-1': 25.0,
                    'eigenvector': [
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0],
                    ],
                },
            ],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(
                tmpdir,
                'ase_fixed_atom_modes.json',
            )

            with open(
                filename,
                'w',
            ) as ase_fh:
                json.dump(
                    ase_data,
                    ase_fh,
                    indent=2,
                )

            parser = vasp_raman.build_parser()

            options = parser.parse_args([
                '--ase-modes',
                filename,
            ])

            (
                eigvals,
                eigvecs,
                norms,
                canonical_modes,
            ) = vasp_raman.load_modes(
                options=options,
                nat=2,
                frozen=[
                    True,
                    False,
                ],
                poscar_symbols=[
                    'Si',
                    'O',
                ],
                direct_ase_settings=None,
            )

        self.assertEqual(
            eigvals,
            [
                100.0,
                50.0,
                25.0,
            ],
        )

        self.assertEqual(
            len(canonical_modes),
            3,
        )

        self.assertEqual(
            len(eigvecs),
            3,
        )

        self.assertEqual(
            len(eigvecs[0]),
            2,
        )

        self.assertEqual(
            eigvecs[0][0],
            [0.0, 0.0, 0.0],
        )

        self.assertEqual(
            norms,
            [
                1.0,
                1.0,
                1.0,
            ],
        )

    def testValidateAseModeMetadataRejectsInconsistencies(self):
        valid_metadata = {
            'schema_version': 1,
            'number_of_atoms': 2,
            'vibrated_atom_indices': [
                1,
            ],
            'expected_mode_count': 3,
        }

        result = (
            vasp_raman.validate_ase_mode_metadata(
                metadata=valid_metadata,
                nat=2,
                frozen=[
                    True,
                    False,
                ],
            )
        )

        self.assertEqual(
            result,
            3,
        )

        invalid_atom_count = dict(
            valid_metadata
        )
        invalid_atom_count[
            'number_of_atoms'
        ] = 3

        with self.assertRaisesRegex(
            ValueError,
            'number_of_atoms',
        ):
            vasp_raman.validate_ase_mode_metadata(
                metadata=invalid_atom_count,
                nat=2,
                frozen=[
                    True,
                    False,
                ],
            )

        invalid_mode_count = dict(
            valid_metadata
        )
        invalid_mode_count[
            'expected_mode_count'
        ] = 6

        with self.assertRaisesRegex(
            ValueError,
            'expected_mode_count',
        ):
            vasp_raman.validate_ase_mode_metadata(
                metadata=invalid_mode_count,
                nat=2,
                frozen=[
                    True,
                    False,
                ],
            )

        invalid_indices = dict(
            valid_metadata
        )
        invalid_indices[
            'vibrated_atom_indices'
        ] = [
            0,
        ]

        with self.assertRaisesRegex(
            ValueError,
            'do not match',
        ):
            vasp_raman.validate_ase_mode_metadata(
                metadata=invalid_indices,
                nat=2,
                frozen=[
                    True,
                    False,
                ],
            )


    def testReadPhonopyQpointsRejectsNonMappingYaml(self):
        invalid_documents = {
            'empty': '',
            'list': '[]\n',
            'scalar': '42\n',
            'string': 'not-a-mapping\n',
        }

        for label, contents in invalid_documents.items():
            with self.subTest(label=label):
                with self.assertRaisesRegex(
                    ValueError,
                    'mapping',
                ):
                    vasp_raman.read_phonopy_qpoints(
                        StringIO(contents)
                    )


    def testValidateOptionsAcceptsAllAseVaspCommandVariables(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            settings_filename = os.path.join(
                tmpdir,
                'ase_vasp_settings.json',
            )

            with open(
                settings_filename,
                'w',
            ) as settings_fh:
                json.dump(
                    {},
                    settings_fh,
                )

            parser = vasp_raman.build_parser()

            args = parser.parse_args([
                '--run-ase-vibrations',
                '--ase-vasp-settings',
                settings_filename,
                '--write-external-modes',
                os.path.join(
                    tmpdir,
                    'output',
                ),
            ])

            for variable_name in (
                'ASE_VASP_COMMAND',
                'VASP_COMMAND',
                'VASP_SCRIPT',
            ):
                with self.subTest(
                    variable_name=variable_name
                ):
                    configuration = (
                        vasp_raman.validate_options(
                            parser,
                            args,
                            environ={
                                variable_name: 'false',
                            },
                        )
                    )

                    self.assertIsNotNone(
                        configuration[
                            'direct_ase_settings'
                        ]
                    )

    def testRunAseVibrationsPassesMovableIndicesExplicitly(self):
        captured = {}

        class FakeFixAtoms:
            def __init__(self, indices):
                self.indices = list(indices)

            def get_indices(self):
                return list(self.indices)

        class FakeAtoms:
            def __init__(self):
                self.constraints = [
                    FakeFixAtoms([
                        0,
                    ]),
                ]
                self.calc = None

            def __len__(self):
                return 2

        class FakeVasp:
            def __init__(
                self,
                *args,
                **kwargs,
            ):
                pass

        class FakeVibrations:
            def __init__(
                self,
                atoms,
                name,
                delta,
                nfree,
                indices=None,
            ):
                captured['indices'] = (
                    None
                    if indices is None
                    else [
                        int(index)
                        for index in indices
                    ]
                )

                # Simulate the old ASE default, where all atoms
                # were selected if indices was not supplied.
                self.indices = (
                    [0, 1]
                    if indices is None
                    else list(indices)
                )

        fake_ase = types.ModuleType(
            'ase'
        )
        fake_ase.__version__ = '3.22.1'
        fake_ase.__path__ = []

        fake_calculators = types.ModuleType(
            'ase.calculators'
        )
        fake_calculators.__path__ = []

        fake_vasp_module = types.ModuleType(
            'ase.calculators.vasp'
        )
        fake_vasp_module.Vasp = FakeVasp

        fake_io_module = types.ModuleType(
            'ase.io'
        )
        fake_io_module.read = (
            lambda *args, **kwargs: FakeAtoms()
        )

        fake_vibrations_module = types.ModuleType(
            'ase.vibrations'
        )
        fake_vibrations_module.Vibrations = (
            FakeVibrations
        )

        fake_constraints_module = types.ModuleType(
            'ase.constraints'
        )
        fake_constraints_module.FixAtoms = (
            FakeFixAtoms
        )

        fake_ase.calculators = (
            fake_calculators
        )
        fake_ase.io = fake_io_module
        fake_ase.vibrations = (
            fake_vibrations_module
        )
        fake_ase.constraints = (
            fake_constraints_module
        )

        fake_calculators.vasp = (
            fake_vasp_module
        )

        fake_modules = {
            'ase': fake_ase,
            'ase.calculators': (
                fake_calculators
            ),
            'ase.calculators.vasp': (
                fake_vasp_module
            ),
            'ase.io': fake_io_module,
            'ase.vibrations': (
                fake_vibrations_module
            ),
            'ase.constraints': (
                fake_constraints_module
            ),
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            structure_filename = os.path.join(
                tmpdir,
                'POSCAR.phon',
            )

            with open(
                structure_filename,
                'w',
            ) as structure_fh:
                structure_fh.write(
                    'placeholder structure\n'
                )

            with mock.patch.dict(
                sys.modules,
                fake_modules,
            ):
                with mock.patch.object(
                    vasp_raman,
                    'prepare_ase_cache',
                    side_effect=RuntimeError(
                        'stop-after-constructor'
                    ),
                ):
                    with self.assertRaisesRegex(
                        RuntimeError,
                        'stop-after-constructor',
                    ):
                        vasp_raman.run_ase_vibrations(
                            structure_filename=(
                                structure_filename
                            ),
                            settings={
                                'encut': 500,
                            },
                            workdir=tmpdir,
                            delta=0.01,
                            nfree=2,
                        )

        self.assertEqual(
            captured['indices'],
            [
                1,
            ],
        )

    def testAseMetadataRejectsNonzeroExcludedAtomRows(self):
        ase_data = {
            'metadata': {
                'schema_version': 1,
                'number_of_atoms': 2,
                'vibrated_atom_indices': [
                    1,
                ],
                'expected_mode_count': 3,
            },
            'modes': [
                {
                    'mode': 1,
                    'frequency_cm-1': 100.0,
                    'eigenvector': [
                        # Atom 0 is excluded, but this row
                        # incorrectly contains a displacement.
                        [0.001, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                    ],
                },
                {
                    'mode': 2,
                    'frequency_cm-1': 50.0,
                    'eigenvector': [
                        [0.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ],
                },
                {
                    'mode': 3,
                    'frequency_cm-1': 25.0,
                    'eigenvector': [
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 1.0],
                    ],
                },
            ],
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            filename = os.path.join(
                tmpdir,
                'ase_nonzero_excluded_atom.json',
            )

            with open(
                filename,
                'w',
            ) as ase_fh:
                json.dump(
                    ase_data,
                    ase_fh,
                    indent=2,
                )

            parser = vasp_raman.build_parser()

            options = parser.parse_args([
                '--ase-modes',
                filename,
            ])

            with self.assertRaisesRegex(
                ValueError,
                'excluded atom 0',
            ):
                vasp_raman.load_modes(
                    options=options,
                    nat=2,
                    frozen=[
                        True,
                        False,
                    ],
                    poscar_symbols=[
                        'Si',
                        'O',
                    ],
                    direct_ase_settings=None,
                )

    def testValidateAseCacheReuseRejectsMalformedNonemptyCache(self):
        expected_manifest = {
            'schema_version': 1,
            'structure_sha256': 'structure-a',
            'structure': {
                'atomic_numbers': [
                    1,
                ],
            },
            'calculator_settings': {
                'encut': 500,
            },
            'delta_angstrom': 0.01,
            'nfree': 2,
            'vibrated_atom_indices': [
                0,
            ],
        }

        invalid_cache_contents = {
            'truncated-json': '{',
            'not-an-object': '[]',
            'missing-forces': '{}',
            'bad-force-shape': (
                '{"forces": [[0.0, 0.0]]}'
            ),
        }

        for label, contents in (
            invalid_cache_contents.items()
        ):
            with self.subTest(label=label):
                with tempfile.TemporaryDirectory() as tmpdir:
                    vibration_directory = os.path.join(
                        tmpdir,
                        'vib',
                    )

                    os.makedirs(
                        vibration_directory
                    )

                    cache_filename = os.path.join(
                        vibration_directory,
                        'cache.eq.json',
                    )

                    with open(
                        cache_filename,
                        'w',
                    ) as cache_fh:
                        cache_fh.write(
                            contents
                        )

                    manifest_filename = os.path.join(
                        tmpdir,
                        'ase_vibration_cache_manifest.json',
                    )

                    with open(
                        manifest_filename,
                        'w',
                    ) as manifest_fh:
                        json.dump(
                            expected_manifest,
                            manifest_fh,
                            indent=2,
                        )

                    with self.assertRaisesRegex(
                        ValueError,
                        'ASE vibration cache',
                    ):
                        vasp_raman.validate_ase_cache_reuse(
                            workdir=tmpdir,
                            expected_manifest=(
                                expected_manifest
                            ),
                        )