#!/usr/bin/env python
#
# vasp_raman.py v. 0.6.0
#
# Raman off-resonant activity calculator
# using VASP as a back-end.
#
# Contributors: Alexandr Fonari (Georgia Tech)
# Shannon Stauffer (UT Austin)
#
# URL: http://raman-sc.github.io
#
# MIT license, 2013 - 2016
#

import argparse
import cmath
import datetime
import json
import math
import os
import re
import sys

from math import pi
from shutil import move

def MAT_m_VEC(m, v):
    p = [0.0 for i in range(len(v))]
    for i in range(len(m)):

        assert len(v) == len(
            m[i]
        ), 'Length of the matrix row is not equal to the length of the vector'

        p[i] = sum([m[i][j] * v[j] for j in range(len(v))])
    return p


def T(m):
    p = [[m[i][j] for i in range(len(m[j]))] for j in range(len(m))]
    return p


def parse_poscar(poscar_fh):
    """Parse a POSCAR and return structure and atom-order information."""
    poscar_fh.seek(0)
    lines = poscar_fh.readlines()

    if len(lines) < 8:
        raise ValueError(
            "POSCAR is incomplete: fewer than eight lines were found."
        )

    try:
        scale = float(lines[1])
    except ValueError as err:
        raise ValueError(
            "POSCAR scale factor is not a number."
        ) from err

    if scale <= 0.0:
        raise ValueError(
            "Zero or negative POSCAR scale factors are not supported."
        )

    lattice = []

    for line_index in range(2, 5):
        values = lines[line_index].split()

        if len(values) < 3:
            raise ValueError(
                f"POSCAR lattice line {line_index + 1} "
                "does not contain three components."
            )

        try:
            vector = [
                float(value) * scale
                for value in values[:3]
            ]
        except ValueError as err:
            raise ValueError(
                f"POSCAR lattice line {line_index + 1} "
                "contains a non-numeric value."
            ) from err

        if not all(math.isfinite(value) for value in vector):
            raise ValueError(
                f"POSCAR lattice line {line_index + 1} "
                "contains a non-finite value."
            )

        lattice.append(vector)

    volume = (
        lattice[0][0] * lattice[1][1] * lattice[2][2]
        + lattice[1][0] * lattice[2][1] * lattice[0][2]
        + lattice[2][0] * lattice[0][1] * lattice[1][2]
        - lattice[0][2] * lattice[1][1] * lattice[2][0]
        - lattice[2][1] * lattice[1][2] * lattice[0][0]
        - lattice[2][2] * lattice[0][1] * lattice[1][0]
    )

    atom_symbols = None

    try:
        counts = [
            int(value)
            for value in lines[5].split()
        ]
        line_at = 6
    except ValueError:
        element_symbols = lines[5].split()

        try:
            counts = [
                int(value)
                for value in lines[6].split()
            ]
        except ValueError as err:
            raise ValueError(
                "POSCAR atom-count line contains a non-integer value."
            ) from err

        if len(element_symbols) != len(counts):
            raise ValueError(
                "POSCAR element and atom-count lines are inconsistent."
            )

        atom_symbols = []

        for symbol, count in zip(element_symbols, counts):
            atom_symbols.extend([symbol] * count)

        line_at = 7

    if not counts or any(count < 0 for count in counts):
        raise ValueError(
            "POSCAR atom counts must be non-negative integers."
        )

    nat = sum(counts)

    if nat == 0:
        raise ValueError(
            "POSCAR does not contain any atoms."
        )

    if line_at >= len(lines):
        raise ValueError(
            "POSCAR is missing its coordinate section."
        )

    if lines[line_at].strip().lower().startswith('s'):
        line_at += 1

    if line_at >= len(lines):
        raise ValueError(
            "POSCAR is missing its coordinate type."
        )

    coordinate_type = lines[line_at].strip().lower()
    is_scaled = not (
        coordinate_type.startswith('c')
        or coordinate_type.startswith('k')
    )

    line_at += 1

    if len(lines) < line_at + nat:
        raise ValueError(
            f"POSCAR declares {nat} atoms but does not contain "
            f"{nat} coordinate lines."
        )

    positions = []
    frozen = []

    for atom_index in range(nat):
        values = lines[line_at + atom_index].split()

        if len(values) < 3:
            raise ValueError(
                f"POSCAR atom {atom_index + 1} coordinate line "
                "does not contain three components."
            )

        try:
            position = [
                float(value)
                for value in values[:3]
            ]
        except ValueError as err:
            raise ValueError(
                f"POSCAR atom {atom_index + 1} coordinate line "
                "contains a non-numeric value."
            ) from err

        if not all(math.isfinite(value) for value in position):
            raise ValueError(
                f"POSCAR atom {atom_index + 1} coordinate line "
                "contains a non-finite value."
            )

        if is_scaled:
            position = MAT_m_VEC(
                T(lattice),
                position,
            )

        positions.append(position)

        flags = [
            value.upper()
            for value in values[3:6]
        ]

        frozen.append(
            len(flags) == 3
            and flags == ['F', 'F', 'F']
        )

    poscar_header = ''.join(
        lines[1:line_at - 1]
    )

    return (
        nat,
        volume,
        lattice,
        positions,
        frozen,
        poscar_header,
        atom_symbols,
    )


def is_selective_dynamics(poscar_header):
    return poscar_header.splitlines()[-1].lower().startswith('s')


def parse_env_params(params):
    """Parse FIRST_LAST_NDERIV_STEP from VASP_RAMAN_PARAMS."""
    values = params.strip().split('_')

    if len(values) != 4:
        raise ValueError(
            "VASP_RAMAN_PARAMS must contain exactly four "
            "underscore-separated values."
        )

    try:
        first = int(values[0])
        last = int(values[1])
        nderiv = int(values[2])
        step_size = float(values[3])
    except ValueError as err:
        raise ValueError(
            "VASP_RAMAN_PARAMS contains an invalid numeric value."
        ) from err

    return first, last, nderiv, step_size


#### subs for the output from VTST tools
def parse_freqdat(freqdat_fh, nat):
    freqdat_fh.seek(0)  # just in case
    #
    eigvals = [0.0 for i in range(nat * 3)]
    #
    # all frequencies should be supplied, regardless of requested to calculate
    for i in range(nat * 3):
        tmp = freqdat_fh.readline().split()
        eigvals[i] = float(tmp[0])
    #
    return eigvals


#
def parse_modesdat(modesdat_fh, nat):
    from math import sqrt
    modesdat_fh.seek(0)  # just in case
    #
    eigvecs = [0.0 for i in range(nat * 3)]
    norms = [0.0 for i in range(nat * 3)]
    #
    # all frequencies should be supplied, regardless of requested to calculate
    for i in range(nat * 3):
        eigvec = []
        for j in range(nat):
            tmp = modesdat_fh.readline().split()
            eigvec.append([float(tmp[x]) for x in range(3)])
        #
        modesdat_fh.readline().split()  # empty line
        eigvecs[i] = eigvec
        norms[i] = sqrt(sum([abs(x)**2 for sublist in eigvec
                             for x in sublist]))
    #
    return eigvecs, norms


#### end subs for VTST
#
# External phonon mode helpers (Phonopy and ASE)
THZ_TO_CM = 33.35640951981521


def order_external_modes(modes):
    """Return external modes ordered from highest to lowest frequency."""
    return sorted(
        modes,
        key=lambda mode: mode['frequency_cm-1'],
        reverse=True,
    )


def real_phonopy_eigenvector(eigenvector, tolerance=1e-8):
    """Convert one Phonopy eigenvector to a real-valued representation."""
    complex_values = []

    for atom_vector in eigenvector:
        row = []

        for component in atom_vector:
            if isinstance(component, (list, tuple)):
                if len(component) != 2:
                    raise ValueError(
                        "Phonopy complex components must contain real and "
                        "imaginary parts."
                    )

                row.append(
                    complex(
                        float(component[0]),
                        float(component[1]),
                    )
                )
            else:
                row.append(
                    complex(
                        float(component),
                        0.0,
                    )
                )

        # Do not validate the vector length here.
        # The shared validator will check that it has three components.
        complex_values.append(row)

    flattened = [
        value
        for row in complex_values
        for value in row
    ]

    reference = max(
        flattened,
        key=abs,
        default=0.0,
    )

    phase = (
        cmath.exp(
            -1j * cmath.phase(reference)
        )
        if reference
        else 1.0
    )

    real_values = []

    for row in complex_values:
        rotated_row = [
            value * phase
            for value in row
        ]

        max_imaginary = max(
            (
                abs(value.imag)
                for value in rotated_row
            ),
            default=0.0,
        )

        if max_imaginary > tolerance:
            raise ValueError(
                "Phonopy mode remains complex after phase alignment. "
                "Only real Gamma-point eigenvectors are currently supported."
            )

        real_values.append([
            value.real
            for value in rotated_row
        ])

    return real_values

def get_phonopy_points(data, expected_natom=None):
    """Find a Phonopy points list matching the expected atom count."""
    if not isinstance(data, dict):
        return None

    candidates = []

    for section in [
        data,
        data.get('primitive_cell'),
        data.get('unit_cell'),
        data.get('supercell'),
    ]:
        if not isinstance(section, dict):
            continue

        points = section.get('points')

        if points:
            candidates.append(points)

    if expected_natom is not None:
        for points in candidates:
            if len(points) == expected_natom:
                return points

    if candidates:
        return candidates[0]

    return None

def parse_phonopy_data(
    data,
    expected_symbols=None,
    metadata=None,
):
    """Parse Gamma-point modes from an already-loaded Phonopy YAML object."""
    phonons = data.get('phonon')
    if not phonons:
        raise ValueError(
            "Phonopy data do not contain a 'phonon' section."
        )

    qpoint = phonons[0]
    q_position = qpoint.get(
        'q-position',
        qpoint.get('q_position', [0, 0, 0])
    )

    if any(abs(float(value)) > 1e-8 for value in q_position):
        raise ValueError(
            "The first Phonopy q-point is not Gamma (0, 0, 0)."
        )

    natom_hint = data.get('natom')

    if natom_hint is not None:
        natom_hint = int(natom_hint)

    points = get_phonopy_points(
        data,
        expected_natom=natom_hint,
    )

    if points is None and metadata is not None:
        points = get_phonopy_points(
            metadata,
            expected_natom=natom_hint,
        )

    if not points:
        raise ValueError(
            "Phonopy qpoints.yaml does not contain atom masses, "
            "and no suitable companion phonopy.yaml, "
            "phonopy_params.yaml or phonopy_disp.yaml file "
            "was found."
        )

    symbols = []
    masses = []

    for atom_index, point in enumerate(
        points,
        start=1,
    ):
        symbol = point.get('symbol')
        mass = point.get('mass')

        if symbol is None:
            raise ValueError(
                f"Phonopy atom {atom_index} has no symbol."
            )

        if mass is None:
            raise ValueError(
                f"Phonopy atom {atom_index} ({symbol}) "
                "has no mass in the companion metadata file."
            )

        symbols.append(symbol)
        masses.append(float(mass))

    nat = len(points)

    if (
        natom_hint is not None
        and nat != natom_hint
    ):
        raise ValueError(
            f"Phonopy qpoints.yaml contains {natom_hint} atoms, "
            f"but the companion metadata contains {nat} atoms."
        )

    if expected_symbols is not None and symbols != expected_symbols:
        raise ValueError(
            "Phonopy atom symbols/order do not match the supplied POSCAR."
        )

    bands = qpoint.get('band')
    if bands is None:
        raise ValueError(
            "Phonopy q-point data do not contain a 'band' section."
        )


    modes = []

    for mode_index, band in enumerate(bands, start=1):
        eigenvector = real_phonopy_eigenvector(
            band['eigenvector']
        )


        mass_scaled = []

        for atom_index, vector in enumerate(eigenvector):
            # Keep malformed extra rows so that the shared
            # validator can report the atom-count error.
            if atom_index < len(masses):
                sqrt_mass = math.sqrt(
                    masses[atom_index]
                )
            else:
                sqrt_mass = 1.0

            mass_scaled.append([
                component / sqrt_mass
                for component in vector
            ])

        modes.append({
            'source_mode': mode_index,
            'frequency_cm-1': (
                float(band['frequency']) * THZ_TO_CM
            ),
            'eigenvector': mass_scaled,
        })

    return modes


def read_phonopy_qpoints(
    qpoints_fh,
    expected_symbols=None,
):
    """Read modes and companion structure metadata from Phonopy YAML."""
    try:
        import yaml
    except ImportError as err:
        raise ImportError(
            "Reading Phonopy qpoints.yaml requires PyYAML."
        ) from err

    qpoints_fh.seek(0)
    data = yaml.safe_load(qpoints_fh)

    metadata = None
    natom_hint = data.get('natom')

    if natom_hint is not None:
        natom_hint = int(natom_hint)

    if get_phonopy_points(
        data,
        expected_natom=natom_hint,
    ) is None:
        qpoints_name = getattr(
            qpoints_fh,
            'name',
            None,
        )

        if qpoints_name:
            directory = os.path.dirname(
                os.path.abspath(qpoints_name)
            )

            for filename in [
                'phonopy.yaml',
                'phonopy_params.yaml',
                'phonopy_disp.yaml',
            ]:
                metadata_path = os.path.join(
                    directory,
                    filename,
                )

                if not os.path.isfile(metadata_path):
                    continue

                with open(
                    metadata_path,
                    'r',
                ) as metadata_fh:
                    candidate = yaml.safe_load(
                        metadata_fh
                    )

                points = get_phonopy_points(
                    candidate,
                    expected_natom=natom_hint,
                )

                if not points:
                    continue

                if not all(
                    point.get('mass') is not None
                    for point in points
                ):
                    continue

                metadata = candidate

                print(
                    "[external modes]: Reading Phonopy "
                    "atom symbols and masses from "
                    + metadata_path
                )

                break

    return parse_phonopy_data(
        data,
        expected_symbols=expected_symbols,
        metadata=metadata,
    )


def parse_ase_data(data):
    """Convert ASE JSON data into the canonical mode representation."""
    if isinstance(data, dict):
        data = data.get('modes')

    if not isinstance(data, list) or not data:
        raise ValueError(
            "ASE mode data must be a non-empty list."
        )

    modes = []

    for mode_index, item in enumerate(
        data,
        start=1,
    ):
        if not isinstance(item, dict):
            raise ValueError(
                f"ASE mode {mode_index} must be a JSON object."
            )

        frequency = item.get(
            'frequency_cm-1',
            item.get('frequency_cm_1'),
        )

        if frequency is None:
            raise ValueError(
                f"ASE mode {mode_index} has no frequency in cm-1."
            )

        eigenvector = item.get('eigenvector')

        if eigenvector is None:
            raise ValueError(
                f"ASE mode {mode_index} has no eigenvector."
            )

        # Only convert values here. Atom count and vector
        # dimensions are checked by validate_external_modes().
        parsed_eigenvector = [
            [
                float(value)
                for value in vector
            ]
            for vector in eigenvector
        ]

        modes.append({
            'source_mode': int(
                item.get('mode', mode_index)
            ),
            'frequency_cm-1': float(frequency),
            'eigenvector': parsed_eigenvector,
        })

    return modes


def read_ase_modes(ase_fh):
    """Read ASE modes from a JSON file handle."""
    ase_fh.seek(0)

    return parse_ase_data(
        json.load(ase_fh)
    )

def parse_ase_vasp_settings(data):
    """Validate ASE VASP calculator settings loaded from JSON."""
    if not isinstance(data, dict):
        raise ValueError(
            "ASE VASP settings must be a JSON object."
        )

    settings = dict(data)

    # ASE controls all atomic displacements. VASP must only
    # calculate the energy and forces for each fixed structure.
    if settings.get('ibrion', -1) != -1:
        raise ValueError(
            "ASE vibration calculations require IBRION = -1."
        )

    if settings.get('nsw', 0) != 0:
        raise ValueError(
            "ASE vibration calculations require NSW = 0."
        )

    settings.setdefault('ibrion', -1)
    settings.setdefault('nsw', 0)
    settings.setdefault('lwave', False)
    settings.setdefault('lcharg', False)

    # These are controlled separately by the direct ASE route.
    reserved_keys = {
        'atoms',
        'directory',
        'label',
        'restart',
        'txt',
    }

    invalid_keys = sorted(
        reserved_keys.intersection(settings)
    )

    if invalid_keys:
        raise ValueError(
            "The following ASE VASP settings are controlled "
            "by vasp_raman.py and must not appear in the JSON file: "
            + ", ".join(invalid_keys)
        )

    return settings


def read_ase_vasp_settings(settings_fh):
    """Read ASE VASP calculator settings from a JSON file handle."""
    settings_fh.seek(0)

    return parse_ase_vasp_settings(
        json.load(settings_fh)
    )

def signed_ase_frequency(frequency, tolerance=1.0e-8):
    """Return a real ASE frequency, using negative values for imaginary modes."""
    frequency = complex(frequency)

    if abs(frequency.imag) > tolerance:
        if abs(frequency.real) > tolerance:
            raise ValueError(
                "ASE frequency contains both significant real and "
                "imaginary parts."
            )

        return -abs(float(frequency.imag))

    return float(frequency.real)


def ase_vibrations_to_modes(vibrations):
    """Convert an ASE Vibrations result into the internal mode format."""
    try:
        import numpy as np
    except ImportError as err:
        raise ImportError(
            "Direct ASE vibration calculations require NumPy."
        ) from err

    frequencies = vibrations.get_frequencies()

    modes = []

    for mode_index, frequency in enumerate(
        frequencies,
        start=1,
    ):
        mode = np.asarray(
            vibrations.get_mode(mode_index - 1)
        )


        maximum_imaginary = float(
            np.max(np.abs(mode.imag))
        ) if np.iscomplexobj(mode) else 0.0

        if maximum_imaginary > 1.0e-8:
            raise ValueError(
                f"ASE mode {mode_index} contains a complex eigenvector."
            )

        modes.append({
            'source_mode': mode_index,
            'frequency_cm-1': signed_ase_frequency(
                frequency
            ),
            'eigenvector': mode.real.tolist(),
        })

    return modes

def run_ase_vibrations(
    structure_filename,
    settings,
    workdir='ase_vibrations',
    delta=0.01,
    nfree=2,
):
    """Run or resume ASE Vibrations with VASP and return the modes."""
    if delta <= 0.0:
        raise ValueError(
            "ASE vibration displacement delta must be positive."
        )

    if nfree not in (2, 4):
        raise ValueError(
            "ASE Vibrations supports nfree = 2 or nfree = 4."
        )

    try:
        import ase
        from ase.calculators.vasp import Vasp
        from ase.io import read
        from ase.vibrations import Vibrations
    except ImportError as err:
        raise ImportError(
            "Direct ASE vibration calculations require ASE."
        ) from err

    structure_filename = os.path.abspath(
        structure_filename
    )
    workdir = os.path.abspath(workdir)

    if not os.path.isfile(structure_filename):
        raise ValueError(
            "ASE vibration structure was not found: "
            + structure_filename
        )

    if not os.path.isdir(workdir):
        os.makedirs(workdir)

    calculator_directory = os.path.join(
        workdir,
        'vasp',
    )
    vibration_name = os.path.join(
        workdir,
        'vib',
    )

    print(
        "[direct ASE]: Reading structure from "
        + structure_filename,
        flush=True,
    )

    atoms = read(
        structure_filename,
        format='vasp',
    )

    print(
        f"[direct ASE]: Loaded {len(atoms)} atoms.",
        flush=True,
    )
    print(
        "[direct ASE]: VASP calculation directory: "
        + calculator_directory,
        flush=True,
    )
    print(
        "[direct ASE]: Vibration cache prefix: "
        + vibration_name,
        flush=True,
    )
    print(
        f"[direct ASE]: delta={delta} Angstrom, nfree={nfree}.",
        flush=True,
    )

    calculator = Vasp(
        directory=calculator_directory,
        txt='vasp.out',
        **settings,
    )

    atoms.calc = calculator

    vibrations = Vibrations(
        atoms,
        name=vibration_name,
        delta=delta,
        nfree=nfree,
    )

    print(
        "[direct ASE]: Running or resuming finite-displacement "
        "calculations...",
        flush=True,
    )

    vibrations.run()

    print(
        "[direct ASE]: Force calculations completed. "
        "Extracting frequencies and modes...",
        flush=True,
    )

    modes = ase_vibrations_to_modes(
        vibrations
    )

    summary_path = os.path.join(
        workdir,
        'ase_vibrations_summary.txt',
    )

    vibrations.summary(
        log=summary_path,
    )

    mode_json_path = os.path.join(
        workdir,
        'ase_gamma_modes.json',
    )

    json_modes = []

    for mode in modes:
        json_modes.append({
            'mode': mode['source_mode'],
            'frequency_cm-1': mode['frequency_cm-1'],
            'frequency_is_imaginary': (
                mode['frequency_cm-1'] < 0.0
            ),
            'eigenvector': mode['eigenvector'],
        })

    with open(mode_json_path, 'w') as mode_fh:
        json.dump(
            json_modes,
            mode_fh,
            indent=2,
        )

    run_log_path = os.path.join(
        workdir,
        'ase_vibration_run.json',
    )

    run_information = {
        'ase_version': ase.__version__,
        'structure': structure_filename,
        'number_of_atoms': len(atoms),
        'number_of_modes': len(modes),
        'delta_angstrom': delta,
        'nfree': nfree,
        'calculator_directory': calculator_directory,
        'vibration_cache': vibration_name,
        'calculator_settings': settings,
    }

    with open(run_log_path, 'w') as run_fh:
        json.dump(
            run_information,
            run_fh,
            indent=2,
        )

    print(
        "[direct ASE]: Wrote " + mode_json_path,
        flush=True,
    )
    print(
        "[direct ASE]: Wrote " + summary_path,
        flush=True,
    )
    print(
        "[direct ASE]: Wrote " + run_log_path,
        flush=True,
    )

    return modes, mode_json_path, summary_path

def validate_external_modes(modes, nat):
    """Validate canonical external modes and return Raman input arrays."""
    if not isinstance(nat, int) or nat <= 0:
        raise ValueError(
            "The number of atoms must be a positive integer."
        )

    if not isinstance(modes, list) or not modes:
        raise ValueError(
            "External mode source must contain a non-empty mode list."
        )

    if len(modes) != 3 * nat:
        raise ValueError(
            f"External mode source contains {len(modes)} modes; "
            f"expected {3 * nat}."
        )

    eigvals = []
    eigvecs = []
    norms = []

    for mode_index, mode in enumerate(
        modes,
        start=1,
    ):
        if not isinstance(mode, dict):
            raise ValueError(
                f"External mode {mode_index} is not a mode dictionary."
            )

        if 'frequency_cm-1' not in mode:
            raise ValueError(
                f"External mode {mode_index} has no frequency."
            )

        try:
            frequency = float(
                mode['frequency_cm-1']
            )
        except (TypeError, ValueError) as err:
            raise ValueError(
                f"External mode {mode_index} has an invalid frequency."
            ) from err

        if not math.isfinite(frequency):
            raise ValueError(
                f"External mode {mode_index} has a non-finite frequency."
            )

        eigenvector = mode.get('eigenvector')

        if not isinstance(
            eigenvector,
            (list, tuple),
        ):
            raise ValueError(
                f"External mode {mode_index} has no valid eigenvector."
            )

        if len(eigenvector) != nat:
            raise ValueError(
                f"External mode {mode_index} contains "
                f"{len(eigenvector)} atoms; expected {nat}."
            )

        validated_eigenvector = []

        for vector in eigenvector:
            if not isinstance(
                vector,
                (list, tuple),
            ) or len(vector) != 3:
                raise ValueError(
                    f"External mode {mode_index} contains "
                    "a non-3D vector."
                )

            try:
                values = [
                    float(value)
                    for value in vector
                ]
            except (TypeError, ValueError) as err:
                raise ValueError(
                    f"External mode {mode_index} contains "
                    "a non-numeric vector."
                ) from err

            if not all(
                math.isfinite(value)
                for value in values
            ):
                raise ValueError(
                    f"External mode {mode_index} contains "
                    "a non-finite vector."
                )

            validated_eigenvector.append(
                values
            )

        norm = math.sqrt(
            sum(
                value ** 2
                for vector in validated_eigenvector
                for value in vector
            )
        )

        if norm == 0.0:
            raise ValueError(
                f"External mode {mode_index} has zero norm."
            )

        eigvals.append(frequency)
        eigvecs.append(
            validated_eigenvector
        )
        norms.append(norm)

    return eigvals, eigvecs, norms


def write_external_mode_files(modes, output_dir='.'):
    """Write external modes in the format read by vasp_raman.py."""
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    freq_path = os.path.join(
        output_dir,
        'freq.dat',
    )
    modes_path = os.path.join(
        output_dir,
        'modes_sqrt_amu.dat',
    )
    order_path = os.path.join(
        output_dir,
        'external_mode_order.csv',
    )

    with open(freq_path, 'w') as freq_fh:
        for mode in modes:
            freq_fh.write(
                f"{mode['frequency_cm-1']:.10f}\n"
            )

    with open(modes_path, 'w') as modes_fh:
        for mode in modes:
            for vector in mode['eigenvector']:
                modes_fh.write(
                    f"{vector[0]: .12e} "
                    f"{vector[1]: .12e} "
                    f"{vector[2]: .12e}\n"
                )

            modes_fh.write('\n')

    with open(order_path, 'w') as order_fh:
        order_fh.write(
            'new_mode,source_mode,frequency_cm-1\n'
        )

        for new_index, mode in enumerate(
            modes,
            start=1,
        ):
            order_fh.write(
                f"{new_index},"
                f"{mode.get('source_mode', new_index)},"
                f"{mode['frequency_cm-1']:.10f}\n"
            )

    return freq_path, modes_path, order_path


def get_modes_from_OUTCAR(outcar_fh, nat, frozen):
    from math import sqrt
    #
    nfreq = 3 * (nat - sum(1 for x in frozen if x))
    print(f'Expected number of modes (excluding frozen atoms): {nfreq}.')
    #
    eigvals, eigvecs, norms = [], [], []
    #
    outcar_fh.seek(0)  # just in case
    while True:
        line = outcar_fh.readline()
        if not line:
            break
        #
        if "Eigenvectors after division by SQRT(mass)" in line:
            outcar_fh.readline()  # empty line
            outcar_fh.readline(
            )  # Eigenvectors and eigenvalues of the dynamical matrix
            outcar_fh.readline(
            )  # ----------------------------------------------------
            outcar_fh.readline()  # empty line
            #
            for i in range(nfreq):
                outcar_fh.readline()  # empty line
                p = re.search(r'^\s*(\d+).+?([\.\d]+) cm-1',
                              outcar_fh.readline())
                #
                eigvals.append(float(p.group(2)))
                #
                # X         Y         Z           dx          dy          dz
                outcar_fh.readline()
                #
                eigvec = []
                for j in range(nat):
                    tmp = outcar_fh.readline().split()
                    #
                    eigvec.append([float(tmp[x]) for x in range(3, 6)])
                    #
                eigvecs.append(eigvec)
                #
                norm = sqrt(
                    sum([abs(x)**2 for sublist in eigvec for x in sublist]))
                norms.append(norm)
            #
            assert len(eigvals) == len(eigvecs) == len(norms)
            #
            return eigvals, eigvecs, norms
        #
    print(
        "[get_modes_from_OUTCAR]: ERROR Couldn't find 'Eigenvectors after division by SQRT(mass)' in OUTCAR. Use 'NWRITE=3' in INCAR. Exiting..."
    )
    sys.exit(1)


def get_epsilon_from_OUTCAR(outcar_fh):
    epsilon = []
    #
    outcar_fh.seek(0)  # just in case
    while True:
        line = outcar_fh.readline()
        if not line:
            break
        #
        if "MACROSCOPIC STATIC DIELECTRIC TENSOR" in line:
            outcar_fh.readline()
            epsilon.append([float(x) for x in outcar_fh.readline().split()])
            epsilon.append([float(x) for x in outcar_fh.readline().split()])
            epsilon.append([float(x) for x in outcar_fh.readline().split()])
            return epsilon
    #
    raise RuntimeError(
        "[get_epsilon_from_OUTCAR]: ERROR Couldn't find dielectric tensor in OUTCAR"
    )
    return 1

def build_parser():
    """Create the command-line argument parser."""
    description = (
        "Raman off-resonant activity calculator using VASP "
        "as a back-end.\n\n"
        "For Raman calculations, set:\n"
        "  VASP_RAMAN_RUN='mpirun vasp'\n"
        "  VASP_RAMAN_PARAMS="
        "'FIRST_LAST_NDERIV_STEP'\n\n"
        "Example:\n"
        "  VASP_RAMAN_RUN='mpirun vasp' "
        "VASP_RAMAN_PARAMS='1_2_2_0.01' "
        "python3 vasp_raman.py"
    )

    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        '-g',
        '--gen',
        help='Generate displaced POSCAR files only.',
        action='store_true',
    )

    parser.add_argument(
        '-u',
        '--use_poscar',
        '--use-poscar',
        dest='use_poscar',
        help='Use an existing POSCAR in the working directory.',
        action='store_true',
    )

    parser.add_argument(
        '--phonopy-qpoints',
        help=(
            'Read Gamma-point modes directly from '
            'Phonopy qpoints.yaml.'
        ),
        metavar='FILE',
    )

    parser.add_argument(
        '--ase-modes',
        help=(
            'Read mass-scaled modes from an '
            'ASE-generated JSON file.'
        ),
        metavar='FILE',
    )

    parser.add_argument(
        '--phonon-poscar',
        help=(
            'Structure matching the phonon modes '
            '(default: POSCAR.phon).'
        ),
        default='POSCAR.phon',
        metavar='FILE',
    )

    parser.add_argument(
        '--keep-mode-order',
        help=(
            'Keep source mode order instead of '
            'sorting from high to low frequency.'
        ),
        action='store_true',
    )

    parser.add_argument(
        '--write-external-modes',
        help=(
            'Write freq.dat, modes_sqrt_amu.dat and '
            'external_mode_order.csv to DIR, then exit.'
        ),
        metavar='DIR',
    )

    parser.add_argument(
        '--run-ase-vibrations',
        help=(
            'Run or resume ASE Vibrations with VASP '
            'and use the resulting Gamma-point modes.'
        ),
        action='store_true',
    )

    parser.add_argument(
        '--ase-vasp-settings',
        help='JSON file containing ASE Vasp calculator settings.',
        metavar='FILE',
    )

    parser.add_argument(
        '--ase-vibration-dir',
        help=(
            'Directory for ASE calculations and '
            'the vibration cache.'
        ),
        default='ase_vibrations',
        metavar='DIR',
    )

    parser.add_argument(
        '--ase-delta',
        help='ASE finite-displacement size in Angstrom.',
        type=float,
        default=0.01,
        metavar='FLOAT',
    )

    parser.add_argument(
        '--ase-nfree',
        help='ASE finite-difference order: 2 or 4.',
        type=int,
        default=2,
        metavar='INT',
    )

    return parser

def validate_options(parser, args, environ=None):
    """Validate command-line options and environment variables."""
    if environ is None:
        environ = os.environ

    external_sources = [
        args.phonopy_qpoints,
        args.ase_modes,
        args.run_ase_vibrations,
    ]

    external_source_count = sum(
        bool(source)
        for source in external_sources
    )

    if external_source_count > 1:
        parser.error(
            "Use only one of --phonopy-qpoints, --ase-modes "
            "or --run-ase-vibrations."
        )

    if (
        args.run_ase_vibrations
        and not args.ase_vasp_settings
    ):
        parser.error(
            "--run-ase-vibrations requires "
            "--ase-vasp-settings FILE."
        )

    if (
        args.ase_vasp_settings
        and not args.run_ase_vibrations
    ):
        parser.error(
            "--ase-vasp-settings requires "
            "--run-ase-vibrations."
        )

    if args.ase_delta <= 0.0:
        parser.error(
            "--ase-delta must be positive."
        )

    if args.ase_nfree not in (2, 4):
        parser.error(
            "--ase-nfree must be 2 or 4."
        )

    if (
        args.write_external_modes
        and external_source_count == 0
    ):
        parser.error(
            "--write-external-modes requires "
            "--phonopy-qpoints, --ase-modes or "
            "--run-ase-vibrations."
        )

    if (
        args.keep_mode_order
        and external_source_count == 0
    ):
        parser.error(
            "--keep-mode-order requires an external mode source."
        )

    direct_ase_settings = None

    if args.run_ase_vibrations:
        if not os.path.isfile(
            args.ase_vasp_settings
        ):
            parser.error(
                "ASE VASP settings file was not found: "
                + args.ase_vasp_settings
            )

        try:
            with open(
                args.ase_vasp_settings,
                'r',
            ) as settings_fh:
                direct_ase_settings = (
                    read_ase_vasp_settings(
                        settings_fh
                    )
                )
        except (OSError, ValueError) as err:
            parser.error(
                "Could not read ASE VASP settings: "
                + str(err)
            )

        if (
            not direct_ase_settings.get('command')
            and not environ.get('ASE_VASP_COMMAND')
        ):
            parser.error(
                "Direct ASE vibrations require either a "
                "'command' entry in the settings JSON or "
                "the ASE_VASP_COMMAND environment variable."
            )

    conversion_only = bool(
        args.write_external_modes
    )

    vasp_raman_run = environ.get(
        'VASP_RAMAN_RUN'
    )
    vasp_raman_params = environ.get(
        'VASP_RAMAN_PARAMS'
    )

    if not conversion_only and not vasp_raman_run:
        parser.error(
            "Set environment variable VASP_RAMAN_RUN."
        )

    if not conversion_only and not vasp_raman_params:
        parser.error(
            "Set environment variable VASP_RAMAN_PARAMS."
        )

    first = None
    last = None
    nderiv = None
    step_size = None
    disps = None
    coeffs = None

    if not conversion_only:
        try:
            (
                first,
                last,
                nderiv,
                step_size,
            ) = parse_env_params(
                vasp_raman_params
            )
        except ValueError as err:
            parser.error(str(err))

        if first < 1:
            parser.error(
                "The first mode must be at least 1."
            )

        if last < first:
            parser.error(
                "The last mode must be equal to or "
                "larger than the first mode."
            )

        if args.gen and last != first:
            parser.error(
                "--gen supports exactly one mode."
            )

        if nderiv != 2:
            parser.error(
                "Only NDERIV = 2 is currently supported."
            )

        disps = [-1, 1]
        coeffs = [-0.5, 0.5]

    return {
        'conversion_only': conversion_only,
        'direct_ase_settings': direct_ase_settings,
        'vasp_raman_run': vasp_raman_run,
        'vasp_raman_params': vasp_raman_params,
        'first': first,
        'last': last,
        'nderiv': nderiv,
        'step_size': step_size,
        'disps': disps,
        'coeffs': coeffs,
    }

def load_structure(filename):
    """Open and parse the structure associated with the phonon modes."""
    try:
        with open(
            filename,
            'r',
        ) as poscar_fh:
            return parse_poscar(
                poscar_fh
            )

    except OSError as err:
        raise ValueError(
            "Could not open input structure "
            + filename
            + ": "
            + str(err)
        ) from err

    except ValueError as err:
        raise ValueError(
            "Could not parse input structure "
            + filename
            + ": "
            + str(err)
        ) from err

def load_modes(
    options,
    nat,
    frozen,
    poscar_symbols,
    direct_ase_settings,
):
    """Load phonon modes from the selected or available source."""
    external_modes = None

    if (
        options.phonopy_qpoints
        or options.ase_modes
        or options.run_ase_vibrations
    ):
        if options.phonopy_qpoints:
            if poscar_symbols is None:
                raise ValueError(
                    "Phonopy input requires a VASP 5 style "
                    "POSCAR containing element symbols."
                )

            print(
                "[external modes]: Reading Phonopy modes from "
                + options.phonopy_qpoints
            )

            with open(
                options.phonopy_qpoints,
                'r',
            ) as external_fh:
                external_modes = (
                    read_phonopy_qpoints(
                        external_fh,
                        poscar_symbols,
                    )
                )

        elif options.ase_modes:
            print(
                "[external modes]: Reading ASE modes from "
                + options.ase_modes
            )

            with open(
                options.ase_modes,
                'r',
            ) as external_fh:
                external_modes = read_ase_modes(
                    external_fh
                )

        else:
            print(
                "[direct ASE]: Starting or resuming "
                "ASE Vibrations.",
                flush=True,
            )

            (
                external_modes,
                ase_mode_json_path,
                ase_summary_path,
            ) = run_ase_vibrations(
                structure_filename=options.phonon_poscar,
                settings=direct_ase_settings,
                workdir=options.ase_vibration_dir,
                delta=options.ase_delta,
                nfree=options.ase_nfree,
            )

            print(
                "[direct ASE]: Mode JSON: "
                + ase_mode_json_path,
                flush=True,
            )
            print(
                "[direct ASE]: Vibration summary: "
                + ase_summary_path,
                flush=True,
            )

        if not options.keep_mode_order:
            external_modes = order_external_modes(
                external_modes
            )

        eigvals, eigvecs, norms = (
            validate_external_modes(
                external_modes,
                nat,
            )
        )

        print(
            f"[external modes]: Loaded "
            f"{len(external_modes)} modes for {nat} atoms."
        )

        return (
            eigvals,
            eigvecs,
            norms,
            external_modes,
        )

    if (
        os.path.isfile('freq.dat')
        and os.path.isfile('modes_sqrt_amu.dat')
    ):
        with open(
            'freq.dat',
            'r',
        ) as freqdat_fh:
            eigvals = parse_freqdat(
                freqdat_fh,
                nat,
            )

        with open(
            'modes_sqrt_amu.dat',
            'r',
        ) as modes_fh:
            eigvecs, norms = parse_modesdat(
                modes_fh,
                nat,
            )

        return (
            eigvals,
            eigvecs,
            norms,
            None,
        )

    if os.path.isfile('OUTCAR.phon'):
        with open(
            'OUTCAR.phon',
            'r',
        ) as outcar_fh:
            eigvals, eigvecs, norms = (
                get_modes_from_OUTCAR(
                    outcar_fh,
                    nat,
                    frozen,
                )
            )

        return (
            eigvals,
            eigvecs,
            norms,
            None,
        )

    raise ValueError(
        "Neither an external mode source, OUTCAR.phon, "
        "nor freq.dat/modes_sqrt_amu.dat was found."
    )

def write_displaced_poscar(
    output_filename,
    calculation_label,
    step_size,
    poscar_header,
    positions,
    frozen,
    eigenvector,
    norm,
    displacement,
):
    """Write one finite-displacement POSCAR for a normal mode."""
    nat = len(positions)

    if len(frozen) != nat:
        raise ValueError(
            "Frozen-atom flags do not match the number of positions."
        )

    if len(eigenvector) != nat:
        raise ValueError(
            "The mode eigenvector does not match the number of positions."
        )

    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError(
            "The mode norm must be a positive finite value."
        )

    is_selective = is_selective_dynamics(
        poscar_header
    )

    with open(
        output_filename,
        'w',
    ) as poscar_fh:
        poscar_fh.write(
            "%s %4.1e \n"
            % (
                calculation_label,
                step_size,
            )
        )
        poscar_fh.write(
            poscar_header
        )
        poscar_fh.write(
            "Cartesian\n"
        )

        for atom_index in range(nat):
            position = positions[atom_index]
            mode_vector = eigenvector[atom_index]

            if (
                len(position) != 3
                or len(mode_vector) != 3
            ):
                raise ValueError(
                    f"Atom {atom_index + 1} does not have "
                    "three Cartesian components."
                )

            displacement_vector = [
                mode_vector[component]
                * step_size
                * displacement
                / norm
                for component in range(3)
            ]

            if frozen[atom_index]:
                displacement_vector = [
                    0.0,
                    0.0,
                    0.0,
                ]

            displaced_position = [
                position[component]
                + displacement_vector[component]
                for component in range(3)
            ]

            poscar_fh.write(
                "%15.10f %15.10f %15.10f"
                % (
                    displaced_position[0],
                    displaced_position[1],
                    displaced_position[2],
                )
            )

            if is_selective:
                flags = (
                    " F F F"
                    if frozen[atom_index]
                    else " T T T"
                )
            else:
                flags = ""

            poscar_fh.write(
                flags + "\n"
            )

def get_displaced_epsilon(
    mode_number,
    displacement,
    step_size,
    poscar_header,
    positions,
    frozen,
    eigenvector,
    norm,
    vasp_command,
    generate_only=False,
    use_poscar=False,
):
    """Return the dielectric tensor for one displaced mode structure."""
    disp_filename = (
        f"OUTCAR.{mode_number:04d}."
        f"{displacement:+d}.out"
    )

    if (
        not generate_only
        and os.path.isfile(disp_filename)
    ):
        print(
            "[progress]: Reusing existing "
            + disp_filename
            + " and reading its dielectric tensor.",
            flush=True,
        )

    else:
        if not use_poscar:
            print(
                "[progress]: File "
                + disp_filename
                + " was not found; preparing displaced POSCAR.",
                flush=True,
            )

            write_displaced_poscar(
                output_filename='POSCAR',
                calculation_label=disp_filename,
                step_size=step_size,
                poscar_header=poscar_header,
                positions=positions,
                frozen=frozen,
                eigenvector=eigenvector,
                norm=norm,
                displacement=displacement,
            )

        elif not os.path.isfile('POSCAR'):
            raise ValueError(
                "--use-poscar was selected, but POSCAR "
                "was not found in the working directory."
            )

        else:
            print(
                "[progress]: Using provided POSCAR.",
                flush=True,
            )

        if generate_only:
            poscar_filename = (
                f"POSCAR.{displacement:+d}.out"
            )

            try:
                move(
                    'POSCAR',
                    poscar_filename,
                )
            except OSError as err:
                raise RuntimeError(
                    "Could not save generated structure as "
                    + poscar_filename
                    + ": "
                    + str(err)
                ) from err

            print(
                "[progress]: Generated "
                + poscar_filename,
                flush=True,
            )

            return None

        if not vasp_command:
            raise ValueError(
                "No VASP command was supplied."
            )

        print(
            "[progress]: Running VASP for mode "
            f"{mode_number}, displacement "
            f"{displacement:+d}.",
            flush=True,
        )

        return_code = os.system(
            vasp_command
        )

        if return_code != 0:
            raise RuntimeError(
                "VASP command returned a non-zero status "
                f"for mode {mode_number}, displacement "
                f"{displacement:+d}. Check job.out."
            )

        try:
            move(
                'OUTCAR',
                disp_filename,
            )
        except OSError as err:
            raise RuntimeError(
                "VASP did not create a usable OUTCAR "
                f"for mode {mode_number}, displacement "
                f"{displacement:+d}: {err}. "
                "Check the VASP command and job.out."
            ) from err

    try:
        with open(
            disp_filename,
            'r',
        ) as outcar_fh:
            epsilon = get_epsilon_from_OUTCAR(
                outcar_fh
            )

    except (
        OSError,
        RuntimeError,
        TypeError,
        ValueError,
    ) as err:
        # Preserve the previous workflow behaviour:
        # return an invalid displaced OUTCAR to the standard
        # OUTCAR filename for inspection.
        if os.path.isfile(
            disp_filename
        ):
            try:
                move(
                    disp_filename,
                    'OUTCAR',
                )
            except OSError:
                pass

        raise RuntimeError(
            "Could not read the dielectric tensor from "
            + disp_filename
            + ": "
            + str(err)
            + ". Check OUTCAR and job.out."
        ) from err

    print(
        "[progress]: Dielectric tensor found in "
        + disp_filename
        + ".",
        flush=True,
    )

    return epsilon

def calculate_raman_properties(raman_tensor):
    """Return alpha, beta squared and Raman activity."""
    if (
        not isinstance(raman_tensor, (list, tuple))
        or len(raman_tensor) != 3
    ):
        raise ValueError(
            "The Raman tensor must contain three rows."
        )

    validated_tensor = []

    for row in raman_tensor:
        if (
            not isinstance(row, (list, tuple))
            or len(row) != 3
        ):
            raise ValueError(
                "Each Raman tensor row must contain "
                "three components."
            )

        try:
            values = [
                float(value)
                for value in row
            ]
        except (TypeError, ValueError) as err:
            raise ValueError(
                "The Raman tensor contains a non-numeric value."
            ) from err

        if not all(
            math.isfinite(value)
            for value in values
        ):
            raise ValueError(
                "The Raman tensor contains a non-finite value."
            )

        validated_tensor.append(
            values
        )

    alpha = (
        validated_tensor[0][0]
        + validated_tensor[1][1]
        + validated_tensor[2][2]
    ) / 3.0

    beta2 = (
        (
            validated_tensor[0][0]
            - validated_tensor[1][1]
        ) ** 2
        + (
            validated_tensor[0][0]
            - validated_tensor[2][2]
        ) ** 2
        + (
            validated_tensor[1][1]
            - validated_tensor[2][2]
        ) ** 2
        + 6.0 * (
            validated_tensor[0][1] ** 2
            + validated_tensor[0][2] ** 2
            + validated_tensor[1][2] ** 2
        )
    ) / 2.0

    activity = (
        45.0 * alpha ** 2
        + 7.0 * beta2
    )

    return alpha, beta2, activity

def run_raman_calculation(
    options,
    configuration,
    volume,
    positions,
    frozen,
    poscar_header,
    eigvals,
    eigvecs,
    norms,
    output_filename='vasp_raman.dat',
):
    """Run finite-displacement Raman calculations for selected modes."""
    first = configuration['first']
    last = configuration['last']
    step_size = configuration['step_size']
    displacements = configuration['disps']
    coefficients = configuration['coeffs']
    vasp_command = configuration[
        'vasp_raman_run'
    ]

    if (
        first is None
        or last is None
        or step_size is None
    ):
        raise ValueError(
            "Raman calculation parameters are incomplete."
        )

    if (
        not math.isfinite(step_size)
        or step_size == 0.0
    ):
        raise ValueError(
            "The Raman finite-difference step must be "
            "a non-zero finite value."
        )

    if not (
        len(eigvals)
        == len(eigvecs)
        == len(norms)
    ):
        raise ValueError(
            "Frequency, eigenvector and norm arrays "
            "have different lengths."
        )

    if first < 1 or last < first:
        raise ValueError(
            "The requested Raman mode range is invalid."
        )

    if last > len(eigvals):
        raise ValueError(
            f"The requested last mode is {last}, but only "
            f"{len(eigvals)} modes are available."
        )

    if (
        not displacements
        or not coefficients
        or len(displacements) != len(coefficients)
    ):
        raise ValueError(
            "Finite-difference displacements and coefficients "
            "are inconsistent."
        )

    results = []

    with open(
        output_filename,
        'w',
    ) as output_fh:
        output_fh.write(
            "# mode    freq(cm-1)    alpha    beta2    activity\n"
        )

        for mode_index in range(
            first - 1,
            last,
        ):
            frequency = eigvals[mode_index]
            eigenvector = eigvecs[mode_index]
            norm = norms[mode_index]

            print()
            print(
                "[raman]: Mode #%i: frequency %10.7f cm-1; "
                "norm: %10.7f"
                % (
                    mode_index + 1,
                    frequency,
                    norm,
                )
            )

            raman_tensor = [
                [0.0 for column in range(3)]
                for row in range(3)
            ]

            for displacement_index, (
                displacement,
                coefficient,
            ) in enumerate(
                zip(
                    displacements,
                    coefficients,
                ),
                start=1,
            ):
                print(
                    "[progress]: Mode %i/%i, displacement "
                    "%+d (%i/%i)."
                    % (
                        mode_index + 1,
                        last,
                        displacement,
                        displacement_index,
                        len(displacements),
                    ),
                    flush=True,
                )

                epsilon = get_displaced_epsilon(
                    mode_number=mode_index + 1,
                    displacement=displacement,
                    step_size=step_size,
                    poscar_header=poscar_header,
                    positions=positions,
                    frozen=frozen,
                    eigenvector=eigenvector,
                    norm=norm,
                    vasp_command=vasp_command,
                    generate_only=options.gen,
                    use_poscar=options.use_poscar,
                )

                if epsilon is None:
                    continue

                if (
                    not isinstance(
                        epsilon,
                        (list, tuple),
                    )
                    or len(epsilon) != 3
                    or any(
                        not isinstance(row, (list, tuple))
                        or len(row) != 3
                        for row in epsilon
                    )
                ):
                    raise ValueError(
                        "The dielectric tensor must be 3 x 3."
                    )

                for row in range(3):
                    for column in range(3):
                        try:
                            epsilon_value = float(
                                epsilon[row][column]
                            )
                        except (
                            TypeError,
                            ValueError,
                        ) as err:
                            raise ValueError(
                                "The dielectric tensor contains "
                                "a non-numeric value."
                            ) from err

                        if not math.isfinite(
                            epsilon_value
                        ):
                            raise ValueError(
                                "The dielectric tensor contains "
                                "a non-finite value."
                            )

                        raman_tensor[row][column] += (
                            epsilon_value
                            * coefficient
                            / step_size
                            * norm
                            * volume
                            / (4.0 * pi)
                        )

            if options.gen:
                print(
                    "[raman]: '--gen' completed; positive "
                    "and negative displaced POSCAR files "
                    "were generated."
                )

                return results

            (
                alpha,
                beta2,
                activity,
            ) = calculate_raman_properties(
                raman_tensor
            )

            result = {
                'mode': mode_index + 1,
                'frequency_cm-1': frequency,
                'alpha': alpha,
                'beta2': beta2,
                'activity': activity,
            }

            results.append(
                result
            )

            print()
            print(
                "! %4i  freq: %10.5f  alpha: %10.7f  "
                "beta2: %10.7f  activity: %10.7f"
                % (
                    result['mode'],
                    result['frequency_cm-1'],
                    result['alpha'],
                    result['beta2'],
                    result['activity'],
                )
            )

            output_fh.write(
                "%03i  %10.5f  %10.7f  %10.7f  %10.7f\n"
                % (
                    result['mode'],
                    result['frequency_cm-1'],
                    result['alpha'],
                    result['beta2'],
                    result['activity'],
                )
            )
            output_fh.flush()

    return results

#
def main(argv=None):
    print()
    print("    Raman off-resonant activity calculator,")
    print("    using VASP as a back-end.")
    print("")
    print("    Contributors: Alexandr Fonari  (Georgia Tech)")
    print("                  Shannon Stauffer (UT Austin)")
    print("    MIT License, 2013")
    print("    URL: http://raman-sc.github.io")
    print("    Started at: " +
          datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
    print()
    #
    parser = build_parser()
    namespace = parser.parse_args(argv)
    configuration = validate_options(
        parser,
        namespace,
    )


    conversion_only = configuration[
        'conversion_only'
    ]
    direct_ase_settings = configuration[
        'direct_ase_settings'
    ]

    VASP_RAMAN_RUN = configuration[
        'vasp_raman_run'
    ]
    VASP_RAMAN_PARAMS = configuration[
        'vasp_raman_params'
    ]


    if VASP_RAMAN_RUN:
        print(
            "[__main__]: VASP_RAMAN_RUN='"
            + VASP_RAMAN_RUN
            + "'"
        )

    if VASP_RAMAN_PARAMS:
        print(
            "[__main__]: VASP_RAMAN_PARAMS='"
            + VASP_RAMAN_PARAMS
            + "'"
        )

    try:
        (
            nat,
            vol,
            b,
            pos,
            frozen,
            poscar_header,
            poscar_symbols,
        ) = load_structure(
            namespace.phonon_poscar
        )

    except ValueError as err:
        parser.error(
            str(err)
        )

    print(f'Total atoms {nat}.')
    if any(frozen):
        frozen_str = ', '.join(
            map(str,
                (idx for idx, is_frozen in enumerate(frozen) if is_frozen)))
        print(f'Frozen atoms (0th-based): {frozen_str}.')
    #
    try:
        (
            eigvals,
            eigvecs,
            norms,
            external_modes,
        ) = load_modes(
            namespace,
            nat,
            frozen,
            poscar_symbols,
            direct_ase_settings,
        )

    except (
        OSError,
        ImportError,
        KeyError,
        TypeError,
        ValueError,
        RuntimeError,
    ) as err:
        parser.error(
            str(err)
        )

    if conversion_only:
        (
            freq_path,
            modes_path,
            order_path,
        ) = write_external_mode_files(
            external_modes,
            namespace.write_external_modes,
        )

        print(
            f"[external modes]: Wrote {freq_path}"
        )
        print(
            f"[external modes]: Wrote {modes_path}"
        )
        print(
            f"[external modes]: Wrote {order_path}"
        )

        return 0

    try:
        run_raman_calculation(
            options=namespace,
            configuration=configuration,
            volume=vol,
            positions=pos,
            frozen=frozen,
            poscar_header=poscar_header,
            eigvals=eigvals,
            eigvecs=eigvecs,
            norms=norms,
        )

    except (
        OSError,
        TypeError,
        ValueError,
        RuntimeError,
    ) as err:
        parser.error(
            str(err)
        )

    return 0


if __name__ == '__main__':
    main()