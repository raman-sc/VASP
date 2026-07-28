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

import cmath
import json
import math
import os
import re
import sys


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
    # modified subroutine from phonopy 1.8.3 (New BSD license)
    #
    poscar_fh.seek(0)  # just in case
    lines = poscar_fh.readlines()
    #
    scale = float(lines[1])
    if scale < 0.0:
        print("[parse_poscar]: ERROR negative scale not implemented.")
        sys.exit(1)
    #
    b = []
    for i in range(2, 5):
        b.append([float(x) * scale for x in lines[i].split()[:3]])
    #
    vol = b[0][0]*b[1][1]*b[2][2] + b[1][0]*b[2][1]*b[0][2] + b[2][0]*b[0][1]*b[1][2] - \
          b[0][2]*b[1][1]*b[2][0] - b[2][1]*b[1][2]*b[0][0] - b[2][2]*b[0][1]*b[1][0]
    #
    try:
        num_atoms = [int(x) for x in lines[5].split()]
        line_at = 6
    except ValueError:
        symbols = [x for x in lines[5].split()]
        num_atoms = [int(x) for x in lines[6].split()]
        line_at = 7
    nat = sum(num_atoms)
    #
    if lines[line_at][0].lower() == 's':
        line_at += 1
    #
    if (lines[line_at][0].lower() == 'c' or lines[line_at][0].lower() == 'k'):
        is_scaled = False
    else:
        is_scaled = True
    #
    line_at += 1
    #
    positions, frozen = [], []
    for i in range(line_at, line_at + nat):
        tmp = [x.strip() for x in lines[i].split()]
        pos = [float(x) for x in tmp[:3]]
        #
        if is_scaled:
            pos = MAT_m_VEC(T(b), pos)
        #
        positions.append(pos)
        #
        is_frozen = len(tmp) == 6 and (tmp[3] == tmp[4] == tmp[5] == 'F')
        frozen.append(is_frozen)
        #
    #
    # will add title and 'Cartesian' later
    poscar_header = ''.join(lines[1:line_at - 1])
    return nat, vol, b, positions, frozen, poscar_header


def is_selective_dynamics(poscar_header):
    return poscar_header.splitlines()[-1].lower().startswith('s')


def parse_env_params(params):
    tmp = params.strip().split('_')
    if len(tmp) != 4:
        print(
            "[parse_env_params]: ERROR there should be exactly four parameters"
        )
        sys.exit(1)
    #
    [first, last, nderiv,
     step_size] = [int(tmp[0]),
                   int(tmp[1]),
                   int(tmp[2]),
                   float(tmp[3])]
    #
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


def parse_poscar_symbols(poscar_fh):
    """Return one chemical symbol per atom from a VASP 5 style POSCAR."""
    poscar_fh.seek(0)
    lines = poscar_fh.readlines()

    try:
        [int(x) for x in lines[5].split()]
    except ValueError:
        symbols = lines[5].split()
        counts = [int(x) for x in lines[6].split()]
    else:
        raise ValueError(
            "External phonon modes require a VASP 5 style POSCAR with "
            "element symbols."
        )

    if len(symbols) != len(counts):
        raise ValueError(
            "POSCAR element and atom-count lines are inconsistent."
        )

    expanded_symbols = []
    for symbol, count in zip(symbols, counts):
        expanded_symbols.extend([symbol] * count)

    return expanded_symbols


def order_external_modes(modes):
    """Return external modes ordered from highest to lowest frequency."""
    return sorted(
        modes,
        key=lambda mode: mode['frequency_cm-1'],
        reverse=True,
    )


def _real_phonopy_eigenvector(eigenvector, tolerance=1e-8):
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
                    complex(float(component[0]), float(component[1]))
                )
            else:
                row.append(complex(float(component), 0.0))

        if len(row) != 3:
            raise ValueError(
                "Each eigenvector row must have three components."
            )

        complex_values.append(row)

    flattened = [
        value
        for row in complex_values
        for value in row
    ]

    reference = max(flattened, key=abs, default=0.0)
    phase = (
        cmath.exp(-1j * cmath.phase(reference))
        if reference
        else 1.0
    )

    real_values = []

    for row in complex_values:
        rotated_row = [value * phase for value in row]
        max_imaginary = max(abs(value.imag) for value in rotated_row)

        if max_imaginary > tolerance:
            raise ValueError(
                "Phonopy mode remains complex after phase alignment. "
                "Only real Gamma-point eigenvectors are currently supported."
            )

        real_values.append([value.real for value in rotated_row])

    return real_values#

def _get_phonopy_points(data, expected_natom=None):
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

    points = _get_phonopy_points(
        data,
        expected_natom=natom_hint,
    )

    if points is None and metadata is not None:
        points = _get_phonopy_points(
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

    if len(bands) != 3 * nat:
        raise ValueError(
            f"Phonopy data contain {len(bands)} modes; "
            f"expected {3 * nat}."
        )

    modes = []

    for mode_index, band in enumerate(bands, start=1):
        eigenvector = _real_phonopy_eigenvector(
            band['eigenvector']
        )

        if len(eigenvector) != nat:
            raise ValueError(
                f"Phonopy mode {mode_index} contains "
                f"{len(eigenvector)} atoms; expected {nat}."
            )

        mass_scaled = []

        for atom_index, vector in enumerate(eigenvector):
            sqrt_mass = math.sqrt(masses[atom_index])
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

    if _get_phonopy_points(
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

                points = _get_phonopy_points(
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


def parse_ase_data(data, nat):
    """Parse mass-scaled modes from an ASE-generated JSON object."""
    if isinstance(data, dict):
        data = data.get('modes')

    if not isinstance(data, list) or not data:
        raise ValueError(
            "ASE mode data must be a non-empty list."
        )

    if len(data) != 3 * nat:
        raise ValueError(
            f"ASE data contain {len(data)} modes; "
            f"expected {3 * nat}."
        )

    modes = []

    for mode_index, item in enumerate(data, start=1):
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

        if len(eigenvector) != nat:
            raise ValueError(
                f"ASE mode {mode_index} contains "
                f"{len(eigenvector)} atoms; expected {nat}."
            )

        parsed_eigenvector = []

        for vector in eigenvector:
            if len(vector) != 3:
                raise ValueError(
                    "Each ASE eigenvector row must have "
                    "three components."
                )

            parsed_eigenvector.append([
                float(value)
                for value in vector
            ])

        # ASE Vibrations.get_mode() returns mass-scaled
        # displacement vectors.
        modes.append({
            'source_mode': int(
                item.get('mode', mode_index)
            ),
            'frequency_cm-1': float(frequency),
            'eigenvector': parsed_eigenvector,
        })

    return modes


def read_ase_modes(ase_fh, nat):
    """Read mass-scaled modes from an ASE JSON file handle."""
    ase_fh.seek(0)
    return parse_ase_data(
        json.load(ase_fh),
        nat,
    )



def validate_external_modes(modes, nat):
    """Validate atom and mode counts and calculate mode norms."""
    if len(modes) != 3 * nat:
        raise ValueError(
            f"External mode source contains {len(modes)} modes; "
            f"expected {3 * nat}."
        )

    eigvals = []
    eigvecs = []
    norms = []

    for mode_index, mode in enumerate(modes, start=1):
        eigenvector = mode['eigenvector']

        if len(eigenvector) != nat:
            raise ValueError(
                f"External mode {mode_index} contains "
                f"{len(eigenvector)} atoms; expected {nat}."
            )

        for vector in eigenvector:
            if len(vector) != 3:
                raise ValueError(
                    f"External mode {mode_index} contains "
                    "a non-3D vector."
                )

        norm = math.sqrt(
            sum(
                abs(value) ** 2
                for vector in eigenvector
                for value in vector
            )
        )

        if norm == 0.0:
            raise ValueError(
                f"External mode {mode_index} has zero norm."
            )

        eigvals.append(float(mode['frequency_cm-1']))
        eigvecs.append(eigenvector)
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


#
if __name__ == '__main__':
    from math import pi
    from shutil import move
    import os
    import datetime
    import optparse
    #
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
    description = "Before run, set environment variables:\n"
    description += "    VASP_RAMAN_RUN='mpirun vasp'\n"
    description += "    VASP_RAMAN_PARAMS='[first-mode]_[last-mode]_[nderiv]_[step-size]'\n\n"
    description += "bash one-liner is:\n"
    description += "VASP_RAMAN_RUN='mpirun vasp' VASP_RAMAN_PARAMS='1_2_2_0.01' python vasp_raman.py"
    #
    parser = optparse.OptionParser(description=description)
    parser.add_option('-g',
                      '--gen',
                      help='Generate POSCAR only',
                      action='store_true')
    parser.add_option(
        '-u',
        '--use_poscar',
        help='Use provided POSCAR in the folder, USE WITH CAUTION!!',
        action='store_true')
    parser.add_option(
        '--phonopy-qpoints',
        help='Read Gamma-point modes directly from Phonopy qpoints.yaml')
    parser.add_option(
        '--ase-modes',
        help='Read mass-scaled modes from an ASE-generated JSON file')
    parser.add_option(
        '--phonon-poscar',
        help='Structure matching the external modes (default: POSCAR.phon)',
        default='POSCAR.phon')
    parser.add_option(
        '--keep-mode-order',
        help='Keep the source mode order instead of sorting high-to-low',
        action='store_true')
    parser.add_option(
        '--write-external-modes',
        help='Write freq.dat and modes_sqrt_amu.dat to DIR and exit',
        metavar='DIR')
    (options, args) = parser.parse_args()
    args = vars(options)
    #
    external_sources = [
        args['phonopy_qpoints'],
        args['ase_modes'],
    ]

    if sum(
        source is not None
        for source in external_sources
    ) > 1:
        print(
            "[__main__]: ERROR Use only one of "
            "--phonopy-qpoints or --ase-modes"
        )
        sys.exit(1)

    if (
        args['write_external_modes']
        and not any(external_sources)
    ):
        print(
            "[__main__]: ERROR --write-external-modes requires "
            "--phonopy-qpoints or --ase-modes"
        )
        sys.exit(1)

    conversion_only = bool(
        args['write_external_modes']
    )
    #
    #
    VASP_RAMAN_RUN = os.environ.get(
        'VASP_RAMAN_RUN'
    )

    if not VASP_RAMAN_RUN and not conversion_only:
        print(
            "[__main__]: ERROR Set environment variable "
            "'VASP_RAMAN_RUN'"
        )
        print()
        parser.print_help()
        sys.exit(1)

    if VASP_RAMAN_RUN:
        print(
            "[__main__]: VASP_RAMAN_RUN='"
            + VASP_RAMAN_RUN
            + "'"
        )

    VASP_RAMAN_PARAMS = os.environ.get(
        'VASP_RAMAN_PARAMS'
    )

    if not VASP_RAMAN_PARAMS and not conversion_only:
        print(
            "[__main__]: ERROR Set environment variable "
            "'VASP_RAMAN_PARAMS'"
        )
        print()
        parser.print_help()
        sys.exit(1)

    if VASP_RAMAN_PARAMS:
        print(
            "[__main__]: VASP_RAMAN_PARAMS='"
            + VASP_RAMAN_PARAMS
            + "'"
        )

    if not conversion_only:
        first, last, nderiv, step_size = (
            parse_env_params(
                VASP_RAMAN_PARAMS
            )
        )

        assert first >= 1, (
            '[__main__]: First mode should be '
            'equal or larger than 1'
        )

        assert last >= first, (
            '[__main__]: Last mode should be '
            'equal or larger than first mode'
        )

        if args['gen']:
            assert last == first, (
                "[__main__]: '-gen' mode -> only generation "
                "for one mode makes sense"
            )

        assert nderiv == 2, (
            '[__main__]: At this time, nderiv = 2 '
            'is the only supported value'
        )

        disps = [-1, 1]
        coeffs = [-0.5, 0.5]
    #
    try:
        poscar_fh = open(
            args['phonon_poscar'],
            'r',
        )
    except IOError:
        print(
            "[__main__]: ERROR Couldn't open input file "
            + args['phonon_poscar']
            + ", exiting...\n"
        )
        sys.exit(1)
    #
    nat, vol, b, pos, frozen, poscar_header = parse_poscar(poscar_fh)
    poscar_symbols = None

    if args['phonopy_qpoints']:
        try:
            poscar_symbols = parse_poscar_symbols(
                poscar_fh
            )
        except ValueError as err:
            print(
                f"[external modes]: ERROR {err}"
            )
            sys.exit(1)

    poscar_fh.close()
    print(f'Total atoms {nat}.')
    if any(frozen):
        frozen_str = ', '.join(
            map(str,
                (idx for idx, is_frozen in enumerate(frozen) if is_frozen)))
        print(f'Frozen atoms (0th-based): {frozen_str}.')
    #
    #
    # Use external sources, VTST-format files,
    # or a VASP OUTCAR.phon.
    if (
        args['phonopy_qpoints']
        or args['ase_modes']
    ):
        try:
            if args['phonopy_qpoints']:
                print(
                    "[external modes]: Reading Phonopy "
                    "modes from "
                    + args['phonopy_qpoints']
                )

                with open(
                    args['phonopy_qpoints'],
                    'r',
                ) as external_fh:
                    modes = read_phonopy_qpoints(
                        external_fh,
                        poscar_symbols,
                    )
            else:
                print(
                    "[external modes]: Reading ASE modes "
                    "from "
                    + args['ase_modes']
                )

                with open(
                    args['ase_modes'],
                    'r',
                ) as external_fh:
                    modes = read_ase_modes(
                        external_fh,
                        nat,
                    )

            if not args['keep_mode_order']:
                modes = order_external_modes(
                    modes
                )

            eigvals, eigvecs, norms = (
                validate_external_modes(
                    modes,
                    nat,
                )
            )

        except (
            IOError,
            ImportError,
            KeyError,
            TypeError,
            ValueError,
        ) as err:
            print(
                f"[external modes]: ERROR {err}"
            )
            sys.exit(1)

        print(
            f"[external modes]: Loaded "
            f"{len(modes)} modes for {nat} atoms."
        )

        if conversion_only:
            (
                freq_path,
                modes_path,
                order_path,
            ) = write_external_mode_files(
                modes,
                args['write_external_modes'],
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

            sys.exit(0)

    elif (
        os.path.isfile('freq.dat')
        and os.path.isfile('modes_sqrt_amu.dat')
    ):
        try:
            freqdat_fh = open(
                'freq.dat',
                'r',
            )
        except IOError:
            print(
                "[__main__]: ERROR Couldn't open "
                "freq.dat, exiting...\n"
            )
            sys.exit(1)

        eigvals = parse_freqdat(
            freqdat_fh,
            nat,
        )
        freqdat_fh.close()

        try:
            modes_fh = open(
                'modes_sqrt_amu.dat',
                'r',
            )
        except IOError:
            print(
                "[__main__]: ERROR Couldn't open "
                "modes_sqrt_amu.dat, exiting...\n"
            )
            sys.exit(1)

        eigvecs, norms = parse_modesdat(
            modes_fh,
            nat,
        )
        modes_fh.close()

    elif os.path.isfile('OUTCAR.phon'):
        try:
            outcar_fh = open(
                'OUTCAR.phon',
                'r',
            )
        except IOError:
            print(
                "[__main__]: ERROR Couldn't open "
                "OUTCAR.phon, exiting...\n"
            )
            sys.exit(1)

        eigvals, eigvecs, norms = (
            get_modes_from_OUTCAR(
                outcar_fh,
                nat,
                frozen,
            )
        )
        outcar_fh.close()

    else:
        print(
            "[__main__]: Neither an external mode source, "
            "OUTCAR.phon, nor freq.dat/modes_sqrt_amu.dat "
            "was found, nothing to do, exiting..."
        )
        sys.exit(1)
    #
    output_fh = open('vasp_raman.dat', 'w')
    output_fh.write("# mode    freq(cm-1)    alpha    beta2    activity\n")
    for i in range(first - 1, last):
        eigval = eigvals[i]
        eigvec = eigvecs[i]
        norm = norms[i]
        #
        print()
        print("[__main__]: Mode #%i: frequency %10.7f cm-1; norm: %10.7f" %
              (i + 1, eigval, norm))
        #
        ra = [[0.0 for x in range(3)] for y in range(3)]
        for j in range(len(disps)):
            print(
                "[progress]: Mode %i/%i, displacement %+d (%i/%i)."
                % (
                    i + 1,
                    last,
                    disps[j],
                    j + 1,
                    len(disps),
                ),
                flush=True,
            )
            disp_filename = 'OUTCAR.%04d.%+d.out' % (i + 1, disps[j])
            #
            try:
                if args['gen']:
                    # Force IOError to skip to POSCAR generation
                    raise IOError

                outcar_fh = open(disp_filename, 'r')
                print(
                    "[progress]: Reusing existing "
                    + disp_filename
                    + " and reading its dielectric tensor.",
                    flush=True,
                )
            except IOError:
                if not args['use_poscar']:
                    print("[__main__]: File " + disp_filename +
                          " not found, preparing displaced POSCAR")
                    poscar_fh = open('POSCAR', 'w')
                    poscar_fh.write("%s %4.1e \n" % (disp_filename, step_size))
                    poscar_fh.write(poscar_header)
                    poscar_fh.write("Cartesian\n")
                    is_selective = is_selective_dynamics(poscar_header)
                    #
                    for k in range(nat):
                        is_frozen = frozen[k]
                        #
                        disp_val = [
                            eigvec[k][l] * step_size * disps[j] / norm
                            for l in range(3)
                        ]
                        if is_frozen:
                            disp_val = [0.0, 0.0, 0.0]
                        #
                        pos_disp = [pos[k][l] + disp_val[l] for l in range(3)]
                        #
                        poscar_fh.write(
                            '%15.10f %15.10f %15.10f' %
                            (pos_disp[0], pos_disp[1], pos_disp[2]))
                        #
                        if is_selective:
                            is_frozen_str = ' F F F' if is_frozen else ' T T T'
                        else:
                            is_frozen_str = ''
                        #
                        poscar_fh.write(f'{is_frozen_str}\n')
                    poscar_fh.close()
                else:
                    print("[__main__]: Using provided POSCAR")
                #
                if args['gen']:  # only generate POSCARs
                    poscar_fn = 'POSCAR.%+d.out' % disps[j]
                    move('POSCAR', poscar_fn)
                    print("[__main__]: '-gen' mode -> " + poscar_fn +
                          " with displaced atoms have been generated")
                    #
                    # last iteration for the current displacements list
                    if j + 1 == len(disps):
                        print(
                            "[__main__]: '-gen' mode -> POSCAR files with displaced atoms have been generated, exiting now"
                        )
                        sys.exit(0)

                    # Continue to the next displacement
                    continue

                else:  # run VASP here
                    print(
                        "[progress]: Mode %i/%i, displacement %+d: "
                        "running VASP..."
                        % (
                            i + 1,
                            last,
                            disps[j],
                        ),
                        flush=True,
                    )

                    return_code = os.system(
                        VASP_RAMAN_RUN
                    )

                    if return_code != 0:
                        print(
                            "[__main__]: ERROR VASP command returned "
                            "a non-zero status for mode %i, "
                            "displacement %+d. Check job.out."
                            % (
                                i + 1,
                                disps[j],
                            ),
                            flush=True,
                        )
                        sys.exit(1)
                    try:
                        move('OUTCAR', disp_filename)
                    except IOError as err:
                        print(
                            "[__main__]: ERROR VASP did not create "
                            "OUTCAR for mode %i, displacement %+d: %s. "
                            "Check the VASP command and job.out."
                            % (
                                i + 1,
                                disps[j],
                                err,
                            ),
                            flush=True,
                        )
                        sys.exit(1)
                    #
                    outcar_fh = open(disp_filename, 'r')
            #
            try:
                eps = get_epsilon_from_OUTCAR(outcar_fh)
                outcar_fh.close()
                print(
                    "[progress]: Mode %i/%i, displacement %+d: "
                    "dielectric tensor found in %s."
                    % (
                        i + 1,
                        last,
                        disps[j],
                        disp_filename,
                    ),
                    flush=True,
                )
            except Exception as err:
                print(
                    "[__main__]: ERROR while reading "
                    + disp_filename
                    + ": "
                    + str(err),
                    flush=True,
                )
                print(
                    "[__main__]: Moving "
                    + disp_filename
                    + " back to OUTCAR and exiting. "
                    "Check OUTCAR and job.out.",
                    flush=True,
                )
                move(
                    disp_filename,
                    'OUTCAR',
                )
                sys.exit(1)
            #
            for m in range(3):
                for n in range(3):
                    ra[m][n] += eps[m][n] * coeffs[
                        j] / step_size * norm * vol / (4.0 * pi)
            #units: A^2/amu^1/2 =         dimless   * 1/A         * 1/amu^1/2  * A^3
        #
        alpha = (ra[0][0] + ra[1][1] + ra[2][2]) / 3.0
        beta2 = ((ra[0][0] - ra[1][1])**2 + (ra[0][0] - ra[2][2])**2 +
                 (ra[1][1] - ra[2][2])**2 + 6.0 *
                 (ra[0][1]**2 + ra[0][2]**2 + ra[1][2]**2)) / 2.0
        print()
        print(
            "! %4i  freq: %10.5f  alpha: %10.7f  beta2: %10.7f  activity: %10.7f "
            % (i + 1, eigval, alpha, beta2, 45.0 * alpha**2 + 7.0 * beta2))
        output_fh.write(
            "%03i  %10.5f  %10.7f  %10.7f  %10.7f\n" %
            (i + 1, eigval, alpha, beta2, 45.0 * alpha**2 + 7.0 * beta2))
        output_fh.flush()
    #
    output_fh.close()