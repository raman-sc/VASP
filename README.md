# vasp_raman.py

Raman off-resonant activity calculator using VASP as a back-end.

## TOC

[Theory](#theory)  
[Installation](#installation)  
[Global variables](#global-variables)  
[External phonon inputs](#external-phonon-inputs)  
[Examples](#examples)  
[Changelog](#changelog)  
[How to cite](#how-to-cite)  
[Contributors](#contributors)

## Theory

In order to calculate off-resonance Raman activity of a mode, one needs to compute the derivative of the polarizability (or macroscopic dielectric tensor) with respect to that normal mode coordinate: `dP/dQ` (or `de/dQ`).  
Thus, two ingredients are required:

 1. Phonons at Γ-point
 2. Macroscopic dielectric tensor

### Phonons at Γ-point
In VASP, phonons at Γ-point can be computed using either:

 * finite displacements: `IBRION=5` or `IBRION=6`; or
 * density functional perturbation theory (DFPT): `IBRION=7` or `IBRION=8`.

Only finite displacements are available when hybrid functional is employed.

### Macroscopic dielectric tensor
In VASP, macroscopic dielectric tensor can be computed using either:
 * DFPT: `LEPSILON=.TRUE.`
 * or from frequency dependent dielectric matrix calculation: `LOPTICS=.TRUE.`.

In the latter case, hybrids functionals could be employed.  
For a more formal description of the method see [D. Porezag, M.R. Pederson, PRB, 54, 7830 (1996)](http://doi.org/10.1103/PhysRevB.54.7830).

## Installation

Python 3.6 or later is required. Copy `vasp_raman.py` into a directory in `$PATH`, or run it directly from the repository.

The original `OUTCAR.phon` and `freq.dat` / `modes_sqrt_amu.dat` routes use only the Python standard library.

Additional requirements for the external phonon options are:

- `--phonopy-qpoints` requires [PyYAML](https://pyyaml.org/);
- `--ase-modes` reads a JSON file and does not require ASE at conversion time;
- ASE is required if the vibration modes are generated using `ase.vibrations.Vibrations`.

PyYAML and ASE can be installed using:

```bash
python3 -m pip install pyyaml ase
```

## Global variables

For Raman calculations, `vasp_raman.py` requires two environment variables to be set:

  - `VASP_RAMAN_PARAMS` is defined as `FIRST_LAST_NDERIV_STEP` where:
      - `FIRST` - integer, first mode for which derivative of the polarizability is computed
      - `LAST`  - integer, last mode for which derivative of the polarizability is computed
      - `NDERIV`     - integer, scheme for finite difference, **currently** only value `2` is supported
      - `STEP`   - float, step-size for finite difference, in Angstroms
        
    Example: `VASP_RAMAN_PARAMS=01_10_2_0.01`

  - `VASP_RAMAN_RUN` the command to execute VASP (can contain MPI call):  
Example: `VASP_RAMAN_RUN='aprun -B /u/afonari/vasp.5.3.2/vasp.5.3/vasp &> job.out'`

Both variables should be `exported` (in Bash language) before running `vasp_raman.py`.
The environment variables are not required when using
`--write-external-modes` to convert Phonopy or ASE mode data without
running a Raman calculation.

An example of PBS script:

```bash
#!/bin/bash
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -l walltime=01:00:00
#PBS -q debug
#PBS -j oe
#PBS -N Example
#PBS -V

cd $PBS_O_WORKDIR

ulimit -s unlimited  # remove limit on stack size

export VASP_RAMAN_RUN='aprun -B /u/afonari/vasp.5.3.2/vasp.5.3/vasp &> job.out'
export VASP_RAMAN_PARAMS='01_10_2_0.01'

python3 vasp_raman.py > vasp_raman.out
```

An example of bash script (in case no scheduler is installed):

```bash
#!/bin/bash

# suggested by Ricardo Faccio, Universidad de la República, Montevideo, Uruguay

# OpenMP variables
#export OMP_NUM_THREADS=1
#export MKL_NUM_THREADS=1

# vasp_raman.py variables
export VASP_RAMAN_RUN='mpirun -np 4 vasp5.3.5_par'
export VASP_RAMAN_PARAMS='01_06_2_0.01'

python3 /home/user/bin/vasp_raman.py > vasp_raman.out

```
## External phonon inputs

`vasp_raman.py` can read Γ-point vibrational modes from several sources:

1. a Phonopy `qpoints.yaml` file;
2. an ASE-generated JSON mode file;
3. a direct ASE `Vibrations` calculation using VASP;
4. existing `freq.dat` and `modes_sqrt_amu.dat` files;
5. the original VASP `OUTCAR.phon` file.

Only one external source should be selected at a time.



### General requirements

The structure supplied with `--phonon-poscar` must have the same:

- number of atoms;
- atom ordering;
- chemical symbols;
- Cartesian coordinate convention

as the external phonon calculation.

For a structure containing `N` atoms, the external source must contain
exactly `3N` modes.

Unless `--keep-mode-order` is used, modes from Phonopy or ASE are sorted
from highest to lowest frequency.

### Phonopy input

The first q-point in the Phonopy file must be Γ:

```text
[0, 0, 0]
```

The Phonopy route:

- reads frequencies from `qpoints.yaml`;
- converts frequencies from THz to `cm^-1`;
- checks the atom symbols and atom order against the supplied POSCAR;
- converts Γ-point eigenvectors to a real representation;
- divides the eigenvectors by the square root of the atomic mass in amu.

A VASP 5 style POSCAR containing element symbols is required.

To convert Phonopy modes into the existing external-mode format:

```bash
python3 vasp_raman.py \
  --phonopy-qpoints qpoints.yaml \
  --phonon-poscar POSCAR.phon \
  --write-external-modes raman_inputs
```

This creates:

```text
raman_inputs/freq.dat
raman_inputs/modes_sqrt_amu.dat
raman_inputs/external_mode_order.csv
```

### ASE input

The ASE input is a JSON file containing the frequencies and Cartesian
eigenvectors.

It may be a list of modes:

```json
[
  {
    "mode": 1,
    "frequency_cm-1": 500.0,
    "eigenvector": [
      [0.1, 0.0, 0.0],
      [0.0, 0.1, 0.0]
    ]
  }
]
```

It may also be an object containing a `modes` list.

Each eigenvector must contain one `x y z` row per atom.

Modes returned by `ase.vibrations.Vibrations.get_mode()` are already
mass-scaled. They should be supplied directly and must not be divided
by the square root of the atomic mass again.

To convert ASE modes:

```bash
python3 vasp_raman.py \
  --ase-modes ase_gamma_modes.json \
  --phonon-poscar POSCAR.phon \
  --write-external-modes raman_inputs
```

### Direct ASE vibration calculation

The `--run-ase-vibrations` option runs or resumes an ASE
`Vibrations` calculation using the ASE VASP calculator. The generated
Gamma-point frequencies and mass-scaled Cartesian modes are passed
directly to the same shared validation and Raman workflow used by the
Phonopy and ASE JSON routes.

This option requires:

- ASE;
- a working VASP installation;
- VASP pseudopotentials;
- a JSON file containing the ASE VASP calculator settings;
- either a `command` entry in that JSON file or the
  `ASE_VASP_COMMAND` environment variable.

The reference structure is supplied using `--phonon-poscar`. The atom
order in this file is retained throughout the vibration calculation.

An example settings file is:

```json
{
  "encut": 500,
  "ediff": 1e-8,
  "xc": "PBE",
  "kpts": [1, 1, 1],
  "gamma": true,
  "ismear": 0,
  "sigma": 0.01,
  "ibrion": -1,
  "nsw": 0,
  "isym": 0,
  "prec": "Accurate",
  "lwave": false,
  "lcharg": false,
  "lreal": false
}
```

The VASP command may be supplied through the environment:

```bash
export ASE_VASP_COMMAND='mpirun -np 16 /path/to/vasp_gam'
```

Alternatively, the settings JSON may contain:

```json
{
  "command": "mpirun -np 16 /path/to/vasp_gam"
}
```

`VASP_PP_PATH` must point to the directory containing the VASP
pseudopotential collections, for example:

```bash
export VASP_PP_PATH=/path/to/vasp/potentials
```

A conversion-only Direct ASE run can be started using:

```bash
python3 vasp_raman.py \
  --run-ase-vibrations \
  --ase-vasp-settings ase_vasp_settings.json \
  --phonon-poscar POSCAR.phon \
  --ase-vibration-dir ase_vibrations \
  --write-external-modes raman_inputs
```

The default displacement size is `0.01 Angstrom` and the default
finite-difference order is `nfree=2`. These may be changed using:

```text
--ase-delta FLOAT
--ase-nfree INT
```

The Direct ASE route writes:

```text
ase_vibrations/ase_gamma_modes.json
ase_vibrations/ase_vibrations_summary.txt
ase_vibrations/ase_vibration_run.json
raman_inputs/freq.dat
raman_inputs/modes_sqrt_amu.dat
raman_inputs/external_mode_order.csv
```

ASE stores the force calculations in a vibration cache below the
selected vibration directory. A later invocation using the same
structure, displacement settings and vibration directory reuses all
completed cache entries. Only missing displacements are calculated.

A completed cache can therefore be used to resume an interrupted run.
Users should not manually edit or remove individual cache files unless
they intend to recalculate the corresponding displacement.

The Direct ASE route generates the vibration modes only. The Raman
activity calculation still uses VASP dielectric tensors as the
back-end.

### Generating displaced structures

External modes can be used directly with `--gen`.

For example:

```bash
export VASP_RAMAN_RUN='true'
export VASP_RAMAN_PARAMS='1_1_2_0.01'

python3 vasp_raman.py \
  --ase-modes ase_gamma_modes.json \
  --phonon-poscar POSCAR.phon \
  --gen
```

This creates:

```text
POSCAR.-1.out
POSCAR.+1.out
```

Use a persistent filename such as `POSCAR.phon` for the reference
structure. Do not use `POSCAR` as the persistent input filename because
`vasp_raman.py` uses `POSCAR` internally as a temporary displaced
structure.

### Running a Raman calculation directly

The external source can also be used without first writing
`freq.dat` and `modes_sqrt_amu.dat`.

For ASE modes:

```bash
export VASP_RAMAN_RUN='mpirun -np 16 vasp_std >> job.out 2>&1'
export VASP_RAMAN_PARAMS='01_21_2_0.01'

python3 vasp_raman.py \
  --ase-modes ase_gamma_modes.json \
  --phonon-poscar POSCAR.phon
```

For Phonopy modes:

```bash
python3 vasp_raman.py \
  --phonopy-qpoints qpoints.yaml \
  --phonon-poscar POSCAR.phon
```

### New command-line options

```text
--phonopy-qpoints FILE
    Read Γ-point modes from a Phonopy qpoints.yaml file.

--ase-modes FILE
    Read mass-scaled modes from an ASE-generated JSON file.

--phonon-poscar FILE
    Read the structure corresponding to the external modes.
    The default is POSCAR.phon.

--keep-mode-order
    Keep the ordering supplied by the external source instead
    of sorting the modes from highest to lowest frequency.

--write-external-modes DIR
    Write freq.dat, modes_sqrt_amu.dat and
    external_mode_order.csv to DIR, then exit.

--run-ase-vibrations
    Run or resume an ASE Vibrations calculation using VASP and use
    the resulting Gamma-point modes.

--ase-vasp-settings FILE
    Read ASE Vasp calculator parameters from a JSON file. This option
    is required with --run-ase-vibrations.

--ase-vibration-dir DIR
    Store ASE calculator directories, vibration cache files, mode JSON,
    summary and run metadata in DIR. The default is ase_vibrations.

--ase-delta FLOAT
    Set the ASE finite-displacement size in Angstrom. The default is
    0.01 and the value must be positive.

--ase-nfree INT
    Set the ASE finite-difference order. Supported values are 2 and 4;
    the default is 2.
```

### Validation

The external-mode route checks:

- the number of atoms;
- the number of modes;
- the number of atomic rows in each eigenvector;
- the three-dimensional shape of each atomic vector;
- zero-norm modes;
- Phonopy Γ-point and atom-order consistency.

The Direct ASE route also checks:

- that `--ase-vasp-settings` is supplied;
- that the settings file can be read;
- that `--ase-delta` is positive;
- that `--ase-nfree` is either `2` or `4`;
- that a VASP command is available through the settings file or
  `ASE_VASP_COMMAND`;
- that only one external mode source is selected.

### Bulk-Si benchmark

The external input routes were tested using the same 8-atom bulk-Si cell
as the original example.

The main Raman-active modes were obtained at approximately:

```text
Original OUTCAR.phon route: 504.47 cm^-1
Phonopy-derived route:      503.80 cm^-1
ASE-derived route:          503.82 cm^-1
```

The dominant Raman activities from the Phonopy and ASE routes agreed
with the original workflow to approximately 1%.

### Direct ASE integration test

The Direct ASE route was also tested from an empty cache using a
three-atom water molecule.

With `delta=0.01 Angstrom` and `nfree=2`:

- all 19 equilibrium and displaced-force cache entries were generated;
- 9 vibrational modes were extracted;
- the three main internal modes were obtained near
  `3829.7`, `3718.4` and `1580.7 cm^-1`;
- a second invocation reused the completed cache;
- the cache hashes were unchanged after restart;
- `freq.dat`, `modes_sqrt_amu.dat` and
  `external_mode_order.csv` were identical between the first and
  restarted runs.

### Known limitations

- The Phonopy route currently supports real Γ-point eigenvectors.
  A mode is rejected if it remains complex after phase alignment.
- The ASE JSON input should use a documented convention for imaginary
  frequencies. If `frequency_cm-1` contains only the real part, an
  imaginary mode may otherwise be written as `0.0 cm^-1`.
- The external phonon calculation and `--phonon-poscar` structure must
  use exactly the same atom ordering.

## Examples

* [Raman activity spectrum for Si using VASP](https://github.com/raman-sc/VASP/tree/master/Sibulk-VASP)
* [Raman activity spectrum for Si using VTST tools](https://github.com/raman-sc/VASP/tree/master/Sibulk-VTST)
* [Raman activity spectrum for cyclopentadiene using VASP](https://github.com/raman-sc/VASP/tree/master/Cyclopentadiene)
* [Raman activity spectrum for Si using VTST tools and PW91 functional](https://github.com/raman-sc/VASP/tree/master/Sibulk)

## Changelog

#### Current
* Switched to Python 3.
* Added direct Γ-point mode input from Phonopy `qpoints.yaml`.
* Added direct input from ASE-generated mode JSON files.
* Added conversion-only output to `freq.dat` and `modes_sqrt_amu.dat`.
* Added external-mode validation and mode-order tracking.
* Added clearer progress and error messages.
* Added unit tests for the Phonopy and ASE input routes.
* Added direct ASE/VASP finite-displacement vibration calculations.
* Added restart and reuse of ASE vibration caches.
* Added ASE VASP settings, displacement-size and finite-difference CLI
  options.
* Replaced `optparse` with `argparse`.
* Moved command-line validation into `validate_options()`.
* Refactored structure loading, mode loading, displaced-structure
  generation, dielectric extraction and Raman calculation into
  reusable functions.
* Expanded the unit-test suite to 30 tests.
* Added real regressions for the original VASP, Phonopy, ASE JSON and
  Direct ASE routes.

#### [0.6.0](https://github.com/raman-sc/VASP/blob/0.6.0/vasp_raman.py)
* ADDED: ability to use phonons obtained from the [vtst tools](http://theory.cm.utexas.edu/vtsttools/dymmat.html)
* FIX: cleaned `POSCAR` parsing code

#### [0.5.1](https://raw.github.com/raman-sc/VASP/3cb3cdf0682609365c4b966472ef6eb5be1defc5/vasp_raman.py)
* FIX: contributors and version are now in the output
* FIX: [Cyclopentadiene example](https://github.com/raman-sc/VASP/tree/master/Cyclopentadiene) is now fully consistent with the version

#### [0.5](https://raw.github.com/raman-sc/VASP/3004f2fd455b0f81c28a2e227542b328d5998bbd/vasp_raman.py)
* Basic working functionality

## How to cite

Use [Bibtext](https://raw.githubusercontent.com/raman-sc/VASP/master/vasp_raman_py.bib) or [RIS](https://raw.githubusercontent.com/raman-sc/VASP/master/vasp_raman_py.ris) file for citation.

## Contributors

Alexandr Fonari (Georgia Tech, PIs: J.-L. Bredas/V. Coropceanu): [Email](mailto:alexandr.fonari[nospam]gatech.edu)  
Shannon Stauffer (UT Austin, PI: G. Henkelman): [Email](mailto:stauffers[nospam]utexas.edu).