# External phonon input routes for raman-sc/VASP

This folder contains helper scripts and notes for using external Γ-point phonon modes with `vasp_raman.py`.

The aim is to make it easier to use vibrational modes calculated outside the original VASP `OUTCAR.phon` route, while still using the existing `freq.dat` and `modes_sqrt_amu.dat` input format supported by `vasp_raman.py`.

## Motivation

The standard `raman-sc/VASP` workflow reads vibrational frequencies and eigenvectors from a VASP frequency calculation. This works well for small benchmark systems, but it can be useful to generate Γ-point phonons with external tools such as Phonopy or the Atomic Simulation Environment (ASE), especially when testing different finite-displacement workflows.

This workflow provides two alternative routes:

1. Phonopy → `freq.dat` + `modes_sqrt_amu.dat`
2. ASE Vibrations → `freq.dat` + `modes_sqrt_amu.dat`

These files can then be read by `vasp_raman.py` in the same way as the existing external-mode input route.

## Files

The main helper scripts are:

- `phonopy_to_raman_inputs.py`  
  Converts Phonopy Γ-point modes from `qpoints.yaml` into `freq.dat` and `modes_sqrt_amu.dat`.

- `ase_to_raman_inputs.py`  
  Converts ASE vibration mode output into `freq.dat` and `modes_sqrt_amu.dat`.

- `check_raman_inputs.py`  
  Performs basic checks on the generated external phonon input files, including atom number, mode number, and mode block structure.

## Input format

The generated files follow the external input format used by `vasp_raman.py`:

- `freq.dat`  
  One vibrational frequency per line, in `cm^-1`.

- `modes_sqrt_amu.dat`  
  Mode eigenvector components written in blocks, with one `x y z` row per atom and a blank line between modes.

The atom order must match the structure used by the Raman calculation, usually `POSCAR.phon`.

## Benchmark test

The workflow was tested on the bulk silicon example from `raman-sc/VASP`, using the same 8-atom cell as the original benchmark.

Three routes were compared:

1. Original `raman-sc/VASP` route using `OUTCAR.phon`
2. Phonopy-derived external mode route
3. ASE-derived external mode route

For the main Raman-active Si modes, the three routes gave very similar results:

| Mode | Original frequency / cm^-1 | Phonopy frequency / cm^-1 | ASE frequency / cm^-1 | Original activity | Phonopy activity | ASE activity |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 504.47552 | 503.79656 | 503.82165 | 5476.88 | 5444.85 | 5458.18 |
| 2 | 504.47464 | 503.79656 | 503.82152 | 5425.07 | 5469.15 | 5447.31 |
| 3 | 504.47201 | 503.79656 | 503.82143 | 5512.51 | 5446.40 | 5463.77 |

The main frequencies are within about `1 cm^-1` of the original workflow, and the dominant Raman activities agree well with the baseline result.

## Notes

- The external phonon modes must use the same atom order as `POSCAR.phon`.
- The generated mode files should be checked before running `vasp_raman.py`.
- In one ASE-derived full-mode test, VASP stopped due to a symmetry-related `RHOSYG` error. Adding `SYMPREC = 1E-4` to the INCAR allowed the calculation to complete.
- Large VASP output files such as `OUTCAR`, `vasprun.xml`, `WAVECAR`, `CHGCAR`, and `POTCAR` should not be included in this repository.

## Status

This is intended as a small helper workflow rather than a replacement for the existing `vasp_raman.py` workflow. The aim is to provide a reproducible way to test Phonopy- and ASE-derived Γ-point phonons with the existing external input route.
