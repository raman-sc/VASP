# vasp_raman.py

Raman off-resonant activity calculator using VASP as a back-end.

## Theory

In order to calculate off-resonance Raman activity of a mode, one needs to compute the derivative of the polarizability (or macroscopic dielectric tensor) with respect to that normal mode coordinate: `dP/dQ` (or `de/dQ`).  
Thus, two ingredients are required:

 1. Phonons at Γ-point
 2. Macroscopic dielectric tensor

### Phonons at Γ-point
In VASP, phonons at Γ-point can be computed either using either:

 * finite displacements: `IBRION=5` or `IBRION=6`; or
 * density functional perturbation theory (DFPT): `IBRION=7` or `IBRION=8`.

Only finite displacements are available when hybrid functional is employed.

### Macroscopic dielectric tensor
In VASP, macroscopic dielectric tensor can be computed using either:
 * DFPT: `LEPSILON=.TRUE.`;
 * or from frequency dependent dielectric matrix calculation: `LOPTICS=.TRUE.`.

In the latter case, hybrids functionals could be employed.

## Changelog

#### 0.6 (will be released soon)
* ADDED: ability to use phonons obtained from the [vtst tools](theory.cm.utexas.edu/vtsttools/dynmat/)
* FIX: cleaned `POSCAR` parsing code

#### [0.5.1](https://raw.github.com/raman-sc/VASP/3cb3cdf0682609365c4b966472ef6eb5be1defc5/vasp_raman.py)
* FIX: contributors and version are now in the output
* FIX: [Cyclopentadiene example](https://github.com/raman-sc/VASP/tree/master/Cyclopentadiene) is now fully consistent with the version

#### [0.5](https://raw.github.com/raman-sc/VASP/3004f2fd455b0f81c28a2e227542b328d5998bbd/vasp_raman.py)
* Basic working functionality

## Installation

Python >= 2.6 is required. Just copy `vasp_raman.py` in the `$PATH` and run! No external dependencies.

## Contributors

Alexandr Fonari (Georgia Tech, PIs: J.-L. Bredas/V. Coropceanu): [Email](mailto:alexandr.fonari[nospam]gatech.edu)  
Shannon Stauffer (UT Austin, PI: G. Henkelman): [Email](mailto:stauffers[nospam]utexas.edu).
