# vasp_raman.py

Raman off-resonant activity calculator using VASP as a back-end.

## TOC

[Theory](#theory)  
[Installation](#installation)  
[Global variables](#global-variables)  
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

Python >= 2.6 is required. Just copy `vasp_raman.py` in the `$PATH` and run! No external dependencies.

## Global variables

`vasp_raman.py` **requires** two environmental variables to be set:

  - `VASP_RAMAN_PARAMS` is defined as `FIRST-MODE_LAST-MODE_NDERIV_STEPSIZE` where:
      - `FIRST_MODE` - integer, first mode for which derivative of the polarizability is computed
      - `LAST-MODE`  - integer, last mode for which derivative of the polarizability is computed
      - `NDERIV`     - integer, scheme for finite difference, **currently** only value `2` is supported
      - `STEPSIZE`   - float, step-size for finite difference, in Angstroms
        
    Example: `VASP_RAMAN_PARAMS=01_10_2_0.01`

  - `VASP_RAMAN_RUN` the command to execute VASP (can contain MPI call):  
Example: `VASP_RAMAN_RUN='aprun -B /u/afonari/vasp.5.3.2/vasp.5.3/vasp &> job.out'`

Both variables should be `exported` (in Bash language) before running `vasp_raman.py`.

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

python27 vasp_raman.py > vasp_raman.out
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

python /home/user/bin/vasp_raman.py > vasp_raman.out

```

## Examples

* [Raman activity spectrum for Si using VASP](https://github.com/raman-sc/VASP/tree/master/Sibulk-VASP)
* [Raman activity spectrum for Si using VTST tools](https://github.com/raman-sc/VASP/tree/master/Sibulk-VTST)
* [Raman activity spectrum for cyclopentadiene using VASP](https://github.com/raman-sc/VASP/tree/master/Cyclopentadiene)
* [Raman activity spectrum for Si using VTST tools and PW91 functional](https://github.com/raman-sc/VASP/tree/master/Sibulk)

## Changelog

#### 0.6 (will be released soon)
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
