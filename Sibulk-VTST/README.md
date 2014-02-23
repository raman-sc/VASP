# WORK in PROGRESS

## Bulk silicon example using phonons from VTST tools

This example requires VTST tools version XXX or higher, download [here](http://theory.cm.utexas.edu/vtsttools/downloads/).

Unit cell parameters were taken from the [phonopy](http://phonopy.sourceforge.net) example folder. Silicon unit cell contains 8 atoms, thus 24 frequencies (phonons), three acoustic close to zero.

For this example `vasp_raman.py` script version [XXX] (https://raw.github.com/raman-sc/VASP/3cb3cdf0682609365c4b966472ef6eb5be1defc5/vasp_raman.py) was used.

Raman intensities for all modes will be computed in one go. For the multi-step procedure see [Cyclopentadiene example](https://github.com/raman-sc/VASP/tree/master/Cyclopentadiene).

Working directory should contain the following files (`FAT` file system has different soft linking ways):
```
INCAR        - should contain LEPSILON=.TRUE. or LCALCEPS=.TRUE. because we want 'MACROSCOPIC STATIC DIELECTRIC TENSOR' in the OUTCAR
KPOINTS      - just kpoints (soft link of the KPOINTS file from the PHONON folder, to prevent errors)
freq.dat     - contains frequencies (soft link of the freq.dat file from the PHONON folder)
modes_sqrt_amu.dat - contains eigenvectors/sqrt(mass in amu) (soft link of the modes_sqrt_amu.dat file from the PHONON folder)
POSCAR.phon  - VASP4 or VASP5 format is supported, at this point only positive scales are supported (soft link of the POSCAR file from the PHONON folder)
POTCAR       - `PAW_PBE Si 05Jan2001` PP (soft link to the POTCAR file from the PHONON folder)
raman.sub    - shell script for the scheduler to *do the deed*
```

### Obtaining phonons using vtst tools
1. Run `dymselsph.pl` in the `PHONON` folder to generate `DISPLACECAR`. Need to make sure that all atoms in the unit cell are found, in this case `8`:
```bash
perl ~/vtstscripts/dymselsph.pl POSCAR 8 3.87 0.001
----------------------------------------------------------------------
Using 1 central atom
Central Coordinates: 0.6250000000000000 0.6250000000000000 0.1250000000000000 
Central atom 1: 0.6250000000000000 0.6250000000000000 0.1250000000000000 
----------------------------------------------------------------------
8 atoms were found within a radius of 3.87 of atom 8, 
leading to 24 degrees of freedom selected.
----------------------------------------------------------------------
```

Be sure to check the `INCAR` file for the correct `vtst` settings ([http://theory.cm.utexas.edu/vtsttools/dynmat/](http://theory.cm.utexas.edu/vtsttools/dynmat/)):
```
! phonons
ICHAIN = 1      ! Run the dynamical matrix code
! IMAGES  M     ! Number of parallel images, if desired as in step 2 above; otherwise, do not add.
NSW = 25        ! (DOF/M)+1   Number of ionic steps
IBRION =   3    ! Tell VASP to run dynamics,
POTIM  = 0.0    ! with a time step of zero (ie, do nothing)
ISYM = 0        ! Make sure that symmetry is off,
EDIFFG = -1E-10 ! and that vasp does not quit due to low forces
```

Once phonons are obtained, we need to generate `freq.dat` and `modes_sqrt_amu.dat`, that are equivalent to the information obtained from the `OUTCAR` in the the case of phonons calculation using VASP:
```
python ~/vtstscripts/dymmatrix.py DISPLACECAR OUTCAR 
Reading DISPLACECAR
Number of displacements: 24
Reading OUTCAR
Building dynamical matrix
Diagonalizing matrix
    0.219483 cm^{-1} ... 1 
    0.023661 cm^{-1} ... 1 
    0.005978 cm^{-1} ... 1 
  145.585614 cm^{-1} ... 0 
  145.750931 cm^{-1} ... 0 
  145.783401 cm^{-1} ... 0 
  145.804204 cm^{-1} ... 0 
  145.823087 cm^{-1} ... 0 
  145.832016 cm^{-1} ... 0 
  402.135683 cm^{-1} ... 0 
  402.137959 cm^{-1} ... 0 
  402.145976 cm^{-1} ... 0 
  402.146508 cm^{-1} ... 0 
  402.150628 cm^{-1} ... 0 
  402.154078 cm^{-1} ... 0 
  447.963957 cm^{-1} ... 0 
  448.030348 cm^{-1} ... 0 
  448.143871 cm^{-1} ... 0 
  448.245366 cm^{-1} ... 0 
  448.402489 cm^{-1} ... 0 
  448.497029 cm^{-1} ... 0 
  504.414547 cm^{-1} ... 0 
  504.447782 cm^{-1} ... 0 
  504.483568 cm^{-1} ... 0 
```

Phonons look good, similar to those obtained from VASP ([Sibulk-VASP example](https://github.com/raman-sc/VASP/tree/master/Sibulk-VASP)).

Now we are ready to compute Raman intensities.

Contents of the `raman.sub`, note that modes are giving in the **increasing** order, contrary to VASP:
```bash
#!/bin/bash
#PBS -A ONRDC17403171
#PBS -l select=5:ncpus=32:mpiprocs=32
#PBS -l walltime=01:00:00
#PBS -q debug
#PBS -j oe
#PBS -N Si_bulk-Raman-VTST
#PBS -V

cd $PBS_O_WORKDIR

ulimit -s unlimited  # remove limit on stack size

export VASP_RAMAN_RUN='aprun -B /u/afonari/vasp.5.3.2/vasp.5.3/vasp &> job.out'
export VASP_RAMAN_PARAMS='04_24_2_0.01'

python /u/afonari/vasp_raman.py > vasp_raman.out

```

Submit all calculations to the scheduler:
```
qsub raman.sub
```

**Similar to the Sibulk-VASP example...**

After the job finished, a lot of `OUTCAR` files have been created (strictly speaking 42... *the answer??!*), `vasp_raman.dat` will contain Raman activities for the future processing (for example with `gnuplot`).

### Validation

From experiment, Raman spectrum of the bulk silicon contains only one intense peak at around 520 cm-1 ([J.H. Parker, *et al.*, Phys Rev, 155, 712 (1967)](http://dx.doi.org/10.1103/PhysRev.155.712)), `vasp_raman.dat` contains:
```
# mode    freq(cm-1)    alpha    beta2    activity
004   145.58561  -0.0001635   0.0000028   0.0000206
005   145.75093  -0.0000817   0.0000009   0.0000064
006   145.78340  -0.0000000   0.0000008   0.0000054
007   145.80420   0.0002044   0.0000003   0.0000039
008   145.82309   0.0000409   0.0000002   0.0000018
009   145.83202   0.0002861   0.0000011   0.0000114
010   402.13568  -0.0016350   0.0015959   0.0112915
011   402.13796  -0.0105863   0.0002777   0.0069868
012   402.14598  -0.0084200   0.0000903   0.0038223
013   402.14651  -0.0000817   0.0000081   0.0000572
014   402.15063   0.0021254   0.0011882   0.0085204
015   402.15408  -0.0058450   0.0007709   0.0069335
016   447.96396  -0.0002861   0.0001196   0.0008405
017   448.03035  -0.0018802   0.0008435   0.0060639
018   448.14387  -0.0018393   0.0039953   0.0281196
019   448.24537  -0.0024116   0.0011679   0.0084368
020   448.40249   0.0025751   0.0001548   0.0013818
021   448.49703  -0.0000817   0.0002679   0.0018757
022   504.41455  -0.0006540  780.5241570  5463.6691183
023   504.44778  -0.0008175  780.5880251  5464.1162060
024   504.48357  -0.0215814  779.9106481  5459.3954955
```
Triple degenerate mode at 504 cm-1 has the largest intensity.

[**Download complete example.**](https://github.com/raman-sc/VASP/raw/master/Sibulk-VTST/Si-VTST-pre0.6.1.tar.gz)

## Contributors

Shannon Stauffer (UT Austin, PI: G. Henkelman): [Email](mailto:stauffers[nospam]utexas.edu).  
Alexandr Fonari (Georgia Tech, PIs: J.-L. Bredas/V. Coropceanu): [Email](mailto:alexandr.fonari[nospam]gatech.edu)

