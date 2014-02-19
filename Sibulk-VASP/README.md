## Bulk silicon example

Unit cell parameters were taken from the [phonopy](http://phonopy.sourceforge.net) example folder. Silicon unit cell contains 8 atoms, thus 24 frequencies (phonons), three acoustic close to zero.

For this example `vasp_raman.py` script version [0.5.1] (https://raw.github.com/raman-sc/VASP/3cb3cdf0682609365c4b966472ef6eb5be1defc5/vasp_raman.py) was used.

Raman intensities for all modes will be computed in one go. For the multi-step procedure see [Cyclopentadiene example](https://github.com/raman-sc/VASP/tree/master/Cyclopentadiene).

Working directory should contain the following files (`FAT` file system has different soft linking ways):
```
INCAR        - should contain LEPSILON=.TRUE. or LCALCEPS=.TRUE. because we want 'MACROSCOPIC STATIC DIELECTRIC TENSOR' in the OUTCAR
KPOINTS      - just kpoints (soft link of the KPOINTS file from the PHONON folder, to prevent errors)
OUTCAR.phon  - should contain 'Eigenvectors after division by SQRT(mass)' either from VASP (soft link of the OUTCAR file from the PHONON folder)
POSCAR.phon  - VASP5 format is *required* (atomic symbols and numbers), at this point only positive scales are supported (soft link of the POSCAR file from the PHONON folder)
POTCAR       - `PAW_PBE Si 05Jan2001` PP (soft link to the POTCAR file from the PHONON folder)
raman.sub    - shell script for the scheduler to *do the deed*
```

Contents of the `raman.sub`:
```bash
#!/bin/bash
#PBS -A ONRDC17403171
#PBS -l select=5:ncpus=32:mpiprocs=32
#PBS -l walltime=01:00:00
#PBS -q debug
#PBS -j oe
#PBS -N Si_bulk-R
#PBS -V

cd $PBS_O_WORKDIR

ulimit -s unlimited  # remove limit on stack size

export VASP_RAMAN_RUN='aprun -B /u/afonari/vasp.5.3.2/vasp.5.3/vasp &> job.out'
export VASP_RAMAN_PARAMS='01_21_2_0.01'

python /u/afonari/vasp_raman.py > vasp_raman.out

```

Submit all calculations to the scheduler:
```
qsub raman.sub
```

After the job finished, a lot of `OUTCAR` files have been created (strictly speaking 42... *the answer??!*), `vasp_raman.dat` will contain Raman activities for the future processing (for example with `gnuplot`).

### Validation

From experiment, Raman spectrum of the bulk silicon contains only one intense peak at around 520 cm-1 ([J.H. Parker, *et al.*, Phys Rev, 155, 712 (1967)](http://dx.doi.org/10.1103/PhysRev.155.712)), `vasp_raman.dat` contains:
```
# mode    freq(cm-1)    alpha    beta2    activity
1   504.47552   0.0000409  780.7952797  5465.5669580
2   504.47464   0.0045779  779.0264832  5453.1863253
3   504.47201  -0.0031882  779.4324411  5456.0275454
4   447.84208   0.0006131   0.0005437   0.0038231
5   447.84107  -0.0017167   0.0000277   0.0003268
6   447.84076   0.0049049   0.0000643   0.0015326
7   447.83928   0.0004905   0.0001417   0.0010030
8   447.83865  -0.0380539   0.0203246   0.2074370
9   447.82907  -0.0295112   0.0906046   0.6734232
10   402.16258  -0.0002044   0.0001276   0.0008950
11   402.16054   0.0001226   0.0000161   0.0001131
12   402.15924  -0.0008992   0.0000119   0.0001199
13   402.15921  -0.0012262   0.0000106   0.0001422
14   402.15746   0.0008584   0.0002334   0.0016672
15   402.15492   0.0005722   0.0000023   0.0000306
16   145.82303   0.0002044   0.0000004   0.0000045
17   145.82113  -0.0000817   0.0000001   0.0000010
18   145.81855  -0.0002044   0.0000002   0.0000029
19   145.81681  -0.0001226   0.0000007   0.0000054
20   145.81480   0.0000817   0.0000004   0.0000029
21   145.81429   0.0000409   0.0000002   0.0000011
```
Triple degenerate mode at 504 cm-1 has the largest intensity.

[**Download complete example.**](https://github.com/raman-sc/VASP/raw/master/Sibulk-VASP/Sibulk-VASP-vasp_raman-0.5.1.tar.gz)

## Contributors

Shannon Stauffer (UT Austin, PI: G. Henkelman): [Email](mailto:stauffers[nospam]utexas.edu).  
Alexandr Fonari (Georgia Tech, PIs: J.-L. Bredas/V. Coropceanu): [Email](mailto:alexandr.fonari[nospam]gatech.edu)

