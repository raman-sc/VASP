## Bulk silicon example

Unit cell parameters were taken from the phonopy example folder. Silicon unit cell contains 8 atoms, thus 24 frequencies (phonons), three acoustic close to zero.

For this example `vasp_raman.py` script version [0.5.1] (https://raw.github.com/raman-sc/VASP/3cb3cdf0682609365c4b966472ef6eb5be1defc5/vasp_raman.py) was used.

In this example, Raman intensities for all modes will be computed in one go. For the multi-step procedure see Cyclopentadiene example.

Working directory should contain the following files (`FAT` file system has different soft linking ways):
```
INCAR        - should contain LEPSILON=.TRUE. or LCALCEPS=.TRUE. because we want 'MACROSCOPIC STATIC DIELECTRIC TENSOR' in the OUTCAR
KPOINTS      - just kpoints (soft link of the KPOINTS file from the PHONON folder, to prevent errors)
OUTCAR.phon  - should contain 'Eigenvectors after division by SQRT(mass)' either from VASP (soft link of the OUTCAR file from the PHONON folder)
POSCAR.phon  - VASP5 format is *required* (atomic symbols and numbers), at this point only positive scales are supported (soft link of the POSCAR file from the PHONON folder)
POTCAR       - hard PAWs usually give better agreement with the experiment (`ATOMTYPE_h` in the potpaw folders) (soft link to the POTCAR file from the PHONON folder)
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

[**Download complete example.**](https://github.com/raman-sc/VASP/raw/master/Cyclopentadiene/Cyclopentadiene-vasp_raman-0.5.1.tar.gz)

## Contributors

Alexandr Fonari (Georgia Tech, PIs: J.-L. Bredas/V. Coropceanu): [Email](mailto:alexandr.fonari[nospam]gatech.edu)  
Shannon Stauffer (UT Austin, PI: G. Henkelman): [Email](mailto:stauffers[nospam]utexas.edu).
