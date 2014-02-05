## Cyclopentadiene example

Cyclopentadiene molecule has 11 atoms: C5H6. Thus, there are 33 modes (lowest three acoustic should be close to zero). Note that VASP prints the modes in the *decreasing* eigenvalue (frequency) order. In the current example, the lowest *four* modes are imaginary (probably due to the optimization issues), these will be skipped in the Raman calculation.

For this example `vasp_raman.py` script version [0.5.1] (https://raw.github.com/raman-sc/VASP/3cb3cdf0682609365c4b966472ef6eb5be1defc5/vasp_raman.py) was used.

The following procedure can be used to obtain Raman activities for the G-point phonons. `POSCAR.phon` and `OUTCAR.phon` should come from a phonon calculation. `INCAR` should have `NWRITE=3` variable set.

### Cyclopentadiene example: One calculation

TODO

### Cyclopentadiene example: Split modes

Create a directory with the following files:
```
INCAR.raman  - should contain LEPSILON=.TRUE. or LCALCEPS=.TRUE. because we want 'MACROSCOPIC STATIC DIELECTRIC TENSOR' in the OUTCAR
KPOINTS      - just kpoints
OUTCAR.phon  - should contain 'Eigenvectors after division by SQRT(mass)' either from VASP itself or from Henkelman's script
POSCAR.phon  - at this point only positive scales are supported.
POTCAR       - hard PAWs usually give better agreement with the experiment (`ATOMTYPE_h` in the potpaw folders)
raman.sub    - see below
raman_sub.sh - see below
```

Contents of the `raman_sub.sh`:
```bash
#!/bin/bash

nderiv_step_size='2_0.01'

export VASP_RAMAN_RUN='aprun -B /u/afonari/vasp.5.3.2/vasp.5.3/vasp &> job.out'

# manually split modes
Modes[0]='01_10'
Modes[1]='11_20'
Modes[2]='21_29'

for m in ${Modes[*]}
do
    VASP_RAMAN_PARAMS="${m}_${nderiv_step_size}"
    export VASP_RAMAN_PARAMS=${VASP_RAMAN_PARAMS}
    #
    FOLDER="Modes_${VASP_RAMAN_PARAMS}"
    #
    echo "Running VASP in ${FOLDER}"
    #
    if [ -d "${FOLDER}" ]; then
        cd "${FOLDER}"
        qsub raman.sub
        cd ..
    else
        mkdir "${FOLDER}"
        cd "${FOLDER}"
        ln -s ../OUTCAR.phon ./OUTCAR.phon
        ln -s ../POSCAR.phon ./POSCAR.phon
        ln -s ../POTCAR ./POTCAR
        ln -s ../INCAR.raman ./INCAR
        ln -s ../raman.sub ./raman.sub
        ln -s ../KPOINTS ./KPOINTS
        qsub raman.sub
        cd ..
    fi
done

```

Contents of the `raman.sub`:
```bash
#!/bin/bash
#PBS -l select=1:ncpus=32:mpiprocs=32
#PBS -l walltime=01:00:00
#PBS -q debug
#PBS -j oe
#PBS -N CPD-RT
#PBS -V

cd $PBS_O_WORKDIR

ulimit -s unlimited  # remove limit on stack size

python27 vasp_raman.py > vasp_raman.out

```

Submit all calculations using:
```
bash ./raman_sub.sh
```

Three folders will be created where calculations will run: `Modes_01_10_2_0.01`, `Modes_11_20_2_0.01` and `Modes_21_29_2_0.01`.

After all calculations are done, it may be a good idea to save all `OUTCAR.*` files from all the `Modes_*` folders for the future reference. Also, `vasp_raman.dat` will contain Raman activities for the future processing (for example with `broaden.py`).

This script will gather all `OUTCAR.*` files in a back up folder, and also concatenate all `vasp_raman.dat` files in one:
```bash
#!/bin/bash

backup_dir="Raman_Backup"
mkdir -p ${backup_dir}

rm -f vasp_raman.dat

for dir in ./Modes*
do
    cat "${dir}"/vasp_raman.dat >> vasp_raman.tmp
    cp "${dir}"/OUTCAR.* "${backup_dir}"
done

head -n 1 vasp_raman.tmp > vasp_raman.dat
sed '/^#/d' vasp_raman.tmp > vasp_raman.unsorted
sort -k 2 -n vasp_raman.unsorted >> vasp_raman.dat

rm -f vasp_raman.tmp
rm -f vasp_raman.unsorted

```

Finally, broaden and plot Raman spectrum (scripts are in the archive):

```bash
python ./broaden.py vasp_raman.dat
bash ./vasp_raman.gnuplot.sh
```

Enjoy your Raman spectrum in `Raman.ps`!

[**Download complete example.**](https://github.com/raman-sc/VASP/raw/master/Cyclopentadiene/Cyclopentadiene-vasp_raman-0.5.1.tar.gz)

## Contributors

Alexandr Fonari (Georgia Tech, PIs: J.-L. Bredas/V. Coropceanu): [Email](mailto:alexandr.fonari[nospam]gatech.edu)  
Shannon Stauffer (UT Austin, PI: G. Henkelman): [Email](mailto:stauffers[nospam]utexas.edu).
