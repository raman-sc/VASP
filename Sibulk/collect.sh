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


