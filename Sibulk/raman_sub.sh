#!/bin/bash

DIR=$(pwd)

nderiv_step_size='2_0.01'

export VASP_RAMAN_RUN='mpirun -n $NSLOTS vasp &> job.out'

# manually split modes
Modes[0]='10_24'
#Modes[1]='33_64'

#Modes[2]='65_96'

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
        
        #move OUTCAR.phon if present, freq.dat and modes* if not
        if [ -a "$DIR/OUTCAR.phon" ]; then
#            echo "VASP CALC"
            ln -s ../OUTCAR.phon ./OUTCAR.phon
        elif [ -a "$DIR/freq.dat" ]; then
            ln -s $DIR/freq.dat 
            ln -s $DIR/modes_sqrt_amu.dat 
            ln -s $DIR/DMFDM.py 
 #           echo "DM CALC"
        fi

        ln -s ../POSCAR.phon ./POSCAR.phon
        ln -s ../POTCAR ./POTCAR
        ln -s ../INCAR.raman ./INCAR
        ln -s ../raman.sub ./raman.sub
        ln -s ../KPOINTS ./KPOINTS
        ln -s ../vasp_raman.py 
        qsub raman.sub
        cd ..
    fi
done
