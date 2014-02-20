## Bulk silicon example using Dynamical Matrix Finite Difference Method (DMFDM)

Unit cell parameters were found from geometric optimization in VASP using PAW data sets with PW-91 xc functional.  Silicon unit cell contains 8 atoms, thus 24 frequencies (phonons), three acoustic close to zero.

For this example `vasp_raman.py` script version [0.6.0] (https://raw.github.com/raman-sc/VASP/51627f0c9ce24526ef757c1d8212cfa2dbd84648/vasp_raman.py) was used.

This example calculates the Raman activity tensor of bulk Si.  
The phonon modes are calculated by a finite difference method (http://theory.cm.utexas.edu/vtsttools/dynmat/) which is implemented through VASP.    

The Raman activity tensor is calculated following the completion of the dynamical matrix calculation using the script vasp_raman.py. as explained in the CPD example [Cyclopentadiene example](https://github.com/  raman-sc/VASP/tree/master/Cyclopentadiene).


Several proceedural changes are necessary to complete setup the calculation from DMFDM output as outlined below:

##Generating normal mode frequency files and eigenvector files
Necessary scripts from vtstscripts package:
dymmatrix.py

Run dymmatrix.py which creates the hessian matrix and gives output files: freq.dat, modes_sqrt_amu.   dat which will be used here, others output files are not relevant for this calcualtion.

Include freq.dat and modes_sqrt_amu.dat in your working directory

With the post-processing of the DM VASP calculation complete,

The order of these modes is from smallest to largest wavenumbers in the DMFDM freq.dat file


Create a directory which contains the following:
```
POSCAR.phon  -used in DM VASP run
OUTCAR       -from DM VASP run
INCAR.raman  -should contain LEPSILON=.TRUE. or LCALCEPS=.TRUE. because we want 'MACROSCOPIC STATIC  DIELECTRIC TENSOR' in the OUTCAR
KPOINTS
POTCAR
vasp_raman.py      -version 0.6.0
DMFDM.py           -support script to vasp_raman.py which parses DMFDM output
raman_sub.sh       -changed to link DMFDM output files
raman.sub   
collect.sh      
freq.dat
modes_sqrt_amu.dat
```
Run raman_sub.sh to setup and exectute vasp_raman.py through local queueing system

Run collect.sh after all VASP calculations are complete.
Contents of 'collect.sh':
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
The remaining steps are the same as in the CPD example [Cyclopentadiene example](https://github.com/raman-sc/VASP/tree/master/Cyclopentadiene).


### Results Show Wrong Raman Active Mode: To be fixed*

Resulting freq.dat output from DMFDM calculation with dr = 0.005 Angstoms, PW-91 xc, 4x4x4 kpt sampling:
```
#freq ... 1=imaginary, 0 = real
0.013841 cm^{-1} ... 1 
0.004294 cm^{-1} ... 1 
0.000438 cm^{-1} ... 1 
141.905626 cm^{-1} ... 0 
141.908828 cm^{-1} ... 0 
141.911775 cm^{-1} ... 0 
141.932114 cm^{-1} ... 0 
141.933816 cm^{-1} ... 0 
141.939454 cm^{-1} ... 0 
412.092481 cm^{-1} ... 0 
412.102361 cm^{-1} ... 0 
412.104556 cm^{-1} ... 0 
412.105947 cm^{-1} ... 0 
412.108592 cm^{-1} ... 0 
412.125218 cm^{-1} ... 0 
463.614422 cm^{-1} ... 0 
463.630034 cm^{-1} ... 0 
463.639223 cm^{-1} ... 0 
463.645767 cm^{-1} ... 0 
463.647886 cm^{-1} ... 0 
463.660010 cm^{-1} ... 0 
515.536841 cm^{-1} ... 0 
515.549200 cm^{-1} ... 0 
515.550550 cm^{-1} ... 0
```

output from 'vasp_raman.dat':
```
# mode    freq(cm-1)    alpha    beta2    activity
4   141.90563   0.0003618   5.6200775  39.3405484
5   141.90883   0.0052292  54.6599140  382.6206283
6   141.91178  -0.0037038  160.0769321  1120.5391422
7   141.93211   0.0000371  78.9567767  552.6974371
8   141.93382   0.0085383  83.6444174  585.5142023
9   141.93945   0.0270537  57.5587791  402.9443895
10   412.09248   0.0000720  171.8463528  1202.9244699
11   412.10236  -0.0020808  66.9794704  468.8564879
12   412.10456   0.0066282  25.5925725  179.1499847
13   412.10595  -0.0114676   4.1176669  28.8295864
14   412.10859  -0.0527437   9.8348916  68.9694266
15   412.12522  -0.0054311  18.4996014  129.4985374
16   463.61442   0.0069102  89.3339891  625.3400725
17   463.63003  -0.0030733  109.7296067  768.1076720
18   463.63922   0.0007015  135.9888766  951.9221585
19   463.64577  -0.0027378  62.7335161  439.1349496
20   463.64789  -0.0022364  38.0224668  266.1574928
21   463.66001  -0.0006996  79.5765437  557.0358282
22   515.53684  -0.0005787  21.8868936  153.2082706
23   515.54920   0.0390079  30.0167765  210.1859082
24   515.55055  -0.0580559  19.9332233  139.6842350
```
An intense peak is expected at 520 cm-1 from experiment and calculation from VASP-phonon calculation wiht PAW@PBE, dr = 0.01, 4x4x4 kpts shows triple degenerate mode at 504 cm-1 has the largest intensity.

From experiment, Raman spectrum of the bulk silicon contains only one intense peak at around 520 cm-1 ([J.H. Parker, *et al.*, Phys Rev, 155, 712 (1967)](http://dx.doi.org/10.1103/PhysRev.155.712))

[**Download complete example.**](https://github.com/raman-sc/VASP/blob/master/Sibulk/Sibulk-PW91-DMFDM.tar.gz)

## Contributors

Shannon Stauffer (UT Austin, PI: G. Henkelman): [Email](mailto:stauffers[nospam]utexas.edu).  
Alexandr Fonari (Georgia Tech, PIs: J.-L. Bredas/V. Coropceanu): [Email](mailto:alexandr.fonari[nospam]gatech.edu)

