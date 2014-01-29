#!/usr/bin/env python
#PBS -l select=6:ncpus=32:mpiprocs=32
#PBS -l walltime=140:00:00
#PBS -q standard
#PBS -j oe
#PBS -N Cp

#PBS -M firstname.lastname@gatech.edu
#PBS -m abe

#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -N test-raman-sibulk-calc2
#$ -m es
#$ -V
#$ -pe mpi8 8 
#$ -o  $JOB_NAME.$JOB_ID
#$ -S /bin/bash


STEP_SIZE = 0.01 # Angstrom
NDERIV = 2

PBS=False
VASP_RAMAN_INP_FN='freqs'

import os
if PBS: #PBS enabled environment

    VASP_ENV = """
    export MPI_DSM_DISTRIBUTE=1
    """
    os.system(VASP_ENV)

    VASP_RUN = "aprun -B /u/afonari/vasp.5.3.2/vasp.5.3/vasp >& job.out"
    os.chdir(os.environ["PBS_O_WORKDIR"])

else:
    VASP_RUN = "mpirun -n $NSLOTS vasp"

############################

def parse_inpcar(inp_fh):
    import sys
    from math import sqrt

    poscar_header = ""
    vol = 0.0
    b = []
    atom_numbers = []

    [first, last] = [int(s) for s in inp_fh.readline().split()]

    scale = float(inp_fh.readline()) # scale
    for i in range(3): b.append( [float(s) for s in inp_fh.readline().split()] )

    if scale > 0.0:
        b = [[ b[i][j]*scale for i in range(3)] for j in range(3) ]
        scale = 1.0

        vol = b[0][0]*b[1][1]*b[2][2] + b[1][0]*b[2][1]*b[0][2] + b[2][0]*b[0][1]*b[1][2] - \
              b[0][2]*b[1][1]*b[2][0] - b[2][1]*b[1][2]*b[0][0] - b[2][2]*b[0][1]*b[1][0]
    else:
        vol = scale

    atomic_labels = inp_fh.readline() # yes, it is hardcoded for VASP5
    atom_numbers = [int(s) for s in inp_fh.readline().split()]
    nat = sum(atom_numbers)

    poscar_header += "%15.12f\n" % scale
    poscar_header += "%15.12f %15.12f %15.12f\n" % (b[0][0], b[0][1], b[0][2])
    poscar_header += "%15.12f %15.12f %15.12f\n" % (b[1][0], b[1][1], b[1][2])
    poscar_header += "%15.12f %15.12f %15.12f\n" % (b[2][0], b[2][1], b[2][2])
    poscar_header += atomic_labels
    poscar_header += ' '.join(str(x) for x in atom_numbers)+"\n"

    eigvals = [ 0.0 for i in range(nat*3) ]
    eigvecs = [ 0.0 for i in range(nat*3) ]
    norms =   [ 0.0 for i in range(nat*3) ]
    pos =     [ 0.0 for i in range(nat) ]
    for i in range(nat*3): # all frequencies should be supplied, regardless of requested to calculate
        inp_fh.readline() # empty line
        p = re.search(r'^\s*(\d+).+?([\.\d]+) cm-1', inp_fh.readline())
        eigvals[i] = float(p.group(2))
        #
        inp_fh.readline() # X         Y         Z           dx          dy          dz
        eigvec = []
        #
        for j in range(nat):
            tmp = inp_fh.readline().split()
            if i == 0: pos[j] = [ float(tmp[x]) for x in range(3) ] # get atomic positions only once
            #
            eigvec.append([ float(tmp[x]) for x in range(3,6) ])
            #
        eigvecs[i] = eigvec
        norms[i] = sqrt( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )
        #
    return first, last, nat, vol, pos, eigvals, eigvecs, norms, poscar_header

def parse_outcar(outcar_fn):
    import re
    import sys
    eps = []

    outcar_fh = open(outcar_fn, 'r')
    while True:
        line = outcar_fh.readline()
        if not line:
            break

        p = re.search(r'\s*MACROSCOPIC STATIC DIELECTRIC TENSOR', line)
        if p:
            outcar_fh.readline()
            eps.append( [float(x) for x in outcar_fh.readline().split()] )
            eps.append( [float(x) for x in outcar_fh.readline().split()] )
            eps.append( [float(x) for x in outcar_fh.readline().split()] )
            print "Found eps in "+outcar_fn+":"
            print eps
            return eps

    print "Couldn't find dielectric tensor in "+outcar_fn+", exiting..."
    sys.exit(1)

if __name__ == '__main__':
    import sys
    from math import pi, sqrt
    from shutil import move
    import re
    #
    if NDERIV == 2:
        disps = [-1, 1]
        coeffs = [-0.5, 0.5]
    else:
        print "Unknown value for NDERIV (use 2), exiting..."
        sys.exit(0)
    #
    # ELECTRONMASS_SI  = 9.10938215E-31      # Kg
    # AMU_SI           = 1.660538782E-27# Kg
    # AMU_AU           = AMU_SI / ELECTRONMASS_SI
    # sqrt_amu_au = sqrt(AMU_AU) # VASP divides by sqrt(AMU), we convert to sqrt(AU); AU aka m_e
    #
    inp_fh = open(VASP_RAMAN_INP_FN, 'r')
    first, last, nat, vol, pos, eigvals, eigvecs, norms, poscar_header = parse_inpcar(inp_fh)
    print "# atoms: %i, volume: %7.5f\n" % (nat, vol)
    #
    for i in range(first-1, last):
        eigval = eigvals[i]
        eigvec = eigvecs[i]
        norm = norms[i]
    
        print "Mode #: %i, frequency: %10.7f cm-1, norm: %10.7f" % ( i+1, eigval, norm )
    
        ra = [[0.0 for x in range(3)] for y in range(3)]
        for j in range(len(disps)):
            disp_filename = 'OUTCAR.%04d.%+d.out' % (i+1, disps[j])
    
            try:
                with open(disp_filename, 'r'):
                    file_exists =True
                    pass
            except IOError:
                file_exists=False
    
            if file_exists:
                print "File "+disp_filename+" exists, parsing..."
            else:
                print "File "+disp_filename+" not found, running VASP..."
                poscar_fh = open('POSCAR', 'w')
                poscar_fh.write("%s %4.1e \n" % (disp_filename, STEP_SIZE))
                poscar_fh.write(poscar_header)
                poscar_fh.write("Cartesian\n")
    
                for k in range(nat):
                    pos_disp = [ pos[k][l] + eigvec[k][l]*STEP_SIZE*disps[j]/norm for l in range(3)]
                    poscar_fh.write( '%15.10f %15.10f %15.10f\n' % (pos_disp[0], pos_disp[1], pos_disp[2]) )
                    #print '%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f' % (pos[k][0], pos[k][1], pos[k][2], dis[k][0], dis[k][1], dis[k][2])
                poscar_fh.close()
    
                #run VASP here
                os.system(VASP_RUN)
                move('OUTCAR', disp_filename)
    
            eps = parse_outcar(disp_filename)
            for m in range(3):
                for n in range(3):
                    ra[m][n]   += eps[m][n] * coeffs[j]/STEP_SIZE * norm * vol/(4.0*pi)
            #units: A^2/amu^1/2 =         dimless   * 1/A         * 1/amu^1/2  * A^3
    
        alpha = (ra[0][0] + ra[1][1] + ra[2][2])/3.0
        beta2 = ( (ra[0][0] - ra[1][1])**2 + (ra[0][0] - ra[2][2])**2 + (ra[1][1] - ra[2][2])**2 + 6.0 * (ra[0][1]**2 + ra[0][2]**2 + ra[1][2]**2) )/2.0
        print 'for mode %i (%10.5f): alpha: %10.7f; beta2: %10.7f; activity: %10.7f ' % (i+1, eigval, alpha, beta2, 45.0*alpha**2 + 7.0*beta2)
