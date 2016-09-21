#!/usr/bin/env python
#
# vasp_raman.py v. 0.6.0
#
# Raman off-resonant activity calculator
# using VASP as a back-end.
#
# Contributors: Alexandr Fonari (Georgia Tech)
# Shannon Stauffer (UT Austin)
#
# URL: http://raman-sc.github.io
#
# MIT license, 2013 - 2016
#


import re
import sys


def MAT_m_VEC(m, v):
    p = [ 0.0 for i in range(len(v)) ]
    for i in range(len(m)):
        assert len(v) == len(m[i]), 'Length of the matrix row is not equal to the length of the vector'
        p[i] = sum( [ m[i][j]*v[j] for j in range(len(v)) ] )
    return p


def T(m):
    p = [[ m[i][j] for i in range(len( m[j] )) ] for j in range(len( m )) ]
    return p


def parse_poscar(poscar_fh):
    # modified subroutine from phonopy 1.8.3 (New BSD license)
    #
    poscar_fh.seek(0) # just in case
    lines = poscar_fh.readlines()
    #
    scale = float(lines[1])
    if scale < 0.0:
        print "[parse_poscar]: ERROR negative scale not implemented."
        sys.exit(1)
    #
    b = []
    for i in range(2, 5):
        b.append([float(x)*scale for x in lines[i].split()[:3]])
    #
    vol = b[0][0]*b[1][1]*b[2][2] + b[1][0]*b[2][1]*b[0][2] + b[2][0]*b[0][1]*b[1][2] - \
          b[0][2]*b[1][1]*b[2][0] - b[2][1]*b[1][2]*b[0][0] - b[2][2]*b[0][1]*b[1][0]
    #
    try:
        num_atoms = [int(x) for x in lines[5].split()]
        line_at = 6
    except ValueError:
        symbols = [x for x in lines[5].split()]
        num_atoms = [int(x) for x in lines[6].split()]
        line_at = 7
    nat = sum(num_atoms)
    #
    if lines[line_at][0].lower() == 's':
        line_at += 1
    #
    if (lines[line_at][0].lower() == 'c' or lines[line_at][0].lower() == 'k'):
        is_scaled = False
    else:
        is_scaled = True
    #
    line_at += 1
    #
    positions = []
    for i in range(line_at, line_at + nat):
        pos = [float(x) for x in lines[i].split()[:3]]
        #
        if is_scaled:
            pos = MAT_m_VEC(T(b), pos)
        #
        positions.append(pos)
    #
    poscar_header = ''.join(lines[1:line_at-1]) # will add title and 'Cartesian' later
    return nat, vol, b, positions, poscar_header


def parse_env_params(params):
    tmp = params.strip().split('_')
    if len(tmp) != 4:
        print "[parse_env_params]: ERROR there should be exactly four parameters"
        sys.exit(1)
    #
    [first, last, nderiv, step_size] = [int(tmp[0]), int(tmp[1]), int(tmp[2]), float(tmp[3])]
    #
    return first, last, nderiv, step_size


#### subs for the output from VTST tools
def parse_freqdat(freqdat_fh, nat):
    freqdat_fh.seek(0) # just in case
    #
    eigvals = [ 0.0 for i in range(nat*3) ]
    #
    for i in range(nat*3): # all frequencies should be supplied, regardless of requested to calculate
        tmp = freqdat_fh.readline().split()
        eigvals[i] = float(tmp[0])
    #
    return eigvals
#
def parse_modesdat(modesdat_fh, nat):
    from math import sqrt
    modesdat_fh.seek(0) # just in case
    #
    eigvecs = [ 0.0 for i in range(nat*3) ]
    norms =   [ 0.0 for i in range(nat*3) ]
    #
    for i in range(nat*3): # all frequencies should be supplied, regardless of requested to calculate
        eigvec = []
        for j in range(nat):
            tmp = modesdat_fh.readline().split()
            eigvec.append([ float(tmp[x]) for x in range(3) ])
        #
        modesdat_fh.readline().split() # empty line
        eigvecs[i] = eigvec
        norms[i] = sqrt( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )
    #
    return eigvecs, norms
#### end subs for VTST
#
def get_modes_from_OUTCAR(outcar_fh, nat):
    from math import sqrt
    eigvals = [ 0.0 for i in range(nat*3) ]
    eigvecs = [ 0.0 for i in range(nat*3) ]
    norms   = [ 0.0 for i in range(nat*3) ]
    #
    outcar_fh.seek(0) # just in case
    while True:
        line = outcar_fh.readline()
        if not line:
            break
        #
        if "Eigenvectors after division by SQRT(mass)" in line:
            outcar_fh.readline() # empty line
            outcar_fh.readline() # Eigenvectors and eigenvalues of the dynamical matrix
            outcar_fh.readline() # ----------------------------------------------------
            outcar_fh.readline() # empty line
            #
            for i in range(nat*3): # all frequencies should be supplied, regardless of those requested to calculate
                outcar_fh.readline() # empty line
                p = re.search(r'^\s*(\d+).+?([\.\d]+) cm-1', outcar_fh.readline())
                eigvals[i] = float(p.group(2))
                #
                outcar_fh.readline() # X         Y         Z           dx          dy          dz
                eigvec = []
                #
                for j in range(nat):
                    tmp = outcar_fh.readline().split()
                    #
                    eigvec.append([ float(tmp[x]) for x in range(3,6) ])
                    #
                eigvecs[i] = eigvec
                norms[i] = sqrt( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )
            #
            return eigvals, eigvecs, norms
        #
    print "[get_modes_from_OUTCAR]: ERROR Couldn't find 'Eigenvectors after division by SQRT(mass)' in OUTCAR. Use 'NWRITE=3' in INCAR. Exiting..."
    sys.exit(1)
#
def get_epsilon_from_OUTCAR(outcar_fh):
    epsilon = []
    #
    outcar_fh.seek(0) # just in case
    while True:
        line = outcar_fh.readline()
        if not line:
            break
        #
        if "MACROSCOPIC STATIC DIELECTRIC TENSOR" in line:
            outcar_fh.readline()
            epsilon.append([float(x) for x in outcar_fh.readline().split()])
            epsilon.append([float(x) for x in outcar_fh.readline().split()])
            epsilon.append([float(x) for x in outcar_fh.readline().split()])
            return epsilon
    #
    raise RuntimeError("[get_epsilon_from_OUTCAR]: ERROR Couldn't find dielectric tensor in OUTCAR")
    return 1
#
if __name__ == '__main__':
    from math import pi
    from shutil import move
    import os
    import datetime
    import time
    #import argparse
    import optparse
    #
    print ""
    print "    Raman off-resonant activity calculator,"
    print "    using VASP as a back-end."
    print ""
    print "    Contributors: Alexandr Fonari  (Georgia Tech)"
    print "                  Shannon Stauffer (UT Austin)"
    print "    MIT License, 2013"
    print "    URL: http://raman-sc.github.io"
    print "    Started at: "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    print ""
    #
    description  = "Before run, set environment variables:\n"
    description += "    VASP_RAMAN_RUN='mpirun vasp'\n"
    description += "    VASP_RAMAN_PARAMS='[first-mode]_[last-mode]_[nderiv]_[step-size]'\n\n"
    description += "bash one-liner is:\n"
    description += "VASP_RAMAN_RUN='mpirun vasp' VASP_RAMAN_PARAMS='1_2_2_0.01' python vasp_raman.py"
    #
    parser = optparse.OptionParser(description=description)
    parser.add_option('-g', '--gen', help='Generate POSCAR only', action='store_true')
    parser.add_option('-u', '--use_poscar', help='Use provided POSCAR in the folder, USE WITH CAUTION!!', action='store_true')
    (options, args) = parser.parse_args()
    #args = vars(parser.parse_args())
    args = vars(options)
    #
    VASP_RAMAN_RUN = os.environ.get('VASP_RAMAN_RUN')
    if VASP_RAMAN_RUN == None:
        print "[__main__]: ERROR Set environment variable 'VASP_RAMAN_RUN'"
        print ""
        parser.print_help()
        sys.exit(1)
    print "[__main__]: VASP_RAMAN_RUN='"+VASP_RAMAN_RUN+"'"
    #
    VASP_RAMAN_PARAMS = os.environ.get('VASP_RAMAN_PARAMS')
    if VASP_RAMAN_PARAMS == None:
        print "[__main__]: ERROR Set environment variable 'VASP_RAMAN_PARAMS'"
        print ""
        parser.print_help()
        sys.exit(1)
    print "[__main__]: VASP_RAMAN_PARAMS='"+VASP_RAMAN_PARAMS+"'"
    #
    first, last, nderiv, step_size = parse_env_params(VASP_RAMAN_PARAMS)
    assert first >= 1,    '[__main__]: First mode should be equal or larger than 1'
    assert last >= first, '[__main__]: Last mode should be equal or larger than first mode'
    if args['gen']: assert last == first, "[__main__]: '-gen' mode -> only generation for the one mode makes sense"
    assert nderiv == 2,   '[__main__]: At this time, nderiv = 2 is the only supported'
    disps = [-1, 1]      # hardcoded for
    coeffs = [-0.5, 0.5] # three point stencil (nderiv=2)
    #
    try:
        poscar_fh = open('POSCAR.phon', 'r')
    except IOError:
        print "[__main__]: ERROR Couldn't open input file POSCAR.phon, exiting...\n"
        sys.exit(1)
    #
    # nat, vol, b, poscar_header = parse_poscar_header(poscar_fh)
    nat, vol, b, pos, poscar_header = parse_poscar(poscar_fh)
    print pos
    #print poscar_header
    #sys.exit(0)
    #
    # either use modes from vtst tools or VASP
    if os.path.isfile('freq.dat') and os.path.isfile('modes_sqrt_amu.dat'):
        try:
            freqdat_fh = open('freq.dat', 'r')
        except IOError:
            print "[__main__]: ERROR Couldn't open freq.dat, exiting...\n"
            sys.exit(1)
        #
        eigvals = parse_freqdat(freqdat_fh, nat)
        freqdat_fh.close()
        #
        try: 
            modes_fh = open('modes_sqrt_amu.dat' , 'r')
        except IOError:
            print "[__main__]: ERROR Couldn't open modes_sqrt_amu.dat, exiting...\n"
            sys.exit(1)
        #
        eigvecs, norms = parse_modesdat(modes_fh, nat)
        modes_fh.close()
    #
    elif os.path.isfile('OUTCAR.phon'):
        try:
            outcar_fh = open('OUTCAR.phon', 'r')
        except IOError:
            print "[__main__]: ERROR Couldn't open OUTCAR.phon, exiting...\n"
            sys.exit(1)
        #
        eigvals, eigvecs, norms = get_modes_from_OUTCAR(outcar_fh, nat)
        outcar_fh.close()
    #
    else:
        print "[__main__]: Neither OUTCAR.phon nor freq.dat/modes_sqrt_amu.dat were found, nothing to do, exiting..."
        sys.exit(1)
    #
    output_fh = open('vasp_raman.dat', 'w')
    output_fh.write("# mode    freq(cm-1)    alpha    beta2    activity\n")
    for i in range(first-1, last):
        eigval = eigvals[i]
        eigvec = eigvecs[i]
        norm = norms[i]
        #
        print ""
        print "[__main__]: Mode #%i: frequency %10.7f cm-1; norm: %10.7f" % ( i+1, eigval, norm )
        #
        ra = [[0.0 for x in range(3)] for y in range(3)]
        for j in range(len(disps)):
            disp_filename = 'OUTCAR.%04d.%+d.out' % (i+1, disps[j])
            #
            try:
                outcar_fh = open(disp_filename, 'r')
                print "[__main__]: File "+disp_filename+" exists, parsing..."
            except IOError:
                if args['use_poscar'] != True:
                    print "[__main__]: File "+disp_filename+" not found, preparing displaced POSCAR"
                    poscar_fh = open('POSCAR', 'w')
                    poscar_fh.write("%s %4.1e \n" % (disp_filename, step_size))
                    poscar_fh.write(poscar_header)
                    poscar_fh.write("Cartesian\n")
                    #
                    for k in range(nat):
                        pos_disp = [ pos[k][l] + eigvec[k][l]*step_size*disps[j]/norm for l in range(3)]
                        poscar_fh.write( '%15.10f %15.10f %15.10f\n' % (pos_disp[0], pos_disp[1], pos_disp[2]) )
                        #print '%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f' % (pos[k][0], pos[k][1], pos[k][2], dis[k][0], dis[k][1], dis[k][2])
                    poscar_fh.close()
                else:
                    print "[__main__]: Using provided POSCAR"
                #
                if args['gen']: # only generate POSCARs
                    poscar_fn = 'POSCAR.%+d.out' % disps[j]
                    move('POSCAR', poscar_fn)
                    print "[__main__]: '-gen' mode -> "+poscar_fn+" with displaced atoms have been generated"
                    #
                    if j+1 == len(disps): # last iteration for the current displacements list
                        print "[__main__]: '-gen' mode -> POSCAR files with displaced atoms have been generated, exiting now"
                        sys.exit(0)
                else: # run VASP here
                    print "[__main__]: Running VASP..."
                    os.system(VASP_RAMAN_RUN)
                    try:
                        move('OUTCAR', disp_filename)
                    except IOError:
                        print "[__main__]: ERROR Couldn't find OUTCAR file, exiting..."
                        sys.exit(1)
                    #
                    outcar_fh = open(disp_filename, 'r')
            #
            try:
                eps = get_epsilon_from_OUTCAR(outcar_fh)
                outcar_fh.close()
            except Exception, err:
                print err
                print "[__main__]: Moving "+disp_filename+" back to 'OUTCAR' and exiting..."
                move(disp_filename, 'OUTCAR')
                sys.exit(1)
            #
            for m in range(3):
                for n in range(3):
                    ra[m][n]   += eps[m][n] * coeffs[j]/step_size * norm * vol/(4.0*pi)
            #units: A^2/amu^1/2 =         dimless   * 1/A         * 1/amu^1/2  * A^3
        #
        alpha = (ra[0][0] + ra[1][1] + ra[2][2])/3.0
        beta2 = ( (ra[0][0] - ra[1][1])**2 + (ra[0][0] - ra[2][2])**2 + (ra[1][1] - ra[2][2])**2 + 6.0 * (ra[0][1]**2 + ra[0][2]**2 + ra[1][2]**2) )/2.0
        print ""
        print "! %4i  freq: %10.5f  alpha: %10.7f  beta2: %10.7f  activity: %10.7f " % (i+1, eigval, alpha, beta2, 45.0*alpha**2 + 7.0*beta2)
        output_fh.write("%03i  %10.5f  %10.7f  %10.7f  %10.7f\n" % (i+1, eigval, alpha, beta2, 45.0*alpha**2 + 7.0*beta2))
        output_fh.flush()
    #
    output_fh.close()
