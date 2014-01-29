#!/usr/bin/env python
#############################
#The MIT License (MIT)
#
#Copyright (c) 2013 Alexandr Fonari
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in
#the Software without restriction, including without limitation the rights to
#use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
#the Software, and to permit persons to whom the Software is furnished to do so,
#subject to the following conditions:
#
#    The above copyright notice and this permission notice shall be included in all
#    copies or substantial portions of the Software.
#
#    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#    FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
#    COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
#    IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#    CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#############################
##Modified by Shannon Stauffer.  University of Texas at Austin
##January 2014
#############################
####This script reads freq.dat and modes_sqrt_amu.dat files from a DYMMAT run using Henkelman group scripts
####This script will convert POSCAR to Cartesian coordinates if needed
#############################

STEP_SIZE = 0.01 # Angstrom
NDERIV = 2

VASP_RAMAN_INP_FN='freqs'

import os
import numpy as np

#################
##Write POSCAR
#################
#def write_poscar(poscar_header, pos):
#    poscar_fh = open('POSCAR', 'w')
#    poscar_fh.write("VASP 5 format\n" )
#    poscar_fh.write(poscar_header)
#    poscar_fh.write("Cartesian")
#    np.savetxt(poscar_fh, pos)
#    poscar_fh.close()

#############################
########Parse POSCAR
#############################
def parse_inphead(inp_head):
    import sys

    poscar_header = ""
    vol = 0.0
    b = []
    atom_numbers = []

    inp_head.readline() #top line
    scale = float(inp_head.readline()) # scale
    for i in range(3): b.append( [float(s) for s in inp_head.readline().split()] )

    if scale > 0.0:
        b = [[ b[i][j]*scale for i in range(3)] for j in range(3) ]
        scale = 1.0

        vol = b[0][0]*b[1][1]*b[2][2] + b[1][0]*b[2][1]*b[0][2] + b[2][0]*b[0][1]*b[1][2] - \
              b[0][2]*b[1][1]*b[2][0] - b[2][1]*b[1][2]*b[0][0] - b[2][2]*b[0][1]*b[1][0]
    else:
        vol = scale

    atomic_labels = inp_head.readline() # yes, it is hardcoded for VASP5
    atom_numbers = [int(s) for s in inp_head.readline().split()]
    nat = sum(atom_numbers)
    conf_type = inp_head.readline().split() # CART or DIRECT

    poscar_header += "%15.12f\n" % scale
    poscar_header += "%15.12f %15.12f %15.12f\n" % (b[0][0], b[0][1], b[0][2])
    poscar_header += "%15.12f %15.12f %15.12f\n" % (b[1][0], b[1][1], b[1][2])
    poscar_header += "%15.12f %15.12f %15.12f\n" % (b[2][0], b[2][1], b[2][2])
    poscar_header += atomic_labels
    poscar_header += ' '.join(str(x) for x in atom_numbers)+"\n"

#    inp_head.readline() #dir line !!! Change for diff VASP versions
    pos =     [ 0.0 for i in range(nat) ]
    for j in range(nat):
    	tmp = inp_head.readline().split()
    #    print tmp
      	pos[j] = [ float(tmp[x]) for x in range(3) ]

    if conf_type[0][0] == "D":
        for i in range(nat):
            pos[i] = np.dot(pos[i], b)
    #    write_poscar(poscar_header, pos)

    return nat, vol, pos, poscar_header

#############################
########Parse freq.dat
#############################
def parse_inpeval(inp_eval, nat):

    eigvals = [ 0.0 for i in range(nat*3) ]
    real = [ 0.0 for i in range(nat*3) ]
    last = len(eigvals)
    first = 0

    for i in range(nat*3): # all frequencies should be supplied, regardless of requested to calculate
        tmp = inp_eval.readline().split()
        eigvals[i] = float(tmp[0])
        real[i] = int(tmp[3])
        if int(tmp[3]) == 0  and first == 0:
        	first = i+1
    return first, last, eigvals, real

#############################
########Parse modes.dat
#############################
def parse_inpevec(inp_evec, nat):

    eigvecs = [ 0.0 for i in range(nat*3) ]
    norms =   [ 0.0 for i in range(nat*3) ]

    for i in range(nat*3): # all frequencies should be supplied, regardless of requested to calculate
        tmp = inp_evec.readline().split()#vecs for 1 mode
        eigvec = []

        for j in range(nat):
            eigvec.append([ float(tmp[x]) for x in range(j,j+3) ])

        eigvecs[i] = eigvec
        norms[i] = sqrt( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )

    return eigvecs, norms 

if __name__ == '__main__':
    import sys
    from math import pi, sqrt
    from shutil import move
    from shutil import copy
    import re
    #
    if NDERIV == 2:
        disps = [-1, 1]
        coeffs = [-0.5, 0.5]
    else:
        print "Unknown value for NDERIV (use 2), exiting..."
        sys.exit(0)

#################
##Read POSCAR, freq.dat, and modes.dat
#################
    inp_head = open('POSCAR', 'r')
    nat, vol, pos, poscar_header = parse_inphead(inp_head)
    inp_eval = open('freq.dat' , 'r')
    first, last, eigvals, real = parse_inpeval(inp_eval, nat)
    inp_evec = open('modes_sqrt_amu.dat' , 'r')
    eigvecs, norms = parse_inpevec(inp_evec, nat)
    print "# atoms: %i, volume: %7.5f\n" % (nat, vol)

#################
##Write "VASP_RAMAN_INP_FN"
################
    inp_fh = open(VASP_RAMAN_INP_FN , 'w')
    inp_fh.write("%i %i \n" % (first, last))
    inp_fh.write(poscar_header)
    inp_fh.write("\n")
    for i in range(nat*3):
    	fir = real[i]
    	eigval = eigvals[i]
    	eigvec = eigvecs[i]
    	norm = norms[i]
        inp_fh.write("%i %i =   %10.7f cm-1 \n" % ( i+1, fir, eigval ))
        #inp_fh.write("%i %i =   %10.7f cm-1 \n" % ( i+1, fir, eigval, norm ))
        inp_fh.write("           X                Y                 Z              dx              dy              dz \n" )
        for k in range(nat):
            disp = eigvec[k]
            inp_fh.write( '%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n' % (pos[k][0], pos[k][1], pos[k][2], disp[0], disp[1], disp[2]))
    
        inp_fh.write("\n")
    inp_fh.close()



