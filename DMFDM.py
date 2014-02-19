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
##February 2014
#############################
####This script reads freq.dat and modes_sqrt_amu.dat files from a DYMMAT run using Henkelman group scripts
####This script will convert POSCAR to Cartesian coordinates if needed
#############################

#################
##Write POSCAR
#################
#def write_poscar(poscar_header, pos):
#    poscar_fh = open('POSCAR', 'w')
#    poscar_fh.write("VASP 5 format\n" )
#    poscar_fh.write(poscar_header)
#    np.savetxt(poscar_fh, pos)
#    poscar_fh.close()

#############################
########Parse POSCAR
#############################
def parse_pos(inp_head, nat, b):
    import sys

    conf_type=inp_head.readline().split() #dir line !!! Change for diff VASP versions
    pos =     [ 0.0 for i in range(nat) ]
    for j in range(nat):
    	tmp = inp_head.readline().split()
    #    print tmp
      	pos[j] = [ float(tmp[x]) for x in range(3) ]

    if conf_type[0][0] == "D":
        for i in range(nat):
            pos[i] = np.dot(pos[i], b)
        write_poscar(poscar_header, pos)

    return pos

#############################
########Parse freq.dat
#############################
def parse_eval(inp_eval, nat):

    eigvals = [ 0.0 for i in range(nat*3) ]

    for i in range(nat*3): # all frequencies should be supplied, regardless of requested to calculate
        tmp = inp_eval.readline().split()
        eigvals[i] = float(tmp[0])

    return eigvals

#############################
########Parse modes.dat
#############################
def parse_evec(inp_evec, nat):
    from math import sqrt

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

