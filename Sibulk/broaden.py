#!/usr/bin/env python

# Copyright (c) "2012, by Alexandr Fonari
#                License: MIT License"

import numpy as np
from pylab import plot
from math import sqrt, log, fabs, fsum, exp
import sys
import argparse


def to_plot(hw,ab,gam=0.001,type='Gaussian'):
    fmin = min(hw)
    fmax = max(hw)
    erange = np.arange(fmin-40*gam,fmax+40*gam,gam/10)#np.arange(fmin-40*gam,fmax+40*gam,gam/10)
    spectrum = 0.0*erange
    for i in range(len(hw)):
        if type=='Gaussian':
            spectrum += (2*np.pi)**(-.5)/gam*np.exp(np.clip(-1.0*(hw[i]-erange)**2/(2*gam**2),-300,300))
        elif type=='Lorentzian':
            spectrum += ab[i]*1/np.pi*gam/((hw[i]-erange)**2+gam**2)

    return erange, spectrum

#    Emin = min(E for E,ab in data)
#    Emax = max(E for E,ab in data)
#    abmin = min([0] + [ab for E,ab in data if ab > 0])
#    abmin = max(abmin,abmin_floor)
#    Edel = Emax-Emin
#    Es = np.linspace(Emin-Edel/10.,Emax+Edel/10.,2000)
#    Spectrum = 0.0*Es
#    for E,ab in data:
#        if ab > 0:
#            Spectrum += ab*np.exp(-alpha*(Es-E)**2)
#        elif plot_zeros:
#            Spectrum -= abmin*exp(-alpha*(Es-E)**2)
#    return Es,Spectrum


hw=np.genfromtxt(sys.argv[1], dtype=float)
cm1 = hw[:,0]
int1 = hw[:,1]
int1 /= np.max(np.abs(int1),axis=0)
#print cm1
Es1,Spectrum1 = to_plot(cm1, int1, 15.0, 'Lorentzian')

#cm1 = hw[:,2]
#int1 = hw[:,3]
#int1 /= np.max(np.abs(int1),axis=0)
#print cm1
#Es2,Spectrum2 = to_plot(cm1, int1, 15.0, 'Lorentzian')

filename = sys.argv[1]+"-a-pw91.out"
print 'Phonons.WriteVibDOSFile: Writing', filename
f = open(filename,'w')
f.write('# energy/meV  VibDOS \n')
for i in range(len(Es1)):
    f.write('%.5e   %.5e\n' % (Es1[i],Spectrum1[i]))
f.close()

#filename = sys.argv[1]+"-paw.out"
#print 'Phonons.WriteVibDOSFile: Writing', filename
#f = open(filename,'w')
#f.write('# energy/meV  VibDOS \n')
#for i in range(len(Es2)):
#    f.write('%.5e   %.5e\n' % (Es2[i],Spectrum2[i]))
#f.close()
