"""
Script sets the location and characteristics
of the Wind Turbines for the Lillgrund offshore wind farm
for incompact3d
Author : Georgios (Yorgos) Deskos 2017
"""
import os
from math import pi, sin, cos
import argparse
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import math


H=1000.
Lx=4*pi*H
Ly=H
Lz=pi*H
nx=9
nz=6
Nt=nx*nz
sx=pi*H/6.
sz=pi*H/6.

TurbineName=[]
TurbinesCoords=np.zeros((Nt,2))

for j in range(0,nz):
    for i in range(0,nx):
        TurbineName.append('WT_'+str(i)+str(j))
        TurbinesCoords[i+nx*j,0]=2*pi*H+sx/2.+i*sx
        TurbinesCoords[i+nx*j,1]=sz/2.+j*sz

with open('Turbines.ad','w') as fturout:
    for i in range(len(TurbineName)):
        fturout.write(str(TurbinesCoords[i,0])+' '+str(100.)+' '+str(TurbinesCoords[i,1])+'  1. 0. 0. 100.\n')

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

plt.figure(1)
plt.plot(TurbinesCoords[:,0]/(pi*1000),TurbinesCoords[:,1]/(1000*pi),'ko')
plt.xlim(0,4)
plt.savefig('fig_HornsRev_accurate.eps')
