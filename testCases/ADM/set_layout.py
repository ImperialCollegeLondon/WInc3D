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
Lx=2*pi*H
Ly=H
Lz=pi*H
Nt=4*6
sx=pi*H/4.
sz=pi*H/6.

TurbineName=[]
TurbinesCoords=np.zeros((Nt,2))

for j in range(0,6):
    for i in range(0,4):
        TurbineName.append('WT_'+str(i)+str(j))
        TurbinesCoords[i+4*j,0]=sx/2.+i*sx
        TurbinesCoords[i+4*j,1]=sz/2.+j*sz

with open('Turbines.ad','w') as fturout:
    for i in range(len(TurbineName)):
        fturout.write(str(TurbinesCoords[i,0])+' '+str(100.)+' '+str(TurbinesCoords[i,1])+'  1. 0. 0. 100.\n')

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

plt.figure(1)
plt.plot(TurbinesCoords[:,0]/1000,TurbinesCoords[:,1]/1000,'ko')
plt.xlabel("x [km]",fontsize=16)
plt.ylabel("y [km]",fontsize=16)
plt.arrow(1,5, 6, 0, head_width=0.25, head_length=0.4, fc='b', ec='b')
plt.text(7.5,5,r"$\theta=270^\circ$",color='blue',fontsize=14)
plt.arrow(1,5,5.196,3 , head_width=0.25, head_length=0.4, fc='b', ec='b')
plt.text(6.5,8.5,r"$\theta=240^\circ$",color='blue',fontsize=14)
plt.arrow(1,5,5.196,-3 , head_width=0.25, head_length=0.4, fc='b', ec='b')
plt.text(6.5,1.0,r"$\theta=300^\circ$",color='blue',fontsize=14)
plt.xlim((0,6.16))
plt.ylim((-1.12,0))
plt.savefig('fig_HornsRev_accurate.eps')

