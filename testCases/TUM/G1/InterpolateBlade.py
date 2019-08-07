#!/usr/bin/env python3

import math
import argparse
import csv
import f90nml
import matplotlib
import numpy as np
from scipy import interpolate
from pylab import *
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description="Script to extract Boundary Layer Flow statistics from *.dat files")
parser.add_argument("-v","--verbose",action="store_true",help="Print location")
parser.add_argument("-p","--plot",action="store_true",help="Plots the wake profiles")
parser.add_argument("-w","--write",action="store_true",help="Write results in a .csv file")
parser.add_argument("NumElem", type=int, help="number of blade elements")
parser.add_argument("NumElemTower", type=int, help="number of blade elements")

args = parser.parse_args()
NElem = args.NumElem
NElemT = args.NumElemTower

def interpolate_to_nearest(values_to_correct, valid_values):
    values_to_correct = np.asarray(values_to_correct)
    n = len(values_to_correct)
    corrected_values = np.zeros((n))
    for i in range(n):
        pos = (np.abs(valid_values - values_to_correct[i])).argmin()
        corrected_values[i] = valid_values[pos]
    return corrected_values

# PARAMETERS
case = 'G1'
R = 0.55 # Tip radius in meters
tower_height = 0.756 # in meters

fname_blade = case + '_Blade.txt'
fname_structure = case + '_Structure.txt'
fname_tower = case + '_Tower.txt'
delimiter='\t' # For the input files

# NumElem=108 NumElemTower=14
# rR_New=np.linspace(rR_ref[0],rR_ref[-1],NElem)
rR_New = np.array([0.0981818182,
                    0.1066100182,
                    0.1150382364,
                    0.1234663636,
                    0.1318947273,
                    0.1403229091,
                    0.1487510909,
                    0.1571792727,
                    0.1656074545,
                    0.1740356364,
                    0.1824638182,
                    0.1908921818,
                    0.1993203636,
                    0.2077485455,
                    0.2161767273,
                    0.2246049091,
                    0.2330330909,
                    0.2414612727,
                    0.2498896364,
                    0.2583178182,
                    0.266746,
                    0.2751741818,
                    0.2836018182,
                    0.2920309091,
                    0.3004581818,
                    0.3088872727,
                    0.3173145455,
                    0.3257436364,
                    0.3341709091,
                    0.3426,
                    0.3510272727,
                    0.3594563636,
                    0.3678836364,
                    0.3763127273,
                    0.38474,
                    0.3931690909,
                    0.4015981818,
                    0.4100254545,
                    0.4184545455,
                    0.4268818182,
                    0.4353109091,
                    0.4437381818,
                    0.4521672727,
                    0.4605945455,
                    0.4690236364,
                    0.4774509091,
                    0.48588,
                    0.4943072727,
                    0.5027363636,
                    0.5111636364,
                    0.5195927273,
                    0.52802,
                    0.5364490909,
                    0.5448763636,
                    0.5533054545,
                    0.5617327273,
                    0.5701618182,
                    0.5785890909,
                    0.5870181818,
                    0.5954454545,
                    0.6038745455,
                    0.6123018182,
                    0.6207309091,
                    0.6291581818,
                    0.6375872727,
                    0.6460145455,
                    0.6544436364,
                    0.6628709091,
                    0.6713,
                    0.6797272727,
                    0.6881563636,
                    0.6965836364,
                    0.7050127273,
                    0.7134418182,
                    0.7218690909,
                    0.7302981818,
                    0.7387254545,
                    0.7471545455,
                    0.7555818182,
                    0.7640109091,
                    0.7724381818,
                    0.7808672727,
                    0.7892945455,
                    0.7977236364,
                    0.8061509091,
                    0.81458,
                    0.8230072727,
                    0.8314363636,
                    0.8398636364,
                    0.8482927273,
                    0.85672,
                    0.8651490909,
                    0.8735763636,
                    0.8820054545,
                    0.8904327273,
                    0.8988618182,
                    0.9072890909,
                    0.9157181818,
                    0.9241454545,
                    0.9325745455,
                    0.9410018182,
                    0.9494309091,
                    0.9578581818,
                    0.9662872727,
                    0.9747145455,
                    0.9831436364,
                    0.9915709091,
                    1])

last = rR_New[-1]
rR_New = rR_New[::4]
if not rR_New[-1] == last:
    rR_New = np.concatenate((rR_New, np.array([last])), axis=0)

valid_thickness = np.array([1,
                            0.9031775485,
                            0.7931077835,
                            0.6789127693,
                            0.5393601293,
                            0.3998074893,
                            0.2884782239,
                            0.2194482997,
                            0.1504198647,
                            0.1168218079,
                            0.1014359629,
                            0.0860501179,
                            0.0848476727,
                            0.0848657801,
                            0.0848762848,
                            0.0848467198,
                            0.0848171548,
                            0.0848])

NElem = len(rR_New)
# Original Blades
A=np.genfromtxt(fname_blade,delimiter=delimiter,skip_header=0)
rR_ref=A[:,0]/R
cR_ref=A[:,1]/R
pitch_ref=A[:,2]
t2c_ref=A[:,3]

cR_New=np.interp(rR_New,rR_ref,cR_ref)
pitch_New=np.interp(rR_New,rR_ref,pitch_ref)
t2c_New=np.interp(rR_New,rR_ref,t2c_ref)
# Round t2c_New to the closest value in the airfoil list
t2c_New = interpolate_to_nearest(t2c_New, valid_thickness)

# Read the strucrure part and make sure that the structural and aerodynamic models are the same
Astruct=np.genfromtxt(fname_structure,delimiter=delimiter,skip_header=2)
rR_struct=Astruct[:,0]
AeroCent_struct=Astruct[:,1]  #Aerodynamic center
StrTwist_struct=Astruct[:,2]  #Structural Twist
BMassDen_struct=Astruct[:,3]  #B
FlpStff_struct=Astruct[:,4]    # Flapwise stiffness (EI)
EdgStff_struct=Astruct[:,5]   # Edge-wise stiffness (EI)
GJStff_struct=Astruct[:,6]    # Torsional stiffness (GJ)
EAStff_struct=Astruct[:,7]    #
Alpha_struct=Astruct[:,8]
FlpInert_struct=Astruct[:,9]
EdgInert_struct=Astruct[:,10]
Precrv_struct=Astruct[:,11]
Preswp_struct=Astruct[:,12]
FlpcgOf_struct=Astruct[:,13]
EdgcgOf_struct=Astruct[:,14]
FlpEAOf_struct=Astruct[:,15]
EdgEAOf_struct=Astruct[:,16]
rR_struct_New=np.linspace(rR_ref[0],rR_ref[-1],NElem)
AeroCent_struct_New=np.interp(rR_struct_New,rR_struct,AeroCent_struct)
StrTwist_struct_New=np.interp(rR_struct_New,rR_struct,StrTwist_struct)
BMassDen_struct_New=np.interp(rR_struct_New,rR_struct,BMassDen_struct)
FlpStff_struct_New=np.interp(rR_struct_New,rR_struct,FlpStff_struct)
EdgStff_struct_New=np.interp(rR_struct_New,rR_struct,EdgStff_struct)
GJStff_struct_New=np.interp(rR_struct_New,rR_struct,GJStff_struct)
EAStff_struct_New=np.interp(rR_struct_New,rR_struct,EAStff_struct)
Alpha_struct_New=np.interp(rR_struct_New,rR_struct,Alpha_struct)
FlpInert_struct_New=np.interp(rR_struct_New,rR_struct,FlpInert_struct)
EdgInert_struct_New=np.interp(rR_struct_New,rR_struct,EdgInert_struct)
Precrv_struct_New=np.interp(rR_struct_New,rR_struct,Precrv_struct)
Preswp_struct_New=np.interp(rR_struct_New,rR_struct,Preswp_struct)
FlpcgOf_struct_New=np.interp(rR_struct_New,rR_struct,FlpcgOf_struct)
EdgcgOf_struct_New=np.interp(rR_struct_New,rR_struct,EdgcgOf_struct)
FlpEAOf_struct_New=np.interp(rR_struct_New,rR_struct,FlpEAOf_struct)
EdgEAOf_struct_New=np.interp(rR_struct_New,rR_struct,EdgEAOf_struct)

with open(case + 'Blade_N'+str(NElem)+'.al','w') as fout:
    fout.write('R  : '+str(R)+' \n')
    fout.write('Spanwise  : 0.0 0.0 1.0 \n')
    fout.write('NStations : '+str(NElem)+'\n')
    for j in range(0,NElem):
        fout.write(str(rR_New[j])+'\t'+str(cR_New[j])+'\t'+str(pitch_New[j])+'\t'+str(t2c_New[j])+'\n')

with open(case + 'Blade_N'+str(NElem)+'.struct','w') as fout:
    fout.write('NStations : '+str(NElem)+'\n')
    for j in range(0,NElem):
        fout.write(str(rR_struct_New[j])+'\t'+str(AeroCent_struct_New[j])+'\t'+str(StrTwist_struct_New[j])+'\t'+str(BMassDen_struct_New[j])+'\t'+str(FlpStff_struct_New[j])+'\t'+str(EdgStff_struct_New[j])+'\t'+str(GJStff_struct_New[j])+'\t'+str(EAStff_struct_New[j])+'\t'+str(Alpha_struct_New[j])+'\t'+str(FlpInert_struct_New[j])+'\t'+str(EdgInert_struct_New[j])+'\t'+str(Precrv_struct_New[j])+'\t'+str(Preswp_struct_New[j])+'\t'+str(FlpcgOf_struct_New[j])+'\t'+str(EdgcgOf_struct_New[j])+'\t'+str(FlpEAOf_struct_New[j])+'\t'+str(EdgEAOf_struct_New[j])+'\n')

B=np.genfromtxt(fname_tower,delimiter=delimiter)
rR_Tref=B[:,0]/tower_height
cR_Tref=B[:,1]/tower_height
pitch_Tref=B[:,2]
t2c_Tref=B[:,3]
rR_TNew=np.linspace(rR_Tref[0],rR_Tref[-1],NElemT)
cR_TNew=np.interp(rR_TNew,rR_Tref,cR_Tref)
pitch_TNew=np.interp(rR_TNew,rR_Tref,pitch_Tref)
t2c_TNew=np.interp(rR_TNew,rR_Tref,t2c_Tref)

with open(case +'Tower_N'+str(NElemT)+'.al','w') as fout:
    fout.write('Length  : '+str(tower_height) +'\n')
    fout.write('Spanwise  : 0.0 -1.0 0.0 \n')
    fout.write('NStations : '+str(NElemT)+'\n')
    for j in range(0,NElemT):
        fout.write(str(rR_TNew[j])+'\t'+str(cR_TNew[j])+'\t'+str(pitch_TNew[j])+'\t'+str(t2c_TNew[j])+'\n')
