#===================================================
#Python script used to visualise the elastic ALM
#Copyright G.Deskos 2018
#===================================================
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import time

class al:
    def __init__(self,name):
        self.name=name

    def read(self,filename):
        self.DATA=np.genfromtxt(filename,skip_header=1,delimiter=',')
        self.NElem=len(self.DATA)
        self.ielem=self.DATA[:,0]
        self.x=self.DATA[:,1]
        self.y=self.DATA[:,2]
        self.z=self.DATA[:,3]
        self.FN=self.DATA[:,17]
        self.FT=self.DATA[:,18]

itime=0
# Compute Averages
blade1=al('Blade1')
blade1.read('ALM/'+str(itime)+'/NREL-5MW_blade_1.load')
blade2=al('Blade2')
blade2.read('ALM/'+str(itime)+'/NREL-5MW_blade_2.load')
blade3=al('Blade3')
blade3.read('ALM/'+str(itime)+'/NREL-5MW_blade_3.load')
assert(blade1.Nelem==blade2.Nelem)
assert(blade2.Nelem==blade3.Nelem)


CP=[]
Power=[]
CT=[]
Thrust=[]
FN_ave=np.zeros((blade1.Nelem))
FT_ave=np.zeros((blade1.Nelem))
for itime in range(14):
    A=np.genfromtxt('ALM/'+str(itime)+'/NREL-5MW.perf',skip_header=1,delimiter=',')
    Power.append(A[-1])
    CP.append(A[6])
    Thrust.append(A[-3])
    CT.append(A[4])
    #Compute Time-average loads
    blade1=al('Blade1')
    blade1.read('ALM/'+str(itime)+'/NREL-5MW_blade_1.load')
    blade2=al('Blade2')
    blade2.read('ALM/'+str(itime)+'/NREL-5MW_blade_2.load')
    blade3=al('Blade3')
    blade3.read('ALM/'+str(itime)+'/NREL-5MW_blade_3.load')
    FN_ave=FN_ave+(blade1.FN+blade2.FN+blade3.FN)/3.0
    FT_ave=FT_ave+(blade1.FT+blaed2.FT+blade3.FT)/3.0


plt.figure(1)
plt.plot(CP)
plt.show()

