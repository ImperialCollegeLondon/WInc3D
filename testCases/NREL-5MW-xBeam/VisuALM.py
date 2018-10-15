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


# Compute Averages
blade1=al('Blade1')
blade1.read('ALM/'+str(itime)+'/NREL-5MW_blade_1.load')
blade2=al('Blade2')
blade2.read('ALM/'+str(itime)+'/NREL-5MW_blade_2.load')
blade3=al('Blade3')
blade3.read('ALM/'+str(itime)+'/NREL-5MW_blade_3.load')
tower=al('tower')
tower.read('ALM/'+str(itime)+'/NREL-5MW_tower.load')


