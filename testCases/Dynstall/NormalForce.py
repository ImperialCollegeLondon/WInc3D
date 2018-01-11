import numpy as np
import matplotlib.pyplot as plt
from math import pi

StaticData1=np.genfromtxt('S809.air',skip_header=8,skip_footer=153)
StaticData2=np.genfromtxt('S826_Ostavan.air',skip_header=8,skip_footer=29)

Static_AOA1=StaticData1[:,0]
Static_CL1=StaticData1[:,1]
Static_CD1=StaticData1[:,2]

Static_AOA2=StaticData2[:,0]
Static_CL2=StaticData2[:,1]
Static_CD2=StaticData2[:,2]

Static_CN1=Static_CL1*np.cos(Static_AOA1)+Static_CD1*np.sin(Static_AOA1)
Static_CN2=Static_CL2*np.cos(Static_AOA2)+Static_CD2*np.sin(Static_AOA2)

Dynamic_AOA=[];
Dynamic_CN=[];
time=[];

for i in range (1,26):
    DynamicData=np.genfromtxt('ALM/'+str(i)+'/Sheng.load',skip_header=1,delimiter=',')
    Dynamic_AOA.append(DynamicData[11,7])
    Dynamic_CN.append(DynamicData[11,15])

plt.figure(1)
plt.plot(Static_AOA1,Static_CL1,'-ko')
plt.plot(Static_AOA2,Static_CL2,'-ks')
plt.plot(Dynamic_AOA,Dynamic_CN,'-r+')
plt.xlim(-5.,35)
plt.ylim(-0.2,1.5)
plt.show()
