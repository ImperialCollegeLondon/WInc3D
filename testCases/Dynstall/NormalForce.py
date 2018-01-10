import numpy as np
import matplotlib.pyplot as plt
from math import pi

StaticData=np.genfromtxt('S809.air',skip_header=8,skip_footer=153)

Static_AOA=StaticData[:,0]
Static_CL=StaticData[:,1]
Static_CD=StaticData[:,2]

Static_CN=Static_CL*np.cos(Static_AOA)+Static_CD*np.sin(Static_AOA)

Dynamic_AOA=[];
Dynamic_CN=[];
time=[];

for i in range (0,101):
    DynamicData=np.genfromtxt('ALM/'+str(i)+'/Sheng.load',skip_header=1,delimiter=',')
    Dynamic_AOA.append(DynamicData[11,7])
    Dynamic_CN.append(DynamicData[11,15])

plt.figure(1)
plt.plot(Static_AOA,Static_CL,'-ko')
plt.plot(Dynamic_AOA,Dynamic_CN,'-r+')
plt.xlim(-5.,35)
plt.ylim(-0.2,1.5)
plt.show()
