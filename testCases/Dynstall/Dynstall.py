import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from math import pi, exp, cos, sin, sqrt

conrad=pi/180.
condeg=180./pi

def read_data(fileID):
	# Read from .dat file, data should be formated in columns AOA, CL, CD, CM
	A=np.genfromtxt(fileID,delimiter='',skip_header=1)
	AOA=A[:,0]*conrad;
	CL=A[:,1];
	CD=A[:,2];
	CM=A[:,2];

	# Compute the Normal and Tangential forces
	CN=CL*np.cos(AOA)+CD*np.sin(AOA)
	CT=-CL*np.sin(AOA)+CD*np.cos(AOA)
	
	# Compute the static stall angle
	# determine the indices of the local maxima
	maxInd = argrelextrema(CL, np.greater)
	# get the actual values using these indices
	AlphaSS = 19.16*conrad #AOA[maxInd]
	
	# Compute the slope up to alphaSS
	Alpha0=AOA[np.nonzero(abs(CL)==min(abs(CL)))]
	Alpha11=7.*conrad
	Alpha1=17.4*conrad  # These two are hard-coded at the moment
	Alpha2=17.4*conrad
	# Region near 0

	Alphareg11=AOA[np.all([AOA > float(Alpha0), AOA < float(Alpha11)], axis=0)]
	CNreg11=CN[np.all([AOA > float(Alpha0), AOA < float(Alpha11)], axis=0)]
	
	# Do curve fitting from these data to find CNAlpha
	[A,B]=np.polyfit(Alphareg11,CNreg11,1)

	CNAlpha=A

	
	# Region 1
	Alphareg1=AOA[np.all([AOA > float(Alpha0), AOA < float(Alpha1)], axis=0)]
	CNreg1=CN[np.all([AOA > float(Alpha0), AOA < float(Alpha1)], axis=0)]
	
	# Region 2
	Alphareg2=AOA[AOA > float(Alpha2)]
	CNreg2=CN[AOA > float(Alpha2)]
	
	# Do the curve fitting for S1 and S2
	SMALL=1e-7
	# Start with region 1
	f1=(2.0*np.sqrt(CNreg1/(CNAlpha*(Alphareg1-Alpha0)))-1)**2.0
	f1[f1>=1]=1.0-SMALL
	f1[f1<=0.0]=SMALL
	
	[A1,B1]=np.polyfit(Alphareg1-Alpha1, np.log((1.0-f1)/0.4),1)
	S1=-1/B1
	
	f2=(2.0*np.sqrt(CNreg2/(CNAlpha*(Alphareg2-Alpha0)))-1)**2.0
	f2[f2>=1]=1.0-SMALL
	f2[f2<SMALL]=SMALL
	
	[A2,B2]=np.polyfit(Alpha1-Alphareg2, np.log((f2-0.02)/0.58),1)
	S2=-1/B2
	
	# Before you leave report the static data findings
	print "-------------------------------------------------------------"
	print " Static Data analysis results: "
	print " CNAlpha ratio : ", float(CNAlpha)
	print " Alpha0 -- zero Lift AOA : ", float(Alpha0)
	print " Alpha1 -- Maximum Lift AOA : ", float(Alpha1)
	print " AlphaSS -- Static Stall AOA : ", float(AlphaSS)
	print " "
	print " --- Fitting parameters for the the separation function -----"
	print " S1 = ", S1
	print " S2 = ", S2	
	print "-------------------------------------------------------------"

	return AOA,CL,CD,CN,CT,CM,Alpha0,AlphaSS,Alpha1,CNAlpha,S1,S2;

def separation_location(Alpha,alpha1,alpha2,S1,S2,S3):
	if(abs(Alpha)<alpha1):
		f=1-0.4*np.exp((abs(Alpha)-alpha1)/S1)
	else:
		f=0.02+0.58*np.exp((alpha1-abs(Alpha))/S2)

	if (f>1): 
		f=1.
	if (f<0):
		f=0. 
	return f;


# READ DATA
AOA_data,CL_data,CD_data,CN_data,CT_data,CM_data,alpha0,alphaSS,alpha1,CNAlpha,S1,S2=read_data('Data/S809.dat')

N=100;
Alpha=np.linspace(-10,50,N)*conrad;
f=np.linspace(1,1,N);
E0=0.16
eta=-0.6
alpha1=10.*conrad
alpha2=16.*conrad
S1=0.02
S2=0.03
S3=0.001
# Calculate the separation function based on Kirchoff flow equation
for i in range(N):
	f[i]=separation_location(Alpha[i],alpha1,alpha2,S1,S2,S3)

# Compute normal, tangential and pitching moment
#CN=CNAlpha*(Alpha-alpha0)*((1.0+np.sqrt(f))/2.0)**2.0
#CT=eta*CNAlpha*(Alpha-alpha0)**2.0*(np.sqrt(f)-E0)

#plt.figure(1)
#plt.plot(AOA_data*condeg,CN_data,'ko',markevery=1,label='Cn exp.')
#plt.plot(Alpha*condeg,CN,'-k',label='Cn recons')
#plt.plot(AOA_data*condeg,CT_data,'ro',markevery=1,label='Ct exp')
#plt.plot(Alpha*condeg,CT,'-r',label='Cc recons')
#plt.rc('text',usetex=True)
#plt.rc('font',family='serif')
#plt.ylabel(r"$C_n, C_c$",fontsize=12)
#plt.xlabel(r"$\alpha$ [deg]",fontsize=12)
#plt.legend()
# STARTING THE DYNAMIC STALL MODEL
# PARAMETERS
Nsteps=2000 #Number of time steps
U=1.;
chord=0.14;
dt=0.0005;
alpha_init=17.5
alpha_amp=10.
k=0.05
ds=2.0*U*dt/chord
omega=2.*k*U/chord
speed_of_sound=343.
Tp=1.7
Tf=3.0
TAlpha=1.
alphaDS0DiffDeg=0.
r0=0.01
Tv=7.0
Tvl=9.0
B1=0.5
B2=0.5
alphaSS=19.16*conrad
CD0=0.01

# Init Variables
alpha=0.
alpha_prev=0.
ds=0.
lambdaL=0.
lambdaL_prev=0.
H=0.
H_prev=0.
lambdaM=0.
lambdaM_prev=0.
J=0.
J_prev=0.
CNP=0.
CNP_prev=0.
DP=0.
DP_prev=0.
CNprime=0.
deltaalpha=0.
Dalpha=0.
Dalpha_prev=0.
StallFlag=False
StallFlag_prev=False
tau=0.
tau_prev=0.
DF=0.
DF_prev=0.
fprime=0.
fprime_prev=0.
fDoublePrime=0.
fDoublePrime_prev=0.
alphaPrime=0.
alphaPrime_prev=0.
X=0.
Y=0.
Z=0.
X_prev=0.
Y_prev=0.
Z_prev=0.
etaL=0.
etaL_prev=0.
#Output variables
AOADyn=[];
CNDyn=[];
CTDyn=[];
CLDyn=[];
CDDyn=[];
CMDyn=[];
print '--------------------------------------'
print 'Starting the dynamic stall simulations'
print '--------------------------------------' 



for i in range(Nsteps):
	t=i*dt
	alpha=alpha_init+alpha_amp*sin(omega*t)
	# Convert from degrees to rads
	alpha=alpha*conrad
	mach=U/speed_of_sound
	deltaalpha=alpha-alpha_prev
	ds=2.0*U*dt/chord
	s=i*ds
	##################################
	# Starting with the attached flow
	##################################

	beta=1-mach**2
        A1=0.165
	A2=0.335
	A3=0.5
	T1=20
	T2=4.5
	T3=1.25*mach
		
	X=X_prev*exp(-beta*ds/T1)+A1*(etaL-etaL_prev)*exp(-beta*ds/(2.0*T1))
	Y=Y_prev*exp(-beta*ds/T2)+A2*(etaL-etaL_prev)*exp(-beta*ds/(2.0*T2))
	Z=Z_prev*exp(-beta*ds/T3)+A3*(etaL-etaL_prev)*exp(-beta*ds/(2.0*T3))
	
	alphaEquiv=alpha#*(1.-X-Y-Z)

	CNC=CNAlpha*(alphaEquiv-alpha0)
       
	# Calculate the impulsive normal force coefficient
        lambdaL=(pi/4.0)*(alpha+chord/(4.0*U)*deltaalpha/dt)

        TI=chord/speed_of_sound*(1.0+3.0*mach)/4.0

        H=H_prev*exp(-dt/TI)+(lambdaL-lambdaL_prev)*exp(-dt/(2.0*TI))

        CNI=4.0/mach*H

        # Calculate the impulsive moment coefficient
        lambdaM=3*pi/16.0*(alpha+chord/(4.0*U)*deltaalpha/dt)+pi/16*chord/U*(deltaalpha/dt)

        J=J_prev*exp(-dt/TI)+(lambdaM-lambdaM_prev)*exp(-dt/(2.0*TI))

        CMI=-4.0/mach*J	
	
    	# Calculate lagged angle of attack
        Dalpha=Dalpha_prev*exp(-ds/TAlpha)+(alpha-alpha_prev)*exp(-ds/(2.0*TAlpha)) 
	alphaPrime=alpha-Dalpha
    	
	# Calculate reduced pitch rate
	r=deltaalpha/dt*chord/(2.*U)*conrad

    	# Calculate alphaDS0
    	dAlphaDS=alphaDS0DiffDeg*conrad
    	alphaDS0=alphaSS + dAlphaDS
    
    	if (abs(r)>=r0):
        	alphaCrit=alphaDS0
		Dalpha1=alphaDS0-alphaSS
    	else:
        	alphaCrit=alphaSS+(alphaDS0-alphaSS)*abs(r)/r0
		Dalpha1=(alphaDS0-alphaSS)*abs(r)/r0 

    	if(abs(alphaPrime)>alphaCrit):
        	StallFlag=True
		print 'Airfoil is stalled'
	
	##################################
	# Starting with the seperated flow
	##################################
	# Calculate trailing-edge separation point
	fprime=separation_location(alphaPrime,alpha1,alpha2,S1,S2,S3)
	
	# Calculate vortex tracking time
	if(StallFlag_prev==False): 
		tau=0.0
	else:
	    	tau=tau_prev+ds
	
	if(tau>Tf):
	    	TF=0.5*Tf
	else:
	    	TF=Tf
	
	# Calculate dynamic separation point
	DF=DF_prev*exp(-ds/TF)+(fprime-fprime_prev)*exp(-ds/(2.*TF))
	fDoublePrime=fprime-DF#*exp(-ds/TF)+(fprime-fprime_prev)*exp(-ds/(2.*TF)) 

	# Calculate vortex modulation parameter
	if(tau>=0. and tau<=Tv):
	    	Vx=sin(pi*tau/(2.0*Tv))**1.5
	elif (tau>Tv): 
	    	Vx=cos(pi*(tau-Tv)/Tvl)**2.0
	
	if (tau>Tvl):
	    	Vx=0.0

	# Calculate static trailing-edge separation point
	floc=separation_location(alpha,alpha1,alpha2,S1,S2,S3)

	# Calculate normal force coefficient including dynamic separation point
	CNF=CNAlpha*(alphaEquiv-alpha0)*((1.0+sqrt(fDoublePrime))/2.0)**2.+CNI
	
	print alphaPrime, fDoublePrime, CNF
	# Calculate tangential force coefficient
	CT=eta*CNAlpha*(alphaEquiv-alpha0)**2.*(sqrt(fDoublePrime)-E0)
	
	# Calculate vortex lift contribution
	CNV=B1*(fDoublePrime-floc)*Vx

	# Total normal force coefficient is the combination of that from 
	# circulatory effects, impulsive effects, dynamic separation, and vortex
	# lift
	CN=CNF+CNV

	# Calculate moment coefficient
	#cmf=(ds%K0+ds%K1*(1-ds%fDoublePrime)+ds%K2*sin(pi*ds%fDoublePrime**2))*ds%CNC
	#! + moment coefficient at Zero lift angle of attack
	#cmv = ds%B2*(1.0-cos(pi*ds%tau/ds%Tv))*ds%CNV
	#ds%CM=cmf+cmv+ds%CMI
	
	############################
	# Calculate final loads
	############################
	AOADyn.append(alpha*condeg);
	CNDyn.append(CN);
	CTDyn.append(CT);
	CLDyn.append(CN*cos(alpha)-CT*sin(alpha))
	CDDyn.append(CN*sin(alpha)+CT*cos(alpha)+CD0)
	print ''
	print ' Current time :', i*dt
	print ' Current AOA  :', alpha*condeg
	print ' Reduced pitch rate :', r
	print ' Separation location :', fprime
	print ' non-dimensional parameter : ', ds
	print 'TAlpha :', TAlpha
	print ''
	
	############################
	# Update previous variables 
	############################
	X_prev=X
	Y_prev=Y
	Z_prev=Z
	etaL_prev=etaL
	Alpha_prev=alpha
	DP_prev=DP
	CNP_prev=CNP
	DF_prev=DF
	fprime_prev=fprime
	CNV_prev=CNV
	StallFlag_prev=StallFlag
	tau_prev=tau
	H_prev=H
	lambdaL_prev=lambdaL
	J_prev=J
	lambdaM_prev=lambdaM
	alphaPrime_prev=alphaPrime
	fDoublePrime_prev=fDoublePrime
	Dalpha_prev=Dalpha



DynLiftData=np.genfromtxt('Data//DynLiftS809.dat',delimiter='',skip_header=1)
AOAL=DynLiftData[:,0]
Lift=DynLiftData[:,1]
DynDragData=np.genfromtxt('Data/DynDragS809.dat',delimiter='',skip_header=1)
AOAD=DynDragData[:,0]
Drag=DynDragData[:,1]
DynMomData=np.genfromtxt('Data/DynMomS809.dat',delimiter='',skip_header=1)
AOAM=DynMomData[:,0]
Mom=DynMomData[:,1]

Dynamic_AOA=[];
Dynamic_CL=[];
Dynamic_CD=[];
Dynamic_CM=[];

for i in range (1,18):
    DynamicData=np.genfromtxt('ALM/'+str(i)+'/Sheng.load',skip_header=1,delimiter=',')
    Dynamic_AOA.append(DynamicData[11,6]-90)
    Dynamic_CL.append(DynamicData[11,12])
    Dynamic_CD.append(DynamicData[11,13])
    Dynamic_CM.append(DynamicData[11,14])

plt.rc('text',usetex=True)
plt.rc('font',family='serif')

plt.figure(2)
plt.plot(Dynamic_AOA,Dynamic_CL,'-ro',markevery=1,label='CL dyn. recon.')
plt.plot(AOAL,Lift,'k*',label='DS S809 exp.')
plt.ylabel(r"$C_L$",fontsize=12)
plt.xlabel(r"$\alpha$ [deg]",fontsize=12)
plt.legend()
plt.grid()
plt.xlim(0,30)
#plt.ylim(0.2,1.8)
plt.savefig('CLDyn.pdf')

plt.figure(3)
#plt.plot(AOADyn[100:2000],CLDyn[100:2000],'-k',markevery=5,label='Cn dyn. recon.')
plt.plot(Dynamic_AOA,Dynamic_CD,'-ro',markevery=1,label='CL dyn. recon.')
plt.plot(AOAD,Drag,'k*',label='DS S809 exp.')
plt.ylabel(r"$C_D$",fontsize=12)
plt.xlabel(r"$\alpha$ [deg]",fontsize=12)
plt.legend()
plt.xlim(0,30)
plt.grid()
plt.savefig('CDDyn.pdf')

plt.figure(4)
#plt.plot(AOADyn[100:2000],CLDyn[100:2000],'-k',markevery=5,label='Cn dyn. recon.')
plt.plot(AOAM,Mom,'k*',label='DS S809 exp.')
plt.ylabel(r"$C_M$",fontsize=12)
plt.xlabel(r"$\alpha$ [deg]",fontsize=12)
plt.legend()
plt.xlim(0,30)
plt.grid()
plt.savefig('CMDyn.pdf')

