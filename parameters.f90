!################################################################################
!This file is part of Incompact3d.
!
!Incompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Incompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Incompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Incompact3d in your publications and 
!    presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for 
!    incompressible flows: a simple and efficient method with the quasi-spectral 
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence 
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical 
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################
!
!********************************************************************
!
subroutine parameter(InputFN)
!
!********************************************************************
  
USE IBM 
USE param
USE decomp_2d

implicit none

character(len=*),intent(in):: InputFN 
real(mytype) :: theta, cfl,cf2
integer :: longueur ,impi,j
character :: a*80

! Have you heard of NAMELISTs ?
NAMELIST/FlowParam/itype,iin,NEddies,sem_file,xlx,yly,zlz,re,sc,u1,u2,noise,noise1,itripping,ibuoyancy,icoriolis,Pr,TempRef,CoriolisFreq 
NAMELIST/NumConfig/nx,ny,nz,p_row,p_col,nclx,ncly,nclz,TurbRadius,ifirst,ilast,nscheme,dt,istret, &
    beta,iskew,iscalar,jles,FSGS,jadv,smagcst,SmagWallDamp,nSmag,iwallmodel,walecst,rxxnu,cnu,dynhypvisc  
NAMELIST/FileParam/ilit,isave,imodulo,ioutflow,InflowPath, NTimeSteps
NAMELIST/IBMParam/ivirt,ibmshape,cex,cey,cez,ra
NAMELIST/ALMParam/ialm,ialmrestart, ialmoutput,NTurbines,TurbinesPath,NActuatorlines,ActuatorlinesPath,eps_factor
NAMELIST/ADMParam/iadm,Ndiscs,ADMcoords,iverifyadm,iadmmode,CT,aind,fileADM
NAMELIST/StatParam/spinup_time,nstat,nvisu,iprobe,Probelistfile,nsampling, y_loc_pencil,& 
                  z_loc_pencil,isnapshot,simin,simax,sjmin,sjmax,skmin,skmax,sfreq 
NAMELIST/ABLParam/iabl,z_zero,k_roughness,PsiM,ustar,IPressureGradient,Ug,dBL,idampingzone,ifringeregion,FLS,FLE,Imassconserve 

#ifdef DOUBLE_PREC 
pi=dacos(-1.d0) 
#else
pi=acos(-1.)
#endif

twopi=2.*pi

! IF variables are not set we will need to give them some default values
nx=64
ny=64
nz=64
p_col=2
p_row=2
nstat=1
nvisu=1
xlx=1.0
yly=1.0
zlz=1.0
re=1000
nclx=0
ncly=2
nclz=0
itype=1
iin=1
TurbRadius=0.0
ifirst=1
ilast=1000
nscheme=1
istret=0
beta=0.28
iskew=1
iscalar=0
jles=0
FSGS=2.0
smagcst=0.1
walecst=0.5
ilit=0
isave=100
imodulo=100
ivirt=0
ialm=0
ialmoutput=50
ialmrestart=0.
eps_factor=2.0
rxxnu=1.0
cnu=0.44
dynhypvisc=0
nSmag=1.0
SmagWallDamp=0
dBL=500 !delta of the boundary layer
UG=[10.0d0,0.0d0,0.0d0]
IPressureGradient=0
ioutflow=0
InflowPath='./'
idampingzone=0
ifringeregion=0
Imassconserve=0
itripping=0
iadm=0
Ndiscs=0
iverifyadm=0
iadmmode=0
CT=0.75
aind=0.25
fileADM='/./'
simin=1
simax=nx
sjmin=1
sjmax=ny
skmin=1
skmax=nz
sfreq=100

! READ PARAMETERS FROM FILE
open(10,file=InputFN) 
read(10,nml=FlowParam)
read(10,nml=NumConfig)
read(10,nml=StatParam)
read(10,nml=FileParam)
read(10,nml=IBMParam)
read(10,nml=ALMParam)
read(10,nml=ADMParam)
read(10,nml=ABLParam)

close(10) 
 
call init_module_parameters()

return  
end subroutine parameter



