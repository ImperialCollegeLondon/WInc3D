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
NAMELIST/FlowParam/itype,iin,xlx,yly,zlz,re,sc,u1,u2,noise,noise1,ibuoyancy,icoriolis,Pr,TempRef,CoriolisFreq
NAMELIST/NumConfig/nx,ny,nz,p_row,p_col,nclx,ncly,nclz,TurbRadius,ifirst,ilast,nscheme,dt,istret, &
    beta,iskew,iscalar,jles,FSGS,jadv,smagcst,SmagWallDamp,nSmag,walecst,rxxnu 
NAMELIST/FileParam/ilit,isave,imodulo
NAMELIST/IBMParam/ivirt,ibmshape,cex,cey,cez,ra
NAMELIST/ALMParam/ialm,NTurbines,TurbinesPath,NActuatorlines,ActuatorlinesPath,eps_factor
NAMELIST/StatParam/spinup_time,nstat,nvisu,iprobe,Probelistfile,nsampling 
NAMELIST/ABLParam/iabl,z_zero,k_roughness,PsiM,ustar,IPressureGradient,Ug,epsilon_pert 
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
eps_factor=2.0
rxxnu=3.0
nSmag=1.0
SmagWallDamp=0
IPressureGradient=0

! READ PARAMETERS FROM FILE
open(10,file=InputFN) 
read(10,nml=FlowParam)
read(10,nml=NumConfig)
read(10,nml=StatParam)
read(10,nml=FileParam)
read(10,nml=IBMParam)
read(10,nml=ALMParam)
read(10,nml=ABLParam)

close(10) 

if (nrank==0) then
print *,'==========================================================='
print *,'==========================================================='
print *,'==========================================================='
print *,'======================Incompact3d=========================='
print *,'===Copyright (c) 2012 Eric Lamballais and Sylvain Laizet==='
print *,'eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com'
print *,'==========================================================='
print *,'==========================================================='
print *,'==========================================================='
print *,''
print *,''
print *,''
if (itype.eq.1) print *,'Constant flow field'
if (itype.eq.2) print *,'Channel flow'
if (itype.eq.3) print *,'Wake flow'
if (itype.eq.4) print *,'Mixing layer with splitter plate'
if (itype.eq.5) print *,'Channel flow'
if (itype.eq.6) print *,'Taylor Green vortices'
if (itype.eq.7) print *,'Cavity flow'
if (itype.eq.8) print *,'Atmospheric boundary layer'
if (itype.eq.9) print *,'Water tank'
write(*,1101) nx,ny,nz
write(*,1103) xlx,yly,zlz 
write(*,1102) nclx,ncly,nclz 
write(*,1104) u1,u2 
write(*,1105) re
write(*,1106) dt
if (nscheme.eq.1) print *,'Temporal scheme   : Adams-bashforth 2'
if (nscheme.eq.2) print *,'Temporal scheme   : Runge-Kutta 3'
if (nscheme.eq.3) print *,'Temporal scheme   : Runge-Kutta 4'
if (iscalar.eq.0) print *,'Passive scalar    : off'
if (iscalar.eq.1) then
   print *,'Passive scalar : on'
   write (*,1113) sc
endif
if (ivirt.eq.0) print *,'Immersed boundary : off'
if (ivirt.eq.1) then
   print *,'Immersed boundary : on old school'
   write(*,1107) cex,cey,cez
   write(*,1110) ra
endif
if (ivirt.eq.2) then
   print *,'Immersed boundary : on with Lagrangian Poly'
endif


 1101 format(' Spatial Resolution: (nx,ny,nz)=(',I4,',',I4,',',I4,')')
 1102 format(' Boundary condition: (nclx,ncly,nclz)=(',I1,',',I1,',',I1,')')
 1103 format(' Domain dimension  : (lx,ly,lz)=(',F6.1,',',F6.1,',',F6.1,')')
 1104 format(' High and low speed: u1=',F6.2,' and u2=',F6.2)
 1105 format(' Reynolds number Re: ',F15.8)
 1106 format(' Time step dt      : ',F15.8)
 1107 format(' Object centred at : (',F6.2,',',F6.2,',',F6.2,')')
 1110 format(' Object length     : ',F6.2)
 1113 format(' Schmidt number    : ',F6.2)
endif
 
call init_module_parameters()

return  
end subroutine parameter



