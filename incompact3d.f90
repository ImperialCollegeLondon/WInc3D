!###############################################################################
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

PROGRAM incompact3d

USE decomp_2d
USE decomp_2d_poisson
use decomp_2d_io
USE param
USE var
USE MPI
USE IBM
USE derivX
USE derivZ

!>>GD INTRODUCE THE actuator_line_modules
use actuator_line_model 
use actuator_line_source

implicit none

integer :: code,nlock,i,j,k,ii,bcx,bcy,bcz,fh,ierror
real(mytype) :: x,y,z,tmp1
double precision :: t1,t2
integer :: ErrFlag, nargin, FNLength, status, DecInd
logical :: back
character(len=80) :: InputFN, FNBase
character(len=20) :: filename

TYPE(DECOMP_INFO) :: phG,ph1,ph2,ph3,ph4

CALL MPI_INIT(code)

!==========================================================================
! Handle Input file
nargin=command_argument_count()
if (nargin <1) then
    write(6,*) 'Please call the program with the name of the input file on the command line Ex. incompact3d input.in'
    stop
endif

call get_command_argument(1,InputFN,FNLength,status)
back=.true.
FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
DecInd=index(FNBase,'.',back)
if (DecInd >1) then
    FNBase=FNBase(1:(DecInd-1))
end if
!===========================================================================

call parameter(InputFN)

! Decompose the problem into rows and columns 
call decomp_2d_init(nx,ny,nz,p_row,p_col)

call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)
call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)

call init_variables

! Print out the parameters  
if (nrank==0) then
print *,'==========================================================='
print *,'==========================================================='
print *,'==========================================================='
print *,'======================Incompact3d=========================='
print *,'===Copyright (c) 2012 Eric Lamballais and Sylvain Laizet==='
print *,'eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com'
print *,'Imperial College London Version: '
print *,'Currently in development by G. Deskos'
print *,'==========================================================='
print *,'==========================================================='
print *,'==========================================================='
print *,''
print *,''
print *,''
if (itype.eq.1) print *,'Fluid flow problem : Constant flow field/wake in a confined domain'
if (itype.eq.2) print *,'Fluid flow problem : Smooth turbulent channel flow -- periodic domain'
if (itype.eq.3) print *,'Fluid flow problem : Various configurations: Flow field is imported from a precursor simulation through planes'
if (itype.eq.4) print *,'Fluid flow problem : Smooth/rough developing boundary layer'
if (itype.eq.5) print *,'Fluid flow problem : Not defined yet'
if (itype.eq.6) print *,'Fluid flow problem : Taylor Green vortices'
if (itype.eq.7) print *,'Fluid flow problem : 3D cavity flow'
if (itype.eq.8) print *,'Fluid flow problem : Atmospheric Boundary Layer'
if (itype.eq.9) print *,'Fluid flow problem : Not defined yet'
write(*,1101) nx,ny,nz
write(*,1103) xlx,yly,zlz 
write(*,1102) nclx,ncly,nclz 
write(*,1104) u1,u2 
write(*,1105) re
write(*,1106) dt
if (nscheme.eq.1) print *,'Temporal scheme   : Adams-Bashforth 2'
if (nscheme.eq.2) print *,'Temporal scheme   : Runge-Kutta 3'
if (nscheme.eq.3) print *,'Temporal scheme   : Runge-Kutta 4'
if (nscheme.eq.2) print *,'Temporal scheme   : Adams-Bashforth 4'

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

! ==========================================================================================================
! Initialize the schemes (Currently deciding between the optimal sixth-order Pade and the sixth-order iSVV
! ==========================================================================================================
if (jLES==0) then
    call schemes_dns()
    if(nrank==0) then
        write(*,*) 'DNS'
    endif
else if (jLES==1) then
    call schemes_iles()
    if(nrank==0) then
        write(*,*) 'Implicit LES with xxnu = 1 / ', rxxnu
    endif
else if (jLES==2.OR.jLES==3.OR.jLES==4.OR.jLES==5) then 
    call init_explicit_les() 
    call schemes_dns()
endif
! ================================================
! Dealing with pressure
! ================================================
if (nclx==0) then
   bcx=0
else
   bcx=1
endif
if (ncly==0) then
   bcy=0
else
   bcy=1
endif
if (nclz==0) then
   bcz=0
else
   bcz=1
endif

call decomp_2d_poisson_init(bcx,bcy,bcz)

call decomp_info_init(nxm,nym,nzm,phG)

! ======================================================
! Initialise flow in the domain -- init or restart 
! ======================================================
if (ilit==0) call init(ux1,uy1,uz1,ep1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)  
if (ilit==1) call restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,&
        px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,0)
call test_speed_min_max(ux1,uy1,uz1)
if (iscalar==1) call test_scalar_min_max(phi1)

!array for stat to zero
umean=0.;vmean=0.;wmean=0.
uumean=0.;vvmean=0.;wwmean=0.
uvmean=0.;uwmean=0.;vwmean=0.
phimean=0.;phiphimean=0.

t1 = MPI_WTIME()

!div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
call decomp_info_init(nxm, nym, nzm, ph1)
call decomp_info_init(nxm, ny, nz, ph4)

!gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
call decomp_info_init(nxm, ny, nz, ph2)  
call decomp_info_init(nxm, nym, nz, ph3) 

! Initialise the Turbine Model
if (ialm==1) then
  call actuator_line_model_init(Nturbines,Nactuatorlines,TurbinesPath,ActuatorlinesPath,dt)  
  call initialize_actuator_source 
endif

! Initialise the Probe inside the domain
if (iprobe==1) then
    call init_probe
endif

! Initialise inflow file
if (iin==3) then
    call read_inflow(ux_inflow,uy_inflow,uz_inflow)
endif

do itime=ifirst,ilast
   t=(itime-1)*dt

   if (nrank==0) then
      write(*,*) '========================================'
      write(*,1001) itime,t
1001  format(' Time step =',i7,', Time unit =',F9.3)
      write(*,*) '========================================'
   endif

   if(ialm==1) then
      if (nrank==0) then
          write(6,*) '' 
          write(6,*) 'Unsteady ACtuator Line Model INFO:'
      endif
      call Compute_Momentum_Source_Term_pointwise            
      call actuator_line_model_update(t,dt)
      if (nrank==0) then
          write(6,*) '' 
      endif
   endif
   
   do itr=1,iadvance_time

      if (nclx.eq.2) then
         call inflow (ux1,uy1,uz1,phi1) !X PENCILS
         call outflow(ux1,uy1,uz1,phi1) !X PENCILS 
      endif 
        
      call convdiff(ux1,uy1,uz1,phi1,uxt,uyt,uzt,ep1,divdiva,curldiva,ta1,tb1,tc1,&
      td1,te1,tf1,tg1,th1,ti1,di1,ux2,uy2,uz2,phi2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,&
      ti2,tj2,di2,ux3,uy3,uz3,phi3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,nut1,ucx1,&
      ucy1,ucz1,tmean,sgszmean,sgsxmean,sgsymean)

      ! Potential Temperature -- to be computed after the convdiff
      if (ibuoyancy==1) then
          ! Transport equation for potential temperature
          call PotentialTemperature(ux1,uy1,uz1,phi1,phis1,phiss1,di1,tg1,th1,ti1,td1,&
              uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,ep1)  
      endif
      
      !X PENCILS
      call intt (ux1,uy1,uz1,gx1,gy1,gz1,hx1,hy1,hz1,ta1,tb1,tc1) 
      
      call pre_correc(ux1,uy1,uz1,phi1,ta1)
      
      if (ivirt==1) then !solid body old school
         !we are in X-pencil
         call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
         call body(ux1,uy1,uz1,ep1)
         call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
      endif

      !X-->Y-->Z
      call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
           td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,pp3,&
           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,1)       

      !POISSON Z-->Z 
      call decomp_2d_poisson_stg(pp3,bcx,bcy,bcz)

      !Z-->Y-->X
      call gradp(px1,py1,pz1,di1,td2,tf2,ta2,tb2,tc2,di2,&
           ta3,tc3,di3,pp3,nxmsize,nymsize,nzmsize,ph2,ph3)

      !X PENCILS
      call corgp(ux1,ux2,uy1,uz1,px1,py1,pz1) 
      
     !does not matter -->output=DIV U=0 (in dv3)
      call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
           td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,dv3,&
           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,2)


      call test_speed_min_max(ux1,uy1,uz1)
      if (iscalar==1) call test_scalar_min_max(phi1)

   enddo
   
   if (t>=spinup_time) then
       call STATISTIC(ux1,uy1,uz1,phi1,ta1,umean,vmean,wmean,phimean,uumean,vvmean,wwmean,&
           uvmean,uwmean,vwmean,phiphimean,tmean)
        
       if(ioutflow==1) call write_outflow(ux1,uy1,uz1) 

       if(ialm==1) call actuator_line_statistics()

       if(iprobe==1) then
           if (mod(itime,nsampling)==0) then
               call probe(ux1,uy1,uz1,phi1)
               call write_probe(itime/nsampling) 
           endif
       endif
   endif

   if (mod(itime,isave)==0) call restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,&
        px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,1)
     
   if (mod(itime,imodulo)==0) then
      call VISU_INSTA(ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
           ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
           ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
      call VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
           ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu) 
   endif

   if (ialm==1) then
    if (nrank==0) then
       call actuator_line_model_write_output(itime/imodulo) ! Write the Rotor output
    end if
   endif 
enddo

t2=MPI_WTIME()-t1
call MPI_ALLREDUCE(t2,t1,1,MPI_REAL8,MPI_SUM, &
                   MPI_COMM_WORLD,code)
if (nrank==0) print *,'time per time_step: ', &
     t1/float(nproc)/(ilast-ifirst+1),' seconds'
if (nrank==0) print *,'simulation with nx*ny*nz=',nx,ny,nz,'mesh nodes'
if (nrank==0) print *,'Mapping p_row*p_col=',p_row,p_col


!call decomp_2d_poisson_finalize
call decomp_2d_finalize
CALL MPI_FINALIZE(code)

end PROGRAM incompact3d
