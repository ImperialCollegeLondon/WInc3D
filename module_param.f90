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

module param 

  use decomp_2d, only : mytype

! Boundary conditions : ncl = 2 --> Dirichlet
! Boundary conditions : ncl = 1 --> Free-slip
! Boundary conditions : ncl = 0 --> Periodic
! l: power of 2,3,4,5 and 6
! if ncl = 1 or 2, --> n  = 2l+ 1 
!                  --> nm = n - 1 
!                  --> m  = n + 1
! If ncl = 0,      --> n  = 2*l
!                  --> nm = n  
!                  --> m  = n + 2
!nstat = size arrays for statistic collection
!2-->every 2 mesh nodes
!4-->every 4 mesh nodes
!nvisu = size for visualization collection
integer,save :: nx,ny,nz
integer,save :: nstat,nvisu
integer,save :: p_row,p_col
integer,save :: nxm,nym,nzm
!end module variables
  
integer :: nclx,ncly,nclz
integer :: ifft, ivirt,istret,iforc_entree,iturb, ialm, Nturbines, NActuatorlines
integer :: itype, iskew, iin, nscheme, ifirst, ilast, iles, jLES, jADV
integer :: isave,ilit,srestart,idebmod, imodulo, idemarre, icommence, irecord
integer :: iscalar,ilag,npif,izap
character :: inflow_file*80, sem_file*80
integer :: iprobe, Nprobes, Nsampling, ioutflow, iinflow, OutflowOnsetIndex, NTimeSteps
integer :: y_loc_pencil(4), z_loc_pencil(4)
integer :: NEddies  !For syntetic eddy method
integer :: iabl, ibmshape, SmagWallDamp, iwallmodel
character :: Probelistfile*80, inflowdir*80
integer :: nxboite, istat,iread,iadvance_time, ibuoyancy, icoriolis
real(mytype) :: xlx,yly,zlz,dx,dy,dz,dx2,dy2,dz2
real(mytype) :: dt,xnu,noise,noise1,pi,twopi,u1,u2,re,sc,Pr,TempRef,CoriolisFreq
real(mytype) :: t,xxk1,xxk2, spinup_time
real(mytype) :: smagcst, nSmag, walecst,dys, FSGS, rxxnu, cnu, dynhypvisc
real(mytype) :: eps_factor ! Smoothing factor 
real(mytype) :: TurbRadius,z_zero,k_roughness,PsiM,ustar,u_shear, dBL
integer :: IPressureGradient,idampingzone, Imassconserve
real(mytype),dimension(3) :: UG
integer :: itr,itime
character :: dirname*80
character :: filesauve*80, filenoise*80, &
     nchamp*80,filepath*80, fileturb*80, filevisu*80 
character,dimension(100) :: turbinesPath*80, ActuatorlinesPath*80 ! Assign a maximum number of 100 turbines and alms
real(mytype), dimension(5) :: adt,bdt,cdt,gdt
! Actuator disc model
integer :: iadm ! Actuator disc flag
integer :: Ndiscs ! number of actuator discs
character(len=100) :: admCoords
integer :: iadmmode ! 0: constnat thrust, 1: from a list
real(mytype) :: CT, aind
character(len=100) :: fileADM

!module filter
real(mytype),dimension(200) :: idata
real(mytype), save, allocatable, dimension(:) :: fiffx, fifcx, fifbx, fisfx, fiscx, fisbx,fifsx,fifwx,fissx,fiswx
real(mytype), save, allocatable, dimension(:) :: fiffxp,fifsxp,fifwxp,fisfxp,fissxp,fiswxp
real(mytype), save, allocatable, dimension(:) :: fiffy, fifcy, fifby, fisfy, fiscy, fisby,fifsy,fifwy,fissy,fiswy
real(mytype), save, allocatable, dimension(:) :: fiffyp,fifsyp,fifwyp,fisfyp,fissyp,fiswyp
real(mytype), save, allocatable, dimension(:) :: fiffz, fifcz, fifbz, fisfz, fiscz, fisbz,fifsz,fifwz,fissz,fiswz
real(mytype), save, allocatable, dimension(:) :: fiffzp,fifszp,fifwzp,fisfzp,fisszp,fiswzp
real(mytype), save, allocatable, dimension(:,:) :: fisx,fivx
real(mytype), save, allocatable, dimension(:,:) :: fisy,fivy
real(mytype), save, allocatable, dimension(:,:) :: fisz,fivz

!module derivative
real(mytype), save, allocatable, dimension(:) :: ffx,fcx,fbx,sfx,scx,sbx,fsx,fwx,ssx,swx
real(mytype), save, allocatable, dimension(:) :: ffxp,fsxp,fwxp,sfxp,ssxp,swxp
real(mytype), save, allocatable, dimension(:) :: ffy,fcy,fby,sfy,scy,sby,fsy,fwy,ssy,swy
real(mytype), save, allocatable, dimension(:) :: ffyp,fsyp,fwyp,sfyp,ssyp,swyp
real(mytype), save, allocatable, dimension(:) :: ffz,fcz,fbz,sfz,scz,sbz,fsz,fwz,ssz,swz
real(mytype), save, allocatable, dimension(:) :: ffzp,fszp,fwzp,sfzp,sszp,swzp
real(mytype), save, allocatable, dimension(:,:) :: sx,vx
real(mytype), save, allocatable, dimension(:,:) :: sy,vy
real(mytype), save, allocatable, dimension(:,:) :: sz,vz

!module pressure
real(mytype), save, allocatable, dimension(:,:) :: dpdyx1,dpdyxn,dpdzx1,dpdzxn
real(mytype), save, allocatable, dimension(:,:) :: dpdxy1,dpdxyn,dpdzy1,dpdzyn
real(mytype), save, allocatable, dimension(:,:) :: dpdxz1,dpdxzn,dpdyz1,dpdyzn

!module solid_body
integer :: nxfin,nyfin,nzfin

!module inflow
real(mytype), save, allocatable, dimension(:,:) :: bxx1,bxy1,bxz1,bxxn,bxyn,bxzn,bxo,byo,bzo
real(mytype), save, allocatable, dimension(:,:) :: byx1,byy1,byz1,byxn,byyn,byzn
real(mytype), save, allocatable, dimension(:,:) :: bzx1,bzy1,bzz1,bzxn,bzyn,bzzn
  
!module tripping
integer ::  z_modes, nxt_itr, itrip
real(mytype) :: x0_tr, xs_tr, ys_tr, ts_tr, zs_param, zs_tr, randomseed, A_trip
real(mytype), allocatable, dimension(:) :: h_coeff, h_nxt,h_i

!module enspec
real(mytype), dimension(3,7000) :: uensp,vensp,wensp

!module derpres
real(mytype), save, allocatable, dimension(:) :: cfx6,ccx6,cbx6,cfxp6,ciwxp6,csxp6,&
     cwxp6,csx6,cwx6,cifx6,cicx6,cisx6   
real(mytype), save, allocatable, dimension(:) :: cibx6,cifxp6,cisxp6,ciwx6
real(mytype), save, allocatable, dimension(:) :: cfi6,cci6,cbi6,cfip6,csip6,cwip6,csi6,&
    cwi6,cifi6,cici6,cibi6,cifip6  
real(mytype), save, allocatable, dimension(:) :: cisip6,ciwip6,cisi6,ciwi6 
real(mytype),save, allocatable, dimension(:) :: cfy6,ccy6,cby6,cfyp6,csyp6,cwyp6,csy6 
real(mytype),save, allocatable, dimension(:) :: cwy6,cify6,cicy6,ciby6,cifyp6,cisyp6,&
     ciwyp6,cisy6,ciwy6 
real(mytype),save, allocatable, dimension(:) :: cfi6y,cci6y,cbi6y,cfip6y,csip6y,cwip6y,&
     csi6y,cwi6y,cifi6y,cici6y  
real(mytype),save, allocatable, dimension(:) :: cibi6y,cifip6y,cisip6y,ciwip6y,cisi6y,ciwi6y  
real(mytype),save, allocatable, dimension(:) :: cfz6,ccz6,cbz6,cfzp6,cszp6,cwzp6,csz6 
real(mytype),save, allocatable, dimension(:) :: cwz6,cifz6,cicz6,cibz6,cifzp6,ciszp6,&
     ciwzp6,cisz6,ciwz6 
real(mytype),save, allocatable, dimension(:) :: cfi6z,cci6z,cbi6z,cfip6z,csip6z,cwip6z,&
     csi6z,cwi6z,cifi6z,cici6z  
real(mytype),save, allocatable, dimension(:) :: cibi6z,cifip6z,cisip6z,ciwip6z,cisi6z,ciwi6z 

!module waves
complex(mytype), save, allocatable, dimension(:) :: zkz,zk2,ezs
complex(mytype), save, allocatable, dimension(:) :: yky,yk2,eys
complex(mytype), save, allocatable, dimension(:) :: xkx,xk2,exs

!module mesh
real(mytype), save, allocatable, dimension(:) :: ppy,pp2y,pp4y
real(mytype), save, allocatable, dimension(:) :: ppyi,pp2yi,pp4yi
real(mytype), save, allocatable, dimension(:) :: yp,ypi,del
real(mytype), save, allocatable, dimension(:) :: yeta,yetai
real(mytype) :: alpha,beta

contains
    
    subroutine init_module_parameters()
    
    use decomp_2d, only: nrank
    implicit none
    ! Allocate parameters
    integer :: j 

    ! Compute the sizing parameters 
        ! X-DIRECTION
    if (nclx==0) then
        if(mod(nx,2)==0) then 
            !if (nrank==0) write(*,*) 'nx ok !'
        else
            if (nrank==0) write(*,*) 'nx should be an even number!'
            stop
        endif
        dx=xlx/nx
        nxm=nx 
    else if (nclx==1 .or. nclx==2) then 
        if(mod(nx,2)==1) then 
            !if (nrank==0) write(*,*) 'nx ok !'
        else
            if (nrank==0) write(*,*) 'nx should be an odd number!'
            stop
        endif
        dx=xlx/(nx-1.)
        nxm=nx-1
    else
        if (nrank==0) write(*,*) 'please select between nclx=0,1 and 2'
        stop
    endif
        ! Y-DIRECTION 
    if (ncly==0) then 
        if(mod(ny,2)==0) then 
            !if (nrank==0) write(*,*) 'ny ok !'
        else
            if (nrank==0) write(*,*) 'ny should be an even number!'
            stop
        endif
        dy=yly/ny 
        nym=ny
    else if (ncly==1.or.ncly==2) then 
        if(mod(ny,2)==1) then 
            !if (nrank==0) write(*,*) 'ny ok !'
        else
            if (nrank==0) write(*,*) 'ny should be an odd number!'
            stop
        endif
        dy=yly/(ny-1.) 
        nym=ny-1
    else
        if (nrank==0) write(*,*) 'please select between ncly=0,1 and 2'
        stop
    endif
        ! Z-DIRECTION  
    if (nclz==0) then 
        if(mod(nz,2)==0) then 
            !if (nrank==0) write(*,*) 'nz ok !'
        else
            if (nrank==0) write(*,*) 'nz should be an even number!'
            stop
        endif
        dz=zlz/nz 
        nzm=nz
    else if (nclz==1.or.nclz==2) then 
        if(mod(nz,2)==1) then 
            !if (nrank==0) write(*,*) 'nz ok !'
        else
            if (nrank==0) write(*,*) 'nz should be an odd number!'
            stop
        endif
        dz=zlz/(nz-1.) 
        nzm=nz-1
    else
        if (nrank==0) write(*,*) 'please select between nclz=0,1 and 2'
        stop
    endif
    
    ! Compute the squares of the grid size
    dx2=dx*dx
    dy2=dy*dy
    dz2=dz*dz
    

    ! Allocate filter and interpolation parameters 
    ! X-direction
    !allocate(fifx(nx),ficx(nx),fibx(nx),fiffx(nx),fibbx(nx),fiz1x(nx),fiz2x(nx))
    !allocate(filax(nx,2),filaxp(nx,2))
    !allocate(fifxp(nx),ficxp(nx),fibxp(nx),fiffxp(nx),fibbxp(nx))
    !! Y-direction 
    !allocate(fify(ny),ficy(ny),fiby(ny),fiffy(ny),fibby(ny),fiz1y(ny),fiz2y(ny))
    !allocate(filay(ny,2),filayp(ny,2))
    !allocate(fifyp(ny),ficyp(ny),fibyp(ny),fiffyp(ny),fibbyp(ny))
    !! Z-direction  
    !allocate(fifz(nz),ficz(nz),fibz(nz),fiffz(nz),fibbz(nz),fiz1z(nz),fiz2z(nz))
    !allocate(filaz(nz,2),filazp(nz,2))
    !allocate(fifzp(nz),ficzp(nz),fibzp(nz),fiffzp(nz),fibbzp(nz))
    allocate(fiffx(nx), fifcx(nx), fifbx(nx), fisfx(nx), fiscx(nx), fisbx(nx),fifsx(nx),fifwx(nx),fissx(nx),fiswx(nx))
    allocate(fiffxp(nx),fifsxp(nx),fifwxp(nx),fisfxp(nx),fissxp(nx),fiswxp(nx))
    allocate(fiffy(ny), fifcy(ny), fifby(ny), fisfy(ny), fiscy(ny), fisby(ny),fifsy(ny),fifwy(ny),fissy(ny),fiswy(ny))
    allocate(fiffyp(ny),fifsyp(ny),fifwyp(ny),fisfyp(ny),fissyp(ny),fiswyp(ny))
    allocate(fiffz(nz), fifcz(nz), fifbz(nz), fisfz(nz), fiscz(nz), fisbz(nz),fifsz(nz),fifwz(nz),fissz(nz),fiswz(nz))
    allocate(fiffzp(nz),fifszp(nz),fifwzp(nz),fisfzp(nz),fisszp(nz),fiswzp(nz))

    
    ! Allocate the derivatives coefficients
    allocate(ffx(nx),fcx(nx),fbx(nx),sfx(nx),scx(nx),sbx(nx),fsx(nx),fwx(nx),ssx(nx),swx(nx))
    allocate(ffxp(nx),fsxp(nx),fwxp(nx),sfxp(nx),ssxp(nx),swxp(nx))
    allocate(ffy(ny),fcy(ny),fby(ny),sfy(ny),scy(ny),sby(ny),fsy(ny),fwy(ny),ssy(ny),swy(ny))
    allocate(ffyp(ny),fsyp(ny),fwyp(ny),sfyp(ny),ssyp(ny),swyp(ny))
    allocate(ffz(nz),fcz(nz),fbz(nz),sfz(nz),scz(nz),sbz(nz),fsz(nz),fwz(nz),ssz(nz),swz(nz))
    allocate(ffzp(nz),fszp(nz),fwzp(nz),sfzp(nz),sszp(nz),swzp(nz))

    ! Allocate Pressure derivatives
    ! X-direction
    allocate(cfx6(nxm),ccx6(nxm),cbx6(nxm),cfxp6(nxm),ciwxp6(nxm),csxp6(nxm),&
     cwxp6(nxm),csx6(nxm),cwx6(nxm),cifx6(nxm),cicx6(nxm),cisx6(nxm))
    allocate(cibx6(nxm),cifxp6(nxm),cisxp6(nxm),ciwx6(nxm))
    allocate(cfi6(nx),cci6(nx),cbi6(nx),cfip6(nx),csip6(nx),cwip6(nx),csi6(nx),cwi6(nx),&
        cifi6(nx),cici6(nx),cibi6(nx),cifip6(nx)) 
    allocate(cisip6(nx),ciwip6(nx),cisi6(nx),ciwi6(nx))
    ! Y-direction
    allocate(cfy6(nym),ccy6(nym),cby6(nym),cfyp6(nym),csyp6(nym),cwyp6(nym),csy6(nym))
    allocate(cwy6(nym),cify6(nym),cicy6(nym),ciby6(nym),cifyp6(nym),cisyp6(nym),&
     ciwyp6(nym),cisy6(nym),ciwy6(nym)) 
    allocate(cfi6y(ny),cci6y(ny),cbi6y(ny),cfip6y(ny),csip6y(ny),cwip6y(ny),&
     csi6y(ny),cwi6y(ny),cifi6y(ny),cici6y(ny)) 
    allocate(cibi6y(ny),cifip6y(ny),cisip6y(ny),ciwip6y(ny),cisi6y(ny),ciwi6y(ny)) 
    ! Z-direction
    allocate(cfz6(nzm),ccz6(nzm),cbz6(nzm),cfzp6(nzm),cszp6(nzm),cwzp6(nzm),csz6(nzm))
    allocate(cwz6(nzm),cifz6(nzm),cicz6(nzm),cibz6(nzm),cifzp6(nzm),ciszp6(nzm),ciwzp6(nzm),cisz6(nzm),ciwz6(nzm))
    allocate(cfi6z(nz),cci6z(nz),cbi6z(nz),cfip6z(nz),csip6z(nz),cwip6z(nz),& 
        csi6z(nz),cwi6z(nz),cifi6z(nz),cici6z(nz))  
    allocate(cibi6z(nz),cifip6z(nz),cisip6z(nz),ciwip6z(nz),cisi6z(nz),ciwi6z(nz)) 

    ! Allocate waves
    allocate(zkz(nz/2+1),zk2(nz/2+1),ezs(nz/2+1))
    allocate(yky(ny),yk2(ny),eys(ny))
    allocate(xkx(nx),xk2(nx),exs(nx))

    ! Allocate the stretch-mesh parameters
    allocate(ppy(ny),pp2y(ny),pp4y(ny),ppyi(ny),pp2yi(ny),pp4yi(ny),yp(ny),ypi(ny),del(ny),yeta(ny),yetai(ny))
    
    ! Do the stretching
    if (istret.eq.0) then
       do j=1,ny
          yp(j)=(j-1.)*dy
          ypi(j)=(j-0.5)*dy
       enddo
    else
       call stretching()
    endif
    
    ! Define non-dimensional parameters    
    xnu=1./re
   
    ! Define the coefficients for time-marching
    adt(:)=0. ; bdt(:)=0. ; cdt(:)=0. ; gdt(:)=0.
    if (nscheme==1) then!AB2
       iadvance_time=1 
       adt(1)=1.5*dt
       bdt(1)=-0.5*dt
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)
    endif
    if (nscheme==2) then !RK3
       iadvance_time=3 
       adt(1)=(8./15.)*dt
       bdt(1)=0.
       gdt(1)=adt(1)
       adt(2)=(5./12.)*dt
       bdt(2)=(-17./60.)*dt
       gdt(2)=adt(2)+bdt(2)
       adt(3)=(3./4.)*dt
       bdt(3)=(-5./12.)*dt
       gdt(3)=adt(3)+bdt(3)
    endif
    if (nscheme==3) then !RK4 Carpenter and Kennedy  
       iadvance_time=5 
       adt(1)=0.
       adt(2)=-0.4178904745
       adt(3)=-1.192151694643
       adt(4)=-1.697784692471
       adt(5)=-1.514183444257
       bdt(1)=0.1496590219993
       bdt(2)=0.3792103129999
       bdt(3)=0.8229550293869
       bdt(4)=0.6994504559488
       bdt(5)=0.1530572479681
       gdt(1)=0.1496590219993*dt
       gdt(2)=0.220741935365*dt
       gdt(3)=0.25185480577*dt
       gdt(4)=0.33602636754*dt
       gdt(5)=0.041717869325*dt
    endif
    
    if (nscheme==4) then!AB3
       iadvance_time=1
       adt(1)= (23./12.)*dt
       bdt(1)=-(16./12.)*dt
       cdt(1)= ( 5./12.)*dt
       gdt(1)=adt(1)+bdt(1)+cdt(1)
       gdt(3)=gdt(1)
    endif


    end subroutine init_module_parameters

end module param 

!module param
!
!use decomp_2d, only : mytype
!
!  integer :: nclx,ncly,nclz
!  integer :: ifft, ivirt,istret,iforc_entree,iturb, ialm, Nturbines, NActuatorlines
!  integer :: itype, iskew, iin, nscheme, ifirst, ilast, iles, jLES, jADV
!  integer :: isave,ilit,srestart,idebmod, imodulo, idemarre, icommence, irecord
!  integer :: iscalar,ilag,npif,izap
!  integer :: iprobe, Nprobes, Nsampling
!  integer :: iabl, ibmshape, SmagWallDamp, nSmag
!  character :: Probelistfile*80
!  integer :: nxboite, istat,iread,iadvance_time, ibuoyancy, icoriolis
!  real(mytype) :: xlx,yly,zlz,dx,dy,dz,dx2,dy2,dz2
!  real(mytype) :: dt,xnu,noise,noise1,pi,twopi,u1,u2,re,sc,Pr,TempRef,CoriolisFreq
!  real(mytype) :: t,xxk1,xxk2, spinup_time
!  real(mytype) :: smagcst, walecst,dys, FSGS, rxxnu
!  real(mytype) :: eps_factor ! Smoothing factor 
!  real(mytype) :: TurbRadius,z_zero,k_roughness,PsiM,ustar,u_shear,IPressureGradient,epsilon_pert
!  real(mytype),dimension(3) :: Ug
!  integer :: itr,itime
!  character :: dirname*80
!  character :: filesauve*80, filenoise*80, &
!       nchamp*80,filepath*80, fileturb*80, filevisu*80 
!  character, dimension(100) :: turbinesPath*80, ActuatorlinesPath*80 ! Assign a maximum number of 100 turbines and alms
!   real(mytype), dimension(5) :: adt,bdt,cdt,gdt
!end module param

module IBM

use decomp_2d, only : mytype

  real(mytype) :: cex,cey,cez,ra
end module IBM

module complex_geometry

use decomp_2d, only : mytype
use param, only : nx,ny,nz

  integer     ,parameter                  :: nobjmax=4
  integer     ,save, allocatable, dimension          (:,:) :: nobjx
  integer     ,save, allocatable, dimension          (:,:) :: nobjy
  integer     ,save, allocatable, dimension          (:,:) :: nobjz
  integer     ,save, allocatable, dimension(:,:,:) :: nxipif,nxfpif
  integer     ,save, allocatable, dimension(:,:,:) :: nyipif,nyfpif
  integer     ,save, allocatable, dimension(:,:,:) :: nzipif,nzfpif
  real(mytype),save, allocatable, dimension(:,:,:) :: xi,xf
  real(mytype),save, allocatable, dimension(:,:,:) :: yi,yf
  real(mytype),save, allocatable, dimension(:,:,:) :: zi,zf
    
  contains
      
        subroutine init_complex_geometry(nx1,ny1,nz1)
            implicit none
            integer nx1, ny1, nz1        
            
            allocate(nobjx(ny1,nz1),nobjy(nx1,nz1),nobjz(nx1,ny1)) 
            allocate(nxipif(0:nobjmax,ny1,nz1),nxfpif(0:nobjmax,ny1,nz1))
            allocate(nyipif(0:nobjmax,nx1,nz1),nyfpif(0:nobjmax,nx1,nz1))
            allocate(nzipif(0:nobjmax,nx1,ny1),nzfpif(0:nobjmax,nx1,ny1))
            allocate(xi(nobjmax,ny1,nz1),xf(nobjmax,ny1,nz1))
            allocate(yi(nobjmax,nx1,nz1),yf(nobjmax,nx1,nz1))
            allocate(zi(nobjmax,nx1,ny1),zf(nobjmax,nx1,ny1)) 

        end subroutine init_complex_geometry

end module complex_geometry


module derivX

use decomp_2d, only : mytype

  real(mytype) :: alcaix6,acix6,bcix6
  real(mytype) :: ailcaix6,aicix6,bicix6,cicix6
  real(mytype) :: alfa1x,af1x,bf1x,cf1x,df1x,alfa2x,af2x,alfanx,afnx,bfnx
  real(mytype) :: cfnx,dfnx,alfamx,afmx,alfaix,afix,bfix,alsa1x,as1x,bs1x
  real(mytype) :: cs1x,ds1x,alsa2x,as2x,alsanx,asnx,bsnx,csnx,dsnx,alsamx
  real(mytype) :: asmx,alsa3x,as3x,bs3x,alsatx,astx,bstx 
  real(mytype) :: alsaixt,asixt,bsixt,csixt,dsixt
  !O6SVV
  real(mytype) :: alsa4x,as4x,bs4x,cs4x
  real(mytype) :: alsattx,asttx,bsttx,csttx
  real(mytype) :: alsaix,asix,bsix,csix,dsix
  !
end module derivX

module derivY

use decomp_2d, only : mytype

  real(mytype) :: alcaiy6,aciy6,bciy6
  real(mytype) :: ailcaiy6,aiciy6,biciy6,ciciy6 
  real(mytype) :: alfa1y,af1y,bf1y,cf1y,df1y,alfa2y,af2y,alfany,afny,bfny
  real(mytype) :: cfny,dfny,alfamy,afmy,alfajy,afjy,bfjy,alsa1y,as1y,bs1y
  real(mytype) :: cs1y,ds1y,alsa2y,as2y,alsany,asny,bsny,csny,dsny,alsamy
  real(mytype) :: asmy,alsa3y,as3y,bs3y,alsaty,asty,bsty 
  real(mytype) :: alsajyt,asjyt,bsjyt,csjyt,dsjyt
  !O6SVV
  real(mytype) :: alsa4y,as4y,bs4y,cs4y
  real(mytype) :: alsatty,astty,bstty,cstty
  real(mytype) :: alsajy,asjy,bsjy,csjy,dsjy
  !
end module derivY

module derivZ

use decomp_2d, only : mytype

  real(mytype) :: alcaiz6,aciz6,bciz6
  real(mytype) :: ailcaiz6,aiciz6,biciz6,ciciz6
  real(mytype) :: alfa1z,af1z,bf1z,cf1z,df1z,alfa2z,af2z,alfanz,afnz,bfnz
  real(mytype) :: cfnz,dfnz,alfamz,afmz,alfakz,afkz,bfkz,alsa1z,as1z,bs1z
  real(mytype) :: cs1z,ds1z,alsa2z,as2z,alsanz,asnz,bsnz,csnz,dsnz,alsamz
  real(mytype) :: asmz,alsa3z,as3z,bs3z,alsatz,astz,bstz
  real(mytype) :: alsakzt,askzt,bskzt,cskzt,dskzt
  !O6SVV
  real(mytype) :: alsa4z,as4z,bs4z,cs4z
  real(mytype) :: alsattz,asttz,bsttz,csttz
  real(mytype) :: alsakz,askz,bskz,cskz,dskz
  !
end module derivZ

! Describes the parameters for the discrete filters in X-Pencil
module parfiX
use decomp_2d, only : mytype
  real(mytype) :: fia1x, fib1x, fic1x, fid1x, fie1x, fif1x, fig1x ! Coefficients for filter at boundary point 1  
  real(mytype) :: fia2x, fib2x, fic2x, fid2x, fie2x, fif2x, fig2x ! Coefficients for filter at boundary point 2
  real(mytype) :: fia3x, fib3x, fic3x, fid3x, fie3x, fif3x, fig3x ! Coefficients for filter at boundary point 3
  real(mytype) :: fialx, fiaix, fibix, ficix, fidix               ! Coefficient for filter at interior points 
  real(mytype) :: fianx, fibnx, ficnx, fidnx, fienx, fifnx, fignx ! Coefficient for filter at boundary point n 
  real(mytype) :: fiamx, fibmx, ficmx, fidmx, fiemx, fifmx, figmx ! Coefficient for filter at boundary point m=n-1 
  real(mytype) :: fiapx, fibpx, ficpx, fidpx, fiepx, fifpx, figpx ! Coefficient for filter at boundary point p=n-2
end module parfiX
!
module parfiY

use decomp_2d, only : mytype
  real(mytype) :: fia1y, fib1y, fic1y, fid1y, fie1y, fif1y, fig1y ! Coefficients for filter at boundary point 1  
  real(mytype) :: fia2y, fib2y, fic2y, fid2y, fie2y, fif2y, fig2y ! Coefficients for filter at boundary point 2
  real(mytype) :: fia3y, fib3y, fic3y, fid3y, fie3y, fif3y, fig3y ! Coefficients for filter at boundary point 3
  real(mytype) :: fialy, fiajy, fibjy, ficjy, fidjy               ! Coefficient for filter at interior points 
  real(mytype) :: fiany, fibny, ficny, fidny, fieny, fifny, figny ! Coefficient for filter at boundary point n 
  real(mytype) :: fiamy, fibmy, ficmy, fidmy, fiemy, fifmy, figmy ! Coefficient for filter at boundary point m=n-1 
  real(mytype) :: fiapy, fibpy, ficpy, fidpy, fiepy, fifpy, figpy ! Coefficient for filter at boundary point p=n-2
end module parfiY

module parfiZ

use decomp_2d, only : mytype
  real(mytype) :: fia1z, fib1z, fic1z, fid1z, fie1z, fif1z, fig1z ! Coefficients for filter at boundary point 1  
  real(mytype) :: fia2z, fib2z, fic2z, fid2z, fie2z, fif2z, fig2z ! Coefficients for filter at boundary point 2
  real(mytype) :: fia3z, fib3z, fic3z, fid3z, fie3z, fif3z, fig3z ! Coefficients for filter at boundary point 3
  real(mytype) :: fialz, fiakz, fibkz, fickz, fidkz               ! Coefficient for filter at interior points 
  real(mytype) :: fianz, fibnz, ficnz, fidnz, fienz, fifnz, fignz ! Coefficient for filter at boundary point n 
  real(mytype) :: fiamz, fibmz, ficmz, fidmz, fiemz, fifmz, figmz ! Coefficient for filter at boundary point m=n-1 
  real(mytype) :: fiapz, fibpz, ficpz, fidpz, fiepz, fifpz, figpz ! Coefficient for filter at boundary point p=n-2
end module parfiZ



