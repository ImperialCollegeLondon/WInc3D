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

module variables

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

!module filter
real(mytype), save, allocatable, dimension(:) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
real(mytype), save, allocatable, dimension(:,:) ::filax,filaxp
real(mytype), save, allocatable, dimension(:) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
real(mytype), save, allocatable, dimension(:) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y
real(mytype), save, allocatable, dimension(:,:) ::filay,filayp
real(mytype), save, allocatable, dimension(:) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
real(mytype), save, allocatable, dimension(:) :: fifz,ficz,fibz,fiffz,fibbz,fiz1z,fiz2z
real(mytype), save, allocatable, dimension(:,:) ::filaz,filazp
real(mytype), save, allocatable, dimension(:) :: fifzp,ficzp,fibzp,fiffzp,fibbzp
integer, dimension(200) :: idata
real(mytype), save, allocatable, dimension(:) :: Cs


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
    
    subroutine init_module_parameters(nx1,ny1,nz1,nxm1,nym1,nzm1,p_row1,p_col1)

    implicit none
    ! Allocate variables
    integer :: nx1,ny1,nz1,nxm1,nym1,nzm1,p_row1,p_col1
     
    ! Allocate filter and interpolation parameters 
    ! X-direction
    allocate(fifx(nx1),ficx(nx1),fibx(nx1),fiffx(nx1),fibbx(nx1),fiz1x(nx1),fiz2x(nx1))
    allocate(filax(nx1,2),filaxp(nx1,2))
    allocate(fifxp(nx1),ficxp(nx1),fibxp(nx1),fiffxp(nx1),fibbxp(nx1))
    ! Y-direction 
    allocate(fify(ny1),ficy(ny1),fiby(ny1),fiffy(ny1),fibby(ny1),fiz1y(ny1),fiz2y(ny1))
    allocate(filay(ny1,2),filayp(ny1,2))
    allocate(fifyp(ny1),ficyp(ny1),fibyp(ny1),fiffyp(ny1),fibbyp(ny1))
    ! Z-direction  
    allocate(fifz(nz1),ficz(nz1),fibz(nz1),fiffz(nz1),fibbz(nz1),fiz1z(nz1),fiz2z(nz1))
    allocate(filaz(nz1,2),filazp(nz1,2))
    allocate(fifzp(nz1),ficzp(nz1),fibzp(nz1),fiffzp(nz1),fibbzp(nz1))
    
    allocate(Cs(p_col1*p_row1))

    ! Allocate the derivatives coefficients
    allocate(ffx(nx1),fcx(nx1),fbx(nx1),sfx(nx1),scx(nx1),sbx(nx1),fsx(nx1),fwx(nx1),ssx(nx1),swx(nx1))
    allocate(ffxp(nx1),fsxp(nx1),fwxp(nx1),sfxp(nx1),ssxp(nx1),swxp(nx1))
    allocate(ffy(ny1),fcy(ny1),fby(ny1),sfy(ny1),scy(ny1),sby(ny1),fsy(ny1),fwy(ny1),ssy(ny1),swy(ny1))
    allocate(ffyp(ny1),fsyp(ny1),fwyp(ny1),sfyp(ny1),ssyp(ny1),swyp(ny1))
    allocate(ffz(nz1),fcz(nz1),fbz(nz1),sfz(nz1),scz(nz1),sbz(nz1),fsz(nz1),fwz(nz1),ssz(nz1),swz(nz1))
    allocate(ffzp(nz1),fszp(nz1),fwzp(nz1),sfzp(nz1),sszp(nz1),swzp(nz1))

    ! Allocate Pressure derivatives
    ! X-direction
    allocate(cfx6(nxm1),ccx6(nxm1),cbx6(nxm1),cfxp6(nxm1),ciwxp6(nxm1),csxp6(nxm1),&
     cwxp6(nxm1),csx6(nxm1),cwx6(nxm1),cifx6(nxm1),cicx6(nxm1),cisx6(nxm1))
    allocate(cibx6(nxm1),cifxp6(nxm1),cisxp6(nxm1),ciwx6(nxm1))
    allocate(cfi6(nx1),cci6(nx1),cbi6(nx1),cfip6(nx1),csip6(nx1),cwip6(nx1),csi6(nx1),cwi6(nx1),&
        cifi6(nx1),cici6(nx1),cibi6(nx1),cifip6(nx1)) 
    allocate(cisip6(nx1),ciwip6(nx1),cisi6(nx1),ciwi6(nx1))
    ! Y-direction
    allocate(cfy6(nym1),ccy6(nym1),cby6(nym1),cfyp6(nym1),csyp6(nym1),cwyp6(nym1),csy6(nym1))
    allocate(cwy6(nym1),cify6(nym1),cicy6(nym1),ciby6(nym1),cifyp6(nym1),cisyp6(nym1),&
     ciwyp6(nym1),cisy6(nym1),ciwy6(nym1)) 
    allocate(cfi6y(ny1),cci6y(ny1),cbi6y(ny1),cfip6y(ny1),csip6y(ny1),cwip6y(ny1),&
     csi6y(ny1),cwi6y(ny1),cifi6y(ny1),cici6y(ny1)) 
    allocate(cibi6y(ny1),cifip6y(ny1),cisip6y(ny1),ciwip6y(ny1),cisi6y(ny1),ciwi6y(ny1)) 
    ! Z-direction
    allocate(cfz6(nzm1),ccz6(nzm1),cbz6(nzm1),cfzp6(nzm1),cszp6(nzm1),cwzp6(nzm1),csz6(nzm1))
    allocate(cwz6(nzm1),cifz6(nzm1),cicz6(nzm1),cibz6(nzm1),cifzp6(nzm1),ciszp6(nzm1),ciwzp6(nzm1),cisz6(nzm1),ciwz6(nzm1))
    allocate(cfi6z(nz1),cci6z(nz1),cbi6z(nz1),cfip6z(nz1),csip6z(nz1),cwip6z(nz1),& 
        csi6z(nz1),cwi6z(nz1),cifi6z(nz1),cici6z(nz1))  
    allocate(cibi6z(nz1),cifip6z(nz1),cisip6z(nz1),ciwip6z(nz1),cisi6z(nz1),ciwi6z(nz1)) 

    ! Allocate waves
    allocate(zkz(nz1/2+1),zk2(nz1/2+1),ezs(nz1/2+1))
    allocate(yky(ny1),yk2(ny1),eys(ny1))
    allocate(xkx(nx1),xk2(nx1),exs(nx1))

    ! Allocate the stretch-mesh parameters
    allocate(ppy(ny1),pp2y(ny1),pp4y(ny1),ppyi(ny1),pp2yi(ny1),pp4yi(ny1),yp(ny1),ypi(ny1),del(ny1),yeta(ny1),yetai(ny1))

    end subroutine init_module_parameters

end module variables

module param

use decomp_2d, only : mytype

  integer :: nclx,ncly,nclz
  integer :: ifft, ivirt,istret,iforc_entree,iturb, ialm, Nturbines, NActuatorlines
  integer :: itype, iskew, iin, nscheme, ifirst, ilast, iles, jLES, jADV
  integer :: isave,ilit,srestart,idebmod, imodulo, idemarre, icommence, irecord
  integer :: iscalar,ilag,npif,izap
  integer :: iprobe, Nprobes, Nsampling
  integer :: iabl, ibmshape, SmagWallDamp, nSmag
  character :: Probelistfile*80
  integer :: nxboite, istat,iread,iadvance_time, ibuoyancy, icoriolis
  real(mytype) :: xlx,yly,zlz,dx,dy,dz,dx2,dy2,dz2
  real(mytype) :: dt,xnu,noise,noise1,pi,twopi,u1,u2,sc,Pr,TempRef,CoriolisFreq
  real(mytype) :: t,xxk1,xxk2, spinup_time
  real(mytype) :: smagcst, walecst,dys, FSGS, rxxnu
  real(mytype) :: eps_factor ! Smoothing factor 
  real(mytype) :: TurbRadius,z_zero,k_roughness,PsiM,ustar,u_shear,IPressureGradient,epsilon_pert
  real(mytype),dimension(3) :: Ug
  integer :: itr,itime
  character :: dirname*80
  character :: filesauve*80, filenoise*80, &
       nchamp*80,filepath*80, fileturb*80, filevisu*80 
  character, dimension(100) :: turbinesPath*80, ActuatorlinesPath*80 ! Assign a maximum number of 100 turbines and alms
   real(mytype), dimension(5) :: adt,bdt,cdt,gdt
end module param

module IBM

use decomp_2d, only : mytype

  real(mytype) :: cex,cey,cez,ra
end module IBM

module complex_geometry

use decomp_2d, only : mytype
use variables, only : nx,ny,nz

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


module parfiX

use decomp_2d, only : mytype

  real(mytype) :: fia1x, fib1x, fic1x, fid1x, fie1x, fia2x, fib2x, fic2x, fid2x
  real(mytype) :: fie2x, fia3x, fib3x, fic3x, fid3x, fie3x, fianx, fibnx, ficnx, fidnx
  real(mytype) :: fienx, fiamx, fibmx, ficmx, fidmx, fiemx, fiapx, fibpx, ficpx, fidpx
  real(mytype) :: fiepx, fiaix, fibix, ficix, fidix, fialx, fibex, fih1x, fih2x, fih3x,fih4x 
end module parfiX
!
module parfiY

use decomp_2d, only : mytype

  real(mytype) :: fia1y, fib1y, fic1y, fid1y, fie1y, fia2y, fib2y, fic2y, fid2y
  real(mytype) :: fie2y, fia3y, fib3y, fic3y, fid3y, fie3y, fiany, fibny, ficny, fidny
  real(mytype) :: fieny, fiamy, fibmy, ficmy, fidmy, fiemy, fiapy, fibpy, ficpy, fidpy
  real(mytype) :: fiepy, fiaiy, fibiy, ficiy, fidiy, fialy, fibey, fih1y, fih2y, fih3y,fih4y 
end module parfiY

module parfiZ

use decomp_2d, only : mytype

  real(mytype) :: fia1z, fib1z, fic1z, fid1z, fie1z, fia2z, fib2z, fic2z, fid2z
  real(mytype) :: fie2z, fia3z, fib3z, fic3z, fid3z, fie3z, fianz, fibnz, ficnz, fidnz
  real(mytype) :: fienz, fiamz, fibmz, ficmz, fidmz, fiemz, fiapz, fibpz, ficpz, fidpz
  real(mytype) :: fiepz, fiaiz, fibiz, ficiz, fidiz, fialz, fibez, fih1z, fih2z, fih3z,fih4z 
end module parfiZ



