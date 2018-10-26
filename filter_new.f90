subroutine filter(af)

USE param 
USE parfiX 
USE parfiY 
USE parfiZ 
!=================================================
! Discrete low-pass filter according to 
!=================================================	
implicit none
real(mytype),intent(in) :: af 

! Set the coefficient for the discrete filter following 
! the tridiagonal filtering of Motheau and Abraham, JCP 2016 
! Filter should be -0.5<filax<0.5

! General Case (entire points)
! alpha*fhat(i-1)+fhat(i)+alpha*fhat(i+1)=af(i)+b/2*[f(i+1)+f(i-1)] + ...

! Coefficients are calculated according to the report of Gaitonde & Visbal, 1998,
! "High-order schemes for Navier-Stokes equations: Algorithm and implementation into FDL3DI"

!========================================
! Define filter coefficients for X-pencil
!========================================
fialx=af                         ! alpha_f
!Interior points
fiaix=(11. + 10.*af)/16.         ! a
fibix=0.5*(15. + 34.*af)/32.     ! b/2 
ficix=0.5*(-3. + 6.*af)/16.      ! c/2
fidix=0.5*(1. - 2.*af)/32.       ! d/2
!Boundary point 1
fia1x=0.5*(63./64.+af/64.)       ! a1
fib1x=0.5*(3./32.+29.*af/32.)    ! b1/2
fic1x=0.5*(-15./64.+ 15.*af/64.) ! c1/2
fid1x=0.5*(5./16.-5.*af/16.)     ! d1/2
fie1x=0.5*(-15./16.+15.*af/64.)  ! e1/2
fif1x=0.5*(3./32.-3.*af/32.)     ! f1/2
fig1x=0.5*(-1./64.+af/64.)       ! g1/2
!Boundary point 2
fia2x=0.5*(1./64.+31.*af/32.)    ! a2
fib2x=0.5*(29./32.+3.*af/16.)    ! b2/2
fic2x=0.5*(15./64.+ 17.*af/32.)  ! c2/2
fid2x=0.5*(-5./16.+5.*af/8.)     ! d2/2
fie2x=0.5*(15./64.-15.*af/32.)   ! e2/2
fif2x=0.5*(-3./32.+3.*af/16.)    ! f2/2
fig2x=0.5*(1./64.-af/32.)        ! g2/2
!Boundary point 3
fia3x=0.5*(-1./64.+af/32.)       ! a3
fib3x=0.5*(3./32.+13.*af/16.)    ! b3/2
fic3x=0.5*(49./64.+15.*af/32.)   ! c3/2
fid3x=0.5*(5./16.+3.*af/8.)      ! d3/2
fie3x=0.5*(-15./64.+15.*af/32.)  ! e3/2
fif3x=0.5*(3./32.-3.*af/16.)     ! f3/2
fig3x=0.5*(-1./64.+af/32.)       ! g3/2
!Boundary point n
fianx=0.5*(63./64.+af/64.)       ! an
fibnx=0.5*(3./32.+29.*af/32.)    ! bn/2
ficnx=0.5*(-15./64.+ 15.*af/64.) ! cn/2
fidnx=0.5*(5./16.-5.*af/16.)     ! dn/2
fienx=0.5*(-15./16.+15.*af/64.)  ! en/2
fifnx=0.5*(3./32.-3.*af/32.)     ! fn/2
fignx=0.5*(-1./64.+af/64.)       ! gn/2
!Boundary point m=n-1
fiamx=0.5*(1./64.+31.*af/32.)    ! am 
fibmx=0.5*(29./32.+3.*af/16.)    ! bm/2
ficmx=0.5*(15./64.+ 17.*af/32.)  ! cm/2
fidmx=0.5*(-5./16.+5.*af/8.)     ! dm/2
fiemx=0.5*(15./64.-15.*af/32.)   ! em/2
fifmx=0.5*(-3./32.+3.*af/16.)    ! fm/2
figmx=0.5*(1./64.-af/32.)        ! gm/2
!Boundary point p=n-2
fiapx=0.5*(-1./64.+af/32.)       ! ap 
fibpx=0.5*(3./32.+13.*af/16.)    ! bp/2
ficpx=0.5*(49./64.+15.*af/32.)   ! cp/2
fidpx=0.5*(5./16.+3.*af/8.)      ! dp/2
fiepx=0.5*(-15./64.+15.*af/32.)  ! ep/2
fifpx=0.5*(3./32.-3.*af/16.)     ! fp/2
figpx=0.5*(-1./64.+af/32.)       ! gp/2

! Set the coefficients for the matrix A
if (nclx.eq.0) then
   fiffx(1)   =af
   fiffx(2)   =af
   fiffx(nx-2)=af
   fiffx(nx-1)=af
   fiffx(nx)  =0.
   fifcx(1)   =2.
   fifcx(2)   =1.
   fifcx(nx-2)=1.
   fifcx(nx-1)=1.
   fifcx(nx  )=1.+af*af
   fifbx(1)   =af
   fifbx(2)   =af
   fifbx(nx-2)=af
   fifbx(nx-1)=af
   fifbx(nx  )=0.
   do i=3,nx-3
      fiffx(i)=af
      fifcx(i)=1.
      fifbx(i)=af
   enddo   
endif

call prepare (fifbx,fifcx,fiffx,fifsx,fifwx,nx)


!======================================
! Define for Y-pencil
!======================================
 
return 
end subroutine filter

!*********************************************************************
!subroutine filx(tx,ux,rx,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,&
!     fiz2x,nx,ny,nz,npaire) 
!*********************************************************************
subroutine filx(tx,ux,rx,fisx,fiffx,fifsx,fifwx,nx,ny,nz,npaire) 
  
USE param, only: nclx, ncly, nclz  
USE parfiX 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx 
real(mytype), dimension(ny,nz) :: fisx
real(mytype), dimension(nx) :: fiffx,fifsx,fifwx

! This is for periodic boundary conditions
if (nclx==0) then
   do k=1,nz
   do j=1,ny
   ! Set coefficients


   enddo
   enddo
endif

return  
end subroutine filx
