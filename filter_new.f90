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
integer ::  i,j,k
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
! Periodic case
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
! Prepare coefficients to be used in the Thomas Algorithm
call prepare (fifbx,fifcx,fiffx,fifsx,fifwx,nx)
!========================================
! Define filter coefficients for Y-pencil
!========================================
fialy=af                         ! alpha_f
!Interior points
fiajy=(11. + 10.*af)/16.         ! a
fibjy=0.5*(15. + 34.*af)/32.     ! b/2 
ficjy=0.5*(-3. + 6.*af)/16.      ! c/2
fidjy=0.5*(1. - 2.*af)/32.       ! d/2
!Boundary point 1
fia1y=0.5*(63./64.+af/64.)       ! a1
fib1y=0.5*(3./32.+29.*af/32.)    ! b1/2
fic1y=0.5*(-15./64.+ 15.*af/64.) ! c1/2
fid1y=0.5*(5./16.-5.*af/16.)     ! d1/2
fie1y=0.5*(-15./16.+15.*af/64.)  ! e1/2
fif1y=0.5*(3./32.-3.*af/32.)     ! f1/2
fig1y=0.5*(-1./64.+af/64.)       ! g1/2
!Boundary point 2
fia2y=0.5*(1./64.+31.*af/32.)    ! a2
fib2y=0.5*(29./32.+3.*af/16.)    ! b2/2
fic2y=0.5*(15./64.+ 17.*af/32.)  ! c2/2
fid2y=0.5*(-5./16.+5.*af/8.)     ! d2/2
fie2y=0.5*(15./64.-15.*af/32.)   ! e2/2
fif2y=0.5*(-3./32.+3.*af/16.)    ! f2/2
fig2y=0.5*(1./64.-af/32.)        ! g2/2
!Boundary point 3
fia3y=0.5*(-1./64.+af/32.)       ! a3
fib3y=0.5*(3./32.+13.*af/16.)    ! b3/2
fic3y=0.5*(49./64.+15.*af/32.)   ! c3/2
fid3y=0.5*(5./16.+3.*af/8.)      ! d3/2
fie3y=0.5*(-15./64.+15.*af/32.)  ! e3/2
fif3y=0.5*(3./32.-3.*af/16.)     ! f3/2
fig3y=0.5*(-1./64.+af/32.)       ! g3/2
!Boundary point n
fiany=0.5*(63./64.+af/64.)       ! an
fibny=0.5*(3./32.+29.*af/32.)    ! bn/2
ficny=0.5*(-15./64.+ 15.*af/64.) ! cn/2
fidny=0.5*(5./16.-5.*af/16.)     ! dn/2
fieny=0.5*(-15./16.+15.*af/64.)  ! en/2
fifny=0.5*(3./32.-3.*af/32.)     ! fn/2
figny=0.5*(-1./64.+af/64.)       ! gn/2
!Boundary point m=n-1
fiamy=0.5*(1./64.+31.*af/32.)    ! am 
fibmy=0.5*(29./32.+3.*af/16.)    ! bm/2
ficmy=0.5*(15./64.+ 17.*af/32.)  ! cm/2
fidmy=0.5*(-5./16.+5.*af/8.)     ! dm/2
fiemy=0.5*(15./64.-15.*af/32.)   ! em/2
fifmy=0.5*(-3./32.+3.*af/16.)    ! fm/2
figmy=0.5*(1./64.-af/32.)        ! gm/2
!Boundary point p=n-2
fiapy=0.5*(-1./64.+af/32.)       ! ap 
fibpy=0.5*(3./32.+13.*af/16.)    ! bp/2
ficpy=0.5*(49./64.+15.*af/32.)   ! cp/2
fidpy=0.5*(5./16.+3.*af/8.)      ! dp/2
fiepy=0.5*(-15./64.+15.*af/32.)  ! ep/2
fifpy=0.5*(3./32.-3.*af/16.)     ! fp/2
figpy=0.5*(-1./64.+af/32.)       ! gp/2
! Define coefficients
if (ncly.eq.1) then
   fiffy(1)   =af+af
   fiffy(2)   =af
   fiffy(ny-2)=af
   fiffy(ny-1)=af
   fiffy(ny)  =0.
   fifcy(1)   =1.
   fifcy(2)   =1.
   fifcy(ny-2)=1.
   fifcy(ny-1)=1.
   fifcy(ny  )=1.
   fifby(1)   =af 
   fifby(2)   =af
   fifby(ny-2)=af
   fifby(ny-1)=af+af
   fifby(ny  )=0.
   do j=3,ny-3
      fiffy(j)=af
      fifcy(j)=1.
      fifby(j)=af
   enddo
endif
call prepare(fifby,fifcy,fiffy,fifsy,fifwy,ny)
!========================================
! Define filter coefficients for Z-pencil
!========================================
fialz=af                         ! alpha_f
!Interior points
fiakz=(11. + 10.*af)/16.         ! a
fibkz=0.5*(15. + 34.*af)/32.     ! b/2 
fickz=0.5*(-3. + 6.*af)/16.      ! c/2
fidkz=0.5*(1. - 2.*af)/32.       ! d/2
!Boundary point 1
fia1z=0.5*(63./64.+af/64.)       ! a1
fib1z=0.5*(3./32.+29.*af/32.)    ! b1/2
fic1z=0.5*(-15./64.+ 15.*af/64.) ! c1/2
fid1z=0.5*(5./16.-5.*af/16.)     ! d1/2
fie1z=0.5*(-15./16.+15.*af/64.)  ! e1/2
fif1z=0.5*(3./32.-3.*af/32.)     ! f1/2
fig1z=0.5*(-1./64.+af/64.)       ! g1/2
!Boundary point 2
fia2z=0.5*(1./64.+31.*af/32.)    ! a2
fib2z=0.5*(29./32.+3.*af/16.)    ! b2/2
fic2z=0.5*(15./64.+ 17.*af/32.)  ! c2/2
fid2z=0.5*(-5./16.+5.*af/8.)     ! d2/2
fie2z=0.5*(15./64.-15.*af/32.)   ! e2/2
fif2z=0.5*(-3./32.+3.*af/16.)    ! f2/2
fig2z=0.5*(1./64.-af/32.)        ! g2/2
!Boundary point 3
fia3z=0.5*(-1./64.+af/32.)       ! a3
fib3z=0.5*(3./32.+13.*af/16.)    ! b3/2
fic3z=0.5*(49./64.+15.*af/32.)   ! c3/2
fid3z=0.5*(5./16.+3.*af/8.)      ! d3/2
fie3z=0.5*(-15./64.+15.*af/32.)  ! e3/2
fif3z=0.5*(3./32.-3.*af/16.)     ! f3/2
fig3z=0.5*(-1./64.+af/32.)       ! g3/2
!Boundary point n
fianz=0.5*(63./64.+af/64.)       ! an
fibnz=0.5*(3./32.+29.*af/32.)    ! bn/2
ficnz=0.5*(-15./64.+ 15.*af/64.) ! cn/2
fidnz=0.5*(5./16.-5.*af/16.)     ! dn/2
fienz=0.5*(-15./16.+15.*af/64.)  ! en/2
fifnz=0.5*(3./32.-3.*af/32.)     ! fn/2
fignz=0.5*(-1./64.+af/64.)       ! gn/2
!Boundary point m=n-1
fiamz=0.5*(1./64.+31.*af/32.)    ! am 
fibmz=0.5*(29./32.+3.*af/16.)    ! bm/2
ficmz=0.5*(15./64.+ 17.*af/32.)  ! cm/2
fidmz=0.5*(-5./16.+5.*af/8.)     ! dm/2
fiemz=0.5*(15./64.-15.*af/32.)   ! em/2
fifmz=0.5*(-3./32.+3.*af/16.)    ! fm/2
figmz=0.5*(1./64.-af/32.)        ! gm/2
!Boundary point p=n-2
fiapz=0.5*(-1./64.+af/32.)       ! ap 
fibpz=0.5*(3./32.+13.*af/16.)    ! bp/2
ficpz=0.5*(49./64.+15.*af/32.)   ! cp/2
fidpz=0.5*(5./16.+3.*af/8.)      ! dp/2
fiepz=0.5*(-15./64.+15.*af/32.)  ! ep/2
fifpz=0.5*(3./32.-3.*af/16.)     ! fp/2
figpz=0.5*(-1./64.+af/32.)       ! gp/2

if (nclz.eq.0) then
      fiffz(1)   =af
      fiffz(2)   =af
      fiffz(nz-2)=af
      fiffz(nz-1)=af
      fiffz(nz)  =0.
      fifcz(1)   =2.
      fifcz(2)   =1.
      fifcz(nz-2)=1.
      fifcz(nz-1)=1.
      fifcz(nz  )=1.+af*af
      fifbz(1)   =af
      fifbz(2)   =af
      fifbz(nz-2)=af
      fifbz(nz-1)=af
      fifbz(nz  )=0.
      do k=3,nz-3
         fiffz(k)=af
         fifcz(k)=1.
         fifbz(k)=af
      enddo
endif
call prepare (fifbz,fifcz,fiffz,fifsz,fifwz,nz)

return 
end subroutine filter


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
      tx(1,j,k)=fiaix*ux(1,j,k)+fibix*(ux(2,j,k)+ux(nx,j,k))& 
                               +ficix*(ux(3,j,k)+ux(nx-1,j,k))&
                               +fidix*(ux(4,j,k)+ux(nx-2,j,k)) 
      rx(1,j,k)=-1.
      tx(2,j,k)=fiaix*ux(2,j,k)+fibix*(ux(3,j,k)+ux(1,j,k))&
                               +ficix*(ux(4,j,k)+ux(nx,j,k))& 
                               +fidix*(ux(5,j,k)+ux(nx-1,j,k)) 
      rx(2,j,k)=0. 
      tx(3,j,k)=fiaix*ux(3,j,k)+fibix*(ux(4,j,k)+ux(2,j,k))&
          
                               +ficix*(ux(5,j,k)+ux(1,j,k))& 
                               +fidix*(ux(6,j,k)+ux(nx,j,k)) 
      rx(3,j,k)=0. 
      do i=4,nx-3
         tx(i,j,k)=fiaix*ux(i,j,k)+fibix*(ux(i+1,j,k)+ux(i-1,j,k))& 
                                  +ficix*(ux(i+2,j,k)+ux(i-2,j,k))&
                                  +fidix*(ux(i+3,j,k)+ux(i-3,j,k)) 
         rx(i,j,k)=0. 
      enddo
      tx(nx-2,j,k)=fiaix*ux(nx-2,j,k)+fibix*(ux(nx-3,j,k)+ux(nx-1,j,k))&
                                     +ficix*(ux(nx-4,j,k)+ux(nx,j,k))& 
                                     +fidix*(ux(nx-5,j,k)+ux(1,j,k)) 
      rx(nx-2,j,k)=0. 
      tx(nx-1,j,k)=fiaix*ux(nx-1,j,k)+fibix*(ux(nx-2,j,k)+ux(nx,j,k))&
                                     +ficix*(ux(nx-3,j,k)+ux(1,j,k))& 
                                     +fidix*(ux(nx-4,j,k)+ux(2,j,k)) 
      rx(nx-1,j,k)=0. 
      tx(nx,j,k)=fiaix*ux(nx,j,k)+fibix*(ux(nx-1,j,k)+ux(1,j,k))&
                                 +ficix*(ux(nx-2,j,k)+ux(2,j,k))& 
                                 +fidix*(ux(nx-3,j,k)+ux(3,j,k)) 
      rx(nx,j,k)=fialx           
      do i=2, nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fifsx(i) 
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*fifsx(i) 
      enddo
      tx(nx,j,k)=tx(nx,j,k)*fifwx(nx) 
      rx(nx,j,k)=rx(nx,j,k)*fifwx(nx) 
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-fiffx(i)*tx(i+1,j,k))*fifwx(i) 
         rx(i,j,k)=(rx(i,j,k)-fiffx(i)*rx(i+1,j,k))*fifwx(i) 
      enddo
        fisx(j,k)=(tx(1,j,k)-fialx*tx(nx,j,k))&
           /(1.+rx(1,j,k)-fialx*rx(nx,j,k)) 
      do i=1,nx 
         tx(i,j,k)=tx(i,j,k)-fisx(j,k)*rx(i,j,k) 
      enddo
   enddo
   enddo
endif

return  
end subroutine filx


!********************************************************************
!
subroutine fily(ty,uy,ry,fisy,fiffy,fifsy,fifwy,ppy,nx,ny,nz,npaire) 
!
!********************************************************************
  
USE param, only: ncly, istret
USE parfiY

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(mytype), dimension(nx,ny,nz) :: ty,uy 
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz)  :: fisy
real(mytype), dimension(ny) :: fiffy,fifsy,fifwy,ppy


if (ncly==1) then 
   if (npaire==1) then 
      do k=1,nz 
      do i=1,nx 
         ty(i,1,k)=0. 
         ty(i,2,k)=fiajy*uy(i,2,k)+fibjy*(uy(i,3,k)+uy(i,1,k))& 
                                  +ficjy*(uy(i,4,k)+uy(i,2,k))&
                                  +fidjy*(uy(i,5,k)+uy(i,3,k)) 
         ty(i,3,k)=fiajy*uy(i,3,k)+fibjy*(uy(i,4,k)+uy(i,2,k))& 
                                  +ficjy*(uy(i,5,k)+uy(i,3,k))&
                                  +fidjy*(uy(i,6,k)+uy(i,4,k)) 
      enddo
      enddo 
      do k=1,nz 
      do j=4,ny-3 
      do i=1,nx 
         ty(i,j,k)=fiajy*uy(i,j,k)+fibjy*(uy(i,j+1,k)+uy(i,j-1,k))& 
                                  +ficjy*(uy(i,j+2,k)+uy(i,j-2,k))&
                                  +fidjy*(uy(i,j+3,k)+uy(i,j-3,k)) 
      enddo
      enddo 
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,ny,k)=0. 
         ty(i,ny-1,k)=fiajy*uy(i,ny-1,k)+fibjy*(uy(i,ny,k)  +uy(i,ny-2,k))& 
                                        +ficjy*(uy(i,ny-1,k)+uy(i,ny-3,k))&
                                        +fidjy*(uy(i,ny-2,k)+uy(i,ny-4,k)) 
         ty(i,ny-2,k)=fiajy*uy(i,ny-2,k)+fibjy*(uy(i,ny-1,k)+uy(i,ny-3,k))& 
                                        +ficjy*(uy(i,ny-2,k)+uy(i,ny-4,k))&
                                        +fidjy*(uy(i,ny-3,k)+uy(i,ny-5,k)) 
      enddo
      enddo 
      do k=1,nz
      do j=2,ny  
      do i=1,nx 
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fifsy(j) 
      enddo
      enddo
      enddo 
      do k=1,nz 
      do i=1,nx 
         ty(i,ny,k)=ty(i,ny,k)*fifwy(ny) 
      enddo 
      enddo 
      do k=1,nz
      do j=ny-1,1,-1  
      do i=1,nx 
      ty(i,j,k)=(ty(i,j,k)-fiffy(j)*ty(i,j+1,k))*fifwy(j) 
      enddo 
      enddo 
      enddo 
   endif
!   if (npaire==0) then 
!      do k=1,nz 
!      do i=1,nx 
!         ty(i,1,k)=afjy*(uy(i,2,k)+uy(i,2,k))&
!              +bfjy*(uy(i,3,k)+uy(i,3,k)) 
!         ty(i,2,k)=afjy*(uy(i,3,k)-uy(i,1,k))&
!              +bfjy*(uy(i,4,k)+uy(i,2,k)) 
!      enddo
!      enddo 
!      do k=1,nz 
!      do j=3,ny-2 
!      do i=1,nx 
!         ty(i,j,k)=afjy*(uy(i,j+1,k)-uy(i,j-1,k))&
!              +bfjy*(uy(i,j+2,k)-uy(i,j-2,k)) 
!      enddo
!      enddo 
!      enddo 
!      do k=1,nz 
!      do i=1,nx 
!         ty(i,ny-1,k)=afjy*(uy(i,ny,k)-uy(i,ny-2,k))&
!              +bfjy*((-uy(i,ny-1,k))-uy(i,ny-3,k)) 
!         ty(i,ny,k)=afjy*((-uy(i,ny-1,k))-uy(i,ny-1,k))&
!              +bfjy*((-uy(i,ny-2,k))-uy(i,ny-2,k)) 
!      enddo
!      enddo 
!      do k=1,nz 
!      do j=2,ny 
!      do i=1,nx 
!         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fsy(j) 
!      enddo 
!      enddo 
!      enddo 
!      do k=1,nz 
!      do i=1,nx 
!         ty(i,ny,k)=ty(i,ny,k)*fwy(ny) 
!      enddo
!      enddo 
!      do k=1,nz
!      do j=ny-1,1,-1  
!      do i=1,nx 
!         ty(i,j,k)=(ty(i,j,k)-ffy(j)*ty(i,j+1,k))*fwy(j) 
!      enddo
!      enddo 
!      enddo 
!   endif
endif

end subroutine fily


subroutine filz(tz,uz,rz,fisz,fiffz,fifsz,fifwz,nx,ny,nz,npaire) 
  
USE param, only: nclx, ncly, nclz  
USE parfiZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: fisz
real(mytype), dimension(nz) :: fiffz,fifsz,fifwz

! This is for periodic boundary conditions
if (nclz==0) then 
   do j=1,ny 
   do i=1,nx 
      tz(i,j,1)=fiakz*uz(i,j,1)+fibkz*(uz(i,j,2)+uz(i,j,nz))& 
                               +fickz*(uz(i,j,3)+uz(i,j,nz-1))&
                               +fidkz*(uz(i,j,4)+uz(i,j,nz-2)) 
      rz(i,j,1)=-1.
      tz(i,j,2)=fiakz*uz(i,j,2)+fibkz*(uz(i,j,3)+uz(i,j,1))&
                               +fickz*(uz(i,j,4)+uz(i,j,nz))& 
                               +fidkz*(uz(i,j,5)+uz(i,j,nz-1)) 
      rz(i,j,2)=0. 
      tz(i,j,3)=fiakz*uz(i,j,3)+fibkz*(uz(i,j,4)+uz(i,j,2))&
                               +fickz*(uz(j,j,5)+uz(i,j,1))& 
                               +fidkz*(uz(k,j,6)+uz(i,j,nz)) 
      rz(i,j,3)=0. 
      do k=4,nz-3
         tz(i,j,k)=fiakz*uz(i,j,k)+fibkz*(uz(i,j,k+1)+uz(i,j,k-1))& 
                                  +fickz*(uz(i,j,k+2)+uz(i,j,k-2))&
                                  +fidkz*(uz(i,j,k+3)+uz(i,j,k-3)) 
         rz(i,j,k)=0. 
      enddo
      tz(i,j,nz-2)=fiakz*uz(i,j,nz-2)+fibkz*(uz(i,j,nz-3)+uz(i,j,nz-1))&
                                     +fickz*(uz(i,j,nz-4)+uz(i,j,nz))& 
                                     +fidkz*(uz(i,j,nz-5)+uz(i,j,1)) 
      rz(i,j,nz-2)=0. 
      tz(i,j,nz-1)=fiakz*uz(i,j,nz-1)+fibkz*(uz(i,j,nz-2)+uz(i,j,nz))&
                                     +fickz*(uz(i,j,nz-3)+uz(i,j,1))& 
                                     +fidkz*(uz(i,j,nz-4)+uz(i,j,2)) 
      rz(i,j,nz-1)=0. 
      tz(i,j,nz)=fiakz*uz(i,j,nz)+fibkz*(uz(i,j,nz-1)+uz(i,j,1))&
                                 +fickz*(uz(i,j,nz-2)+uz(i,j,2))& 
                                 +fidkz*(uz(i,j,nz-3)+uz(i,j,3)) 
      rz(i,j,nz)=fialz           
      do k=2, nz
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fifsz(k) 
         rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*fifsz(k) 
      enddo
      tz(i,j,nz)=tz(i,j,nz)*fifwz(nz) 
      rz(i,j,nz)=rz(i,j,nz)*fifwz(nz) 
      do k=nz-1,1,-1
         tz(i,j,k)=(tz(i,j,k)-fiffz(k)*tz(i,j,k+1))*fifwz(k) 
         rz(i,j,k)=(rz(i,j,k)-fiffz(k)*rz(i,j,k+1))*fifwz(k) 
      enddo
        fisz(i,j)=(tz(i,j,1)-fialz*tz(i,j,nz))&
           /(1.+rz(i,j,1)-fialz*rz(i,j,nz)) 
      do k=1,nz 
         tz(i,j,k)=tz(i,j,k)-fisz(i,j)*rz(i,j,k) 
      enddo
   enddo
   enddo
endif

return  
end subroutine filz

