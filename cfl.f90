!===================================================
! Subroutine for computing the local and global CFL
! number, according to Lele 1992.
!===================================================
subroutine cfl_compute(uxmax,uymax,uzmax)

use param
use variables
use var
    
implicit none

real(mytype),intent(in) :: uxmax,uymax,uzmax
real(mytype) :: cfl_x_adv,cfl_x_diff,cfl_y_adv,cfl_y_diff,cfl_z_adv,cfl_z_diff
real(mytype) :: cfl_conv_lim, cfl_diff_lim
real(mytype) :: sigma_conv(3), sigma_diff(3)
real(mytype) :: visc

! Set the constants (this is true for periodic boundaries)
sigma_conv=[0.0, sqrt(3.0), 2.85]
sigma_diff=[2.0, 2.5, 2.9]

if(jles==0) then
visc=xnu
elseif (jles==1) then
visc=20*rxxnu*xnu
endif

! This is considering 1D peridic boundaries
! Do x-direction
cfl_x_adv=abs(uxmax)*dt/dx
cfl_x_diff=visc*dt/dx**2.0
! Do y-direction
cfl_y_adv=abs(uymax)*dt/dy
cfl_y_diff=visc*dt/dy**2.0
! Do z-direction
cfl_z_adv=abs(uzmax)*dt/dz
cfl_z_diff=visc*dt/dz**2.0

! So far we will focus on uniform grids
if(nrank==0) then
    write(*,*) ' '
    write(*,1002) cfl_x_adv, cfl_x_diff
1002  format('CFL x-direction (Adv and Diff) =',F9.4,',',F9.4)
    write(*,1003) cfl_y_adv, cfl_y_diff
1003  format('CFL y-direction (Adv and Diff) =',F9.4,',',F9.4)
    write(*,1004) cfl_z_adv, cfl_z_diff
1004  format('CFL z-direction (Adv and Diff) =',F9.4,',',F9.4)
    cfl_conv_lim=sigma_conv(nscheme)/sqrt(3.0)
    cfl_diff_lim=sigma_diff(nscheme)/6.0
    write(*,1005) cfl_conv_lim, cfl_diff_lim
    write(*,*) ' '
1005  format('CFL limits (Adv and Diff) : ',F9.4,',',F9.4)
endif

end subroutine
