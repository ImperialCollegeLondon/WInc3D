!==========================================================
! Subroutines for computing the SFS from the ABL
!==========================================================
subroutine abl_boundary_stresses(ux,uy,uz,taux,tauz,delta)

use param
use variables
use var

implicit none
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,phi
integer :: i,j,k,code
real(mytype) :: ut,ut1,utt,ut11, abl_vel, ABLtaux, ABLtauz, delta
real(mytype) :: ux_HAve_local, uz_HAve_local,S_HAve_local
real(mytype) :: ux_HAve, uz_HAve,S_HAve
real(mytype),dimension(xsize(1),xsize(3)) :: taux,tauz 

! Determine the shear stress using Moeng's formulation
!*****************************************************************************************
    ux_HAve_local=0.
    uz_HAve_local=0.
    S_HAve_local=0.
    do k=1,xsize(3)
    do i=1,xsize(1)
        ux_HAve_local=ux_HAve_local+0.5*(ux(i,1,k)+ux(i,2,k))
        uz_HAve_local=uz_HAve_local+0.5*(uz(i,1,k)+uz(i,2,k))
        S_HAve_local=S_HAve_local+sqrt((0.5*(ux(i,1,k)+ux(i,2,k)))**2.+ (0.5*(uz(i,1,k)+uz(i,2,k)))**2.) 
    enddo
    enddo
    
    ux_HAve_local=ux_HAve_local/xsize(3)/xsize(1)
    uz_HAve_local=uz_HAve_local/xsize(3)/xsize(1)
     S_HAve_local= S_HAve_local/xsize(3)/xsize(1)
    
    !call MPI_ALLREDUCE(ux_HAve_local,ux_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    !call MPI_ALLREDUCE(uz_HAve_local,uz_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    !call MPI_ALLREDUCE(S_HAve_local,S_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    
    ux_HAve=ux_HAve_local!/p_col
    uz_HAve=uz_HAve_local!/p_col
     S_HAve= S_HAve_local!/p_col
    
    if (istret.ne.0) delta=(yp(2)-yp(1))/2.0
    if (istret.eq.0) delta=dy/2.0   
    
    ! Compute the friction velocity u_shear
    u_shear=k_roughness*sqrt(ux_HAve**2+uz_HAve**2)/log(delta/z_zero)
    !Compute the shear stresses
    do k=1,xsize(3)
    do i=1,xsize(1)
    taux(i,k)=-u_shear**2.0*(sqrt((0.5*(ux(i,1,k)+ux(i,2,k)))**2.+(0.5*(ux(i,1,k)+ux(i,2,k)))**2.)*ux_HAve+&
                             S_HAve*(0.5*(ux(i,1,k)+ux(i,2,k))-ux_HAve))/(S_Have*sqrt(ux_HAve**2.+uz_HAve**2.))
    tauz(i,k)=-u_shear**2.0*(sqrt((0.5*(uz(i,1,k)+uz(i,2,k)))**2.+(0.5*(uz(i,1,k)+uz(i,2,k)))**2.)*uz_HAve+&
                             S_HAve*(0.5*(uz(i,1,k)+uz(i,2,k))-uz_HAve))/(S_Have*sqrt(ux_HAve**2.+uz_HAve**2.))
    enddo
    enddo

return

end subroutine abl_boundary_stresses

!subroutine abl_turbulent_flux(ux1,uy1,uz1,tfluxx1,tfluxy1,tfluxz1)
!
!    USE param
!    USE variables
!    USE decomp_2d
!    
!    implicit none
!    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
!    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: taux1,tauy1,tauz1
!    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: di1
!    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
!    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: taux2,tauy2,tauz2
!    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: di2
!    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3
!    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: taux3,tauy3,tauz3
!    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: di3
!
!    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tfx1,tfy1,tfz1
!    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tfx2,tfy2,tfz2
!    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: tfx3,tfy3,tfz3
!    
!    !Final turbulent fluxes
!    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tfluxx1,tfluxy1,tfluxz1
!    real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tfluxx2,tfluxy2,tfluxz2
!    real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: tfluxx3,tfluxy3,tfluxz3
!    
!    integer :: i,j,k,code
!
!    ! First compute the boundary stresses
!    call abl_boundary_stresses(ux1,uy1,uz1,taux1,tauy1,tauz1)
!   
!    ! Compute the derivatives in the x-direction 
!    call derx (tfx1,taux1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!    
!    tfluxx1(:,:,:)=tfx1(:,:,:)
!    
!    !WORK Y-PENCILS
!    call transpose_x_to_y(taux1,taux2)
!    call transpose_x_to_y(tauz1,tauz2)
!    
!    ! Compute the derivatives in the x-direction 
!    call dery(tfx2,taux2,di2,sy,ffyp,fsyp,fwyp,ysize(1),ysize(2),ysize(3),1)
!    call dery(tfz2,tauz2,di2,sy,ffyp,fsyp,fwyp,ysize(1),ysize(2),ysize(3),1)
!    
!    tfluxy2(:,:,:)=tfx2(:,:,:)+tfz2(:,:,:)
!
!    !WORK Z-PENCILS
!    call transpose_y_to_z(tauy2,tauy3)
!    
!    call derz(tfz3,tauz3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
!    
!    tfluxz3(:,:,:)=tfz3(:,:,:)
!    
!    call transpose_z_to_y(tfluxz3,tfluxz2)
!    call transpose_y_to_x(tfluxz2,tfluxz1)
!    call transpose_y_to_x(tfluxy2,tfluxy1)
!
!    
!    if(xstart(2)==1) then
!	write(*,*) tfluxx1(:,1,:), tfluxy1(:,1,:), tfluxz1(:,1,:)
!	tfluxx1(:,1,:)=0.01
!	tfluxy1(:,1,:)=0
!	tfluxz1(:,1,:)=0
!    endif
!	stop
!    
!    return
!end subroutine abl_turbulent_flux
