!==========================================================
! Subroutines for computing the SFS from the ABL
!==========================================================
subroutine wall_shear_flux(ux,uy,uz,tauwallxy1,tauwallzy1,wallfluxx,wallfluxy,wallfluxz)

USE param
USE variables
USE decomp_2d
USE MPI

implicit none
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: wallfluxx,wallfluxy,wallfluxz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tauwallxy1, tauwallzy1 
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: tauwallxy2, tauwallzy2 
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: tauwallxy3, tauwallzy3
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxy1,gyx1,gyz1,gzy1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gxy2,gzy2,gyz2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: gyz3,di3
integer :: i,j,k,code
real(mytype) :: ut,ut1,utt,ut11, abl_vel, ABLtaux, ABLtauz, delta
real(mytype) :: ux_HAve_local, uz_HAve_local,S_HAve_local
real(mytype) :: ux_HAve, uz_HAve,S_HAve


! First step we need to zero the coeffcients
tauwallxy1=0.; tauwallzy1=0.;
tauwallxy2=0.; tauwallzy2=0.;

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
    
    call MPI_ALLREDUCE(ux_HAve_local,ux_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uz_HAve_local,uz_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(S_HAve_local,S_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    
    ux_HAve=ux_HAve_local/p_col
    uz_HAve=uz_HAve_local/p_col
     S_HAve= S_HAve_local/p_col
   
    if (istret.ne.0) delta=(yp(2)-yp(1))/2.0
    if (istret.eq.0) delta=dy/2.0   
    
    ! Compute the friction velocity u_shear
    u_shear=k_roughness*sqrt(ux_HAve**2.+uz_HAve**2.)/log(delta/z_zero)
    if (nrank==0) write(*,*) "Friction velocity ... ", u_shear 
    !Compute the shear stresses -- only on the wall

    if (xstart(2)==1) then
    do k=1,xsize(3)
    do i=1,xsize(1)
    tauwallxy1(i,1,k)=-u_shear**2.0*(sqrt((0.5*(ux(i,1,k)+ux(i,2,k)))**2.+(0.5*(ux(i,1,k)+ux(i,2,k)))**2.)*ux_HAve+&
                             S_HAve*(0.5*(ux(i,1,k)+ux(i,2,k))-ux_HAve))/(S_Have*sqrt(ux_HAve**2.+uz_HAve**2.))
    tauwallzy1(i,1,k)=-u_shear**2.0*(sqrt((0.5*(uz(i,1,k)+uz(i,2,k)))**2.+(0.5*(uz(i,1,k)+uz(i,2,k)))**2.)*uz_HAve+&
                             S_HAve*(0.5*(uz(i,1,k)+uz(i,2,k))-uz_HAve))/(S_Have*sqrt(ux_HAve**2.+uz_HAve**2.))
    enddo
    enddo
    endif
!*********************************************************************************************************

! Computing the wall fluxes 
! Derivates for x 
call derx (gyx1,tauwallxy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)

! Transpose X --> Y
call transpose_x_to_y(tauwallxy1,tauwallxy2)
call transpose_x_to_y(tauwallzy1,tauwallzy2)

! Differentiate for y
!call dery (gyx2,tauwallxy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
!call dery (gyz2,tauwallzy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (gxy2,tauwallxy2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (gzy2,tauwallzy2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)

! Transpose Y --> Z
call transpose_y_to_z(tauwallzy2,tauwallzy3)

! Differentiate for z
call derz(gyz3,tauwallzy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)

!Transpose Z --> Y
call transpose_z_to_y(gyz3,gyz2)

! Transpose Y --> X
call transpose_y_to_x(gyz2,gyz1)
call transpose_y_to_x(gxy2,gxy1)
call transpose_y_to_x(gzy2,gzy1)

wallfluxx(:,:,:) = -gxy1(:,:,:)
wallfluxy(:,:,:) = -(gyx1(:,:,:)+gyz1(:,:,:))
wallfluxz(:,:,:) = -gzy1(:,:,:)
if (nrank==0) write(*,*)  'Maximum wallflux for x, y and z', maxval(wallfluxx), maxval(wallfluxy), maxval(wallfluxz)
if (nrank==0) write(*,*)  'Minimum wallflux for x, y and z', minval(wallfluxx), minval(wallfluxy), minval(wallfluxz)
return
end subroutine wall_shear_flux
