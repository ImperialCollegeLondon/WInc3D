!==========================================================
! Subroutines for computing the SFS from the ABL
!==========================================================
subroutine wall_shear_stress(ux,uy,uz,nut1,sxy1,syz1,tauwallxy,tauwallzy,wallfluxx,wallfluxy,wallfluxz)

    USE MPI
    USE decomp_2d
    USE param

implicit none
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,nut1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxf,uyf,uzf
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: wallfluxx,wallfluxy,wallfluxz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxy1, syz1 
real(mytype),dimension(xsize(1),xsize(3)) :: tauwallxy, tauwallzy 
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxy1,gyx1,gyz1,gzy1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gxy2,gzy2,gyz2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: gyz3,di3
integer :: i,j,k,code
real(mytype) :: ut,ut1,utt,ut11, abl_vel, ABLtaux, ABLtauz, delta
real(mytype) :: ux_HAve_local, uz_HAve_local,S_HAve_local
real(mytype) :: ux_HAve, uz_HAve,S_HAve



call filter()

! Filter the velocity with twice the grid scale according to Bou-zeid et al 2005
call filx(uxf,ux,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(uzf,uz,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)

! Determine the shear stress using Moeng's formulation
!*****************************************************************************************
    ux_HAve_local=0.
    uz_HAve_local=0.
    S_HAve_local=0.
    
    if (xstart(2)==1) then
    
    do k=1,xsize(3)
    do i=1,xsize(1)
        ux_HAve_local=ux_HAve_local+0.5*(uxf(i,1,k)+uxf(i,2,k))
        uz_HAve_local=uz_HAve_local+0.5*(uzf(i,1,k)+uzf(i,2,k))
    enddo
    enddo
    
    ux_HAve_local=ux_HAve_local/xsize(3)/xsize(1)
    uz_HAve_local=uz_HAve_local/xsize(3)/xsize(1)
   
    else 
    
    ux_HAve_local=0.  
    uz_HAve_local=0.
    S_HAve_local =0.
   
    endif
    
    call MPI_ALLREDUCE(ux_HAve_local,ux_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uz_HAve_local,uz_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    
    ux_HAve=ux_HAve/p_col
    uz_HAve=uz_HAve/p_col
     !S_HAve= S_HAve/p_col
   
    if (istret.ne.0) delta=(yp(2)-yp(1))/2.0
    if (istret.eq.0) delta=dy/2.0   
    
    ! Compute the friction velocity u_shear
    u_shear=k_roughness*sqrt(ux_HAve**2.+uz_HAve**2.)/log(delta/z_zero)
    if (nrank==0) write(*,*) "Horizontally-averaged velocity at y=1/2... ", ux_HAve,0,uz_Have
    if (nrank==0) write(*,*) "Friction velocity ... ", u_shear 
    !Compute the shear stresses -- only on the wall
    !u_shear=ustar
    wallfluxx = 0. 
    wallfluxy = 0.
    wallfluxz = 0.
    
    if (xstart(2)==1) then
    do k=1,xsize(3)
    do i=1,xsize(1)                        
    tauwallxy(i,k)=-u_shear**2.0*0.5*(uxf(i,1,k)+uxf(i,2,k))/sqrt(ux_HAve**2.+uz_HAve**2.)
    tauwallzy(i,k)=-u_shear**2.0*0.5*(uzf(i,1,k)+uzf(i,2,k))/sqrt(ux_HAve**2.+uz_Have**2.)
    
    if(jLES.ge.2) then ! Apply third order one-sided finite difference 
    wallfluxx(i,1,k) = -(-1./2.*(-2.*nut1(i,3,k)*sxy1(i,3,k))+&
        2.*(-2.*nut1(i,2,k)*sxy1(i,2,k))-3./2.*tauwallxy(i,k))/(2.*delta)
    wallfluxy(i,1,k) = -(tauwallxy(i,k)-tauwallxy(i-1,k))/dx-(tauwallzy(i,k)-tauwallzy(i,k-1))/dz
    wallfluxz(i,1,k) = -(-1./2.*(-2.*nut1(i,3,k)*syz1(i,3,k))+&
        2.*(-2.*nut1(i,2,k)*syz1(i,2,k))-3./2.*tauwallzy(i,k))/(2.*delta)
    else
    wallfluxx(i,1,k) = 0.
    wallfluxy(i,1,k) = 0.
    wallfluxz(i,1,k) = 0.
    endif
    
    enddo
    enddo
     
    endif
!*********************************************************************************************************

if (nrank==0) write(*,*)  'Maximum wall shear stress for x and z', maxval(tauwallxy), maxval(tauwallzy)
if (nrank==0) write(*,*)  'Minimum wall shear stress for x and z', minval(tauwallxy), minval(tauwallzy)

return

end subroutine wall_shear_stress
