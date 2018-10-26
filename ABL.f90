!==========================================================
! Subroutines for computing the SFS from the ABL
!==========================================================
subroutine wall_sgs(ux,uy,uz,nut1,sxy1,syz1,tauwallxy,tauwallzy,wallfluxx,wallfluxy,wallfluxz)

    USE MPI
    USE decomp_2d
    USE param

implicit none
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,nut1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxf,uyf,uzf
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: wallfluxx,wallfluxy,wallfluxz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: wallfluxxf,wallfluxyf,wallfluxzf
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxy1, syz1 
real(mytype),dimension(xsize(1),xsize(3)) :: tauwallxy, tauwallzy 
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxy1,gyx1,gyz1,gzy1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gxy2,gzy2,gyz2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: gyz3,di3
integer :: i,j,k,code
real(mytype) :: ut,ut1,utt,ut11, abl_vel, ABLtaux, ABLtauz, delta
real(mytype) :: ux_HAve_local, uz_HAve_local,S_HAve_local
real(mytype) :: ux_HAve, uz_HAve,S_HAve,ux12,uz12,S12
real(mytype) :: sxy_HAve_local, szy_HAve_local, nut_HAve_local, nutsxy_HAve_local, nutszy_HAve_local
real(mytype) :: sxy_HAve, szy_HAve, nut_HAve, nutsxy_HAve, nutszy_HAve
real(mytype) :: nutprimes
real(mytype) :: nuLESBar, scriptR, xi1, nuLES, TS1, TR1
real(mytype) :: CD ! drag coefficient

call filter(0.)
!
! Filter the velocity with twice the grid scale according to Bou-zeid et al 2005
call filx(uxf,ux,di1,fisx,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0) 

!call filx(uzf,uz,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
!fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)

! Determine the shear stress using Moeng's formulation
!*****************************************************************************************
    ux_HAve_local=0.
    uz_HAve_local=0.
    S_HAve_local=0.
    
    if (xstart(2)==1) then 
    !if (nrank==0) print *, 'Max of the filtered velocity', maxval(ux), 'Max of the unfiltered velocity', maxval(ux) 
    do k=1,xsize(3)
    do i=1,xsize(1)
        ux_HAve_local=ux_HAve_local+0.5*(ux(i,1,k)+ux(i,2,k))
        uz_HAve_local=uz_HAve_local+0.5*(uz(i,1,k)+uz(i,2,k))
        S_HAve_local=S_HAve_local+sqrt((0.5*(ux(i,1,k)+ux(i,2,k)))**2.+(0.5*(uz(i,1,k)+uz(i,2,k)))**2.)
    enddo
    enddo
        ux_HAve_local=ux_HAve_local/xsize(3)/xsize(1)
        uz_HAve_local=uz_HAve_local/xsize(3)/xsize(1)
        S_HAve_local=S_HAve_local/xsize(3)/xsize(1) 
    else 
    ux_HAve_local=0.  
    uz_HAve_local=0.
    S_HAve_local =0.
    endif
    
    call MPI_ALLREDUCE(ux_HAve_local,ux_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uz_HAve_local,uz_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(S_HAve_local,S_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    
    ux_HAve=ux_HAve/p_col
    uz_HAve=uz_HAve/p_col
    S_HAve= S_HAve/p_col
   
    if (istret.ne.0) delta=(yp(2)-yp(1))/2.0
    if (istret.eq.0) delta=dy/2.
    
    ! Compute the friction velocity u_shear and the drag coefficient

    u_shear=k_roughness*sqrt(ux_HAve**2.+uz_HAve**2.)/log(delta/z_zero)
    CD=k_roughness**2./log(delta/z_zero)**2.

    if (nrank==0) then 
        print *, ' '
        print *, ' ABL:'
        print *, ' Horizontally-averaged velocity at y=1/2... ', ux_HAve,uz_Have
        print *, ' Friction velocity : ', u_shear 
        print *, ' Drag Coefficient : ', CD 
        print *, ' xi1 : ', xi1
        print *, ' scriptR ', scriptR
    endif
    !Compute the shear stresses -- only on the wall
    !u_shear=ustar
    wallfluxx = 0. 
    wallfluxy = 0.
    wallfluxz = 0.
   
    ! Apply BCs locally
    if (xstart(2)==1) then
    do k=1,xsize(3)
    do i=1,xsize(1)                       
    ux12=0.5*(ux(i,1,k)+ux(i,2,k)) 
    uz12=0.5*(uz(i,1,k)+uz(i,2,k))
    S12=sqrt(ux12**2.+uz12**2.)
    if(iwallmodel==1) then ! MOENG    
    tauwallxy(i,k)=-(k_roughness/log(delta/z_zero))**2.0*ux_HAve*S_HAve
    tauwallzy(i,k)=-(k_roughness/log(delta/z_zero))**2.0*uz_HAve*S_HAve
    else !Parlage model 
    tauwallxy(i,k)=-(k_roughness/log(delta/z_zero))**2.0*ux12*S_HAve
    tauwallzy(i,k)=-(k_roughness/log(delta/z_zero))**2.0*uz12*S_HAve
    endif
    if(jLES.ge.2) then ! Apply third order one-sided finite difference 
    wallfluxx(i,1,k) = -(-1./2.*(-2.*nut1(i,3,k)*sxy1(i,3,k))+&
        2.*(-2.*nut1(i,2,k)*sxy1(i,2,k))-3./2.*tauwallxy(i,k))/(2.*delta)
    wallfluxy(i,1,k) = 0.
    wallfluxz(i,1,k) = -(-1./2.*(-2.*nut1(i,3,k)*syz1(i,3,k))+&
        2.*(-2.*nut1(i,2,k)*syz1(i,2,k))-3./2.*tauwallzy(i,k))/(2.*delta)
    else
    wallfluxx(i,1,k) = -4*CD*ux12*u_shear/delta 
    wallfluxy(i,1,k) = 0.!
    wallfluxz(i,1,k) = -4*CD*uz12*u_shear/delta
    endif
    enddo
    enddo
    endif
!*********************************************************************************************************
! Filter wallfluxes
!call filx(wallfluxxf,wallfluxx,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
!fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
!call filx(wallfluxyf,wallfluxy,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
!fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
!call filx(wallfluxzf,wallfluxz,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
!fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
!wallfluxx(:,:,:) = wallfluxxf(:,:,:)
!wallfluxy(:,:,:) = wallfluxyf(:,:,:)
!wallfluxz(:,:,:) = wallfluxzf(:,:,:)



if (nrank==0) write(*,*)  'Maximum wall shear stress for x and z', maxval(tauwallxy), maxval(tauwallzy)
if (nrank==0) write(*,*)  'Minimum wall shear stress for x and z', minval(tauwallxy), minval(tauwallzy)

if (nrank==0) write(*,*)  'Wall flux x ', maxval(wallfluxx), maxval(wallfluxz)
if (nrank==0) write(*,*)  'Wall flux z ', minval(wallfluxx), minval(wallfluxz)

return

end subroutine wall_sgs

subroutine wall_shear_stress(ux,uy,uz,di1,tauwallxy,tauwallzy)

    USE MPI
    USE decomp_2d
    USE param

implicit none
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,nut1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxf,uyf,uzf
real(mytype),dimension(xsize(1),xsize(3)) :: tauwallxy, tauwallzy 
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: di1
integer :: i,j,k,code
real(mytype) :: delta
real(mytype) :: ux_HAve_local, uz_HAve_local,S_HAve_local
real(mytype) :: ux_HAve, uz_HAve,S_HAve,ux12,uz12,S12
real(mytype) :: sxy_HAve_local, szy_HAve_local, nut_HAve_local, nutsxy_HAve_local, nutszy_HAve_local
real(mytype) :: sxy_HAve, szy_HAve, nut_HAve, nutsxy_HAve, nutszy_HAve
real(mytype) :: nutprimes
real(mytype) :: nuLESBar, scriptR, xi1, nuLES, TS1, TR1
real(mytype) :: CD ! drag coefficient

call filter(-0.49)
call filx(uxf,ux,di1,fisx,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0) 
call filx(uzf,uz,di1,fisx,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0) 

! Filter the velocity with twice the grid scale according to Bou-zeid et al 2005

! Determine the shear stress using Moeng's formulation
!*****************************************************************************************
    ux_HAve_local=0.
    uz_HAve_local=0.
    S_HAve_local=0.

    if (xstart(2)==1) then 
    if (nrank==0) print *, 'Max of the filtered velocity', maxval(uxf), 'Max of the unfiltered velocity', maxval(ux) 
    do k=1,xsize(3)
    do i=1,xsize(1)
        ux_HAve_local=ux_HAve_local+0.5*(uxf(i,1,k)+uxf(i,2,k))
        uz_HAve_local=uz_HAve_local+0.5*(uzf(i,1,k)+uzf(i,2,k))
        S_HAve_local=S_HAve_local+sqrt((0.5*(uxf(i,1,k)+uxf(i,2,k)))**2.+(0.5*(uzf(i,1,k)+uzf(i,2,k)))**2.)
    enddo
    enddo
        ux_HAve_local=ux_HAve_local/xsize(3)/xsize(1)
        uz_HAve_local=uz_HAve_local/xsize(3)/xsize(1)
        S_HAve_local=S_HAve_local/xsize(3)/xsize(1) 
    else 
    ux_HAve_local=0.  
    uz_HAve_local=0.
    S_HAve_local =0.
    endif
    
    call MPI_ALLREDUCE(ux_HAve_local,ux_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uz_HAve_local,uz_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(S_HAve_local,S_HAve,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    
    ux_HAve=ux_HAve/p_col
    uz_HAve=uz_HAve/p_col
    S_HAve= S_HAve/p_col
   
    if (istret.ne.0) delta=(yp(2)-yp(1))/2.0
    if (istret.eq.0) delta=dy/2.
    
    ! Compute the friction velocity u_shear and the drag coefficient

    u_shear=k_roughness*sqrt(ux_HAve**2.+uz_HAve**2.)/log(delta/z_zero)
    CD=k_roughness**2./log(delta/z_zero)**2.

    if (nrank==0) then 
        print *, ' '
        print *, ' ABL:'
        print *, ' Horizontally-averaged velocity at y=1/2... ', ux_HAve,uz_Have
        print *, ' Friction velocity : ', u_shear 
        print *, ' Drag Coefficient : ', CD 
        print *, ' xi1 : ', xi1
        print *, ' scriptR ', scriptR
    endif
    !Compute the shear stresses -- only on the wall
   
    if (xstart(2)==1) then
    do k=1,xsize(3)
    do i=1,xsize(1)                       
    ux12=0.5*(uxf(i,1,k)+uxf(i,2,k)) 
    uz12=0.5*(uzf(i,1,k)+uzf(i,2,k))
    S12=sqrt(ux12**2.+uz12**2.)

    if(iwallmodel==1) then ! MOENG    
    tauwallxy(i,k)=-u_shear**2.0*(S12*ux_HAve+S_HAve*(ux12-ux_HAve))/(S_HAve*sqrt(ux_HAve**2.+uz_HAve**2.))
    tauwallzy(i,k)=-u_shear**2.0*(S12*uz_HAve+S_HAve*(uz12-uz_HAve))/(S_HAve*sqrt(ux_HAve**2.+uz_Have**2.))
    else !Parlage model 
    tauwallxy(i,k)=-(k_roughness/log(delta/z_zero))**2.0*ux12*S12
    tauwallzy(i,k)=-(k_roughness/log(delta/z_zero))**2.0*uz12*S12
    endif
 
    enddo
    enddo
     
    endif
!*********************************************************************************************************

if (nrank==0) write(*,*)  'Maximum wall shear stress for x and z', maxval(tauwallxy), maxval(tauwallzy)
if (nrank==0) write(*,*)  'Minimum wall shear stress for x and z', minval(tauwallxy), minval(tauwallzy)

return

end subroutine wall_shear_stress


subroutine numerical_tripping(Ftrip,x0,lx,ly,time,ts)
    
   USE decomp_2d
   USE param

implicit none
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: Ftrip 
real(mytype),intent(in) :: x0, lx,ly,time,ts
real(mytype) :: x,y,z
integer :: i,j,k

do k=1,xsize(3)
    z=(k+xstart(3)-1-1)*dz
do j=1,xsize(2)
   if (istret.eq.0) y=(j+xstart(2)-1-1)*dy
   if (istret.ne.0) y=yp(j)
do i=1,xsize(1)
    x=(i+xstart(1)-1-1)*dx
    Ftrip=exp(((x-x0)/lx)**2.-(y/ly)**2.)
enddo
enddo
enddo

end subroutine numerical_tripping
