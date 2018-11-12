!==========================================================
! Subroutines for computing the SFS from the ABL
!==========================================================
subroutine wall_sgs(ux,uy,uz,nut1,sxy1,syz1,tauwallxy,tauwallzy,wallfluxx,wallfluxy,wallfluxz)

    USE MPI
    USE decomp_2d
    USE param 
    USE var, only: uxf1,uzf1, di1
implicit none
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz, nut1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: wallfluxx,wallfluxy,wallfluxz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: wallfluxxf,wallfluxyf,wallfluxzf
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxy1, syz1 
real(mytype),dimension(xsize(1),xsize(3)) :: tauwallxy, tauwallzy 
integer :: i,j,k,code
real(mytype) :: ut,ut1,utt,ut11, abl_vel, ABLtaux, ABLtauz, delta
real(mytype) :: ux_HAve_local, uz_HAve_local,S_HAve_local
real(mytype) :: ux_HAve, uz_HAve,S_HAve,ux12,uz12,S12
real(mytype) :: sxy_HAve_local, szy_HAve_local, nut_HAve_local, nutsxy_HAve_local, nutszy_HAve_local
real(mytype) :: sxy_HAve, szy_HAve, nut_HAve, nutsxy_HAve, nutszy_HAve
real(mytype) :: nutprimes
real(mytype) :: nuLESBar, scriptR, xi1, nuLES, TS1, TR1
real(mytype) :: CD ! drag coefficient

call filter(0.0d0)
! Filter the velocity with twice the grid scale according to Bou-zeid et al 2005
call filx(uxf1,ux,di1,fisx,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0) 
call filx(uzf1,uz,di1,fisx,fiffx,fifsx,fifwx,xsize(1),xsize(2),xsize(3),0) 
if (nrank==0) then
open(30,file='filter.dat')
do i=1,xsize(1)
write(30,*) ux(i,1,1), uxf1(i,1,1)
enddo
close(30)
endif

! Determine the shear stress using Moeng's formulation
!*****************************************************************************************
    ux_HAve_local=0.
    uz_HAve_local=0.
    S_HAve_local=0.
    
    if (xstart(2)==1) then 
    !if (nrank==0) print *, 'Max of the filtered velocity', maxval(ux), 'Max of the unfiltered velocity', maxval(ux) 
    do k=1,xsize(3)
    do i=1,xsize(1)
        ux_HAve_local=ux_HAve_local+0.5*(uxf1(i,1,k)+uxf1(i,2,k))
        uz_HAve_local=uz_HAve_local+0.5*(uzf1(i,1,k)+uzf1(i,2,k))
        S_HAve_local=S_HAve_local+sqrt((0.5*(uxf1(i,1,k)+uxf1(i,2,k)))**2.+(0.5*(uzf1(i,1,k)+uzf1(i,2,k)))**2.)
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
    ux12=0.5*(uxf1(i,1,k)+uxf1(i,2,k)) 
    uz12=0.5*(uzf1(i,1,k)+uzf1(i,2,k))
    S12=sqrt(ux12**2.+uz12**2.)
    if(iwallmodel==1) then ! MOENG    
    tauwallxy(i,k)=-(k_roughness/log(delta/z_zero))**2.0*ux_HAve*S_HAve
    tauwallzy(i,k)=-(k_roughness/log(delta/z_zero))**2.0*uz_HAve*S_HAve
    else !Parlage model 
    tauwallxy(i,k)=-(k_roughness/log(delta/z_zero))**2.0*ux12*S12
    tauwallzy(i,k)=-(k_roughness/log(delta/z_zero))**2.0*uz12*S12
    endif
    if(jLES.ge.2) then ! Apply second order central finite difference schemes for the near wall 
    wallfluxx(i,1,k) = -(-1./2.*(-2.*nut1(i,3,k)*sxy1(i,3,k))+&
        2.*(-2.*nut1(i,2,k)*sxy1(i,2,k))-3./2.*tauwallxy(i,k))/(2.*delta)
    wallfluxy(i,1,k) = 0.
    wallfluxz(i,1,k) = -(-1./2.*(-2.*nut1(i,3,k)*syz1(i,3,k))+&
        2.*(-2.*nut1(i,2,k)*syz1(i,2,k))-3./2.*tauwallzy(i,k))/(2.*delta)
    else
    wallfluxx(i,1,k) = -CD*S12*delta 
    wallfluxy(i,1,k) = 0.!
    wallfluxz(i,1,k) = -CD*S12*delta
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

call filter(0.)
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

subroutine tripping(tb,ta) !TRIPPING SUBROUTINE FOR ATMOSPHERIC BOUNDARY LAYERS

USE param
USE var
USE decomp_2d
USE MPI

implicit none

integer :: i,j,k
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb, ta
integer :: seed0, ii, code
real(mytype) :: z_pos, randx, p_tr, b_tr, x_pos, y_pos, A_tr

!Done in X-Pencils
seed0=randomseed !Seed for random number
!A_tr=A_trip*min(1.0,0.8+real(itime)/200.0)
!xs_tr=4.0/2.853
!ys_tr=2.0/2.853
!ts_tr=4.0/2.853
!x0_tr=40.0/2.853
A_tr = 0.1*dt

if ((itime.eq.ifirst).and.(nrank.eq.0)) then
call random_seed(SIZE=ii)
call random_seed(PUT=seed0*(/ (1, i = 1, ii) /))

!DEBUG:
!call random_number(randx)
!call MPI_BCAST(randx,1,real_type,0,MPI_COMM_WORLD,code)
!write(*,*) 'RANDOM:', nrank, randx, ii
!First random generation of h_nxt


  do j=1,z_modes

    call random_number(randx)
    h_coeff(j)=1.0*(randx-0.5)
  enddo
    h_coeff=h_coeff/sqrt(DBLE(z_modes)) 
endif

!Initialization h_nxt  (always bounded by xsize(3)^2 operations)
if (itime.eq.ifirst) then
  call MPI_BCAST(h_coeff,z_modes,real_type,0,MPI_COMM_WORLD,code)
  nxt_itr=0
  do k=1,xsize(3)
     h_nxt(k)=0.0
     z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
   do j=1,z_modes
      h_nxt(k)= h_nxt(k)+h_coeff(j)*sin(2.0*pi*j*z_pos/zlz)
   enddo
  enddo
end if



!Time-loop
 i=int(t/ts_tr)
    if (i.ge.nxt_itr) then  !Nxt_itr is a global variable
        nxt_itr=i+1
        
        !First random generation of h
        h_i(:)=h_nxt(:)
        if (nrank .eq. 0) then
         do j=1,z_modes
            call random_number(randx)
            h_coeff(j)=1.0*(randx-0.5)
         enddo
        h_coeff=h_coeff/sqrt(DBLE(z_modes)) !Non-dimensionalization
        end if
        
        call MPI_BCAST(h_coeff,z_modes,real_type,0,MPI_COMM_WORLD,code)
          
        
        !Initialization h_nxt  (always bounded by z_steps^2 operations)
       do k=1,xsize(3)
          h_nxt(k)=0.0
          z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
          do j=1,z_modes
            h_nxt(k)= h_nxt(k)+h_coeff(j)*sin(2.0*pi*j*z_pos/zlz)
          enddo
        enddo
     endif


 !Time coefficient
  p_tr=t/ts_tr-i
  b_tr=3.0*p_tr**2-2.0*p_tr**3
  
  !Creation of tripping velocity
  do i=1,xsize(1)
    x_pos=(xstart(1)+(i-1)-1)*dx
    do j=1,xsize(2) 
      !y_pos=(xstart(2)+(j-1)-1)*dy   
      y_pos=yp(xstart(2)+(j-1))
      do k=1,xsize(3)
       !g(z)*EXP_F(X,Y)
       ta(i,j,k)=((1.0-b_tr)*h_i(k)+b_tr*h_nxt(k))
       !ta(i,j,k)=A_tr*exp(-((x_pos-x0_tr)/xs_tr)**2-(y_pos/ys_tr)**2)*ta(i,j,k)
       ta(i,j,k)=A_tr*exp(-((x_pos-x0_tr)/xs_tr)**2-((y_pos-0.5)/ys_tr)**2)*ta(i,j,k)  
       tb(i,j,k)=tb(i,j,k)+ta(i,j,k)
       
       z_pos=-zlz/2.0+(xstart(3)+(k-1)-1)*dz
      ! if ((((x_pos-x0_tr)**2).le.9.0e-3).and.(y_pos.le.0.0001).and.((z_pos).le.0.03))then
      !       open(442,file='tripping.dat',form='formatted',position='APPEND')
      !  write(442,*) t,ta(i,j,k)
      !  close(442)   
      ! end if

      enddo
    enddo
  enddo
    
return   
end subroutine tripping

