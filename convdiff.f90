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

!*************************************************************************************
subroutine convdiff(ux1,uy1,uz1,phi1,uxt,uyt,uzt,ep1,divdiva,curldiva,ta1,tb1,tc1,&
     td1,te1,tf1,tg1,th1,ti1,di1,ux2,uy2,uz2,phi2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,&
     ti2,tj2,di2,ux3,uy3,phi3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,nut1,shrt_coeff, &
     ucx1,ucy1,ucz1,tmean,sgszmean,sgsxmean,sgsymean)
!*************************************************************************************

USE param
USE var, only: FTx, FTy, FTz
USE decomp_2d
USE decomp_2d_io
USE MPI


implicit none

! GLOBAL ARRAYS
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1,ep1,ucx1,ucy1,ucz1,deltaphi1
real(mytype),dimension(xsize(1),xsize(2),xsize(3),25) :: uxt,uyt,uzt
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tf1_abl,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2,phi2,deltaphi2 
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tf2_abl,ta2_abl,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3,phi3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,ta3_abl

!LES Arrays
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxx1,syy1,szz1,&
sxy1,sxz1,syz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: asxx1,asyy1,aszz1,&
asxy1,asxz1,asyz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: nut1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sgsx1,sgsy1,sgsz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxx1,gyx1,gzx1,&
gxy1,gyy1,gzy1,gxz1,gyz1,gzz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: srt_smag, shrt2, shrt_coeff
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dsmagcst

!STATISTICS Arrays
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: tmean,sgszmean,sgsxmean,sgsymean
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: divdiva,curldiva
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu

! Buoyancy 
real(mytype),dimension(ysize(2)) :: tmpphi, phiPlaneAve !Horizontally-averaged potential temperature

!ABL boundary conditions
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tablx1, tably1, tablz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: wallfluxx1,wallfluxy1,wallfluxz1,wallfluxphi1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tauwallxy1, tauwallzy1 

integer, dimension(2) :: dims, dummy_coords
logical, dimension(2) :: dummy_periods

integer :: ijk,nvect1,nvect2,nvect3,i,j,k,code
character(len=20) :: filename
real(mytype) :: x,y,z, lambda

ta1=0.;tb1=0.;tc1=0.

nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)

if (iskew==0) then !UROTU!
!WORK X-PENCILS
   call derx (ta1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call derx (tb1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call transpose_x_to_y(ux1,ux2)
   call transpose_x_to_y(uy1,uy2)
   call transpose_x_to_y(uz1,uz2)
   call transpose_x_to_y(ta1,ta2)
   call transpose_x_to_y(tb1,tb2)
!WORK Y-PENCILS
   call dery (tc2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
   call dery (td2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
   call transpose_y_to_z(ux2,ux3)
   call transpose_y_to_z(uy2,uy3)
   call transpose_y_to_z(uz2,uz3)
   call transpose_y_to_z(ta2,ta3)
   call transpose_y_to_z(tb2,tb3)
   call transpose_y_to_z(tc2,tc3)
   call transpose_y_to_z(td2,td3)
!WORK Z-PENCILS
   call derz (te3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (tf3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   do ijk=1,nvect3
      ta3(ijk,1,1)=uz3(ijk,1,1)*(te3(ijk,1,1)-tb3(ijk,1,1))-&
           uy3(ijk,1,1)*(ta3(ijk,1,1)-tc3(ijk,1,1))
      tb3(ijk,1,1)=ux3(ijk,1,1)*(ta3(ijk,1,1)-tc3(ijk,1,1))-&
           uz3(ijk,1,1)*(td3(ijk,1,1)-tf3(ijk,1,1))
      tc3(ijk,1,1)=uy3(ijk,1,1)*(td3(ijk,1,1)-tf3(ijk,1,1))-&
           ux3(ijk,1,1)*(te3(ijk,1,1)-tb3(ijk,1,1))
   enddo
else !SKEW!
!############################## STARTING LES MODELLING HERE #######################

if (jLES==1.and.dynhypvisc==1) then
call shear_rate_coeff(ux1,uy1,uz1,gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,&
sxx1,syy1,szz1,sxy1,sxz1,syz1,shrt2,shrt_coeff,ta2,ta3,di1,di2,di3)
endif

!CLASSIC SMAGORINSKY (plus required rates)
if (jLES==2) then
sgsx1=0.;sgsy1=0.;sgsz1=0.
dsmagcst=0.

call smag(ux1,uy1,uz1,gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,&
sxx1,syy1,szz1,sxy1,sxz1,syz1,srt_smag,nut1,ta2,ta3,di1,di2,di3)

call lesdiff(ux1,uy1,uz1,phi1,gxx1,gyy1,gzz1,gxy1,gxz1,gyz1,gyx1,gzx1,gzy1,nut1,&
    sgsx1,sgsy1,sgsz1,ep1,ta1,td1,te1,tf1,di1,ta2,td2,te2,tf2,tj2,di2,&
    ta3,td3,te3,tf3,di3)
!call compute_sgs(ux1, uy1, uz1, ep1, sxx1, syy1, szz1, sxy1, sxz1, syz1, nut1, &
!     sgsx1, sgsy1, sgsz1, ta1, tb1, tc1, td1, te1, tf1, di1, tb2, tc2, td2, te2, tf2, di2, &
!     tc3, te3, tf3, di3)

elseif (jLES == 3) then !WALE
sgsx1=0.;sgsy1=0.;sgsz1=0.
dsmagcst=0.
    
! First Calculating classic Smagorinsky everywhere in the domain
call smag(ux1,uy1,uz1,gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,&
sxx1,syy1,szz1,sxy1,sxz1,syz1,srt_smag,nut1,ta2,ta3,di1,di2,di3)

! Then adapt the smagorinsky coefficient near the boundaries
call wale(gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,srt_smag,nut1)

call lesdiff(ux1,uy1,uz1,phi1,gxx1,gyy1,gzz1,gxy1,gxz1,gyz1,gyx1,gzx1,gzy1,nut1,&
    sgsx1,sgsy1,sgsz1,ep1,ta1,td1,te1,tf1,di1,ta2,td2,te2,tf2,tj2,di2,&
    ta3,td3,te3,tf3,di3)

elseif (jLES == 4) then !DYNAMIC SMAGORINSKY
sgsx1=0.;sgsy1=0.;sgsz1=0.
dsmagcst=0.

!call dynsmag(ux1,uy1,uz1,ep1,sxx1,syy1,szz1,sxy1,sxz1,syz1,&
!srt_smag,dsmagcst,nut1,di1,ta1,tb1,tc1,td1,ta2,tb2,tc2,td2,te2,tf2,&
!tg2,th2,ti2,di2,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)

call lesdiff(ux1,uy1,uz1,phi1,gxx1,gyy1,gzz1,gxy1,gxz1,gyz1,gyx1,gzx1,gzy1,nut1,&
    sgsx1,sgsy1,sgsz1,ep1,ta1,td1,te1,tf1,di1,ta2,td2,te2,tf2,tj2,di2,&
    ta3,td3,te3,tf3,di3)

endif


!SGS Calculation Over!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ta1=0.;tb1=0.;tc1=0.
td1=0.;te1=0.;tf1=0.

ta2=0.;tb2=0.;tc2=0.
td2=0.;te2=0.;tf2=0.
tj2=0.

ta3=0.;tb3=0.;tc3=0.
td3=0.;te3=0.;tf3=0.
!########################## ENDING LES TERMS ##################################
!SKEW CONVECTIVE TERMS!

   !WORK X-PENCILS
   ta1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)
   tb1(:,:,:)=ux1(:,:,:)*uy1(:,:,:)
   tc1(:,:,:)=ux1(:,:,:)*uz1(:,:,:)
   
   call derx (td1,ta1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call derx (te1,tb1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx (tf1,tc1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
   call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)

   ta1(:,:,:)=0.5*td1(:,:,:)+0.5*ux1(:,:,:)*ta1(:,:,:)
   tb1(:,:,:)=0.5*te1(:,:,:)+0.5*ux1(:,:,:)*tb1(:,:,:)
   tc1(:,:,:)=0.5*tf1(:,:,:)+0.5*ux1(:,:,:)*tc1(:,:,:)

   call transpose_x_to_y(ux1,ux2)
   call transpose_x_to_y(uy1,uy2)
   call transpose_x_to_y(uz1,uz2)
   call transpose_x_to_y(ta1,ta2)
   call transpose_x_to_y(tb1,tb2)
   call transpose_x_to_y(tc1,tc2)

   !WORK Y-PENCILS
   td2(:,:,:)=ux2(:,:,:)*uy2(:,:,:)
   te2(:,:,:)=uy2(:,:,:)*uy2(:,:,:)
   tf2(:,:,:)=uz2(:,:,:)*uy2(:,:,:)
   
   call dery (tg2,td2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   call dery (th2,te2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   call dery (ti2,tf2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   
   ta2(:,:,:)=ta2(:,:,:)+0.5*tg2(:,:,:)+0.5*uy2(:,:,:)*td2(:,:,:)
   tb2(:,:,:)=tb2(:,:,:)+0.5*th2(:,:,:)+0.5*uy2(:,:,:)*te2(:,:,:)
   tc2(:,:,:)=tc2(:,:,:)+0.5*ti2(:,:,:)+0.5*uy2(:,:,:)*tf2(:,:,:)
   
   call transpose_y_to_z(ux2,ux3)
   call transpose_y_to_z(uy2,uy3)
   call transpose_y_to_z(uz2,uz3)
   call transpose_y_to_z(ta2,ta3)
   call transpose_y_to_z(tb2,tb3)
   call transpose_y_to_z(tc2,tc3)
!WORK Z-PENCILS
   
   td3(:,:,:)=ux3(:,:,:)*uz3(:,:,:)
   te3(:,:,:)=uy3(:,:,:)*uz3(:,:,:)
   tf3(:,:,:)=uz3(:,:,:)*uz3(:,:,:)
   
   call derz (tg3,td3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
   call derz (th3,te3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
   call derz (ti3,tf3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (td3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (te3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
   call derz (tf3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

   ta3(:,:,:)=ta3(:,:,:)+0.5*tg3(:,:,:)+0.5*uz3(:,:,:)*td3(:,:,:)
   tb3(:,:,:)=tb3(:,:,:)+0.5*th3(:,:,:)+0.5*uz3(:,:,:)*te3(:,:,:)
   tc3(:,:,:)=tc3(:,:,:)+0.5*ti3(:,:,:)+0.5*uz3(:,:,:)*tf3(:,:,:)

endif
!ALL THE CONVECTIVE TERMS ARE IN TA3, TB3 and TC3

td3(:,:,:)=ta3(:,:,:)
te3(:,:,:)=tb3(:,:,:)
tf3(:,:,:)=tc3(:,:,:)

!DIFFUSIVE TERMS IN Z
if (jLES==1) then ! IMPLICIT LES
    call derzz_iles (ta3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    call derzz_iles (tb3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    call derzz_iles (tc3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)
else
    call derzz (ta3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    call derzz (tb3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    call derzz (tc3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)
endif

!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2)
call transpose_z_to_y(tb3,tb2)
call transpose_z_to_y(tc3,tc2)
call transpose_z_to_y(td3,td2)
call transpose_z_to_y(te3,te2)
call transpose_z_to_y(tf3,tf2)

tg2(:,:,:)=td2(:,:,:)
th2(:,:,:)=te2(:,:,:)
ti2(:,:,:)=tf2(:,:,:)

!DIFFUSIVE TERMS IN Y
!-->for ux
if (istret.ne.0) then
    if(jLES==1) then
        call deryy_iles (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    else
        call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    endif
   call dery (te2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      td2(i,j,k)=td2(i,j,k)*pp2y(j)-pp4y(j)*te2(i,j,k)
   enddo
   enddo
   enddo
else
    if(jLES==1) then
        call deryy_iles (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    else
        call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    endif
   !call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
endif
!-->for uy
if (istret.ne.0) then
    if(jLES==1) then
   call deryy_iles (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
   else
   call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
    endif
   call dery (tf2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      te2(i,j,k)=te2(i,j,k)*pp2y(j)-pp4y(j)*tf2(i,j,k)
   enddo
   enddo
   enddo
else
    if(jLES==1) then
        call deryy_iles (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
    else
        call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
    endif
endif
!-->for uz
if (istret.ne.0) then
    if(jLES==1) then
        call deryy_iles (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    else
        call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    endif
    call dery (tj2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      tf2(i,j,k)=tf2(i,j,k)*pp2y(j)-pp4y(j)*tj2(i,j,k)
   enddo
   enddo
   enddo
else
    if(jLES==1) then
        call deryy_iles (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
    else 
        call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    endif
endif

ta2(:,:,:)=ta2(:,:,:)+td2(:,:,:)
tb2(:,:,:)=tb2(:,:,:)+te2(:,:,:)
tc2(:,:,:)=tc2(:,:,:)+tf2(:,:,:)

!WORK X-PENCILS
call transpose_y_to_x(ta2,ta1)
call transpose_y_to_x(tb2,tb1)
call transpose_y_to_x(tc2,tc1) !diff
call transpose_y_to_x(tg2,td1)
call transpose_y_to_x(th2,te1)
call transpose_y_to_x(ti2,tf1) !conv

tg1(:,:,:)=td1(:,:,:)
th1(:,:,:)=te1(:,:,:)
ti1(:,:,:)=tf1(:,:,:)

!DIFFUSIVE TERMS IN X
if(jLES==1) then
call derxx_iles (td1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
call derxx_iles (te1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
call derxx_iles (tf1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
else
call derxx (td1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
call derxx (te1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
call derxx (tf1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
endif


ta1(:,:,:)=ta1(:,:,:)+td1(:,:,:)
tb1(:,:,:)=tb1(:,:,:)+te1(:,:,:)
tc1(:,:,:)=tc1(:,:,:)+tf1(:,:,:)

if (iabl==1) then
    ! In case of ABL set to zero the SGS model at level 1 (This will be computed later by the 
    ! SGS wall stress model    
    if (xstart(2)==1) then
    sgsx1(:,1,:)=0.
    sgsy1(:,1,:)=0.
    sgsz1(:,1,:)=0.
    endif
endif

!FINAL SUM: DIFF TERMS + CONV TERMS
if(jLES==0) then ! DNS 
    ta1(:,:,:)=xnu*ta1(:,:,:)-tg1(:,:,:)
    tb1(:,:,:)=xnu*tb1(:,:,:)-th1(:,:,:)
    tc1(:,:,:)=xnu*tc1(:,:,:)-ti1(:,:,:)
elseif (jLES==1) then ! implicit LES 
    if(dynhypvisc==1) then ! use the dynamic definition of nu0/nu 
    ta1(:,:,:)=xnu*shrt_coeff(:,:,:)*ta1(:,:,:)-tg1(:,:,:)
    tb1(:,:,:)=xnu*shrt_coeff(:,:,:)*tb1(:,:,:)-th1(:,:,:)
    tc1(:,:,:)=xnu*shrt_coeff(:,:,:)*tc1(:,:,:)-ti1(:,:,:)
    else
    ta1(:,:,:)=xnu*ta1(:,:,:)-tg1(:,:,:)
    tb1(:,:,:)=xnu*tb1(:,:,:)-th1(:,:,:)
    tc1(:,:,:)=xnu*tc1(:,:,:)-ti1(:,:,:)
    endif
elseif (jLES==2) then ! Classic Smagorisnky Model
    ta1(:,:,:)=xnu*ta1(:,:,:)-tg1(:,:,:)+sgsx1(:,:,:)
    tb1(:,:,:)=xnu*tb1(:,:,:)-th1(:,:,:)+sgsy1(:,:,:)
    tc1(:,:,:)=xnu*tc1(:,:,:)-ti1(:,:,:)+sgsz1(:,:,:)
elseif (jLES==3) then ! WALE 
    ta1(:,:,:)=xnu*ta1(:,:,:)-tg1(:,:,:)+sgsx1(:,:,:)
    tb1(:,:,:)=xnu*tb1(:,:,:)-th1(:,:,:)+sgsy1(:,:,:)
    tc1(:,:,:)=xnu*tc1(:,:,:)-ti1(:,:,:)+sgsz1(:,:,:) 
elseif(jLES==4) then ! dynamic Smagorinsky 
    ta1(:,:,:)=xnu*ta1(:,:,:)-tg1(:,:,:)+sgsx1(:,:,:)
    tb1(:,:,:)=xnu*tb1(:,:,:)-th1(:,:,:)+sgsy1(:,:,:)
    tc1(:,:,:)=xnu*tc1(:,:,:)-ti1(:,:,:)+sgsz1(:,:,:)
else
    if(nrank==0) then
    write(*,*) 'Dont know what to do. This LES model is not defined'
    write(*,*) 'Choose between : 0--> DNS'
    write(*,*) '               : 1--> Implicit SVV'
    write(*,*) '               : 2--> standard Smagorinsky'
    write(*,*) '               : 3--> Wall-Adaptive LES'
    write(*,*) '               : 4--> Scale-invariant dynamic smagorinsky model'
    endif
    stop
endif
!***************************************
! Compute additional Models
!***************************************
! 
if (iabl==1) then
    call wall_sgs(ux1,uy1,uz1,nut1,sxy1,syz1,tauwallxy1,tauwallzy1,wallfluxx1,wallfluxy1,wallfluxz1)
    ta1(:,:,:)=ta1(:,:,:)+wallfluxx1(:,:,:)
    tb1(:,:,:)=tb1(:,:,:)+wallfluxy1(:,:,:)
    tc1(:,:,:)=tc1(:,:,:)+wallfluxz1(:,:,:)
endif

! Buoyancy Effects
if (ibuoyancy==1) then     
    ! Average quantities over the x-z plane 
    deltaphi1(:,:,:)=(phi1(:,:,:)-TempRef)/TempRef
    tb1(:,:,:)=tb1(:,:,:) + 9.81*deltaphi1(:,:,:)
endif

if (IPressureGradient==1) then
    ta1(:,:,:)=ta1(:,:,:)+ustar**2./yly ! Apply a pressure gradient in the stream-wise direction
endif

! Coriolis Effects
if (icoriolis==1) then    
    ta1(:,:,:)=ta1(:,:,:)+CoriolisFreq*uz1(:,:,:) ! This is the stream-wise direction
    tc1(:,:,:)=tc1(:,:,:)-CoriolisFreq*ux1(:,:,:) ! This is not the vertical direction but the lateral horizontal
endif

if (idampingzone==1) then
do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)

if (istret.eq.0) y=(j+xstart(2)-1-1)*dy
if (istret.ne.0) y=yp(j+xstart(2)-1)

if (y>1.1*dBL) then
    lambda=1.0
elseif (y>0.9*dBL.and.y<1.1*dBL) then
    lambda=0.5*(1-cos(pi*(y-0.9*dBL)/(0.2*dBL)))
else 
    lambda=0.
endif
ta1(i,j,k)=ta1(i,j,k)-5.*ustar/yly*lambda*(ux1(i,j,k)-UG(1))
tb1(i,j,k)=tb1(i,j,k)-5.*ustar/yly*lambda*(uy1(i,j,k)-UG(2))
tc1(i,j,k)=tc1(i,j,k)-5.*ustar/yly*lambda*(uz1(i,j,k)-UG(3))
enddo
enddo
enddo
endif

! Turbine forcing through an Actuator Line Model
if (ialm==1) then
    ta1(:,:,:)=ta1(:,:,:)+FTx(:,:,:)
    tb1(:,:,:)=tb1(:,:,:)+FTy(:,:,:)
    tc1(:,:,:)=tc1(:,:,:)+FTz(:,:,:)
endif


end subroutine convdiff

!************************************************************
subroutine PotentialTemperature(ux1,uy1,uz1,nut1,phi1,phis1,phiss1,di1,ta1,tb1,tc1,td1,&
     uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,tc3,epsi)
!************************************************************

USE param
USE decomp_2d

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,nut1,phi1,phis1,&
                                              phiss1,di1,ta1,tb1,tc1,td1,epsi,sgsphi1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uy2,uz2,phi2,nut2,di2,ta2,tb2,tc2,td2,sgsphi2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,phi3,nut3,di3,ta3,tb3,tc3,td3,sgsphi3

integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz
real(mytype) :: x,y,z, PsiH, delta

!X PENCILS
ta1(:,:,:)=ux1(:,:,:)*phi1(:,:,:)

call derx (tb1,ta1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tc1,nut1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (td1,phi1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)

if (jles==1) then
call derxx_iles (ta1,phi1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
else

call derxx (ta1,phi1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
sgsphi1=tc1/Pr*td1+nut1/Pr*ta1
endif

call transpose_x_to_y(phi1,phi2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)
!For explicit LES
call transpose_x_to_y(nut1,nut2)
call transpose_x_to_y(sgsphi1,sgsphi2)

!Y PENCILS
ta2(:,:,:)=uy2(:,:,:)*phi2(:,:,:)

call dery(tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery(td2,nut2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)

if (istret.ne.0) then         
    call dery(tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    
    if (jles==1) then
        call deryy_iles (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    else
        call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
    endif
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      ta2(i,j,k)=ta2(i,j,k)*pp2y(j)-pp4y(j)*tc2(i,j,k)
   enddo
   enddo
   enddo
        sgsphi2=sgsphi2+td2/Pr+nut1/Pr*ta2
else
    if (jles==1) then
        call deryy_iles (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
    else
        call deryy (ta2,phi2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1) 
        sgsphi2=sgsphi2+td2/Pr*tc2+nut2/Pr*ta2
    endif
endif
!Before applying the SGS model -- do the correction for ABL
if(iabl==1.and.jles.ge.2.and.ysize(2).eq.1) then
if (istret.ne.0) delta=(yp(2)-yp(1))/2.0
if (istret.eq.0) delta=dy/2.
do k=1,ysize(3)
do i=1,ysize(1)
PsiH=0!-7.8*delta/ObukhovL
sgsphi2(i,1,k) = -(-1./2.*(-2.*nut2(i,3,k)*tc2(i,3,k))/Pr+&
        2.*(-2.*nut2(i,2,k)*tc2(i,2,k))/Pr-3./2.*(u_shear*k_roughness*(phi2(i,1,k)-0.5*(phi2(i,1,k)+phi2(i,2,k)))/(log(delta/z_zero)-PsiH)))/(2.*delta)
    enddo
    enddo
endif

call transpose_y_to_z(phi2,phi3)
call transpose_y_to_z(uz2,uz3)
!For explicit LES
call transpose_y_to_z(nut2,nut3)
call transpose_y_to_z(sgsphi2,sgsphi3)

!Z PENCILS
ta3(:,:,:)=uz3(:,:,:)*phi3(:,:,:)

call derz(tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
call derz(tc3,nut3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
call derz(td3,phi3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

if (jles==1) then
    call derzz_iles (ta3,phi3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
else
    call derzz (ta3,phi3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
    sgsphi3=sgsphi3+tc3/Pr*td3+nut3/Pr*ta3
endif

call transpose_z_to_y(ta3,tc2)
call transpose_z_to_y(tb3,td2)

!Y PENCILS ADD TERMS
tc2(:,:,:)=tc2(:,:,:)+ta2(:,:,:)
td2(:,:,:)=td2(:,:,:)+tb2(:,:,:)

call transpose_y_to_x(tc2,tc1)
call transpose_y_to_x(td2,td1)

call transpose_z_to_y(sgsphi3,sgsphi2)
call transpose_y_to_x(sgsphi2,sgsphi1)

!X PENCILS ADD TERMS
ta1(:,:,:)=ta1(:,:,:)+tc1(:,:,:) !SECOND DERIVATIVE
tb1(:,:,:)=tb1(:,:,:)+td1(:,:,:) !FIRST DERIVATIVE
 
if(jles.le.1) then
ta1(:,:,:)=xnu/Pr*ta1(:,:,:)-tb1(:,:,:) 
else
ta1(:,:,:)=xnu/Pr*ta1(:,:,:)-tb1(:,:,:)+sgsphi1(:,:,:) 
endif



!TIME ADVANCEMENT (Doing it locally)
if ((nscheme.eq.1).or.(nscheme.eq.2)) then
   if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
        (nscheme.eq.2.and.itr.eq.1)) then
        phi1(:,:,:)=gdt(itr)*ta1(:,:,:)+phi1(:,:,:)
        phis1(:,:,:)=ta1(:,:,:)          
   else
         phi1(:,:,:)=adt(itr)*ta1(:,:,:)+bdt(itr)*phis1(:,:,:)+phi1(:,:,:)
         phis1(:,:,:)=ta1(:,:,:)          
endif
endif

if (nscheme.eq.3) then 
endif

if (nscheme==4) then
   if ((itime.eq.1).and.(ilit.eq.0)) then
      if (nrank==0) print *,'start with Euler',itime
      do ijk=1,nxyz !start with Euler
         phi1(ijk,1,1)=dt*ta1(ijk,1,1)+phi1(ijk,1,1)
         phis1(ijk,1,1)=ta1(ijk,1,1)          
      enddo
   else
      if  ((itime.eq.2).and.(ilit.eq.0)) then
         if (nrank==0) print *,'then with AB2',itime
         do ijk=1,nxyz
            phi1(ijk,1,1)=1.5*dt*ta1(ijk,1,1)-0.5*dt*phis1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo 
      else
         do ijk=1,nxyz
            phi1(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*phis1(ijk,1,1)+&
                 cdt(itr)*phiss1(ijk,1,1)+phi1(ijk,1,1)
            phiss1(ijk,1,1)=phis1(ijk,1,1)
            phis1(ijk,1,1)=ta1(ijk,1,1)
         enddo
      endif
   endif
endif
end subroutine PotentialTemperature

