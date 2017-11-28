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

subroutine init_explicit_les

    USE param
    USE variables
    USE decomp_2d

    implicit none
    integer:: j
    if(nrank==0) then
    write(*,*) ' '
    write(*,*) '++++++++++++++++++++++++++++++++'
    write(*,*) 'Initializing explicit LES Filter'
        if(jLES==2) then
        write(*,*) ' Classic Smagorinsky is used ... '
        write(*,*) ' Smagorinsky constant = ', smagcst
        write(*,*) ' Filter Size / Grid Size = ', FSGS
        else if (jLES==3) then
        write(*,*) ' Wall-adaptive LES (WALES) is used ... '
        else if (jLES==4) then
        call filter() ! Apply the first filter to the equations
        write(*,*) ' Scale-invariant Dynamic Smagorinsky is used ... '
        else if (jLES==5) then
        call filter() ! Apply the first filter to the equations
        write(*,*) ' Scale-dependent Dynamic Smagorinsky is used ... '
        endif
    write(*,*) '++++++++++++++++++++++++++++++++'
    write(*,*) ' '
    endif 

    if (istret.eq.0) then 
    del(:)=FSGS*(dx*dy*dz)**(1.0/3.0)
    else
        do j=1,ny-1
            del(j) = FSGS*(dx*(yp(j+1)-yp(j))*dz)**(1.0/3.0)
        enddo
        del(ny) = del(ny-1)
    endif
    
if (nrank==0) write(*,*) maxval(del)

end subroutine


!************************************************************
!
subroutine smag(ux1,uy1,uz1,gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,&
    sxx1,syy1,szz1,sxy1,sxz1,syz1,srt_smag,nut1,ta2,ta3,di1,di2,di3)
!
!************************************************************
USE param
USE variables
USE decomp_2d

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,di3

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxx1,syy1,szz1,&
sxy1,sxz1,syz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxx1,gyx1,gzx1,&
gxy1,gyy1,gzy1,gxz1,gyz1,gzz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: srt_smag,nut1

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: syy2,szz2,&
sxy2,syz2,srt_smag2,nut2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gxy2,gyy2,gzy2,&
gxz2,gyz2,gzz2

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: szz3,syz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: gxz3,gyz3,gzz3

real(mytype) :: smag_constant, y, yplus, delta, length

integer :: i,j,k

srt_smag(:,:,:) = 0.

call derx (gxx1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (gyx1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (gzx1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)

sxx1(:,:,:) = gxx1(:,:,:)

!WORK Y-PENCILS
call transpose_x_to_y(ux1,ux2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)
call transpose_x_to_y(gyx1,ta2)

call dery (gxy2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (gyy2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (gzy2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)

sxy2(:,:,:)=0.5*(gxy2(:,:,:)+ta2(:,:,:))
syy2(:,:,:)=gyy2(:,:,:)

!WORK Z-PENCILS
call transpose_y_to_z(ux2,ux3)
call transpose_y_to_z(uy2,uy3)
call transpose_y_to_z(uz2,uz3)
call transpose_y_to_z(gzy2,ta3)

call derz(gxz3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz(gyz3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz(gzz3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

szz3(:,:,:)=gzz3(:,:,:)
syz3(:,:,:)=0.5*(gyz3(:,:,:)+ta3(:,:,:))

!WORK Y-PENCILS
call transpose_z_to_y(syz3,syz2)
call transpose_z_to_y(szz3,szz2)

call transpose_z_to_y(gxz3,gxz2)
call transpose_z_to_y(gyz3,gyz2)
call transpose_z_to_y(gzz3,gzz2)

!WORK X-PENCILS
call transpose_y_to_x(sxy2,sxy1)
call transpose_y_to_x(syy2,syy1)
call transpose_y_to_x(syz2,syz1)
call transpose_y_to_x(szz2,szz1)

call transpose_y_to_x(gxy2,gxy1)
call transpose_y_to_x(gyy2,gyy1)
call transpose_y_to_x(gzy2,gzy1)
call transpose_y_to_x(gxz2,gxz1)
call transpose_y_to_x(gyz2,gyz1)
call transpose_y_to_x(gzz2,gzz1)

sxz1(:,:,:)=0.5*(gzx1(:,:,:)+gxz1(:,:,:))

srt_smag(:,:,:) = sxx1(:,:,:)*sxx1(:,:,:)+syy1(:,:,:)*syy1(:,:,:) &
+szz1(:,:,:)*szz1(:,:,:)+2.*sxy1(:,:,:)*sxy1(:,:,:) &
+2.*sxz1(:,:,:)*sxz1(:,:,:)+2.*syz1(:,:,:)*syz1(:,:,:)

!if (nrank==0) print *, "srt_smag = ",maxval(srt_smag)

call transpose_x_to_y(srt_smag,srt_smag2)
do k=1,ysize(3)
do j=1,ysize(2)
do i=1,ysize(1)

if(SmagWallDamp.eq.1) then
    ! Mason & Thomson damping coefficients
    if (istret.eq.0) y=(j+ystart(2)-1-1)*dy
    if (istret.ne.0) y=yp(j+ystart(2)-1)
    smag_constant=(smagcst**(-nSmag)+(k_roughness*(y/del(j)+z_zero/del(j)))**(-nSmag))**(-1./nSmag)
    length=smag_constant*del(j)
else if(SmagWallDamp.eq.2) then
    ! van Driest damping coefficient 
    if (istret.eq.0) y=(j+ystart(2)-1-1)*dy-yly/2.
    if (istret.ne.0) y=yp(j+ystart(2)-1)-yly/2.
    delta=u_shear/xnu
    yplus=y/delta
    length=smagcst*del(j)*sqrt((1.-exp(yplus/25.))**3.)
else
    length=smagcst*del(j)
endif
    nut2(i,j,k) = length**2.0*sqrt(2.*srt_smag2(i,j,k))
enddo
enddo
enddo
call transpose_y_to_x(nut2,nut1)

end subroutine smag

!************************************************************
!
subroutine wale(gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,srt_smag,nut1)
!
!************************************************************
USE param
USE variables
USE decomp_2d


implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxx1,gyx1,gzx1,&
gxy1,gyy1,gzy1,gxz1,gyz1,gzz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: srt_wale,srt_smag
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: srt_wale2,srt_smag2,nut2
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxxd1,syyd1,szzd1,&
sxyd1,sxzd1,syzd1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: nut1

integer::i,j,k

do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)

sxxd1(i,j,k)=gxx1(i,j,k)*gxx1(i,j,k)+gxy1(i,j,k)*gyx1(i,j,k)+&
gxz1(i,j,k)*gzx1(i,j,k)-(1./3.)*(gxx1(i,j,k)*gxx1(i,j,k)+&
gyy1(i,j,k)*gyy1(i,j,k)+gzz1(i,j,k)*gzz1(i,j,k)+&
2.*gxy1(i,j,k)*gyx1(i,j,k)+2.*gxz1(i,j,k)*gzx1(i,j,k)+&
2.*gzy1(i,j,k)*gyz1(i,j,k))

sxyd1(i,j,k)=0.5*(gxx1(i,j,k)*gxy1(i,j,k)+gxy1(i,j,k)*gyy1(i,j,k)+&
gxz1(i,j,k)*gzy1(i,j,k)+gyx1(i,j,k)*gxx1(i,j,k)+&
gyy1(i,j,k)*gyx1(i,j,k)+gyz1(i,j,k)*gzx1(i,j,k))

sxzd1(i,j,k)=0.5*(gxx1(i,j,k)*gxz1(i,j,k)+gxy1(i,j,k)*gyz1(i,j,k)+&
gxz1(i,j,k)*gzz1(i,j,k)+gzx1(i,j,k)*gxx1(i,j,k)+&
gzy1(i,j,k)*gyx1(i,j,k)+gzz1(i,j,k)*gzx1(i,j,k))

syyd1(i,j,k)=gyx1(i,j,k)*gxy1(i,j,k)+gyy1(i,j,k)*gyy1(i,j,k)+&
gyz1(i,j,k)*gzy1(i,j,k)-(1./3.)*(gxx1(i,j,k)*gxx1(i,j,k)+&
gyy1(i,j,k)*gyy1(i,j,k)+gzz1(i,j,k)*gzz1(i,j,k)+&
2.*gxy1(i,j,k)*gyx1(i,j,k)+2.*gxz1(i,j,k)*gzx1(i,j,k)+&
2.*gzy1(i,j,k)*gyz1(i,j,k))

syzd1(i,j,k)=0.5*(gyx1(i,j,k)*gxz1(i,j,k)+gyy1(i,j,k)*gyz1(i,j,k)+&
gyz1(i,j,k)*gzz1(i,j,k)+gzx1(i,j,k)*gxy1(i,j,k)+&
gzy1(i,j,k)*gyy1(i,j,k)+gzz1(i,j,k)*gzy1(i,j,k))

szzd1(i,j,k)=gzx1(i,j,k)*gxz1(i,j,k)+gzy1(i,j,k)*gyz1(i,j,k)+&
gzz1(i,j,k)*gzz1(i,j,k)-(1./3.)*(gxx1(i,j,k)*gxx1(i,j,k)+&
gyy1(i,j,k)*gyy1(i,j,k)+gzz1(i,j,k)*gzz1(i,j,k)+&
2.*gxy1(i,j,k)*gyx1(i,j,k)+2.*gxz1(i,j,k)*gzx1(i,j,k)+&
2.*gzy1(i,j,k)*gyz1(i,j,k))

srt_wale(i,j,k) = sxxd1(i,j,k)*sxxd1(i,j,k)+syyd1(i,j,k)*syyd1(i,j,k)+&
szzd1(i,j,k)*szzd1(i,j,k)+ 2.*sxyd1(i,j,k)*sxyd1(i,j,k)+ &
2.*sxzd1(i,j,k)*sxzd1(i,j,k)+ 2.*syzd1(i,j,k)*syzd1(i,j,k)

enddo
enddo
enddo

call transpose_x_to_y(srt_smag,srt_smag2)
call transpose_x_to_y(srt_wale,srt_wale2)

do k=1,ysize(3)
do j=1,ysize(2)
do i=1,ysize(1)
nut2(i,j,k)=((walecst*del(j))**(2.0)*srt_wale2(i,j,k)**(3./2.)) &
/(srt_smag2(i,j,k)**(5./2.)+srt_wale2(i,j,k)**(5./4.))
enddo
enddo
enddo

call transpose_y_to_x(nut2,nut1)

end subroutine wale

!************************************************************
!
subroutine dynsmag(ux1,uy1,uz1,ep1,sxx1,syy1,szz1,sxy1,sxz1,syz1,&
srt_smag,dsmagcst,nut1,di1,ta1,tb1,tc1,td1,ta2,tb2,tc2,td2,te2,tf2,&
tg2,th2,ti2,di2,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)
!
!************************************************************
USE param
USE variables
USE decomp_2d
USE decomp_2d_io
USE MPI


implicit none

TYPE(DECOMP_INFO) :: ph3
!real(mytype),dimension(1:nx,1:ny,1:nz) :: smagCgbl
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: di1,ta1,tb1,tc1,td1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxx1,syy1,szz1,&
sxy1,sxz1,syz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: nut1,srt_smag
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: nut2,srt_smag2
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxx1,uyy1,uzz1,&
uxy1,uxz1,uyz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: uxx1f,uyy1f,uzz1f,&
uxy1f,uxz1f,uyz1f,ux1f,uy1f,uz1f
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: lxx1,lyy1,lzz1,&
lxy1,lxz1,lyz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxx1f,gyx1f,gzx1f,&
gxz1f
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxx1f,syy1f,szz1f,&
sxy1f,sxz1f,syz1f
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: bbxx1,bbyy1,bbzz1,&
bbxy1,bbxz1,bbyz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: axx1,ayy1,azz1,&
axy1,axz1,ayz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: axx1f,ayy1f,azz1f,&
axy1f,axz1f,ayz1f
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: mxx1,myy1,mzz1,&
mxy1,mxz1,myz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: smagC,smagC1f

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uxx2f,uyy2f,uzz2f,&
uxy2f,uxz2f,uyz2f,ux2f,uy2f,uz2f
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gxy2f,gyy2f,gzy2f,&
gxz2f
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: syy2f,szz2f,&
sxy2f,syz2f
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: axx2f,ayy2f,azz2f,&
axy2f,axz2f,ayz2f
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: smagC2f

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uxx3f,uyy3f,uzz3f,&
uxy3f,uxz3f,uyz3f,ux3f,uy3f,uz3f
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: gxz3f,gyz3f,gzz3f
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: szz3f,syz3f
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: axx3f,ayy3f,azz3f,&
axy3f,axz3f,ayz3f
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: smagC3f
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: dsmagcst3
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: dsmagcst2
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: dsmagcst
real(mytype),dimension(xsize(2)) :: tmpa1, dsmagHP1
real(mytype) :: smagC1,smagC2,dsmaggbl
integer::i,j,k,ierror,i2,j2,code
real(mytype) :: denom ! Denominator in dynamic Smagorinsky


uxx1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)
uyy1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)
uzz1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)
uxy1(:,:,:)=ux1(:,:,:)*uy1(:,:,:)
uxz1(:,:,:)=ux1(:,:,:)*uz1(:,:,:)
uyz1(:,:,:)=uy1(:,:,:)*uz1(:,:,:)

call filx(ux1f,ux1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(uy1f,uy1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(uz1f,uz1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(uxx1f,uxx1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(uyy1f,uyy1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(uzz1f,uzz1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(uxy1f,uxy1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(uxz1f,uxz1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(uyz1f,uyz1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)

call transpose_x_to_y(ux1f,ta2)
call transpose_x_to_y(uy1f,tb2)
call transpose_x_to_y(uz1f,tc2)
call transpose_x_to_y(uxx1f,td2)
call transpose_x_to_y(uyy1f,te2)
call transpose_x_to_y(uzz1f,tf2)
call transpose_x_to_y(uxy1f,tg2)
call transpose_x_to_y(uxz1f,th2)
call transpose_x_to_y(uyz1f,ti2)

call fily(ux2f,ta2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(uy2f,tb2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(uz2f,tc2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(uxx2f,td2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(uyy2f,te2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(uzz2f,tf2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(uxy2f,tg2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(uxz2f,th2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(uyz2f,ti2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)

ta2=0.;tb2=0.;tc2=0.
td2=0.;te2=0.;tf2=0.
tg2=0.;th2=0.;ti2=0.

!
call transpose_y_to_z(ux2f,ta3)
call transpose_y_to_z(uy2f,tb3)
call transpose_y_to_z(uz2f,tc3)
call transpose_y_to_z(uxx2f,td3)
call transpose_y_to_z(uyy2f,te3)
call transpose_y_to_z(uzz2f,tf3)
call transpose_y_to_z(uxy2f,tg3)
call transpose_y_to_z(uxz2f,th3)
call transpose_y_to_z(uyz2f,ti3)


call filz(ux3f,ta3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(uy3f,tb3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(uz3f,tc3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(uxx3f,td3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(uyy3f,te3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(uzz3f,tf3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(uxy3f,tg3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(uxz3f,th3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(uyz3f,ti3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)

ta3=0.;tb3=0.;tc3=0.
td3=0.;te3=0.;tf3=0.
tg3=0.;th3=0.;ti3=0.

ux2f=0.;uy2f=0.;uz2f=0.
call transpose_z_to_y(ux3f,ux2f)
call transpose_z_to_y(uy3f,uy2f)
call transpose_z_to_y(uz3f,uz2f)
call transpose_z_to_y(uxx3f,uxx2f)
call transpose_z_to_y(uyy3f,uyy2f)
call transpose_z_to_y(uzz3f,uzz2f)
call transpose_z_to_y(uxy3f,uxy2f)
call transpose_z_to_y(uxz3f,uxz2f)
call transpose_z_to_y(uyz3f,uyz2f)

ux1f=0.;uy1f=0.;uz1f=0.
call transpose_y_to_x(ux2f,ux1f)
call transpose_y_to_x(uy2f,uy1f)
call transpose_y_to_x(uz2f,uz1f)
call transpose_y_to_x(uxx2f,uxx1f)
call transpose_y_to_x(uyy2f,uyy1f)
call transpose_y_to_x(uzz2f,uzz1f)
call transpose_y_to_x(uxy2f,uxy1f)
call transpose_y_to_x(uxz2f,uxz1f)
call transpose_y_to_x(uyz2f,uyz1f)

!Lij tensor OK
lxx1(:,:,:)=uxx1f(:,:,:)-ux1f(:,:,:)*ux1f(:,:,:)
lyy1(:,:,:)=uyy1f(:,:,:)-uy1f(:,:,:)*uy1f(:,:,:)
lzz1(:,:,:)=uzz1f(:,:,:)-uz1f(:,:,:)*uz1f(:,:,:)
lxy1(:,:,:)=uxy1f(:,:,:)-ux1f(:,:,:)*uy1f(:,:,:)
lxz1(:,:,:)=uxz1f(:,:,:)-ux1f(:,:,:)*uz1f(:,:,:)
lyz1(:,:,:)=uyz1f(:,:,:)-uy1f(:,:,:)*uz1f(:,:,:)
!

call derx (gxx1f,ux1f,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (gyx1f,uy1f,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (gzx1f,uz1f,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)

sxx1f(:,:,:) = gxx1f(:,:,:)

!WORK Y-PENCILS
call transpose_x_to_y(ux1f,ux2f)
call transpose_x_to_y(uy1f,uy2f)
call transpose_x_to_y(uz1f,uz2f)
call transpose_x_to_y(gyx1f,ta2)

call dery (gxy2f,ux2f,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (gyy2f,uy2f,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (gzy2f,uz2f,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)

sxy2f(:,:,:)=0.5*(gxy2f(:,:,:)+ta2(:,:,:))
syy2f(:,:,:)=gyy2f(:,:,:)

!WORK Z-PENCILS
call transpose_y_to_z(ux2f,ux3f)
call transpose_y_to_z(uy2f,uy3f)
call transpose_y_to_z(uz2f,uz3f)
call transpose_y_to_z(gzy2f,ta3)

call derz(gxz3f,ux3f,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
call derz(gyz3f,uy3f,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
call derz(gzz3f,uz3f,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

szz3f(:,:,:)=gzz3f(:,:,:)
syz3f(:,:,:)=0.5*(gyz3f(:,:,:)+ta3(:,:,:))

!WORK Y-PENCILS
call transpose_z_to_y(syz3f,syz2f)
call transpose_z_to_y(szz3f,szz2f)

call transpose_z_to_y(gxz3f,gxz2f)

!WORK X-PENCILS
call transpose_y_to_x(sxy2f,sxy1f)
call transpose_y_to_x(syy2f,syy1f)
call transpose_y_to_x(syz2f,syz1f)
call transpose_y_to_x(szz2f,szz1f)

call transpose_y_to_x(gxz2f,gxz1f)

sxz1f(:,:,:)=0.5*(gzx1f(:,:,:)+gxz1f(:,:,:))

do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)

!Bij tensor with u test filtered OK
bbxx1(i,j,k)=-2*sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
syy1f(i,j,k)*syy1f(i,j,k)+&
szz1f(i,j,k)*szz1f(i,j,k)+&
2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
2.*syz1f(i,j,k)*syz1f(i,j,k)))&
*sxx1f(i,j,k)

bbyy1(i,j,k)=-2*sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
syy1f(i,j,k)*syy1f(i,j,k)+&
szz1f(i,j,k)*szz1f(i,j,k)+&
2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
2.*syz1f(i,j,k)*syz1f(i,j,k)))&
*syy1f(i,j,k)

bbzz1(i,j,k)=-2*sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
syy1f(i,j,k)*syy1f(i,j,k)+&
szz1f(i,j,k)*szz1f(i,j,k)+&
2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
2.*syz1f(i,j,k)*syz1f(i,j,k)))&
*szz1f(i,j,k)

bbxy1(i,j,k)=-2*sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
syy1f(i,j,k)*syy1f(i,j,k)+&
szz1f(i,j,k)*szz1f(i,j,k)+&
2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
2.*syz1f(i,j,k)*syz1f(i,j,k)))&
*sxy1f(i,j,k)

bbxz1(i,j,k)=-2*sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
syy1f(i,j,k)*syy1f(i,j,k)+&
szz1f(i,j,k)*szz1f(i,j,k)+&
2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
2.*syz1f(i,j,k)*syz1f(i,j,k)))&
*sxz1f(i,j,k)

bbyz1(i,j,k)=-2*sqrt(2.*(sxx1f(i,j,k)*sxx1f(i,j,k)+&
syy1f(i,j,k)*syy1f(i,j,k)+&
szz1f(i,j,k)*szz1f(i,j,k)+&
2.*sxy1f(i,j,k)*sxy1f(i,j,k)+&
2.*sxz1f(i,j,k)*sxz1f(i,j,k)+&
2.*syz1f(i,j,k)*syz1f(i,j,k)))&
*syz1f(i,j,k)
!

!Aij tensor with u
axx1(i,j,k)=-2*sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
syy1(i,j,k)*syy1(i,j,k)+&
szz1(i,j,k)*szz1(i,j,k)+&
2.*sxy1(i,j,k)*sxy1(i,j,k)+&
2.*sxz1(i,j,k)*sxz1(i,j,k)+&
2.*syz1(i,j,k)*syz1(i,j,k)))&
*sxx1(i,j,k)

ayy1(i,j,k)=-2*sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
syy1(i,j,k)*syy1(i,j,k)+&
szz1(i,j,k)*szz1(i,j,k)+&
2.*sxy1(i,j,k)*sxy1(i,j,k)+&
2.*sxz1(i,j,k)*sxz1(i,j,k)+&
2.*syz1(i,j,k)*syz1(i,j,k)))&
*syy1(i,j,k)

azz1(i,j,k)=-2*sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
syy1(i,j,k)*syy1(i,j,k)+&
szz1(i,j,k)*szz1(i,j,k)+&
2.*sxy1(i,j,k)*sxy1(i,j,k)+&
2.*sxz1(i,j,k)*sxz1(i,j,k)+&
2.*syz1(i,j,k)*syz1(i,j,k)))&
*szz1(i,j,k)

axy1(i,j,k)=-2*sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
syy1(i,j,k)*syy1(i,j,k)+&
szz1(i,j,k)*szz1(i,j,k)+&
2.*sxy1(i,j,k)*sxy1(i,j,k)+&
2.*sxz1(i,j,k)*sxz1(i,j,k)+&
2.*syz1(i,j,k)*syz1(i,j,k)))&
*sxy1(i,j,k)

axz1(i,j,k)=-2*sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
syy1(i,j,k)*syy1(i,j,k)+&
szz1(i,j,k)*szz1(i,j,k)+&
2.*sxy1(i,j,k)*sxy1(i,j,k)+&
2.*sxz1(i,j,k)*sxz1(i,j,k)+&
2.*syz1(i,j,k)*syz1(i,j,k)))&
*sxz1(i,j,k)

ayz1(i,j,k)=-2*sqrt(2.*(sxx1(i,j,k)*sxx1(i,j,k)+&
syy1(i,j,k)*syy1(i,j,k)+&
szz1(i,j,k)*szz1(i,j,k)+&
2.*sxy1(i,j,k)*sxy1(i,j,k)+&
2.*sxz1(i,j,k)*sxz1(i,j,k)+&
2.*syz1(i,j,k)*syz1(i,j,k)))&
*syz1(i,j,k)
!
enddo
enddo
enddo

call transpose_x_to_y(axx1,ta2)
call transpose_x_to_y(ayy1,tb2)
call transpose_x_to_y(azz1,tc2)
call transpose_x_to_y(axy1,td2)
call transpose_x_to_y(axz1,te2)
call transpose_x_to_y(ayz1,tf2)

do j=1,ysize(2)
ta2(:,j,:) = ta2(:,j,:)*(del(j))**2.
tb2(:,j,:) = tb2(:,j,:)*(del(j))**2.
tc2(:,j,:) = tc2(:,j,:)*(del(j))**2.
td2(:,j,:) = td2(:,j,:)*(del(j))**2.
te2(:,j,:) = te2(:,j,:)*(del(j))**2.
tf2(:,j,:) = tf2(:,j,:)*(del(j))**2.
enddo

call transpose_y_to_x(ta2,axx1)
call transpose_y_to_x(tb2,ayy1)
call transpose_y_to_x(tc2,azz1)
call transpose_y_to_x(td2,axy1)
call transpose_y_to_x(te2,axz1)
call transpose_y_to_x(tf2,ayz1)

call transpose_x_to_y(bbxx1,ta2)
call transpose_x_to_y(bbyy1,tb2)
call transpose_x_to_y(bbzz1,tc2)
call transpose_x_to_y(bbxy1,td2)
call transpose_x_to_y(bbxz1,te2)
call transpose_x_to_y(bbyz1,tf2)

! Multiply the test filter
do j=1,ysize(2)
ta2(:,j,:) = ta2(:,j,:)*(2.*del(j))**2.
tb2(:,j,:) = tb2(:,j,:)*(2.*del(j))**2.
tc2(:,j,:) = tc2(:,j,:)*(2.*del(j))**2.
td2(:,j,:) = td2(:,j,:)*(2.*del(j))**2.
te2(:,j,:) = te2(:,j,:)*(2.*del(j))**2.
tf2(:,j,:) = tf2(:,j,:)*(2.*del(j))**2.
enddo

call transpose_y_to_x(ta2,bbxx1)
call transpose_y_to_x(tb2,bbyy1)
call transpose_y_to_x(tc2,bbzz1)
call transpose_y_to_x(td2,bbxy1)
call transpose_y_to_x(te2,bbxz1)
call transpose_y_to_x(tf2,bbyz1)

!Need to filter Aij components
call filx(axx1f,axx1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(ayy1f,ayy1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(azz1f,azz1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(axy1f,axy1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(axz1f,axz1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)
call filx(ayz1f,ayz1,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)

call transpose_x_to_y(axx1f,ta2)
call transpose_x_to_y(ayy1f,tb2)
call transpose_x_to_y(azz1f,tc2)
call transpose_x_to_y(axy1f,td2)
call transpose_x_to_y(axz1f,te2)
call transpose_x_to_y(ayz1f,tf2)

call fily(axx2f,ta2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(ayy2f,tb2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(azz2f,tc2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(axy2f,td2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(axz2f,te2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)
call fily(ayz2f,tf2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)

ta2=0.;tb2=0.;tc2=0.
td2=0.;te2=0.;tf2=0.

call transpose_y_to_z(axx2f,ta3)
call transpose_y_to_z(ayy2f,tb3)
call transpose_y_to_z(azz2f,tc3)
call transpose_y_to_z(axy2f,td3)
call transpose_y_to_z(axz2f,te3)
call transpose_y_to_z(ayz2f,tf3)

call filz(axx3f,ta3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(ayy3f,tb3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(azz3f,tc3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(axy3f,td3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(axz3f,te3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
call filz(ayz3f,tf3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)

ta3=0.;tb3=0.;tc3=0.
td3=0.;te3=0.;tf3=0.

call transpose_z_to_y(axx3f,axx2f)
call transpose_z_to_y(ayy3f,ayy2f)
call transpose_z_to_y(azz3f,azz2f)
call transpose_z_to_y(axy3f,axy2f)
call transpose_z_to_y(axz3f,axz2f)
call transpose_z_to_y(ayz3f,ayz2f)

call transpose_y_to_x(axx2f,axx1f)
call transpose_y_to_x(ayy2f,ayy1f)
call transpose_y_to_x(azz2f,azz1f)
call transpose_y_to_x(axy2f,axy1f)
call transpose_y_to_x(axz2f,axz1f)
call transpose_y_to_x(ayz2f,ayz1f)

!Mij tensor OK
mxx1(:,:,:)=bbxx1(:,:,:)-axx1f(:,:,:)
myy1(:,:,:)=bbyy1(:,:,:)-ayy1f(:,:,:)
mzz1(:,:,:)=bbzz1(:,:,:)-azz1f(:,:,:)
mxy1(:,:,:)=bbxy1(:,:,:)-axy1f(:,:,:)
mxz1(:,:,:)=bbxz1(:,:,:)-axz1f(:,:,:)
myz1(:,:,:)=bbyz1(:,:,:)-ayz1f(:,:,:)
!

!Lij deviator
lxx1(:,:,:)=lxx1(:,:,:)-(lxx1(:,:,:)+lyy1(:,:,:)+lzz1(:,:,:))/3.
lyy1(:,:,:)=lyy1(:,:,:)-(lxx1(:,:,:)+lyy1(:,:,:)+lzz1(:,:,:))/3.
lzz1(:,:,:)=lzz1(:,:,:)-(lxx1(:,:,:)+lyy1(:,:,:)+lzz1(:,:,:))/3.

if(itime==1) then
dsmaggbl=0.01
else
smagC = 0.
do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)

Denom=(mxx1(i,j,k)*mxx1(i,j,k)+&
myy1(i,j,k)*myy1(i,j,k)+&
mzz1(i,j,k)*mzz1(i,j,k)+&
2.*mxy1(i,j,k)*mxy1(i,j,k)+&
2.*mxz1(i,j,k)*mxz1(i,j,k)+&
2.*myz1(i,j,k)*myz1(i,j,k))
if(abs(Denom).le.1e-6) then
    smagC(i,j,k)=0.01
else
smagC(i,j,k)= (lxx1(i,j,k)*mxx1(i,j,k)+&
lyy1(i,j,k)*myy1(i,j,k)+&
lzz1(i,j,k)*mzz1(i,j,k)+&
2.*lxy1(i,j,k)*mxy1(i,j,k)+&
2.*lxz1(i,j,k)*mxz1(i,j,k)+&
2.*lyz1(i,j,k)*myz1(i,j,k))/&
(mxx1(i,j,k)*mxx1(i,j,k)+&
myy1(i,j,k)*myy1(i,j,k)+&
mzz1(i,j,k)*mzz1(i,j,k)+&
2.*mxy1(i,j,k)*mxy1(i,j,k)+&
2.*mxz1(i,j,k)*mxz1(i,j,k)+&
2.*myz1(i,j,k)*myz1(i,j,k))
endif
!write(*,*) smagC(i,j,k)
!ERIC LIMITEUR SI BESOIN
!if(smagC(i,j,k).lt.0) then
!   smagC(i,j,k)=0.
!endif

enddo
enddo
enddo

!FILTERING THE NON-CONSTANT CONSTANT
call filx(smagC1f,smagC,di1,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,&
fiz1x,fiz2x,xsize(1),xsize(2),xsize(3),0)

call transpose_x_to_y(smagC1f,ta2)

call fily(smagC2f,ta2,di2,sy,vy,fiffy,fify,ficy,&
fiby,fibby,filay,fiz1y,fiz2y,ysize(1),ysize(2),ysize(3),1)

ta2=0.
call transpose_y_to_z(smagC2f,ta3)

call filz(smagC3f,ta3,di3,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,&
fiz1z,fiz2z,zsize(1),zsize(2),zsize(3),1)
ta3=0.

call transpose_z_to_y(smagC3f,smagC2f)
call transpose_y_to_x(smagC2f,smagC)

! Average coefficient within a horizontal plane
tmpa1=0.
do j=1,xsize(2)
do k=1,xsize(3)
do i=1,xsize(1)
tmpa1(j)=tmpa1(j)+smagC(i,j,k)
enddo
enddo
! DO the averaging
tmpa1(j)=tmpa1(j)/xsize(1)/xsize(3)
enddo

endif
call MPI_ALLREDUCE(tmpa1,dsmagHP1,xsize(2),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,code)

dsmagHP1(:) = dsmagHP1(:)/p_col

!if (nrank==0) print*,"Cst = ",maxval(dsmagcst),minval(dsmagcst)
!if (mod(itime,50)==0) print*,"Cst = ",maxval(dsmagcst),minval(dsmagcst)

call transpose_x_to_y(srt_smag,srt_smag2)
!!call transpose_x_to_y(dsmagcst,dsmagcst2)
do k=1,ysize(3)
do j=1,ysize(2)
do i=1,ysize(1)
nut2(i,j,k)=dsmagHP1(j)*(del(j)**2.0)*sqrt(2.*srt_smag2(i,j,k))
enddo
enddo
enddo
call transpose_y_to_x(nut2,nut1)
if (nrank==0) print*,"Max and Min of the  Constant = ",maxval(dsmagHP1), minval(dsmagHP1)
!!call test_sgs_min_max(dsmagcst,dsmagcst,dsmagcst,4)

return 

end subroutine dynsmag

!************************************************************
!
subroutine lesdiff(ux1,uy1,uz1,gxx1,gyy1,gzz1,gxy1,gxz1,gyz1,gyx1,gzx1,gzy1,nut1,&
    sgsx1,sgsy1,sgsz1,ep1,ta1,td1,te1,tf1,di1,ta2,td2,te2,tf2,tj2,di2,&
    ta3,td3,te3,tf3,di3)
!
!************************************************************
USE param
USE variables
USE decomp_2d


implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,nut1,ep1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,td1,te1,tf1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2,nut2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,td2,te2,tf2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3,nut3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,td3,te3,tf3,di3

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sgsx1,sgsy1,sgsz1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: sgsx2,sgsy2,sgsz2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: sgsx3,sgsy3,sgsz3

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxx1,gyx1,gzx1,&
gxy1,gyy1,gzy1,gxz1,gyz1,gzz1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gxy2,gyy2,gzy2,&
gxz2,gyz2,gzz2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: gxz3,gyz3,gzz3

integer :: i,j,k

ta1=0.;ta2=0.;ta3=0.;
sgsx1=0.;sgsy1=0.;sgsz1=0.
sgsx2=0.;sgsy2=0.;sgsz2=0.
sgsx3=0.;sgsy3=0.;sgsz3=0.
!WORK X-PENCILS
call derx (ta1,nut1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)

call derxx (td1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
call derxx (te1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
call derxx (tf1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)

sgsx1(:,:,:)=td1(:,:,:)*nut1(:,:,:)+gxx1(:,:,:)*ta1(:,:,:)
sgsy1(:,:,:)=te1(:,:,:)*nut1(:,:,:)+gyx1(:,:,:)*ta1(:,:,:)
sgsz1(:,:,:)=tf1(:,:,:)*nut1(:,:,:)+gzx1(:,:,:)*ta1(:,:,:)

!WORK Y-PENCILS
call transpose_x_to_y(sgsx1,sgsx2)
call transpose_x_to_y(sgsy1,sgsy2)
call transpose_x_to_y(sgsz1,sgsz2)
call transpose_x_to_y(gxy1,gxy2)
call transpose_x_to_y(gyy1,gyy2)
call transpose_x_to_y(gzy1,gzy2)
call transpose_x_to_y(gxz1,gxz2)
call transpose_x_to_y(gyz1,gyz2)
call transpose_x_to_y(gzz1,gzz2)
call transpose_x_to_y(nut1,nut2)
call transpose_x_to_y(ux1,ux2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)

call dery (ta2,nut2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
td2=0.
if (istret.ne.0) then
call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
call dery (te2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
do k=1,ysize(3)
do j=1,ysize(2)
do i=1,ysize(1)
td2(i,j,k)=td2(i,j,k)*pp2y(j)-pp4y(j)*te2(i,j,k)
enddo
enddo
enddo
else
call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
endif
te2=0.
!-->for uy
if (istret.ne.0) then
call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
call dery (tf2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
do k=1,ysize(3)
do j=1,ysize(2)
do i=1,ysize(1)
te2(i,j,k)=te2(i,j,k)*pp2y(j)-pp4y(j)*tf2(i,j,k)
enddo
enddo
enddo
else
call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
endif
tf2=0.
!-->for uz
if (istret.ne.0) then
call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
call dery (tj2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
do k=1,ysize(3)
do j=1,ysize(2)
do i=1,ysize(1)
tf2(i,j,k)=tf2(i,j,k)*pp2y(j)-pp4y(j)*tj2(i,j,k)
enddo
enddo
enddo
else
call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
endif

!      call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
!      call deryy (te2,uy2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
!      call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)

sgsx2(:,:,:)=sgsx2(:,:,:)+nut2(:,:,:)*td2(:,:,:)+gxy2(:,:,:)*ta2(:,:,:)
sgsy2(:,:,:)=sgsy2(:,:,:)+nut2(:,:,:)*te2(:,:,:)+gyy2(:,:,:)*ta2(:,:,:)
sgsz2(:,:,:)=sgsz2(:,:,:)+nut2(:,:,:)*tf2(:,:,:)+gzy2(:,:,:)*ta2(:,:,:)

!WORK Z-PENCILS
call transpose_y_to_z(sgsx2,sgsx3)
call transpose_y_to_z(sgsy2,sgsy3)
call transpose_y_to_z(sgsz2,sgsz3)
call transpose_y_to_z(gxz2,gxz3)
call transpose_y_to_z(gyz2,gyz3)
call transpose_y_to_z(gzz2,gzz3)
call transpose_y_to_z(nut2,nut3)
call transpose_y_to_z(ux2,ux3)
call transpose_y_to_z(uy2,uy3)
call transpose_y_to_z(uz2,uz3)

call derz(ta3,nut3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derzz (td3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
call derzz (te3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
call derzz (tf3,uz3,di3,sz,sfz,ssz,swz,zsize(1),zsize(2),zsize(3),0)

sgsx3(:,:,:)=sgsx3(:,:,:)+nut3(:,:,:)*td3(:,:,:)+gxz3(:,:,:)*ta3(:,:,:)
sgsy3(:,:,:)=sgsy3(:,:,:)+nut3(:,:,:)*te3(:,:,:)+gyz3(:,:,:)*ta3(:,:,:)
sgsz3(:,:,:)=sgsz3(:,:,:)+nut3(:,:,:)*tf3(:,:,:)+gzz3(:,:,:)*ta3(:,:,:)

call transpose_z_to_y(sgsx3,sgsx2)
call transpose_z_to_y(sgsy3,sgsy2)
call transpose_z_to_y(sgsz3,sgsz2)

call transpose_y_to_x(sgsx2,sgsx1)
call transpose_y_to_x(sgsy2,sgsy1)
call transpose_y_to_x(sgsz2,sgsz1)

!do k=1,xsize(3)
!do j=1,xsize(2)
!do i=1,xsize(1)
!if(ep1(i,j,k).eq.1) then
!sgsx1(i,j,k)=0.
!sgsy1(i,j,k)=0.
!sgsz1(i,j,k)=0.
!endif
!enddo
!enddo
!enddo

!call transpose_x_to_y(sgsx1,sgsx2)
!call transpose_x_to_y(sgsy1,sgsy2)
!call transpose_x_to_y(sgsz1,sgsz2)

!sgsx2(:,1,:) = 0.0; sgsx2(:,ysize(2),:) = 0.0
!sgsy2(:,1,:) = 0.0; sgsy2(:,ysize(2),:) = 0.0
!sgsz2(:,1,:) = 0.0; sgsz2(:,ysize(2),:) = 0.0

!call transpose_y_to_x(sgsx2,sgsx1)
!call transpose_y_to_x(sgsy2,sgsy1)
!call transpose_y_to_x(sgsz2,sgsz1)

!call test_sgs_min_max(sgsx1,sgsy1,sgsz1,2)

!call transpose_x_to_y(ux1,ux2)
!call transpose_x_to_y(uy1,uy2)
!call transpose_x_to_y(uz1,uz2)

!ux2(:,1,:)=1.;ux2(:,ysize(2),:)=1.
!uy2(:,1,:)=0.;uy2(:,ysize(2),:)=0.
!uz2(:,1,:)=0.;uz2(:,ysize(2),:)=0.


end subroutine lesdiff


