subroutine shear_rate_coeff(ux1,uy1,uz1,gxx1,gyx1,gzx1,gxy1,gyy1,gzy1,gxz1,gyz1,gzz1,&
    sxx1,syy1,szz1,sxy1,sxz1,syz1,shrt2,shrt_coeff,ta2,ta3,di1,di2,di3)
!
!************************************************************
USE decomp_2d
USE param
USE MPI

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,di3

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: sxx1,syy1,szz1,&
sxy1,sxz1,syz1,nut1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gxx1,gyx1,gzx1,&
gxy1,gyy1,gzy1,gxz1,gyz1,gzz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: shrt2, shrt_coeff

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: syy2,szz2,&
sxy2,syz2,srt_smag2,nut2, shrt_coeff2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: gxy2,gyy2,gzy2,&
gxz2,gyz2,gzz2

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: szz3,syz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: gxz3,gyz3,gzz3

real(mytype) :: shearAve, shearAve_loc

integer::i,j,k,ierror,i2,j2,code

shrt2(:,:,:) = 0.
shrt_coeff(:,:,:) = 0.

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

shrt2(:,:,:) = sxx1(:,:,:)*sxx1(:,:,:)+syy1(:,:,:)*syy1(:,:,:) &
+szz1(:,:,:)*szz1(:,:,:)+2.*sxy1(:,:,:)*sxy1(:,:,:) &
+2.*sxz1(:,:,:)*sxz1(:,:,:)+2.*syz1(:,:,:)*syz1(:,:,:)

!================= Find the max |S|oo ================ !
shearAve_loc=maxval(sqrt(shrt2))

call MPI_ALLREDUCE(shearAve_loc,shearAve,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)

shrt_coeff=sqrt(shrt2)/shearAve

do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)
if (shrt_coeff(i,j,k)*rxxnu<10) shrt_coeff(i,j,k)=10./rxxnu ! Limits the value to 10 (just to be consistent)
enddo
enddo
enddo

if (nrank==0) print*, "Dynamic Hyper eddy viscosity with max and min of shrt_coeff=", maxval(shrt_coeff), minval(shrt_coeff)

return

end subroutine shear_rate_coeff


