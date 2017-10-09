!#########################################################
!
program fuz_sondes
!
!#########################################################
!!!  ifort -shared-intel -mcmodel=large -assume byterecl
implicit none

integer,parameter :: nt=8000,nx=48,nt_total=80000
real(8),dimension(nx) :: ux,uy,uz
real(8),dimension(nt_total,nx) :: ux1,uy1,uz1
real(4) :: dx,rmin,rmax,xl,dx1,rmin1,rmax1,xl1,dt
character(len=40) :: name0,name1,name2,name3,name4,name5,name6
integer :: ii,itime,i,j,k,l,nrank,nfile,iinter,ix,iy,icount
real(8),dimension(nt_total) :: p1,p2,p3

dt=0.00075*5.
nfile=10
nrank=1000


ux1=0.;uy1=0.;uz1=0.

do ii=1,nfile

990 format(I3.3,'/probes0341')
write(name0, 990) ii
open(nrank,file=name0,form='unformatted')
print *,'open DONE  ',name0
ux=0.;uy=0.;uz=0.
do itime=1,nt
read(nrank) ux,uy,uz
ux1(itime+nt*(ii-1),:)=ux(:)
uy1(itime+nt*(ii-1),:)=uy(:)
uz1(itime+nt*(ii-1),:)=uz(:)
enddo
print *,'read DONE',ii
close(nrank)

enddo

print *,'start writing file from different restart files'


open(10,file='FUZ_SONDES',form='unformatted')
write(10) ((ux1(i,j),i=1,nt_total),j=1,nx),((uy1(i,j),i=1,nt_total),j=1,nx),&
          ((uz1(i,j),i=1,nt_total),j=1,nx)
close(10)

open(10,file='test.dat',form='formatted')
do i=1,nt_total
   write(10,*) (i-1)*dt,ux1(i,24),uy1(i,36),uz1(i,12)
enddo
close(10)


open(10,file='test1.dat',form='formatted')
do i=1,nt_total
   write(10,*) (i-1)*dt,ux1(i,4),uy1(i,4),uz1(i,4)
enddo
close(10)

!open(96,file='moy.dat',form='formatted',status='unknown')
!do i=1,nt_total
!   p1(:)=0.
!   p2(:)=0.
!   p3(:)=0.
!   do j=1,i
!      p1(i)=p1(i)+ux1(j,24)
!      p2(i)=p2(i)+uy1(j,24)
!      p3(i)=p3(i)+uz1(j,24)
!   enddo
!   p1(i)=p1(i)/i
!   p2(i)=p2(i)/i
!   p3(i)=p3(i)/i
!   write(96,*) i*0.00075*5.,p1(i),p2(i),p3(i)
!enddo
!close(96)


end program fuz_sondes

