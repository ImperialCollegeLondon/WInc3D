!==================================================
! Subroutine for sampling the velocity from a point
!==================================================

! Routine that initializes the probe locations
subroutine init_probe

use param
!use variables
use var

    implicit none
    character(1000) :: ReadLine
    integer :: i

    open(70,file=Probelistfile)
    ! Read the Number of Blades
    read(70,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Nprobes

    ! Allocate the probe locations
    allocate(xprobe(Nprobes),yprobe(Nprobes),zprobe(Nprobes),uprobe(Nprobes),vprobe(Nprobes),wprobe(Nprobes))
    allocate(uprobe_part(Nprobes),vprobe_part(Nprobes),wprobe_part(Nprobes))

    do i=1,Nprobes

    read(70,'(A)') ReadLine ! Probe point ...

    read(ReadLine,*) xprobe(i), yprobe(i), zprobe(i)

    end do

end subroutine init_probe

subroutine init_probe_pencil

    USE param
    USE decomp_2d
    use var
    use constants
    implicit none

    NProbes=nx*MaxNPencils ! or xsize(1) it should be the same
    ! Allocate the probe locations
    allocate(xprobe(Nprobes),yprobe(Nprobes),zprobe(Nprobes),uprobe(Nprobes),vprobe(Nprobes),wprobe(Nprobes))
    allocate(uprobe_part(Nprobes),vprobe_part(Nprobes),wprobe_part(Nprobes))


end subroutine init_probe_pencil

subroutine probe_pencil(ux,uy,uz,phi)

    USE param
    USE decomp_2d
    use var
    use MPI

    implicit none

    real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,phi
    real(mytype) :: xmesh, ymesh, zmesh, dist, min_dist
    integer :: ipr, iy, ix, i, j, k, min_i, min_j, min_k, ierr
    real(mytype) :: ymin, ymax, zmin,zmax

    if (istret.eq.0) then
    ymin=(xstart(2)-1)*dy-dy/2.0 ! Add -dy/2.0 overlap
    ymax=(xend(2)-1)*dy+dy/2.0   ! Add +dy/2.0 overlap
    else
    ymin=yp(xstart(2))
    ymax=yp(xend(2))
    endif

    zmin=(xstart(3)-1)*dz-dz/2.0 ! Add a -dz/2.0 overlap
    zmax=(xend(3)-1)*dz+dz/2.0   ! Add a +dz/2.0 overlap

    do ipr=1,NProbes

        ix=mod(ipr,xsize(1))
        if (ix==0) ix=xsize(1)
        xprobe(ipr)=(ix-1)*dx

        iy=ipr/xsize(1)+1
        if (mod(ipr,xsize(1))==0) iy=ipr/xsize(1)
        yprobe(ipr)=y_loc_pencil(iy)*dy
        zprobe(ipr)=z_loc_pencil(iy)*dz

    if((yprobe(ipr)>=ymin).and.(yprobe(ipr)<=ymax).and.(zprobe(ipr)>=zmin).and.(zprobe(ipr)<=zmax)) then
        min_dist=1e6
        do k=1,xsize(3)
        zmesh=(k-1)*dz
        do j=1,xsize(2)


        if (istret.eq.0) ymesh=(j-1)*dy
        if (istret.ne.0) ymesh=yp(j)

        do i=1,xsize(1)
        xmesh=(i-1)*dx
        dist = sqrt((xprobe(ipr)-xmesh)**2.+(yprobe(ipr)-ymesh)**2.+(zprobe(ipr)-zmesh)**2.)

        if (dist<min_dist) then
            min_dist=dist
            min_i=i
            min_j=j
            min_k=k
        endif

        enddo
        enddo
        enddo
        uprobe_part(ipr)=ux(min_i,min_j,min_k)
        vprobe_part(ipr)=uy(min_i,min_j,min_k)
        wprobe_part(ipr)=uz(min_i,min_j,min_k)
    else
        uprobe_part(ipr)=0.0
        vprobe_part(ipr)=0.0
        wprobe_part(ipr)=0.0
        !write(*,*) 'Warning: I do not own this node'
    endif

    enddo

        call MPI_ALLREDUCE(uprobe_part,uprobe,Nprobes,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(vprobe_part,vprobe,Nprobes,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(wprobe_part,wprobe,Nprobes,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)

    return

end subroutine probe_pencil

subroutine probe(ux,uy,uz,phi)

use actuator_line_model_utils ! used only for the trilinear interpolation
USE param
USE decomp_2d
!USE variables
use var
use MPI

implicit none

real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,phi
integer :: i,j,k, ipr, ierr
real(mytype) :: ymin,ymax,zmin,zmax
real(mytype) :: xmesh,ymesh,zmesh
real(mytype) :: dist, min_dist
real(mytype) :: x0,y0,z0,x1,y1,z1,x,y,z,u000,u100,u001,u101,u010,u110,u011,u111
integer :: min_i,min_j,min_k
integer :: i_lower, j_lower, k_lower, i_upper, j_upper, k_upper

        if (istret.eq.0) then
        ymin=(xstart(2)-1)*dy-dy/2.0 ! Add -dy/2.0 overlap
        ymax=(xend(2)-1)*dy+dy/2.0   ! Add +dy/2.0 overlap
        else
        ymin=yp(xstart(2))
        ymax=yp(xend(2))
        endif

        zmin=(xstart(3)-1)*dz-dz/2.0 ! Add a -dz/2.0 overlap
        zmax=(xend(3)-1)*dz+dz/2.0   ! Add a +dz/2.0 overlap

        do ipr=1,Nprobes

        min_dist=1e6
        if((yprobe(ipr)>=ymin).and.(yprobe(ipr)<=ymax).and.(zprobe(ipr)>=zmin).and.(zprobe(ipr)<=zmax)) then
            !write(*,*) 'Warning: I own this node'
            do k=xstart(3),xend(3)
            zmesh=(k-1)*dz
            do j=xstart(2),xend(2)
            if (istret.eq.0) ymesh=(j-1)*dy
            if (istret.ne.0) ymesh=yp(j)
            do i=xstart(1),xend(1)
            xmesh=(i-1)*dx
            dist = sqrt((xprobe(ipr)-xmesh)**2.+(yprobe(ipr)-ymesh)**2.+(zprobe(ipr)-zmesh)**2.)

            if (dist<min_dist) then
                min_dist=dist
                min_i=i
                min_j=j
                min_k=k
            endif

            enddo
            enddo
            enddo

            if(yprobe(ipr)>ymax.or.yprobe(ipr)<ymin) then
            write(*,*) 'In processor ', nrank
            write(*,*) 'yprobe =', yprobe(ipr),'is not within the', ymin, ymax, 'limits'
            stop
            endif

            if(xprobe(ipr)>(min_i-1)*dx) then
                i_lower=min_i
                i_upper=min_i+1
            else if(xprobe(ipr)<(min_i-1)*dx) then
                i_lower=min_i-1
                i_upper=min_i
            else if(xprobe(ipr)==(min_i-1)*dx) then
                i_lower=min_i
                i_upper=min_i
            endif

            if (istret.eq.0) then
            if(yprobe(ipr)>(min_j-1)*dy.and.yprobe(ipr)<(xend(2)-1)*dy) then
                j_lower=min_j
                j_upper=min_j+1
            else if(yprobe(ipr)>(min_j-1)*dy.and.yprobe(ipr)>(xend(2)-1)*dy) then
                j_lower=min_j
                j_upper=min_j
            else if(yprobe(ipr)<(min_j-1)*dy.and.yprobe(ipr)>(xstart(2)-1)*dy) then
                j_lower=min_j-1
                j_upper=min_j
            else if(yprobe(ipr)<(min_j-1)*dy.and.yprobe(ipr)<(xstart(2)-1)*dy) then
                j_lower=min_j
                j_upper=min_j
            else if (yprobe(ipr)==(min_j-1)*dy) then
                j_lower=min_j
                j_upper=min_j
            endif
            else

            if(yprobe(ipr)>yp(min_j)) then
                j_lower=min_j
                j_upper=min_j+1
            else if(yprobe(ipr)<yp(min_j)) then
                j_lower=min_j-1
                j_upper=min_j
            else if (yprobe(ipr)==yp(min_j)) then
                j_lower=min_j
                j_upper=min_j
            endif
            endif

            if(zprobe(ipr)>(min_k-1)*dz.and.zprobe(ipr)<(xend(3)-1)*dz) then
                k_lower=min_k
                k_upper=min_k+1
            elseif(zprobe(ipr)>(min_k-1)*dz.and.zprobe(ipr)>(xend(3)-1)*dz) then
                k_lower=min_k
                k_upper=min_k
            else if(zprobe(ipr)<(min_k-1)*dz.and.zprobe(ipr)>(xstart(3)-1)*dz) then
                k_lower=min_k-1
                k_upper=min_k
            else if(zprobe(ipr)<(min_k-1)*dz.and.zprobe(ipr)<(xstart(3)-1)*dz) then
                k_lower=min_k
                k_upper=min_k
            else if (zprobe(ipr)==(min_k-1)*dz) then
                k_lower=min_k
                k_upper=min_k
            endif


            ! Prepare for interpolation
            x0=(i_lower-1)*dx
            x1=(i_upper-1)*dx
            if (istret.eq.0) then
            y0=(j_lower-1)*dy
            y1=(j_upper-1)*dy
            else
            y0=yp(j_lower)
            y1=yp(j_upper)
            endif
            z0=(k_lower-1)*dz
            z1=(k_upper-1)*dz

            x=xprobe(ipr)
            y=yprobe(ipr)
            z=zprobe(ipr)

            ! Apply interpolation kernels from 8 neighboring nodes
            uprobe_part(ipr)= trilinear_interpolation(x0,y0,z0, &
                                                  x1,y1,z1, &
                                                  x,y,z, &
                                                  ux1(i_lower,j_lower,k_lower), &
                                                  ux1(i_upper,j_lower,k_lower), &
                                                  ux1(i_lower,j_lower,k_upper), &
                                                  ux1(i_upper,j_lower,k_upper), &
                                                  ux1(i_lower,j_upper,k_lower), &
                                                  ux1(i_upper,j_upper,k_lower), &
                                                  ux1(i_lower,j_upper,k_upper), &
                                                  ux1(i_upper,j_upper,k_upper))

              vprobe_part(ipr)= trilinear_interpolation(x0,y0,z0, &
                                                  x1,y1,z1, &
                                                  x,y,z, &
                                                  uy1(i_lower,j_lower,k_lower), &
                                                  uy1(i_upper,j_lower,k_lower), &
                                                  uy1(i_lower,j_lower,k_upper), &
                                                  uy1(i_upper,j_lower,k_upper), &
                                                  uy1(i_lower,j_upper,k_lower), &
                                                  uy1(i_upper,j_upper,k_lower), &
                                                  uy1(i_lower,j_upper,k_upper), &
                                                  uy1(i_upper,j_upper,k_upper))

              wprobe_part(ipr)= trilinear_interpolation(x0,y0,z0, &
                                                  x1,y1,z1, &
                                                  x,y,z, &
                                                  uz1(i_lower,j_lower,k_lower), &
                                                  uz1(i_upper,j_lower,k_lower), &
                                                  uz1(i_lower,j_lower,k_upper), &
                                                  uz1(i_upper,j_lower,k_upper), &
                                                  uz1(i_lower,j_upper,k_lower), &
                                                  uz1(i_upper,j_upper,k_lower), &
                                                  uz1(i_lower,j_upper,k_upper), &
                                                  uz1(i_upper,j_upper,k_upper))

        else
            uprobe_part(ipr)=0.0
            vprobe_part(ipr)=0.0
            wprobe_part(ipr)=0.0
            !write(*,*) 'Warning: I do not own this node'
        endif
        enddo

        call MPI_ALLREDUCE(uprobe_part,uprobe,Nprobes,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(vprobe_part,vprobe,Nprobes,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(wprobe_part,wprobe,Nprobes,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)

end subroutine probe

subroutine write_probe(isample)
USE param
USE decomp_2d
!USE variables
use var

    implicit none
    integer,intent(in) :: isample
    integer :: ipr
    character(len=20) :: filename, Format1

990 format('probe',I6.6)
write(filename, 990) isample

    if (nrank==0) then
        open(2019,File=filename)
        write(2019,*) 'Probe ID, xprobe, yprobe, zprobe, u_probe, v_probe, z_probe'
        Format1="(I5,A,6(E14.7,A))"
        do ipr=1,Nprobes
        write(2019,Format1) ipr,',',xprobe(ipr), ',', yprobe(ipr),',', zprobe(ipr),',', uprobe(ipr),',', vprobe(ipr),',', wprobe(ipr)
        end do
        close(2019)
    endif
end subroutine
