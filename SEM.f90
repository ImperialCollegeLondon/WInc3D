! Routine that makes a synthetic eddy method for simple inflow-outflow
! simulations, based on the work of Jarrin et al. 2006
module sem_module 

USE decomp_2d
USE param
    
    type eddytype
        real(mytype) :: eddy_len  ! Eddy length scale
        real(mytype) :: x_pos     ! Eddy's x position
        real(mytype) :: y_pos     ! Eddy's y position
        real(mytype) :: z_pos     ! Eddy's z position
        real(mytype) :: x_int     ! Eddy's z position
        real(mytype) :: y_int     ! Eddy's z position
        real(mytype) :: z_int     ! Eddy's z position
    end type eddytype

    type(eddytype), allocatable, dimension (:) :: eddy
    real(mytype), allocatable :: U(:,:), RS(:,:)
contains 

    subroutine initialise_sem

        implicit none
        integer :: i, Nprofile
        character(1000) :: ReadLine

        ! Allocate the number of eddies
        allocate(eddy(NEddies))

        
        ! Open file and READ the inlet flow statistics (MEAN VEL and REYNOLDS STRESSES)
        open(15,file=SEM_File)
        ! Read the Number of DATA 
        read(15,'(A)') ReadLine
        read(ReadLine(index(ReadLine,':')+1:),*) Nprofile 
    
        allocate(U(3,Nprofile),RS(6,Nprofile))
        ! Read the stations specs
        do i=1,Nprofile
    
        read(15,'(A)') ReadLine ! Profile

        read(ReadLine,*) U(1,i), U(2,i), U(3,i), RS(1,i), RS(2,i), RS(3,i), RS(4,i), RS(5,i), RS(6,i)

        end do
    
        close(15)

        return

    end subroutine initialise_sem

    subroutine add_synthetic_turbulence()
    
    implicit none
    integer :: i 
    call random_number(bxo)
    call random_number(byo)
    call random_number(bzo)
    

    ! Zero inlet
    bxx1=0; bxy1=0; bxz1=0;
   
    return
    end subroutine add_synthetic_turbulence
   
    subroutine fluctuation_generation
        
        implicit none
    
        real(mytype) :: time_sta, time_end, x0, y0, z0, f
        real(mytype) :: xmesh,ymesh,zmesh
        integer :: ieddy,j,k

    ! First generate fluctuations without combining mean dna rms data
    do k=1,xsize(3)
    zmesh=(k+xstart(3)-1-1)*dz
    do j=1,xsize(2)
    if (istret.eq.0) ymesh=(j+xstart(2)-1-1)*dy
    if (istret.ne.0) ymesh=yp(j+xstart(2)-1)
        do ieddy=1,NEddies
        x0=(0.-eddy(ieddy)%x_pos)/eddy(ieddy)%eddy_len
        y0=(ymesh-eddy(ieddy)%x_pos)/eddy(ieddy)%eddy_len
        z0=(zmesh-eddy(ieddy)%x_pos)/eddy(ieddy)%eddy_len
        enddo
    enddo
    enddo
    
    
    end subroutine fluctuation_generation

end module sem_module

