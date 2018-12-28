module actuator_disc_model

    ! Use the actuator_line Modules
    use decomp_2d, only: mytype, nrank
    use actuator_line_model_utils 
    use airfoils

    implicit none
    
    type ActuatorDiscType
        integer :: ID                       ! Actuator disk ID
        real(mytype) :: D                   ! Actuator disk diameter 
        real(mytype) :: COR(3)              ! Center of Rotation
        real(mytype) :: RotN(3)             ! axis of rotation
        real(mytype) :: CT                  ! Thrust coefficient
        real(mytype) :: alpha               ! Induction coefficient
        real(mytype) :: Udisc               ! Disc-averaged velocity
        real(mytype) :: Udisc_prev          ! Disc-averaged velocity
        real(mytype) :: Power
        real(mytype) :: Udisc_ave=0.0_mytype
        real(mytype) :: Power_ave=0.0_mytype
    end type ActuatorDiscType

    type(ActuatorDiscType), allocatable, save :: ActuatorDisc(:)
    integer,save :: Nad ! Number of the actuator disk turbines 

contains

    subroutine actuator_disc_model_init(Ndiscs,admCoords,iadmmode,CT,aind,fileADM)
        
        USE var, only: Fdiscx, Fdiscy, Fdiscz, GammaDisc
        USE decomp_2d
        USE decomp_2d_io
        USE MPI

        implicit none
        integer,intent(in) :: Ndiscs
        character(len=100),intent(in) :: admCoords, fileADM
        character(1000) :: ReadLine
        integer,intent(in) :: iadmmode
        real(mytype) :: CT, aind
        integer :: idisc
        !### Specify the actuator discs
        Nad=Ndiscs
        if (nrank==0) then
        print *, "============================================================================================================="
        print *, " The actuator disc model is enabled"
        print *, 'Number of Actuator discs : ', Nad
        endif
        if (Nad>0) then 
            allocate(ActuatorDisc(Nad))
            open(15,file=admCoords)
            do idisc=1,Nad
            actuatordisc(idisc)%ID=idisc
            read(15,'(A)') ReadLine
            read(Readline,*) ActuatorDisc(idisc)%COR(1),ActuatorDisc(idisc)%COR(2),ActuatorDisc(idisc)%COR(3),ActuatorDisc(idisc)%RotN(1),ActuatorDisc(idisc)%RotN(2),ActuatorDisc(idisc)%RotN(3),ActuatorDisc(idisc)%D 
            !if (nrank==0) then
            !    print *, 'actuator ', idisc, ' --> (X,Y,Z)=  ',ActuatorDisc(idisc)%COR(1), ActuatorDisc(idisc)%COR(2),ActuatorDisc(idisc)%COR(3) 
            !    print *, '                          Axis =   ',ActuatorDisc(idisc)%RotN(1), ActuatorDisc(idisc)%RotN(2),ActuatorDisc(idisc)%RotN(3) 
            !    print *, '                          Diameter =   ',ActuatorDisc(idisc)%D 
            !endif
            if(iadmmode==0) then
                ActuatorDisc(idisc)%CT=CT
                ActuatorDisc(idisc)%alpha=aind
            else if(iadmmode==1) then
                if(nrank==0) print *, "not available yet"
                stop
            endif
            enddo
            close(15)
        endif
        !### Create the source term

        return
    end subroutine actuator_disc_model_init

    subroutine actuator_disc_model_compute_source(ux1,uy1,uz1)
        
        use decomp_2d, only: mytype, nproc, xstart, xend, xsize, update_halo
        use MPI
        use param, only: dx,dy,dz,eps_factor,xnu,yp,istret,xlx,yly,zlz,dt,u1,u2,iverifyadm,itime,ustar,dBL, spinup_time
        use var, only: FDiscx, FDiscy, FDiscz, GammaDisc
        
        implicit none
        real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1       
        real(mytype) :: xmesh, ymesh,zmesh,deltax,deltay,deltaz,deltar,dr,gamma_disc_partial
        real(mytype) :: uave,CTprime, T_relax, alpha_relax, Ratio
        real(mytype), allocatable, dimension(:) :: Udisc_partial
        integer,allocatable, dimension(:) :: counter, counter_total
        integer :: i,j,k, idisc, ierr

        ! First compute Gamma
        GammaDisc=0.
        do idisc=1,Nad
        do k=1,xsize(3)
        zmesh=(xstart(3)+k-1-1)*dz 
        do j=1,xsize(2)
        if (istret.eq.0) ymesh=(xstart(2)+j-1-2)*dy
        if (istret.ne.0) ymesh=yp(xstart(2)+j)
        do i=1,xsize(1)
        xmesh=(xstart(1)+i-1-1)*dx
        deltax=abs(xmesh-actuatordisc(idisc)%COR(1))
        deltay=abs(ymesh-actuatordisc(idisc)%COR(2))
        deltaz=abs(zmesh-actuatordisc(idisc)%COR(3))
        deltar=sqrt(deltay**2.+deltaz**2.)
        dr=sqrt(dx**2.+dz**2.)
        if(deltax>dx/2.) then 
            gamma_disc_partial=0.
        elseif (deltax<=dx/2.) then
            if(deltar<=actuatordisc(idisc)%D/2.) then
                gamma_disc_partial=1.
            elseif(deltar>actuatordisc(idisc)%D/2..and.deltar<=actuatordisc(idisc)%D/2.+2.*dr) then    
                gamma_disc_partial=1.-(deltar-actuatordisc(idisc)%D/2.)/(2*dr)    
            else
                gamma_disc_partial=0.
            endif
    
        endif
        GammaDisc(i,j,k)=GammaDisc(i,j,k)+gamma_disc_partial
        enddo
        enddo
        enddo 
        enddo
        
        ! Then compute disc-averaged velocity
        allocate(Udisc_partial(Nad))
        allocate(counter(Nad))
        allocate(counter_total(Nad))
     
        do idisc=1,Nad
        uave=0.
        counter(idisc)=0
        counter_total(idisc)=0.
        do k=1,xsize(3)
        zmesh=(xstart(3)+k-1-1)*dz 
        do j=1,xsize(2)
        if (istret.eq.0) ymesh=(xstart(2)+j-1-2)*dy
        if (istret.ne.0) ymesh=yp(xstart(2)+j)
        do i=1,xsize(1)
        xmesh=(i-1)*dx
        deltax=abs(xmesh-actuatordisc(idisc)%COR(1))
        deltay=abs(ymesh-actuatordisc(idisc)%COR(2))
        deltaz=abs(zmesh-actuatordisc(idisc)%COR(3))
        deltar=sqrt(deltay**2.+deltaz**2.)
        if(deltar<=actuatordisc(idisc)%D/2..and.deltax<dx/2.) then
            ! Take the inner product to compute the rotor-normal velocity
            uave=uave+ux1(i,j,k)*actuatordisc(idisc)%RotN(1)+&
                 uy1(i,j,k)*actuatordisc(idisc)%RotN(2)+&
                 uz1(i,j,k)*actuatordisc(idisc)%RotN(3)
                 counter(idisc)=counter(idisc)+1
        endif
        enddo
        enddo
        enddo
            Udisc_partial(idisc)=uave
        enddo

        call MPI_ALLREDUCE(Udisc_partial,actuatordisc%Udisc,Nad,MPI_REAL8,MPI_SUM, &
            MPI_COMM_WORLD,ierr)
        call MPI_ALLREDUCE(counter,counter_total,Nad,MPI_INTEGER,MPI_SUM, &
            MPI_COMM_WORLD,ierr)

        do idisc=1,Nad
        if (counter_total(idisc)==0) then
            print *, 'counter=0 for the disc --- something is wrong with the disc number: ', idisc
            stop
        endif
        actuatordisc(idisc)%Udisc=actuatordisc(idisc)%Udisc/counter_total(idisc)
        enddo
        
        deallocate(Udisc_partial,counter,counter_total)

        ! Time relaxation -- low pass filter
        if (itime==0) then
            do idisc=1,Nad
            actuatordisc(idisc)%Udisc_prev=actuatordisc(idisc)%Udisc
            enddo
        else
            do idisc=1,Nad
            T_relax=10 !0.27*dBL/ustar
            alpha_relax=(dt/T_relax)/(1.+dt/T_relax)
            actuatordisc(idisc)%Udisc=alpha_relax*actuatordisc(idisc)%Udisc+(1.-alpha_relax)*actuatordisc(idisc)%Udisc_prev
            actuatordisc(idisc)%Udisc_prev=actuatordisc(idisc)%Udisc
            enddo
        endif

        ! Compute the average values after spin_up time
        if(itime*dt>=spinup_time) then
            do idisc=1,Nad
             actuatordisc(idisc)%Udisc_ave=actuatordisc(idisc)%Udisc_ave+actuatordisc(idisc)%Udisc  
             actuatordisc(idisc)%Power_ave=actuatordisc(idisc)%Power_ave+actuatordisc(idisc)%Power   
            enddo
        endif
        ! Compute the forces
        do idisc=1,Nad

        CTprime=actuatordisc(idisc)%CT/(1-actuatordisc(idisc)%alpha)**2.
        ! Compute power
        actuatordisc(idisc)%Power=0.5*CTprime*actuatordisc(idisc)%Udisc**3.*pi*actuatordisc(idisc)%D**2./4.*0.432/0.56
        Fdiscx(:,:,:)=-0.5*CTprime*actuatordisc(idisc)%Udisc**2.*GammaDisc(:,:,:)/dx
        Fdiscy(:,:,:)=0.  
        Fdiscz(:,:,:)=0.
        ! Check if the total disc actuator disc F_t is equal to 
        if(iverifyadm.eq.1) then
            Ratio=sum(Fdiscx(:,:,:))/(-0.5*actuatordisc(idisc)%CT*((u1+u2)/2.0_mytype)**2.0_mytype*pi/4.0_mytype*actuatordisc(idisc)%D**2.0_mytype)
            if (nrank==0) print *, Ratio
        endif
        enddo

        return
    end subroutine actuator_disc_model_compute_source

 
    subroutine actuator_disc_model_write_output(dump_no)

        implicit none
        integer,intent(in) :: dump_no
        integer :: idisc
        character(len=100) :: dir, Format

        !call system('mkdir -p ADM/')
        
        if (Nad>0) then
        open(2020,File='discs_time'//trim(int2str(dump_no))//'.adm')
        write(2020,*) 'Udisc, CT, Power'
            Format="(3(E14.7,A))"
            do idisc=1,Nad
            write(2020,Format) actuatordisc(idisc)%Udisc,',',actuatordisc(idisc)%CT,',',actuatordisc(idisc)%Power,actuatordisc(idisc)%Udisc_ave,',',actuatordisc(idisc)%Power_ave
            end do
        close(2020)
        endif

        return 
    end subroutine actuator_disc_model_write_output
end module actuator_disc_model
