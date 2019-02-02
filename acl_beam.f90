module actuator_line_beam_model

    ! Use the actuator_line Modules
    use decomp_2d, only: mytype, nrank
    use actuator_line_element
    use xbeam_shared

    type BeamType
    real(mytype), allocatable :: pos_ini(:,:)                       ! positions
    real(mytype), allocatable :: psi_ini(:,:,:)                     ! CRV of the nodes in the elements
    real(mytype), allocatable :: pos_def(:,:)                       ! positions
    real(mytype), allocatable :: psi_def(:,:,:)                     ! CRV of the nodes in the elements
    real(mytype), allocatable :: pos_dot_def(:,:)                   ! positions
    real(mytype), allocatable :: psi_dot_def(:,:,:)                 ! CRV of the nodes in the elements
    real(mytype), allocatable :: static_forces(:,:)                 ! Static forces
    real(mytype), allocatable :: dynamic_forces(:,:)                ! Dynamic forces
    real(mytype), allocatable :: gravity_forces(:,:)                ! Gravity forces
     
    real(mytype), allocatable :: StructuralTwist(:)
    real(mytype), allocatable :: frame_of_reference_delta(:)

    integer :: Nnodes                                               ! Number of nodes
    integer :: Ndofs                                                ! Number of degrees of freedom=6*(num_nodes-3) (blades are clamped at the root)
    integer :: NElems                                               ! Number of elements
    type(xbelem), allocatable  :: elem(:)                           ! Element information.
    type(xbnode), allocatable  :: node(:)                           ! Nodal information.

    ! Beam Characteristics 
    real(mytype), allocatable :: rR(:),AeroCent(:), StrcTwist(:), BMassDen(:)  
    real(mytype), allocatable :: FlpStff(:),EdgStff(:), GJStff(:), EAStff(:)  
    real(mytype), allocatable :: Alpha(:), FlpInert(:), EdgInert(:), PrecrvRef(:), PreswpRef(:) 
    real(mytype), allocatable :: FlpcgOf(:),EdgcgOf(:),FlpEAOf(:),EdgEAOf(:)

    end type BeamType

contains

    subroutine actuator_line_beam_model_init(beam,acl,Nblades,FN)

        implicit none
    type(BeamType) :: beam
    integer, intent(in) :: Nblades
    type(ActuatorLineType),intent(in),dimension(3) :: acl(3)
    integer :: iblade,jElem,knode
    character(len=100),intent(in)  :: FN ! FileName of the geometry file
    integer :: Nstations
    integer :: i
    character(1000) :: ReadLine
    
    ! OPEN and READ the structural characteristics of the blades
    
    open(15,file=FN)
    ! Read the Number of Nstations  
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Nstations

    allocate(beam%rR(Nstations),     beam%AeroCent(Nstations), beam%StrcTwist(Nstations), beam%BMassDen(Nstations))
    allocate(beam%FlpStff(Nstations),beam%EdgStff(Nstations),  beam%GJStff(Nstations),    beam%EAStff(Nstations))
    allocate(beam%Alpha(Nstations),  beam%FlpInert(Nstations), beam%EdgInert(Nstations),  beam%PrecrvRef(Nstations), beam%PreswpRef(Nstations))
    allocate(beam%FlpcgOf(Nstations),beam%EdgcgOf(Nstations),  beam%FlpEAOf(Nstations),   beam%EdgEAOf(Nstations))
    
    ! Read the station specs
    do i=1,NStations
    
    read(15,'(A)') ReadLine ! Structure ....

    read(ReadLine,*) beam%rR(i),    beam%AeroCent(i), beam%StrcTwist(i), beam%BMassDen(i),  beam%FlpStff(i),   beam%EdgStff(i), beam%GJStff(i),  beam%EAStff(i), & 
                     beam%Alpha(i), beam%FlpInert(i), beam%EdgInert(i),  beam%PrecrvRef(i), beam%PreswpRef(i), beam%FlpcgOf(i), beam%EdgcgOf(i), beam%FlpEAOf(i), beam%EdgEAOf(i)   

    end do
    
    close(15)


    ! Degrees of freedom for the beam. It should be equal to the number of blades times twice the number of elements plus one (midpoints + edges)
    beam%Nnodes=Nblades*(2*acl(1)%Nelem+1)
    beam%NElems=Nblades*acl(1)%Nelem
    beam%Ndofs=6*(beam%Nnodes-Nblades)

    ! First allocate
    allocate(beam%pos_ini(beam%Nnodes,3))
    allocate(beam%elem(beam%NElems))
    allocate(beam%node(beam%Nnodes))
    
    do iblade=1,Nblades
    
        ! Init the coordinates from the actuator line model
        do jelem=1,acl(iblade)%Nelem
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem-1,1)=acl(iblade)%QCx(jelem)  ! First-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem,1)=acl(iblade)%PEx(jelem)      ! Mid-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem+1,1)=acl(iblade)%QCx(jelem+1) ! Last-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem-1,2)=acl(iblade)%QCy(jelem)  ! First-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem,2)=acl(iblade)%PEy(jelem)      ! Mid-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem+1,2)=acl(iblade)%QCy(jelem+1) ! Last-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem-1,3)=acl(iblade)%QCz(jelem)  ! First-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem,3)=acl(iblade)%PEz(jelem)      ! Mid-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem+1,3)=acl(iblade)%QCz(jelem+1) ! Last-point of the element
        enddo
    
        !Define Element information (beam%elem) based on the beam information
        do jelem=1,acl(iblade)%NElem    
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%NumNodes=3
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%MemNo=iblade
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%Conn=(/1.0d0,1.0d0,1.0d0/)*((jelem-1)*(3-1)+(iblade-1))+(/0.0d0,2.0d0,1.0d0/)
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%Master(3,1)=0 ! DONT KNOW THIS ONE
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%Length=acl(iblade)%EDS(jelem) ! DONT KNOW THIS ONE
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%PreCurv=(/0.d0,0.d0,0.d0/) ! DONT KNOW THIS ONE
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%Psi=(/0.d0,0.d0,0.d0/) ! DONT KNOW THIS ONE
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%Vector=(/0.d0,0.d0,0.d0/) ! DONT KNOW THIS ONE  
        enddo
    enddo 
    
     
    return
    end subroutine actuator_line_beam_model_init


    subroutine actuator_line_beam_solve(beam,dt)

    use cbeam3_solv ! This loads the module from xbeam (WInc3D needs to be compiled against the xbeam library)

        implicit none
    type(BeamType),intent(inout) :: beam
	real(mytype),intent(in) :: dt

    ! Update the position of the beams
    
    ! Update the velocities


    ! Call the cbeam solver (timestep variation -- not sure if that is the correct thing to do)
    !call cbeam3_solv_nlndyn_step(beam%Ndofs, &
        !                            beam%NElems,&
        !                            beam%Nnodes,&
        !                            dt,&
        !                            elem,&
        !                            node,&
        !                            static_forces,&
        !                            dynamic_forces,&
        !                            gravity_forces,&
        !                            quat,&
        !                            for_vel,&
        !                            for_acc,&
        !                            pos_ini,&
        !                            psi_ini,&
        !                            pos_def,&
        !                            psi_def,&
        !                            pos_dot_def,&
        !                            psi_dot_def,&
        !                            options)

    return

    end subroutine actuator_line_beam_solve

    subroutine  read_beam_characteristics(FN,rR,AeroCent,StrcTwist,BMassDen,FlpStff,EdgStff,GJStff,EAStff,Alpha,&
                       FlpInert,EdgInert,PrecrvRef,PreswpRef,FlpcgOf,EdgcgOf,FlpEAOf,EdgEAOf,Nstations)
    
    implicit none
    character(len=100),intent(in)  :: FN ! FileName of the geometry file
    real(mytype), allocatable,intent(out) :: rR(:),AeroCent(:), StrcTwist(:), BMassDen(:)  
    real(mytype), allocatable,intent(out) :: FlpStff(:),EdgStff(:), GJStff(:), EAStff(:)  
    real(mytype), allocatable,intent(out) :: Alpha(:), FlpInert(:), EdgInert(:), PrecrvRef(:), PreswpRef(:) 
    real(mytype), allocatable,intent(out) :: FlpcgOf(:),EdgcgOf(:),FlpEAOf(:),EdgEAOf(:)
    integer, intent(out) :: Nstations
    integer :: i
    character(1000) :: ReadLine
    
    open(15,file=FN)
    ! Read the Number of Nstations  
    read(15,'(A)') ReadLine
    read(ReadLine(index(ReadLine,':')+1:),*) Nstations

    allocate(rR(Nstations),AeroCent(Nstations),StrcTwist(Nstations),BMassDen(Nstations))
    allocate(FlpStff(Nstations),EdgStff(Nstations),GJStff(Nstations),EAStff(Nstations))
    allocate(Alpha(Nstations),FlpInert(Nstations),EdgInert(Nstations),PrecrvRef(Nstations),PreswpRef(Nstations))
    allocate(FlpcgOf(Nstations),EdgcgOf(Nstations),FlpEAOf(Nstations),EdgEAOf(Nstations))
    
    ! Read the station specs
    do i=1,NStations
    
    read(15,'(A)') ReadLine ! Structure ....

    read(ReadLine,*) rR(i), AeroCent(i), StrcTwist(i), BMassDen(i), FlpStff(i), EdgStff(i), GJStff(i), EAStff(i), & 
                     Alpha(i), FlpInert(i), EdgInert(i), PrecrvRef(i), PreswpRef(i), FlpcgOf(i), EdgcgOf(i), FlpEAOf(i), EdgEAOf(i)   

    end do
    
    close(15)


    end subroutine read_beam_characteristics 




end module actuator_line_beam_model
