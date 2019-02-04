module actuator_line_beam_model

    ! Use the actuator_line Modules
    use decomp_2d, only: mytype, nrank
    use airfoils 
    use actuator_line_element
    use xbeam_shared

    type BeamType
    ! Global variables
    real(mytype), allocatable :: pos_ini(:,:)                       ! positions
    real(mytype), allocatable :: psi_ini(:,:,:)                     ! CRV of the nodes in the elements
    real(mytype), allocatable :: pos_def(:,:)                       ! positions
    real(mytype), allocatable :: psi_def(:,:,:)                     ! CRV of the nodes in the elements
    real(mytype), allocatable :: pos_dot_def(:,:)                   ! positions
    real(mytype), allocatable :: psi_dot_def(:,:,:)                 ! CRV of the nodes in the elements
    real(mytype), allocatable :: static_forces(:,:)                 ! Static forces
    real(mytype), allocatable :: dynamic_forces(:,:)                ! Dynamic forces
    real(mytype), allocatable :: gravity_forces(:,:)                ! Gravity forces     
    integer :: Nnodes                                               ! Number of nodes
    integer :: Ndofs                                                ! Number of degrees of freedom=6*(num_nodes-3) (blades are clamped at the root)
    integer :: NElems                                               ! Number of elements
    type(xbelem), allocatable  :: elem(:)                           ! Element information.
    type(xbnode), allocatable  :: node(:)                           ! Nodal information.

    ! Beam Characteristics from File Read 
    real(mytype), allocatable :: rR(:),AeroCent(:), StrcTwist(:), BMassDen(:)  
    real(mytype), allocatable :: FlpStff(:),EdgStff(:), GJStff(:), EAStff(:)  
    real(mytype), allocatable :: Alpha(:), FlpInert(:), EdgInert(:), PrecrvRef(:), PreswpRef(:) 
    real(mytype), allocatable :: FlpcgOf(:),EdgcgOf(:),FlpEAOf(:),EdgEAOf(:)

    ! Beam working variables
    real(mytype), allocatable :: prebending(:),presweept(:),node_struct_twist(:),EA(:),EIy(:),EIz(:)
    real(mytype), allocatable :: GJ(:),GAy(:),GAz(:),mass(:),inertia_xb(:),inertia_yb(:),inertia_zb(:),pos_cg_B(:,:)
    real(mytype), allocatable :: blade_y_BFoR(:,:), mass_db(:,:,:),stiffness_matrix_db(:,:,:)
    ! Constants
    real(mytype) :: poisson_coeff=0.3_mytype

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
    
    ! OPEN and READ the structural characteristics of single blades 
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
   
    ! Make sure that the structural and aerodynamic blade elements are equal in number
    if(NStations/=acl(1)%Nelem+1) then
        stop
    endif
   
    ! Define some variables
    beam%StrcTwist(:)=-beam%StrcTwist(:)*conrad
    
    ! Allocate working arrays (for a single blade)
    allocate(beam%prebending(Nstations),beam%presweept(Nstations),beam%node_struct_twist(Nstations))
    allocate(beam%EA(Nstations-1), beam%EIy(Nstations-1),beam%EIz(Nstations-1),beam%GJ(Nstations-1))
    allocate(beam%GAy(Nstations-1),beam%GAz(Nstations-1),beam%mass(Nstations-1),beam%inertia_xb(Nstations-1),beam%inertia_yb(Nstations-1),beam%inertia_zb(Nstations-1))
    allocate(beam%pos_cg_B(Nstations-1,3),beam%blade_y_BFoR(Nstations-1,3),beam%mass_db(NStations-1,6,6),beam%stiffness_matrix_db(Nstations-1,6,6))

    ! Compute working arrays
        ! x    z
        ! ^   ^    Local co-ordinate system for the working arrays 
        ! |  /
        ! | / 
        ! |/  
        ! O--------> y
    do i=1,Nstations-1
        beam%prebending(i)=beam%PrecrvRef(i)
        beam%presweept(i)=beam%PreswpRef(i)
        beam%node_struct_twist(i)=beam%StrcTwist(i)
        beam%EA(i) =0.5_mytype*(beam%EAStff(i)+beam%EAStff(i+1))
        beam%EIz(i)=0.5_mytype*(beam%FlpStff(i)+beam%FlpStff(i+1)) ! Flap and edge-wise directions are in a local co-ordinate system
        beam%EIy(i)=0.5_mytype*(beam%FlpStff(i)+beam%FlpStff(i+1)) 
        beam%GJ(i)=0.5_mytype*(beam%GJStff(i)+beam%GJStff(i+1)) 
        beam%GAy(i)=beam%EA(i)/2.0_mytype*(1.0_mytype+beam%poisson_coeff)  
        beam%GAz(i)=beam%EA(i)/2.0_mytype*(1.0_mytype+beam%poisson_coeff)
        beam%mass(i)=0.5_mytype*(beam%BMassDen(i)+beam%BMassDen(i+1))
        beam%inertia_yb(i)=0.5_mytype*(beam%EdgInert(i)+beam%EdgInert(i+1))
        beam%inertia_zb(i)=0.5_mytype*(beam%FlpInert(i)+beam%FlpInert(i+1))
        beam%inertia_xb(i)=beam%inertia_yb(i)+beam%inertia_zb(i) 
        beam%pos_cg_B(i,2)=0.5_mytype*(beam%FlpcgOf(i)+beam%FlpcgOf(i+1))
        beam%pos_cg_B(i,3)=0.5_mytype*(beam%EdgcgOf(i)+beam%EdgcgOf(i+1))
        ! Frame of Reference Delta
        beam%blade_y_BFoR(i,:)=(/1.0_mytype,0.0_mytype,0.0_mytype/)
    enddo
    ! Compute last value (only for node-arrays)
    beam%prebending(Nstations)=beam%PrecrvRef(Nstations)
    beam%presweept(Nstations)=beam%PreswpRef(Nstations)
    beam%node_struct_twist(Nstations)=beam%StrcTwist(Nstations)
    beam%blade_y_BFoR(Nstations,:)=(/1.0_mytype,0.0_mytype,0.0_mytype/)

    ! Compute the element Mass matrix
    do i=1,Nstations-1
        beam%mass_db(i,1:3,1:3)=beam%mass(i)
        beam%mass_db(i,1:3,4:6)=reshape((/0.0_mytype,beam%mass(i)*beam%pos_cg_B(i,3),-1.0_mytype*beam%mass(i)*beam%pos_cg_B(i,2),&
                                        -1.0_mytype*beam%mass(i)*beam%pos_cg_B(i,3),0.0_mytype,beam%mass(i)*beam%pos_cg_B(i,1),&
                                        beam%mass(i)*beam%pos_cg_B(i,2),-1.0_mytype*beam%mass(i)*beam%pos_cg_B(i,1),0.0_mytype/),(/3,3/)) 
        beam%mass_db(i,4:6,1:3)=-1.0_mytype*beam%mass_db(i,1:3,4:6)  ! Mass anti-symmetry
        beam%mass_db(i,4:6,4:6)=reshape((/beam%inertia_xb(i),0.0_mytype,0.0_mytype,0.0_mytype,beam%inertia_yb(i),0.0_mytype,0.0_mytype,0.0_mytype,beam%inertia_zb(i)/),(/3,3/))
    enddo
   
    ! Compute the elemebt stiffness matrix
    do i=1,Nstations-1
        beam%stiffness_matrix_db(i,1:6,1:6)=reshape((/beam%EA(i),0.0_mytype,0.0_mytype,0.0_mytype,0.0_mytype,0.0_mytype,&
                                                     0.0_mytype,beam%GAy(i),0.0_mytype,0.0_mytype,0.0_mytype,0.0_mytype,&
                                                     0.0_mytype,0.0_mytype,beam%GAz(i),0.0_mytype,0.0_mytype,0.0_mytype,&
                                                     0.0_mytype, 0.0_mytype,0.0_mytype,beam%GJ(i),0.0_mytype,0.0_mytype,&
                                                     0.0_mytype, 0.0_mytype,0.0_mytype,0.0_mytype,beam%EIy(i),0.0_mytype,&
                                                     0.0_mytype, 0.0_mytype,0.0_mytype,0.0_mytype,0.0_mytype ,beam%EIz(i)/),(/6,6/))
    enddo
    
    ! Degrees of freedom for the multibeam. It should be equal to the number of blades times twice the number of elements plus one (midpoints + edges)
    beam%Nnodes=Nblades*(2*acl(1)%Nelem+1)
    beam%NElems=Nblades*acl(1)%Nelem
    beam%Ndofs=6*(beam%Nnodes-Nblades)

    ! First allocate the global variables
    allocate(beam%pos_ini(beam%Nnodes,3))
    allocate(beam%psi_ini(beam%Nnodes,3,3))           
    allocate(beam%pos_def(beam%Nnodes,3))           
    allocate(beam%psi_def(beam%Nnodes,3,3))           
    allocate(beam%pos_dot_def(beam%Nnodes,3))                   
    allocate(beam%psi_dot_def(beam%Nnodes,3,3))                 
    allocate(beam%static_forces(beam%Nnodes,6))                
    allocate(beam%dynamic_forces(beam%Nnodes,6))                
    allocate(beam%gravity_forces(beam%Nnodes,6))               
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

end module actuator_line_beam_model
