module actuator_line_beam_model

    ! Use the actuator_line Modules
    use decomp_2d, only: mytype, nrank
    use actuator_line_model_utils
    use airfoils, only : conrad 
    use actuator_line_element
    use xbeam_shared, only: xbelem, xbnode, xbopts
    use xbeam_undef, only: xbeam_undef_relatree
    use cbeam3_interface
    use lib_rotvect 
    
    type BeamType
    ! Global variables
    real(mytype), allocatable :: pos_ini(:,:)                       ! positions
    real(mytype), allocatable :: psi_ini(:,:,:)                     ! CRV of the nodes in the elements
    real(mytype),  allocatable :: pos_def(:,:)                       ! positions
    real(mytype),  allocatable :: psi_def(:,:,:)                     ! CRV of the nodes in the elements
    real(mytype),  allocatable :: pos_dot_def(:,:)                   ! positions
    real(mytype),  allocatable :: psi_dot_def(:,:,:)                 ! CRV of the nodes in the elements
    real(mytype),  allocatable :: static_forces(:,:)                 ! Static forces
    real(mytype),  allocatable :: dynamic_forces(:,:)                ! Dynamic forces
    real(mytype),  allocatable :: gravity_forces(:,:)                ! Gravity forces     
    real(mytype), dimension(4) :: quat
    integer :: Nnodes                                               ! Number of nodes
    integer :: Ndofs                                                ! Number of degrees of freedom=6*(num_nodes-3) (blades are clamped at the root)
    integer :: NElems                                               ! Number of elements
    type(xbelem), allocatable  :: elem(:)                           ! Element information.
    type(xbnode), allocatable  :: node(:)                           ! Nodal information.
    type(xbopts) :: options

    ! Beam Characteristics from File Read 
    real(mytype), allocatable :: rR(:),AeroCent(:), StrcTwist(:), BMassDen(:)  
    real(mytype), allocatable :: FlpStff(:),EdgStff(:), GJStff(:), EAStff(:)  
    real(mytype), allocatable :: Alpha(:), FlpInert(:), EdgInert(:), PrecrvRef(:), PreswpRef(:) 
    real(mytype), allocatable :: FlpcgOf(:),EdgcgOf(:),FlpEAOf(:),EdgEAOf(:)

    ! Beam working variables
    real(mytype), allocatable :: prebending(:),presweept(:),node_struct_twist(:),EA(:),EIy(:),EIz(:)
    real(mytype), allocatable :: GJ(:),GAy(:),GAz(:),mass(:),inertia_xb(:),inertia_yb(:),inertia_zb(:),pos_cg_B(:,:)
    real(mytype), allocatable :: blade_y_BFoR(:,:), mass_db(:,:,:),stiffness_matrix_db(:,:,:)
    real(mytype), allocatable :: Force(:,:)
    ! Constants
    real(mytype) :: poisson_coeff=0.3_mytype
   
    end type BeamType

contains

    subroutine actuator_line_beam_model_init(beam,acl,Nblades,FN)

        implicit none
    type(BeamType) :: beam
    integer, intent(in) :: Nblades
    type(ActuatorLineType),intent(in),dimension(3) :: acl(3)
    real(mytype) :: BFoRx(3), BFory(3),BForz(3), RotMat(3,3)
    integer :: iblade,jElem,knode, vcounter, fcounter
    character(len=100),intent(in)  :: FN ! FileName of the geometry file
    integer :: Nstations
    integer :: i,j, i_elem, i_node_local, master_elem, master_node, ielem
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
        beam%EIy(i)=0.5_mytype*(beam%EdgStff(i)+beam%EdgStff(i+1)) 
        beam%GJ(i)=0.5_mytype*(beam%GJStff(i)+beam%GJStff(i+1)) 
        beam%GAy(i)=beam%EA(i)/2.0_mytype*(1.0_mytype+beam%poisson_coeff)  
        beam%GAz(i)=beam%EA(i)/2.0_mytype*(1.0_mytype+beam%poisson_coeff)
        beam%mass(i)=0.5_mytype*(beam%BMassDen(i)+beam%BMassDen(i+1))
        beam%inertia_yb(i)=0.5_mytype*(beam%EdgInert(i)+beam%EdgInert(i+1))
        beam%inertia_zb(i)=0.5_mytype*(beam%FlpInert(i)+beam%FlpInert(i+1))
        beam%inertia_xb(i)=beam%inertia_yb(i)+beam%inertia_zb(i) 
        beam%pos_cg_B(i,2)=0.5_mytype*(beam%FlpcgOf(i)+beam%FlpcgOf(i+1))
        beam%pos_cg_B(i,3)=0.5_mytype*(beam%EdgcgOf(i)+beam%EdgcgOf(i+1))
    enddo
    ! Compute last value (only for node-arrays)
    beam%prebending(Nstations)=beam%PrecrvRef(Nstations)
    beam%presweept(Nstations)=beam%PreswpRef(Nstations)
    beam%node_struct_twist(Nstations)=beam%StrcTwist(Nstations)

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
    
    allocate(beam%Force(beam%Nnodes,6))

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
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem-1,1)=acl(iblade)%QCx(jelem)  -acl(iblade)%COR(1)       ! First-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem,  1)=acl(iblade)%PEx(jelem)  -acl(iblade)%COR(1)       ! Mid-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem+1,1)=acl(iblade)%QCx(jelem+1)-acl(iblade)%COR(1)     ! Last-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem-1,2)=acl(iblade)%QCy(jelem)  -acl(iblade)%COR(2)       ! First-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem,  2)=acl(iblade)%PEy(jelem)  -acl(iblade)%COR(2)       ! Mid-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem+1,2)=acl(iblade)%QCy(jelem+1)-acl(iblade)%COR(2)     ! Last-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem-1,3)=acl(iblade)%QCz(jelem)  -acl(iblade)%COR(3)       ! First-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem,  3)=acl(iblade)%PEz(jelem)  -acl(iblade)%COR(3)       ! Mid-point of the element
        beam%pos_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem+1,3)=acl(iblade)%QCz(jelem+1)-acl(iblade)%COR(3)     ! Last-point of the element
        
        ! BFoRx --> spanwise
        ! BFoRy --> tangential (negative)
        ! BFoRz --> normal (positive)
        BFoRx(1)= acl(iblade)%sEx(jelem)
        BFoRx(2)= acl(iblade)%sEy(jelem)
        BFoRx(3)= acl(iblade)%sEz(jelem)
        BFoRy(1)=-acl(iblade)%tEx(jelem)
        BFoRy(2)=-acl(iblade)%tEy(jelem)
        BFoRy(3)=-acl(iblade)%tEz(jelem)
        BFoRz(1)= acl(iblade)%nEx(jelem)
        BFoRz(2)= acl(iblade)%nEy(jelem)
        BFoRz(3)= acl(iblade)%nEz(jelem)

        RotMat=reshape((/BForx(1), BForx(2), BForx(3),&
                         BFory(1), BFory(2), BFory(3),&
                         BForz(1), BForz(2), BForz(3)/),(/3,3/))
         
        beam%psi_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem-1,1,:)=rotvect_mat2psi(RotMat)     ! First-point of the element
        beam%psi_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem,  1,:)=rotvect_mat2psi(RotMat)     ! Mid-point of the element
        beam%psi_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem+1,1,:)=rotvect_mat2psi(RotMat)     ! Last-point of the element
        beam%psi_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem-1,2,:)=rotvect_mat2psi(RotMat)     ! First-point of the element
        beam%psi_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem,  2,:)=rotvect_mat2psi(RotMat)     ! Mid-point of the element
        beam%psi_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem+1,2,:)=rotvect_mat2psi(RotMat)     ! Last-point of the element
        beam%psi_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem-1,3,:)=rotvect_mat2psi(RotMat)     ! First-point of the element
        beam%psi_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem,  3,:)=rotvect_mat2psi(RotMat)     ! Mid-point of the element
        beam%psi_ini((iblade-1)*(2*acl(iblade)%NElem+1)+2*jelem+1,3,:)=rotvect_mat2psi(RotMat)     ! Last-point of the element
        enddo
    
        beam%pos_def=beam%pos_ini
        beam%psi_def=beam%psi_ini
        beam%pos_dot_def=0.0_mytype  
        beam%psi_dot_def=0.0_mytype


        !Define Element information (beam%elem) based on the beam information
        do jelem=1,acl(iblade)%NElem    
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%NumNodes=3
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%MemNo=iblade
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%Conn=(/1,1,1/)*(2*(jelem-1)+2*(iblade-1)*acl(iblade)%NElem)+(/1,1,1/)*iblade+(/0,2,1/)
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%Length=acl(iblade)%EDS(jelem) ! Does not to be computed directly 
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%PreCurv=(/0.d0,0.d0,0.d0/) !  This is not used and therefore is set to
                                                                                     !  zero
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%Psi=(/0.d0,0.d0,0.d0/) ! Not used
            beam%elem(jelem+(iblade-1)*acl(iblade)%NElem)%Vector=(/0.d0,0.d0,0.d0/) ! Not used  
        enddo
    enddo 
        ! Define the master (tree of nodes shared between elements)
        call xbeam_undef_relatree(beam%elem)

        
        ! Define node (beam%node) params 
     do i_elem=1,beam%NElems
         do i_node_local=1,3
             if(beam%elem(i_elem)%master(i_node_local,1)==0) then
                  if(beam%node(beam%elem(i_elem)%Conn(i_node_local))%master(1)<1) then
                     beam%node(beam%elem(i_elem)%Conn(i_node_local))%master(1)=i_elem
                     beam%node(beam%elem(i_elem)%Conn(i_node_local))%master(2)=i_node_local
                 endif
             else
                 master_elem=beam%elem(i_elem)%master(i_node_local,1)
                 master_node=beam%elem(i_elem)%master(i_node_local,2)
                  if(beam%node(beam%elem(i_elem)%Conn(i_node_local))%master(1)<1) then
                     beam%node(beam%elem(i_elem)%Conn(i_node_local))%master(1)=master_elem
                     beam%node(beam%elem(i_elem)%Conn(i_node_local))%master(2)=master_node
                 endif
             endif
         enddo
     enddo
      

    do iblade=1,Nblades 
        fcounter=0
        vcounter=0
        do jelem=1,acl(iblade)%NElem
            if(jelem>1.and.jelem<acl(iblade)%NElem) then
            vcounter=vcounter+1
            fcounter=fcounter+1
            beam%node(jelem+(iblade-1)*acl(iblade)%NElem)%Vdof=vcounter
            beam%node(jelem+(iblade-1)*acl(iblade)%NElem)%Fdof=fcounter
            elseif(jelem==1) then
            fcounter=fcounter+1
            beam%node(jelem+(iblade-1)*acl(iblade)%NElem)%Fdof=fcounter
            elseif(jelem==acl(iblade)%NElem)  then
            vcounter=vcounter+1
            beam%node(jelem+(iblade-1)*acl(iblade)%NElem)%Vdof=vcounter
            endif
        enddo
    enddo
    beam%quat=(/1.0_mytype,0.0_mytype,0.0_mytype,0.0_mytype/)

    ! Set options 
    beam%options%deltaCurved=1d-2

    return
    end subroutine actuator_line_beam_model_init


    subroutine actuator_line_beam_solve(blade,beam,RotN,Nblades,angularVel,dt)

    use cbeam3_solv ! This loads the module from xbeam (WInc3D needs to be compiled against the xbeam library)

    implicit none
    
    type(BeamType) :: beam
    integer, intent(in) :: Nblades
    type(ActuatorLineType),intent(in) :: blade(3)
	real(mytype), intent(in) :: dt, angularVel, RotN(3)
    real(mytype) :: quat(4), for_vel(6), for_acc(6)
    real(mytype), allocatable :: CBA(:,:,:), q(:), dqdt(:) 
    integer :: inode, imaster_elem, imaster_node,i_index, iblade, jelem
    
    allocate(CBA(beam%Nnodes,3,3))
    allocate(q(beam%Ndofs+10),dqdt(beam%Ndofs+10))

    ! Translate the beam forces from the global to the local CRF
    beam%static_forces=0.0
    for_vel=(/0.0_mytype,0.0_mytype,0.0_mytype,RotN(1)*angularVel,RotN(2)*angularVel,RotN(3)*angularVel/)
    for_acc=0.0_mytype ! No acceleration
    ! Translate the dynamic forces
  
    do iblade=1,Nblades 
        do jelem=1,blade(iblade)%Nelem
            do inode=1,3
                i_index=(iblade-1)*(2*blade(iblade)%NElem+1)+2*jelem-1+inode-1
                beam%Force(i_index,1)=blade(iblade)%EFx(jelem)        
                beam%Force(i_index,2)=blade(iblade)%EFy(jelem)        
                beam%Force(i_index,3)=blade(iblade)%EFz(jelem)        
                beam%Force(i_index,4)=0.0_mytype       
                beam%Force(i_index,5)=0.0_mytype       
                beam%Force(i_index,6)=0.0_mytype        
            enddo
        enddo
    enddo
    
    do inode=1,beam%Nnodes
        CBA(inode, :, :)=rotvect_psi2Mat(beam%psi_def(beam%node(inode)%master(1),beam%node(inode)%master(2),:))
        beam%dynamic_forces(inode, 1:3) = MATMUL(CBA(inode, :, :), beam%FORCE(inode, 1:3)) ! FORCES IN A FoR
        beam%static_forces=0._mytype
    enddo


    !do iblade=1,Nblades 
    !    do jelem=1,blade(iblade)%Nelem
    !        do inode=1,3
    !            i_index=(iblade-1)*(2*blade(iblade)%NElem+1)+2*jelem-1+inode-1
    !        enddo
    !    enddo
    !enddo


    call cbeam3_solv_nlndyn_step(beam%Ndofs, &
                                 beam%NElems,&
                                 beam%Nnodes,&
                                 dt,&
                                 beam%elem,&
                                 beam%node,&
                                 beam%static_forces,&
                                 beam%dynamic_forces,&
                                 beam%gravity_forces,&
                                 beam%quat,&
                                 for_vel,&              ! Frame of reference velocity
                                 for_acc,&              ! Frame of reference acceleration
                                 beam%pos_ini,&         ! Initial -- unloaded position of the nodes
                                 beam%psi_ini,&         ! Initial -- unloaded cartertisian rotation vector
                                 beam%pos_def,&
                                 beam%psi_def,&
                                 beam%pos_dot_def,&
                                 beam%psi_dot_def,&
                                 q,&
                                 dqdt,&
                                 beam%options)

    return

    end subroutine actuator_line_beam_solve

end module actuator_line_beam_model
