module actuator_line_beam_model

    ! Use the actuator_line Modules
    use decomp_2d, only: mytype, nrank

    type BeamType
        character(len=100):: name                                ! Beam model name (the same as the turbine name)
        real(mytype), allocatable :: X(:)                        ! Element X coordinate (Local frame of reference) relative to the axis of rotation
        real(mytype), allocatable :: Y(:)                        ! Element Y coordinate (Local frame of reference) relative to the axis of rotation
        real(mytype), allocatable :: Z(:)                        ! Element Z coordinate (Local frame of reference) relative to the axis of rotation
        real(mytype), allocatable :: StructuralTwist(:)          ! Structural twist 
        real(mytype), allocatable :: fram_of_reference_delta(:)  ! Frame of reference delta 
        integer, allocatable :: Conn(:,:)                        ! Connectivity matrix 
    end type BeamType


contains

    subroutine actuator_line_beam_model_init()

        !=================================================
        ! Reading the beam information from .fem.h5 sharpy
        !=================================================
        implicit none
        
        
    end subroutine actuator_line_beam_model_init


    subroutine actuator_line_beam_solve()

        implicit none

    end subroutine actuator_line_beam_solve

end module actuator_line_beam_model
