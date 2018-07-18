module actuator_disk_model

    ! Use the actuator_line Modules
    use decomp_2d, only: mytype, nrank
    use actuator_line_model_utils 
    use airfoils

    implicit none
    
    type ActuatorDiskType
        character(len=100):: name           ! Actuator disk name
        real(mytype) :: D                   ! Actuator disk diameter 
        real(mytype) :: COR(3)              ! Center of Rotation
        real(mytype) :: RotN(3)             ! axis of rotation
        real(mytype) :: CT                  ! Thrust coefficient
    end type ActuatorDiskType

    type(ActuatorDiskType), allocatable, save :: ActuatorDisk(:)
    integer,save :: Nad ! Number of the actuator disk turbines 

contains

    subroutine actuator_disk_model_init(Nactuatordisks)

        implicit none
        integer :: Nactuatordisks

    end subroutine actuator_disk_model_init

end module actuator_disk_model
