module acl_control
    
use decomp_2d, only: mytype, nrank
use actuator_line_turbine

contains 

subroutine Controller(turbine,time)
    
    ! This is bladed-style controller is used to implement a variable-speed
    ! generator-torque controller and PI collective blade pitch contoller for
    ! the NREL Offshore 5MW baseline wind turbine. 

    implicit none
    ! Passed variables:
    type(TurbineType), intent(inout) :: turbine
    real(mytype) :: time 
    
    ! Local variables
    real(mytype) :: Alpha                               ! Current coefficient n the recursive, signal
    real(mytype), dimension(3) :: BlPitch(3)            ! Current values of the blade pitch angles
    real(mytype) :: ElapTime                            ! Elapsed time since the last call to the con
    real(mytype), parameter :: CornerFreq = 1.570796    ! Corner frequency (-3 dB point) in the recursion
    real(mytype) :: GenSpeed                            ! Current HSS (generator) speed, rad/s
    real(mytype), save :: GenSpeedF                     ! Filtered HSS (generator) speed, rad/s
    real(mytype) :: GenTrq                              ! Electrical generator torque, N-m
    real(mytype) :: GK                                  ! Current value of the gain correction factor
    real(mytype) :: HorWindV                            ! Horizontal hub_height wind speed
    real(mytype), save :: IntSpdErr                     ! Current integral of speed error w.r.t time
    real(mytype), save :: LastGenTrq                    ! Commanded electrical generator torque the last time it was computed
    real(mytype), save :: LastTime                      ! LastTime the controller was used 
    real(mytype), save :: LastTimePC                    ! LastTime pitch controller was called
    real(mytype), save :: LastTimeVS                    ! LastTime torque controller was used

    ! Check validity of the input parameters

end subroutine Controller

end module acl_control
