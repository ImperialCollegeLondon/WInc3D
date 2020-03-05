module constants

    use decomp_2d, only: mytype

    ! Mathematical constants
    real(mytype), parameter :: pi = 3.14159265358979323846_mytype
    real(mytype), parameter :: conrad = pi / 180.0
    real(mytype), parameter :: condeg = 180.0 / pi

    ! Definition of maximum size of arrays
    integer, parameter :: MaxNPencils = 90 ! Maximum number of pencils (in each XY direction) to be read from the input file

    integer, parameter :: MaxNAirfoils = 30 ! Maximum number of airfoils to be read
    integer, parameter :: MaxReVals = 10 ! Maximum number of tables (one per Reynolds number) that will be read
    integer, parameter :: MaxAOAVals = 1000 ! Maximum number of angles of attack in each polar that will be read

    integer, parameter :: MaxReadLine = 1000 ! I think this is never used

end module constants
