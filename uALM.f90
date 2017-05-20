program uALM

    use actuator_line_model

    implicit none

    integer :: ErrFlag, nargin, FNLength, status, DecInd
    logical :: back
    character(len=80) :: InputFN, FNBase
    character(len=80),dimension(100) :: TurbinesPath, ActuatorlinesPath
    character(len=20) :: filename
    integer :: ilit, isave, imodulo
    real(mytype) :: xlx, yly, zlz, re, sc, u1, u2, noise, noise1, dt, eps_factor, t
    integer :: itime, istart
    integer :: ialm, NTurbines, NActuatorlines
    
    NAMELIST/FlowParam/xlx,yly,zlz,re,sc,u1,u2,noise,noise1,dt
    NAMELIST/FileParam/ilit,isave,imodulo
    NAMELIST/ALMParam/ialm,NTurbines,TurbinesPath,NActuatorlines,ActuatorlinesPath,eps_factor

    !==========================================================================
    ! Handle Input file
    nargin=command_argument_count()
    if (nargin <1) then
        write(6,*) 'Please call the program with the name of the input file on the command line Ex. uALM input.in'
        stop
    endif
    
    call get_command_argument(1,InputFN,FNLength,status)
    back=.true.
    FNBase=inputFN((index(InputFN,'/',back)+1):len(InputFN))
    DecInd=index(FNBase,'.',back)
    if (DecInd >1) then
        FNBase=FNBase(1:(DecInd-1))
    end if
    !===========================================================================
    open(10,file=InputFN) 
    read(10,nml=FlowParam)
    read(10,nml=ALMParam)
    close(10) 

    ! Initialize the Actuator line model
    call actuator_line_model_init(Nturbines,Nactuatorlines,TurbinesPath,ActuatorlinesPath,dt)  
   
    Visc=1e-6

    ! Update
    do itime=1,10
        t=(itime-1)*dt
        call actuator_line_model_compute_forces
        call actuator_line_model_update(t,dt)
    end do

end program uALM
