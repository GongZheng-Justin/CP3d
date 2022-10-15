module m_Parameters
  use m_LogInfo
  use m_TypeDef
  implicit none  
  private

  ! Decomp2d options
  integer,public:: p_row,p_col
    
  ! Flow type option
  integer,parameter,public:: FT_CH=1  ! Channel
  integer,parameter,public:: FT_HC=2  ! Half Channel
  integer,parameter,public:: FT_TG=3  ! Taylor-Green vortex
  integer,parameter,public:: FT_HI=4  ! Homogenerous isotropic turbulence
  integer,parameter,public:: FT_AN=5  ! Added new
  integer, public:: FlowType
  real(RK),public:: PrGradAver=zero
  real(RK),public:: PrGradNow=zero
  
  ! mesh options
  real(RK),public::xlx, yly, zlz  ! domain length in three directions
  integer,public::nxp,nyp,nzp     ! grid point number in three directions
  integer,public::nxc,nyc,nzc     ! grid center number in three directions. nxc = nxp-1
  
  integer,parameter,public:: x_pencil=1
  integer,parameter,public:: y_pencil=2
  integer,parameter,public:: z_pencil=3
  integer,parameter,public:: xm_dir=1  ! x- direction
  integer,parameter,public:: xp_dir=2  ! x+ direction
  integer,parameter,public:: ym_dir=3  ! y- direction
  integer,parameter,public:: yp_dir=4  ! y+ direction
  integer,parameter,public:: zm_dir=5  ! z- direction
  integer,parameter,public:: zp_dir=6  ! z+ direction  
  
  ! Physical properties
  real(RK),public:: xnu                    ! Kinematic viscosity
  real(RK),public:: FluidDensity           ! Fluid density 
  real(RK),dimension(3),public:: gravity   ! Gravity or  other constant body forces (if any)
  
  ! Time stepping scheme and Projection method options
  integer,parameter,public:: FI_AB2 =1          !  AB2 for convective term
  integer,parameter,public:: FI_RK2 =2          !  RK2 for convective term
  integer,parameter,public:: FI_RK3 =3          !  RK3 for convective term
  real(RK),public:: dt         ! current time step
  real(RK),public:: dtMax      ! Maxium time step
  real(RK),public:: SimTime    ! Real simulation time
  integer,public :: iCFL       ! Use CFL condition to change time step dynamically( 1: yes, 2:no ).
  real(RK),public:: CFLc       ! CFL parameter
  integer,public::  itime      ! current time step
  integer,public::  ifirst     ! First iteration
  integer,public::  ilast      ! Last iteration 
  integer,public::  iadvance
  integer,public::  ischeme
  integer,public::  IsImplicit  !(0=full explicit, 1=partial implicit, 2=full implicit )
  real(RK),dimension(3),public :: pmGammaConst
  real(RK),dimension(3),public :: pmThetaConst
  real(RK),dimension(3),public :: pmAlphaConst
  real(RK),public::pmGamma
  real(RK),public::pmTheta
  real(RK),public::pmAlpha
  real(RK),public::pmAlphaC
  real(RK),public::pmBeta
  real(RK),public::pmBetaT

  ! FFT_option
  integer,public:: FFTW_plan_type
  
  ! Boundary conditions
  !  0: Periodic
  ! -1: NoSlip
  ! -2: Slip
  ! -3: Convective (ONLY AVAILABLE for x+ and y+)
  integer,parameter,public::BC_PERIOD  =  0 
  integer,parameter,public::BC_NoSlip  = -1
  integer,parameter,public::BC_FreeSlip= -2
  integer,parameter,public::BC_OutFlow = -3
  integer, public,dimension(6):: BcOption
  real(RK),public,dimension(6):: uxBcValue
  real(RK),public,dimension(6):: uyBcValue
  real(RK),public,dimension(6):: uzBcValue

  ! I/O, Statistics
  integer,public:: ivstats            ! time step interval for statistics calculation 
  integer,public:: SaveVisu           ! Output visulizing file frequency
  integer,public:: BackupFreq         ! Output Restarting file frequency
  integer,public:: SaveStat           ! Output Statistics file frequency

  logical,public::       RestartFlag  ! restart or not
  character(64),public:: RunName      ! Run name
  character(64),public:: ResultsDir   ! Result directory
  character(64),public:: RestartDir   ! Restart directory
  integer,public:: Cmd_LFile_Freq= 1  ! report frequency in the terminal 
  integer,public:: LF_file_lvl   = 5  ! logfile report level      
  integer,public:: LF_cmdw_lvl   = 3  ! terminal report level

  ! limited velocity and div
  real(RK),public:: vel_limit
  real(RK),public:: div_limit

  ! LES options
  ! 0:none; 1:Smagorinsky Model; 2:constant Smagorinsky Model; 3:Dynamic Smagorinsky Model; 4:MTS Model;
  integer,public:: LES_type
  integer,public:: FilterType   ! 0:trapezoidal type, 1:Simpson type

  public:: ReadAndInitParameters,DumpReadedParam,PMcoeUpdate
contains
    
  !******************************************************************
  ! InitParameters
  !****************************************************************** 
  subroutine ReadAndInitParameters(chFile)
    implicit none 
    character(*),intent(in)::chFile
    
    ! locals
    integer:: nUnitFile,ierror
    NAMELIST/BasicParam/FlowType,xlx,yly,zlz,nxc,nyc,nzc,xnu,dtMax,iCFL,CFLc,ifirst,ilast,ischeme, &
                        IsImplicit,FFTW_plan_type,BcOption,gravity,uxBcValue,uyBcValue,uzBcValue,  &
                        ivstats,BackupFreq,SaveStat,SaveVisu,RestartFlag,RunName,Cmd_LFile_Freq,   &
                        ResultsDir,RestartDir,LF_file_lvl,LF_cmdw_lvl,p_row,p_col,vel_limit,div_limit,FluidDensity
    NAMELIST/LesOptions/LES_type, FilterType
 
    open(newunit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0) then
      print*,"Cannot open file: "//trim(chFile); STOP
    endif
    read(nUnitFile, nml=BasicParam)
    rewind(nUnitFile)
    read(nUnitFile, nml=LesOptions)
    close(nUnitFile,IOSTAT=ierror)  

    nxp= nxc+1
    nyp= nyc+1
    nzp= nzc+1
    if(ischeme==FI_AB2) then
      iadvance = 1
      pmGammaConst = (/ three/two, zero, zero/)
      pmThetaConst = (/     -half, zero, zero/)
    elseif(ischeme==FI_RK2) then
      iadvance = 2
      pmGammaConst = (/ half,  one, zero/)
      pmThetaConst = (/ zero,-half, zero/)
    elseif(ischeme==FI_RK3) then
      iadvance = 3
      pmGammaConst = (/ eight/fifteen,       five/twelve,    three/four/)
      pmThetaConst = (/          zero,  -seventeen/sixty,  -five/twelve/)
    else
      print*,"Time scheme WRONG!!! ischeme=",ischeme; STOP      
    endif
    pmAlphaConst = pmGammaConst + pmThetaConst   
    
  end subroutine ReadAndInitParameters

  !******************************************************************
  ! DumpReadedParam
  !****************************************************************** 
  subroutine DumpReadedParam()
    implicit none

    ! locals
    NAMELIST/BasicParam/FlowType,xlx,yly,zlz,nxc,nyc,nzc,xnu,dtMax,iCFL,CFLc,ifirst,ilast,ischeme, &
                        IsImplicit,FFTW_plan_type,BcOption,gravity,uxBcValue,uyBcValue,uzBcValue,  &
                        ivstats,BackupFreq,SaveStat,SaveVisu,RestartFlag,RunName,Cmd_LFile_Freq,   &
                        ResultsDir,RestartDir,LF_file_lvl,LF_cmdw_lvl,p_row,p_col,vel_limit,div_limit,FluidDensity
    NAMELIST/LesOptions/LES_type, FilterType
    write(MainLog%nUnit, nml=BasicParam)
    write(MainLog%nUnit, nml=LesOptions)
  end subroutine DumpReadedParam

  !******************************************************************
  ! PMcoeUpdate
  !******************************************************************
  subroutine PMcoeUpdate(ns)
    implicit none
    integer,intent(in)::ns
      
    pmGamma = pmGammaConst(ns) * dt
    pmTheta = pmThetaConst(ns) * dt
    pmAlpha = pmAlphaConst(ns) * dt
    pmBeta  = half * pmAlphaConst(ns) *xnu * dt
    SimTime=SimTime+dt*real(pmAlphaConst(ns),RK)

    pmAlphaC= pmAlphaConst(ns)              ! for pressure gradient purpose only
    pmBetaT = half * pmAlphaConst(ns) * dt  ! for LES purpose only
    if(ns==1) PrGradAver=zero
  end subroutine PMcoeUpdate
  
end module m_Parameters
