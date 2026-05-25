#include "definitions_inc.f90"
module f2_Parameters
  use mc_LogInfo
  use mc_TypeDef
  implicit none  
  private

  ! Log 
  type(LogType),public::MainLog
  
  ! Decomp2d options
  integer,public:: p_row,p_col
    
  ! Flow type option
  integer,parameter,public:: FT_AN  =0  ! Added new  
  integer,parameter,public:: FT_CH  =1  ! Channel
  integer,parameter,public:: FT_HC  =2  ! Open-Channel
  integer,parameter,public:: FT_Duct=3  ! Open-/Closed-Duct
  integer,parameter,public:: FT_TG  =4  ! Taylor-Green vortex
  integer,parameter,public:: FT_HIT =5  ! Homogenerous isotropic turbulence
  logical,public::IsUxConst=.false.
  real(RK),public::uBulk=0.0
  integer,public:: FlowType
  real(RK),dimension(4),public:: PrGradData=0.0_RK ! PrGradAve, PrGradAveOld, PrGradNow, ForcedOld
      
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
  real(RK),public:: SimTime=0.0_RK  ! Real simulation time
  integer,public :: iCFL       ! Use CFL condition to change time step dynamically( 1: yes, 2:no ).
  real(RK),public:: CFLc       ! CFL parameter
  integer,public::  itime      ! current time step
  integer,public::  ifirst     ! First iteration
  integer,public::  ilast      ! Last iteration 
  integer,public::  iadvance
  integer,public::  ischeme
  integer,public::  IsImplicit  !(0=full explicit, 1=partial implicit, 2=full implicit )
  real(RK),public:: pmGamma
  real(RK),public:: pmTheta
  real(RK),public:: pmAlpha
  real(RK),public:: pmBeta

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

  logical,public:: RestartFlag                       ! restart or not
  character(len=:),allocatable,public:: RunName_     ! Run name
  character(len=:),allocatable,public:: ResultsDir_  ! Result directory
  character(len=:),allocatable,public:: RestartDir_  ! Restart directory
  integer,public:: Cmd_LFile_Freq= 1                 ! report frequency in the terminal 
  integer,public:: LF_file_lvl   = 5                 ! logfile report level      
  integer,public:: LF_cmdw_lvl   = 3                 ! terminal report level
  logical, public :: delete_prev_restart_flag = .true. ! Delete previous restart or not

  ! limited velocity and div
  real(RK),public:: vel_limit
  real(RK),public:: div_limit

  public:: ReadAndInitParameters,DumpReadedParam,PMcoeUpdate
contains
    
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitParameters
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90 
  subroutine ReadAndInitParameters(chFile); implicit none 
    character(len=*),intent(in)::chFile
    
    ! locals
    integer:: iUnit,ierr, len_str
    character(512):: RunName, ResultsDir, RestartDir
    NAMELIST/BasicParam/FlowType,IsUxConst,uBulk,xlx,yly,zlz,nxc,nyc,nzc,xnu,dtMax,iCFL,CFLc,ifirst,    &
                        ilast,ischeme,IsImplicit,FFTW_plan_type,BcOption,gravity,uxBcValue,uyBcValue,   &
                        uzBcValue,ivstats,BackupFreq,SaveStat,SaveVisu,RestartFlag,Cmd_LFile_Freq,RunName,   &
                        ResultsDir,RestartDir,LF_file_lvl,LF_cmdw_lvl,p_row,p_col,vel_limit,div_limit,FluidDensity, &
                        delete_prev_restart_flag
 
    open(newunit=iUnit, file=chFile, status='old',form='formatted',IOSTAT=ierr )
    if(ierr/=0) then
      print*,"Cannot open file: "//strip(chFile); STOP
    endif
    read(iUnit, nml=BasicParam)
    close(iUnit,IOSTAT=ierr)
    RunName_    = strip(RunName)

    ResultsDir_ = strip(ResultsDir);  len_str = len(ResultsDir_)
    if(ResultsDir_(len_str : len_str) /= '/') ResultsDir_ = ResultsDir_ // '/'
    RestartDir_ = strip(RestartDir);  len_str = len(RestartDir_)
    if(RestartDir_(len_str : len_str) /= '/') RestartDir_ = RestartDir_ // '/'

    nxp= nxc+1
    nyp= nyc+1
    nzp= nzc+1
    if(ischeme==FI_AB2) then
      iadvance = 1
    elseif(ischeme==FI_RK2) then
      iadvance = 2
    elseif(ischeme==FI_RK3) then
      iadvance = 3
    else
      print*,"Time scheme WRONG!!! ischeme=",ischeme; STOP      
    endif   
  end subroutine ReadAndInitParameters

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! DumpReadedParam
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90 
  subroutine DumpReadedParam(); implicit none

    ! locals
    NAMELIST/BasicParam/FlowType,IsUxConst,uBulk,xlx,yly,zlz,nxc,nyc,nzc,xnu,dtMax,iCFL,CFLc,ifirst,    &
                        ilast,ischeme,IsImplicit,FFTW_plan_type,BcOption,gravity,uxBcValue,uyBcValue,   &
                        uzBcValue,ivstats,BackupFreq,SaveStat,SaveVisu,RestartFlag,Cmd_LFile_Freq,RunName_,  &
                        ResultsDir_,RestartDir_, LF_file_lvl,LF_cmdw_lvl,p_row,p_col,vel_limit,div_limit,FluidDensity, &
                        delete_prev_restart_flag
    write(MainLog%iUnit, nml=BasicParam)
  end subroutine DumpReadedParam

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PMcoeUpdate
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PMcoeUpdate(ns); implicit none
    integer,intent(in)::ns

    ! locals
    real(RK),dimension(3):: pmGammaConst,pmThetaConst,pmAlphaConst

    if(ischeme==FI_AB2) then
      pmGammaConst= [ 1.5_RK, 0.0_RK, 0.0_RK]
      pmThetaConst= [-0.5_RK, 0.0_RK, 0.0_RK]
      if((.not. RestartFlag) .and. itime==ifirst) then
        pmGammaConst= [ 1.0_RK, 0.0_RK, 0.0_RK]
        pmThetaConst= [ 0.0_RK, 0.0_RK, 0.0_RK]
      endif
    elseif(ischeme==FI_RK2) then
      pmGammaConst= [ 0.5_RK, 1.0_RK, 0.0_RK]
      pmThetaConst= [ 0.0_RK,-0.5_RK, 0.0_RK]
    else
      pmGammaConst= [8.0_RK/15.0_RK, 5.0_RK/12.0_RK, 0.75_RK]
      pmThetaConst= [0.0_RK, -17.0_RK/60.0_RK, -5.0_RK/12.0_RK]   
    endif
    pmAlphaConst = pmGammaConst + pmThetaConst
       
    pmGamma = pmGammaConst(ns) *dt
    pmTheta = pmThetaConst(ns) *dt
    pmAlpha = pmAlphaConst(ns) *dt
    pmBeta  = 0.5_RK*xnu *pmAlpha
    SimTime = SimTime +pmAlpha
    if(ns==1) then
      PrGradData(2)=PrGradData(1)
      PrGradData(1)=0.0_RK
    endif
  end subroutine PMcoeUpdate
end module f2_Parameters
