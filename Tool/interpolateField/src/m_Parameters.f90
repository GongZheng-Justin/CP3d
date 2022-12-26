module m_Parameters
  use m_TypeDef
  use m_LogInfo
  implicit none  
  private

  ! Decomp2d options
  integer,public:: p_row,p_col
  
  integer,parameter,public:: xm_dir=1  ! x- direction
  integer,parameter,public:: xp_dir=2  ! x+ direction
  integer,parameter,public:: ym_dir=3  ! y- direction
  integer,parameter,public:: yp_dir=4  ! y+ direction
  integer,parameter,public:: zm_dir=5  ! z- direction
  integer,parameter,public:: zp_dir=6  ! z+ direction

  logical, public:: IsUxConst
  real(RK),public:: ubulk

  ! mesh options
  real(RK),public::yly
  integer,public::nxpOld,nypOld,nzpOld     ! Old grid point number in three directions
  integer,public::nxcOld,nycOld,nzcOld     ! Old grid center number in three directions. nxc = nxp-1
  integer,public::nxpNew,nypNew,nzpNew     ! New grid point number in three directions
  integer,public::nxcNew,nycNew,nzcNew     ! New grid center number in three directions. nxc = nxp-1
  
  ! Boundary conditions
  !  0 : periodic conditions for all variables
  ! -1: no slip for three velocity components and ZERO GRADIENT for pressure field
  ! -2: free slip for three velocity components and ZERO GRADIENT for pressure field
  integer,parameter,public::BC_PERIOD =  0  
  integer,parameter,public::BC_NSLIP  = -1
  integer,parameter,public::BC_FSLIP  = -2
  integer, public,dimension(6):: BcOption
  real(RK),public,dimension(6):: uxBcValue
  real(RK),public,dimension(6):: uyBcValue
  real(RK),public,dimension(6):: uzBcValue

  character(64),public:: RunName      ! Run name
  character(64),public:: ResultsDir   ! Result directory
  character(128),public:: OldRestartName
  character(128),public:: NewRestartName
  integer,public:: Cmd_LFile_Freq= 1  ! report frequency in the terminal 
  integer,public:: LF_file_lvl   = 5  ! logfile report level      
  integer,public:: LF_cmdw_lvl   = 3  ! terminal report level

#ifdef ScalarFlow
  !=== Scalar Options ===!
  ! Scalar boundary conditions
  ! -1: Dirichlet Bc. C     = F
  ! -2: Neumann Bc.   dC/dy = S
  ! -3: Robin Bc.     dC/dy + S*C = F
  ! y-,  y+
  logical, public:: IsScalarConst
  real(RK),public:: ScalarMean
  integer, public,dimension(2):: ScalarBcOption
  real(RK),public,dimension(4):: ScalarBcValues !(F_b,F_t,S_b,S_t)
#endif

  public:: ReadAndInitParameters,DumpReadedParam
contains
    
  !******************************************************************
  ! InitParameters
  !****************************************************************** 
  subroutine ReadAndInitParameters(chFile)
    implicit none 
    character(*),intent(in)::chFile
    
    ! locals
    integer:: nUnit,myistat
    NAMELIST/BasicParam/nxcOld,nycOld,nzcOld,nxcNew,nycNew,nzcNew,BcOption,uxBcValue,uyBcValue,uzBcValue,    &
                        RunName,ResultsDir,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,p_row,p_col,yly,IsUxConst, &
                        ubulk,OldRestartName,NewRestartName
#ifdef ScalarFlow
    NAMELIST/ScalarFlowOptions/IsScalarConst,ScalarMean,ScalarBcOption,ScalarBcValues
#endif

    open(newunit=nUnit, file=chFile, status='old',form='formatted',IOSTAT=myistat )
    if(myistat/=0) call MainLog%CheckForError(ErrT_Abort,"ReadAndInitParameters","Open file wrong: "//trim(chFile))
    read(nUnit,nml=BasicParam)
#ifdef ScalarFlow
    rewind(nUnit)
    read(nUnit, nml=ScalarFlowOptions)
#endif
    close(nUnit)  

    nxpOld =nxcOld+1; nypOld=nycOld+1; nzpOld=nzcOld+1
    nxpNew =nxcNew+1; nypNew=nycNew+1; nzpNew=nzcNew+1
  end subroutine ReadAndInitParameters

  !******************************************************************
  ! DumpReadedParam
  !****************************************************************** 
  subroutine DumpReadedParam()
    implicit none

    ! locals
    NAMELIST/BasicParam/nxpOld,nypOld,nzpOld,nxpNew,nypNew,nzpNew,BcOption,uxBcValue,uyBcValue,uzBcValue,    &
                        RunName,ResultsDir,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,p_row,p_col,yly,IsUxConst, &
                        ubulk,OldRestartName,NewRestartName
#ifdef ScalarFlow
    NAMELIST/ScalarFlowOptions/IsScalarConst,ScalarMean,ScalarBcOption,ScalarBcValues
#endif

    write(MainLog%nUnit, nml=BasicParam)
#ifdef ScalarFlow
    write(MainLog%nUnit, nml=ScalarFlowOptions)
#endif
  end subroutine DumpReadedParam
  
end module m_Parameters
