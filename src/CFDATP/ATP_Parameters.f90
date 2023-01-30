module ATP_Parameters
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d,only:nrank
  use m_Parameters,only: dtMax, ifirst, ilast, BackupFreq,xlx,yly,zlz,SaveVisu
  implicit none
  private
  
  ! Log
  type(LogType),public::ATPLogInfo
  
  integer,parameter,public:: x_axis = 1
  integer,parameter,public:: y_axis = 2
  integer,parameter,public:: z_axis = 3 
    
  integer,parameter,public:: PIM_AB2 = 2
  integer,parameter,public:: PIM_AB3 = 3
   
  ! default values  
  type ATP_Options
    logical:: RestartFlag=.false.
    integer:: numPrtcl    = 8000     ! total particle number
    integer:: np_InDomain            ! particle in domain
    integer:: ifirst                 ! first time step
    integer:: ilast                  ! last time step
    real(RK):: dt   =  1.0E-5_RK     ! time step 
    type(real3):: SimDomain_min
    type(real3):: SimDomain_max
    logical,dimension(3):: IsPeriodic = .false.

    integer:: PI_Method = PIM_AB2      ! integration scheme for translational motion   
    integer:: numPrtcl_Type=1          ! number of particle type 
    character(64)::RunName  = "ATPRun" ! run name
    character(64)::ResultsDir  = "."   ! result directory 
    character(64)::RestartDir="."      ! restart directory
    integer:: SaveVisu      = 1000     ! save frequency for visulizing file
    integer:: BackupFreq    = 100000   ! save frequency for restarting file
    integer:: Cmd_LFile_Freq= 500      ! report frequency in the terminal 
    integer:: LF_file_lvl   = 5        ! logfile report level      
    integer:: LF_cmdw_lvl   = 3        ! terminal report level
  contains 
    procedure :: ReadATPOption => LO_ReadATPOption
  end type ATP_Options
  type(ATP_Options),public::  ATP_opt
    
contains

  !**********************************************************************
  ! LO_ReadATPOption
  !**********************************************************************
  subroutine LO_ReadATPOption(this, chFile)
    implicit none
    class(ATP_Options):: this
    character(*),intent(in)::chFile
           
    ! locals
    logical::RestartFlag
    character(64):: RunName, ResultsDir,RestartDir
    integer::numPrtcl,SaveVisuATP,PI_Method,numPrtcl_Type,Cmd_LFile_Freq,LF_file_lvl, &
             LF_cmdw_lvl,nUnitFile, ierror
    NAMELIST/ATPOptions/ RestartFlag,numPrtcl,PI_Method,numPrtcl_Type,RunName,ResultsDir, &
                          RestartDir,Cmd_LFile_Freq,LF_file_lvl, LF_cmdw_lvl
               
    open(newunit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) then
       print*, "Cannot open file: "//trim(chFile); STOP
    endif
    read(nUnitFile, nml=ATPOptions)
    close(nUnitFile,IOSTAT=ierror)
    
    SaveVisuATP=SaveVisu
    this%RestartFlag = RestartFlag
    this%numPrtcl    = numPrtcl
    this%dt       = dtMax
    this%ifirst   = ifirst
    this%ilast    = ilast
    this%SimDomain_min = zero_r3
    this%SimDomain_max = real3(xlx,yly,zlz)
    this%IsPeriodic(1)=.true.
    this%IsPeriodic(2)=.false.
    this%IsPeriodic(3)=.true.
           
    this%PI_Method = PI_Method
    this%numPrtcl_Type = numPrtcl_Type
           
    write(this%RunName,"(A)") RunName
    write(this%ResultsDir,"(A)") ResultsDir
    write(this%RestartDir,"(A)") RestartDir
    this%SaveVisu = SaveVisuATP
    this%BackupFreq = BackupFreq
    this%Cmd_LFile_Freq = Cmd_LFile_Freq
    this%LF_file_lvl = LF_file_lvl
    this%LF_cmdw_lvl = LF_cmdw_lvl

  end subroutine LO_ReadATPOption
end module ATP_Parameters
