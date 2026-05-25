#include "definitions_inc.f90"
module ap_Parameters
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d,only:nrank
#ifdef CFDSecondOrder
  use f2_Parameters,only: BcOption
  use f2_Parameters,only: dtMax, ifirst, ilast, BackupFreq,xlx,yly,zlz,SaveVisu
#else
  use f4_Parameters,only: dtMax, ifirst, ilast, BackupFreq,xlx,yly,zlz,SaveVisu
#endif
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
    character(len=:),allocatable::RunName       ! run name
    character(len=:),allocatable::ResultsDir    ! result directory 
    character(len=:),allocatable::RestartDir    ! restart directory
    integer:: SaveVisu      = 1000     ! save frequency for Visulizing file
    integer:: BackupFreq    = 100000   ! save frequency for restarting file
    integer:: Cmd_LFile_Freq= 500      ! report frequency in the terminal 
    integer:: LF_file_lvl   = 5        ! logfile report level      
    integer:: LF_cmdw_lvl   = 3        ! terminal report level
  contains 
    procedure :: ReadATPOption => LO_ReadATPOption
  end type ATP_Options
  type(ATP_Options),public::  ATP_opt
    
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! LO_ReadATPOption
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine LO_ReadATPOption(this, chFile); implicit none
    class(ATP_Options):: this
    character(len=*),intent(in)::chFile
           
    ! locals
    logical::RestartFlag
    character(512):: RunName, ResultsDir,RestartDir
    integer::numPrtcl,SaveVisuATP,PI_Method,numPrtcl_Type,Cmd_LFile_Freq,LF_file_lvl, &
             LF_cmdw_lvl,iUnit, ierr,len_str
    NAMELIST/ATPOptions/ RestartFlag,numPrtcl,PI_Method,numPrtcl_Type,RunName,ResultsDir, &
                         RestartDir,Cmd_LFile_Freq,LF_file_lvl, LF_cmdw_lvl
               
    open(newunit=iUnit, file=chFile, status='old',form='formatted',IOSTAT=ierr)
    if(ierr/=0) then
       print*, "Cannot open file: "//strip(chFile); STOP
    endif
    read(iUnit, nml=ATPOptions)
    close(iUnit,IOSTAT=ierr)
    
    SaveVisuATP=SaveVisu
    this%RestartFlag = RestartFlag
    this%numPrtcl    = numPrtcl
    this%dt       = dtMax
    this%ifirst   = ifirst
    this%ilast    = ilast
    this%SimDomain_min = zero_r3
    this%SimDomain_max = real3(xlx,yly,zlz)

#ifdef CFDSecondOrder
    if(BcOption(1)==0) this%IsPeriodic(1)=.true.
    if(BcOption(3)==0) this%IsPeriodic(2)=.true.
    if(BcOption(5)==0) this%IsPeriodic(3)=.true.
#elif CFDFourthOrder
    this%IsPeriodic(1)=.true.
    this%IsPeriodic(2)=.false.
    this%IsPeriodic(3)=.true.
#endif
           
    this%PI_Method = PI_Method
    this%numPrtcl_Type = numPrtcl_Type
           
    this%RunName    = strip(RunName)

    this%ResultsDir = strip(ResultsDir);  len_str = len(this%ResultsDir)
    if(this%ResultsDir(len_str : len_str) /= '/') this%ResultsDir = this%ResultsDir // '/'
    this%RestartDir = strip(RestartDir);  len_str = len(this%RestartDir)
    if(this%RestartDir(len_str : len_str) /= '/') this%RestartDir = this%RestartDir // '/'

    this%SaveVisu = SaveVisuATP
    this%BackupFreq = BackupFreq
    this%Cmd_LFile_Freq = Cmd_LFile_Freq
    this%LF_file_lvl = LF_file_lvl
    this%LF_cmdw_lvl = LF_cmdw_lvl

  end subroutine LO_ReadATPOption
end module ap_Parameters
