module LPT_Parameters
  use LPT_TypeDef
  use m_Decomp2d,only:nrank
#ifdef CFDSecondOrder
  use m_Parameters,only: BcOption
#endif
  use m_Parameters,only: dtMax, ifirst, ilast, BackupFreq,xlx,yly,zlz,SaveVisu
  implicit none
  private
  
  integer,parameter,public:: x_axis = 1
  integer,parameter,public:: y_axis = 2
  integer,parameter,public:: z_axis = 3 
    
  integer,parameter,public:: PIM_AB2 = 2
  integer,parameter,public:: PIM_AB3 = 3
   
  ! default values  
  type LPT_Options
    logical:: RestartFlag=.false.
    integer:: numPrtcl    = 8000     ! total particle number
    integer:: np_InDomain            ! particle in domain
    integer:: ifirst                 ! first time step
    integer:: ilast                  ! last time step
    real(RK):: dt   =  1.0E-5_RK     ! time step 
    type(real3):: SimDomain_min
    type(real3):: SimDomain_max
    type(real3):: gravity = real3(zero,-9.81_RK,zero) ! gravity or other constant body forces if any
    logical,dimension(3):: IsPeriodic = .false.

    integer:: PI_Method = PIM_AB2      ! integration scheme for translational motion   
    integer:: numPrtcl_Type=1          ! number of particle type 
    character(64)::RunName  = "LPTRun" ! run name
    character(64)::ResultsDir  = "."   ! result directory 
    character(64)::RestartDir="."      ! restart directory
    integer:: SaveVisu      = 1000     ! save frequency for visulizing file
    integer:: BackupFreq    = 100000   ! save frequency for restarting file
    integer:: Cmd_LFile_Freq= 500      ! report frequency in the terminal 
    integer:: LF_file_lvl   = 5        ! logfile report level      
    integer:: LF_cmdw_lvl   = 3        ! terminal report level
  contains 
    procedure :: ReadLPTOption => LO_ReadLPTOption
  end type LPT_Options
  type(LPT_Options),public::  LPT_opt
    
contains

  !**********************************************************************
  ! LO_ReadLPTOption
  !**********************************************************************
  subroutine LO_ReadLPTOption(this, chFile)
    implicit none
    class(LPT_Options):: this
    character(*),intent(in)::chFile
           
    ! locals
    real(RK)::gravity(3)
    logical::RestartFlag
    character(64):: RunName, ResultsDir,RestartDir
    integer::numPrtcl,SaveVisuLPT,PI_Method,numPrtcl_Type,Cmd_LFile_Freq,LF_file_lvl, &
             LF_cmdw_lvl,nUnitFile, ierror
    NAMELIST /LPTOptions/ RestartFlag,numPrtcl,gravity,PI_Method,numPrtcl_Type,RunName,ResultsDir, &
                          RestartDir,Cmd_LFile_Freq,LF_file_lvl, LF_cmdw_lvl
               
    open(newunit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) then
       print*, "Cannot open file: "//trim(chFile); STOP
    endif
    read(nUnitFile, nml=LPTOptions)
    close(nUnitFile,IOSTAT=ierror)
    
    SaveVisuLPT=SaveVisu
    this%RestartFlag = RestartFlag
    this%numPrtcl    = numPrtcl
    this%dt       = dtMax
    this%ifirst   = ifirst
    this%ilast    = ilast
    this%gravity  = gravity
    this%SimDomain_min = zero_r3
    this%SimDomain_max = real3(xlx,yly,zlz)

#ifdef CFDSecondOrder
    if(BcOption(1)==0)this%IsPeriodic(1)=.true.
    if(BcOption(3)==0)this%IsPeriodic(2)=.true.
    if(BcOption(5)==0)this%IsPeriodic(3)=.true.
#elif CFDFourthOrder
    this%IsPeriodic(1)=.true.
    this%IsPeriodic(2)=.false.
    this%IsPeriodic(3)=.true.
#endif
           
    this%PI_Method = PI_Method
    this%numPrtcl_Type = numPrtcl_Type
           
    write(this%RunName,"(A)") RunName
    write(this%ResultsDir,"(A)") ResultsDir
    write(this%RestartDir,"(A)") RestartDir
    this%SaveVisu = SaveVisuLPT
    this%BackupFreq = BackupFreq
    this%Cmd_LFile_Freq = Cmd_LFile_Freq
    this%LF_file_lvl = LF_file_lvl
    this%LF_cmdw_lvl = LF_cmdw_lvl

  end subroutine LO_ReadLPTOption
end module LPT_Parameters
