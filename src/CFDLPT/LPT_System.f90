module LPT_System
  use MPI
  use m_Timer
  use LPT_Comm
  use m_TypeDef
  use LPT_Property
  use LPT_Geometry
  use LPT_decomp_2d
  use LPT_Variables
  use LPT_IOAndVisu
  use LPT_Parameters
  use LPT_Statistics
  use LPT_Integration
  use LPT_ContactSearchPW
  use m_Decomp2d,only:nrank
  use m_Parameters,only:ivstats
  implicit none
  private
    
  !// LPTSystem class 
  type LPTSystem
    integer :: iterNumber   = 0  ! iteration number 
        
    !// timers
    type(timer):: m_total_timer
    type(timer):: m_integration_timer
    type(timer):: m_write_prtcl_timer
    type(timer):: m_comm_exchange_timer
  contains
    procedure:: Initialize => LPT_Initialize
    procedure:: iterate    => LPT_iterate
  end type LPTSystem
  type(LPTSystem),public::LPT
  
  integer::iCountLPT
contains

!********************************************************************************
!   Initializing LPTSystem object with particles which are inserted from a 
!   predefined plane 
!********************************************************************************
  subroutine  LPT_Initialize(this,chLPTPrm)
    implicit none
    class(LPTSystem)::this
    character(*),intent(in)::chLPTPrm
    
    ! locals
    integer::ierror
    character(256)::chStr
    real(RK)::t_restart1,t_restart2,t_res_tot
    
    !// Initializing main log info and visu
    iCountLPT=0
    if(LPT_Opt%RestartFlag) iCountLPT=10
    this%IterNumber=LPT_Opt%ifirst-1
    write(chStr,"(A)") 'mkdir -p '//LPT_opt%ResultsDir//' '//LPT_opt%RestartDir//' 2> /dev/null'
    if (nrank==0) call system(trim(adjustl(chStr)))
    call LPTLogInfo%InitLog(LPT_opt%ResultsDir,LPT_opt%RunName,LPT_opt%LF_file_lvl,LPT_opt%LF_cmdw_lvl)
    if(nrank==0) call LPTLogInfo%CreateFile(LPT_opt%RunName)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call LPTLogInfo%OpenFile()
    if(nrank==0)call Write_LPT_Opt_to_Log()
    call InitLPTStatistics(chLPTPrm)

    ! Step1: Physical property
    call LPTProperty%InitPrtclProperty(chLPTPrm)
    if(nrank==0) then
      call LPTLogInfo%OutInfo("Step1: Physical properties of particels and walls are set.",1)
      call LPTLogInfo%OutInfo("Physical properties contains "// trim( num2str(LPT_opt%numPrtcl_Type ) ) //" particle types ",2)
    endif

    ! Step2: set the geometry
    call LPTGeometry%MakeGeometry()
    if(nrank==0) then
      call LPTLogInfo%OutInfo("Step2: Geometry is set", 1 )
      call LPTLogInfo%OutInfo("Geometry Contains "//trim(num2str(LPTGeometry%num_pWall))//" Plane walls.", 2)
    endif

    ! Step3: initilize all the particle variables
    call GPrtcl_list%AllocateAllVar()
    call LPT_IO%Init_visu(chLPTPrm,1)
    t_restart1=MPI_WTIME()
    if(.not.LPT_Opt%RestartFlag) then
      if(LPT_Opt%numPrtcl>0) call GPrtcl_list%MakingAllPrtcl(chLPTPrm)
      if(nrank==0) then
        call LPTLogInfo%OutInfo("Step3: Initial Particle coordinates are MAKING into LPTSystem ...", 1 )
        call LPTLogInfo%OutInfo("Number of particles avaiable in the system:"//trim(num2str(LPT_opt%numPrtcl)),2)
      endif
      LPT_opt%np_InDomain = LPT_opt%numPrtcl
    else
      if(LPT_Opt%numPrtcl>0) call LPT_IO%Read_Restart()
      if(nrank==0) then
        call LPTLogInfo%OutInfo("Step3: Particles are READING from the Resarting file ...", 1 )
        call LPTLogInfo%OutInfo("Number of particles avaiable in domain:"//trim(num2str(LPT_opt%np_InDomain)),2)
      endif
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1
    if(nrank==0 .and. LPT_Opt%RestartFlag) call LPTLogInfo%OutInfo("Restart time [sec] :"//trim(num2str(t_res_tot)),2)

    ! Step4: Initializing visu
    call LPT_IO%Init_visu(chLPTPrm,2)
    
    ! Step5: initialize the inter-processors communication
    call LPTComm%InitComm()
    if(nrank==0) call LPTLogInfo%OutInfo("Step4: Initializing the inter-processors communication . . . ", 1 )
    
    ! Step6: timers for recording the execution time of different parts of program
    if(nrank==0) call LPTLogInfo%OutInfo("Step5: Initializing timers . . . ", 1 )
    call this%m_total_timer%reset()
    call this%m_integration_timer%reset()
    call this%m_comm_exchange_timer%reset()
    call this%m_write_prtcl_timer%reset()
  end subroutine LPT_Initialize

  !********************************************************************************
  !   iterating over time 
  !   calls all the required methods to do numIter iterations in the LPT system
  !********************************************************************************
  subroutine LPT_iterate(this,itime)
    implicit none
    class(LPTSystem) this
    integer,intent(in)::itime

    ! locals
    integer::ierror
    
    ! body
    call this%m_total_timer%start()

    ! correcting position and velocities 
    call this%m_integration_timer%start()
    iCountLPT=iCountLPT+1
    call Prtcl_Integrate(iCountLPT)
    call LPTContactSearchPW%FindContactsPW()
    call this%m_integration_timer%finish()

    ! inter-processor commucation for exchange
    call this%m_comm_exchange_timer%start()
    call LPTComm%Comm_For_Exchange()
    call this%m_comm_exchange_timer%finish()
    this%iterNumber = this%iterNumber + 1

    ! writing results to the output file and Restart file
    call this%m_write_prtcl_timer%start()
    call MPI_ALLREDUCE(GPrtcl_list%nlocal, LPT_Opt%np_InDomain, 1, int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    if( mod(itime,LPT_opt%SaveVisu)== 0)   call LPT_IO%dump_visu(itime)
    if( mod(itime,LPT_opt%BackupFreq)== 0 .or. itime==LPT_opt%ilast) then
      call LPT_IO%Write_Restart(itime)
      call LPT_IO%Delete_Prev_Restart(itime)
    endif
    call this%m_write_prtcl_timer%finish()
    call this%m_total_timer%finish()
    
    if(mod(itime,ivstats)==0) call ClcLPTStatistics()
    
    ! output to log file and terminal/command window
    IF((this%IterNumber==1 .or. mod(itime,LPT_opt%Cmd_LFile_Freq)==0) ) THEN
      if(nrank/=0) return
    
      ! command window and log file output
      call LPTLogInfo%OutInfo("LPT performed "//trim(num2str(itime))//" iterations up to here!",1)
      call LPTLogInfo%OutInfo("Execution time [tot, last, ave] [sec]: "//trim(num2str(this%m_total_timer%tot_time))//", "// &
      trim(num2str(this%m_total_timer%last_time ))//", "//trim(num2str(this%m_total_timer%average())),2)

      call LPTLogInfo%OutInfo("Integration time [tot, ave]        : "//trim(num2str(this%m_integration_timer%tot_time))//", "// &
      trim(num2str(this%m_integration_timer%average())), 3)

      call LPTLogInfo%OutInfo("Comm_For_Exchange [tot, ave]       : "//trim(num2str(this%m_comm_exchange_timer%tot_time))//", "// &
      trim(num2str(this%m_comm_exchange_timer%average())), 3)

      call LPTLogInfo%OutInfo("Write to file time [tot, ave]      : "//trim(num2str(this%m_write_prtcl_timer%tot_time))//", "// &
      trim(num2str(this%m_write_prtcl_timer%average())), 3)
     
      call LPTLogInfo%OutInfo("Particle number in  domain:  "//trim(num2str(LPT_Opt%np_InDomain)), 2)        
    ENDIF

  end subroutine LPT_iterate

  !**********************************************************************
  ! Write_LPT_Opt_to_Log
  !**********************************************************************
  subroutine Write_LPT_Opt_to_Log()
    implicit none

    ! locals
    real(RK)::dtLPT
    logical::RestartFlag,IsPeriodic(3)
    type(real3):: gravity, minpoint, maxpoint
    character(64):: RunName, ResultsDir,RestartDir
    integer::SaveVisuLPT,BackupFreqLPT, Cmd_LFile_Freq, LF_file_lvl, LF_cmdw_lvl
    integer::numPrtcl,PI_Method,numPrtcl_Type,ifirstLPT,ilastLPT
    NAMELIST /LPTOptions/ RestartFlag,numPrtcl,dtLPT,gravity,minpoint,maxpoint,PI_Method,numPrtcl_Type,RunName,RestartDir,   &
                          ResultsDir,BackupFreqLPT,SaveVisuLPT,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,ifirstLPT,ilastLPT,IsPeriodic

    RestartFlag = LPT_opt%RestartFlag 
    numPrtcl    = LPT_opt%numPrtcl 
    dtLPT       = LPT_opt%dt       
    ifirstLPT   = LPT_opt%ifirst   
    ilastLPT    = LPT_opt%ilast    
    gravity     = LPT_opt%gravity  
    minpoint    = LPT_opt%SimDomain_min 
    maxpoint    = LPT_opt%SimDomain_max 
    IsPeriodic  = LPT_opt%IsPeriodic 
    PI_Method   = LPT_opt%PI_Method
    write(RunName,"(A)")LPT_opt%RunName 
    write(ResultsDir,"(A)")LPT_opt%ResultsDir 

    numPrtcl_Type = LPT_opt%numPrtcl_Type
    write(RestartDir,"(A)") LPT_opt%RestartDir 
    BackupFreqLPT = LPT_opt%BackupFreq 
    SaveVisuLPT   = LPT_opt%SaveVisu 
    Cmd_LFile_Freq= LPT_opt%Cmd_LFile_Freq 
    LF_file_lvl   = LPT_opt%LF_file_lvl 
    LF_cmdw_lvl   = LPT_opt%LF_cmdw_lvl
    write(LPTLogInfo%nUnit, nml=LPTOptions)
  end subroutine Write_LPT_Opt_to_Log

end module LPT_System
