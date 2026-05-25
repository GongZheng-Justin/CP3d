#include "definitions_inc.f90"
module lp_System
  use MPI
  use mc_Timer
  use mc_TypeDef
  use mc_Decomp2d,only:nrank
  use mc_FileOperator,only:mkdir
#ifdef CFDSecondOrder
  use f2_Parameters,only:ivstats
#else
  use f4_Parameters,only:ivstats
#endif
  use lp_Comm
  use lp_Property
  use lp_Geometry
  use lp_decomp_2d
  use lp_Variables
  use lp_IOAndVisu
  use lp_Parameters
  use lp_Statistics
  use lp_Integration
  use lp_ContactSearchPW
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

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
!   Initializing LPTSystem object with particles which are inserted from a 
!   predefined plane 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine LPT_Initialize(this,chLPTPrm); implicit none
    class(LPTSystem)::this
    character(len=*),intent(in)::chLPTPrm
    
    ! locals
    integer::ierr
    real(RK)::t_restart1,t_restart2,t_res_tot
    
    ! Initializing main log info and visu
    iCountLPT=0
    if(LPT_Opt%RestartFlag) iCountLPT=10
    this%IterNumber=LPT_Opt%ifirst-1
    if(nrank==0) then
      call mkdir(LPT_Opt%ResultsDir, ierr); if(ierr<0) then; print*,'Cannot create folder ResultsDir'; stop; endif
      call mkdir(LPT_Opt%RestartDir, ierr); if(ierr<0) then; print*,'Cannot create folder RestartDir'; stop; endif
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call LPTLogInfo%InitLog(LPT_Opt%ResultsDir,LPT_Opt%RunName,LPT_Opt%LF_file_lvl,LPT_Opt%LF_cmdw_lvl)
    if(nrank==0) call LPTLogInfo%CreateFile(LPT_Opt%RunName)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call LPTLogInfo%OpenFile()
    if(nrank==0)call Write_LPT_Opt_to_Log()
    call InitLPTStatistics(chLPTPrm)

    ! Step1: Physical property
    call LPTProperty%InitPrtclProperty(chLPTPrm)
    if(nrank==0) then
      call LPTLogInfo%OutInfo("Step1: Physical properties of particels and walls are set.",1)
      call LPTLogInfo%OutInfo("Physical properties contains "// strip( num2str(LPT_Opt%numPrtcl_Type ) ) //" particle types ",2)
    endif

    ! Step2: set the geometry
    call LPTGeometry%MakeGeometry()
    if(nrank==0) then
      call LPTLogInfo%OutInfo("Step2: Geometry is set", 1 )
      call LPTLogInfo%OutInfo("Geometry Contains "//strip(num2str(LPTGeometry%num_pWall))//" Plane walls.", 2)
    endif

    ! Step3: initilize all the particle variables
    call GPrtcl_list%AllocateAllVar()
    call LPT_IO%Init_visu(chLPTPrm,1)
    t_restart1=MPI_WTIME()
    if(.not.LPT_Opt%RestartFlag) then
      if(LPT_Opt%numPrtcl>0) call GPrtcl_list%MakingAllPrtcl(chLPTPrm)
      if(nrank==0) then
        call LPTLogInfo%OutInfo("Step3: Initial Particle coordinates are MAKING into LPTSystem ...", 1 )
        call LPTLogInfo%OutInfo("Number of particles avaiable in the system:"//strip(num2str(LPT_Opt%numPrtcl)),2)
      endif
      LPT_Opt%np_InDomain = LPT_Opt%numPrtcl
    else
      if(LPT_Opt%numPrtcl>0) call LPT_IO%Read_Restart()
      if(nrank==0) then
        call LPTLogInfo%OutInfo("Step3: Particles are READING from the Resarting file ...", 1 )
        call LPTLogInfo%OutInfo("Number of particles avaiable in domain:"//strip(num2str(LPT_Opt%np_InDomain)),2)
      endif
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1
    if(nrank==0 .and. LPT_Opt%RestartFlag) call LPTLogInfo%OutInfo("Restart time [sec] :"//strip(num2str(t_res_tot)),2)

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

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  !   iterating over time 
  !   calls all the required methods to do numIter iterations in the LPT system
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine LPT_iterate(this,itime); implicit none
    class(LPTSystem):: this
    integer,intent(in)::itime

    ! locals
    integer::ierr
    
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
    call MPI_ALLREDUCE(GPrtcl_list%nlocal, LPT_Opt%np_InDomain, 1, int_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    if( mod(itime,LPT_Opt%SaveVisu)== 0)   call LPT_IO%dump_visu(itime)
    if( mod(itime,LPT_Opt%BackupFreq)== 0 .or. itime==LPT_Opt%ilast) then
      call LPT_IO%Write_Restart(itime)
      call LPT_IO%Delete_Prev_Restart(itime)
    endif
    call this%m_write_prtcl_timer%finish()
    call this%m_total_timer%finish()
    
    if(mod(itime,ivstats)==0) call ClcLPTStatistics()
    
    ! output to log file and terminal/command window
    IF((this%IterNumber==1 .or. mod(itime,LPT_Opt%Cmd_LFile_Freq)==0) ) THEN
      if(nrank/=0) return
    
      ! command window and log file output
      call LPTLogInfo%OutInfo("LPT performed "//strip(num2str(itime))//" iterations up to here!",1)
      call LPTLogInfo%OutInfo("Execution time [tot, last, ave] [sec]: "//strip(num2str(this%m_total_timer%tot_time))//", "// &
      strip(num2str(this%m_total_timer%last_time ))//", "//strip(num2str(this%m_total_timer%average())),2)

      call LPTLogInfo%OutInfo("Integration time [tot, ave]        : "//strip(num2str(this%m_integration_timer%tot_time))//", "// &
      strip(num2str(this%m_integration_timer%average())), 3)

      call LPTLogInfo%OutInfo("Comm_For_Exchange [tot, ave]       : "//strip(num2str(this%m_comm_exchange_timer%tot_time))//", "// &
      strip(num2str(this%m_comm_exchange_timer%average())), 3)

      call LPTLogInfo%OutInfo("Write to file time [tot, ave]      : "//strip(num2str(this%m_write_prtcl_timer%tot_time))//", "// &
      strip(num2str(this%m_write_prtcl_timer%average())), 3)
     
      call LPTLogInfo%OutInfo("Particle number in  domain:  "//strip(num2str(LPT_Opt%np_InDomain)), 2)        
    ENDIF
  end subroutine LPT_iterate

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Write_LPT_Opt_to_Log
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Write_LPT_Opt_to_Log(); implicit none

    ! locals
    real(RK)::dtLPT
    logical::RestartFlag,IsPeriodic(3)
    type(real3):: gravity, minpoint, maxpoint
    character(len=:),allocatable:: RunName, ResultsDir,RestartDir
    integer::SaveVisuLPT,BackupFreqLPT, Cmd_LFile_Freq, LF_file_lvl, LF_cmdw_lvl
    integer::numPrtcl,PI_Method,numPrtcl_Type,ifirstLPT,ilastLPT
    NAMELIST /LPTOptions/ RestartFlag,numPrtcl,dtLPT,gravity,minpoint,maxpoint,PI_Method,numPrtcl_Type,RunName,RestartDir,   &
                          ResultsDir,BackupFreqLPT,SaveVisuLPT,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,ifirstLPT,ilastLPT,IsPeriodic

    RestartFlag= LPT_Opt%RestartFlag 
    numPrtcl   = LPT_Opt%numPrtcl 
    dtLPT      = LPT_Opt%dt       
    ifirstLPT  = LPT_Opt%ifirst   
    ilastLPT   = LPT_Opt%ilast    
    gravity    = LPT_Opt%gravity  
    minpoint   = LPT_Opt%SimDomain_min 
    maxpoint   = LPT_Opt%SimDomain_max 
    IsPeriodic = LPT_Opt%IsPeriodic 
    PI_Method  = LPT_Opt%PI_Method
    RunName    = LPT_Opt%RunName
    ResultsDir = LPT_Opt%ResultsDir
    RestartDir = LPT_Opt%RestartDir
    
    numPrtcl_Type = LPT_Opt%numPrtcl_Type
    BackupFreqLPT = LPT_Opt%BackupFreq 
    SaveVisuLPT   = LPT_Opt%SaveVisu 
    Cmd_LFile_Freq= LPT_Opt%Cmd_LFile_Freq 
    LF_file_lvl   = LPT_Opt%LF_file_lvl 
    LF_cmdw_lvl   = LPT_Opt%LF_cmdw_lvl
    write(LPTLogInfo%iUnit, nml=LPTOptions)
  end subroutine Write_LPT_Opt_to_Log

end module lp_System
