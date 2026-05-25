#include "definitions_inc.f90"
module ap_System
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
  use ap_Comm
  use ap_Property
  use ap_decomp_2d
  use ap_Variables
  use ap_IOAndVisu
  use ap_Parameters
  use ap_Statistics
  use ap_Integration
  use ap_ContactSearchPW
  implicit none
  private
    
  !// ATPSystem class 
  type ATPSystem
    integer :: iterNumber   = 0  ! iteration number 
        
    !// timers
    type(timer):: m_total_timer
    type(timer):: m_integration_timer
    type(timer):: m_write_prtcl_timer
    type(timer):: m_comm_exchange_timer
  contains
    procedure:: Initialize => ATP_Initialize
    procedure:: iterate    => ATP_iterate
  end type ATPSystem
  type(ATPSystem),public::ATP
  
  integer::iCountATP
contains

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
!   Initializing ATPSystem object with particles which are inserted from a 
!   predefined plane 
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ATP_Initialize(this,chATPPrm); implicit none
    class(ATPSystem)::this
    character(len=*),intent(in)::chATPPrm
    
    ! locals
    integer::ierr
    real(RK)::t_restart1,t_restart2,t_res_tot
    
    ! Initializing main log info and Visu
    iCountATP=0
    if(ATP_Opt%RestartFlag) iCountATP=10
    this%IterNumber=ATP_Opt%ifirst-1
    if(nrank==0) then
      call mkdir(ATP_Opt%ResultsDir, ierr); if(ierr<0) then; print*,'Cannot create folder ResultsDir'; stop; endif
      call mkdir(ATP_Opt%RestartDir, ierr); if(ierr<0) then; print*,'Cannot create folder RestartDir'; stop; endif    
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call ATPLogInfo%InitLog(ATP_opt%ResultsDir,ATP_opt%RunName,ATP_opt%LF_file_lvl,ATP_opt%LF_cmdw_lvl)
    if(nrank==0) call ATPLogInfo%CreateFile(ATP_opt%RunName)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call ATPLogInfo%OpenFile()
    if(nrank==0)call Write_ATP_Opt_to_Log()
    call InitATPStatistics(chATPPrm)

    ! Step1: Physical property
    call ATPProperty%InitPrtclProperty(chATPPrm)
    if(nrank==0) then
      call ATPLogInfo%OutInfo("Step1: Physical properties of particels and walls are set.",1)
      call ATPLogInfo%OutInfo("Physical properties contains "// strip( num2str(ATP_opt%numPrtcl_Type ) ) //" particle types ",2)
    endif

    ! Step2: initilize all the particle variables
    call GPrtcl_list%AllocateAllVar()
    call ATP_IO%Init_Visu(chATPPrm,1)
    t_restart1=MPI_WTIME()
    if(.not.ATP_Opt%RestartFlag) then
      if(ATP_Opt%numPrtcl>0) call GPrtcl_list%MakingAllPrtcl()
      if(nrank==0) then
        call ATPLogInfo%OutInfo("Step2: Initial Particle coordinates are MAKING into ATPSystem ...", 1 )
        call ATPLogInfo%OutInfo("Number of particles avaiable in the system:"//strip(num2str(ATP_opt%numPrtcl)),2)
      endif
      ATP_opt%np_InDomain = ATP_opt%numPrtcl
    else
      if(ATP_Opt%numPrtcl>0) call ATP_IO%Read_Restart()
      if(nrank==0) then
        call ATPLogInfo%OutInfo("Step3: Particles are READING from the Resarting file ...", 1 )
        call ATPLogInfo%OutInfo("Number of particles avaiable in domain:"//strip(num2str(ATP_opt%np_InDomain)),2)
      endif
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1
    if(nrank==0 .and. ATP_Opt%RestartFlag) call ATPLogInfo%OutInfo("Restart time [sec] :"//strip(num2str(t_res_tot)),2)

    ! Step3: Initializing Visu
    call ATP_IO%Init_Visu(chATPPrm,2)
    
    ! Step4: initialize the inter-processors communication
    call ATPComm%InitComm()
    if(nrank==0) call ATPLogInfo%OutInfo("Step4: Initializing the inter-processors communication . . . ", 1 )
    
    ! Step5: timers for recording the execution time of different parts of program
    if(nrank==0) call ATPLogInfo%OutInfo("Step5: Initializing timers . . . ", 1 )
    call this%m_total_timer%reset()
    call this%m_integration_timer%reset()
    call this%m_comm_exchange_timer%reset()
    call this%m_write_prtcl_timer%reset()
  end subroutine ATP_Initialize

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  !   iterating over time 
  !   calls all the required methods to do numIter iterations in the ATP system
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ATP_iterate(this,itime); implicit none
    class(ATPSystem) this
    integer,intent(in)::itime

    ! locals
    integer::ierr

    ! body
    call this%m_total_timer%start()

    ! correcting position and velocities 
    call this%m_integration_timer%start()
    iCountATP=iCountATP+1
    call Prtcl_Integrate(iCountATP)
    call ATPContactSearchPW%FindContactsPW()
    call this%m_integration_timer%finish()

    ! inter-processor commucation for exchange
    call this%m_comm_exchange_timer%start()
    call ATPComm%Comm_For_Exchange()
    call this%m_comm_exchange_timer%finish()
    this%iterNumber = this%iterNumber + 1

    ! writing results to the output file and Restart file
    call this%m_write_prtcl_timer%start()
    call MPI_ALLREDUCE(GPrtcl_list%nlocal, ATP_Opt%np_InDomain, 1, int_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    if( mod(itime,ATP_opt%SaveVisu)== 0)   call ATP_IO%dump_Visu(itime)
    if( mod(itime,ATP_opt%BackupFreq)== 0 .or. itime==ATP_opt%ilast) then
      call ATP_IO%Write_Restart(itime)
      call ATP_IO%Delete_Prev_Restart(itime)
    endif
    call this%m_write_prtcl_timer%finish()
    call this%m_total_timer%finish()
    
    if(mod(itime,ivstats)==0) call ClcATPStatistics()
    
    ! output to log file and terminal/command window
    IF((this%IterNumber==1 .or. mod(itime,ATP_opt%Cmd_LFile_Freq)==0) ) THEN
      if(nrank/=0) return
    
      ! command window and log file output
      call ATPLogInfo%OutInfo("ATP performed "//strip(num2str(itime))//" iterations up to here!",1)
      call ATPLogInfo%OutInfo("Execution time [tot, last, ave] [sec]: "//strip(num2str(this%m_total_timer%tot_time))//", "// &
      strip(num2str(this%m_total_timer%last_time ))//", "//strip(num2str(this%m_total_timer%average())),2)

      call ATPLogInfo%OutInfo("Integration time [tot, ave]        : "//strip(num2str(this%m_integration_timer%tot_time))//", "// &
      strip(num2str(this%m_integration_timer%average())), 3)

      call ATPLogInfo%OutInfo("Comm_For_Exchange [tot, ave]       : "//strip(num2str(this%m_comm_exchange_timer%tot_time))//", "// &
      strip(num2str(this%m_comm_exchange_timer%average())), 3)

      call ATPLogInfo%OutInfo("Write to file time [tot, ave]      : "//strip(num2str(this%m_write_prtcl_timer%tot_time))//", "// &
      strip(num2str(this%m_write_prtcl_timer%average())), 3)
     
      call ATPLogInfo%OutInfo("Particle number in  domain:  "//strip(num2str(ATP_Opt%np_InDomain)), 2)        
    ENDIF

  end subroutine ATP_iterate

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Write_ATP_Opt_to_Log
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Write_ATP_Opt_to_Log(); implicit none

    ! locals
    real(RK)::dtATP
    type(real3):: minpoint, maxpoint
    logical::RestartFlag,IsPeriodic(3)
    character(len=:),allocatable:: RunName, ResultsDir,RestartDir
    integer::SaveVisuATP,BackupFreqATP, Cmd_LFile_Freq, LF_file_lvl, LF_cmdw_lvl
    integer::numPrtcl,PI_Method,numPrtcl_Type,ifirstATP,ilastATP
    NAMELIST /ATPOptions/ RestartFlag,numPrtcl,dtATP,minpoint,maxpoint,PI_Method,numPrtcl_Type,RunName,RestartDir,   &
                          ResultsDir,BackupFreqATP,SaveVisuATP,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,ifirstATP,ilastATP,IsPeriodic

    RestartFlag = ATP_opt%RestartFlag
    numPrtcl    = ATP_opt%numPrtcl 
    dtATP       = ATP_opt%dt       
    ifirstATP   = ATP_opt%ifirst   
    ilastATP    = ATP_opt%ilast
    minpoint    = ATP_opt%SimDomain_min 
    maxpoint    = ATP_opt%SimDomain_max 
    IsPeriodic  = ATP_opt%IsPeriodic 
    PI_Method   = ATP_opt%PI_Method
    RunName     = ATP_opt%RunName 
    ResultsDir  = ATP_opt%ResultsDir 
    RestartDir  = ATP_opt%RestartDir
    numPrtcl_Type = ATP_opt%numPrtcl_Type
    BackupFreqATP = ATP_opt%BackupFreq 
    SaveVisuATP   = ATP_opt%SaveVisu 
    Cmd_LFile_Freq= ATP_opt%Cmd_LFile_Freq 
    LF_file_lvl   = ATP_opt%LF_file_lvl 
    LF_cmdw_lvl   = ATP_opt%LF_cmdw_lvl
    write(ATPLogInfo%iUnit, nml=ATPOptions)
  end subroutine Write_ATP_Opt_to_Log

end module ap_System
