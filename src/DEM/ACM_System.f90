module Prtcl_System
  use MPI
  use m_Timer
  use m_TypeDef
  use Prtcl_Comm
  use Prtcl_Property
  use Prtcl_Geometry
  use Prtcl_decomp_2d
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_IOAndVisu
  use Prtcl_DumpPrtcl
  use Prtcl_Parameters
  use Prtcl_Integration
  use Prtcl_ContactSearch
  use Prtcl_ContactSearchPW
  use m_Decomp2d,only:nrank
  implicit none
  private
    
  !// DEMSystem class 
  type DEMSystem
    integer :: iterNumber   = 0  ! iteration number 
        
    !// timers
    type(timer):: m_total_timer
    type(timer):: m_pre_iter_timer
    type(timer):: m_comm_cs_timer
    type(timer):: m_CSCF_PP_timer
    type(timer):: m_CSCF_PW_timer
    type(timer):: m_Acceleration_timer
    type(timer):: m_integration_timer
    type(timer):: m_write_prtcl_timer
    type(timer):: m_comm_exchange_timer
  contains
    procedure:: Initialize => DEMS_Initialize
    
    ! iterating simulation for n time steps
    procedure:: iterate     => DEMS_iterate
    
    ! performing pre-iterations 
    procedure:: preIteration    => DEMS_preIteration
    
  end type DEMSystem
  type(DEMSystem),public::DEM
  
  integer::iCountACM
contains

!********************************************************************************
!   Initializing DEMSystem object with particles which are inserted from a 
!   predefined plane 
!********************************************************************************
  subroutine  DEMS_Initialize(this,chDEMPrm)
    implicit none
    class(DEMSystem)::this
    character(*),intent(in)::chDEMPrm
    
    ! locals
    integer::ierror
    character(256)::chStr
    real(RK)::t_restart1,t_restart2,t_res_tot
    
    !// Initializing main log info
    iCountACM=0
    if(DEM_Opt%RestartFlag) iCountACM=10
    this%IterNumber=DEM_Opt%ifirst-1
    write(chStr,"(A)") 'mkdir -p '//DEM_Opt%ResultsDir//' '//DEM_Opt%RestartDir//' 2> /dev/null'
    if(nrank==0) call system(trim(adjustl(chStr)))
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call DEMLogInfo%InitLog(DEM_Opt%ResultsDir,DEM_Opt%RunName,DEM_Opt%LF_file_lvl,DEM_Opt%LF_cmdw_lvl)
    if(nrank==0) call DEMLogInfo%CreateFile(DEM_Opt%RunName)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call DEMLogInfo%OpenFile()
    if(nrank==0) call Write_DEM_Opt_to_Log()

    ! Step1: Physical property
    call DEMProperty%InitPrtclProperty(chDEMPrm)
    call DEMProperty%InitWallProperty(chDEMPrm)
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step1: Physical properties of particels and walls are set.",1)
      call DEMLogInfo%OutInfo("Physical properties contains "// trim( num2str(DEM_Opt%numPrtcl_Type ) ) // &
                              " particle types and "//trim( num2str(DEM_Opt%numWall_type ) )// " wall types.",2)
    endif

    ! Step2: set the geometry
    call DEMGeometry%MakeGeometry(chDEMPrm)
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step2: Geometry is set", 1 )
      call DEMLogInfo%OutInfo("Geometry Contains "//trim(num2str(DEMGeometry%num_pWall))//" Plane walls.", 2)
    endif

    ! Step3: initilize all the particle variables
    call GPrtcl_list%AllocateAllVar()
    call DEM_IO%Init_visu(chDEMPrm,1)
    t_restart1=MPI_WTIME()
    if(.not.DEM_Opt%RestartFlag) then
      if(DEM_Opt%numPrtcl>0) call DEM_IO%ReadInitialCoord()
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are MAKING into DEMSystem ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in the system:"//trim(num2str(DEM_Opt%numPrtcl)),2)
      endif
      DEM_Opt%np_InDomain = DEM_Opt%numPrtcl   
    else
      if(DEM_Opt%numPrtcl>0) then
        call DEM_IO%Read_Restart()
      else
        DEM_Opt%np_InDomain=0
      endif
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are READING from the Restarting file ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in domain:"//trim(num2str(DEM_Opt%np_InDomain)),2)
      endif
    endif
    if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedCoord()
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1

    ! Step4: Initializing visu
    call DEM_IO%Init_visu(chDEMPrm,2)
    call Initialize_DumpPrtcl(chDEMPrm)

    ! Step5: initialize the inter-processors communication
    call DEM_Comm%InitComm()
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step4: Initializing the inter-processors communication . . . ", 1 )
    endif

    ! Step6: Initializing contact list and contact force
    t_restart1=MPI_WTIME()
    call GPPW_CntctList%InitContactList()
    if(DEM_Opt%RestartFlag) then
      if(DEM_Opt%np_InDomain>0)call DEM_IO%RestartCL()
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1+t_res_tot
    if(nrank==0 .and. DEM_Opt%RestartFlag) call DEMLogInfo%OutInfo("Restart time [sec] :"//trim(num2str(t_res_tot)),2)
    
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step5: Initializing contact list and contact force models . . . ", 1 )
      if(DEM_Opt%CF_Type == ACM_LSD ) then
        write(chStr,"(A)") "Adaptive linear spring-dashpot model"
      elseif(DEM_Opt%CF_Type == ACM_nLin ) then
        write(chStr,"(A)") "Adaptive non-linear visco-elastic model"
      elseif(DEM_Opt%CF_Type == DEM_LSD ) then
        write(chStr,"(A)") "Typical linear spring-dashpot model"
      elseif(DEM_Opt%CF_Type == DEM_nLin ) then
        write(chStr,"(A)") "Typical non-linear visco-elastic model"
      endif
      call DEMLogInfo%OutInfo("Contact force model is "//trim(chStr), 2)

      if(DEM_Opt%PI_Method==PIM_FE) then
        write(chStr,"(A)") "Forward Euler             "
      elseif(DEM_Opt%PI_Method==PIM_AB2) then
        write(chStr,"(A)") "Adams Bashforth: 2nd Order"
      elseif(DEM_Opt%PI_Method==PIM_AB3) then
        write(chStr,"(A)") "Adams Bashforth: 3nd Order"
      endif
      call DEMLogInfo%OutInfo("Linear   movement Integration scheme is : "//trim(chStr),2)

      if(DEM_Opt%PRI_Method==PIM_FE) then
        write(chStr,"(A)") "Forward Euler             "
      elseif(DEM_Opt%PRI_Method==PIM_AB2) then
        write(chStr,"(A)") "Adams Bashforth: 2nd Order"
      elseif(DEM_Opt%PRI_Method==PIM_AB3) then
        write(chStr,"(A)") "Adams Bashforth: 3nd Order"
      endif
      call DEMLogInfo%OutInfo("Rotating movement Integration scheme is : "//trim(chStr),2)
    endif

    ! Step7: Initializing contact search method
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step6: Initializing contact search method . . . ", 1)
      call DEMLogInfo%OutInfo("Particle-Particle contact search intialization...",2)
      call DEMLogInfo%OutInfo("Particle-Wall contact search intialization...",2)
    endif
    call DEMContactSearch%InitContactSearch()
    call DEMContactSearchPW%InitContactSearchPW()
    
    ! Step8: timers for recording the execution time of different parts of program
    if(nrank==0) call DEMLogInfo%OutInfo("Step7: Initializing timers . . . ", 1 )
    call this%m_total_timer%reset()
    call this%m_pre_iter_timer%reset()
    call this%m_comm_cs_timer%reset()
    call this%m_CSCF_PP_timer%reset()
    call this%m_CSCF_PW_timer%reset()
    call this%m_Acceleration_timer%reset()
    call this%m_integration_timer%reset()
    call this%m_comm_exchange_timer%reset()
    call this%m_write_prtcl_timer%reset()
    
    call DEM_IO%dump_visu((DEM_Opt%ifirst-1)/icouple)            
  end subroutine DEMS_Initialize

  !********************************************************************************
  !   iterating over time 
  !   calls all the required methods to do numIter iterations in the DEM system
  !********************************************************************************
  subroutine DEMS_iterate(this,itime)
    implicit none
    class(DEMSystem) this
    integer,intent(in)::itime

    ! locals
    character(256)::chLine
    integer::Consv_Cont(2),Consv_Cont1(2),ierror,npwcs(4)

  IF(UpdateACMflag) THEN
    ! body
    call this%m_total_timer%start()

    ! pre-iteration adjustments 
    call this%m_pre_iter_timer%start()
    call this%preIteration()
    call this%m_pre_iter_timer%finish()

    ! inter-processor commucation for contact search ( ghost particle )
    call this%m_comm_cs_timer%start()
    call DEM_Comm%Comm_For_Cntct()
    call this%m_comm_cs_timer%finish()

    ! finding contacts among particels, and then calculating contact forces
    call this%m_CSCF_PP_timer%start()
    call DEMContactSearch%FindContacts()
    call this%m_CSCF_PP_timer%finish()

    ! finding contacts between particles and walls, and then calculating contact forces
    call this%m_CSCF_PW_timer%start()
    call DEMContactSearchPW%FindContactsPW()
    call this%m_CSCF_PW_timer%finish()

    ! correcting position and velocities 
    call this%m_integration_timer%start()
    iCountACM=iCountACM+1
    call Prtcl_Integrate(iCountACM)
    call this%m_integration_timer%finish()
   
    ! inter-processor commucation for exchange
    call this%m_comm_exchange_timer%start()
    call DEM_Comm%Comm_For_Exchange()
    call GPPW_CntctList%RemvReleased()
    call this%m_comm_exchange_timer%finish()
  ENDIF
    this%iterNumber = this%iterNumber + 1

    ! writing results to the output file and Restart file
    call this%m_write_prtcl_timer%start()
    call MPI_ALLREDUCE(GPrtcl_list%nlocal, DEM_Opt%np_InDomain, 1, int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    if( mod(this%IterNumber,DEM_Opt%SaveVisu)== 0)   call DEM_IO%dump_visu(itime/icouple)
    if( mod(this%IterNumber,DumpPrtclFreq)== 0)      call WriteDumpCache(itime)
    if( mod(this%IterNumber,DEM_Opt%BackupFreq)== 0 .or. itime==DEM_Opt%ilast) then
      call DEM_IO%Write_Restart(itime)
      call DEM_IO%Delete_Prev_Restart(itime)
      call PrtclVarDump(itime)
    endif
    call this%m_write_prtcl_timer%finish()
    call this%m_total_timer%finish()
            
  IF(UpdateACMflag) THEN    
    ! output to log file and terminal/command window
    IF((this%IterNumber==DEM_Opt%ifirst .or. mod(this%IterNumber,DEM_Opt%Cmd_LFile_Freq)==0) ) THEN
      Consv_Cont1 =  DEMContactSearch%get_numContact()
      call MPI_REDUCE(Consv_Cont1, Consv_Cont,       2,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
      call MPI_REDUCE(GPPW_CntctList%numCntcts,npwcs,4,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
      if(nrank/=0) return
    
      ! command window and log file output
      call DEMLogInfo%OutInfo("DEM performed "//trim(num2str(this%IterNumber))//" iterations up to here!",1)
      call DEMLogInfo%OutInfo("Execution time [tot, last, ave] [sec]: "//trim(num2str(this%m_total_timer%tot_time))//", "// &
      trim(num2str(this%m_total_timer%last_time ))//", "//trim(num2str(this%m_total_timer%average())),2)

      call DEMLogInfo%OutInfo("PreItertion time [tot, ave]        : "//trim(num2str(this%m_pre_iter_timer%tot_time))//", "// &
      trim(num2str(this%m_pre_iter_timer%average())),3)

      call DEMLogInfo%OutInfo("Comm_For_Contact [tot, ave]        : "//trim(num2str(this%m_comm_cs_timer%tot_time))//", "// &
      trim(num2str(this%m_comm_cs_timer%average())), 3)

      call DEMLogInfo%OutInfo("CS and CF P-P time [tot, ave]      : "//trim(num2str(this%m_CSCF_PP_timer%tot_time))//", "// &
      trim(num2str(this%m_CSCF_PP_timer%average())), 3)

      call DEMLogInfo%OutInfo("CS and CF P-W time [tot, ave]      : "//trim(num2str(this%m_CSCF_PW_timer%tot_time))//", "// &
      trim(num2str(this%m_CSCF_PW_timer%average())), 3)
      call DEMLogInfo%OutInfo("Integration time [tot, ave]        : "//trim(num2str(this%m_integration_timer%tot_time))//", "// &
      trim(num2str(this%m_integration_timer%average())), 3)
      call DEMLogInfo%OutInfo("Comm_For_Exchange [tot, ave]       : "//trim(num2str(this%m_comm_exchange_timer%tot_time))//", "// &
      trim(num2str(this%m_comm_exchange_timer%average())), 3)
      call DEMLogInfo%OutInfo("Write to file time [tot, ave]      : "//trim(num2str(this%m_write_prtcl_timer%tot_time))//", "// &
      trim(num2str(this%m_write_prtcl_timer%average())), 3)
      write(chLine,"(A)") "Particle number in  domain:  "//trim(num2str(DEM_Opt%np_InDomain))
      call DEMLogInfo%OutInfo(chLine, 2)
      call DEMLogInfo%OutInfo("Contact information", 2)
      write(chLine,"(A)") "No. consrvtv. contacts, same level | cross level: "// trim(num2str(Consv_Cont(1)))//" | "//trim(num2str(Consv_Cont(2)))
      call DEMLogInfo%OutInfo(chLine, 3)
      write(chLine,"(A)") "No. exact contacts P-P | P-GP | P-FP | P-W : "//trim(num2str(npwcs(1)))//" | "//trim(num2str(npwcs(2)))//" | "//trim(num2str(npwcs(3)))//" | "//trim(num2str(npwcs(4)))
      call DEMLogInfo%OutInfo(chLine, 3)
    ENDIF
  ENDIF
  end subroutine DEMS_iterate
    
  !**********************************************************************
  ! DEMS_preIteration
  !**********************************************************************
  subroutine DEMS_preIteration(this)
    implicit none
    class(DEMSystem)::this
        
    ! update wall neighbor list if necessary
    call DEMContactSearchPW%UpdateNearPrtclsPW(this%iterNumber)
    GPrtcl_cntctForce =zero_r3
    GPrtcl_torque= zero_r3
    GPrtcl_HighSt= "N"
    call GPPW_CntctList%PreIteration()
  end subroutine DEMS_preIteration

  !**********************************************************************
  ! Write_DEM_Opt_to_Log
  !**********************************************************************
  subroutine Write_DEM_Opt_to_Log()
    implicit none

    ! locals
    logical::RestartFlag
    real(RK)::dtDEM,Wall_neighbor_ratio,Prtcl_cs_ratio
    character(64):: RunName,ResultsDir,RestartDir,Geom_Dir
    integer::Wall_max_update_iter,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl
    integer::numPrtcl,numPrtclFix,CS_Method,CF_Type,PI_Method,PRI_Method,GeometrySource
    integer::numPrtcl_Type,numWall_type,ncvAllowed,CntctList_Size,CS_numlvls,Base_wall_id
    
    logical,dimension(3)::IsPeriodic
    type(real3)::gravity,minpoint,maxpoint
    integer::ifirstDEM,ilastDEM,BackupFreqDEM,SaveVisuDEM
    NAMELIST /DEMOptions/ RestartFlag,numPrtcl,numPrtclFix,dtDEM,gravity,minpoint,maxpoint,CS_Method,CF_Type,   &
                          PI_Method,PRI_Method,numPrtcl_Type,numWall_type,CS_numlvls,ncvAllowed,CntctList_Size, &
                          RunName,Wall_max_update_iter,Wall_neighbor_ratio,ResultsDir,RestartDir,BackupFreqDEM, &
                          SaveVisuDEM,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,GeometrySource,Geom_Dir,ifirstDEM, &
                          ilastDEM,Prtcl_cs_ratio,IsPeriodic,Base_wall_id
#ifdef CFDACM
    NAMELIST/CFDACMCoupling/UpdateACMflag,icouple,nForcingExtra,IBM_Scheme,Klub_pp,Klub_pw,Lub_ratio, &
                            Ndt_coll,IsDryColl,St_Crit,IsAddFluidPressureGradient
#endif

    RestartFlag = DEM_Opt%RestartFlag 
    numPrtcl    = DEM_Opt%numPrtcl   
    numPrtclFix = DEM_Opt%numPrtclFix  
    dtDEM       = DEM_Opt%dt       
    ifirstDEM   = DEM_Opt%ifirst   
    ilastDEM    = DEM_Opt%ilast    
    gravity     = DEM_Opt%gravity  
    minpoint    = DEM_Opt%SimDomain_min 
    maxpoint    = DEM_Opt%SimDomain_max 
    IsPeriodic  = DEM_Opt%IsPeriodic 
           
    Prtcl_cs_ratio=  DEM_Opt%Prtcl_cs_ratio 
    CS_Method  = DEM_Opt%CS_Method 
    CF_Type    = DEM_Opt%CF_Type
    PI_Method  = DEM_Opt%PI_Method  
    PRI_Method = DEM_Opt%PRI_Method
           
    numPrtcl_Type = DEM_Opt%numPrtcl_Type
    numWall_type  = DEM_Opt%numWall_type 
    
    ncvAllowed     = DEM_Opt%ncvAllowed
    CntctList_Size = DEM_Opt%CntctList_Size 
    CS_numlvls     = DEM_Opt%CS_numlvls    
           
    Base_wall_id  =  DEM_Opt%Base_wall_id 
    Wall_max_update_iter = DEM_Opt%Wall_max_update_iter
    Wall_neighbor_ratio  = DEM_Opt%Wall_neighbor_ratio 
           
    write(RunName,"(A)") DEM_Opt%RunName 
    write(ResultsDir,"(A)") DEM_Opt%ResultsDir 
    write(RestartDir,"(A)")DEM_Opt%RestartDir 
    BackupFreqDEM = DEM_Opt%BackupFreq 
    SaveVisuDEM   = DEM_Opt%SaveVisu 
    Cmd_LFile_Freq= DEM_Opt%Cmd_LFile_Freq 
    LF_file_lvl   = DEM_Opt%LF_file_lvl 
    LF_cmdw_lvl   = DEM_Opt%LF_cmdw_lvl 
           
    GeometrySource=  DEM_Opt%GeometrySource 
    write(Geom_Dir,"(A)") DEM_Opt%Geom_Dir 
    write(DEMLogInfo%nUnit, nml=DEMOptions)
#ifdef CFDACM
    write(DEMLogInfo%nUnit, nml=CFDACMCoupling)
#endif
  end subroutine Write_DEM_Opt_to_Log
end module Prtcl_System
