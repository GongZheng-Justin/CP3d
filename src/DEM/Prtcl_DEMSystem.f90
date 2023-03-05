module Prtcl_DEMSystem
  use MPI
  use m_Timer
  use m_TypeDef
  use Prtcl_Comm
  use Prtcl_Property
  use Prtcl_Geometry
  use Prtcl_IOAndVisu
  use Prtcl_Variables
  use Prtcl_decomp_2d
  use Prtcl_CL_and_CF
  use Prtcl_Parameters
  use Prtcl_Integration
  use Prtcl_ContactSearch
  use Prtcl_ContactSearchPW
#ifdef CFDDEM
  use Prtcl_DumpPrtcl
  use m_Decomp2d,only: nrank
#endif  
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
   
  integer::iCountDEM
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
    integer:: ierror
    character(256):: chStr
    real(RK)::t_restart1,t_restart2,t_res_tot
    
    !// Initializing main log info and visu
    iCountDEM=0
    if(DEM_Opt%RestartFlag) iCountDEM=10
    this%IterNumber=DEM_Opt%ifirst-1
    write(chStr,"(A)") 'mkdir -p '//DEM_opt%ResultsDir//' '//DEM_opt%RestartDir//' 2> /dev/null'
    if(nrank==0) call system(trim(adjustl(chStr)))
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call DEMLogInfo%InitLog(DEM_opt%ResultsDir,DEM_opt%RunName,DEM_opt%LF_file_lvl,DEM_opt%LF_cmdw_lvl)
    if(nrank==0) call DEMLogInfo%CreateFile(DEM_opt%RunName)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call DEMLogInfo%OpenFile()
    if(nrank==0) call Write_DEM_Opt_to_Log()

    ! Step1: Physical property
    call DEMProperty%InitPrtclProperty(chDEMPrm)
    call DEMProperty%InitWallProperty(chDEMPrm)
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step1: Physical properties of particels and walls are set.",1)
      call DEMLogInfo%OutInfo("Physical properties contains "// trim( num2str(DEM_opt%numPrtcl_Type ) ) // &
                              " particle types and "//trim( num2str(DEM_opt%numWall_type ) )// " wall types.",2)
    endif

    ! Step2: set the geometry
    call DEMGeometry%MakeGeometry(chDEMPrm)
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step2: Geometry is set", 1 )
      call DEMLogInfo%OutInfo("Geometry Contains "//trim(num2str(DEMGeometry%num_pWall))//" Plane walls.", 2)
    endif

    ! Step3: initilize all the particle variables and IO
    call GPrtcl_list%AllocateAllVar()
    call DEM_IO%Init_visu(chDEMPrm,1)
    t_restart1=MPI_WTIME()
#ifdef CFDDEM
    if(.not.DEM_Opt%RestartFlag) then
      if(DEM_Opt%numPrtcl>0) call DEM_IO%ReadInitialCoord()
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Initial Particle coordinates are READING into DEMSystem ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in the system:"//trim(num2str(DEM_opt%numPrtcl)),2)
      endif
      DEM_opt%np_InDomain = DEM_opt%numPrtcl
      if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedCoord()
    else
      if(DEM_Opt%numPrtcl>0) then
        call DEM_IO%Read_Restart()
      else
        DEM_opt%np_InDomain=0
      endif
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are READING from the Resarting file ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in domain:"//trim(num2str(DEM_opt%np_InDomain)),2)
      endif
      if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedRestart()
    endif
#else
    if(.not.DEM_Opt%RestartFlag) then
      call GPrtcl_list%MakingAllPrtcl(chDEMPrm)
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are MAKING into DEMSystem ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in the system:"//trim(num2str(DEM_opt%numPrtcl)),2)
      endif
      DEM_opt%np_InDomain = DEM_opt%numPrtcl
    else
      call DEM_IO%Read_Restart()
      if(nrank==0) then
        call DEMLogInfo%OutInfo("Step3: Particles are READING from the Resarting file ...", 1 )
        call DEMLogInfo%OutInfo("Number of particles avaiable in domain:"//trim(num2str(DEM_opt%np_InDomain)),2)
      endif
    endif
    if(DEM_Opt%numPrtclFix>0) call DEM_IO%ReadFixedCoord()
#endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1

    call DEM_IO%Init_visu(chDEMPrm,2)
#ifdef CFDDEM
    call DEM_PDump%Initialize(chDEMPrm)
#endif

    ! Step4: initialize the inter-processors communication
    call DEM_Comm%InitComm()
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step4: Initializing the inter-processors communication . . . ", 1 )
    endif

    ! Step5: Initializing contact list and contact force 
    call  GPPW_CntctList%InitContactList()
    t_restart1=MPI_WTIME()
    if(DEM_Opt%RestartFlag) then
      if(DEM_opt%np_InDomain>0)call DEM_IO%RestartCL()
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    t_restart2=MPI_WTIME(); t_res_tot=t_restart2-t_restart1+t_res_tot
    if(nrank==0 .and. DEM_Opt%RestartFlag) call DEMLogInfo%OutInfo("Restart time [sec] :"//trim(num2str(t_res_tot)),2)
    
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step5: Initializing contact list and contact force models . . . ", 1 )
      if( DEM_opt%CF_Type == DEM_LSD ) then
        write(chStr,"(A)") "linear spring-dashpot with limited tangential displacement"
      elseif( DEM_opt%CF_Type == DEM_nLin ) then
        write(chStr,"(A)") "non-linear visco-elastic model with limited tangential displacement"
      endif
      call DEMLogInfo%OutInfo("Contact force model is "//trim(chStr), 2 )

      if(DEM_opt%PI_Method==PIM_FE) then
        write(chStr,"(A)") "Forward Euler             "
      elseif(DEM_opt%PI_Method==PIM_AB2) then
        write(chStr,"(A)") "Adams Bashforth: 2nd Order"
      elseif(DEM_opt%PI_Method==PIM_AB3) then
        write(chStr,"(A)") "Adams Bashforth: 3nd Order"
      endif
      call DEMLogInfo%OutInfo("Linear   movement Integration scheme is : "//trim(chStr),2)

      if(DEM_opt%PRI_Method==PIM_FE) then
        write(chStr,"(A)") "Forward Euler             "
      elseif(DEM_opt%PRI_Method==PIM_AB2) then
        write(chStr,"(A)") "Adams Bashforth: 2nd Order"
      elseif(DEM_opt%PRI_Method==PIM_AB3) then
        write(chStr,"(A)") "Adams Bashforth: 3nd Order"
      endif
      call DEMLogInfo%OutInfo("Rotating movement Integration scheme is : "//trim(chStr),2)
    endif

    ! Step6: Initializing contact search method
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Step6: Initializing contact search method . . . ", 1)
      call DEMLogInfo%OutInfo("Particle-Particle contact search intialization...",2)
      call DEMLogInfo%OutInfo("Particle-Wall contact search intialization...",2)
    endif
    call DEMContactSearch%InitContactSearch()
    call DEMContactSearchPW%InitContactSearchPW()
    
    ! Step7: timers for recording the execution time of different parts of program
    if(nrank==0) call DEMLogInfo%OutInfo("Step7: Initializing timers . . . ", 1 )
    call this%m_total_timer%reset()
    call this%m_pre_iter_timer%reset()
    call this%m_comm_cs_timer%reset()
    call this%m_CSCF_PP_timer%reset()
    call this%m_CSCF_PW_timer%reset()
    call this%m_integration_timer%reset()
    call this%m_write_prtcl_timer%reset()
    call this%m_comm_exchange_timer%reset()
    
#ifdef CFDDEM
    call DEM_IO%dump_visu((DEM_opt%ifirst-1)/icouple)            
#else
    call DEM_IO%dump_visu(DEM_opt%ifirst-1)
#endif
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
    character(256):: chLine
    integer:: Consv_Cont(2),Consv_Cont1(2),ierror,npwcs(4)

#ifdef CFDDEM
  IF(UpdateDEMflag) THEN
#endif
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

    ! calculate linear and angular accelerations, position and velocities 
    call this%m_integration_timer%start()
    iCountDEM=iCountDEM+1
    call Prtcl_Integrate(iCountDEM)
    call this%m_integration_timer%finish()
   
    ! inter-processor commucation for exchange
    call this%m_comm_exchange_timer%start()
    call DEM_Comm%Comm_For_Exchange()
    call GPPW_CntctList%RemvReleased()
    call this%m_comm_exchange_timer%finish()
#ifdef CFDDEM
  ENDIF
#endif
    this%iterNumber = this%iterNumber + 1

    ! writing results to the output file and Restart file
    call this%m_write_prtcl_timer%start()
    call MPI_ALLREDUCE(GPrtcl_list%nlocal, DEM_Opt%np_InDomain, 1, int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
#ifdef CFDDEM
    if( mod(this%IterNumber,DEM_opt%SaveVisu)== 0)   call DEM_IO%dump_visu(itime/icouple)
    if( mod(this%IterNumber,DEM_PDump%WriteCacheFreq)== 0)  call DEM_PDump%WriteCache(itime)
#else
    if( mod(this%IterNumber,DEM_opt%SaveVisu)== 0)   call DEM_IO%dump_visu(itime)
#endif
    if( mod(this%IterNumber,DEM_opt%BackupFreq)== 0 .or. itime==DEM_opt%ilast) then
      call DEM_IO%Write_Restart(itime)
#ifdef CFDDEM
      call DEM_IO%WriteFixedRestart(itime)
      call DEM_PDump%PrtclVarDump(itime)
#endif
      call DEM_IO%Delete_Prev_Restart(itime)
    endif
    call this%m_write_prtcl_timer%finish()
    call this%m_total_timer%finish()

#ifdef CFDDEM
  IF(UpdateDEMflag) THEN
#endif    
    ! output to log file and terminal/command window
    IF((this%IterNumber==DEM_Opt%ifirst .or. mod(this%IterNumber,DEM_opt%Cmd_LFile_Freq)==0) ) THEN
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
#ifdef CFDDEM
  ENDIF
#endif
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
    integer::numPrtcl_Type,numWall_type,CntctList_Size,CS_numlvls,Base_wall_id
    integer::numPrtcl,numPrtclFix,CS_Method,CF_Type,PI_Method,PRI_Method,GeometrySource
    
    logical,dimension(3)::IsPeriodic
    type(real3)::gravity,minpoint,maxpoint
    integer::ifirstDEM,ilastDEM,BackupFreqDEM,SaveVisuDEM
    NAMELIST /DEMOptions/ RestartFlag,numPrtcl,numPrtclFix,dtDEM,gravity,minpoint,maxpoint,CS_Method,CF_Type,     &
                          PI_Method,PRI_Method,numPrtcl_Type,numWall_type,CS_numlvls,CntctList_Size,Base_wall_id, &
                          Wall_max_update_iter,Wall_neighbor_ratio,RunName,ResultsDir,RestartDir,BackupFreqDEM,   &
                          SaveVisuDEM,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,GeometrySource,Geom_Dir,ifirstDEM,   &
                          ilastDEM,Prtcl_cs_ratio,IsPeriodic
#ifdef CFDDEM
    NAMELIST/CFDDEMCoupling/icouple,UpdateDEMflag,is_clc_Lift,is_clc_Basset,is_clc_Basset_fixed,is_clc_ViscousForce,&
                            is_clc_PressureGradient,is_clc_ViscousForce,is_clc_FluidAcc,FluidAccCoe,SaffmanConst,   &
                            RatioSR,IsAddFluidPressureGradient
    NAMELIST/BassetOptions/ mWinBasset, mTailBasset, BassetAccuracy, BassetTailType
#endif

    RestartFlag = DEM_opt%RestartFlag 
    numPrtcl    = DEM_opt%numPrtcl   
    numPrtclFix = DEM_opt%numPrtclFix  
    dtDEM       = DEM_opt%dt       
    ifirstDEM   = DEM_opt%ifirst   
    ilastDEM    = DEM_opt%ilast    
    gravity     = DEM_opt%gravity  
    minpoint    = DEM_opt%SimDomain_min 
    maxpoint    = DEM_opt%SimDomain_max 
    IsPeriodic  = DEM_opt%IsPeriodic 
           
    Prtcl_cs_ratio=  DEM_opt%Prtcl_cs_ratio 
    CS_Method  = DEM_opt%CS_Method 
    CF_Type    = DEM_opt%CF_Type
    PI_Method  = DEM_opt%PI_Method  
    PRI_Method = DEM_opt%PRI_Method
           
    numPrtcl_Type = DEM_opt%numPrtcl_Type
    numWall_type  = DEM_opt%numWall_type 
          
    CntctList_Size = DEM_opt%CntctList_Size 
    CS_numlvls     = DEM_opt%CS_numlvls    
           
    Base_wall_id   = DEM_opt%Base_wall_id 
    Wall_max_update_iter = DEM_opt%Wall_max_update_iter
    Wall_neighbor_ratio  = DEM_opt%Wall_neighbor_ratio 
           
    write(RunName,"(A)") DEM_opt%RunName 
    write(ResultsDir,"(A)") DEM_opt%ResultsDir 
    write(RestartDir,"(A)")DEM_opt%RestartDir 
    BackupFreqDEM = DEM_opt%BackupFreq 
    SaveVisuDEM   = DEM_opt%SaveVisu 
    Cmd_LFile_Freq= DEM_opt%Cmd_LFile_Freq 
    LF_file_lvl   = DEM_opt%LF_file_lvl 
    LF_cmdw_lvl   = DEM_opt%LF_cmdw_lvl 
           
    GeometrySource=  DEM_opt%GeometrySource 
    write(Geom_Dir,"(A)") DEM_opt%Geom_Dir 
    write(DEMLogInfo%nUnit, nml=DEMOptions)
#ifdef CFDDEM
    write(DEMLogInfo%nUnit, nml=CFDDEMCoupling)
    write(DEMLogInfo%nUnit, nml=BassetOptions)
#endif
  end subroutine Write_DEM_Opt_to_Log

end module Prtcl_DEMSystem
