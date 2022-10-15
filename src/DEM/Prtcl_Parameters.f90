module Prtcl_Parameters
  use Prtcl_TypeDef
#ifdef CFDDEM
  use m_Parameters,only: dtMax,ifirst,ilast,BackupFreq,SaveVisu,xlx,yly,zlz,BcOption
#endif
#ifdef CFDACM
  use m_Decomp2d,only: nrank
  use m_Parameters,only: dtMax,ifirst,ilast,BackupFreq,SaveVisu,xlx,yly,zlz,BcOption,ischeme,FI_AB2,FI_RK2,FI_RK3
#endif
  implicit none
  private

  integer,parameter,public:: x_axis = 1
  integer,parameter,public:: y_axis = 2
  integer,parameter,public:: z_axis = 3    
    
  integer,parameter,public:: CSM_NBS_Munjiza       = 1
  integer,parameter,public:: CSM_NBS_Munjiza_Hrchl = 2
  integer,parameter,public:: DEM_LSD = 1
  integer,parameter,public:: DEM_nLin= 2
#ifdef CFDACM
  integer,parameter,public:: ACM_LSD = 3
  integer,parameter,public:: ACM_nLin= 4
#endif
  integer,parameter,public:: PIM_FE  = 1
  integer,parameter,public:: PIM_AB2 = 2
  integer,parameter,public:: PIM_AB3 = 3

#ifdef CFDACM
  real(RK),public::Klub_pp, Klub_pw,Lub_ratio,Ndt_coll,St_Crit
  logical,public::UpdateACMflag,IsDryColl,IsAddFluidPressureGradient
  integer,public::icouple,nForcingExtra,IBM_Scheme,idem_advance_start(3),idem_advance_end(3)
#endif
#ifdef CFDDEM
  integer,public::  icouple
  real(RK),public:: RatioSR                            ! search radius ratio
  real(RK),public:: FluidAccCoe,SaffmanConst
  logical,public::  IsAddFluidPressureGradient,UpdateDEMflag
  logical,public::  is_clc_Lift, is_clc_Basset,is_clc_Basset_fixed, is_clc_ViscousForce, is_clc_PressureGradient,is_clc_FluidAcc
  integer,public::  mWinBasset, mTailBasset, BassetAccuracy, BassetTailType 
  type BassetDataSeq
    integer:: HistoryStage=0
    integer:: HistStageFix=0
    integer:: nDataLen
    integer:: iWindowStart
    integer:: iWindowEnd
    integer:: iTailStart
    integer:: iTailEnd
  end type BassetDataSeq
  type(BassetDataSeq),public:: GPrtcl_BassetSeq
#endif
   
  ! default values  
  type DEMS_Options
    logical:: RestartFlag=.false.
    integer:: numPrtcl    = 8000     ! total particle number
    integer:: numPrtclFix = 1000     ! total fixed particle number
    integer:: np_InDomain            ! particle in domain
    integer:: ifirst                 ! first time step
    integer:: ilast                  ! last time step
    real(RK):: dt   =  1.0E-5_RK     ! time step       
    type(real3):: gravity = real3(zero,-9.81_RK,zero) ! gravity or other constant body forces if any
    type(real3):: SimDomain_min
    type(real3):: SimDomain_max
    logical,dimension(3):: IsPeriodic = .false.
        
    real(RK):: Prtcl_cs_ratio
    integer:: CS_Method = CSM_NBS_Munjiza      ! contact search method        
    integer:: CF_Type   = DEM_LSD              ! contact force type
    integer:: PI_Method = PIM_FE   ! integration scheme for translational motion
    integer:: PRI_Method = PIM_FE  ! integration scheme for rotational motion        
        
    integer:: numPrtcl_Type=1      ! number of particle type
    integer:: numWall_type =1      ! number of wall type
        
    ! contact list size, 6 means every particle can contact with 12 neighbour particles/walls in average 
    integer:: CntctList_Size = 6
    ! means default behavior, number of levels in multi-level contact search 
    integer:: CS_numlvls = 0  
    
    ! Global wall id will start at (Base_wall_id+1). Please make sure that  Base_wall_id≥numPrtcl    
    integer:: Base_wall_id = 10000000  
    ! Near wall list will be updated no more than every 100 iterations
    integer:: Wall_max_update_iter = 100
    ! The particle withthin 2*MaxRadius, will be considered into the NEAR WALL LIST 
    real(RK):: Wall_neighbor_ratio = two           
        
    character(64)::RunName  = "DEMRun" ! run name
    character(64)::ResultsDir= "."     ! result directory 
    character(64)::RestartDir="."      ! restart directory
    integer:: SaveVisu      = 1000     ! save frequency for visulizing file
    integer:: BackupFreq    = 100000   ! save frequency for restarting file
    integer:: Cmd_LFile_Freq= 500      ! report frequency in the terminal 
    integer:: LF_file_lvl   = 5        ! logfile report level      
    integer:: LF_cmdw_lvl   = 3        ! terminal report level
        
    ! Where is the geometry from? 0: Added directly in the program, 1: From DEM.prm, 2: From external STL file
    integer:: GeometrySource =0
    ! If GeometrySource =2, please give a STL file routine, if not, just ignore this variable      
    character(64):: Geom_Dir ="DEMGeom.stl"  
  contains 
    procedure :: ReadDEMOption => DO_ReadDEMOption
  end type DEMS_Options
  type(DEMS_Options),public::  DEM_opt
contains

  !**********************************************************************
  ! DO_ReadDEMOption
  !**********************************************************************
  subroutine DO_ReadDEMOption(this, chFile)
    implicit none
    class(DEMS_Options):: this
    character(*),intent(in)::chFile
           
    ! locals
    logical::RestartFlag
    character(64)::RunName,ResultsDir,RestartDir,Geom_Dir
    real(RK)::Wall_neighbor_ratio,Prtcl_cs_ratio,gravity(3)
    integer::Wall_max_update_iter,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl
    integer::numPrtcl,numPrtclFix,CS_Method,CF_Type,PI_Method,PRI_Method,GeometrySource
    integer::numPrtcl_Type,numWall_type,CntctList_Size,CS_numlvls,Base_wall_id,nUnitFile,ierror
    
    logical::IsPeriodic(3)=.false.
    real(RK)::dtDEM,minpoint(3),maxpoint(3)
    integer:: ifirstDEM,ilastDEM,BackupFreqDEM,SaveVisuDEM
#if defined(CFDDEM) || defined(CFDACM)
    NAMELIST /DEMOptions/ RestartFlag,numPrtcl,numPrtclFix,gravity,CS_Method,CF_Type,PI_Method,PRI_Method,        &
                          numPrtcl_Type,numWall_type,CS_numlvls,CntctList_Size,Base_wall_id,Wall_max_update_iter, &
                          Wall_neighbor_ratio,RunName,ResultsDir,RestartDir,Geom_Dir,Cmd_LFile_Freq,LF_file_lvl,  &
                          LF_cmdw_lvl,GeometrySource,Prtcl_cs_ratio
#else
    NAMELIST /DEMOptions/ RestartFlag,numPrtcl,numPrtclFix,dtDEM,gravity,minpoint,maxpoint,CS_Method,CF_Type,     &
                          PI_Method,PRI_Method,numPrtcl_Type,numWall_type,CS_numlvls,CntctList_Size,Base_wall_id, &
                          Wall_max_update_iter,Wall_neighbor_ratio,RunName,ResultsDir,RestartDir,BackupFreqDEM,   &
                          SaveVisuDEM,Cmd_LFile_Freq,LF_file_lvl,LF_cmdw_lvl,GeometrySource,Geom_Dir,ifirstDEM,   &
                          ilastDEM,Prtcl_cs_ratio,IsPeriodic
#endif

#ifdef CFDDEM
    NAMELIST/CFDDEMCoupling/icouple,UpdateDEMflag,is_clc_Lift,is_clc_Basset,is_clc_Basset_fixed,is_clc_ViscousForce,&
                            is_clc_PressureGradient,is_clc_ViscousForce,is_clc_FluidAcc,FluidAccCoe,SaffmanConst,   &
                            RatioSR,IsAddFluidPressureGradient
    NAMELIST/BassetOptions/ mWinBasset, mTailBasset, BassetAccuracy, BassetTailType 
#endif
#ifdef CFDACM
    integer::icouple_sub
    NAMELIST/CFDACMCoupling/UpdateACMflag,icouple,nForcingExtra,IBM_Scheme,Klub_pp,Klub_pw,Lub_ratio, &
                            Ndt_coll,IsDryColl,St_Crit,IsAddFluidPressureGradient
#endif
              
    open(newunit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) then
       print*, "Cannot open file: "//trim(adjustl(chFile)); STOP
    endif
    read(nUnitFile, nml=DEMOptions)

#ifdef CFDDEM 
    rewind(nUnitFile)
    read(nUnitFile, nml=CFDDEMCoupling)
    rewind(nUnitFile)
    read(nUnitFile, nml=BassetOptions )
    dtDEM     =  dtMax/real(icouple,kind=RK)
    ifirstDEM =  icouple*(ifirst-1)+1
    ilastDEM  =  icouple* ilast
    BackupFreqDEM = icouple* BackupFreq
    SaveVisuDEM   = icouple* SaveVisu
    minpoint= (/zero, zero,zero /)
    maxpoint= (/xlx, yly, zlz/)
    if(BcOption(1)==0)IsPeriodic(1)=.true.
    if(BcOption(3)==0)IsPeriodic(2)=.true.
    if(BcOption(5)==0)IsPeriodic(3)=.true.

    ! Basset history part
    GPrtcl_BassetSeq%HistoryStage= 0
    GPrtcl_BassetSeq%nDataLen    = 1+ (mWinBasset+2) + mTailBasset
    GPrtcl_BassetSeq%iWindowStart= 2
    GPrtcl_BassetSeq%iWindowEnd  = GPrtcl_BassetSeq%iWindowStart+ (mWinBasset+2)-1
    GPrtcl_BassetSeq%iTailStart  = GPrtcl_BassetSeq%iWindowEnd  + 1
    GPrtcl_BassetSeq%iTailEnd    = GPrtcl_BassetSeq%nDataLen
#endif
#ifdef CFDACM
    rewind(nUnitFile)
    read(nUnitFile, nml=CFDACMCoupling)
    if(ischeme==FI_RK2 .and. mod(icouple, 2)/=0 .and. nrank==0) then
      print*,"icouple WRONG!!! icouple must be multiple of  2 for FI_RK2 scheme"; STOP
    endif
    if(ischeme==FI_RK3 .and. mod(icouple,15)/=0 .and. nrank==0) then
      print*,"icouple WRONG!!! icouple must be multiple of 15 for FI_RK3 scheme"; STOP
    endif    
    if(ischeme==FI_AB2) then
      idem_advance_start(1)=1
      idem_advance_end(1)  =icouple
    elseif(ischeme==FI_RK2) then
      icouple_sub= icouple/2
      idem_advance_start(1)=1
      idem_advance_start(2)=icouple_sub+1 
      idem_advance_end(1)  =icouple_sub
      idem_advance_end(2)  =icouple_sub*2
    elseif(ischeme==FI_RK3) then
      icouple_sub= icouple/15
      idem_advance_start(1)=1
      idem_advance_start(2)=icouple_sub*8  + 1 
      idem_advance_start(3)=icouple_sub*10 + 1
      idem_advance_end(1)  =icouple_sub*8
      idem_advance_end(2)  =icouple_sub*10
      idem_advance_end(3)  =icouple_sub*15
    endif
    dtDEM     =  dtMax/real(icouple,kind=RK)
    ifirstDEM =  icouple*(ifirst-1)+1
    ilastDEM  =  icouple* ilast
    BackupFreqDEM = icouple* BackupFreq
    SaveVisuDEM   = icouple* SaveVisu
    minpoint= (/zero, zero,zero /)
    maxpoint= (/xlx, yly, zlz/)
    if(BcOption(1)==0)IsPeriodic(1)=.true.
    if(BcOption(3)==0)IsPeriodic(2)=.true.
    if(BcOption(5)==0)IsPeriodic(3)=.true.
#endif
    close(nUnitFile,IOSTAT=ierror)
           
    this%RestartFlag = RestartFlag
    this%numPrtcl    = numPrtcl
    this%numPrtclFix = numPrtclFix
    this%dt       = dtDEM
    this%ifirst   = ifirstDEM
    this%ilast    = ilastDEM
    this%gravity  = gravity
    this%SimDomain_min = minpoint
    this%SimDomain_max = maxpoint
    this%IsPeriodic = IsPeriodic
           
    this%Prtcl_cs_ratio = Prtcl_cs_ratio
    this%CS_Method = CS_Method
    this%CF_Type   = CF_Type
    this%PI_Method = PI_Method  
    this%PRI_Method= PRI_Method
           
    this%numPrtcl_Type= numPrtcl_Type
    this%numWall_type = numWall_type
           
    this%CntctList_Size= CntctList_Size
    this%CS_numlvls    = CS_numlvls
           
    this%Base_wall_id =  Base_wall_id 
    this%Wall_max_update_iter= Wall_max_update_iter
    this%Wall_neighbor_ratio = Wall_neighbor_ratio 
           
    write(this%RunName,"(A)") RunName
    write(this%ResultsDir,"(A)") ResultsDir
    write(this%RestartDir,"(A)") RestartDir
    this%BackupFreq = BackupFreqDEM
    this%SaveVisu = SaveVisuDEM
    this%Cmd_LFile_Freq = Cmd_LFile_Freq
    this%LF_file_lvl = LF_file_lvl
    this%LF_cmdw_lvl = LF_cmdw_lvl
    this%GeometrySource = GeometrySource
    if( this%GeometrySource==2 ) write(this%Geom_Dir,"(A)")Geom_Dir

  end subroutine DO_ReadDEMOption
end module Prtcl_Parameters
