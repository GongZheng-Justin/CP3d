module Prtcl_DumpPrtcl
  use MPI
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_Parameters
  use m_Decomp2d,only: nrank
  implicit none
  private
#define RKP_Dump 4
  
  character(128)::DumpPrtclDir
  logical::DumpPrtclFlag,ResetDumpFlag
  integer,  dimension(:,:),allocatable::DumpInteMat
  real(RKP_Dump),dimension(:,:),allocatable::DumpRealMat
  integer::nDumpPrtclSize,mDumpPrtclSize,DumpPrtclFreq,iDump
  
  public::DumpPrtclFreq,Initialize_DumpPrtcl,WriteDumpCache,PrtclVarDump
contains
!#define OnlyDumpFpForce
#define nDumpPrtclInte 4

#ifdef OnlyDumpFpForce
#define nDumpPrtclReal 9
#else
#ifdef CFDACM
#define nDumpPrtclReal 21
#else
#define nDumpPrtclReal 9
#endif
#endif

#define Prtcl_Dump_Flag 93

  !**********************************************************************
  ! Initialize_DumpPrtcl
  !**********************************************************************
  subroutine Initialize_DumpPrtcl(chFile)
    implicit none
    character(*),intent(in)::chFile

    ! locals
    real(RK)::yDump
    character(256)::chStr
    integer:: pid,nUnitFile,ierror,ierrTmp
    namelist/DumpPrtclOptions/DumpPrtclFlag,ResetDumpFlag,yDump,mDumpPrtclSize,DumpPrtclDir,DumpPrtclFreq
  
    yDump=0.0_RK
    open(newunit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=ierror)
    if(ierror /= 0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize_DumpPrtcl","Cannot open file:"//trim(adjustl(chFile)))
    read(nUnitFile, nml=DumpPrtclOptions)
    close(nUnitFile,IOSTAT=ierror)

    if(.not. DumpPrtclFlag) return
    if(ResetDumpFlag) then
      do pid=1,GPrtcl_list%nlocal
        if(GPrtcl_PosR(pid)%y >= yDump) then
          GPrtcl_usrMark(pid)=Prtcl_Dump_Flag
        else
          GPrtcl_usrMark(pid)=1
        endif
      enddo
    endif

    nDumpPrtclSize=0; iDump=0; ierror=0
    if(mDumpPrtclSize<10000 .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize_DumpPrtcl","So small mDumpPrtclSize:"//trim(num2str(mDumpPrtclSize)))
    endif
    allocate(DumpInteMat(nDumpPrtclInte,mDumpPrtclSize),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(DumpRealMat(nDumpPrtclReal,mDumpPrtclSize),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize_DumpPrtcl","Allocation failed")

    write(chStr,"(A)") 'mkdir -p '//trim(adjustl(DumpPrtclDir))//' 2> /dev/null'
    if(nrank==0) call system(trim(adjustl(chStr)))

#ifdef OnlyDumpFpForce
    if(nrank==0) call DEMLogInfo%OutInfo("Choose to only dump Fluid-particle force",2)
#else
    if(nrank==0) call DEMLogInfo%OutInfo("Choose to dump full particle information",2)
#endif
  end subroutine Initialize_DumpPrtcl

  !**********************************************************************
  ! WriteDumpCache
  !**********************************************************************
  subroutine WriteDumpCache(itime)
    implicit none
    integer,intent(in)::itime

    ! locals
    integer::k,pid,nlocal

    if(.not. DumpPrtclFlag) return
    nlocal=GPrtcl_list%nlocal
    do pid=1,nlocal
      if(GPrtcl_usrMark(pid)/=Prtcl_Dump_Flag) cycle
      k=1
      nDumpPrtclSize=nDumpPrtclSize+1
      DumpInteMat(1,nDumpPrtclSize)=itime
      DumpInteMat(2,nDumpPrtclSize)=GPrtcl_id(pid)
      DumpInteMat(3,nDumpPrtclSize)=GPrtcl_pType(pid)
      DumpInteMat(4,nDumpPrtclSize)=GPPW_CntctList%IsCntct(pid)
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_PosR(pid)%x, RKP_Dump);       k=k+1 ! 01
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_PosR(pid)%y, RKP_Dump);       k=k+1 ! 02         
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_PosR(pid)%z, RKP_Dump);       k=k+1 ! 03         
#ifdef OnlyDumpFpForce
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%x, RKP_Dump);    k=k+1 ! 04 
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%y, RKP_Dump);    k=k+1 ! 05 
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%z, RKP_Dump);    k=k+1 ! 06       
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%x, RKP_Dump);   k=k+1 ! 07
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%y, RKP_Dump);   k=k+1 ! 08
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%z, RKP_Dump);   k=k+1 ! 09
#else
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_linVel(1,pid)%x, RKP_Dump);   k=k+1 ! 04 
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_linVel(1,pid)%y, RKP_Dump);   k=k+1 ! 05      
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_linVel(1,pid)%z, RKP_Dump);   k=k+1 ! 06      
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_RotVel(1,pid)%x, RKP_Dump);   k=k+1 ! 07     
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_RotVel(1,pid)%y, RKP_Dump);   k=k+1 ! 08    
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_RotVel(1,pid)%z, RKP_Dump);   k=k+1 ! 09    
#ifdef CFDACM     
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%x, RKP_Dump);    k=k+1 ! 10     
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%y, RKP_Dump);    k=k+1 ! 11 
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%z, RKP_Dump);    k=k+1 ! 12       
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%x, RKP_Dump);   k=k+1 ! 13
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%y, RKP_Dump);   k=k+1 ! 14
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%z, RKP_Dump);   k=k+1 ! 15
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_CntctForce(pid)%x, RKP_Dump); k=k+1 ! 16
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_CntctForce(pid)%y, RKP_Dump); k=k+1 ! 17
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_CntctForce(pid)%z, RKP_Dump); k=k+1 ! 18
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_Torque(pid)%x, RKP_Dump);     k=k+1 ! 19
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_Torque(pid)%y, RKP_Dump);     k=k+1 ! 20
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_Torque(pid)%z, RKP_Dump);     k=k+1 ! 21
#endif
#endif
      if(nDumpPrtclSize==mDumpPrtclSize) then
        call PrtclVarDump(itime)
        iDump=iDump+1
      endif
    enddo
    iDump=0
  end subroutine WriteDumpCache

  !**********************************************************************
  ! PrtclVarDump
  !**********************************************************************
  subroutine PrtclVarDump(itime)
    implicit none
    integer,intent(in)::itime

    ! locals
    integer::ierror,nUnit
    character(128)::chFile

    if(.not.DumpPrtclFlag)return
    if(nDumpPrtclSize==0) return
    write(chFile,'(A,I5.5,A,I10.10,A,I2.2)')trim(DumpPrtclDir)//'rank',nrank,'_',itime/icouple,'_',iDump

    open(newunit=nUnit,file=trim(chFile),status='replace',form='unformatted',access='stream',IOSTAT=ierror)
    IF(ierror/=0) THEN
      call DEMLogInfo%CheckForError(ErrT_Pass,"PrtclVarDump","Cannot open file: "//trim(chFile))
    ELSE
      write(nUnit)nDumpPrtclSize ! Added by Zheng Gong, 2023-05-04
      write(nUnit)DumpInteMat(:,1:nDumpPrtclSize)
      write(nUnit)DumpRealMat(:,1:nDumpPrtclSize)
    ENDIF
    close(nUnit,IOSTAT=ierror)

    nDumpPrtclSize=0
  end subroutine PrtclVarDump
end module Prtcl_DumpPrtcl

#undef nDumpPrtclInte
#undef nDumpPrtclReal
#undef Prtcl_Dump_Flag

#ifdef OnlyDumpFpForce
#undef OnlyDumpFpForce
#endif

#undef RKP_Dump
