module Prtcl_DumpPrtcl
  use MPI
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_Parameters
  use m_Decomp2d,only: nrank
  implicit none
  private
  
  integer,parameter::RKP=4
  character(128)::DumpPrtclDir
  logical::DumpPrtclFlag,ResetDumpFlag
  integer,  dimension(:,:),allocatable::DumpInteMat
  real(RKP),dimension(:,:),allocatable::DumpRealMat
  integer::nDumpPrtclSize,mDumpPrtclSize,WriteCacheFreq,iDump
  
  type Prtcl_Dump
    integer:: WriteCacheFreq=100000
  contains
    procedure:: Initialize
    procedure:: WriteCache   
    procedure:: PrtclVarDump   
  end type Prtcl_Dump
  type(Prtcl_Dump),public:: DEM_PDump
contains

#define nDumpPrtclInte 4
#ifdef CFDACM
#define nDumpPrtclReal 21
#else
#define nDumpPrtclReal 9
#endif

#define Prtcl_Dump_Flag 93

  !**********************************************************************
  ! Initialize
  !**********************************************************************
  subroutine Initialize(this,chFile)
    implicit none
    class(Prtcl_Dump)::this
    character(*),intent(in)::chFile

    ! locals
    real(RK)::yDump=zero
    character(256)::chStr
    integer:: pid,nUnitFile,ierror,ierrTmp
    namelist/DumpPrtclOptions/DumpPrtclFlag,ResetDumpFlag,yDump,mDumpPrtclSize,DumpPrtclDir,WriteCacheFreq
  
    open(newunit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=ierror)
    if(ierror /= 0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitDumpPrtclVar","Cannot open file:"//trim(adjustl(chFile)))
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
    DEM_PDump%WriteCacheFreq =  WriteCacheFreq

    nDumpPrtclSize=0; iDump=0; ierror=0
    if(mDumpPrtclSize<10000 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitDumpPrtclVar","So small mDumpPrtclSize:"//trim(num2str(mDumpPrtclSize)) )
    allocate(DumpInteMat(nDumpPrtclInte,mDumpPrtclSize),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(DumpRealMat(nDumpPrtclReal,mDumpPrtclSize),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitDumpPrtclVar","Allocation failed")

    write(chStr,"(A)") 'mkdir -p '//trim(adjustl(DumpPrtclDir))//' 2> /dev/null'
    if(nrank==0) call system(trim(adjustl(chStr)))
  end subroutine Initialize

  !**********************************************************************
  ! WriteCache
  !**********************************************************************
  subroutine WriteCache(this,itime)
    implicit none
    class(Prtcl_Dump)::this
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
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_PosR(pid)%x, RKP);       k=k+1 ! 01
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_PosR(pid)%y, RKP);       k=k+1 ! 02         
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_PosR(pid)%z, RKP);       k=k+1 ! 03         
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_linVel(1,pid)%x, RKP);   k=k+1 ! 04 
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_linVel(1,pid)%y, RKP);   k=k+1 ! 05      
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_linVel(1,pid)%z, RKP);   k=k+1 ! 06      
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_RotVel(1,pid)%x, RKP);   k=k+1 ! 07     
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_RotVel(1,pid)%y, RKP);   k=k+1 ! 08    
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_RotVel(1,pid)%z, RKP);   k=k+1 ! 09    
#ifdef CFDACM     
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%x, RKP);    k=k+1 ! 10     
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%y, RKP);    k=k+1 ! 11 
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpForce(pid)%z, RKP);    k=k+1 ! 12       
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_CntctForce(pid)%x, RKP); k=k+1 ! 13
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_CntctForce(pid)%y, RKP); k=k+1 ! 14
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_CntctForce(pid)%z, RKP); k=k+1 ! 15
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%x, RKP);   k=k+1 ! 16
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%y, RKP);   k=k+1 ! 17
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_FpTorque(pid)%z, RKP);   k=k+1 ! 18
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_Torque(pid)%x, RKP);     k=k+1 ! 19
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_Torque(pid)%y, RKP);     k=k+1 ! 20
      DumpRealMat(k,nDumpPrtclSize)=real(GPrtcl_Torque(pid)%z, RKP);     k=k+1 ! 21
#endif
      if(nDumpPrtclSize==mDumpPrtclSize) then
        call this%PrtclVarDump(itime)
        iDump=iDump+1
      endif
    enddo
    iDump=0
  end subroutine WriteCache

  !**********************************************************************
  ! PrtclVarDump
  !**********************************************************************
  subroutine PrtclVarDump(this,itime)
    implicit none
    class(Prtcl_Dump)::this
    integer,intent(in)::itime

    ! locals
    integer::ierror,nUnit
    character(128)::chFile

    if(.not.DumpPrtclFlag)return
    if(nDumpPrtclSize==0) return
    write(chFile,'(A,I5.5,A,I10.10,A,I2.2)')trim(DumpPrtclDir)//'rank',nrank,'_',itime,'_',iDump

    open(newunit=nUnit,file=trim(chFile),status='replace',form='unformatted',access='stream',IOSTAT=ierror)
    IF(ierror/=0) THEN
      call DEMLogInfo%CheckForError(ErrT_Pass,"PrtclVarDump","Cannot open file: "//trim(chFile))
    ELSE
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
