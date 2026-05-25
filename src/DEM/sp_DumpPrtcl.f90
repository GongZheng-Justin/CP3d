#include "definitions_inc.f90"
module sp_DumpPrtcl
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d,only: nrank
  use mc_FileOperator,only:mkdir
  use sp_Variables
  use sp_CL_and_CF
  use sp_Parameters
  implicit none
  private
#define RKP_Dump 4
  
  logical::DumpPrtclFlag,ResetDumpFlag
  character(len=:),allocatable::DumpPrtclDir_
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

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Initialize_DumpPrtcl
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Initialize_DumpPrtcl(chFile); implicit none
    character(len=*),intent(in)::chFile

    ! locals
    real(RK)::yDump
    character(len=512)::DumpPrtclDir
    integer:: pid,iUnit,ierr,ierrTmp
    namelist/DumpPrtclOptions/DumpPrtclFlag,ResetDumpFlag,yDump,mDumpPrtclSize,DumpPrtclDir,DumpPrtclFreq
  
    yDump=0.0_RK
    open(newunit=iUnit, file=chFile, status='old', form='formatted', IOSTAT=ierr)
    if(ierr /= 0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize_DumpPrtcl","Cannot open file:"//strip(chFile))
    read(iUnit, nml=DumpPrtclOptions)
    close(iUnit,IOSTAT=ierr)
    DumpPrtclDir_ = strip(DumpPrtclDir)

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

    nDumpPrtclSize=0; iDump=0; ierr=0
    if(mDumpPrtclSize<10000 .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize_DumpPrtcl","So small mDumpPrtclSize:"//strip(num2str(mDumpPrtclSize)))
    endif
    allocate(DumpInteMat(nDumpPrtclInte,mDumpPrtclSize),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
    allocate(DumpRealMat(nDumpPrtclReal,mDumpPrtclSize),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
    if(ierr/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize_DumpPrtcl","Allocation failed")

    if(nrank==0) then
      call mkdir(DumpPrtclDir_, ierr)
      if(ierr < 0) call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize_DumpPrtcl","DumpPrtclDir failed")
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

#ifdef OnlyDumpFpForce
    if(nrank==0) call DEMLogInfo%OutInfo("Choose to only dump Fluid-particle force",2)
#else
    if(nrank==0) call DEMLogInfo%OutInfo("Choose to dump full particle information",2)
#endif
  end subroutine Initialize_DumpPrtcl

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! WriteDumpCache
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine WriteDumpCache(itime); implicit none
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
      DumpInteMat(4,nDumpPrtclSize)=IsPrtclContact(pid)
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

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PrtclVarDump
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PrtclVarDump(itime); implicit none
    integer,intent(in)::itime

    ! locals
    integer::ierr,iUnit
    character(len=:),allocatable::chFile

    if(.not.DumpPrtclFlag)return
    if(nDumpPrtclSize==0) return
    chFile = strip(DumpPrtclDir_) // 'rank' // int2str(nrank,5) // '_' // int2str(itime/icouple,10) // '_' // int2str(iDump,2)

    open(newunit=iUnit,file=chFile,status='replace',form='unformatted',access='stream',IOSTAT=ierr)
    IF(ierr/=0) THEN
      call DEMLogInfo%CheckForError(ErrT_Pass,"PrtclVarDump","Cannot open file: "//chFile)
    ELSE
      write(iUnit)nDumpPrtclSize ! Added by Zheng Gong, 2023-05-04
      write(iUnit)DumpInteMat(:,1:nDumpPrtclSize)
      write(iUnit)DumpRealMat(:,1:nDumpPrtclSize)
    ENDIF
    close(iUnit,IOSTAT=ierr)

    nDumpPrtclSize=0
  end subroutine PrtclVarDump
end module sp_DumpPrtcl

#undef nDumpPrtclInte
#undef nDumpPrtclReal
#undef Prtcl_Dump_Flag

#ifdef OnlyDumpFpForce
#undef OnlyDumpFpForce
#endif

#undef RKP_Dump
