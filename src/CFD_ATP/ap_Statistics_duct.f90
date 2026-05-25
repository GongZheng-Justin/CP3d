#include "definitions_inc.f90"
module ap_Statistics
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d !,only:nrank
#ifdef CFDSecondOrder
  use f2_MeshAndMetries !,only:yp
  use f2_Parameters
#else
  use f4_MeshAndMetries !,only:yp
  use f4_Parameters
#endif
  use ap_Property
  use ap_Decomp_2d
  use ap_Variables
  use ap_Parameters
  implicit none
  private

  integer::nslab,nShannon
  integer::npType,npstime
  real(RK),dimension(:),allocatable:: ypForPs  ! y point for particle statistics 
  integer,dimension(:,:,:),allocatable::npsum, npsumR

  public::InitATPStatistics,ClcATPStatistics
contains
  
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
!   InitATPStatistics
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitATPStatistics(chFile); implicit none
    character(len=*),intent(in)::chFile
    
    ! locals
    integer:: iUnit,ierr,j,k
    integer:: iErr01,iErr02,iErr03,iErr04,iErr05,iErrSum
    NAMELIST/ParticleStatisticOption/nslab,nShannon

    open(newunit=iUnit, file=chFile,status='old',form='formatted',IOSTAT=ierr)
    if(ierr/=0) call ATPLogInfo%CheckForError(ErrT_Abort,"InitATPStatistics", "Cannot open file: "//strip(chFile))
    read(iUnit, nml=ParticleStatisticOption)
    if(nrank==0)write(ATPLogInfo%iUnit, nml=ParticleStatisticOption)
    close(iUnit,IOSTAT=ierr)
    if(nrank==0 .and. mod(nyc,nslab)/=0) call ATPLogInfo%CheckForError(ErrT_Abort,"InitATPStatistics","mod(nyc,nslab)/=0")
    if(nrank==0 .and. (nShannon>nyc .or. nShannon<0)) call ATPLogInfo%CheckForError(ErrT_Abort,"InitATPStatistics","nShannon wrong")    
    
    npType= ATP_opt%numPrtcl_Type
    allocate(ypForPs(nslab+1),      npsum(nslab,nzc,npType),  npsumR(nslab,nzc,npType),  Stat=iErr01)
    iErrSum=abs(iErr01)
    if(iErrSum/=0) call ATPLogInfo%CheckForError(ErrT_Abort,"ATP_InitStat ","Allocation failed")
    
    k=nyc/nslab
    do j=0,nslab
      ypForPs(j+1)=yp(j*k+1)
    enddo
    call ResetStatVar()
    
    if(nrank/=0) return
  end subroutine InitATPStatistics

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ResetStatVar
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ResetStatVar(); implicit none

    npstime=0;       npsum=0
  end subroutine ResetStatVar

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ClcATPStatistics
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ClcATPStatistics(); implicit none

    ! locals
    type(real3)::pos
    character(len=:),allocatable::filename
    integer::pid,nlocal,itype,js,je,jc,ierr,iUnit,kc
    real(RK),allocatable,dimension(:,:,:)::sumTmp
    real(RK)::inpstime
    
    allocate(sumTmp(npType,nslab,nzc))
    nlocal= GPrtcl_list%nlocal
    do pid=1,nlocal
      pos  = GPrtcl_PosR(pid)
      itype= GPrtcl_pType(pid)
      kc= floor(pos%z*rdz)+1; kc=min(kc,y1end(3)); kc=max(kc,y1start(3));
      
      ! if pos%y is within [0,yly), jc will be within [1,nslab]
      js=0
      je=nslab+2
      do
        jc=(js+je)/2
        if(je-js==1) exit
        if(pos%y< ypForPs(jc)) then
          je =jc
        else
          js =jc
        endif
      enddo
      npsum(jc,kc,itype)     = npsum(jc,kc,itype) + 1
    enddo
    call MPI_REDUCE(npsum,npsumR, npType*nslab*nzc,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    npstime = npstime + 1
    if(mod(itime,SaveStat)/=0) return
    if(nrank==0) then
      inpstime = 1.0_RK/real(npstime,RK)
      sumTmp=inpstime*real(npsumR,RK)
      do itype=1,npType
        filename = strip(ATP_Opt%ResultsDir) // 'pstats' // int2str(itime,10) // '_set' // int2str(itype,4)
        open(newunit=iUnit,file=filename,status='replace',form='formatted',IOSTAT=ierr)
        if(ierr /= 0) then
          call ATPLogInfo%CheckForError(ErrT_Pass,"ClcATPStatistics","Cannot open file: "//filename)
        else              
        Block 
          character(len=128)::FormatStr
          write(FormatStr,'(A,I3,A)')'(',nslab,'ES24.15)'
          do kc=1,nzc
            write(iUnit,FormatStr) sumTmp(1:nslab,kc,itype)
          enddo
        End Block
        endif
        close(iUnit,IOSTAT=ierr)
      enddo
    endif

    call ResetStatVar()
  end subroutine ClcATPStatistics

end module ap_Statistics
