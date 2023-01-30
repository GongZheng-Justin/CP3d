module LPT_Statistics
  use MPI
  use m_TypeDef
  use m_LogInfo
  use LPT_Property
  use LPT_Decomp_2d
  use LPT_Variables
  use LPT_Parameters
  use m_Decomp2d,only:nrank
  use m_MeshAndMetries,only:yp
  use m_Parameters,only: ivstats,saveStat,nyc,itime,yly,ilast 
  implicit none
  private

  integer::nslab,nShannon
  integer::npType,npstime
  real(RK),dimension(:),allocatable:: ypForPs  ! y point for particle statistics 
  integer,dimension(:,:),allocatable::npsum
  real(RK),dimension(:,:),allocatable::upsum,vpsum,wpsum,upupsum,vpvpsum,wpwpsum,upvpsum
  real(RK),dimension(:,:),allocatable::ufsum,vfsum,wfsum,ufufsum,vfvfsum,wfwfsum,ufvfsum

  public::InitLPTStatistics,ClcLPTStatistics
contains
  
!********************************************************************************
!   InitLPTStatistics
!********************************************************************************
  subroutine InitLPTStatistics(chFile)
    implicit none
    character(*),intent(in)::chFile
    
    ! locals
    character(len=128)::filename
    integer:: nUnit,ierror,j,k
    integer:: iErr01,iErr02,iErr03,iErr04,iErr05,iErrSum
    NAMELIST/ParticleStatisticOption/nslab,nShannon

    open(newunit=nUnit, file=chFile,status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"InitLPTStatistics", "Cannot open file: "//trim(chFile))
    read(nUnit, nml=ParticleStatisticOption)
    if(nrank==0)write(LPTLogInfo%nUnit, nml=ParticleStatisticOption)
    close(nUnit,IOSTAT=ierror)
    if(nrank==0 .and. mod(nyc,nslab)/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"InitLPTStatistics","mod(nyc,nslab)/=0")
    if(nrank==0 .and. (nShannon>nyc .or. nShannon<0)) call LPTLogInfo%CheckForError(ErrT_Abort,"InitLPTStatistics","nShannon wrong")    
    
    npType= LPT_opt%numPrtcl_Type
    allocate(ypForPs(nslab+1),      npsum(npType,nslab),   Stat=iErr01)
    allocate(upsum(npType,nslab),   vpsum(npType,nslab),   wpsum(npType,nslab),   Stat=iErr02)
    allocate(ufsum(npType,nslab),   vfsum(npType,nslab),   wfsum(npType,nslab),   Stat=iErr03)
    allocate(upupsum(npType,nslab), vpvpsum(npType,nslab), wpwpsum(npType,nslab), upvpsum(npType,nslab), Stat=iErr04)
    allocate(ufufsum(npType,nslab), vfvfsum(npType,nslab), wfwfsum(npType,nslab), ufvfsum(npType,nslab), Stat=iErr05)
    iErrSum=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)+abs(iErr05)
    if(iErrSum/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"LPT_InitStat ","Allocation failed")
    
    k=nyc/nslab
    do j=0,nslab
      ypForPs(j+1)=yp(j*k+1)
    enddo
    call ResetStatVar()
    
    if(nrank/=0) return
    write(filename,'(A,I10.10)') trim(LPT_Opt%ResultsDir)//'shannon',ilast
    open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror /= 0) call LPTLogInfo%CheckForError(ErrT_Abort,"InitLPTStatistics","Cannot open file: "//trim(filename))
    close(nUnit,IOSTAT=ierror)    
  end subroutine InitLPTStatistics

  !******************************************************************
  ! ResetStatVar
  !******************************************************************
  subroutine ResetStatVar()
    implicit none

    npstime=0;       npsum=0
    upsum=0.0_RK;      vpsum=0.0_RK;      wpsum=0.0_RK    
    ufsum=0.0_RK;      vfsum=0.0_RK;      wfsum=0.0_RK
    upupsum=0.0_RK;    vpvpsum=0.0_RK;    wpwpsum=0.0_RK;    upvpsum=0.0_RK;
    ufufsum=0.0_RK;    vfvfsum=0.0_RK;    wfwfsum=0.0_RK;    ufvfsum=0.0_RK;
  end subroutine ResetStatVar

  !********************************************************************************
  ! ClcLPTStatistics
  !********************************************************************************  
  subroutine ClcLPTStatistics()
    implicit none

    ! locals
    character(len=128)::filename
    type(real3)::pos,VPrtcl,VFluid
    real(RK)::irnpsum,inpstime,dyShannon,pse
    real(RK),dimension(npType)::shannon_entropy
    integer::pid,nlocal,itype,js,je,jc,ierror,nUnit
    real(RK),dimension(14,npType,nslab)::sumStat
    real(RK),dimension(:,:,:),allocatable::sumStatR
    integer,dimension(:,:),allocatable::npsslot,npsslotR
    
    allocate(npsslot(npType,nShannon),npsslotR(npType,nShannon))
    sumStat=0.0_RK; npsslot=0
    nlocal= GPrtcl_list%nlocal
    dyShannon=yly/real(nShannon,RK)
    do pid=1,nlocal
      pos  = GPrtcl_posR(pid)
      itype= GPrtcl_pType(pid)
      VPrtcl=GPrtcl_linVel(1,pid)
      VFluid=GPrtcl_VFluid(pid)

      ! Shannon entropy
      jc=int(pos%y/dyShannon)+1
      npsslot(itype,jc)=npsslot(itype,jc)+1
      
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
      npsum(itype,jc)     = npsum(itype,jc) + 1
      sumStat(1,itype,jc) = sumStat(1,itype,jc) + VPrtcl%x
      sumStat(2,itype,jc) = sumStat(2,itype,jc) + VPrtcl%y
      sumStat(3,itype,jc) = sumStat(3,itype,jc) + VPrtcl%z
      sumStat(4,itype,jc) = sumStat(4,itype,jc) + VPrtcl%x *VPrtcl%x
      sumStat(5,itype,jc) = sumStat(5,itype,jc) + VPrtcl%y *VPrtcl%y
      sumStat(6,itype,jc) = sumStat(6,itype,jc) + VPrtcl%z *VPrtcl%z
      sumStat(7,itype,jc) = sumStat(7,itype,jc) + VPrtcl%x *VPrtcl%y
      sumStat(8,itype,jc) = sumStat(8,itype,jc) + VFluid%x
      sumStat(9,itype,jc) = sumStat(9,itype,jc) + VFluid%y
      sumStat(10,itype,jc)= sumStat(10,itype,jc)+ VFluid%z
      sumStat(11,itype,jc)= sumStat(11,itype,jc)+ VFluid%x *VFluid%x
      sumStat(12,itype,jc)= sumStat(12,itype,jc)+ VFluid%y *VFluid%y
      sumStat(13,itype,jc)= sumStat(13,itype,jc)+ VFluid%z *VFluid%z
      sumStat(14,itype,jc)= sumStat(14,itype,jc)+ VFluid%x *VFluid%y
    enddo
    do jc=1,nslab
      do itype=1,npType
        upsum(itype,jc)  = upsum(itype,jc)  + sumStat( 1,itype,jc)
        vpsum(itype,jc)  = vpsum(itype,jc)  + sumStat( 2,itype,jc)
        wpsum(itype,jc)  = wpsum(itype,jc)  + sumStat( 3,itype,jc)
        upupsum(itype,jc)= upupsum(itype,jc)+ sumStat( 4,itype,jc)
        vpvpsum(itype,jc)= vpvpsum(itype,jc)+ sumStat( 5,itype,jc)
        wpwpsum(itype,jc)= wpwpsum(itype,jc)+ sumStat( 6,itype,jc)
        upvpsum(itype,jc)= upvpsum(itype,jc)+ sumStat( 7,itype,jc)
        ufsum(itype,jc)  = ufsum(itype,jc)  + sumStat( 8,itype,jc)
        vfsum(itype,jc)  = vfsum(itype,jc)  + sumStat( 9,itype,jc)
        wfsum(itype,jc)  = wfsum(itype,jc)  + sumStat(10,itype,jc)
        ufufsum(itype,jc)= ufufsum(itype,jc)+ sumStat(11,itype,jc)
        vfvfsum(itype,jc)= vfvfsum(itype,jc)+ sumStat(12,itype,jc)
        wfwfsum(itype,jc)= wfwfsum(itype,jc)+ sumStat(13,itype,jc)
        ufvfsum(itype,jc)= ufvfsum(itype,jc)+ sumStat(14,itype,jc)
      enddo
    enddo
    call MPI_REDUCE(npsslot,npsslotR, npType*nShannon,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
       
    ! Shannon entropy part
    shannon_entropy = 0.0_RK    
    IF(nrank==0) THEN
      DO itype=1,npType
        nlocal=LPTProperty%nPrtcl_in_Bin(itype)
        do jc=1,nShannon
          pse = real(npsslotR(itype,jc),RK)/real(nlocal,RK)
          if(pse > 0.0_RK) then
            shannon_entropy(itype) = shannon_entropy(itype) -pse*log(pse)
          endif
        enddo
        shannon_entropy(itype) = shannon_entropy(itype)/log(real(nShannon,RK))
      ENDDO
      
      write(filename,'(A,I10.10)') trim(LPT_Opt%ResultsDir)//'shannon',ilast
      open(newunit=nUnit,file=filename,status='old',position='append',form='formatted',IOSTAT=ierror)
      if(ierror /= 0) then
        call LPTLogInfo%CheckForError(ErrT_Pass,"ClcLPTStatistics","Cannot open file: "//trim(filename))
      else  
        write(nUnit,*)itime,shannon_entropy
      endif
      close(nUnit,IOSTAT=ierror)
    ENDIF
    deallocate(npsslot,npsslotR)
    
    npstime = npstime + 1
    if(mod(itime,SaveStat)/=0) return
    
    sumStat( 1,:,:)=upsum;    sumStat( 2,:,:)=vpsum;    sumStat( 3,:,:)=wpsum
    sumStat( 4,:,:)=upupsum;  sumStat( 5,:,:)=vpvpsum;  sumStat( 6,:,:)=wpwpsum;  sumStat( 7,:,:)=upvpsum
    sumStat( 8,:,:)=ufsum;    sumStat( 9,:,:)=vfsum;    sumStat(10,:,:)=wfsum
    sumStat(11,:,:)=ufufsum;  sumStat(12,:,:)=vfvfsum;  sumStat(13,:,:)=wfwfsum;  sumStat(14,:,:)=ufvfsum
    allocate(npsslot(npType,nslab),sumStatR(14,npType,nslab))
    call MPI_REDUCE(npsum,  npsslot,    npType*nslab, int_type, MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(sumStat,sumStatR,14*npType*nslab, real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)

    if(nrank==0) then
      inpstime = 1.0_RK/real(npstime,RK)
      do itype=1,npType
        write(filename,"(A,I10.10,A,I4.4)") trim(LPT_Opt%ResultsDir)//'pstats',itime,'_set',itype
        open(newunit=nUnit,file=trim(filename),status='replace',form='formatted',IOSTAT=ierror)
        if(ierror /= 0) then
          call LPTLogInfo%CheckForError(ErrT_Pass,"ClcLPTStatistics","Cannot open file: "//trim(filename))
        else              
          write(nUnit,'(A,I7,A,I7,A,I7)')'    The time step range for this particle statistics is ', &
                                       itime-(npstime-1)*ivstats, ':', ivstats, ':', itime
          write(nUnit,*)
          write(nUnit,'(A)')'  yp, np, up, vp, wp, upup, vpvp, wpwp, upvp, uf, vf, wf, ufuf, vfvf, wfwf, ufvf' 
          do jc=1,nslab
            if(npsslot(itype,jc)==0) then
              write(nUnit,'(20ES24.15)')(ypForPs(jc)+ypForPs(jc+1))*0.5_RK,   & ! 1
                                     0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,   &
                                     0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK
            else
              irnpsum=1.0_RK/real(npsslot(itype,jc),RK)
              write(nUnit,'(20ES24.15)') (ypForPs(jc)+ypForPs(jc+1))*0.5_RK, & ! 1
                                        real(npsslot(itype,jc),RK)*inpstime, & ! 2
                                              sumStatR( 1,itype,jc)*irnpsum, & ! 3
                                              sumStatR( 2,itype,jc)*irnpsum, & ! 4
                                              sumStatR( 3,itype,jc)*irnpsum, & ! 5
                                              sumStatR( 4,itype,jc)*irnpsum, & ! 6
                                              sumStatR( 5,itype,jc)*irnpsum, & ! 7
                                              sumStatR( 6,itype,jc)*irnpsum, & ! 8
                                              sumStatR( 7,itype,jc)*irnpsum, & ! 9
                                              sumStatR( 8,itype,jc)*irnpsum, & ! 10
                                              sumStatR( 9,itype,jc)*irnpsum, & ! 11
                                              sumStatR(10,itype,jc)*irnpsum, & ! 12
                                              sumStatR(11,itype,jc)*irnpsum, & ! 13
                                              sumStatR(12,itype,jc)*irnpsum, & ! 14
                                              sumStatR(13,itype,jc)*irnpsum, & ! 15
                                              sumStatR(14,itype,jc)*irnpsum    ! 16
            endif
          enddo
        endif
        close(nUnit,IOSTAT=ierror)
      enddo
    endif

    call ResetStatVar()
  end subroutine ClcLPTStatistics

end module LPT_Statistics
