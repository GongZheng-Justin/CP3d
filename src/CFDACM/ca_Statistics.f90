module ca_Statistics
  use MPI
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_Decomp_2d
  use Prtcl_Parameters
  use m_Decomp2d,only:nrank
  use m_MeshAndMetries,only:yp
  use m_Parameters,only: ivstats,saveStat,nyc,itime,xlx,zlz,ilast
  implicit none
  private

  integer::npType,nslab,npstime
  integer,dimension(:,:),allocatable::npsum
  real(RK),dimension(:),allocatable:: ypForPs  ! y point for particle statistics
  real(RK),dimension(:,:),allocatable::wxsum,wysum,wzsum,wxwxsum,wywysum,wzwzsum
  real(RK),dimension(:,:),allocatable::upsum,vpsum,wpsum,upupsum,vpvpsum,wpwpsum,upvpsum

  public::InitCAStatistics,ClcCAStatistics
#define ClcTransportRate
contains
  
  !********************************************************************************
  !   InitCAStatistics
  !********************************************************************************
  subroutine InitCAStatistics(chFile)
    implicit none
    character(*),intent(in)::chFile
        
    ! locals
    character(len=128)::filename
    NAMELIST/ParticleStatisticOption/nslab
    integer::j,k,iErr01,iErr02,iErr03,ierror,nUnit

    npType= DEM_Opt%numPrtcl_Type
    open(newunit=nUnit, file=chFile,status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitDEMStatistics", "Cannot open file: "//trim(chFile))
    read(nUnit, nml=ParticleStatisticOption)
    if(nrank==0)write(DEMLogInfo%nUnit, nml=ParticleStatisticOption)
    close(nUnit,IOSTAT=ierror)
    if(nrank==0 .and. mod(nyc,nslab)/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitCAStatistics","mod(nyc,nslab)/=0")
    
    allocate(ypForPs(nslab+1),Stat=iErr01)
    allocate(npsum(npType,nslab), upsum(npType,nslab),   vpsum(npType,nslab),   wpsum(npType,nslab),   Stat=iErr02)
    allocate(upupsum(npType,nslab), vpvpsum(npType,nslab), wpwpsum(npType,nslab),upvpsum(npType,nslab),wxsum(npType,nslab), &
             wysum(npType,nslab),wzsum(npType,nslab),wxwxsum(npType,nslab),wywysum(npType,nslab),wzwzsum(npType,nslab),Stat=iErr03)
    ierror=abs(iErr01)+abs(iErr02)+abs(iErr03)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitCAStatistics ","Allocation failed")
    k=nyc/nslab
    do j=0,nslab
      ypForPs(j+1)=yp(j*k+1)
    enddo
    call ResetStatVar()

#ifdef ClcTransportRate
    if(nrank==0) then
      write(filename,'(A,I10.10)')trim(DEM_opt%ResultsDir)//"TransportRate",ilast
      open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
      if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitCAStatistics","Cannot open file: "//trim(filename))
      close(nUnit,IOSTAT=ierror)
    endif
#endif
  end subroutine InitCAStatistics

  !******************************************************************
  ! ResetStatVar
  !******************************************************************
  subroutine ResetStatVar()
    implicit none

    npstime=0;       npsum=0
    upsum=zero;      vpsum=zero;      wpsum=zero
    upupsum=zero;    vpvpsum=zero;    wpwpsum=zero;    upvpsum=zero
    wxsum=zero;      wysum=zero;      wzsum=zero
    wxwxsum=zero;    wywysum=zero;    wzwzsum=zero
  end subroutine ResetStatVar

  !********************************************************************************
  ! ClcCAStatistics
  !********************************************************************************  
  subroutine ClcCAStatistics()
    implicit none

    ! locals
    character(len=128)::filename
    type(real3)::pos,VPrtcl,OmePrtcl
    integer,dimension(npType,nslab)::npsslot
    integer::pid,nlocal,itype,js,je,jc,ierror,nUnit
    real(RK)::irnpsum,inpstime,SumVel,SumVelR
    real(RK),dimension(13,npType,nslab)::sumStat,sumStatR

    sumStat=0
    nlocal= GPrtcl_list%nlocal
    do pid=1,nlocal
      pos  = GPrtcl_posR(pid)
      itype= GPrtcl_pType(pid)
      VPrtcl=GPrtcl_linVel(1,pid)
      OmePrtcl=GPrtcl_rotVel(1,pid)
      
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
      sumStat( 1,itype,jc) = sumStat( 1,itype,jc) + VPrtcl%x
      sumStat( 2,itype,jc) = sumStat( 2,itype,jc) + VPrtcl%y
      sumStat( 3,itype,jc) = sumStat( 3,itype,jc) + VPrtcl%z
      sumStat( 4,itype,jc) = sumStat( 4,itype,jc) + VPrtcl%x *VPrtcl%x
      sumStat( 5,itype,jc) = sumStat( 5,itype,jc) + VPrtcl%y *VPrtcl%y
      sumStat( 6,itype,jc) = sumStat( 6,itype,jc) + VPrtcl%z *VPrtcl%z
      sumStat( 7,itype,jc) = sumStat( 7,itype,jc) + VPrtcl%x *VPrtcl%y
      sumStat( 8,itype,jc) = sumStat( 8,itype,jc) + OmePrtcl%x
      sumStat( 9,itype,jc) = sumStat( 9,itype,jc) + OmePrtcl%y
      sumStat(10,itype,jc) = sumStat(10,itype,jc) + OmePrtcl%z
      sumStat(11,itype,jc) = sumStat(11,itype,jc) + OmePrtcl%x *OmePrtcl%x    
      sumStat(12,itype,jc) = sumStat(12,itype,jc) + OmePrtcl%y *OmePrtcl%y 
      sumStat(13,itype,jc) = sumStat(13,itype,jc) + OmePrtcl%z *OmePrtcl%z 
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
        wxsum(itype,jc)  = wxsum(itype,jc)  + sumStat( 8,itype,jc)
        wysum(itype,jc)  = wysum(itype,jc)  + sumStat( 9,itype,jc)
        wzsum(itype,jc)  = wzsum(itype,jc)  + sumStat(10,itype,jc)
        wxwxsum(itype,jc)= wxwxsum(itype,jc)+ sumStat(11,itype,jc)
        wywysum(itype,jc)= wywysum(itype,jc)+ sumStat(12,itype,jc)
        wzwzsum(itype,jc)= wzwzsum(itype,jc)+ sumStat(13,itype,jc)
      enddo
    enddo

#ifdef ClcTransportRate
    SumVel=zero
    do pid=1,nlocal
      itype = GPrtcl_pType(pid)
      SumVel= SumVel+ GPrtcl_linVel(1,pid)%x* DEMProperty%Prtcl_PureProp(itype)%Volume
    enddo
    call MPI_REDUCE(SumVel,SumVelR,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      write(filename,'(A,I10.10)')trim(DEM_opt%ResultsDir)//"TransportRate",ilast
      open(newunit=nUnit,file=filename,status='old',form='formatted',position='append',IOSTAT=ierror)
      if(ierror /= 0) then
        call DEMLogInfo%CheckForError(ErrT_Pass,"ClcLCAStatistics: ","Cannot open file: "//trim(filename))
      else
        write(nUnit,'(I7,ES24.15)')itime,SumVelR/(xlx*zlz)
        close(nUnit,IOSTAT=ierror)
      endif
    endif
#endif

    npstime = npstime + 1
    if(mod(itime,SaveStat)/=0) return
    sumStat( 1,:,:)=upsum;    sumStat( 2,:,:)=vpsum;    sumStat(3,:,:)=wpsum
    sumStat( 4,:,:)=upupsum;  sumStat( 5,:,:)=vpvpsum;  sumStat(6,:,:)=wpwpsum;  sumStat(7,:,:)=upvpsum
    sumStat( 8,:,:)=wxsum;    sumStat( 9,:,:)=wysum;    sumStat(10,:,:)=wzsum
    sumStat(11,:,:)=wxwxsum;  sumStat(12,:,:)=wywysum;  sumStat(13,:,:)=wzwzsum    
    call MPI_REDUCE(npsum,  npsslot,    npType*nslab, int_type, MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(sumStat,sumStatR,13*npType*nslab, real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    
    if(nrank==0) then
      inpstime = one/real(npstime,RK)
      do itype=1,npType
        write(filename,"(A,I10.10,A,I4.4)") trim(DEM_Opt%ResultsDir)//'pstats',itime,'_set',itype
        open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
        if(ierror /= 0) then
          call DEMLogInfo%CheckForError(ErrT_Pass,"ClcCAStatistics","Cannot open file: "//trim(filename))
        else              
          write(nUnit,'(a,I7,a,I7,a,I7)')'    The time step range for this particle statistics is ', &
                                       itime-(npstime-1)*ivstats, ':', ivstats, ':', itime
          write(nUnit,*)
          write(nUnit,'(A)')'  yp, np, up, vp, wp, upup, vpvp, wpwp, upvp, wx, wy, wz, wxwx, wywy, wzwz' 
          do jc=1,nslab
            if(npsslot(itype,jc)==0) then
              write(nUnit,'(15E24.15)') (ypForPs(jc)+ypForPs(jc+1))*half, & ! 1
                zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero
            else
              irnpsum=one/real(npsslot(itype,jc),RK)
              write(nUnit,'(15ES24.15)') (ypForPs(jc)+ypForPs(jc+1))*half, & ! 1
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
                                          sumStatR(13,itype,jc)*irnpsum    ! 15
            endif
          enddo
          close(nUnit,IOSTAT=ierror)
        endif
      enddo
    endif
    call ResetStatVar()
  end subroutine ClcCAStatistics
end module ca_Statistics
