#include "definitions_inc.f90"
module cd_Statistics
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d,only:nrank
  use f2_MeshAndMetries,only: yp
  use f2_Parameters,only: ivstats,saveStat,nyc,itime,xlx,zlz,ilast
  use sp_Property
  use sp_Variables
  use sp_Decomp_2d
  use sp_Parameters
  implicit none
  private

  integer::npType,nslab,npstime
  integer,dimension(:,:),allocatable::npsum
  real(RK),dimension(:),allocatable:: ypForPs  ! y point for particle statistics
  real(RK),dimension(:,:),allocatable::upsum,vpsum,wpsum,upupsum,vpvpsum,wpwpsum,upvpsum

  public::InitCDStatistics,ClcCDStatistics
#define ClcTransportRate
contains
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  !   InitCDStatistics
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitCDStatistics(chFile); implicit none
    character(len=*),intent(in)::chFile
        
    ! locals
    character(len=:),allocatable::filename
    NAMELIST/ParticleStatisticOption/nslab
    integer::j,k,iErr01,iErr02,iErr03,ierr,iUnit

    npType= DEM_Opt%numPrtcl_Type
    open(newunit=iUnit, file=chFile,status='old',form='formatted',IOSTAT=ierr)
    if(ierr/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitDEMStatistics", "Cannot open file: "//strip(chFile))
    read(iUnit, nml=ParticleStatisticOption)
    if(nrank==0)write(DEMLogInfo%iUnit, nml=ParticleStatisticOption)
    close(iUnit,IOSTAT=ierr)
    if(nrank==0 .and. mod(nyc,nslab)/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitCDStatistics","mod(nyc,nslab)/=0")
    
    allocate(ypForPs(nslab+1),Stat=iErr01)
    allocate(npsum(npType,nslab), upsum(npType,nslab),   vpsum(npType,nslab),   wpsum(npType,nslab),   Stat=iErr02)
    allocate(upupsum(npType,nslab), vpvpsum(npType,nslab), wpwpsum(npType,nslab), upvpsum(npType,nslab), Stat=iErr03)
    ierr=abs(iErr01)+abs(iErr02)+abs(iErr03)
    if(ierr/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitCDStatistics ","Allocation failed")
    k=nyc/nslab
    do j=0,nslab
      ypForPs(j+1)=yp(j*k+1)
    enddo
    call ResetStatVar()

#ifdef ClcTransportRate
    if(nrank==0) then
      filename = strip(DEM_opt%ResultsDir)//"TransportRate"//int2str(ilast,10)
      open(newunit=iUnit,file=filename,status='replace',form='formatted',IOSTAT=ierr)
      if(ierr/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"InitCDStatistics","Cannot open file: "//filename)
      close(iUnit,IOSTAT=ierr)
    endif
#endif
  end subroutine InitCDStatistics

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ResetStatVar
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ResetStatVar(); implicit none

    npstime=0;         npsum=0
    upsum=0.0_RK;      vpsum=0.0_RK;      wpsum=0.0_RK
    upupsum=0.0_RK;    vpvpsum=0.0_RK;    wpwpsum=0.0_RK;    upvpsum=0.0_RK
  end subroutine ResetStatVar

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ClcCDStatistics
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90  
  subroutine ClcCDStatistics(); implicit none

    ! locals
    type(real3)::pos,VPrtcl
    character(len=:),allocatable:: filename
    integer,dimension(npType,nslab)::npsslot
    integer::pid,nlocal,itype,js,je,jc,ierr,iUnit
    real(RK)::irnpsum,inpstime,SumVel,SumVelR
    real(RK),dimension(7,npType,nslab)::sumStat,sumStatR

    sumStat=0;  filename= ' '
    nlocal= GPrtcl_list%nlocal
    do pid=1,nlocal
      pos  = GPrtcl_posR(pid)
      itype= GPrtcl_pType(pid)
      VPrtcl=GPrtcl_linVel(1,pid)

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
    enddo
    do jc=1,nslab
      do itype=1,npType
        upsum(itype,jc)  = upsum(itype,jc)  + sumStat(1,itype,jc)
        vpsum(itype,jc)  = vpsum(itype,jc)  + sumStat(2,itype,jc)
        wpsum(itype,jc)  = wpsum(itype,jc)  + sumStat(3,itype,jc)
        upupsum(itype,jc)= upupsum(itype,jc)+ sumStat(4,itype,jc)
        vpvpsum(itype,jc)= vpvpsum(itype,jc)+ sumStat(5,itype,jc)
        wpwpsum(itype,jc)= wpwpsum(itype,jc)+ sumStat(6,itype,jc)
        upvpsum(itype,jc)= upvpsum(itype,jc)+ sumStat(7,itype,jc)
      enddo
    enddo

#ifdef ClcTransportRate
    SumVel=0.0_RK
    do pid=1,nlocal
      itype = GPrtcl_pType(pid)
      SumVel= SumVel+ GPrtcl_linVel(1,pid)%x* DEMProperty%Prtcl_PureProp(itype)%Volume
    enddo
    call MPI_REDUCE(SumVel,SumVelR,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(nrank==0) then
      filename = strip(DEM_opt%ResultsDir)//"TransportRate"//int2str(ilast,10)
      open(newunit=iUnit,file=filename,status='old',form='formatted',position='append',IOSTAT=ierr)
      if(ierr /= 0) then
        call DEMLogInfo%CheckForError(ErrT_Pass,"ClcLCDStatistics: ","Cannot open file: "//filename)
      else
        write(iUnit,'(I7,ES24.15)')itime,SumVelR/(xlx*zlz)
        close(iUnit,IOSTAT=ierr)
      endif
    endif
#endif

    npstime = npstime + 1
    if(mod(itime,SaveStat)/=0) return
    sumStat(1,:,:)=upsum;    sumStat(2,:,:)=vpsum;    sumStat(3,:,:)=wpsum
    sumStat(4,:,:)=upupsum;  sumStat(5,:,:)=vpvpsum;  sumStat(6,:,:)=wpwpsum;  sumStat(7,:,:)=upvpsum
    call MPI_REDUCE(npsum,  npsslot,   npType*nslab, int_type, MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(sumStat,sumStatR,7*npType*nslab, real_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    if(nrank==0) then
      inpstime = 1.0_RK/real(npstime,RK)
      do itype=1,npType
        filename = strip(DEM_Opt%ResultsDir) //'pstats' // int2str(itime,10) // '_set' //int2str(itype,4)
        open(newunit=iUnit,file=filename,status='replace',form='formatted',IOSTAT=ierr)
        if(ierr /= 0) then
          call DEMLogInfo%CheckForError(ErrT_Pass,"ClcCDStatistics","Cannot open file: "//filename)
        else              
          write(iUnit,'(a,I7,a,I7,a,I7)')'    The time step range for this particle statistics is ', &
                                       itime-(npstime-1)*ivstats, ':', ivstats, ':', itime
          write(iUnit,*)
          write(iUnit,'(A)')'  yp, np, up, vp, wp, upup, vpvp, wpwp, upvp' 
          do jc=1,nslab
            if(npsslot(itype,jc)==0) then
              write(iUnit,'(15E24.15)') (ypForPs(jc)+ypForPs(jc+1))*0.5_RK, & ! 1
                                        0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK,0.0_RK
            else
              irnpsum=1.0_RK/real(npsslot(itype,jc),RK)
              write(iUnit,'(15ES24.15)') (ypForPs(jc)+ypForPs(jc+1))*0.5_RK, & ! 1
                                   real(npsslot(itype,jc),RK)*inpstime, & ! 2
                                          sumStatR(1,itype,jc)*irnpsum, & ! 3
                                          sumStatR(2,itype,jc)*irnpsum, & ! 4
                                          sumStatR(3,itype,jc)*irnpsum, & ! 5
                                          sumStatR(4,itype,jc)*irnpsum, & ! 6
                                          sumStatR(5,itype,jc)*irnpsum, & ! 7
                                          sumStatR(6,itype,jc)*irnpsum, & ! 8
                                          sumStatR(7,itype,jc)*irnpsum    ! 9
            endif
          enddo
          close(iUnit,IOSTAT=ierr)
        endif
      enddo
    endif
    call ResetStatVar()
  end subroutine ClcCDStatistics
end module cd_Statistics
