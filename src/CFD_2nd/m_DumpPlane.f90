module m_DumpPlane
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use m_Variables,only: mb1
  implicit none
  private
#define RKP_Dump 4

  logical,parameter::DumpPlaneFlag=.false.
  
  integer,parameter::nyPlane=3
  integer,dimension(nyPlane),parameter::iyPlane=[1, 2, 3]
  
  integer::nTimePlane,iTimePlane
  real(RKP_Dump),dimension(:,:,:,:),allocatable::uxcPlane_y, uycPlane_y, uypPlane_y, uzcPlane_y, prcPlane_y
  
  public::Initialize_DumpPlane, dump_plane
contains

  !******************************************************************
  ! Initialize_DumpPlane
  !******************************************************************
  subroutine Initialize_DumpPlane()
    implicit none
    
    ! locals
    integer::ierrTmp,ierror=0
    
    if(.not. DumpPlaneFlag) return
    if(mod(BackupFreq, ivstats)/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"Initialize_DumpPlane","nTimePlane WRONG")
    endif
    nTimePlane=BackupFreq/ivstats; iTimePlane=0
    allocate(uxcPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane,nyPlane),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(uycPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane,nyPlane),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)    
    allocate(uypPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane,nyPlane),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(uzcPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane,nyPlane),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)    
    allocate(prcPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane,nyPlane),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)          
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"Initialize_DumpPlane","Allocation failed")
  end subroutine Initialize_DumpPlane
  
  !******************************************************************
  ! dump_plane
  !******************************************************************
  subroutine dump_plane(ntime,ux,uy,uz,pr)
    implicit none
    integer,intent(in)::ntime
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pr
    
    ! locals
    character(len=128)::filename
    integer,dimension(3)::sizes,subsizes,starts
    integer::ic,jc,kc,jplane1,jplane2,data_type,ierror,newtype,nUnit
        
    if(.not. DumpPlaneFlag) return
    iTimePlane=iTimePlane+1
    
    do kc=y1start(3),y1end(3)
      do jc=1,nyPlane
        jplane1=iyPlane(jc)
        jplane2=jplane1+1
        do ic=y1start(1),y1end(1)
          uxcPlane_y(ic,kc,iTimePlane,jc)=real(ux(ic,jplane1,kc), RKP_Dump)
          uycPlane_y(ic,kc,iTimePlane,jc)=real(uy(ic,jplane1,kc), RKP_Dump)
          uzcPlane_y(ic,kc,iTimePlane,jc)=real(uz(ic,jplane1,kc), RKP_Dump)
          prcPlane_y(ic,kc,iTimePlane,jc)=real(pr(ic,jplane1,kc), RKP_Dump)          
          uypPlane_y(ic,kc,iTimePlane,jc)=real(uy(ic,jplane2,kc), RKP_Dump)
        enddo
      enddo
    enddo
    if(iTimePlane/=nTimePlane) return

    ! Write plane info ========================
    if((RKP_Dump-4)==0) then
      data_type= MPI_REAL
    else    
      data_type= MPI_DOUBLE_PRECISION    
    endif
    sizes(1)= nxc
    sizes(2)= nzc
    sizes(3)= nTimePlane
    subsizes(1)= y1size(1)
    subsizes(2)= y1size(3)
    subsizes(3)= nTimePlane
    starts(1)= y1start(1)-1
    starts(2)= y1start(3)-1
    starts(3)= 0
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
        
    ! uxcPlane_y
    do jc=1,nyPlane
      write(filename,"(A,A,I2.2,A,I10.10)") trim(ResultsDir),'uxcPlane_y',jc,'_',ntime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FILE_SET_VIEW(nUnit,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(nUnit,uxcPlane_y(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    enddo
    
    ! uycPlane_y
    do jc=1,nyPlane
      write(filename,"(A,A,I2.2,A,I10.10)") trim(ResultsDir),'uycPlane_y',jc,'_',ntime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FILE_SET_VIEW(nUnit,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(nUnit,uycPlane_y(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    enddo    

    ! uypPlane_y
    do jc=1,nyPlane
      write(filename,"(A,A,I2.2,A,I10.10)") trim(ResultsDir),'uypPlane_y',jc,'_',ntime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FILE_SET_VIEW(nUnit,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(nUnit,uypPlane_y(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    enddo  

    ! uzcPlane_y
    do jc=1,nyPlane
      write(filename,"(A,A,I2.2,A,I10.10)") trim(ResultsDir),'uzcPlane_y',jc,'_',ntime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FILE_SET_VIEW(nUnit,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(nUnit,uzcPlane_y(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    enddo  

    ! prcPlane_y
    do jc=1,nyPlane
      write(filename,"(A,A,I2.2,A,I10.10)") trim(ResultsDir),'prcPlane_y',jc,'_',ntime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_FILE_SET_VIEW(nUnit,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(nUnit,prcPlane_y(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
      call MPI_FILE_CLOSE(nUnit,ierror)
    enddo  
    
    call MPI_TYPE_FREE(newtype,ierror)    
    iTimePlane=0
  end subroutine dump_plane
end module m_DumpPlane

#undef RKP_Dump
