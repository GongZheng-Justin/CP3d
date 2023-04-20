program main_channel3d_2nd
  use MPI
  use m_Decomp2d
  use m_Variables
  use m_IOAndVisu
  use m_Parameters
  use m_ChannelSystem
  use m_Timer,only:time2str  
  use m_Poisson,only:Destory_Poisson_FFT_Plan
  implicit none
  integer::intT,ierror
  character(len=128)::chPrm
  character(len=10)::RowColStr
  
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! read Channel options
  intT=command_argument_count()
  if((intT/=1 .and. intT/=3) .and. nrank==0) then
    write(*,*)'command argument wrong!'; stop
  endif
  call get_command_argument(1,chPrm)
  call ReadAndInitParameters(chPrm)
  if(intT==3) then
    call get_command_argument(2,RowColStr)
    read(RowColStr,*) p_row
    call get_command_argument(3,RowColStr)
    read(RowColStr,*) p_col
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  ! Initialize Decomp-2d
  call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,y_pencil,BcOption)
 
  call ChannelInitialize(chPrm)    ! Topest level initialing for Channel body
  do itime=ifirst, ilast
    call ChannelIterate()
  enddo
  if(nrank==0)call MainLog%OutInfo("Good job! Channel3d_2nd finished successfully at "//time2str(),1)

  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize()
  call MPI_FINALIZE(ierror)
end program main_channel3d_2nd
