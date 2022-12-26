program main_channelACM
  use MPI
  use ca_IBM
  use ca_system
  use m_Decomp2d
  use m_Variables
  use m_IOAndVisu
  use m_Parameters
  use Prtcl_System
  use ca_Statistics
  use m_ChannelSystem
  use Prtcl_decomp_2d
  use Prtcl_IOAndVisu
  use Prtcl_Variables
  use ca_IBM_implicit
  use Prtcl_Parameters
  use m_Timer,only:time2str
  use m_Poisson,only:Destory_Poisson_FFT_Plan
  implicit none
  integer::intT,ierror
  character(len=128)::chPrm
  character(len=10)::RowColStr

  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! ================== initialize channel options ==================
  intT=command_argument_count()
  if((intT/=2 .and. intT/=4) .and. nrank==0) then
    write(*,*)'command argument wrong!'; stop
  endif
  call get_command_argument(1,chPrm)
  call ReadAndInitParameters(chPrm)
  if(intT==4) then
    call get_command_argument(3,RowColStr)
    read(RowColStr,*) p_row
    call get_command_argument(4,RowColStr)
    read(RowColStr,*) p_col
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    
  call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,y_pencil,BcOption)
  call ChannelInitialize(chPrm)    ! Topest level initialing for Channel body

  ! ================== initialize DEM options ==================
  call get_command_argument(2,chPrm)
  call DEM_opt%ReadDEMOption( chPrm)
  call DEM_decomp%Init_DECOMP()
  call DEM%Initialize(chPrm) ! Topest level initialing for DEM body

  print*, nrank,GPrtcl_list%nlocal,GPrtcl_list%mlocalFix
  if(IBM_Scheme<2) then
    call Init_IBM(chPrm)
  else
    call Init_IBM_imp(chPrm)
  endif
  call InitCAStatistics(chPrm)

  call ChannelACM_Iterate()
  call DEM_IO%Final_visu()

  if(nrank==0)call MainLog%OutInfo("Good job! ChannelACM finished successfully at "//time2str(),1)
  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
end program main_channelACM
