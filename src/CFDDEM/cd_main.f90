program main_channelDEM
  use MPI
  use cd_System
  use m_LogInfo
  use m_Decomp2d
  use cd_FpForce
  use m_Variables
  use m_IOAndVisu
  use m_Parameters
  use cd_Statistics
  use m_ChannelSystem
  use Prtcl_decomp_2d
  use Prtcl_DEMSystem
  use Prtcl_IOAndVisu
  use Prtcl_Variables
  use Prtcl_Parameters
  use m_Timer,only:time2str
  use Prtcl_TypeDef,only:zero_r3
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
  call DEM_Opt%ReadDEMOption( chPrm)
  call DEM_decomp%Init_DECOMP()
  call DEM%Initialize(chPrm) ! Topest level initialing for DEM body

  ! ================== initialize CFD-DEM coupling part ==================
  call InitDistribute()
  call InitCDStatistics(chPrm) 

  print*, nrank,GPrtcl_list%nlocal,GPrtcl_list%mlocalFix
  call ChannelDEM_Iterate()
  call DEM_IO%Final_visu()

  if(nrank==0)call MainLog%OutInfo("Good job! ChannelDEM finished successfully at "//time2str(),1)
  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
end program main_channelDEM
