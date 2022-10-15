program main_channelLPT
  use MPI
  use m_Timer
  use m_Tools
  use m_LogInfo
  use LPT_System
  use m_Decomp2d
  use LPT_Fpforce
  use m_Variables
  use m_IOAndVisu  
  use m_Parameters
  use LPT_decomp_2d
  use LPT_IOAndVisu
  use LPT_Variables
  use LPT_Parameters
  use m_ChannelSystem
  use m_MeshAndMetries
  use m_TypeDef,only:num2str
  use m_Poisson,only:Destory_Poisson_FFT_Plan
  implicit none
  integer::intT,pid,ierror
  type(timer)::CoupleTimer
  character(len=128)::chPrm
  character(len=10)::RowColStr
#ifdef CFDFourthOrder
  integer::BcOption(6)
#endif

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
  
#ifdef CFDFourthOrder
  if(FlowType==FT_CH) then
    BcOption=(/0,0,-1,-1,0,0/)
  elseif(FlowType==FT_HC) then
    BcOption=(/0,0,-1,-2,0,0/)
  endif
#endif
  call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,y_pencil,BcOption)
  call ChannelInitialize(chPrm)    ! Topest level initialing for Channel body

  ! ================== initialize LPT options ==================
  call get_command_argument(2,chPrm)
  call LPT_opt%ReadLPTOption( chPrm)
  call LPT_decomp%Init_DECOMP()
  call LPT%Initialize(chPrm)       ! Topest level initialing for LPT body

  ! ================== initialize CFD-LPT coupling part ==================
  call InitFpForce(chPrm)

  ! =============== dump initial visulizing files ===============
  asso_Q: associate(Q_vor =>RealArr1)
  call  dump_visu(ifirst-1,ux,uy,uz,pressure,Q_vor)  ! channel3d
  end associate asso_Q
  if(LPT_Opt%RestartFlag)call LPT_IO%dump_visu(ifirst-1)
  print*, nrank,GPrtcl_list%nlocal

  call CoupleTimer%reset()
  do itime=ifirst,ilast

    ! CFD-LPT coupling part
    call CoupleTimer%start()
    call PrepareInterpolation()
    call clc_VelInterpolation(ux,uy,uz)
    if(itime==ifirst .and. (.not. LPT_Opt%RestartFlag)) then
      do pid=1,GPrtcl_list%nlocal
        GPrtcl_linVel(:,pid)=GPrtcl_VFluid(pid)
      enddo
      call LPT_IO%dump_visu(ifirst-1)                    ! LPT
    endif
    call clc_FpForce() !clc_FpForce(ux,uy,uz,pressure)
#ifdef CFDLPT_TwoWay
    call distribute_FpForce()
#endif

    call FinalFpForce()
    call CoupleTimer%finish()

    ! Largrangian Particle Trackiing part
    call LPT%iterate(itime)
    
    ! CFD Iterate
    call ChannelIterate()
    
    if(nrank==0 .and. mod(itime, Cmd_LFile_Freq)==0) then
      call MainLog%OutInfo("Coupling time  [tot, last, ave] [sec]: "//trim(num2str(CoupleTimer%tot_time))//", "// &
          trim(num2str(CoupleTimer%last_time ))//", "//trim(num2str(CoupleTimer%average())),2)
    endif  
  enddo
  call LPT_IO%Final_visu()

  if(nrank==0)call MainLog%OutInfo("Good job! ChannelLPM finished successfully at "//time2str(),1)
  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)
end program main_channelLPT
