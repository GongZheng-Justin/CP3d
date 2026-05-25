#include "definitions_inc.f90"
program main_CFDLPT
  use MPI
  use mc_Timer
  use mc_Decomp2d
  use mc_TypeDef,only:strip
  use mc_TypeDef,only:num2str
#ifdef CFDSecondOrder
  use f2_Tools
  use f2_Variables
  use f2_IOAndVisu
  use f2_CFDSystem
  use f2_Parameters
  use f2_MeshAndMetries
  use f2_Poisson,only:Destory_Poisson_FFT_Plan
#else
  use f4_Tools
  use f4_Variables
  use f4_IOAndVisu
  use f4_CFDSystem
  use f4_Parameters
  use f4_MeshAndMetries
  use f4_Poisson,only:Destory_Poisson_FFT_Plan
#endif
  use lp_System
  use lp_Fpforce
  use lp_decomp_2d
  use lp_IOAndVisu
  use lp_Variables
  use lp_Parameters
  implicit none
  integer::ierr
  type(timer)::CoupleTimer

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! ================== initialize CFD options ==================  
  block
    character(len=10)::RowColStr
    character(len=512)::LPTPrmTmp
    character(len=:),allocatable::LPTPrm
#ifdef CFDFourthOrder
    integer::BcOption(6)
#endif
    
    ierr=command_argument_count()
    if((ierr/=2 .and. ierr/=4) .and. nrank==0) then
      write(*,*)'command argument wrong!'; stop
    endif
    call get_command_argument(1,LPTPrmTmp)
    LPTPrm = strip(LPTPrmTmp)
    call ReadAndInitParameters(LPTPrm)
    if(ierr==4) then
      call get_command_argument(3,RowColStr)
      read(RowColStr,*) p_row
      call get_command_argument(4,RowColStr)
      read(RowColStr,*) p_col
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
#ifdef CFDFourthOrder
    if(FlowType==FT_CH) then
      BcOption=(/0,0,-1,-1,0,0/)
    elseif(FlowType==FT_HC) then
      BcOption=(/0,0,-1,-2,0,0/)
    endif
#endif
    call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,y_pencil,BcOption)
    call CFDInitialize(LPTPrm)    ! Topest level initialing for CFD body

    ! ================== initialize LPT options ==================
    call get_command_argument(2,LPTPrmTmp)
    LPTPrm = strip(LPTPrmTmp)
    call LPT_Opt%ReadLPTOption(LPTPrm)
    call LPT_decomp%Init_DECOMP()
    call LPT%Initialize(LPTPrm)       ! Topest level initialing for LPT body

    ! ================== initialize CFD-LPT coupling part ==================
    call InitFpForce(LPTPrm)
  end block

  ! =============== dump initial visulizing files ===============
  asso_Q: associate(Q_vor =>RealArr1)
  call  dump_visu(ifirst-1,ux,uy,uz,pressure,Q_vor)  ! CFD3d
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
      do ierr=1,GPrtcl_list%nlocal
        GPrtcl_linVel(:,ierr)=GPrtcl_VFluid(ierr)
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
    call CFDIterate()
    
    if(nrank==0 .and. mod(itime, Cmd_LFile_Freq)==0) then
      call MainLog%OutInfo("Coupling time  [tot, last, ave] [sec]: "//strip(num2str(CoupleTimer%tot_time))//", "// &
          strip(num2str(CoupleTimer%last_time ))//", "//strip(num2str(CoupleTimer%average())),2)
    endif  
  enddo
  call LPT_IO%Final_visu()

  if(nrank==0)call MainLog%OutInfo("Good job! CFD_LPT finished successfully at "//time2str(),1)
  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
end program main_CFDLPT
