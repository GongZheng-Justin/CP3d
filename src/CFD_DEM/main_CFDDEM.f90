#include "definitions_inc.f90"
program main_CFDDEM
  use MPI
  use mc_Decomp2d
  use mc_TypeDef,only:strip
  use mc_Timer,only:time2str
  use f2_Variables
  use f2_IOAndVisu
  use f2_CFDSystem
  use f2_Parameters
  use f2_Poisson,only:Destory_Poisson_FFT_Plan
  use sp_System
  use sp_decomp_2d
  use sp_IOAndVisu
  use sp_Variables
  use sp_Parameters
  use cd_System
  use cd_FpForce
  use cd_Statistics
  implicit none
  integer::ierr

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! ================== initialize CFD options ==================
  block
    character(len=10)::RowColStr
    character(len=512)::cdPrmTmp
    character(len=:),allocatable::cdPrm
    ierr=command_argument_count()
    if((ierr/=2 .and. ierr/=4) .and. nrank==0) then
      write(*,*)'command argument wrong!'; stop
    endif
    call get_command_argument(1,cdPrmTmp)
    cdPrm = strip(cdPrmTmp)
    call ReadAndInitParameters(cdPrm)
    if(ierr==4) then
      call get_command_argument(3,RowColStr)
      read(RowColStr,*) p_row
      call get_command_argument(4,RowColStr)
      read(RowColStr,*) p_col
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  
    call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,y_pencil,BcOption)
    call CFDInitialize(cdPrm)    ! Topest level initialing for CFD body

    ! ================== initialize DEM options ==================
    call get_command_argument(2,cdPrmTmp)
    cdPrm = strip(cdPrmTmp)
    call DEM_Opt%ReadDEMOption(cdPrm)
    call DEM_decomp%Init_DECOMP()
    call DEM%Initialize(cdPrm) ! Topest level initialing for DEM body

    ! ================== initialize CFD-DEM coupling part ==================
    call InitDistribute()
    call InitCDStatistics(cdPrm) 
  end block
  
  print*, nrank,GPrtcl_list%nlocal,GPrtcl_list%mlocalFix
  call CFDDEM_Iterate()
  call Prtcl_Final_Visu()

  if(nrank==0) call MainLog%OutInfo("Good job! CFD_DEM finished successfully at "//time2str(),1)
  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
end program main_CFDDEM
