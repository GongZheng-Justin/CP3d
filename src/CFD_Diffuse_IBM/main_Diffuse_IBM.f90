#include "definitions_inc.f90"
program main_Diffuse_IBM
  use MPI
  use mc_Decomp2d
  use mc_TypeDef,only:strip
  use mc_Timer,only:time2str
  use f2_Variables
  use f2_IOAndVisu
  use f2_CFDSystem
  use f2_Parameters
  use f2_Poisson,only:Destory_Poisson_FFT_Plan
  use dIBM_IBM
  use dIBM_system
  implicit none
  integer::ierr

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! ================== initialize CFD options ==================
  block 
    character(len=10)::RowColStr
    character(len=512)::caPrmTmp
    character(len=:),allocatable::caPrm
  
    ierr=command_argument_count()
    if((ierr/=1 .and. ierr/=3) .and. nrank==0) then
      write(*,*)'command argument wrong!'; stop
    endif
    call get_command_argument(1,caPrmTmp)
    caPrm = strip(caPrmTmp)
    call ReadAndInitParameters(caPrm)
    if(ierr==3) then
      call get_command_argument(2,RowColStr)
      read(RowColStr,*) p_row
      call get_command_argument(3,RowColStr)
      read(RowColStr,*) p_col
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,y_pencil,BcOption)
    call CFDInitialize(caPrm)    ! Topest level initialing for CFD body

    call get_command_argument(1,caPrmTmp)
    caPrm = strip(caPrmTmp)
    call Init_IBM(caPrm)
#ifdef IBMDistributeLinear
    if(nrank==0) call MainLog%OutInfo("Choose to distribute IBM force using linear function",2)
#else
    if(nrank==0) call MainLog%OutInfo("Choose to distribute IBM force using 3-point Dirac function",2)
#endif
  end block

  call CFDACM_Iterate()

  if(nrank==0) call MainLog%OutInfo("Good job! CFD_ACM finished successfully at "//time2str(),1)
  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
end program main_Diffuse_IBM
