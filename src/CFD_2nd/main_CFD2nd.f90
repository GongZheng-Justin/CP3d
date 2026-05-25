#include "definitions_inc.f90"
program main_CFD_2nd
  use MPI
  use mc_Decomp2d
  use mc_TypeDef,only:strip
  use mc_Timer,only:time2str
  use f2_Variables
  use f2_IOAndVisu
  use f2_Parameters
  use f2_CFDSystem 
  use f2_Poisson,only:Destory_Poisson_FFT_Plan
  implicit none
  integer::ierr
  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! read CFD options
  ierr=command_argument_count()
  if((ierr/=1 .and. ierr/=3) .and. nrank==0) then
    write(*,*)'command argument wrong!'; stop
  endif
  block
    character(len=10)::RowColStr
    character(len=512)::chPrmTmp
    character(len=:),allocatable::chPrm
    
    call get_command_argument(1,chPrmTmp)
    chPrm = strip(chPrmTmp)
    call ReadAndInitParameters(chPrm)
    if(ierr==3) then
      call get_command_argument(2,RowColStr)
      read(RowColStr,*) p_row
      call get_command_argument(3,RowColStr)
      read(RowColStr,*) p_col
    endif
    
    ! Initialize Decomp-2d
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,y_pencil,BcOption)
    call CFDInitialize(chPrm)    ! Topest level initialing for CFD body
  end block

  do itime=ifirst, ilast
    call CFDIterate()
  enddo
  if(nrank==0) call MainLog%OutInfo("Good job! CFD_2nd finished successfully at "//time2str(),1)

  call Destory_Poisson_FFT_Plan()
  call decomp_2d_finalize()
  call MPI_FINALIZE(ierr)
end program main_CFD_2nd
