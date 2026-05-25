#include "definitions_inc.f90"
program main_DEM
  use MPI
  use mc_TypeDef,only:strip
  use mc_Timer,only:time2str
  use sp_System
  use sp_decomp_2d
  use sp_IOAndVisu
  use sp_Variables
#ifdef TestDEMRestart
  use sp_CL_and_CF
#endif
  use sp_Parameters
  implicit none
  integer::ierr
  
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! read DEM options
  ierr=command_argument_count()
  if((ierr/=1 .and. ierr/=3) .and. nrank==0) then
    write(*,*)'command argument wrong!'; stop
  endif
  
  block
    integer::RowCol(2)
    character(len=10)::RowColStr
    character(len=512)::DEMPrmTmp
    character(len=:),allocatable::DEMPrm
    
    call get_command_argument(1,DEMPrmTmp)
    DEMPrm = strip(DEMPrmTmp)
    call DEM_opt%ReadDEMOption(DEMPrm)
    if(ierr==3) then
      call get_command_argument(2,RowColStr)
      read(RowColStr,*) RowCol(1)
      call get_command_argument(3,RowColStr)
      read(RowColStr,*) RowCol(2)
      call DEM_decomp%Init_DECOMP(DEMPrm,RowCol)
    else
      call DEM_decomp%Init_DECOMP(DEMPrm)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call DEM%Initialize(DEMPrm) ! Topest level initialing for DEM body
  end block
#ifdef TestDEMRestart
  call GPPW_CntctList%printCL(DEM_opt%ifirst-1)
#endif

  print*, nrank,GPrtcl_list%nlocal,GPrtcl_list%mlocalFix
  do ierr= DEM_opt%ifirst, DEM_opt%ilast
    call DEM%iterate(ierr)
  enddo
#ifdef TestDEMRestart
  call GPPW_CntctList%printCL(DEM_opt%ilast)
#endif
  call Prtcl_Final_Visu()

  if(nrank==0)call DEMLogInfo%OutInfo("Good job! DEM finished successfully at "//time2str(),1)
  call MPI_FINALIZE(ierr)
end program main_DEM
