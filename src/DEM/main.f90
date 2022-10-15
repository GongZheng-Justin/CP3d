program main_DEM
  use MPI
  use Prtcl_LogInfo
  use Prtcl_decomp_2d
  use Prtcl_DEMSystem
  use Prtcl_IOAndVisu
  use Prtcl_Variables
#ifdef TestDEMRestart
  use Prtcl_CL_and_CF
#endif
  use Prtcl_Parameters
  use Prtcl_Timer,only:time2str
  implicit none
  character(len=10)::RowColStr
  character(len=128)::chDEMPrm
  integer::ierror,intT,RowCol(2)
  
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! read DEM options
  intT=command_argument_count()
  if((intT/=1 .and. intT/=3) .and. nrank==0) then
    write(*,*)'command argument wrong!'; stop
  endif
  call get_command_argument(1,chDEMPrm)
  call DEM_opt%ReadDEMOption( chDEMPrm)
  if(intT==3) then
    call get_command_argument(2,RowColStr)
    read(RowColStr,*) RowCol(1)
    call get_command_argument(3,RowColStr)
    read(RowColStr,*) RowCol(2)
    call DEM_decomp%Init_DECOMP(chDEMPrm,RowCol)
  else
    call DEM_decomp%Init_DECOMP(chDEMPrm)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  
  call DEM%Initialize(chDEMPrm) ! Topest level initialing for DEM body
#ifdef TestDEMRestart
  call GPPW_CntctList%printCL(DEM_opt%ifirst-1)
#endif

  print*, nrank,GPrtcl_list%nlocal,GPrtcl_list%mlocalFix
  do intT= DEM_opt%ifirst, DEM_opt%ilast
    call DEM%iterate(intT)
  enddo
#ifdef TestDEMRestart
  call GPPW_CntctList%printCL(DEM_opt%ilast)
#endif
  call DEM_IO%Final_visu()

  if(nrank==0)call DEMLogInfo%OutInfo("Good job! DEM finished successfully at "//time2str(),1)
  call MPI_FINALIZE(ierror)
end program main_DEM
