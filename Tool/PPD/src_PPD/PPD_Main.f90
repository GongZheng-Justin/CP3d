program main
  use MPI
  use PPD_Tools
  use PPD_TypeDef
  use PPD_Parameters
  use PPD_Statistics
  use PPD_ReadAndSortDump
  implicit none
  integer::ierror
  real(RK)::SetTime,EndTime
  character(len=256)::PrmStr

  if(nrank==0)SetTime=MPI_WTIME()
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)

  ! initialize 
  if(command_argument_count()/=1 .and. nrank==0) write(*,*)'command argument wrong!'
  call get_command_argument(1,PrmStr)
  call ReadParameters(PrmStr)
  call FileSystem%Initialize()

  PrmStr= 'mkdir -p '//Stat_Dir//' 2> /dev/null'
  if(nrank==0) call system(trim(adjustl(PrmStr)))

  ! Step01: Read particle dump file
  if(nrank==0 .and. IsReadPrtclDumpOut) call ReadPrtclDump()
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  ! Step02: Sort dump file in time sequence.
  if(IsSortPrtclDumpOut) call SortPrtclDump()
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  ! Step03: Shrink dump file if necessary
  if(IsShrinkPrtclDump) call ShrinkPrtclDump()
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  ! Step04: Calculate the necessary Hop statistic data
  if(IsClcHopStat) call ClcHopStat()
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  ! Step05: Calculate the necessary Cross-diffusion statistic data
  if(IsClcCrossDiffusion) call ClcCrossDiffusion()
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  if(nrank==0) then
    EndTime=MPI_WTIME()
    print*,'Consuming time is [s]:',EndTime-SetTime
  endif
  call MPI_FINALIZE(ierror)
end program main
