!******************************************************************
! Main program
!******************************************************************
program main
  implicit none
  integer,parameter::RK=8
  integer,parameter::nxc=8640
  integer,parameter::nyc=256
  integer,parameter::nzc=1680
  character(len=128)::FileIn,FileOut
  real(RK),allocatable,dimension(:,:,:)::ArrTmp
  integer::nUnitIn,nUnitOut,it,ierror
  
  allocate(ArrTmp(nxc,nyc,nzc))
  FileIn='RestartForOCT1000L_090000150000'
  FileOut='Restart'
  open(newunit=nUnitIn,file=FileIn,form='unformatted',access='stream',status='old',action='read',position='rewind',iostat=ierror)
  open(newunit=nUnitOut,file=FileOut,form='unformatted',access='stream',status='replace',action='write',iostat=ierror)
  do it=1,4
    print*,it
    read(unit=nUnitIn)ArrTmp
    write(unit=nUnitOut)ArrTmp(1:2160,1:256,1:1120)
  enddo
  close(unit=nUnitIn,iostat=ierror)
  close(unit=nUnitOut,iostat=ierror)    
end program main
