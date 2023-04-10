program main
  implicit none
  integer,parameter::RKR=8
  integer,parameter::RKW=4
  character(128),parameter::StrRead ='PartVisuForFixedTotal_P0000000000_'
  character(128),parameter::StrWrite='PartVisuForFixedTotal_P0000000000'
    
  ! locals
  integer::ierror,nUnit
  integer(8)::disp,RealNum
  real(RKR),dimension(:),allocatable::VecRead
  real(RKW),dimension(:),allocatable::VecWrite

  open(newunit=nUnit,file=trim(StrRead),form='unformatted',action='read',status='old', access='stream',position='append',iostat=ierror)
  inquire(unit=nUnit,Pos=disp); disp=disp-1_8
  RealNum=disp/RKR
  rewind(unit=nUnit,IOSTAT=ierror)
  print*,' Total real num=',RealNum
  allocate(VecRead(RealNum))
  allocate(VecWrite(RealNum))
  read(unit=nUnit)VecRead
  close(unit=nUnit,iostat=ierror)

  VecWrite=reaL(VecRead,RKW)  
  open(newunit=nUnit,file=trim(StrWrite),form='unformatted',action='write',status='replace',access='stream',iostat=ierror)
  write(unit=nUnit)VecWrite
  close(unit=nUnit,iostat=ierror)
  deallocate(VecRead)
  deallocate(VecWrite)
end program main
