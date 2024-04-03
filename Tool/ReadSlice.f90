! gfortran -cpp -O3 -Wall -ffree-line-length-none ReadSlice.f90 -o reads
program main
  implicit none
  integer,parameter::RKR=4
  integer,parameter::nxc= 18400
  integer,parameter::nyc= 240
  integer,parameter::nzc= 1680
  character(128),parameter::StrRead ='./CFD_D060/VisuForRD060_06_ux_0000036000'
  character(128),parameter::StrWrite='ux_11'
    
  ! locals
  integer::ierror,nUnit
  integer(8)::disp,RealNum
  real(RKR),dimension(:,:,:),allocatable::VecRead
  
  allocate(VecRead(nxc, nyc, nzc))
  open(newunit=nUnit,file=trim(StrRead),form='unformatted',action='read',status='old', access='stream',position='append',iostat=ierror)
  inquire(unit=nUnit,Pos=disp); disp=disp-1_8
  RealNum=disp/int(RKR,8)
  if(RealNum /= int(nxc,8)*int(nyc,8)*int(nzc,8)) then
    print*,'RealNum wrong, stop'; stop
  else
    print*,' Total real num=',RealNum
  endif
  rewind(unit=nUnit,IOSTAT=ierror)
  read(unit=nUnit)VecRead
  close(unit=nUnit,iostat=ierror)

  open(newunit=nUnit,file=trim(StrWrite),form='unformatted',action='write',status='replace',access='stream',iostat=ierror)
  write(unit=nUnit) VecRead(11,:,:)
  close(unit=nUnit,iostat=ierror)
  deallocate(VecRead)
end program main
