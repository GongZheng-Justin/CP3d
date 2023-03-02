program main
  implicit none
  integer,parameter::nxc=9216
  integer,parameter::nyc=144
  integer,parameter::nzc=1400
  real(8),dimension(:,:,:),allocatable::ux1
  real(4),dimension(:,:,:),allocatable::ux2
  integer::nUnit,ierror,ic,jc,kc
  
  
  allocate(ux1(nxc,nyc,nzc))
  allocate(ux2(nxc,nzc,nyc))
  open(newunit=nUnit,file='RestartForB05_010000033000',form='unformatted',access='stream',status='old',action='read',iostat=ierror)
  read(unit=nUnit)ux1
  close(unit=nUnit,iostat=ierror)
  do kc=1,nzc
    do jc=1,nyc
      do ic=1,nxc
        ux2(ic,kc,jc)=real(ux1(ic,jc,kc),4) 
      enddo
    enddo
  enddo  
  open(newunit=nUnit,file='VisuB05',action='write',form='unformatted',status='replace',access='stream',iostat=ierror)
  write(unit=nUnit)ux2
  close(unit=nUnit,iostat=ierror)  
end program main
