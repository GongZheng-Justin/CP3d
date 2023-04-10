!******************************************************************
! Main program
!******************************************************************
program main
  implicit none
  integer(8),parameter::MemoryLimit=21*1024*1024*1024 ! 6 Gb
  character(len=128),parameter::FileIn='RestartForOCT0550M0000000000'
  
  
  ! locals
  character(len=128)::FileOut
  character(len=1)::SingleStr
  integer::ierror,nUnitIn,nUnitOut,k
  character(len=1),dimension(:),allocatable::ReadMat
  integer(8)::FileByteIn,StrByte,DispRead,MatLenght,nLimit,nLeft,nRead
    
  open(newunit=nUnitIn,file=trim(FileIn),form='unformatted',access='stream',status='old',action='read',position='append',iostat=ierror)
  if(ierror/=0) then
    print*,'Cannot openfile:'//trim(FileIn); stop
  endif
  inquire(unit=nUnitIn,Pos=FileByteIn); FileByteIn=FileByteIn-1_8
  rewind(unit=nUnitIn)
  if(FileByteIn==0) then
    print*,'FileByteIn=0, Stop'
    stop
  endif
  StrByte=sizeof(SingleStr)
  if(mod(FileByteIn,StrByte)/=0) then
    print*,'StrByte Wrong: ',StrByte,FileByteIn; stop  
  endif
  MatLenght=FileByteIn/StrByte
  nLimit=MemoryLimit/StrByte
  allocate(ReadMat(nLimit),stat=ierror)
  if(ierror/=0) then
    print*,'too big MemoryLimit. Might need to decrease the MemoryLimit'; stop
  endif
  
  ! Begin to read
  nLeft=MatLenght
  k=0; DispRead=1_8
  do
    k=k+1
    nRead=min(nLeft,nLimit)
    read(unit=nUnitIn,Pos=DispRead)ReadMat(1:nRead)
    ! 
    write(FileOut,'(A,A,I3.3)')trim(adjustl(FileIn)),'_',k
    open(newunit=nUnitOut,file=trim(FileOut),form='unformatted',access='stream',status='replace',action='write',iostat=ierror)
    if(ierror/=0) then
      print*,'Cannot openfile:'//trim(FileOut); stop
    endif
    write(unit=nUnitOut)ReadMat(1:nRead)
    close(nUnitOut,iostat=ierror)
    print*,'write file successfully:',trim(FileOut)
    
    DispRead=DispRead+StrByte*nRead
    nLeft=nLeft-nRead
    if(nLeft==0) exit
  enddo
  
  close(nUnitIn,iostat=ierror)
  deallocate(ReadMat,stat=ierror)
  print*,'************** OK **************'
end program main
