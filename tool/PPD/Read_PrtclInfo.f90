module File_Tools
  implicit none
  private

  public::GetFileByte
contains

  !******************************************************************
  ! GetFileByte
  !******************************************************************
  integer(kind=8) function GetFileByte(FileStr)
    implicit none
    character(*),intent(in)::FileStr

    ! locals
    integer::ierror,nUnit
    open(newunit=nUnit,file=FileStr,form='unformatted',access='stream',status='old',action='read',position='append',iostat=ierror)
    if(ierror/=0) then
       GetFileByte=-1;return
    endif
    inquire(unit=nUnit,Pos=GetFileByte)
    close(unit=nUnit,iostat=ierror)
    GetFileByte=GetFileByte-1
  end function GetFileByte
end module File_Tools

!******************************************************************
! Main program
! gfortran -O3 -cpp -Wall -ffree-line-length-none -fbacktrace -g Read_PrtclInfo.f90 -o Read_PrtclInfo
!******************************************************************
#define nDumpPrtclInte 4
#define nDumpPrtclReal 21

#define iDump_itime    1
#define iDump_id       2
#define iDump_IsCntct  4
program main
  use File_Tools
  implicit none
  integer,parameter::IKF=4
  integer,parameter::RKF=4
  integer,parameter::numPrtcl=305000
  integer,parameter::DumpPrtclFreq=600
  integer,parameter::iTimeSet=120000/20                             ! Need to modify-1
  integer,parameter::iTimeEnd=124000/20                             ! Need to modify-2
  character(128),parameter::FileStrRead="PrtclInfo_0000124000"     ! Need to modify-3
  character(128),parameter::FileStrWrite="PrtclArrange_0000124000" ! Need to modify-4
    
  ! locals
  real(8)::ctime1,ctime2
  integer(8),  allocatable,dimension(:) ::seque
  integer(IKF),dimension(nDumpPrtclInte,2)::InteTmp
  real(RKF),   dimension(nDumpPrtclReal,2)::RealTmp
  integer(IKF),allocatable,dimension(:,:)::InteArrRead
  real(RKF),   allocatable,dimension(:,:)::RealArrRead
  integer(8)::FileByte,nDump,disp,nDumpS,nDumpE,ind,iprev,inext
  integer::OneDumpByte,nUnit,ierror,mproc,iDump,pid,itime,iRearrange
  
  if(mod(nDumpPrtclReal,3)/=0 ) then
    print*,'nDumpPrtclReal wrong'; stop
  endif
  if(trim(adjustl(FileStrWrite))==trim(adjustl(FileStrRead))) then
    print*,'FileStrWrite wrong'; stop
  endif
  FileByte=GetFileByte(FileStrRead)
  print*,'FileByte=',FileByte
  OneDumpByte=IKF*nDumpPrtclInte+RKF*nDumpPrtclReal
  
  ! Prepare
  call cpu_time(ctime1)
  nDump=0_8; disp=1_8; mproc=0
  open(newunit=nUnit,file=trim(FileStrRead),status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
  do 
    mproc=mproc+1
    read(nUnit,pos=disp,IOSTAT=ierror) iDump
    nDump=nDump+int(iDump,8)
    disp=disp+4_8+int(OneDumpByte,8)*int(iDump,8)
    if(disp==FileByte+1_8) exit
  enddo
  disp=int(OneDumpByte,8)*int(nDump,8)+int(mproc,8)*IKF
  if(disp/=FileByte) then
    print*,'File byte wrong-1'; stop
  endif
  allocate(InteArrRead(nDumpPrtclInte,nDump))
  allocate(RealArrRead(nDumpPrtclReal,nDump))
  
  ! Read
  disp=1_8; nDumpE=0_8
  do ind=1,mproc
    read(nUnit,pos=disp,IOSTAT=ierror) iDump
    disp=disp+int(IKF,8)
    nDumpS=nDumpE+1_8
    nDumpE=nDumpS+int(iDump,8)-1_8
   
    ! Integer
    read(nUnit,pos=disp,IOSTAT=ierror) InteArrRead(:,nDumpS:nDumpE)
    disp=disp+int(iDump,8)*int(IKF,8)*int(nDumpPrtclInte,8)
    
    ! Real
    read(nUnit,pos=disp,IOSTAT=ierror) RealArrRead(:,nDumpS:nDumpE)
    disp=disp+int(iDump,8)*int(RKF,8)*int(nDumpPrtclReal,8)
  enddo
  if(disp /= FileByte+1_8) then
    print*,'File byte wrong-2'; stop
  endif
  if(nDumpE /= nDump) then
    print*,'nDumpE wrong',nDumpE,nDump; stop 
  endif
  close(nUnit,IOSTAT=ierror)
  call cpu_time(ctime2)
  print*,'Read File succeccfully: ',trim(FileStrRead),ctime2-ctime1
  
  ! Seque
  call cpu_time(ctime1)
  allocate(seque(nDump))
  do ind=1,nDump
    pid  =InteArrRead(iDump_id,ind)
    itime=InteArrRead(iDump_itime,ind)
    if(pid>numPrtcl .or. pid<1) then
      print*,'pid wrong-1',pid; stop
    endif
    if(mod(itime,DumpPrtclFreq)/=0) then
      print*,'itime wrong-1',itime; stop
    endif
    itime=itime/DumpPrtclFreq
    if(itime<=iTimeSet .or. itime>iTimeEnd) then
      print*,'itime wrong-2',itime; stop
    endif
    itime=itime-iTimeSet
    seque(ind)=(itime-1_8)*int(numPrtcl,8)+int(pid,8)
  enddo
  
  ! Rearrange
  iRearrange=0
  do ind=1,nDump
    if(seque(ind)==ind) cycle
    iprev=ind
    inext=seque(iprev)
    seque(iprev)=-99_8
       
    InteTmp(:,1)=InteArrRead(:,iprev)
    RealTmp(:,1)=RealArrRead(:,iprev)
    do
      iRearrange=iRearrange+1
      
      InteTmp(:,2)=InteArrRead(:,inext)
      InteArrRead(:,inext)=InteTmp(:,1)
      InteTmp(:,1)=InteTmp(:,2)
      
      RealTmp(:,2)=RealArrRead(:,inext)
      RealArrRead(:,inext)=RealTmp(:,1)
      RealTmp(:,1)=RealTmp(:,2)      
            
      iprev=inext
      inext=seque(iprev)
      seque(iprev)=iprev
      if(inext==-99_8) exit
    enddo
  enddo
  call cpu_time(ctime2)
  print*,'iRearrange=',iRearrange,ctime2-ctime1
    
  ! Check
  call cpu_time(ctime1)
  do ind=1,nDump
    if(seque(ind) /= ind) then
      print*,'seque(ind) wrong',ind,seque(ind); stop
    endif
    pid=int(mod(ind,numPrtcl),4)
    if(pid==0) pid=numPrtcl
    if(InteArrRead(iDump_id,ind) /= pid) then
      print*,'pid wrong-2',ind,pid,InteArrRead(iDump_id,ind); stop    
    endif
    
    itime=InteArrRead(iDump_itime,ind)
    if(mod(itime,DumpPrtclFreq)/=0) then
      print*,'itime wrong-3',itime; stop
    endif
    itime=itime/DumpPrtclFreq
    if(itime<=iTimeSet .or. itime>iTimeEnd) then
      print*,'itime wrong-4',itime; stop
    endif        
    itime=itime-iTimeSet
    
    if(pid==numPrtcl) then
      iRearrange=int(ind/numPrtcl, 4)
    else
      iRearrange=int((ind-int(pid,8))/numPrtcl, 4)+1
    endif
    if(itime /= iRearrange) then
      print*,'itime wrong-5',itime; stop   
    endif
    
    iRearrange=InteArrRead(iDump_IsCntct,ind)
    if(iRearrange/=0 .and. iRearrange/=1) then
      print*,'iDump_IsCntct wrong',InteArrRead(iDump_IsCntct,ind); stop      
    endif
  enddo
  deallocate(seque)
  call cpu_time(ctime2)
  print*,'Pass check',ctime2-ctime1
  
  ! Write
  call cpu_time(ctime1)
  iRearrange=int(mod(nDump,numPrtcl),4)
  if(iRearrange/=0) then
    print*,'numPrtcl wrong',nDump,numPrtcl ; stop
  endif
  iRearrange=int(nDump/numPrtcl,4)
  
  nDumpE=0_8
  iDump = nDumpPrtclReal/3
  open(newunit=nUnit,file=trim(FileStrWrite),status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
  do itime=1,iRearrange
    nDumpS=nDumpE+1_8
    nDumpE=nDumpS+int(numPrtcl,8)-1_8
    do pid =1,iDump
      write(nUnit) RealArrRead(3*pid-2: 3*pid, nDumpS:nDumpE)
    enddo
    write(nUnit) InteArrRead( 1,    nDumpS:nDumpE)
    write(nUnit) InteArrRead( 2,    nDumpS:nDumpE)
    write(nUnit) InteArrRead( 3,    nDumpS:nDumpE)
    write(nUnit) InteArrRead( 4,    nDumpS:nDumpE)
  enddo
  close(nUnit,IOSTAT=ierror)
  call cpu_time(ctime2)
  print*,'Write time: ',ctime2-ctime1
 
  deallocate(InteArrRead)
  deallocate(RealArrRead)
  print*,'Finished successfully'
end program main

#undef nDumpPrtclInte
#undef nDumpPrtclReal

#undef iDump_itime
#undef iDump_id
#undef iDump_IsCntct
