module File_Tools
  implicit none
  private
  type FileSystemType
    integer::nFile
    integer::mFile
    logical::IsInitialize=.false.
    character(len=128),dimension(:),allocatable::FileName
  contains
    procedure:: GetFileName
    procedure,private::ReallocateFileName
    procedure,private::InitFileSystem
  end type FileSystemType
  type(FileSystemType),public::FileSystem

  public::GetFileByte,CreateDir,WriteFileNameList,DeleteFileNameList
contains

  !******************************************************************
  ! InitFileSystem
  !******************************************************************
  subroutine InitFileSystem(this)
    implicit none
    class(FileSystemType)::this

    this%nFile=0;this%mFile=10000
    allocate(this%FileName(this%mFile))
  end subroutine InitFileSystem

  !******************************************************************
  ! GetFileName
  !******************************************************************
  subroutine GetFileName(this,DirStr)
    implicit none
    class(FileSystemType)::this
    character(*),intent(in)::DirStr

    ! locals
    character(len=128)::ch
    integer::nUnitRead,ierror,DirLen
    
    if(.not. this%IsInitialize) then
      call this%InitFileSystem()
      this%IsInitialize=.true.
    endif
    this%nFile=0
    DirLen=len(trim(adjustl(DirStr)))
    open(newunit=nUnitRead,file='FileNameList.txt',status='old',form='formatted',iostat=ierror)
    if(ierror/=0) then
      write(*,*)"GetFileName: Cannot open file: FileNameList.txt"; STOP
    endif
    do
      read(unit=nUnitRead,fmt="(A)",iostat=ierror)ch
      if(ierror/=0) exit
      ch=trim(adjustl(ch))
      this%nFile=this%nFile+1
      if(this%nFile>this%mFile) call this%ReallocateFileName()
      this%FileName(this%nFile)=ch(DirLen+1:)
    enddo
    close(nUnitRead)
  end subroutine GetFileName

  !******************************************************************
  ! WriteFileNameList
  !******************************************************************
  subroutine WriteFileNameList(DirStr,FilterStr)
    implicit none
    character(*),intent(in)::DirStr,FilterStr
    call system('ls '//trim(adjustl(DirStr))//trim(adjustl(FilterStr))//' 2>/dev/null >FileNameList.txt')
  end subroutine WriteFileNameList

  !******************************************************************
  ! DeleteFileNameList
  !******************************************************************
  subroutine DeleteFileNameList()
    implicit none
    call system('rm FileNameList.txt 2>/dev/null')
  end subroutine DeleteFileNameList

  !******************************************************************
  ! ReallocateFileName
  !******************************************************************
  subroutine ReallocateFileName(this)
    implicit none
    class(FileSystemType)::this

    ! locals
    integer::sizep,sizen
    character(len=128),dimension(:),allocatable::NameTemp
    
    sizep=this%mFile
    sizen=int(real(sizep)*1.2) +1

    call move_alloc(this%FileName,NameTemp)
    allocate(this%FileName(sizen))
    this%FileName(1:sizep)=NameTemp
    this%mFile=sizen
  end subroutine ReallocateFileName

  !******************************************************************
  ! GetFileByte
  !******************************************************************
  integer(kind=8) function GetFileByte(FileStr)
    implicit none
    character(*),intent(in)::FileStr

    ! locals
    integer::ierror,nUnitRead
    open(newunit=nUnitRead,file=FileStr,form='unformatted',access='stream',status='old',action='read',position='append',iostat=ierror)
    if(ierror/=0) then
       GetFileByte=-1;return
    endif
    inquire(unit=nUnitRead,Pos=GetFileByte)
    close(unit=nUnitRead,iostat=ierror)
    GetFileByte=GetFileByte-1
  end function GetFileByte

  !******************************************************************
  ! CreateDir
  !******************************************************************
  subroutine CreateDir(clearFlag,DirStr)
    implicit none
    logical,intent(in)::clearFlag
    character(*),intent(in)::DirStr

    if(clearFlag) then
      call system("rm -rf "//trim(adjustl(DirStr))//" 2>/dev/null")
      print*,"delete directory: "//trim(adjustl(DirStr))//"  successfully"
    endif
    call system("mkdir "//trim(adjustl(DirStr))//" 2>/dev/null")
  end subroutine CreateDir
end module File_Tools

!******************************************************************
! Main program
! gfortran -O3 -cpp -Wall -ffree-line-length-none -fbacktrace -g AvePlane2Plane.f90 -o AveP2P
!******************************************************************
program main
  use File_Tools
  implicit none
  integer,parameter::RKF=4
  integer,parameter::nxc=22400
  integer,parameter::nzc=1120
  integer,parameter::nTime=200
  integer,parameter::nDelta=14
  character(len=128),parameter::ReadFromDir="./PlaneData/"
  character(len=128),parameter::WriteToDir ="./AvePlane2Info/"
  
  ! locals
  real(8)::rAve
  integer(8)::OneDumpByte,FileByte
  real(RKF),dimension(:,:,:),allocatable::uxp
  integer,dimension(:),allocatable::ipc,iave,kave
  real(8),dimension(:,:,:),allocatable::uxc,uxa
  character(128)::FileStrRead,FileStrWrite,StrTmp
  integer::nxa,nza,i,nUnitRead,ierror,ic,kc,itime,ind,knd,nUnitWrite
  
  if(mod(nxc,nDelta)/=0 .or. mod(nzc,nDelta)/=0) then
    print*,'nxc OR nzc wrong'; stop
  endif
  nxa= nxc/nDelta; nza=nzc/nDelta
  
  allocate(uxp(nxc,nzc,nTime) )
  allocate(uxc(nxc,nzc,nTime) )
  allocate(uxa(nxa,nza,nTime) )
  OneDumpByte=sizeof(uxp)
  
  allocate(ipc(nxc))
  do ic=1,nxc-1
    ipc(ic)=ic+1
  enddo
  ipc(nxc)=1
  
  allocate(iave(nxc))
  do ic=1,nxc
    if(mod(ic,nDelta)==0) then
      iave(ic)= ic/nDelta
    else
      iave(ic)= ic/nDelta+1
    endif
  enddo
  allocate(kave(nzc))
  do kc=1,nzc
    if(mod(kc,nDelta)==0) then
      kave(kc)= kc/nDelta
    else
      kave(kc)= kc/nDelta +1
    endif
  enddo
  rAve=1.0/real(nDelta,8)/real(nDelta,8)
  
  ! Read file name
  call WriteFileNameList(ReadFromDir,"uxcPlane_y02_*")
  call FileSystem%GetFileName(ReadFromDir)
  do i=1,FileSystem%nFile  
    write(FileStrRead,'(A,A)')trim(ReadFromDir),trim(FileSystem%FileName(i))
    FileByte=GetFileByte(trim(FileStrRead))
    if(FileByte /= OneDumpByte) then
      print*,'FileByte Error-1: ',trim(FileStrRead); stop
    endif
    
    open(newunit=nUnitRead, file=trim(FileStrRead), status='old',    form='unformatted',access='stream',action='read', IOSTAT=ierror)
    read(unit=nUnitRead,IOSTAT=ierror) uxp
    close(unit=nUnitRead, IOSTAT=ierror)
    do itime=1,nTime
      do kc=1,nzc
        do ic=1,nxc
          ind = ipc(ic)
          uxc(ic,kc,itime)= (real(uxp(ic, kc,itime), 8)+real(uxp(ind, kc,itime), 8))*0.5_8
        enddo
      enddo
    enddo

    uxa=0.0_8
    StrTmp=FileSystem%FileName(i)
    ind=len_trim(StrTmp)
    write(FileStrWrite,'(A,A,A)')trim(WriteToDir),'AP_uxcPlane_y02_',StrTmp(ind-9:ind)    
    do itime=1,nTime
      do kc=1,nzc
        knd=kave(kc)
        do ic=1,nxc
          ind=iave(ic)
          uxa(ind,knd,itime)= uxa(ind,knd,itime)+ uxc(ic,kc,itime)
        enddo
      enddo
    enddo
    uxa= rAve*uxa
    open(newunit=nUnitWrite,file=trim(FileStrWrite),status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
    write(unit=nUnitWrite,IOSTAT=ierror) real(uxa, RKF)
    close(unit=nUnitWrite,IOSTAT=ierror)
    print*,'Average plane from ',trim(FileStrRead),' to ',trim(FileStrWrite),' successfully!'
  enddo
  call DeleteFileNameList()
  deallocate(ipc,iave,kave,uxp,uxc,uxa)
    
end program main
