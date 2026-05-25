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
    integer::nUnit,ierror,DirLen
    
    if(.not. this%IsInitialize) then
      call this%InitFileSystem()
      this%IsInitialize=.true.
    endif
    this%nFile=0
    DirLen=len(trim(adjustl(DirStr)))
    open(newunit=nUnit,file='FileNameList.txt',status='old',form='formatted',iostat=ierror)
    if(ierror/=0) then
      write(*,*)"GetFileName: Cannot open file: FileNameList.txt"; STOP
    endif
    do
      read(unit=nUnit,fmt="(A)",iostat=ierror)ch
      if(ierror/=0) exit
      ch=trim(adjustl(ch))
      this%nFile=this%nFile+1
      if(this%nFile>this%mFile) call this%ReallocateFileName()
      this%FileName(this%nFile)=ch(DirLen+1:)
    enddo
    close(nUnit)
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
    integer::nUnit,ierror
    open(newunit=nUnit,file="FileNameList.txt",IOSTAT=ierror)
    close(unit=nUnit,status='delete',IOSTAT=ierror)
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
    integer::ierror,nUnit
    open(newunit=nUnit,file=FileStr,form='unformatted',access='stream',status='old',action='read',position='append',iostat=ierror)
    if(ierror/=0) then
       GetFileByte=-1;return
    endif
    inquire(unit=nUnit,Pos=GetFileByte)
    close(unit=nUnit,iostat=ierror)
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
! Distibutes grid points in one dimension
!******************************************************************
subroutine distribute(ndata,nproc,nrank,istart,isize)
  implicit none
  integer(8),intent(in)::ndata
  integer,intent(in)::nproc,nrank
  integer(8),intent(out)::istart,isize

  ! locals
  integer(8)::nl,iend

  isize=ndata/nproc
  nl=nproc-mod(ndata,nproc)
  istart=nrank*isize+1_8+max(nrank-nl,0_8)
  iend=(nrank+1)*isize+max(nrank+1-nl,0_8)  
  isize=iend-istart+1
end subroutine distribute

!******************************************************************
! Main program
! mpif90 -O3 -cpp -Wall -ffree-line-length-none -fbacktrace -g AverageBinaryNew.f90 -o meanBinary
!******************************************************************
program main
  use MPI
  use File_Tools
  implicit none
  integer,parameter::RKF=8
  integer,parameter::nReadTime=8
  character(128),parameter::ReadFromDir="./Tmp/"
  character(128),parameter::FilterStr="SpecX*"
  character(128),parameter::FileWriteStr="AveSpec2D0000015000"
    
  ! locals
  integer,parameter::RKG=8
  character(128)::FileReadStr
  real(RKF),dimension(:),allocatable::DataRead
  real(RKG),dimension(:),allocatable::DataGather
  integer::nrank,nproc,ierror,iFile,FileByteFlag,iErrTmp,nUnit,nUnitW
  integer(8)::FileByte,FileByteSave,nReal,istart,isize,iMat,disp,itime,iRead,ileft
  
  call MPI_INIT(ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
  
  ! 
  if(nrank==0) call WriteFileNameList(ReadFromDir,FilterStr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  call FileSystem%GetFileName(ReadFromDir)
  if(FileSystem%nFile<1) then
    if(nrank==0) then
      print*,"FileSystem%nFile<1, Stop"
      call DeleteFileNameList()
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_FINALIZE(ierror)
    stop
  endif
  FileByteFlag=0
  FileReadStr=trim(adjustl(ReadFromDir))//trim(adjustl(FileSystem%FileName(1)))
  FileByteSave=GetFileByte(FileReadStr)
  if(nrank==0) then
    do iFile=1,FileSystem%nFile
      FileReadStr=trim(adjustl(ReadFromDir))//trim(adjustl(FileSystem%FileName(iFile)))
      print*, FileReadStr
      FileByte=GetFileByte(FileReadStr)
      if(FileByte/=FileByteSave) then
        FileByteFlag=1
      endif
    enddo
    if(FileByteFlag==1) then
      print*,"FileByte not equal, Stop"
    elseif(mod(FileByteSave,RKF)/=0)  then
      print*,"FileByte wrong-1, Stop"
      FileByteFlag=2
    endif
  endif
  call MPI_BCAST(FileByteFlag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  if(FileByteFlag/=0) then
    if(nrank==0) call DeleteFileNameList()
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_FINALIZE(ierror)
    stop
  endif
  nReal=FileByteSave/RKF
  call distribute(nReal,nproc,nrank,istart,isize)
  disp=int(istart-1,8)*int(sizeof(DataRead(1)),8)+1_8
  
  ! Allocate file
  ierror=0
  iMat=isize/int(nReadTime,8)+1_8
  allocate(DataRead(max(iMat,1_8)),   stat=iErrTmp); ierror=ierror+iErrTmp
  allocate(DataGather(max(iMat,1_8)), stat=iErrTmp); ierror=ierror+iErrTmp  
  if(ierror/=0) then
    print*,'Allocate wrong'
    stop
  endif
  call MPI_FILE_OPEN(MPI_COMM_WORLD, FileWriteStr, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
  call MPI_FILE_SET_SIZE(nUnit,0_8,ierror)  ! Guarantee overwriting
  call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
  call MPI_FILE_PREALLOCATE(nUnit,FileByteSave,ierror)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)    
  call MPI_FILE_CLOSE(nUnit,ierror)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  if(nrank==0) then
    if(GetFileByte(trim(FileWriteStr))==FileByteSave) then
      print*,'FileByte right-1'
    else
      print*,"FileByte wrong-2, Stop"
      stop
    endif
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  
  ! Read and write file
  ileft=isize
  open(newunit=nUnitW,file=trim(FileWriteStr),status='old',form='unformatted',access='stream',action='write',IOSTAT=ierror)
  DO itime=1,nReadTime+1
    print*,nrank,itime
    iRead=min(iMat,ileft)
    DataRead=0.0_RKF; DataGather=0.0_RKG
    
    ! Read file
    do iFile=1,FileSystem%nFile
      FileReadStr=trim(adjustl(ReadFromDir))//trim(adjustl(FileSystem%FileName(iFile)))
      open(newunit=nUnit,file=trim(FileReadStr),status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
      read(nUnit,pos=disp,IOSTAT=ierror) DataRead(1:iRead)
      close(nUnit,IOSTAT=ierror)
      DataGather=DataGather+ real(DataRead, RKG)
    enddo
    DataGather=DataGather/real(FileSystem%nFile, RKG)
    DataRead=real(DataGather, RKF)
    
    ! Write file
    write(nUnitW,pos=disp,IOSTAT=ierror) DataRead(1:iRead)
    
    disp=disp+int(iRead,8)*int(sizeof(DataRead(1)),8)
    ileft=ileft-iRead
    if(ileft==0) exit
  ENDDO
  close(nUnitW,IOSTAT=ierror)
  
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  if(nrank==0) then
    if(GetFileByte(trim(FileWriteStr))==FileByteSave) then
      print*,'FileByte right-2'
    else
      print*,"FileByte wrong-3, Stop"
      stop
    endif
    call DeleteFileNameList()
    print*,'AverageBinary ends successfully'
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  call MPI_FINALIZE(ierror)
end program main
