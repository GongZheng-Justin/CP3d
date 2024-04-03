module PPD_Tools
  implicit none
  private

  interface QuickSort
    module procedure QuickSort_int,QuickSort_real4,QuickSort_real8
  end interface QuickSort

  type FileSystemType
    integer::nFile
    integer::mFile
    character(len=128),dimension(:),allocatable::FileName
  contains
    procedure:: Initialize => InitFileSystem
    procedure:: GetFileName
    procedure:: WriteFileNameList
    procedure:: DeleteFileNameList
    procedure,private:: ReallocateFileName
  end type FileSystemType

  public:: QuickSort,GetFileByte,CreateDir
  type(FileSystemType),public::FileSystem
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
    integer::myistat,DirLen

  
    this%nFile=0
    DirLen=len(trim(adjustl(DirStr)))
    open(25,file='FileNameList.txt',status='old',form='formatted')
    do
      read(25,"(A)",IOSTAT=myistat)ch
      if(myistat/=0) exit
      ch=trim(adjustl(ch))
      this%nFile=this%nFile+1
      if(this%nFile>this%mFile) call this%ReallocateFileName()
      this%FileName(this%nFile)=ch(DirLen+1:)
    enddo
    close(25)
  end subroutine GetFileName

  !******************************************************************
  ! WriteFileNameList
  !******************************************************************
  subroutine WriteFileNameList(this,DirStr,FilterStr)
    implicit none
    class(FileSystemType)::this
    character(*),intent(in)::DirStr,FilterStr
    call system('ls '//trim(adjustl(DirStr))//trim(adjustl(FilterStr))//' 2>/dev/null >FileNameList.txt')
  end subroutine WriteFileNameList

  !******************************************************************
  ! DeleteFileNameList
  !******************************************************************
  subroutine DeleteFileNameList(this)
    implicit none
    class(FileSystemType)::this

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
  ! QuickSort_int
  !******************************************************************
  recursive subroutine QuickSort_int(array,seque,first,last)
    implicit none
    integer,intent(in)::first,last
    integer,dimension(:),intent(inout)::array
    integer,dimension(:),intent(inout)::seque

    ! locals
    integer::pivot,tempA
	  integer::i,j,middle,tempS

    if(first>=last) return
    middle=(last+first)/2

	  pivot=array(middle)
    tempA=array(middle);array(middle)=array(last);array(last)=tempA
    tempS=seque(middle);seque(middle)=seque(last);seque(last)=tempS
	  i=first
    do j=first,last-1
      if(array(j)<=pivot) then
        tempA=array(i);array(i)=array(j);array(j)=tempA
        tempS=seque(i);seque(i)=seque(j);seque(j)=tempS
        i=i+1
      endif
    enddo
    tempA=array(i);array(i)=array(last);array(last)=tempA
    tempS=seque(i);seque(i)=seque(last);seque(last)=tempS
 
    if(i>first) call QuickSort_int(array,seque,first,i-1)
    call QuickSort_int(array,seque,i+1,last)
  end subroutine QuickSort_int

  !******************************************************************
  ! QuickSort_real4
  !******************************************************************
  recursive subroutine QuickSort_real4(array,seque,first,last)
    implicit none
    integer,intent(in)::first,last
    real(4),dimension(:),intent(inout)::array
    integer,dimension(:),intent(inout)::seque

    ! locals
    real(4)::pivot,tempA
  	integer::i,j,middle,tempS

    if(first>=last) return
    middle=(last+first)/2

	  pivot=array(middle)
    tempA=array(middle);array(middle)=array(last);array(last)=tempA
    tempS=seque(middle);seque(middle)=seque(last);seque(last)=tempS
	  i=first
    do j=first,last-1
      if(array(j)<=pivot) then
        tempA=array(i);array(i)=array(j);array(j)=tempA
        tempS=seque(i);seque(i)=seque(j);seque(j)=tempS
        i=i+1
      endif
    enddo
    tempA=array(i);array(i)=array(last);array(last)=tempA
    tempS=seque(i);seque(i)=seque(last);seque(last)=tempS
 
    if(i>first) call QuickSort_real4(array,seque,first,i-1)
    call QuickSort_real4(array,seque,i+1,last)
  end subroutine QuickSort_real4

  !******************************************************************
  ! QuickSort_real8
  !******************************************************************
  recursive subroutine QuickSort_real8(array,seque,first,last)
    implicit none
    integer,intent(in)::first,last
    real(8),dimension(:),intent(inout)::array
    integer,dimension(:),intent(inout)::seque

    ! locals
    real(8)::pivot,tempA
	  integer::i,j,middle,tempS

    if(first>=last) return
    middle=(last+first)/2

  	pivot=array(middle)
    tempA=array(middle);array(middle)=array(last);array(last)=tempA
    tempS=seque(middle);seque(middle)=seque(last);seque(last)=tempS
	  i=first
    do j=first,last-1
      if(array(j)<=pivot) then
        tempA=array(i);array(i)=array(j);array(j)=tempA
        tempS=seque(i);seque(i)=seque(j);seque(j)=tempS
        i=i+1
      endif
    enddo
    tempA=array(i);array(i)=array(last);array(last)=tempA
    tempS=seque(i);seque(i)=seque(last);seque(last)=tempS
 
    if(i>first) call QuickSort_real8(array,seque,first,i-1)
    call QuickSort_real8(array,seque,i+1,last)
  end subroutine QuickSort_real8

  !******************************************************************
  ! GetFileByte
  !******************************************************************
  integer function GetFileByte(FileStr)
    implicit none
    character(*),intent(in)::FileStr

    ! locals
    integer::ierror
    open(unit=50,file=FileStr,form='unformatted',access='stream',status='old',action='read',position='append',iostat=ierror)
    if(ierror/=0) then
       GetFileByte=-1;return
    endif
    inquire(unit=50,Pos=GetFileByte)
    close(unit=50)
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

end module PPD_Tools
