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
! Main program
!******************************************************************
program main
  use File_Tools
  implicit none
  integer,parameter::RKS=4
  integer,parameter::RKD=8
  integer,parameter::nxc=9216
  integer,parameter::nyc=144
  integer,parameter::nzc=1400
  real(RKD),parameter::uTau=0.036245_RKD
  logical,parameter::IsOpenChannel=.true.
  character(len=128),parameter::ReadFromDir="./FieldOC0550_48P_4P/"
  character(len=128),parameter::RunName="OCT0550L_02"
  character(len=128),parameter::StatisticIn="./Profile_OC0550_48P_4P.txt"
  
  ! locals
  character(2)::VarStr
  integer(8)::FileByte
  character(128)::FileName
  integer,parameter::NCHASTAT=4
  real(RKS),allocatable,dimension(:,:,:)::prCell
  integer::i,nUnit,ierror,nLen,ic,jc,kc,jt,numAve
  real(RKD)::ReadMeanPr(13,nyc),SumPr(NCHASTAT,nyc),SumVec(NCHASTAT),rTemp,prMean,inxz

  ! Read mean statistic data   
  open(newunit=nUnit,file=StatisticIn,form='formatted',status='old',action='read',position='rewind',iostat=ierror)
  if(ierror/=0) then
    print*, "Cannot open file "//trim(StatisticIn); stop
  endif
  read(nUnit,*)
  if(IsOpenChannel) then
    do jc=1,nyc
      read(nUnit,*,iostat=ierror) ReadMeanPr(:,jc)
      if(ierror/=0) then
        print*, "Read file wrong -1"; stop
      endif
    enddo
  else
    do jc=1,nyc/2
      read(nUnit,*,iostat=ierror) ReadMeanPr(:,jc)
      if(ierror/=0) then
        print*, "Read file wrong -2"; stop
      endif
    enddo
    do jc=nyc/2+1,nyc
      jt=nyc+1-jc
      ReadMeanPr(:,jc)=ReadMeanPr(:,jt)
    enddo
  endif
  close(unit=nUnit,iostat=ierror)

  ! determine the file number for average
  call WriteFileNameList(ReadFromDir,"*"//trim(RunName)//"*")
  call FileSystem%GetFileName(ReadFromDir)
  numAve=0
  do i=1,FileSystem%nFile
    FileName=FileSystem%FileName(i)
    nLen=len(trim(adjustl(FileName)))
    VarStr=FileName(nLen-12:nLen-11)
    if(VarStr=="pr") numAve=numAve+1
  enddo
  print*, 'nyc=   ',nyc
  print*, 'numAve=',numAve
    
  ! calculate statistics
  SumPr=0.0_RKD
  allocate(prCell(nxc,nyc,nzc))
  inxz=1.0_RKD/real(nxc,RKD)/real(nzc,RKD)
  do i=1,FileSystem%nFile
    FileName=FileSystem%FileName(i)
    nLen=len(trim(adjustl(FileName)))
    VarStr=FileName(nLen-12:nLen-11)
    if(VarStr/="pr") cycle
    FileName=trim(adjustl(ReadFromDir))//trim(adjustl(FileName))

    ! Read pr
    FileByte=GetFileByte(FileName)
    if(FileByte/=sizeof(prCell)) then
      print*,"FileByte wrong-1!"; stop
    endif
    open(newunit=nUnit,file=FileName,form='unformatted',access='stream',status='old',action='read',position='rewind',iostat=ierror)
    read(unit=nUnit)prCell
    close(unit=nUnit,iostat=ierror)
    DO jc=1,nyc
      SumVec=0.0_RKD
      prMean=ReadMeanPr(6,jc)
      do kc=1,nzc
        do ic=1,nxc
          rTemp=prCell(ic,jc,kc)/(uTau*uTau)
          SumVec(1)=SumVec(1)+rTemp
          rTemp=rTemp-prMean
          SumVec(2)=SumVec(2)+rTemp*rTemp
          SumVec(3)=SumVec(3)+rTemp*rTemp*rTemp
          SumVec(4)=SumVec(4)+rTemp*rTemp*rTemp*rTemp
        enddo
      enddo
      do kc=1,NCHASTAT
        SumPr(kc,jc)=SumPr(kc,jc)+SumVec(kc)*inxz
      enddo
    ENDDO
    print*, "Read file "//trim(FileName)//" successfully"
  enddo
  call DeleteFileNameList()
  SumPr=SumPr/real(numAve,RKD)
  do jc=1,nyc
    rTemp=ReadMeanPr(13,jc)
    SumPr(2,jc)=sqrt(SumPr(2,jc))
    SumPr(3,jc)=SumPr(3,jc)/(rTemp**3)
    SumPr(4,jc)=SumPr(4,jc)/(rTemp**4)
  enddo
  
  ! Write statistics
  FileName='PrStat_For'//trim(RunName)//".txt"
  open(newunit=nUnit,file=FileName,form='formatted',status='replace',action='write',iostat=ierror)
  if(ierror/=0) then
    print*, "Cannot open file "//trim(StatisticIn); stop
  endif
  write(nUnit,"(A)")" ybar, yplus, pr, prRms, S(pr), F(pr)"
  if(IsOpenChannel) then
    do jc=1,nyc
      write(nUnit,"(14ES24.15)")ReadMeanPr(1,jc),ReadMeanPr(2,jc),SumPr(:,jc)
    enddo
  else
    do jc=1,nyc/2
      jt=nyc+1-jc
      write(nUnit,"(14ES24.15)")ReadMeanPr(1,jc), ReadMeanPr(2,jc), &
                               (SumPr( 1,jc)+SumPr( 1,jt))*0.5_RKD, &
                               (SumPr( 2,jc)+SumPr( 2,jt))*0.5_RKD, &
                               (SumPr( 3,jc)+SumPr( 3,jt))*0.5_RKD, &
                               (SumPr( 4,jc)+SumPr( 4,jt))*0.5_RKD
    enddo  
  endif
  close(unit=nUnit,iostat=ierror)
end program main
