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
  integer,parameter::nyc=288
  integer,parameter::nzc=1400
  real(RKD),parameter::uTau=0.036245_RKD
  logical,parameter::IsOpenChannel=.false.
  character(len=128),parameter::ReadFromDir="./"
  character(len=128),parameter::RunName="CCT0550L_04"
  character(len=128),parameter::StatisticIn="./Profile_CC0550_48P_4P.txt"
  
  ! locals
  logical::ExistFlag
  character(2)::VarStr
  integer(8)::FileByte
  character(128)::FileName 
  integer::i,nUnit,ierror,nLen,ic,jc,kc,jt,numAve
  real(RKS),allocatable,dimension(:,:,:)::uxCell,uyCell,ArrTmp
  real(RKD)::ReadMeanVec(4,nyc),SumVec(12,nyc),SumVecT(12,nyc)

  ! Read mean statistic data   
  open(newunit=nUnit,file=StatisticIn,form='formatted',status='old',action='read',position='rewind',iostat=ierror)
  if(ierror/=0) then
    print*, "Cannot open file "//trim(StatisticIn); stop
  endif
  read(nUnit,*)
  if(IsOpenChannel) then
    do jc=1,nyc
      read(nUnit,*,iostat=ierror) ReadMeanVec(:,jc)
      if(ierror/=0) then
        print*, "Read file wrong -1"; stop
      endif
    enddo
  else
    do jc=1,nyc/2
      read(nUnit,*,iostat=ierror) ReadMeanVec(:,jc)
      if(ierror/=0) then
        print*, "Read file wrong -2"; stop
      endif
    enddo
    do jc=nyc/2+1,nyc
      jt=nyc+1-jc
      ReadMeanVec(:,jc)=ReadMeanVec(:,jt)
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
    if(VarStr=="ux") then    
      FileName(nLen-12:nLen-11)="uy"
      FileName=trim(adjustl(ReadFromDir))//trim(adjustl(FileName))
      inquire(file=FileName,exist=ExistFlag)    
      if(ExistFlag)numAve=numAve+1
    endif 
  enddo
  print*, 'nyc=   ',nyc
  print*, 'numAve=',numAve
    
  ! calculate statistics
  SumVec=0.0_RKD
  allocate(uxCell(nxc,nzc,nyc))
  allocate(uyCell(nxc,nzc,nyc))
  do i=1,FileSystem%nFile
    FileName=FileSystem%FileName(i)
    nLen=len(trim(adjustl(FileName)))
    VarStr=FileName(nLen-12:nLen-11)
    if(VarStr/="ux") cycle
    FileName(nLen-12:nLen-11)="uy"
    FileName=trim(adjustl(ReadFromDir))//trim(adjustl(FileName))
    inquire(file=FileName,exist=ExistFlag)     
    if(.not.ExistFlag) cycle
    allocate(ArrTmp(nxc,nyc,nzc))

    ! Read uy
    FileByte=GetFileByte(FileName)
    if(FileByte/=sizeof(ArrTmp)) then
      print*,"FileByte wrong-1!"; stop
    endif
    open(newunit=nUnit,file=FileName,form='unformatted',access='stream',status='old',action='read',position='rewind',iostat=ierror)
    read(unit=nUnit)ArrTmp
    call clc_uy_cell(ArrTmp,uyCell,nxc,nyc,nzc)
    close(unit=nUnit,iostat=ierror)
    print*, "Read file "//trim(FileName)//" successfully"
    
    ! Read ux
    FileName=trim(adjustl(ReadFromDir))//trim(adjustl(FileSystem%FileName(i)))
    FileByte=GetFileByte(FileName)
    if(FileByte/=sizeof(ArrTmp)) then
      print*,"FileByte wrong-2 !"; stop
    endif
    open(newunit=nUnit,file=FileName,form='unformatted',access='stream',status='old',action='read',position='rewind',iostat=ierror)
    read(unit=nUnit)ArrTmp
    close(unit=nUnit,iostat=ierror)    
    print*, "Read file "//trim(FileName)//" successfully"  
    call clc_ux_cell(ArrTmp,uxCell,nxc,nyc,nzc)
    deallocate(ArrTmp)
      
    call clc_SumVec(SumVecT)
    SumVec=SumVec+SumVecT
  enddo
  call DeleteFileNameList()
  
  SumVec=SumVec/real(numAve,RKD)
  SumVec(1:2,:)=SumVec(1:2,:)/uTau
  SumVec(3:8,:)=SumVec(3:8,:)/(uTau**2)
  SumVec(3,:)=sqrt(SumVec(3,:)-SumVec(1,:)*SumVec(1,:))
  SumVec(4,:)=sqrt(SumVec(4,:)-SumVec(2,:)*SumVec(2,:))
  
  ! Write statistics
  FileName='ReynoldsFor'//trim(RunName)//".txt"
  open(newunit=nUnit,file=FileName,form='formatted',status='replace',action='write',iostat=ierror)
  if(ierror/=0) then
    print*, "Cannot open file "//trim(StatisticIn); stop
  endif
  write(nUnit,"(A)")" ybar, yplus, u, v, uRms, vRms, uv1, uv2, uv3, uv4, Rat1, Rat2, Rat3, Rat4"
  if(IsOpenChannel) then
    do jc=1,nyc
      write(nUnit,"(14ES24.15)")ReadMeanVec(1,jc),ReadMeanVec(2,jc),SumVec(:,jc)
    enddo
  else
    do jc=1,nyc/2
      jt=nyc+1-jc
      write(nUnit,"(14ES24.15)")ReadMeanVec(1,jc), ReadMeanVec(2,jc), &
                               (SumVec( 1,jc)+SumVec( 1,jt))*0.5_RKD, &
                               (SumVec( 2,jc)-SumVec( 2,jt))*0.5_RKD, &
                               (SumVec( 3,jc)+SumVec( 3,jt))*0.5_RKD, &                                                               
                               (SumVec( 4,jc)+SumVec( 4,jt))*0.5_RKD, &
                               (SumVec( 5,jc)-SumVec( 8,jt))*0.5_RKD, &
                               (SumVec( 6,jc)-SumVec( 7,jt))*0.5_RKD, &
                               (SumVec( 7,jc)-SumVec( 6,jt))*0.5_RKD, &
                               (SumVec( 8,jc)-SumVec( 5,jt))*0.5_RKD, &
                               (SumVec( 9,jc)+SumVec(12,jt))*0.5_RKD, &
                               (SumVec(10,jc)+SumVec(11,jt))*0.5_RKD, &
                               (SumVec(11,jc)+SumVec(10,jt))*0.5_RKD, &
                               (SumVec(12,jc)+SumVec( 9,jt))*0.5_RKD
    enddo  
  endif
  close(unit=nUnit,iostat=ierror)
contains

  subroutine clc_SumVec(SumVecOut)
    implicit none
    real(RKD),intent(out)::SumVecOut(12,nyc)
    
    ! locals
    real(RKD)::uxMean,uyMean,uxValue,uyValue,inxz  
    
    SumVecOut=0.0_RKD
    do jc=1,nyc
      uxMean=ReadMeanVec(3,jc)*uTau
      uyMean=ReadMeanVec(4,jc)*uTau
      do kc=1,nzc
       do ic=1,nxc
         uxValue=uxCell(ic,kc,jc)
         uyValue=uyCell(ic,kc,jc)
         SumVecOut(1,jc)=SumVecOut(1,jc)+uxValue
         SumVecOut(2,jc)=SumVecOut(2,jc)+uyValue
         SumVecOut(3,jc)=SumVecOut(3,jc)+uxValue*uxValue
         SumVecOut(4,jc)=SumVecOut(4,jc)+uyValue*uyValue
       enddo
      enddo
      do kc=1,nzc
       do ic=1,nxc
         uxValue=uxCell(ic,kc,jc)-uxMean
         uyValue=uyCell(ic,kc,jc)-uyMean
         if(uxValue>0.0_RKD      .and. uyValue>0.0_RKD) then
           SumVecOut(5,jc)=SumVecOut(5,jc)+uxValue*uyValue
           SumVecOut(9,jc)=SumVecOut(9,jc)+1.0_RKD
         elseif(uxValue<=0.0_RKD .and. uyValue>0.0_RKD) then
           SumVecOut(6,jc) =SumVecOut(6,jc) +uxValue*uyValue
           SumVecOut(10,jc)=SumVecOut(10,jc)+1.0_RKD
         elseif(uxValue<=0.0_RKD .and. uyValue<=0.0_RKD) then
           SumVecOut(7,jc) =SumVecOut(7,jc) +uxValue*uyValue
           SumVecOut(11,jc)=SumVecOut(11,jc)+1.0_RKD
         else
           SumVecOut(8,jc) =SumVecOut(8,jc)+uxValue*uyValue
           SumVecOut(12,jc)=SumVecOut(12,jc)+1.0_RKD                 
         endif
       enddo
      enddo
    enddo
    inxz=1.0_RKD/real(nxc,RKD)/real(nzc,RKD)
    SumVecOut=SumVecOut*inxz
  end subroutine clc_SumVec
end program main

!******************************************************************
! clc_ux_cell
!******************************************************************
subroutine clc_ux_cell(ArrIn,ArrOut,nxc,nyc,nzc)
  implicit none
  integer::nxc,nyc,nzc
  real(4),dimension(nxc,nyc,nzc),intent(in)::ArrIn
  real(4),dimension(nxc,nzc,nyc),intent(out)::ArrOut

  ! locals
  integer::ic,jc,kc
  real(8),allocatable,dimension(:)::VecTmp
  real(8)::RealTmp,InterpCoe1,InterpCoe2,InterpCoe3,InterpCoe4
  
  allocate(VecTmp(0:nxc+2))
  InterpCoe1= -1.0_8/16.0_8; InterpCoe2=  9.0_8/16.0_8; InterpCoe3=  9.0_8/16.0_8; InterpCoe4= -1.0_8/16.0_8
  do kc=1,nzc
    do jc=1,nyc
      do ic=1,nxc
        VecTmp(ic)=real(ArrIn(ic,jc,kc),8)
      enddo
      VecTmp(0) =VecTmp(nxc)
      VecTmp(nxc+1)=VecTmp(1)
      VecTmp(nxc+2)=VecTmp(2)
      do ic=1,nxc
        RealTmp=InterpCoe1*VecTmp(ic-1) +InterpCoe2*VecTmp(ic) +InterpCoe3*VecTmp(ic+1) +InterpCoe4*VecTmp(ic+2)
        ArrOut(ic,kc,jc)=real(RealTmp,4)
      enddo
    enddo
  enddo
  deallocate(VecTmp)
end subroutine clc_ux_cell

!******************************************************************
! clc_uy_cell
!******************************************************************
subroutine clc_uy_cell(ArrIn,ArrOut,nxc,nyc,nzc)
  implicit none
  integer::nxc,nyc,nzc
  real(4),dimension(nxc,nyc,nzc),intent(in)::ArrIn
  real(4),dimension(nxc,nzc,nyc),intent(out)::ArrOut

  ! locals
  real(8)::RealTmp
  integer::ic,jc,kc  
  real(8),allocatable,dimension(:)::VecTmp
  
  allocate(VecTmp(nyc+1))
  do kc=1,nzc
    do ic=1,nxc
     do jc=1,nyc
       VecTmp(jc)=real(ArrIn(ic,jc,kc),8)
     enddo
     VecTmp(nyc+1)=0.0_8
     do jc=1,nyc
       RealTmp=0.5_8*(VecTmp(jc)+VecTmp(jc+1))
       ArrOut(ic,kc,jc)=real(RealTmp,4)
     enddo
    enddo
  enddo
  deallocate(VecTmp)
end subroutine clc_uy_cell
