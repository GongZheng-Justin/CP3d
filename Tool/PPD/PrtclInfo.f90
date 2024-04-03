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
! gfortran -O3 -cpp -Wall -ffree-line-length-none -fbacktrace -g PrtclInfo.f90 -o ClcPrtclInfo
!******************************************************************
#define nDumpPrtclInte 4
#define nDumpPrtclReal 21

#define iDump_itime    1
#define iDump_id       2
#define iDump_IsCntct  4
#define PI (3.141592653589793238462643383279502884_8)
program main
  use File_Tools
  implicit none
  integer,parameter::IKF=4
  integer,parameter::RKF=4
  integer,parameter::numPrtcl=1000000
  integer,parameter::inx=1600
  integer,parameter::inz=80
  integer,parameter::nData=2
  real(8),parameter::Gravity = 9.81_8
  real(8),parameter::Diameter=0.001_8
  real(8),parameter::rho_p  =1700.0_8
  real(8),parameter::rho_f  =1000.0_8
  real(8),parameter::Shield_c= 0.05_8
  real(8),parameter::NetWeight_c=1.0E-5_8
  real(8),parameter::DeltaXZ = 0.00125_8
  character(len=128),parameter::ReadFromDir="./"
  character(len=128),parameter::WriteToDir ="./PrtclInfo/"
  
  ! locals
  logical::BedFlag
  integer(8)::FileByte,disp,OneDumpByte
  real(8),dimension(:,:,:),allocatable::InfoArr
  real(8)::u_th,F_th,y_th,u_mag,VolumeDDelta,PosY
  character(128)::FileStrRead,FileStrWrite1,FileStrWrite2,StrTmp
  integer::i,nUnitRead,ierror,nTime,pid,itime,pi,pk,nUnitWrite1,nUnitWrite2
  real(RKF),dimension(:,:),allocatable::GPrtcl_Pos,GPrtcl_linVel,GPrtcl_CntctForce
  
  y_th=8.0*Diameter
  u_th=sqrt(Shield_c)*sqrt((rho_p/rho_f-1.0_8)*Gravity*Diameter)
  F_th=NetWeight_c*(rho_p-rho_f)*PI*(Diameter**3)*Gravity/6.0_8;
  
  allocate(GPrtcl_Pos(3,numPrtcl))
  allocate(GPrtcl_linVel(3,numPrtcl))
  allocate(GPrtcl_CntctForce(3,numPrtcl))
  allocate(InfoArr(inx,inz,nData))
  OneDumpByte=int(IKF*nDumpPrtclInte+RKF*nDumpPrtclReal,8)*int(numPrtcl,8)
  
  ! Read file name
  VolumeDDelta = (PI*Diameter*Diameter*Diameter/6.0_8)/(DeltaXZ*DeltaXZ)
  call WriteFileNameList(ReadFromDir,"PrtclArrange_*")
  call FileSystem%GetFileName(ReadFromDir)
  do i=1,FileSystem%nFile
    write(FileStrRead,'(A,A)')trim(ReadFromDir),trim(FileSystem%FileName(i))
    FileByte=GetFileByte(trim(FileStrRead))
    if(mod(FileByte,OneDumpByte)/=0) then
      print*,'FileByte Error-1: ',trim(FileStrRead); stop
    endif
    nTime=int(FileByte/OneDumpByte,4)
        
    StrTmp=FileSystem%FileName(i)
    itime=len_trim(StrTmp)
    write(FileStrWrite1,'(A,A,A)')trim(WriteToDir),'PI_Vel_',StrTmp(itime-9:itime)
    write(FileStrWrite2,'(A,A,A)')trim(WriteToDir),'PI_Hei_',StrTmp(itime-9:itime)
    open(newunit=nUnitRead, file=trim(FileStrRead), status='old',    form='unformatted',access='stream',action='read', IOSTAT=ierror)
    open(newunit=nUnitWrite1,file=trim(FileStrWrite1),status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
    open(newunit=nUnitWrite2,file=trim(FileStrWrite2),status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
    do itime=1,nTime
      disp=1_8 +int(itime-1,8)*OneDumpByte
    
      ! GPrtcl_Pos
      read(unit=nUnitRead,pos=disp,IOSTAT=ierror) GPrtcl_Pos
      disp=disp+int(RKF*3,8)*int(numPrtcl,8)
      
      ! GPrtcl_linVel
      read(unit=nUnitRead,pos=disp,IOSTAT=ierror) GPrtcl_linVel
      disp=disp+int(RKF*3,8)*int(numPrtcl,8)
    
      ! CntctForce, Skip RotVel, FpForce, FpTorque, (CntctForce, Torque)
      disp=disp+int(RKF*9,8)*int(numPrtcl,8)
      read(unit=nUnitRead,pos=disp,IOSTAT=ierror) GPrtcl_CntctForce
      
      ! Divide the bed to individual cells
      InfoArr=0.0_8
      do pid=1, numPrtcl
        PosY = GPrtcl_Pos(2,pid)
        pi= floor(GPrtcl_Pos(1,pid)/DeltaXZ) +1
        pk= floor(GPrtcl_Pos(3,pid)/DeltaXZ) +1
        if(pi>inx+1 .or. pi<0 .or. pk>inz+1 .or. pk<0) then
          print*,"pi or pk wrong; stop"; stop
        endif
        if(pi>inx) pi=1
        if(pi<1)   pi=inx
        if(pk>inz) pk=1
        if(pk<1)   pk=inz
        InfoArr(pi,pk,1)= InfoArr(pi,pk,1) +GPrtcl_linVel(1,pid)*VolumeDDelta
        
        ! determine whether this particle is bed particle
        BedFlag=.true.
        u_mag=sqrt(GPrtcl_linVel(1,pid)*GPrtcl_linVel(1,pid)+GPrtcl_linVel(2,pid)*GPrtcl_linVel(2,pid) &
                  +GPrtcl_linVel(3,pid)*GPrtcl_linVel(3,pid)) 
        if(u_mag>u_th) BedFlag=.false.
        if(abs(GPrtcl_CntctForce(2,pid))<=F_th .and. PosY>=y_th) BedFlag=.false.
        if(.not. BedFlag) cycle 
        
        PosY = PosY/Diameter + 0.3_8       
        if(InfoArr(pi,pk,2) < PosY) InfoArr(pi,pk,2)=PosY
      enddo
      write(unit=nUnitWrite1,IOSTAT=ierror) real(InfoArr(:,:,1), RKF)
      write(unit=nUnitWrite2,IOSTAT=ierror) real(InfoArr(:,:,2), RKF)
    enddo
    close(unit=nUnitRead, IOSTAT=ierror)
    close(unit=nUnitWrite1,IOSTAT=ierror)
    close(unit=nUnitWrite2,IOSTAT=ierror)
    FileByte=GetFileByte(trim(FileStrWrite1))
    if(FileByte /= int(inx,8)*int(inz,8)*nTime*int(RKF,8)) then
      print*,'FileByte Error-2: ',trim(FileStrWrite1); stop
    endif
    FileByte=GetFileByte(trim(FileStrWrite2))
    if(FileByte /= int(inx,8)*int(inz,8)*nTime*int(RKF,8)) then
      print*,'FileByte Error-2: ',trim(FileStrWrite2); stop
    endif
    print*,'Extract bed from ',trim(FileStrRead),' to ',trim(FileStrWrite2),' successfully!'
  enddo
  call DeleteFileNameList()
end program main
  
#undef nDumpPrtclInte
#undef nDumpPrtclReal

#undef iDump_itime
#undef iDump_id
#undef iDump_IsCntct
