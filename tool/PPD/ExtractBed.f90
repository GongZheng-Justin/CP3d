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
! gfortran -O3 -cpp -Wall -ffree-line-length-none -fbacktrace -g ExtractBed.f90 -o generateBed
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
  integer,parameter::ifirst=0
  integer,parameter::ilast= 92000
  integer,parameter::DumpFreq=20    
  real(8),parameter::Gravity = 9.81_8
  real(8),parameter::Diameter=0.001_8
  real(8),parameter::rho_p  =1700.0_8
  real(8),parameter::rho_f  =1000.0_8
  real(8),parameter::Shield_c= 0.05_8
  real(8),parameter::NetWeight_c=1.0E-5_8  
  integer,parameter::numPrtcl=1000000
  character(len=128),parameter::ReadFromDir="./PrtclInfo/"
  character(len=128),parameter::WriteToDir ="./PrtclOut/"
  
  ! locals
  logical::BedFlag
  real(8)::u_th,F_th,y_th,u_mag,realTmp
  character(128)::FileStrRead,FileStrWrite,StrTmp,XdmfFile
  integer(8)::FileByte,disp,OneDumpByte,nOutTotal,disp_xdmf
  real(RKF),dimension(:,:),allocatable::GPrtcl_Pos,GPrtcl_linVel,GPrtcl_CntctForce
  integer::i,nUnitRead,nUnitWrite,nUnitXdmf,ierror,nTime,pid,itime,npOut,indent,nflds,iFile,dims,iprec
  
  y_th=8.0*Diameter
  u_th=sqrt(Shield_c)*sqrt((rho_p/rho_f-1.0_8)*Gravity*Diameter)
  F_th=NetWeight_c*(rho_p-rho_f)*PI*(Diameter**3)*Gravity/6.0_8;
  
  ! Xdmf
  indent=2
  nflds = (ilast - ifirst)/DumpFreq
  write(XdmfFile,"(A)") "PrtclBed.xmf"
  open(newunit=nUnitXdmf, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierror)
  write(nUnitXdmf,'(A)') '<?xml version="1.0" ?>'
  write(nUnitXdmf,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(nUnitXdmf,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  write(nUnitXdmf,'(A)') '<Domain>'  
  write(nUnitXdmf,'(A)')repeat(' ',indent)//'<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  indent = indent + 2
  write(nUnitXdmf,'(A)')repeat(' ',indent)//'<Time TimeType="List">'
  indent = indent + 2
  write(nUnitXdmf,'(A,I6,A)')repeat(' ',indent)//'<DataItem Format="XML" NumberType="Int" Dimensions="',nflds,'">' 
  write(nUnitXdmf,'(A)',advance='no') repeat(' ',indent)
  do i = 1,nflds
    if(mod(i,20)/=0) then
      write(nUnitXdmf,'(I9)',advance='no')  i*DumpFreq+ifirst
    else
      write(nUnitXdmf,'(I9)')  i*DumpFreq+ifirst
    endif
  enddo
  write(nUnitXdmf,'(A)')repeat(' ',indent)//'</DataItem>'
  indent = indent- 2
  write(nUnitXdmf,'(A)')repeat(' ',indent)//'</Time>'  
  
  allocate(GPrtcl_Pos(3,numPrtcl))
  allocate(GPrtcl_linVel(3,numPrtcl))
  allocate(GPrtcl_CntctForce(3,numPrtcl))
  OneDumpByte=int(IKF*nDumpPrtclInte+RKF*nDumpPrtclReal,8)*int(numPrtcl,8)
  
  ! Read file name
  call WriteFileNameList(ReadFromDir,"PrtclArrange_*")
  call FileSystem%GetFileName(ReadFromDir)  
  iFile=ifirst
  do i=1,FileSystem%nFile
    write(FileStrRead,'(A,A)')trim(ReadFromDir),trim(FileSystem%FileName(i))
    FileByte=GetFileByte(trim(FileStrRead))
    if(mod(FileByte,OneDumpByte)/=0) then
      print*,'FileByte Error-1: ',trim(FileStrRead); stop
    endif
    nTime=int(FileByte/OneDumpByte,4)
        
    nOutTotal=0; disp_xdmf=0_8
    StrTmp=FileSystem%FileName(i)
    itime=len_trim(StrTmp)
    write(FileStrWrite,'(A,A,A)')trim(WriteToDir),'PrtclBed_',StrTmp(itime-9:itime)
    open(newunit=nUnitRead, file=trim(FileStrRead), status='old',    form='unformatted',access='stream',action='read', IOSTAT=ierror)
    open(newunit=nUnitWrite,file=trim(FileStrWrite),status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
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
      
      npOut=numPrtcl
      do pid=numPrtcl,1,-1
        BedFlag=.true.
        u_mag=sqrt(GPrtcl_linVel(1,pid)*GPrtcl_linVel(1,pid)+GPrtcl_linVel(2,pid)*GPrtcl_linVel(2,pid) &
                  +GPrtcl_linVel(3,pid)*GPrtcl_linVel(3,pid)) 
        if(u_mag>u_th) BedFlag=.false.
        if(abs(GPrtcl_CntctForce(2,pid))<=F_th .and. GPrtcl_Pos(2,pid)>=y_th) BedFlag=.false.
        if(BedFlag) cycle
        if(pid/=npOut) then
          GPrtcl_Pos(:,pid)       =GPrtcl_Pos(:,npOut)
          GPrtcl_linVel(:,pid)    =GPrtcl_linVel(:,npOut)
          GPrtcl_CntctForce(:,pid)=GPrtcl_CntctForce(:,npOut)
        endif
        npOut=npOut-1
      enddo
      nOutTotal=nOutTotal+int(npOut,8)
      write(unit=nUnitWrite,IOSTAT=ierror) npOut
      write(unit=nUnitWrite,IOSTAT=ierror) GPrtcl_Pos(1:3,1:npOut)
      do pid=1,npOut
        realTmp=real(GPrtcl_Pos(2,pid),8)
        GPrtcl_Pos(2,pid)=real(realTmp/Diameter,4)
      enddo
      write(unit=nUnitWrite,IOSTAT=ierror) GPrtcl_Pos(2,  1:npOut)

      ! xdmf
      indent=4; iFile=iFile+DumpFreq; disp_xdmf=disp_xdmf+int(IKF,8) ! Skip "npOut"
      dims=3; iprec=RKF
      write(nUnitXdmf,'(A,I10.10,A)')repeat(' ',indent)//'<Grid Name="T',iFile,'" GridType="Uniform">'
      indent = indent + 2
      write(nUnitXdmf,'(A,I9,A)') repeat(' ',indent)//'<Topology TopologyType="Polyvertex" NodesPerElement="',npOut,'"/>'
      write(nUnitXdmf,'(A)')repeat(' ',indent)//'<Geometry GeometryType="'//"XYZ"//'">'
      indent = indent + 2
      write(nUnitXdmf,'(A,I1,A,I2,I9,A,I15,A)')repeat(' ',indent)// '<DataItem Format="Binary"' // &
            ' DataType="Float" Precision="',RKF,'" Endian="Native"' // &
            ' Dimensions="',dims,npOut,'" Seek="',disp_xdmf,'">'
      disp_xdmf = disp_xdmf + int(npOut,8)*int(dims*iprec,8)
      indent = indent + 2
      write(nUnitXdmf,'(A,I10.10)')repeat(' ',indent)//trim(FileStrWrite)
      indent = indent - 2
      write(nUnitXdmf,'(A)')repeat(' ',indent)//'</DataItem>'
      indent = indent - 2
      write(nUnitXdmf,'(A)')repeat(' ',indent)//'</Geometry>'
      
      dims=1; iprec=RKF
      call Write_XDMF_One(nUnitXdmf,dims,iprec,npOut,FileStrWrite,"PosY/D","Scalar","Float",disp_xdmf)
      write(nUnitXdmf,'(A)')'    </Grid>'
    enddo
    close(unit=nUnitRead, IOSTAT=ierror)
    close(unit=nUnitWrite,IOSTAT=ierror)
    FileByte=GetFileByte(trim(FileStrWrite))
    if(FileByte /= int(RKF*4,8)*nOutTotal+int(IKF,8)*nTime) then
      print*,'FileByte Error-2: ',trim(FileStrWrite); stop
    endif
    print*,'Extract bed from ',trim(FileStrRead),' to ',trim(FileStrWrite),' successfully!'
  enddo
  call DeleteFileNameList()

    ! XDMF/XMF Tail
  write(nUnitXdmf,'(A)') '  </Grid>'
  write(nUnitXdmf,'(A)') '</Domain>'
  write(nUnitXdmf,'(A)') '</Xdmf>'
  close(unit=nUnitXdmf, IOSTAT=ierror)
end program main

  !**********************************************************************
  ! Write_XDMF_One
  !**********************************************************************
  subroutine Write_XDMF_One(nUnitFile,dims,iprec,np,chFile,chName,chAttribute,chDataType,disp)
    implicit none
    integer,intent(in)::nUnitFile,dims,iprec,np
    character(*),intent(in)::chFile,chName,chAttribute,chDataType
    integer(kind=8),intent(inout)::disp
    
    ! locals
    integer:: indent
    
    indent = 6
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Attribute Type="'//trim(chAttribute)//'" Center="Node" Name="'//trim(chName)//'">'
    indent = indent + 2
    write(nUnitFile,'(3A,I1,A,I2,I9,A,I15,A)')repeat(' ',indent)// '<DataItem Format="Binary"' // &
          ' DataType="',trim(chDataType),'" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,np,'" Seek="',disp,'">'
    disp = disp+np*dims*iprec
    indent = indent + 2
    write(nUnitFile,'(A)')repeat(' ',indent)//trim(chFile)
    indent = indent - 2
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 2
    write(nUnitFile,'(A)')repeat(' ',indent)//'</Attribute>'
  end subroutine Write_XDMF_One
  
#undef nDumpPrtclInte
#undef nDumpPrtclReal

#undef iDump_itime
#undef iDump_id
#undef iDump_IsCntct
