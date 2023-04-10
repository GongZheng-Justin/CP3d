module ATP_IOAndVisu
  use MPI
  use ATP_Comm
  use m_TypeDef
  use m_LogInfo
  use ATP_Property
  use ATP_Variables
  use ATP_Decomp_2d
  use ATP_Parameters
  use m_Decomp2d,only: nrank,nproc
  implicit none
  private

  integer,parameter::IK=4
  integer::Prev_BackUp_itime= 53456791
  logical::saveXDMFOnce,save_ID,save_Type,save_UsrMark,save_LinVel,save_SwimDir,save_MoveDist

  type::part_io_size_vec
    integer,dimension(1)::sizes
    integer,dimension(1)::subsizes
    integer,dimension(1)::starts
  end type part_io_size_vec
  type::part_io_size_mat
    integer,dimension(2)::sizes
    integer,dimension(2)::subsizes
    integer,dimension(2)::starts
  end type part_io_size_mat

  type:: Prtcl_IO_Visu
  contains
    procedure:: Init_visu     =>  PIO_Init_visu
    procedure:: Final_visu    =>  PIO_Final_visu
    procedure:: Dump_visu     =>  PIO_Dump_visu
    procedure:: Read_Restart  =>  PIO_Read_Restart
    procedure:: Write_Restart =>  PIO_Write_Restart
    procedure:: Delete_Prev_Restart =>  PIO_Delete_Prev_Restart
    procedure,private:: Write_XDMF  =>  PIO_Write_XDMF
  end type Prtcl_IO_Visu
  type(Prtcl_IO_Visu),public:: ATP_IO

  ! useful interfaces
  interface Prtcl_dump
    module procedure Prtcl_dump_int_vector,  Prtcl_dump_int_matrix
    module procedure Prtcl_dump_real_vector, Prtcl_dump_real3_vector
  end interface Prtcl_dump

contains

  !**********************************************************************
  ! PIO_Init_visu
  !**********************************************************************
  subroutine PIO_Init_visu(this,chFile,iStage)
    implicit none
    class(Prtcl_IO_Visu)::this
    character(*),intent(in)::chFile
    integer,intent(in)::iStage
    
    ! locals
    integer::nUnitFile,ierror,indent,nflds,ifld
    NAMELIST /PrtclVisuOption/saveXDMFOnce,save_ID,save_Type,save_UsrMark,save_LinVel,save_SwimDir,save_MoveDist
    character(128)::XdmfFile
  
    if(iStage==1) then
      open(newunit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=ierror)
      if(ierror/=0)call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Init_visu", "Cannot open file: "//trim(chFile))
      read(nUnitFile, nml=PrtclVisuOption)
      if(nrank==0)write(ATPLogInfo%nUnit, nml=PrtclVisuOption)
      close(nUnitFile,IOSTAT=ierror)
      return
    endif

    ! initialize the XDMF/XDF file
    if(nrank/=0) return
    write(xdmfFile,"(A)") trim(ATP_opt%ResultsDir)//"PartVisuFor"//trim(ATP_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror /= 0) call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Init_visu","Cannot open file: "//trim(XdmfFile))
    ! XDMF/XMF Title
    write(nUnitFile,'(A)') '<?xml version="1.0" ?>'
    write(nUnitFile,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(nUnitFile,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(nUnitFile,'(A)') '<Domain>'

    ! Time series
    indent =  4
    nflds = (ATP_Opt%ilast - ATP_Opt%ifirst +1)/ATP_Opt%SaveVisu  + 1
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    indent = indent + 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Time TimeType="List">'
    indent = indent + 4
    write(nUnitFile,'(A,I6,A)')repeat(' ',indent)//'<DataItem Format="XML" NumberType="Int" Dimensions="',nflds,'">' 
    write(nUnitFile,'(A)',advance='no') repeat(' ',indent)
    do ifld = 1,nflds
      write(nUnitFile,'(I9)',advance='no') ((ifld-1)*ATP_Opt%SaveVisu + ATP_Opt%ifirst-1)
    enddo
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,fmt='(A)')repeat(' ',indent)//'</Time>'
    close(nUnitFile,IOSTAT=ierror)
    if( .not. saveXDMFOnce) return
    
    do ifld = 1,nflds
      call this%Write_XDMF((ifld-1)*ATP_Opt%SaveVisu + ATP_Opt%ifirst-1)
    enddo

    ! XDMF/XMF Tail
    open(newunit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    write(nUnitFile,'(A)') '    </Grid>'
    write(nUnitFile,'(A)') '</Domain>'
    write(nUnitFile,'(A)') '</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)    
  end subroutine PIO_Init_visu

  !**********************************************************************
  ! Purpose:
  !   Create a xdmf/xmf file in order to view the simulation results
  !     by Paraview directly
  ! 
  ! Original Author: 
  !   Pedro Costa
  ! 
  ! Modified by:
  !   Zheng Gong
  ! 
  ! Original Source file is downloaded from ( April 2020 ):
  !   https://github.com/p-costa/gen_xdmf_particles
  !              
  !**********************************************************************
  subroutine PIO_Write_XDMF(this,itime) 
    implicit none 
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime

    !locals
    integer(kind=MPI_OFFSET_KIND)::disp
    integer:: indent,nUnitFile,ierror,np,dims,iprec
    character(128)::XdmfFile

    if(nrank/=0) return 
    np=ATP_Opt%np_InDomain
    write(xdmfFile,"(A)") trim(ATP_opt%ResultsDir)//"PartVisuFor"//trim(ATP_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Write_XDMF","Cannot open file: "//trim(XdmfFile))

    indent = 8; disp = 0_MPI_OFFSET_KIND
    write(xdmfFile,"(A)") "PartVisuFor"//trim(ATP_opt%RunName)
    dims=3; iprec=RK
    write(nUnitFile,'(A,I10.10,A)')repeat(' ',indent)//'<Grid Name="T',itime,'" GridType="Uniform">'
    indent = indent + 4
    write(nUnitFile,'(A,I9,A)')repeat(' ',indent)//'<Topology TopologyType="Polyvertex" NodesPerElement="',np,'"/>'
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Geometry GeometryType="'//"XYZ"//'">'
    indent = indent + 4
    write(nUnitFile,'(A,I1,A,I2,I9,A,I15,A)')repeat(' ',indent)// '<DataItem Format="Binary"' // &
          ' DataType="Float" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,np,'" Seek="',disp,'">'
    disp = disp+np*dims*iprec
    indent = indent + 4
    write(nUnitFile,'(A,I10.10)')repeat(' ',indent)//trim(XdmfFile),itime
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</Geometry>'

    IF(save_ID) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"ID","Scalar","Int",disp)
    ENDIF
    IF(save_Type) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"Type","Scalar","Int",disp)
    ENDIF
    IF(save_UsrMark) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"UsrMark","Scalar","Int",disp)
    ENDIF
    IF(save_LinVel) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"LinVel","Vector","Float",disp)
    ENDIF
    IF(save_SwimDir) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"SwimDir","Vector","Float",disp)
    ENDIF
    IF(save_MoveDist) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,XdmfFile,"MoveDist","Vector","Float",disp)
    ENDIF
    write(nUnitFile,'(A)')'        </Grid>'
    close(nUnitFile)
  end subroutine PIO_Write_XDMF

  !**********************************************************************
  ! Write_XDMF_One
  !**********************************************************************
  subroutine Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,chName,chAttribute,chDataType,disp)
    implicit none
    integer,intent(in)::nUnitFile,dims,iprec,np,itime
    character(*),intent(in)::chFile,chName,chAttribute,chDataType
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    
    ! locals
    integer:: indent
    indent = 12

    write(nUnitFile,'(A)')repeat(' ',indent)//'<Attribute Type="'//trim(chAttribute)//'" Center="Node" Name="'//trim(chName)//'">'
    indent = indent + 4
    write(nUnitFile,'(3A,I1,A,I2,I9,A,I15,A)')repeat(' ',indent)// '<DataItem Format="Binary"' // &
          ' DataType="',trim(chDataType),'" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,np,'" Seek="',disp,'">'
    disp = disp+np*dims*iprec
    indent = indent + 4
    write(nUnitFile,'(A,I10.10)')repeat(' ',indent)//trim(chFile),itime
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</Attribute>'
  end subroutine Write_XDMF_One

  !**********************************************************************
  ! PIO_Delete_Prev_Restart
  !**********************************************************************
  subroutine PIO_Delete_Prev_Restart(this,itime)
    implicit none
    class(Prtcl_IO_Visu)::this
    integer:: itime

    ! locals
    integer::nUnit,ierror
    character(128)::chFile

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if(nrank/=0) return
    !
    write(chFile,"(A,I10.10)") trim(ATP_opt%RestartDir)//"RestartFor"//trim(ATP_opt%RunName),Prev_BackUp_itime 
    open(newunit=nUnit,file=trim(chFile),IOSTAT=ierror)
    close(unit=nUnit,status='delete',IOSTAT=ierror)
    !    
    Prev_BackUp_itime = itime
  end subroutine PIO_Delete_Prev_Restart

  !**********************************************************************
  ! PIO_Read_Restart
  !**********************************************************************
  subroutine PIO_Read_Restart(this)
    implicit none
    class(Prtcl_IO_Visu)::this

    ! locals
    character(128)::chFile
    integer,parameter::NumRead=2000
    real(RK)::xst,xed,yst,yed,zst,zed
    integer,allocatable,dimension(:):: nP_in_bin
    type(real3),allocatable,dimension(:)::real3Vec,PosVec
    integer(kind=8)::disp,disp_pos,disp_int,disp_real3
    integer::itime,nUnit,ierror,nlocal,np,i,k,itype,nlocal_sum,nreal3,tsize,nLeft,nRead,int_t(3)

    itime = ATP_Opt%ifirst - 1
    xst=ATP_decomp%xSt; xed=ATP_decomp%xEd
    yst=ATP_decomp%ySt; yed=ATP_decomp%yEd
    zst=ATP_decomp%zSt; zed=ATP_decomp%zEd

    ! Begin to write Restart file
    write(chFile,"(A,I10.10)") trim(ATP_opt%RestartDir)//"RestartFor"//trim(ATP_opt%RunName),itime
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart","Cannot open file: "//trim(chFile))
    disp =1_8; read(nUnit,pos=disp,IOSTAT=ierror)np; disp=disp+int_byte
    if(np>ATP_Opt%numPrtcl .and. nrank==0) then
      call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart"," np_InDomain > numPrtcl " )
    endif
    ATP_Opt%np_InDomain = np

    tsize=GPrtcl_list%tsize
    nreal3 = 2*tsize+2
    allocate(real3Vec(nreal3),PosVec(NumRead))
    allocate(nP_in_bin(ATP_opt%numPrtcl_Type)); nP_in_bin=0

    nlocal=0; nLeft=np
    disp_pos  = disp
    disp_int  = disp_pos+ real3_byte*np
    disp_real3= disp_int+ int_byte*np*3
    DO
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp_pos,IOSTAT=ierror)PosVec(1:nRead)
      disp_pos=disp_pos+int(real3_byte,8)*int(nRead,8)
      do i=1,nRead
        if(PosVec(i)%x>=xst .and. PosVec(i)%x< xed .and. PosVec(i)%y>=yst .and. &
           PosVec(i)%y< yed .and. PosVec(i)%z>=zst .and. PosVec(i)%z<zed) then
          if(nlocal>=GPrtcl_list%mlocal)  call GPrtcl_list%ReallocatePrtclVar(nlocal)
          nlocal=nlocal+1

          read(nUnit,pos=disp_int,IOSTAT=ierror)int_t(1:3)
          GPrtcl_id(nlocal)=int_t(1)      ! id
          itype=int_t(2)
          GPrtcl_pType(nlocal)=itype      ! pType
          nP_in_bin(itype)= nP_in_bin(itype)+1
          GPrtcl_UsrMark(nlocal)=int_t(3) ! Usr_Mark

          GPrtcl_PosR(nlocal)= PosVec(i)  ! PosR
          k=0;
          read(nUnit,pos=disp_real3,IOSTAT=ierror)real3Vec(1:nreal3)
          GPrtcl_SwimDir(nlocal)=real3Vec(k+1);       k=k+1
          GPrtcl_MoveDistance(nlocal)=real3Vec(k+1);  k=k+1
          GPrtcl_LinVel(1:tsize,nlocal)  = real3Vec(k+1:k+tsize); k=k+tsize ! LinVec
          GPrtcl_SwimAcc(1:tsize,nlocal) = real3Vec(k+1:k+tsize); k=k+tsize ! LinAcc
        endif
        disp_int  = disp_int  + int_byte*3
        disp_real3= disp_real3+ real3_byte*nreal3       
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    deallocate(PosVec,Real3Vec)
    call MPI_ALLREDUCE(nP_in_bin, ATPProperty%nPrtcl_in_Bin,ATP_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    deallocate(nP_in_bin)
    close(nUnit,IOSTAT=ierror)
    GPrtcl_list%nlocal = nlocal
    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nlocal_sum/= np .and. nrank==0) then
      call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart: "," nlocal_sum/= np_InDomain " )
    endif
  end subroutine PIO_Read_Restart

  !**********************************************************************
  ! PIO_Write_Restart
  !**********************************************************************
  subroutine PIO_Write_Restart(this,itime)
    implicit none 
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime
   
    ! locals 
    character(128)::chFile
    type(part_io_size_vec)::pvsize
    type(part_io_size_mat)::pmsize
    integer,parameter::NumRestart=200
    integer,allocatable,dimension(:,:)::IntMat
    type(real3),allocatable,dimension(:)::real3Vec
    integer(kind=MPI_OFFSET_KIND)::disp,bgn_byte,FileSize
    integer::pid,i,k,nlocal,bgn_ind,prank,pProc,ierror,fh
    integer::tsize,color,key,Prtcl_WORLD,nreal3,nRestart,nLeft

    ! Calculate the bgn_ind
    nlocal = GPrtcl_list%nlocal
    bgn_ind= clc_bgn_ind(nlocal)

    ! Create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD,ierror)
    if(color==2) return

    ! Begin to write Restart file
    write(chFile,"(A,I10.10)") trim(ATP_opt%RestartDir)//"RestartFor"//trim(ATP_opt%RunName),itime
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_8,ierror)  ! Guarantee overwriting
    call MPI_BARRIER(Prtcl_WORLD,ierror)
    disp = 0_MPI_OFFSET_KIND

    ! Write ATP_Opt%np_InDomain, in the begining of the Restart file
    call MPI_COMM_RANK(Prtcl_WORLD, prank, ierror)
    call MPI_COMM_SIZE(Prtcl_WORLD, pProc, ierror)
    if(prank==0) then
      call MPI_FILE_WRITE_AT(fh,disp,ATP_Opt%np_InDomain,1,int_type,MPI_STATUS_IGNORE,ierror)
    endif
    disp = disp + int_byte
    call MPI_BARRIER(Prtcl_WORLD,ierror)

    ! Begin to write
    pvsize%sizes(1)   = ATP_Opt%np_InDomain
    pvsize%subsizes(1)= nlocal
    pvsize%starts(1)  = bgn_ind
    allocate(real3Vec(nlocal))
    do pid=1,nlocal
      real3Vec(pid)=GPrtcl_PosR(pid)
    enddo
    call Prtcl_dump(fh,disp, real3Vec(1:nlocal),  pvsize)
    deallocate(real3Vec)
    
    pmsize%sizes(1)   = 3;   pmsize%sizes(2)   = ATP_Opt%np_InDomain
    pmsize%subsizes(1)= 3;   pmsize%subsizes(2)= nlocal
    pmsize%starts(1)  = 0;   pmsize%starts(2)  = bgn_ind
    allocate(IntMat(3,nlocal))
    do pid=1,nlocal
      IntMat(1,pid)= GPrtcl_id(pid)
      IntMat(2,pid)= GPrtcl_pType(pid)
      IntMat(3,pid)= GPrtcl_UsrMark(pid)
    enddo
    call Prtcl_dump(fh,disp,IntMat(1:3,1:nlocal),  pmsize)
    deallocate(IntMat)
    call MPI_FILE_CLOSE(fh,ierror)
    
    tsize=GPrtcl_list%tsize
    nreal3 = 2*tsize+2
    call MPI_FILE_OPEN(Prtcl_WORLD,chFile, MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierror)
    call MPI_FILE_GET_SIZE(fh,FileSize,ierror)
    FileSize=FileSize+int(nreal3*real3_byte,8)*int(ATP_Opt%np_InDomain,8)
    call MPI_BARRIER(Prtcl_WORLD,ierror)
    call MPI_FILE_PREALLOCATE(fh,FileSize,ierror)
    call MPI_BARRIER(Prtcl_WORLD,ierror)
    
    allocate(real3Vec(NumRestart*nreal3))
    nLeft=nlocal; pid=0
    bgn_byte=disp+int(nreal3*real3_byte,8)*int(bgn_ind,8)
    DO
      nRestart=min(nLeft,NumRestart)
      k=0
      do i=1,nRestart
        pid=pid+1
        real3Vec(k+1)=GPrtcl_SwimDir(pid);      k=k+1
        real3Vec(k+1)=GPrtcl_MoveDistance(pid); k=k+1
        real3Vec(k+1:k+tsize)=GPrtcl_LinVel(1:tsize,pid); k=k+tsize
        real3Vec(k+1:k+tsize)=GPrtcl_SwimAcc(1:tsize,pid); k=k+tsize
      enddo
      call MPI_FILE_WRITE_AT(fh,bgn_byte,real3Vec,k,real3_type,MPI_STATUS_IGNORE,ierror)
      bgn_byte=bgn_byte+int(nreal3*real3_byte,8)*int(nRestart,8)
      nLeft=nLeft-nRestart
      if(nLeft==0)exit
    ENDDO
    deallocate(real3Vec)
    call MPI_FILE_CLOSE(fh, ierror)
    call MPI_COMM_FREE(Prtcl_WORLD,ierror)
  end subroutine PIO_Write_Restart

  !**********************************************************************
  ! PIO_Init_visu
  !**********************************************************************
  subroutine PIO_Final_visu(this)
    implicit none 
    class(Prtcl_IO_Visu)::this

    ! locals
    integer::nUnitFile,ierror
    character(128)::XdmfFile

    if(nrank/=0 .or. saveXDMFOnce) return
    write(xdmfFile,"(A)") trim(ATP_opt%ResultsDir)//"PartVisuFor"//trim(ATP_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Final_visu","Cannot open file:  "//trim(XdmfFile))
    ! XDMF/XMF Tail
    write(nUnitFile,'(A)') '    </Grid>'
    write(nUnitFile,'(A)') '</Domain>'
    write(nUnitFile,'(A)') '</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)

  end subroutine PIO_Final_visu

  !**********************************************************************
  ! PIO_Dump_visu
  !**********************************************************************
  subroutine PIO_Dump_visu(this, itime)
    implicit none
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime

    ! locals
    character(128)::chFile
    type(part_io_size_vec)::pvsize
    integer(kind=MPI_OFFSET_KIND)::disp
    integer,allocatable,dimension(:)::intVec
    real(RK),allocatable,dimension(:)::realVec
    type(real3),allocatable,dimension(:)::real3Vec
    integer :: ierror,fh,i,color,key,Prtcl_WORLD,nlocal,bgn_ind

    ! write xdmf file first
    if(.not.saveXDMFOnce) call this%Write_XDMF(itime)

    ! update the bgn_ind
    nlocal = GPrtcl_list%nlocal
    bgn_ind=clc_bgn_ind(nlocal)

    ! create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD,ierror)
    if(color==2) return
    
    ! begin to dump
    write(chFile,"(A,I10.10)") trim(ATP_opt%ResultsDir)//"PartVisuFor"//trim(ATP_opt%RunName),itime
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_8,ierror)  ! guarantee overwriting
    call MPI_BARRIER(Prtcl_WORLD,ierror)
    disp = 0_MPI_OFFSET_KIND
    pvsize%sizes(1)     = ATP_Opt%np_InDomain
    pvsize%subsizes(1)  = nlocal
    pvsize%starts(1)    = bgn_ind
    if(nlocal<=0) return

    allocate(real3Vec(nlocal))
    do i=1,nlocal
      real3Vec(i)=GPrtcl_PosR(i)
    enddo
    call Prtcl_dump(fh,disp, real3Vec(1:nlocal),  pvsize)
    deallocate(real3Vec)
    if(save_ID)       call Prtcl_dump(fh,disp, GPrtcl_id(1:nlocal),       pvsize)
    if(save_Type)     call Prtcl_dump(fh,disp, GPrtcl_pType(1:nlocal),    pvsize)
    if(save_UsrMark)  call Prtcl_dump(fh,disp, GPrtcl_UsrMark(1:nlocal),  pvsize)
    if(save_LinVel)   call Prtcl_dump(fh,disp, GPrtcl_LinVel(1,1:nlocal), pvsize)
    if(save_SwimDir)  call Prtcl_dump(fh,disp, GPrtcl_SwimDir(1:nlocal),  pvsize)
    if(save_MoveDist) call Prtcl_dump(fh,disp, GPrtcl_MoveDistance(1:nlocal),  pvsize)
    call MPI_FILE_CLOSE(fh, ierror) 
    call MPI_COMM_FREE(Prtcl_WORLD, ierror)
  end subroutine PIO_Dump_visu

  !**********************************************************************
  ! clc_bgn_ind
  !**********************************************************************
  function clc_bgn_ind(nlocal) result(bgn_ind)
    implicit none
    integer,intent(in)::nlocal
    integer::bgn_ind

    ! locals
    integer::end_ind,ierror,SRstatus(MPI_STATUS_SIZE)
  
    bgn_ind=0
    if(nproc<=1) return
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    IF(nrank==0)THEN
      end_ind=bgn_ind+nlocal
      call MPI_SEND(end_ind, 1,int_type, nrank+1,0,MPI_COMM_WORLD,ierror)
    ELSEIF(nrank /= nproc-1) THEN
      call MPI_RECV(bgn_ind, 1,int_type, nrank-1,0,MPI_COMM_WORLD,SRstatus,ierror)
      end_ind=bgn_ind+nlocal
      call MPI_SEND(end_ind, 1,int_type, nrank+1,0,MPI_COMM_WORLD,ierror)
    ELSE
      call MPI_RECV(bgn_ind, 1,int_type, nrank-1,0,MPI_COMM_WORLD,SRstatus,ierror)
    ENDIF
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end function clc_bgn_ind

  !**********************************************************************
  ! Prtcl_dump_int_vector
  !**********************************************************************
  subroutine Prtcl_dump_int_vector(fh,disp,var,pvsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_vec),intent(in)::pvsize  
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    integer,dimension(1:pvsize%subsizes(1)),intent(in)::var

    ! locals
    integer:: ierror,newtype
    integer,dimension(1) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes    = pvsize%sizes
    subsizes = pvsize%subsizes
    starts   = pvsize%starts

    ! write the particle revelant integer vector
    call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, int_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,int_type, newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),int_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp = disp + sizes(1) * int_byte
  end subroutine Prtcl_dump_int_vector

  !**********************************************************************
  ! Prtcl_dump_int_matrix
  !**********************************************************************
  subroutine Prtcl_dump_int_matrix(fh,disp,var,pmsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_mat),intent(in)::pmsize    
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    integer,dimension(1:pmsize%subsizes(1),1:pmsize%subsizes(2)),intent(in)::var

    ! locals
    integer :: ierror,newtype
    integer, dimension(2) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes     = pmsize%sizes
    subsizes  = pmsize%subsizes
    starts    = pmsize%starts

    ! write the particle relevant real matrix
    call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, int_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,int_type, newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1)*subsizes(2),int_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp = disp + sizes(1) * sizes(2) * int_byte
  end subroutine Prtcl_dump_int_matrix

  !**********************************************************************
  ! Prtcl_dump_real_vector
  !**********************************************************************
  subroutine Prtcl_dump_real_vector(fh,disp,var,pvsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_vec),intent(in)::pvsize  
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    real(RK),dimension(1:pvsize%subsizes(1)),intent(in)::var

    ! locals
    integer:: ierror,newtype
    integer,dimension(1) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes    = pvsize%sizes
    subsizes = pvsize%subsizes
    starts   = pvsize%starts

    ! write the particle revelant real vector
    call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,real_type, newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp = disp + sizes(1) * real_byte
  end subroutine Prtcl_dump_real_vector

  !**********************************************************************
  ! Prtcl_dump_real3_vector
  !**********************************************************************
  subroutine Prtcl_dump_real3_vector(fh,disp,var,pvsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_vec),intent(in)::pvsize
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    type(real3),dimension(1:pvsize%subsizes(1)),intent(in)::var

    ! locals
    integer:: ierror,newtype
    integer,dimension(1) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes    = pvsize%sizes
    subsizes = pvsize%subsizes
    starts   = pvsize%starts

    ! write the particle revelant real3 vector
    call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real3_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,real3_type, newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),real3_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp = disp + sizes(1) * real3_byte
  end subroutine Prtcl_dump_real3_vector

end module ATP_IOAndVisu
