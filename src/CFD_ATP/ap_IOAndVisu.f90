#include "definitions_inc.f90"
module ap_IOAndVisu
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d,only: nrank,nproc
  use ap_Comm
  use ap_Property
  use ap_Variables
  use ap_Decomp_2d
  use ap_Parameters
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
    procedure:: Init_Visu     => PIO_Init_Visu
    procedure:: Dump_Visu     => PIO_Dump_Visu
    procedure,nopass:: Final_Visu    => PIO_Final_Visu
    procedure,nopass:: Read_Restart  => PIO_Read_Restart
    procedure,nopass:: Write_Restart => PIO_Write_Restart
    procedure,nopass:: Delete_Prev_Restart =>  PIO_Delete_Prev_Restart
    procedure,nopass,private:: Write_XDMF  =>  PIO_Write_XDMF
  end type Prtcl_IO_Visu
  type(Prtcl_IO_Visu),public:: ATP_IO

  ! useful interfaces
  interface Prtcl_dump
    module procedure Prtcl_dump_int_vector,  Prtcl_dump_int_matrix
    module procedure Prtcl_dump_real_vector, Prtcl_dump_real3_vector
  end interface Prtcl_dump
  public:: Prtcl_dump
contains
#include "Prtcl_Dump_MPI_inc.f90"

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PIO_Init_Visu
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PIO_Init_Visu(this,chFile,iStage); implicit none
    class(Prtcl_IO_Visu)::this
    character(len=*),intent(in)::chFile
    integer,intent(in)::iStage
    
    ! locals
    integer::iUnit,ierr,nfld,ifld
    character(len=:),allocatable::XdmfFile
    NAMELIST /PrtclVisuOption/ saveXDMFOnce,save_ID,save_Type,save_UsrMark,save_LinVel,save_SwimDir,save_MoveDist
  
    if(iStage==1) then
      open(newunit=iUnit, file=chFile,status='old',form='formatted',IOSTAT=ierr)
      if(ierr/=0)call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Init_Visu", "Cannot open file: "//strip(chFile))
      read(iUnit, nml=PrtclVisuOption)
      if(nrank==0)write(ATPLogInfo%iUnit, nml=PrtclVisuOption)
      close(iUnit,IOSTAT=ierr)
      return
    endif

    ! Write XDMF file
    if(nrank/=0) return
    XdmfFile = strip(ATP_Opt%ResultsDir)//"PartVisuFor"//strip(ATP_Opt%RunName)//".xmf"
    open(newunit=iUnit, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierr)
    if(ierr /= 0) call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Init_Visu","Cannot open file:  "//XdmfFile)
    write(iUnit,'(A)') '<?xml version="1.0" ?>'
    write(iUnit,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(iUnit,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(iUnit,'(A)') '<Domain>'

    ! Time series
    nfld = (ATP_Opt%ilast - ATP_Opt%ifirst +1)/ATP_Opt%SaveVisu  + 1
    write(iUnit,'(A)')'  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    write(iUnit,'(A)')'    <Time TimeType="List">'
    write(iUnit,'(A,I6,A)')'      <DataItem Format="XML" NumberType="Int" Dimensions="',nfld,'">' 
    write(iUnit,'(A)',advance='no')'        '
    do ifld =1, nfld
      if(mod(ifld,10)==0) then
        write(iUnit,'(I10)') (ifld-1)*ATP_Opt%SaveVisu + ATP_Opt%ifirst-1
        if(ifld < nfld) write(iUnit,'(A)',advance='no') '        '
      else
        write(iUnit,'(I9)',advance='no') (ifld-1)*ATP_Opt%SaveVisu + ATP_Opt%ifirst-1
      endif
    enddo
    if(mod(nfld,10) /=0) write(iUnit,*)' '
    write(iUnit,'(A)') '      </DataItem>'
    write(iUnit,'(A)') '    </Time>'
    close(iUnit, IOSTAT=ierr)
    if(.not. saveXDMFOnce) return

    do ifld = 1,nfld
      call this%Write_XDMF( (ifld-1)*ATP_Opt%SaveVisu + ATP_Opt%ifirst-1 )
    enddo

    ! XDMF/XMF Tail
    open(newunit=iUnit, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierr)
    write(iUnit,'(A)')'  </Grid>'
    write(iUnit,'(A)')'</Domain>'
    write(iUnit,'(A)')'</Xdmf>'
    close(iUnit, IOSTAT=ierr)
  end subroutine PIO_Init_Visu

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
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
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PIO_Write_XDMF(itime); implicit none
    integer,intent(in)::itime

    ! locals
    integer:: iUnit,ierr,np,dims,iprec
    integer(kind=MPI_OFFSET_KIND)::disp
    character(len=:),allocatable::XdmfFile

    if(nrank/=0) return 
    np=ATP_Opt%np_InDomain
    XdmfFile = strip(ATP_Opt%ResultsDir)//"PartVisuFor"//strip(ATP_Opt%RunName)//".xmf"
    open(newunit=iUnit, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierr)
    if(ierr/=0) then
      call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Write_XDMF","Cannot open file: "//XdmfFile)
    endif

    disp = 0_MPI_OFFSET_KIND
    XdmfFile = "PartVisuFor"//strip(ATP_Opt%RunName)
    dims=3; iprec=RK
    write(iUnit,'(A,I10.10,A)') '    <Grid Name="T',itime,'" GridType="Uniform">'
    write(iUnit,'(A,I9,A)') '      <Topology TopologyType="Polyvertex" NodesPerElement="',np,'"/>'
    write(iUnit,'(A)') '      <Geometry GeometryType="'//"XYZ"//'">'
    write(iUnit,'(A,I1,A,I2,I9,A,I15,A)')  '        <DataItem Format="Binary"' // &
          ' DataType="Float" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,np,'" Seek="',disp,'">'
    disp = disp+np*dims*iprec
    write(iUnit,'(A,I10.10)') '          ' // XdmfFile, itime
    write(iUnit,'(A)') '        </DataItem>'
    write(iUnit,'(A)') '      </Geometry>'

    IF(save_ID) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"ID","Scalar","Int",disp)
    ENDIF
    IF(save_Type) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"Type","Scalar","Int",disp)
    ENDIF
    IF(save_UsrMark) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"UsrMark","Scalar","Int",disp)
    ENDIF
    IF(save_LinVel) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"LinVel","Vector","Float",disp)
    ENDIF
    IF(save_SwimDir) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"SwimDir","Vector","Float",disp)
    ENDIF
    IF(save_MoveDist) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"MoveDist","Vector","Float",disp)
    ENDIF
    write(iUnit,'(A)')'    </Grid>'
    close(iUnit)
  end subroutine PIO_Write_XDMF

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Write_XDMF_One
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,chName,chAttribute,chDataType,disp); implicit none
    integer,intent(in)::iUnit,dims,iprec,np,itime
    character(len=*),intent(in)::XdmfFile,chName,chAttribute,chDataType
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

    write(iUnit,'(A)') '      <Attribute Type="'//strip(chAttribute)//'" Center="Node" Name="'//strip(chName)//'">'
    write(iUnit,'(3A,I1,A,I2,I9,A,I15,A)')  '        <DataItem Format="Binary"' // &
          ' DataType="',chDataType,'" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,np,'" Seek="',disp,'">'
    disp = disp+np*dims*iprec
    write(iUnit,'(A,I10.10)') '          ' // XdmfFile, itime
    write(iUnit,'(A)') '        </DataItem>'
    write(iUnit,'(A)') '      </Attribute>'
  end subroutine Write_XDMF_One

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PIO_Delete_Prev_Restart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PIO_Delete_Prev_Restart(itime); implicit none
    integer:: itime

    ! locals
    integer::iUnit,ierr
    character(len=:),allocatable::chFile

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(nrank/=0) return
    !
    chFile = strip(ATP_Opt%RestartDir) // "RestartFor" // strip(ATP_Opt%RunName) //int2str(Prev_BackUp_itime,10)
    open(newunit=iUnit,file=chFile,IOSTAT=ierr)
    close(unit=iUnit,status='delete',IOSTAT=ierr)
    !    
    Prev_BackUp_itime = itime
  end subroutine PIO_Delete_Prev_Restart

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PIO_Read_Restart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
#define ATP_NumRead 2000
  subroutine PIO_Read_Restart(); implicit none

    ! locals   
    real(RK)::xst,xed,yst,yed,zst,zed
    character(len=:),allocatable::chFile
    integer,allocatable,dimension(:):: nP_in_bin
    type(real3),allocatable,dimension(:)::real3Vec,PosVec
    integer(kind=8)::disp,disp_pos,disp_int,disp_real3
    integer::itime,iUnit,ierr,nlocal,np,i,k,itype,nlocal_sum,nreal3,tsize,nLeft,nRead,int_t(3)

    itime = ATP_Opt%ifirst - 1
    xst=ATP_decomp%xSt; xed=ATP_decomp%xEd
    yst=ATP_decomp%ySt; yed=ATP_decomp%yEd
    zst=ATP_decomp%zSt; zed=ATP_decomp%zEd

    ! Begin to write Restart file
    chFile = strip(ATP_Opt%RestartDir)//"RestartFor"//strip(ATP_Opt%RunName) //int2str(itime,10)
    open(newunit=iUnit,file=chFile,status='old',form='unformatted',access='stream',action='read',IOSTAT=ierr)
    if(ierr/=0 .and. nrank==0) call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart","Cannot open file: "//chFile)
    disp =1_8; read(iUnit,pos=disp,IOSTAT=ierr)np; disp=disp+int_byte
    if(np>ATP_Opt%numPrtcl .and. nrank==0) then
      call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart"," np_InDomain > numPrtcl " )
    endif
    ATP_Opt%np_InDomain = np

    tsize=GPrtcl_list%tsize
    nreal3 = 2*tsize+2
    allocate(real3Vec(nreal3),PosVec(ATP_NumRead))
    allocate(nP_in_bin(ATP_Opt%numPrtcl_Type)); nP_in_bin=0

    nlocal=0; nLeft=np
    disp_pos  = disp
    disp_int  = disp_pos+ real3_byte*np
    disp_real3= disp_int+ int_byte*np*3
    DO
      nRead=min(nLeft,ATP_NumRead)
      read(iUnit,pos=disp_pos,IOSTAT=ierr)PosVec(1:nRead)
      disp_pos=disp_pos+int(real3_byte,8)*int(nRead,8)
      do i=1,nRead
        if(PosVec(i)%x>=xst .and. PosVec(i)%x< xed .and. PosVec(i)%y>=yst .and. &
           PosVec(i)%y< yed .and. PosVec(i)%z>=zst .and. PosVec(i)%z<zed) then
          if(nlocal>=GPrtcl_list%mlocal)  call GPrtcl_list%ReallocatePrtclVar(nlocal)
          nlocal=nlocal+1

          read(iUnit,pos=disp_int,IOSTAT=ierr)int_t(1:3)
          GPrtcl_id(nlocal)=int_t(1)      ! id
          itype=int_t(2)
          GPrtcl_pType(nlocal)=itype      ! pType
          nP_in_bin(itype)= nP_in_bin(itype)+1
          GPrtcl_UsrMark(nlocal)=int_t(3) ! Usr_Mark

          GPrtcl_PosR(nlocal)= PosVec(i)  ! PosR
          k=0;
          read(iUnit,pos=disp_real3,IOSTAT=ierr)real3Vec(1:nreal3)
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
    call MPI_ALLREDUCE(nP_in_bin, ATPProperty%nPrtcl_in_Bin,ATP_Opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    deallocate(nP_in_bin)
    close(iUnit,IOSTAT=ierr)
    GPrtcl_list%nlocal = nlocal
    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(nlocal_sum/= np .and. nrank==0) then
      call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart: "," nlocal_sum/= np_InDomain " )
    endif
  end subroutine PIO_Read_Restart
#undef ATP_NumRead

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PIO_Write_Restart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
#define ATP_NumRestart 200
  subroutine PIO_Write_Restart(itime); implicit none 
    integer,intent(in)::itime
   
    ! locals
    type(part_io_size_vec)::pvsize
    type(part_io_size_mat)::pmsize
    character(len=:),allocatable::chFile
    integer,allocatable,dimension(:,:)::IntMat
    type(real3),allocatable,dimension(:)::real3Vec
    integer(kind=MPI_OFFSET_KIND)::disp,bgn_byte,FileSize
    integer::pid,i,k,nlocal,bgn_ind,prank,pProc,ierr,fh
    integer::tsize,color,key,Prtcl_WORLD,nreal3,nRestart,nLeft

    ! Calculate the bgn_ind
    nlocal = GPrtcl_list%nlocal
    bgn_ind= clc_bgn_ind(nlocal)

    ! Create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD,ierr)
    if(color==2) return

    ! Begin to write Restart file
    chFile = strip(ATP_Opt%RestartDir) // "RestartFor" // strip(ATP_Opt%RunName) // int2str(itime,10)
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_BARRIER(Prtcl_WORLD,ierr)
    call my_mpi_file_set_size(fh,0_8,ierr)  ! Guarantee overwriting
    call MPI_BARRIER(Prtcl_WORLD,ierr)
    disp = 0_MPI_OFFSET_KIND

    ! Write ATP_Opt%np_InDomain, in the begining of the Restart file
    call MPI_COMM_RANK(Prtcl_WORLD, prank, ierr)
    call MPI_COMM_SIZE(Prtcl_WORLD, pProc, ierr)
    if(prank==0) then
      call MPI_FILE_WRITE_AT(fh,disp,ATP_Opt%np_InDomain,1,int_type,MPI_STATUS_IGNORE,ierr)
    endif
    disp = disp + int_byte
    call MPI_BARRIER(Prtcl_WORLD,ierr)

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
    call MPI_FILE_CLOSE(fh,ierr)
    
    tsize=GPrtcl_list%tsize
    nreal3 = 2*tsize+2
    call MPI_FILE_OPEN(Prtcl_WORLD,chFile, MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierr)
    call MPI_FILE_GET_SIZE(fh,FileSize,ierr)
    FileSize=FileSize+int(nreal3*real3_byte,8)*int(ATP_Opt%np_InDomain,8)
    call MPI_BARRIER(Prtcl_WORLD,ierr)
    call my_mpi_file_set_size(fh,FileSize,ierr)
    call MPI_BARRIER(Prtcl_WORLD,ierr)
    
    allocate(real3Vec(ATP_NumRestart*nreal3))
    nLeft=nlocal; pid=0
    bgn_byte=disp+int(nreal3*real3_byte,8)*int(bgn_ind,8)
    DO
      nRestart=min(nLeft,ATP_NumRestart)
      k=0
      do i=1,nRestart
        pid=pid+1
        real3Vec(k+1)=GPrtcl_SwimDir(pid);      k=k+1
        real3Vec(k+1)=GPrtcl_MoveDistance(pid); k=k+1
        real3Vec(k+1:k+tsize)=GPrtcl_LinVel(1:tsize,pid); k=k+tsize
        real3Vec(k+1:k+tsize)=GPrtcl_SwimAcc(1:tsize,pid); k=k+tsize
      enddo
      call MPI_FILE_WRITE_AT(fh,bgn_byte,real3Vec,k,real3_type,MPI_STATUS_IGNORE,ierr)
      bgn_byte=bgn_byte+int(nreal3*real3_byte,8)*int(nRestart,8)
      nLeft=nLeft-nRestart
      if(nLeft==0)exit
    ENDDO
    deallocate(real3Vec)
    call MPI_FILE_CLOSE(fh, ierr)
    call MPI_COMM_FREE(Prtcl_WORLD,ierr)
  end subroutine PIO_Write_Restart
#undef ATP_NumRestart

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PIO_Init_Visu
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PIO_Final_Visu(); implicit none 

    ! locals
    integer::iUnit,ierr
    character(len=:),allocatable::XdmfFile

    if(nrank/=0 .or. saveXDMFOnce) return
    xdmfFile = strip(ATP_Opt%ResultsDir) // "PartVisuFor" // strip(ATP_Opt%RunName) // ".xmf"
    open(newunit=iUnit, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierr)
    if(ierr/=0 .and. nrank==0) call ATPLogInfo%CheckForError(ErrT_Abort,"PIO_Final_Visu","Cannot open file:  "//XdmfFile)
    ! XDMF/XMF Tail
    write(iUnit,'(A)') '    </Grid>'
    write(iUnit,'(A)') '</Domain>'
    write(iUnit,'(A)') '</Xdmf>'
    close(iUnit,IOSTAT=ierr)
  end subroutine PIO_Final_Visu

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PIO_Dump_Visu
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PIO_Dump_Visu(this, itime); implicit none
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime

    ! locals
    type(part_io_size_vec)::pvsize 
    integer(kind=MPI_OFFSET_KIND)::disp
    character(len=:),allocatable::chFile
    type(real3),allocatable,dimension(:)::real3Vec
    integer :: ierr,fh,i,color,key,Prtcl_WORLD,nlocal,bgn_ind

    ! write xdmf file first
    if(.not.saveXDMFOnce) call this%Write_XDMF(itime)

    ! update the bgn_ind
    nlocal = GPrtcl_list%nlocal
    bgn_ind=clc_bgn_ind(nlocal)

    ! create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD,ierr)
    if(color==2) return
    
    ! begin to dump
    chFile = strip(ATP_Opt%ResultsDir) // "PartVisuFor" // strip(ATP_Opt%RunName) // int2str(itime,10)
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_BARRIER(Prtcl_WORLD,ierr)
    call my_mpi_file_set_size(fh,0_8,ierr)  ! guarantee overwriting
    call MPI_BARRIER(Prtcl_WORLD,ierr)
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
    call MPI_FILE_CLOSE(fh, ierr) 
    call MPI_COMM_FREE(Prtcl_WORLD, ierr)
  end subroutine PIO_Dump_Visu

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clc_bgn_ind
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  function clc_bgn_ind(nlocal) result(bgn_ind); implicit none
    integer,intent(in)::nlocal
    integer::bgn_ind

    ! locals
    integer::end_ind,ierr,SRstatus(MPI_STATUS_SIZE)
  
    bgn_ind=0
    if(nproc<=1) return
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(nrank==0)THEN
      end_ind=bgn_ind+nlocal
      call MPI_SEND(end_ind, 1,int_type, nrank+1,0,MPI_COMM_WORLD,ierr)
    ELSEIF(nrank /= nproc-1) THEN
      call MPI_RECV(bgn_ind, 1,int_type, nrank-1,0,MPI_COMM_WORLD,SRstatus,ierr)
      end_ind=bgn_ind+nlocal
      call MPI_SEND(end_ind, 1,int_type, nrank+1,0,MPI_COMM_WORLD,ierr)
    ELSE
      call MPI_RECV(bgn_ind, 1,int_type, nrank-1,0,MPI_COMM_WORLD,SRstatus,ierr)
    ENDIF
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end function clc_bgn_ind
end module ap_IOAndVisu
