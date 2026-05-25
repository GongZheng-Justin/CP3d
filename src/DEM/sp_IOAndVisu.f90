#include "definitions_inc.f90"
module sp_IOAndVisu
  use MPI
  use mc_TypeDef
  use mc_LogInfo
#if defined(CFDDEM) || defined(CFDACM)
  use mc_Decomp2d,only: nrank,nproc
#endif
  use sp_Comm
  use sp_Property
  use sp_Variables
  use sp_CL_and_CF
  use sp_Decomp_2d
  use sp_Parameters
  implicit none
  private

  integer,parameter:: IK = 4
  integer::Prev_BackUp_itime = 53456791
  logical::saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel
  logical::save_LinAcc,save_Theta,save_RotVel,save_RotAcc,save_CntctForce,save_Torque
#ifdef CFDACM
  logical::save_HighSt
#endif

  type part_io_size_vec
    integer,dimension(1)::sizes
    integer,dimension(1)::subsizes
    integer,dimension(1)::starts
  end type part_io_size_vec
  type part_io_size_mat
    integer,dimension(2)::sizes
    integer,dimension(2)::subsizes
    integer,dimension(2)::starts
  end type part_io_size_mat

  ! useful interfaces
  interface Prtcl_dump
    module procedure Prtcl_dump_int_vector,  Prtcl_dump_int_matrix
    module procedure Prtcl_dump_real_vector, Prtcl_dump_real3_vector
#ifdef CFDDEM
    module procedure Prtcl_dump_real3_matrix
#endif  
  end interface Prtcl_dump
  
  public:: Prtcl_Dump_Visu, Prtcl_Final_Visu, Prtcl_Delete_Prev_Restart, Prtcl_Write_Restart
  public:: Prtcl_Init_Visu, Prtcl_Read_Restart, Prtcl_ReadFixedCoord, Prtcl_Restart_ContactList
#ifdef CFDDEM
  public:: Prtcl_ReadFixedRestart, Prtcl_WriteFixedRestart
#endif

contains
#include "Prtcl_Dump_MPI_inc.f90"

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_Init_Visu
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Prtcl_Init_Visu(chFile,iStage); implicit none
    character(len=*),intent(in)::chFile
    integer,intent(in)::iStage
    
    ! locals
    integer::iUnit,ierr,nfld,ifld,iDivided
    character(len=:),allocatable::XdmfFile
#ifdef CFDACM
    NAMELIST /PrtclVisuOption/ saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel,  &
                               save_LinAcc,save_Theta,save_RotVel,save_RotAcc,save_CntctForce,save_Torque,save_HighSt
#else
    NAMELIST /PrtclVisuOption/ saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel,  &
                               save_LinAcc,save_Theta,save_RotVel,save_RotAcc,save_CntctForce,save_Torque
#endif
  
    if(iStage==1) then
      open(newunit=iUnit, file=chFile,status='old',form='formatted',IOSTAT=ierr)
      if(ierr/=0)call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_Init_Visu", "Cannot open file: "//strip(chFile))
      read(iUnit, nml=PrtclVisuOption)
      if(nrank==0)write(DEMLogInfo%iUnit, nml=PrtclVisuOption)
      close(iUnit,IOSTAT=ierr)
      return
    endif

    ! Write XDMF file
    if(nrank/=0) return
    XdmfFile = strip(DEM_Opt%ResultsDir)//"PartVisuFor"//strip(DEM_Opt%RunName)//".xmf"
    open(newunit=iUnit, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierr)
    if(ierr /= 0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_Init_Visu","Cannot open file:  "//XdmfFile)
    write(iUnit,'(A)') '<?xml version="1.0" ?>'
    write(iUnit,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(iUnit,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(iUnit,'(A)') '<Domain>'

    ! Time series
#if defined(CFDDEM) || defined(CFDACM)
    iDivided = icouple
#else
    iDivided = 1
#endif
    nfld = (DEM_Opt%ilast - DEM_Opt%ifirst +1)/DEM_Opt%SaveVisu  + 1
    write(iUnit,'(A)')'  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    write(iUnit,'(A)')'    <Time TimeType="List">'
    write(iUnit,'(A,I6,A)')'      <DataItem Format="XML" NumberType="Int" Dimensions="',nfld,'">' 
    write(iUnit,'(A)',advance='no')'        '
    do ifld =1, nfld
      if(mod(ifld,10)==0) then
        write(iUnit,'(I10)') ((ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1)/iDivided
        if(ifld < nfld) write(iUnit,'(A)',advance='no') '        '
      else
        write(iUnit,'(I9)',advance='no') ((ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1)/iDivided
      endif
    enddo
    if(mod(nfld,10) /=0) write(iUnit,*)' '
    write(iUnit,'(A)') '      </DataItem>'
    write(iUnit,'(A)') '    </Time>'
    close(iUnit, IOSTAT=ierr)
    if(.not. saveXDMFOnce) return

    do ifld = 1,nfld
      call Write_XDMF( ((ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1)/iDivided )
    enddo

    ! XDMF/XMF Tail
    open(newunit=iUnit, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierr)
    write(iUnit,'(A)')'  </Grid>'
    write(iUnit,'(A)')'</Domain>'
    write(iUnit,'(A)')'</Xdmf>'
    close(iUnit, IOSTAT=ierr)
  end subroutine Prtcl_Init_Visu

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
  subroutine Write_XDMF(itime); implicit none 
    integer,intent(in)::itime

    ! locals
    integer:: iUnit,ierr,np,dims,iprec
    integer(kind=MPI_OFFSET_KIND)::disp
    character(len=:),allocatable::XdmfFile

    if(nrank/=0) return 
    np=DEM_Opt%np_InDomain
    XdmfFile = strip(DEM_Opt%ResultsDir)//"PartVisuFor"//strip(DEM_Opt%RunName)//".xmf"
    open(newunit=iUnit, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierr)
    if(ierr/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Write_XDMF","Cannot open file: "//XdmfFile)
    endif

    disp = 0_MPI_OFFSET_KIND
    XdmfFile = "PartVisuFor"//strip(DEM_Opt%RunName)
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
    IF(save_Diameter) THEN
      dims=1; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"Diameter","Scalar","Float",disp)
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
    IF(save_LinAcc) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"LinAcc","Vector","Float",disp)
    ENDIF
    IF(save_Theta) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"Theta","Vector","Float",disp)
    ENDIF
    IF(save_RotVel) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"RotVel","Vector","Float",disp)
    ENDIF
    IF(save_RotAcc) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"RotAcc","Vector","Float",disp)
    ENDIF
    IF(save_CntctForce) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"CntctForce","Vector","Float",disp)
    ENDIF
    IF(save_Torque) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"Torque","Vector","Float",disp)
    ENDIF
#ifdef CFDACM
    IF(save_HighSt) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,"IsHighSt","Scalar","Int",disp)
    ENDIF
#endif    
    write(iUnit,'(A)')'    </Grid>'
    close(iUnit,IOSTAT=ierr)
  end subroutine Write_XDMF

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
  ! Prtcl_Delete_Prev_Restart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Prtcl_Delete_Prev_Restart(itime); implicit none
    integer:: itime

    ! locals
    integer::iUnit,ierr
    character(len=:),allocatable::chFile

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(nrank/=0) return
    chFile = strip(DEM_Opt%RestartDir)//"RestartFor"//strip(DEM_Opt%RunName)//int2str(Prev_BackUp_itime,10)
    open(newunit=iUnit,file=chFile,IOSTAT=ierr)
    close(unit=iUnit,status='delete',IOSTAT=ierr)

#ifdef CFDDEM
    if(Is_clc_FluidAcc .or. (is_clc_Basset .and. is_clc_Basset_Fixed)) then
      chFile = strip(DEM_Opt%RestartDir)//"FixedSpheresRestart"//int2str(Prev_BackUp_itime,10)
      open(newunit=iUnit,file=chFile,IOSTAT=ierr)
      close(unit=iUnit,status='delete',IOSTAT=ierr)
    endif
    Prev_BackUp_itime = itime/icouple
#elif CFDACM
    Prev_BackUp_itime = itime/icouple
#else
    Prev_BackUp_itime = itime
#endif
  end subroutine Prtcl_Delete_Prev_Restart

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_ReadFixedCoord
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
#define DEM_NumRead 100
  subroutine Prtcl_ReadFixedCoord(); implicit none
 
    ! locals
    type(real3)::real3t
    integer(kind=8)::byte_total,disp
    character(len=:),allocatable::chFile
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd,radius,diam
    type(real4),allocatable,dimension(:)::real4Vec
    integer,allocatable,dimension(:)::nP_in_bin,nP_in_bin_reduce,IntVec
    integer::i,pid,ierr,nLeft,nRead,numPrtclFix,nfix,nfix_sum,pType,nfixNew,pbyte,iUnit,Init_SpheresCoord_Type
    real(RK),dimension(5,DEM_NumRead)::FixedPrtcl
    logical::lexist

    ! FixedPrtcl_type
    pbyte=real_byte*5

    numPrtclFix = DEM_Opt%numPrtclFix
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd
#ifdef ChanBraunJFM2011
    ySt=-yEd
#endif

    ! The data storage sequence in file "FixedSpheresCoord.dat" is as follow:
    ! Particle1: Position(real3 type), Diameter(real type), Prtcl_Type(real type)
    ! Particle2: Position(real3 type), Diameter(real type), Prtcl_Type(real type) ..
    Init_SpheresCoord_Type = 0
    chFile = strip(DEM_Opt%RestartDir)//"FixedSpheresCoord.bin"
    inquire(file=chFile, exist=lexist)   ! Try to find the binary file firstly.
    if(.not. lexist) then
      chFile = strip(DEM_Opt%RestartDir)//"FixedSpheresCoord.txt"
      inquire(file=chFile, exist=lexist) ! Try to find the plain text file secondly.
      if(.not. lexist) then
        if(nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedCoord","Cannot open particle Coordinate file")
      else
        Init_SpheresCoord_Type = 1
      endif
    endif
    IF(Init_SpheresCoord_Type == 0) THEN
      chFile = strip(DEM_Opt%RestartDir)//"FixedSpheresCoord.bin"
      open(newunit=iUnit,file=chFile,status='old',form='unformatted',access='stream',action='read',position='append',IOSTAT=ierr)
      if(ierr/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedCoord","Cannot open file: "//chFile)
      inquire(unit=iUnit,Pos=disp); disp=disp-1_8
      rewind(unit=iUnit,IOSTAT=ierr)
      byte_total= int(pbyte,8)*int(numPrtclFix,8)
      if(disp/=byte_total .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedCoord","file byte wrong")
    ELSE
      chFile = strip(DEM_Opt%RestartDir)//"FixedSpheresCoord.txt"
      open(newunit=iUnit,file=chFile,status='old',form='formatted',action='read',IOSTAT=ierr)
      if(ierr/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedCoord","Cannot open file: "//chFile)    
    ENDIF

    nfix=0; nLeft=numPrtclFix; pid=0; disp= 1_8
    allocate(nP_in_bin(DEM_Opt%numPrtcl_Type));        nP_in_bin=0
    allocate(nP_in_bin_reduce(DEM_Opt%numPrtcl_Type)); nP_in_bin_reduce=0
    DO 
      nRead=min(nLeft,DEM_NumRead)
      IF(Init_SpheresCoord_Type == 0) THEN
        read(iUnit,pos=disp,IOSTAT=ierr) FixedPrtcl(:,1:nRead)
      ELSE
        do i=1,nRead
          read(unit=iUnit,fmt=*,IOSTAT=ierr) FixedPrtcl(:,i)
          if(ierr /= 0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedCoord", "Read plain text wrong") 
        enddo
      ENDIF
      disp=disp+int(pbyte,8)*int(nRead,8)
      do i=1,nRead
        pid=pid+1
        real3t%x=FixedPrtcl(1,i)
        real3t%y=FixedPrtcl(2,i)
        real3t%z=FixedPrtcl(3,i)
        diam    =FixedPrtcl(4,i)
        pType   =nint(FixedPrtcl(5,i))
        radius= DEMProperty%Prtcl_PureProp(pType)%Radius
        if( abs(2.0_RK*radius/diam -1.0_RK)>1.0E-6 .and. nrank==0 )then
          call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedCoord","Diameter not coordinate")
        endif
        if(real3t%x< xSt .or. real3t%y< ySt .or. real3t%z< zSt .or. &
           real3t%x>=xEd .or. real3t%y>=yEd .or. real3t%z>=zEd) cycle
        if(nfix>=GPrtcl_list%mlocalFix)  then
          nfixNew= int(1.2_RK*real(nfix,kind=RK))
          nfixNew= max(nfixNew, nfix+1)
          nfixNew= min(nfixNew,DEM_Opt%numPrtclFix)
          GPrtcl_list%mlocalFix= nfixNew

          call move_alloc(GPFix_id,IntVec)
          allocate(GPFix_id(nfixNew))
          GPFix_id(1:nfix)=IntVec
          call move_alloc(GPFix_pType,IntVec)
          allocate(GPFix_pType(nfixNew))
          GPFix_pType(1:nfix)=IntVec
          deallocate(IntVec)
          call move_alloc(GPFix_PosR,real4Vec)
          allocate(GPFix_PosR(nfixNew))
          GPFix_PosR(1:nfix)=real4Vec
          deallocate(real4Vec)
        endif
        nfix=nfix+1
         
        GPFix_id(nfix)     = pid  + DEM_Opt%numPrtcl         ! NOTE HERE
        GPFix_pType(nfix)  = pType
        GPFix_PosR(nfix)   = real3t
        GPFix_PosR(nfix)%w = radius
        nP_in_bin(pType)   = nP_in_bin(pType)+1
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    close(iUnit,IOSTAT=ierr)

    call MPI_ALLREDUCE(nP_in_bin, nP_in_bin_reduce,DEM_Opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    DEMProperty%nPrtcl_in_Bin= DEMProperty%nPrtcl_in_Bin+ nP_in_bin_reduce
    deallocate(nP_in_bin,nP_in_bin_reduce)

    if(nfix>0) then
      call move_alloc(GPFix_id,IntVec)
      allocate(GPFix_id(nfix))
      GPFix_id=IntVec(1:nfix)
      call move_alloc(GPFix_pType,IntVec)
      allocate(GPFix_pType(nfix))
      GPFix_pType=IntVec(1:nfix)
      deallocate(IntVec)
      call move_alloc(GPFix_PosR,real4Vec)
      allocate(GPFix_PosR(nfix))
      GPFix_PosR=real4Vec(1:nfix)
      deallocate(real4Vec)
#ifdef CFDDEM
      allocate(GPFix_VFluid(2,nfix));GPFix_VFluid=zero_r3
#endif
    else
      deallocate(GPFix_id,GPFix_pType,GPFix_PosR)
    endif

    GPrtcl_list%mlocalFix = nfix
    call MPI_REDUCE(nfix,nfix_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(nfix_sum/= numPrtclFix .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedCoord"," nfix_sum/= numPrtclFix " )
    endif
  end subroutine Prtcl_ReadFixedCoord
#undef DEM_NumRead

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_Restart_ContactList
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
#define DEM_NumRead 2000
  subroutine Prtcl_Restart_ContactList(); implicit none

    ! locals
    type(real4)::real4t
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd
    character(len=:),allocatable::chFile
    integer,allocatable,dimension(:)::ncvVec
    integer,dimension(:),allocatable::CntctVec
    type(real3),allocatable,dimension(:)::PosVec
    type(real4),dimension(:),allocatable::TanDel_Un
    integer::intvec(2),i,j,ncv,nCntctTotal,nLeft,nRead
    integer::itime,iUnit,ierr,nlocal,np,nreal3,pbyte,ncvMax
    integer(kind=8)::disp,disp_pos,disp_ncv,disp_CL,disp_TanStart,disp_Tan

    itime = DEM_Opt%ifirst - 1
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd
    
    ! Begin to read Restart_Contact_List
#if defined(CFDDEM) || defined(CFDACM)
    chFile = strip(DEM_Opt%RestartDir)//"RestartFor"//strip(DEM_Opt%RunName)//int2str(itime/icouple,10)
#else
    chFile = strip(DEM_Opt%RestartDir)//"RestartFor"//strip(DEM_Opt%RunName)//int2str(itime,10)
#endif
    open(newunit=iUnit,file=chFile,status='old',form='unformatted',access='stream',action='read',IOSTAT=ierr)
    if(ierr/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_Restart_ContactList","Cannot open file: "//chFile)
    np= DEM_Opt%np_InDomain; disp=1_8 + int_byte  ! firstly skip the np_InDomain
    read(iUnit,pos=disp,IOSTAT=ierr) nCntctTotal; disp=disp+int_byte;
    
    nreal3 = 2*(1+GPrtcl_list%tsize+GPrtcl_list%rsize)
#ifdef CFDDEM
    nreal3=nreal3 +3      !(GPrtcl_FpForce, GPrtcl_linVelOld, GPrtcl_VFluid(1,:))
    if(Is_clc_Basset) nreal3= nreal3+ GPrtcl_BassetSeq%nDataLen
    disp=disp+2*int_byte; !skip Is_clc_Basset and HistoryStage
#endif
#ifdef CFDACM
    nreal3=nreal3+4
#endif
    pbyte=int_byte*3 + nreal3*real3_byte  ! corresponding to the subroutine 'Prtcl_Write_Restart'    

    allocate(PosVec(DEM_NumRead),ncvVec(DEM_NumRead))
    
    ! Determine the maxinum of ncvVec
    nLeft=np; ncvMax=0
    disp_ncv = disp + pbyte*np
    DO
      nRead=min(nLeft,DEM_NumRead)
      read(iUnit,pos=disp_ncv,IOSTAT=ierr)ncvVec(1:nRead)
      disp_ncv=disp_ncv+int(int_byte,8)*int(nRead,8)
      do i=1,nRead
        if(ncvMax<ncvVec(i)) ncvMax=ncvVec(i)
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    if(ncvMax==0) ncvMax=1
    allocate(TanDel_Un(ncvMax))
    allocate(CntctVec(ncvMax))
        
    nlocal=0; nLeft=np
    disp_pos = disp
    disp_ncv = disp + pbyte*np
    disp_CL  = disp_ncv+ int_byte*np
    disp_TanStart= disp_CL + 2*nCntctTotal*int_byte
    DO
      nRead=min(nLeft,DEM_NumRead)
      read(iUnit,pos=disp_pos,IOSTAT=ierr)PosVec(1:nRead)
      read(iUnit,pos=disp_ncv,IOSTAT=ierr)ncvVec(1:nRead)
      disp_pos=disp_pos+int(real3_byte,8)*int(nRead,8)
      disp_ncv=disp_ncv+int(int_byte,8)*int(nRead,8)
      do i=1,nRead
        ! ncv: number of particles/walls which have overlap with this particle
        ncv=ncvVec(i)
        if(PosVec(i)%x< xSt .or. PosVec(i)%y< ySt .or. PosVec(i)%z< zSt  .or. &
           PosVec(i)%x>=xEd .or. PosVec(i)%y>=yEd .or. PosVec(i)%z>=zEd) then
          disp_CL=disp_CL+int_byte*2*ncv
        else
          nlocal= nlocal+ 1
          do j=1,ncv
            read(iUnit,pos=disp_CL,IOSTAT=ierr)intvec(1:2); disp_CL=disp_CL+int_byte*2
            disp_Tan = disp_TanStart + real4_byte*(intvec(2)-1)
            read(iUnit,pos=disp_Tan,IOSTAT=ierr)real4t
            CntctVec(j) = intvec(1)
            TanDel_Un(j)= real4t
          enddo
          if(ncv>0) call GPPW_CntctList%Add_RestartCntctlink(nlocal,ncv,CntctVec,TanDel_Un)
        endif
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    deallocate(PosVec, ncvVec)
    deallocate(TanDel_Un, CntctVec)
    close(iUnit,IOSTAT=ierr)
  end subroutine Prtcl_Restart_ContactList
#undef DEM_NumRead

#ifdef CFDDEM
#define DEM_NumWrite 100
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_WriteFixedRestart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Prtcl_WriteFixedRestart(itime); implicit none
    integer,intent(in)::itime

    ! locals
    logical::IsWriteBasset
    type(part_io_size_vec)::pvsize
    type(part_io_size_mat)::pmsize
    character(len=:),allocatable::chFile
    real(RK),dimension(5,DEM_NumWrite)::FixedPrtclOut
    integer(kind=MPI_OFFSET_KIND)::disp,disp_pos,FileSize
    integer::bgn_ind,color,key,ierr,PrtclFix_WORLD,fh,nLeft,nWrite,i,pbyte,pid,intVec(2)

    if(DEM_Opt%numPrtclFix<1) return
    IsWriteBasset=.false.
    if(is_clc_Basset .and. is_clc_Basset_Fixed) IsWriteBasset=.true.
    if((.not.Is_clc_FluidAcc) .and. (.not.IsWriteBasset)) return

    ! Create and initialize file
    if(IsWriteBasset) then
      intVec=[1,GPrtcl_BassetSeq%HistStageFix]
    else
      intVec=[0,0]
    endif
    chFile = strip(DEM_Opt%RestartDir)//"FixedSpheresRestart"//int2str(itime/icouple,10)
    if(nrank==0) then
      open(newunit=fh,file=chFile,status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierr)
      write(fh,pos=1_8,IOSTAT=ierr) intVec(1:2)
      close(fh,IOSTAT=ierr)
    endif
    disp=int_byte*2_8
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    ! Create the Prtcl_GROUP
    bgn_ind= clc_bgn_ind(GPrtcl_list%mlocalFix)
    color = 1; key=nrank
    if(GPrtcl_list%mlocalFix<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,PrtclFix_WORLD,ierr)
    if(color==2) return

    ! PosD, id, pType
    pbyte= real_byte*5
    FileSize=disp+int(pbyte,8)*int(DEM_Opt%numPrtclFix,8)
    call MPI_FILE_OPEN(PrtclFix_WORLD, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_BARRIER(PrtclFix_WORLD,ierr)
    call my_mpi_file_set_size(fh,FileSize,ierr)
    call MPI_BARRIER(PrtclFix_WORLD,ierr)
      
    nLeft=GPrtcl_list%mlocalFix;  pid=0
    disp_pos=disp+int(pbyte,8)*int(bgn_ind,8)
    DO 
      nWrite=min(nLeft,DEM_NumWrite)
      do i=1,nWrite
        pid=pid+1
        FixedPrtclOut(1,i)=GPFix_PosR(pid)%x
        FixedPrtclOut(2,i)=GPFix_PosR(pid)%y
        FixedPrtclOut(3,i)=GPFix_PosR(pid)%z
        FixedPrtclOut(4,i)=real(GPFix_id(pid),RK)
        FixedPrtclOut(5,i)=real(GPFix_pType(pid),RK)
      enddo
      call MPI_FILE_WRITE_AT(fh,disp_pos,FixedPrtclOut,5*nWrite,real_type,MPI_STATUS_IGNORE, ierr)
      disp_pos=disp_pos+int(pbyte,8)*int(nWrite,8)
      nLeft=nLeft-nWrite
      if(nLeft==0)exit
    ENDDO
    disp=disp+int(pbyte,8)*int(DEM_Opt%numPrtclFix,8)
    call MPI_BARRIER(PrtclFix_WORLD,ierr)

    ! Write GPFix_Vfluid and GPFix_BassetData
    if(Is_clc_FluidAcc) then
      pvsize%sizes(1)   = DEM_Opt%numPrtclFix
      pvsize%subsizes(1)= GPrtcl_list%mlocalFix
      pvsize%starts(1)  = bgn_ind
      call Prtcl_dump(fh, disp, GPFix_Vfluid(1,:), pvsize)
    endif
    if(IsWriteBasset) then
      pmsize%sizes(1)   = GPrtcl_BassetSeq%nDataLen;  pmsize%sizes(2)   = DEM_Opt%numPrtclFix
      pmsize%subsizes(1)= GPrtcl_BassetSeq%nDataLen;  pmsize%subsizes(2)= GPrtcl_list%mlocalFix
      pmsize%starts(1)  = 0 ;                         pmsize%starts(2)  = bgn_ind
      call Prtcl_dump(fh, disp, GPFix_BassetData, pmsize)
    endif
    call MPI_FILE_CLOSE(fh, ierr) 
    call MPI_COMM_FREE(PrtclFix_WORLD,ierr)
  end subroutine Prtcl_WriteFixedRestart
#undef DEM_NumWrite

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_ReadFixedRestart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
#define DEM_NumRead 500
  subroutine Prtcl_ReadFixedRestart(); implicit none

    ! locals
    type(real3)::real3t
    logical::IsReadBasset
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd
    character(len=:),allocatable::chFile
    integer(kind=8)::byte_total,disp,disp_Fluid
    type(real3),allocatable,dimension(:)::Real3Vec
    integer,dimension(:),allocatable::nP_in_bin,nP_in_bin_reduce,idVec,IntVec
    integer::i,pbyte,pid,iPos,pType,iUnit,ierr,numPrtclFix,nfix,nfixNew,nRead,nLeft,nfix_sum,int_t(2)
    real(RK)::FixedPrtclIn(5,DEM_NumRead),FixedPrtclOne(5)

    IsReadBasset=.false.
    if(is_clc_Basset .and. is_clc_Basset_Fixed) IsReadBasset=.true.
    if((.not.Is_clc_FluidAcc) .and. (.not.IsReadBasset)) then
      call Prtcl_ReadFixedCoord()
      return
    endif

    numPrtclFix = DEM_Opt%numPrtclFix
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd

    ! PosD, id, pType
    pbyte=real_byte*5

    chFile = strip(DEM_Opt%RestartDir)//"FixedSpheresRestart"//int2str((DEM_Opt%ifirst - 1)/icouple, 10)
    open(newunit=iUnit,file=chFile,status='old',form='unformatted',access='stream',action='read',position='append',IOSTAT=ierr)
    if(ierr/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedRestart","Cannot open file:"//chFile)

    byte_total= int_byte*2+pbyte*numPrtclFix ! 'int_byte*2' corresponds to HistStageFix in  Prtcl_WriteFixedRestart
    if(is_clc_FluidAcc) byte_total= byte_total+ real3_byte*numPrtclFix
    if(IsReadBasset) byte_total= byte_total+ real3_byte*numPrtclFix*GPrtcl_BassetSeq%nDataLen
    inquire(unit=iUnit,Pos=disp); disp=disp-1_8
    rewind(unit=iUnit,IOSTAT=ierr)
    if(disp/=byte_total .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedCoordRestart","file byte wrong")

    disp=1_8
    read(iUnit,pos=disp,IOSTAT=ierr)int_t(1:2); disp=disp+int_byte*2
    if(IsReadBasset) then
      if(int_t(1)/= 1) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedRestart"," Is_clc_Basset_Fixed Wrong 1" )
      GPrtcl_BassetSeq%HistStageFix= int_t(2)
    else
      if(int_t(1)/= 0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedRestart"," Is_clc_Basset_Fixed Wrong 2" )
    endif   

    nfix=0; nLeft=numPrtclFix; pid=0
    allocate(idVec(GPrtcl_list%mlocalFix))           
    allocate(nP_in_bin(DEM_Opt%numPrtcl_Type));        nP_in_bin=0
    allocate(nP_in_bin_reduce(DEM_Opt%numPrtcl_Type)); nP_in_bin_reduce=0    
    DO 
      nRead=min(nLeft,DEM_NumRead)
      read(iUnit,pos=disp,IOSTAT=ierr)FixedPrtclIn(:,1:nRead)
      disp=disp+int(pbyte,8)*int(nRead,8)
      do i=1,nRead
        pid=pid+1
        real3t%x=FixedPrtclIn(1,i)
        real3t%y=FixedPrtclIn(2,i)
        real3t%z=FixedPrtclIn(3,i)
        if(real3t%x< xSt .or. real3t%y< ySt .or. real3t%z< zSt .or. &
           real3t%x>=xEd .or. real3t%y>=yEd .or. real3t%z>=zEd) cycle
        if(nfix>=GPrtcl_list%mlocalFix)  then
          nfixNew= int(1.2_RK*real(nfix,kind=RK))
          nfixNew= max(nfixNew, nfix+1)
          nfixNew= min(nfixNew,numPrtclFix)
          GPrtcl_list%mlocalFix= nfixNew
          call move_alloc(idVec,IntVec)
          allocate(idVec(nfixNew))
          idVec(1:nfix)=IntVec
          deallocate(IntVec)
        endif
        nfix=nfix+1
        idVec(nfix)=pid
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    call MPI_ALLREDUCE(nP_in_bin, nP_in_bin_reduce,DEM_Opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    DEMProperty%nPrtcl_in_Bin= DEMProperty%nPrtcl_in_Bin+ nP_in_bin_reduce
    deallocate(nP_in_bin,nP_in_bin_reduce)

    GPrtcl_list%mlocalFix = nfix
    call MPI_REDUCE(nfix,nfix_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(nfix_sum/= numPrtclFix .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_ReadFixedRestart"," nfix_sum/= numPrtclFix " )
    endif

    if(nfix>0) then
      call move_alloc(idVec,IntVec)
      allocate(idVec(nfix))
      idVec=IntVec(1:nfix)
      deallocate(IntVec)
      deallocate(GPFix_id);   allocate(GPFix_id(nfix))
      deallocate(GPFix_pType);allocate(GPFix_pType(nfix))
      deallocate(GPFix_PosR); allocate(GPFix_PosR(nfix))
      allocate(GPFix_VFluid(2,nfix));GPFix_VFluid=zero_r3
      if(IsReadBasset) then
        allocate(GPFix_BassetData(GPrtcl_BassetSeq%nDataLen, nfix), Stat= ierr)
        if(ierr/=0) call DEMLogInfo%checkForError(ErrT_Abort,"Prtcl_ReadFixedRestart","allocate wrong2")
        GPFix_BassetData=zero_r3
      endif
    else
      if(allocated(IntVec))deallocate(IntVec)
      deallocate(GPFix_id,GPFix_pType,GPFix_PosR)
    endif

    allocate(Real3Vec(GPrtcl_BassetSeq%nDataLen))
    do pid=1,nfix
      iPos=idVec(pid)-1
      disp=1_8+int_byte*2+int(iPos,8)*int(pbyte,8)
      read(iUnit,pos=disp,IOSTAT=ierr)FixedPrtclOne
      GPFix_id(pid)   = nint(FixedPrtclOne(4))
      pType = nint(FixedPrtclOne(5))
      GPFix_pType(pid)= pType
      GPFix_PosR(pid)%x= FixedPrtclOne(1)
      GPFix_PosR(pid)%y= FixedPrtclOne(2)
      GPFix_PosR(pid)%z= FixedPrtclOne(3)
      GPFix_PosR(pid)%w= DEMProperty%Prtcl_PureProp(pType)%Radius
      if(Is_clc_FluidAcc) then
        disp_Fluid= 1_8+int_byte*2+int(numPrtclFix,8)*int(pbyte,8) +int(iPos,8)*int(real3_byte,8)
        read(iUnit,pos=disp_Fluid,IOSTAT=ierr)real3t
        GPFix_Vfluid(1,pid)=real3t
      endif
      if(IsReadBasset) then
        disp_Fluid= 1_8+int_byte*2+int(numPrtclFix,8)*int(pbyte+real3_byte,8) +int(iPos,8)*int(real3_byte*GPrtcl_BassetSeq%nDataLen,8)
        read(iUnit,pos=disp_Fluid,IOSTAT=ierr)Real3Vec(1:GPrtcl_BassetSeq%nDataLen)
        GPFix_BassetData(:,pid)= Real3Vec
      endif
    enddo
    deallocate(idVec,Real3Vec)
    close(iUnit,IOSTAT=ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end subroutine Prtcl_ReadFixedRestart
#undef DEM_NumRead
#endif

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_Read_Restart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
#define DEM_NumRead 2000
  subroutine Prtcl_Read_Restart(); implicit none

    ! locals
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd
    character(len=:),allocatable::chFile
    integer,allocatable,dimension(:):: nP_in_bin
    type(real3),allocatable,dimension(:)::real3Vec,PosVec
    integer::iUnit,tsize,rsize,nlocal_sum,nreal3,itype
    integer::itime,ierr,nlocal,i,k,np,nLeft,nRead,int_t(3)
    integer(kind=MPI_OFFSET_KIND)::disp,disp_pos,disp_int,disp_real3
    
    itime = DEM_Opt%ifirst - 1
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd

    ! Begin to read Restart_file
#if defined(CFDDEM) || defined(CFDACM)
    chFile = strip(DEM_Opt%RestartDir)//"RestartFor"//strip(DEM_Opt%RunName)//int2str(itime/icouple,10)
#else
    chFile = strip(DEM_Opt%RestartDir)//"RestartFor"//strip(DEM_Opt%RunName)//int2str(itime,10)
#endif
    open(newunit=iUnit,file=chFile,status='old',form='unformatted',access='stream',action='read',IOSTAT=ierr)
    if(ierr/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_Read_Restart","Cannot open file: "//chFile)
    disp =1_8; read(iUnit,pos=disp,IOSTAT=ierr)np; disp=disp+int_byte
    if(np>DEM_Opt%numPrtcl .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_Read_Restart: "," np_InDomain > numPrtcl " )
    endif
    DEM_Opt%np_InDomain = np
    disp=disp+int_byte   ! skip nCntctTotal
#ifdef CFDDEM
    read(iUnit,pos=disp,IOSTAT=ierr)int_t(1:2); disp=disp+2*int_byte;
    if(Is_clc_Basset) then
      if(int_t(1)/= 1) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_Read_Restart"," Is_clc_Basset Wrong 1" )
      GPrtcl_BassetSeq%HistoryStage= int_t(2)
    else
      if(int_t(1)/= 0) call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_Read_Restart"," Is_clc_Basset Wrong 2" )
    endif
#endif

    tsize=GPrtcl_list%tsize; rsize=GPrtcl_list%rsize;
    nreal3 = 2*(1+tsize+rsize)
#ifdef CFDDEM
    nreal3=nreal3+3
    if(Is_clc_Basset) nreal3=nreal3+GPrtcl_BassetSeq%nDataLen
#endif
#ifdef CFDACM
    nreal3=nreal3+4
#endif
    nreal3=nreal3-1   ! "-1" corresponds to GPrtcl_Pos
    allocate(real3Vec(nreal3),PosVec(DEM_NumRead))
    allocate(nP_in_bin(DEM_Opt%numPrtcl_Type)); nP_in_bin=0

    nlocal=0; nLeft=np
    disp_pos  = disp
    disp_int  = disp_pos+ real3_byte*np
    disp_real3= disp_int+ int_byte*np*3
    DO
      nRead=min(nLeft,DEM_NumRead)
      read(iUnit,pos=disp_pos,IOSTAT=ierr)PosVec(1:nRead)
      disp_pos=disp_pos+int(real3_byte,8)*int(nRead,8)
      do i=1,nRead
        if(PosVec(i)%x>=xSt .and. PosVec(i)%x< xEd .and. PosVec(i)%y>=ySt .and. &
           PosVec(i)%y< yEd .and. PosVec(i)%z>=zSt .and. PosVec(i)%z<zEd) then
          if(nlocal>=GPrtcl_list%mlocal)  call GPrtcl_list%ReallocatePrtclVar(nlocal)
          nlocal=nlocal+1

          read(iUnit,pos=disp_int,IOSTAT=ierr)int_t(1:3)
          GPrtcl_id(nlocal)=int_t(1)      ! id
          itype=int_t(2)
          GPrtcl_pType(nlocal)=itype      ! pType
          nP_in_bin(itype)= nP_in_bin(itype)+1
          GPrtcl_UsrMark(nlocal)=int_t(3) ! Usr_Mark

          GPrtcl_PosR(nlocal)= PosVec(i)  ! PosR
          GPrtcl_PosR(nlocal)%w=DEMProperty%Prtcl_PureProp(itype)%Radius
          k=0;
          read(iUnit,pos=disp_real3,IOSTAT=ierr)real3Vec(1:nreal3)
          GPrtcl_LinVel(1:tsize,nlocal)  =real3Vec(k+1:k+tsize); k=k+tsize ! LinVec
          GPrtcl_LinAcc(1:tsize,nlocal)  =real3Vec(k+1:k+tsize); k=k+tsize ! LinAcc
          GPrtcl_theta(nlocal)           =real3Vec(k+1);         k=k+1     ! Theta
          GPrtcl_RotVel(1:rsize,nlocal)  =real3Vec(k+1:k+rsize); k=k+rsize ! RotVel
          GPrtcl_RotAcc(1:rsize,nlocal)  =real3Vec(k+1:k+rsize); k=k+rsize ! RotAcc
#ifdef CFDACM
          GPrtcl_FpForce(nlocal)         =real3Vec(k+1);         k=k+1
          GPrtcl_FpTorque(nlocal)        =real3Vec(k+1);         k=k+1 
          GPrtcl_FluidIntOld(1:2,nlocal) =real3Vec(k+1:k+2);     k=k+2     ! FluidIntOld
#endif
#ifdef CFDDEM     
          GPrtcl_FpForce(nlocal)         =real3Vec(k+1);         k=k+1    
          GPrtcl_linVelOld(nlocal)       =real3Vec(k+1);         k=k+1
          GPrtcl_VFluid(1,nlocal)        =real3Vec(k+1);         k=k+1 
          if(Is_clc_Basset) then
            itype=GPrtcl_BassetSeq%nDataLen
            GPrtcl_BassetData(1:itype,nlocal)=real3Vec(k+1:k+itype);     k=k+itype
          endif
#endif
        endif
        disp_int  = disp_int  + int_byte*3
        disp_real3= disp_real3+ real3_byte*nreal3       
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    deallocate(PosVec,Real3Vec)
    call MPI_ALLREDUCE(nP_in_bin, DEMProperty%nPrtcl_in_Bin,DEM_Opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    deallocate(nP_in_bin)

    GPrtcl_list%nlocal = nlocal
    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(nlocal_sum/= np .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_Read_Restart: "," nlocal_sum/= np_InDomain " )
    endif
    close(iUnit,IOSTAT=ierr)
  end subroutine Prtcl_Read_Restart
#undef DEM_NumRead

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_Write_Restart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
#define DEM_NumRestart 100
  subroutine Prtcl_Write_Restart(itime); implicit none 
    integer,intent(in)::itime
   
    ! locals
    type(real4)::TanDel_Un
    type(part_io_size_vec)::pvsize
    type(part_io_size_mat)::pmsize
    character(len=:),allocatable::chFile
    integer,allocatable,dimension(:)::IntVec
    integer,allocatable,dimension(:,:)::IntMat
    type(real3),allocatable,dimension(:)::real3Vec
    type(real4),allocatable,dimension(:)::real4Vec
    integer(kind=MPI_OFFSET_KIND)::disp,bgn_byte,FileSize
    integer::pid,i,j,k,nlocal,bgn_ind,nreal3,rsize,nRestart
    integer::nCntct,nCntctTotal,nCnt_ind,nTanDel,nTan_ind,nLeft
    integer::ierr,fh,color,key,Prtcl_WORLD1,Prtcl_WORLD2,ncv,prev,now,iCntct,tsize
    integer,dimension(:),allocatable::CntctVec

    if(DEM_Opt%np_InDomain<1) return
    nlocal = GPrtcl_list%nlocal
    call get_numContacts(nCntct,nTanDel)
  
    ! Calculate the bgn_ind and nlink_ind
    bgn_ind = clc_bgn_ind(nlocal)
    nCnt_ind= clc_bgn_ind(nCntct)
    nTan_ind= clc_bgn_ind(nTanDel)
    call MPI_ALLREDUCE(nCntct,nCntctTotal,1,int_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    call GPPW_CntctList%Prepare_Restart(nTan_ind)

    ! Create and empty file, Write DEM_Opt%np_InDomain,nCntctTotal, in the begining of the Restart file
    disp = 1_8
#if defined(CFDDEM) || defined(CFDACM)
    chFile = strip(DEM_Opt%RestartDir)//"RestartFor"//strip(DEM_Opt%RunName)//int2str(itime/icouple,10)
#else
    chFile = strip(DEM_Opt%RestartDir)//"RestartFor"//strip(DEM_Opt%RunName)//int2str(itime,10)
#endif
    if(nrank==0) then
      open(newunit=fh,file=chFile,status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierr)
      write(unit=fh,pos=disp,IOSTAT=ierr) DEM_Opt%np_InDomain,nCntctTotal; disp=disp+int_byte*2
#ifdef CFDDEM
      if(Is_clc_Basset) then
        i=1; j=GPrtcl_BassetSeq%HistoryStage
      else
        i=0; j=0
      endif
      write(unit=fh,pos=disp,IOSTAT=ierr) i,j
#endif
      close(fh,IOSTAT=ierr)
    endif
#ifdef CFDDEM
    disp = int_byte*4
#else
    disp = int_byte*2
#endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    ! Create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD1,ierr)
    if(color==2) return

    ! Begin to write Restart file    
    call MPI_FILE_OPEN(Prtcl_WORLD1, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    pvsize%sizes(1)   = DEM_Opt%np_InDomain
    pvsize%subsizes(1)= nlocal
    pvsize%starts(1)  = bgn_ind
    allocate(real3Vec(nlocal))
    do pid=1,nlocal
      real3Vec(pid)=GPrtcl_PosR(pid)
    enddo
    call Prtcl_dump(fh,disp, real3Vec(1:nlocal),  pvsize)
    deallocate(real3Vec)

    pmsize%sizes(1)     = 3;   pmsize%sizes(2)     = DEM_Opt%np_InDomain
    pmsize%subsizes(1)  = 3;   pmsize%subsizes(2)  = nlocal
    pmsize%starts(1)    = 0;   pmsize%starts(2)    = bgn_ind
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
    rsize=GPrtcl_list%rsize
    nreal3 = 2*(1+GPrtcl_list%tsize+GPrtcl_list%rsize)
#ifdef CFDDEM
    nreal3=nreal3+3
    if(Is_clc_Basset) nreal3=nreal3+GPrtcl_BassetSeq%nDataLen
#endif
#ifdef CFDACM
    nreal3=nreal3+4
#endif
    nreal3=nreal3-1   ! "-1" corresponds to GPrtcl_Pos
    FileSize=disp+int(nreal3*real3_byte,8)*int(DEM_Opt%np_InDomain,8)
    call MPI_FILE_OPEN(Prtcl_WORLD1, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_BARRIER(Prtcl_WORLD1,ierr)
    call my_mpi_file_set_size(fh,FileSize,ierr)
    call MPI_BARRIER(Prtcl_WORLD1,ierr)      
      
    allocate(real3Vec(DEM_NumRestart*nreal3))
    nLeft=nlocal; pid=0
    bgn_byte=disp+int(nreal3*real3_byte,8)*int(bgn_ind,8)
    DO
      nRestart=min(nLeft,DEM_NumRestart)
      k=0
      do i=1,nRestart
        pid=pid+1
        real3Vec(k+1:k+tsize)=GPrtcl_LinVel(1:tsize,pid); k=k+tsize
        real3Vec(k+1:k+tsize)=GPrtcl_LinAcc(1:tsize,pid); k=k+tsize
        real3Vec(k+1)        =GPrtcl_theta(pid);          k=k+1
        real3Vec(k+1:k+rsize)=GPrtcl_RotVel(1:rsize,pid); k=k+rsize
        real3Vec(k+1:k+rsize)=GPrtcl_RotAcc(1:rsize,pid); k=k+rsize
#ifdef CFDDEM
        real3Vec(k+1)        =GPrtcl_FpForce(pid);        k=k+1
        real3Vec(k+1)        =GPrtcl_linVelOld(pid);      k=k+1
        real3Vec(k+1)        =GPrtcl_VFluid(1,pid);       k=k+1
        if(Is_clc_Basset) then
          j=GPrtcl_BassetSeq%nDataLen
          real3Vec(k+1:k+j)  =GPrtcl_BassetData(1:j,pid); k=k+j
        endif
#endif
#ifdef CFDACM
        real3Vec(k+1)        =GPrtcl_FpForce(pid);        k=k+1
        real3Vec(k+1)        =GPrtcl_FpTorque(pid);       k=k+1
        real3Vec(k+1:k+2)    =GPrtcl_FluidIntOld(1:2,pid);k=k+2
#endif
      enddo
      call MPI_FILE_WRITE_AT(fh,bgn_byte,real3Vec,k,real3_type,MPI_STATUS_IGNORE,ierr)
      bgn_byte=bgn_byte+int(nreal3*real3_byte,8)*int(nRestart,8)
      nLeft=nLeft-nRestart
      if(nLeft==0)exit
    ENDDO
    deallocate(real3Vec)
    disp=disp+int(nreal3*real3_byte,8)*int(DEM_Opt%np_InDomain,8)
    call MPI_BARRIER(Prtcl_WORLD1,ierr)
      
    ! ncv: number of particles/walls which have overlap with this particle
    iCntct=0
    allocate(IntVec(nlocal))
    DO pid=1,nlocal
      call count_ContactList(pid,ncv)
      IntVec(pid)=ncv
      if(iCntct<ncv) iCntct=ncv
    ENDDO
    if(iCntct==0) iCntct=1
    call Prtcl_dump(fh,disp,IntVec(1:nlocal),pvsize)
    allocate( CntctVec(2*iCntct) )
    deallocate(IntVec)
    call MPI_FILE_CLOSE(fh,ierr)

    ! Begin to write Contact List file
    color = 1; key=nrank
    if(nTanDel<=0) color=2
    call MPI_COMM_SPLIT(Prtcl_WORLD1,color,key,Prtcl_WORLD2,ierr)
    call MPI_COMM_FREE( Prtcl_WORLD1, ierr)
    if(color==2) return
    call MPI_FILE_OPEN(Prtcl_WORLD2, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    pvsize%sizes(1)   = 2*nCntctTotal
    pvsize%subsizes(1)= 2*nCntct     
    pvsize%starts(1)  = 2*nCnt_ind   
    allocate(IntVec(max(2*nCntct,1)))
    iCntct=0
    DO pid=1,nlocal
      call resemble_ContactList(pid,ncv,CntctVec)
      if(ncv>0) then
        IntVec(iCntct+1:iCntct+2*ncv)=CntctVec(1:2*ncv)
        iCntct=iCntct+2*ncv
      endif
    ENDDO
    call Prtcl_dump(fh,disp,IntVec,pvsize)
    deallocate(IntVec)
    deallocate(CntctVec)
    call MPI_FILE_CLOSE(fh,ierr)

    call MPI_ALLREDUCE(nTanDel,nCntctTotal,1,int_type,MPI_SUM,Prtcl_WORLD2,ierr)
    FileSize=disp+int(real4_byte,8)*int(nCntctTotal,8)    
    call MPI_FILE_OPEN(Prtcl_WORLD2, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_BARRIER(Prtcl_WORLD2,ierr)
    call my_mpi_file_set_size(fh,FileSize,ierr)
    call MPI_BARRIER(Prtcl_WORLD2,ierr)
      
    prev=1; j=0; k=0
    disp = disp + nTan_ind * real4_byte
    allocate(real4Vec(nlocal))
    DO
      call GPPW_CntctList%GetNextTanDel_Un(TanDel_Un,prev,now);  prev=now+1
      j=j+1; k=k+1; real4Vec(j)=TanDel_Un
      if(k==nTanDel) then
        call MPI_FILE_WRITE_AT(fh,disp, real4Vec, j, real4_type, MPI_STATUS_IGNORE, ierr)
        disp = disp +j*real4_byte; exit
      endif
      if(j==nlocal) then
        call MPI_FILE_WRITE_AT(fh,disp, real4Vec, j, real4_type, MPI_STATUS_IGNORE, ierr)
        disp = disp +j*real4_byte; j=0
      endif
    ENDDO
    deallocate(real4Vec)
    call MPI_FILE_CLOSE(fh, ierr)
    call MPI_COMM_FREE(Prtcl_WORLD2, ierr)
  end subroutine Prtcl_Write_Restart
#undef DEM_NumRestart

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_Final_Visu
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Prtcl_Final_Visu(); implicit none 

    ! locals
    integer::iUnit,ierr
    character(len=:),allocatable::XdmfFile

    if(nrank/=0 .or. saveXDMFOnce) return
    xdmfFile = strip(DEM_Opt%ResultsDir)//"PartVisuFor"//strip(DEM_Opt%RunName)//".xmf"
    open(newunit=iUnit, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierr)
    if(ierr /= 0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"Prtcl_Final_Visu","Cannot open file: "//XdmfFile)
    endif
    ! XDMF/XMF Tail
    write(iUnit,'(A)') '    </Grid>'
    write(iUnit,'(A)') '</Domain>'
    write(iUnit,'(A)') '</Xdmf>'
    close(iUnit)
  end subroutine Prtcl_Final_Visu

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_Dump_Visu
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Prtcl_Dump_Visu(itime); implicit none
    integer,intent(in)::itime

    ! locals
    type(part_io_size_vec)::pvsize
    integer(kind=MPI_OFFSET_KIND)::disp
    character(len=:),allocatable::chFile
    real(RK),allocatable,dimension(:)::realVec
    type(real3),allocatable,dimension(:)::real3Vec
    integer::ierr,fh,i,color,key,Prtcl_WORLD,nlocal,bgn_ind

    ! write xdmf file first
    if(.not.saveXDMFOnce) call Write_XDMF(itime)

    ! update the bgn_ind
    nlocal = GPrtcl_list%nlocal
    bgn_ind=clc_bgn_ind(nlocal)

    ! Create and empty file
    chFile = strip(DEM_Opt%ResultsDir)//"PartVisuFor"//strip(DEM_Opt%RunName)//int2str(itime,10)
    if(nrank==0) then
      open(newunit=fh,file=chFile,status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierr)
      close(fh,IOSTAT=ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    ! create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD,ierr)
    if(color==2) return
    
    ! begin to dump
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_BARRIER(Prtcl_WORLD,ierr)
    disp = 0_MPI_OFFSET_KIND
    pvsize%sizes(1)     = DEM_Opt%np_InDomain
    pvsize%subsizes(1)  = nlocal
    pvsize%starts(1)    = bgn_ind
    if(nlocal<=0) return

    allocate(real3Vec(nlocal))
    do i=1,nlocal
      real3Vec(i)=GPrtcl_PosR(i)
    enddo
    call Prtcl_dump(fh,disp, real3Vec(1:nlocal),  pvsize)
    deallocate(real3Vec)
    if(save_ID) call Prtcl_dump(fh,disp, GPrtcl_id(1:nlocal),  pvsize)
    if(save_Diameter) then
      allocate(realVec(nlocal))
      do i=1,nlocal
        realVec(i)= 2.0_RK*GPrtcl_PosR(i)%w
      enddo   
      call Prtcl_dump(fh,disp, realVec(1:nlocal),  pvsize)
      deallocate(realVec)
    endif
    if(save_Type)       call Prtcl_dump(fh,disp, GPrtcl_pType(1:nlocal),     pvsize)
    if(save_UsrMark)    call Prtcl_dump(fh,disp, GPrtcl_UsrMark(1:nlocal),   pvsize)
    if(save_LinVel)     call Prtcl_dump(fh,disp, GPrtcl_LinVel(1,1:nlocal),  pvsize)
    if(save_LinAcc)     call Prtcl_dump(fh,disp, GPrtcl_LinAcc(1,1:nlocal),  pvsize)
    if(save_Theta)      call Prtcl_dump(fh,disp, GPrtcl_Theta(1:nlocal),     pvsize)
    if(save_RotVel)     call Prtcl_dump(fh,disp, GPrtcl_RotVel(1,1:nlocal),  pvsize)
    if(save_RotAcc)     call Prtcl_dump(fh,disp, GPrtcl_RotAcc(1,1:nlocal),  pvsize)
    if(save_CntctForce) call Prtcl_dump(fh,disp, GPrtcl_CntctForce(1:nlocal),pvsize)
    if(save_Torque)     call Prtcl_dump(fh,disp, GPrtcl_Torque(1:nlocal),    pvsize)
#ifdef CFDACM
    if(save_HighSt)  then
      block
        integer,allocatable,dimension(:)::intVec
        allocate(IntVec(nlocal))
        do i=1,nlocal
          if(GPrtcl_HighSt(i)=="N") then
            IntVec(i)=0
          else
            IntVec(i)=1
          endif
        enddo   
        call Prtcl_dump(fh,disp, IntVec(1:nlocal),  pvsize)
        deallocate(IntVec)
      end block
    endif
#endif 
    call MPI_FILE_CLOSE(fh, ierr) 
    call MPI_COMM_FREE( Prtcl_WORLD, ierr)
  end subroutine Prtcl_Dump_Visu

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
end module sp_IOAndVisu
