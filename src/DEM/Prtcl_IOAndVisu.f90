module Prtcl_IOAndVisu
  use MPI
  use Prtcl_Comm
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_Decomp_2d
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank,nproc
#endif
  implicit none
  private

  integer,parameter:: IK = 4
  integer::Prev_BackUp_itime= 123456789
  logical::saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel
  logical::save_LinAcc,save_Theta,save_RotVel,save_RotAcc,save_CntctForce,save_Torque
#ifdef CFDACM
  logical::save_HighSt
#endif

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
    procedure:: ReadFixEdCoord=>  PIO_ReadFixEdCoord
    procedure:: RestartCL     =>  PIO_RestartCL
    procedure:: Write_Restart =>  PIO_Write_Restart
    procedure:: Delete_Prev_Restart =>  PIO_Delete_Prev_Restart
    procedure,private:: Write_XDMF  =>  PIO_Write_XDMF
#ifdef CFDDEM
    procedure:: ReadInitialCoord  =>  PIO_ReadInitialCoord
    procedure:: ReadFixEdRestart  =>  PIO_ReadFixEdRestart
    procedure:: WriteFixEdRestart =>  PIO_WriteFixEdRestart
#endif
#ifdef CFDACM
    procedure:: ReadInitialCoord  =>  PIO_ReadInitialCoord
#endif
  end type Prtcl_IO_Visu
  type(Prtcl_IO_Visu),public:: DEM_IO

  ! useful interfaces
  interface Prtcl_dump
    module procedure Prtcl_dump_int_vector,  Prtcl_dump_int_matrix
    module procedure Prtcl_dump_real_vector, Prtcl_dump_real3_vector, Prtcl_dump_real3_matrix
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
    character(128)::XdmfFile
    integer::nUnitFile,ierror,indent,nflds,ifld
#ifdef CFDACM
    NAMELIST /PrtclVisuOption/ saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel,   &
                               save_LinAcc,save_Theta,save_RotVel,save_RotAcc,save_CntctForce,save_Torque,save_HighSt
#else
    NAMELIST /PrtclVisuOption/ saveXDMFOnce,save_ID,save_Diameter,save_Type,save_UsrMark,save_LinVel,   &
                               save_LinAcc,save_Theta,save_RotVel,save_RotAcc,save_CntctForce,save_Torque
#endif
  
    if(iStage==1) then
      open(newunit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=ierror)
      if(ierror/=0)call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Init_visu", "Cannot open file: "//trim(chFile))
      read(nUnitFile, nml=PrtclVisuOption)
      if(nrank==0)write(DEMLogInfo%nUnit, nml=PrtclVisuOption)
      close(nUnitFile,IOSTAT=ierror)
      return
    endif

    ! initialize the XDMF/XDF file
    if(nrank/=0) return
    write(XdmfFile,"(A)") trim(DEM_opt%ResultsDir)//"PartVisuFor"//trim(DEM_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror /= 0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Init_visu","Cannot open file:  "//trim(XdmfFile))
    write(nUnitFile,'(A)') '<?xml version="1.0" ?>'
    write(nUnitFile,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(nUnitFile,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(nUnitFile,'(A)') '<Domain>'

    ! Time series
    indent =  4
    nflds = (DEM_Opt%ilast - DEM_Opt%ifirst +1)/DEM_Opt%SaveVisu  + 1
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    indent = indent + 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'<Time TimeType="List">'
    indent = indent + 4
    write(nUnitFile,'(A,I6,A)')repeat(' ',indent)//'<DataItem Format="XML" NumberType="Int" Dimensions="',nflds,'">' 
    write(nUnitFile,'(A)',advance='no') repeat(' ',indent)
    do ifld = 1,nflds
#if defined(CFDDEM) || defined(CFDACM)
      write(nUnitFile,'(I9)',advance='no') ((ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1)/icouple
#else
      write(nUnitFile,'(I9)',advance='no')  (ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1
#endif
    enddo
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,fmt='(A)')repeat(' ',indent)//'</Time>'
    close(nUnitFile,IOSTAT=ierror)
    if( .not. saveXDMFOnce) return

    do ifld = 1,nflds
#if defined(CFDDEM) || defined(CFDACM)
      call this%Write_XDMF(((ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1)/icouple)
#else
      call this%Write_XDMF((ifld-1)*DEM_Opt%SaveVisu + DEM_Opt%ifirst-1)
#endif
    enddo

    ! XDMF/XMF Tail
    open(newunit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    write(nUnitFile,'(A)') '    </Grid>'
    write(nUnitFile,'(A)') '</Domain>'
    write(nUnitFile,'(A)') '</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)
  end subroutine PIO_Init_visu

  !**********************************************************************
  ! PIO_Delete_Prev_Restart
  !**********************************************************************
  subroutine PIO_Delete_Prev_Restart(this,itime)
    implicit none
    class(Prtcl_IO_Visu)::this
    integer:: itime

    ! locals
    character(128)::chFile

    if(nrank/=0) return
    write(chFile,"(A,I10.10)") trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName),Prev_BackUp_itime 
    call system("rm "//trim(adjustl(chFile))//" 2> /dev/null")

#ifdef CFDDEM
    if(Is_clc_FluidAcc .or. (is_clc_Basset .and. is_clc_Basset_fixEd)) then
      write(chFile,"(A,I10.10)") trim(DEM_opt%RestartDir)//"FixEdSpheresRestart",Prev_BackUp_itime
      call system("rm "//trim(adjustl(chFile))//" 2> /dev/null")
    endif
    Prev_BackUp_itime = itime/icouple
#elif CFDACM
    Prev_BackUp_itime = itime/icouple
#else
    Prev_BackUp_itime = itime
#endif
  end subroutine PIO_Delete_Prev_Restart

  !**********************************************************************
  ! PIO_ReadFixEdCoord
  !**********************************************************************
  subroutine PIO_ReadFixEdCoord(this)
    implicit none
    class(Prtcl_IO_visu)::this
 
    ! locals
    type(real3)::real3t
    character(128)::chFile
    integer,parameter::NumRead=100
    integer(kind=8)::byte_total,disp
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd,radius,diam
    type(real4),allocatable,dimension(:)::real4Vec
    integer,allocatable,dimension(:)::nP_in_bin,nP_in_bin_reduce,IntVec
    integer::i,pid,ierror,nLeft,nRead,numPrtclFix,nfix,nfix_sum,pType,nfixNew,pbyte,nUnit
    real(RK),dimension(5,NumRead)::FixEdPrtcl

    ! FixEdPrtcl_type
    pbyte=real_byte*5

    numPrtclFix = DEM_opt%numPrtclFix
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd
#ifdef ChanBraunJFM2011
    ySt=-yEd
#endif
#ifdef JiaYan
    if(DEM_decomp%ProcNgh(1) == MPI_PROC_NULL)zEd=zEd+2.0E-3
    if(DEM_decomp%ProcNgh(2) == MPI_PROC_NULL)zSt=zSt-2.0E-3    
    if(DEM_decomp%ProcNgh(3) == MPI_PROC_NULL)xEd=xEd+2.0E-3
    if(DEM_decomp%ProcNgh(4) == MPI_PROC_NULL)xSt=xSt-2.0E-3
#endif

    ! The data storage sequence in file "FixedSpheresCoord.dat" is as follow:
    ! Particle1: Position(real3 type), Diameter(real type), Prtcl_Type(real type)
    ! Particle2: Position(real3 type), Diameter(real type), Prtcl_Type(real type) ..
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"FixedSpheresCoord.dat"
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',position='append',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdCoord","Cannot open file: "//trim(chFile))
    inquire(unit=nUnit,Pos=disp); disp=disp-1_8
    rewind(unit=nUnit,IOSTAT=ierror)
    byte_total= int(pbyte,8)*int(numPrtclFix,8)
    if(disp/=byte_total .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdCoord","file byte wrong")

    nfix=0; nLeft=numPrtclFix; pid=0; disp= 1_8
    allocate(nP_in_bin(DEM_opt%numPrtcl_Type));        nP_in_bin=0
    allocate(nP_in_bin_reduce(DEM_opt%numPrtcl_Type)); nP_in_bin_reduce=0
    DO 
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp,IOSTAT=ierror)FixEdPrtcl(:,1:nRead)
      disp=disp+int(pbyte,8)*int(nRead,8)
      do i=1,nRead
        pid=pid+1
        real3t%x=FixEdPrtcl(1,i)
        real3t%y=FixEdPrtcl(2,i)
        real3t%z=FixEdPrtcl(3,i)
        diam    =FixEdPrtcl(4,i)
        pType   =int(FixEdPrtcl(5,i)+0.2)
        radius= DEMProperty%Prtcl_PureProp(pType)%Radius
        if( abs(two*radius/diam -one)>1.0E-6 .and. nrank==0 )then
          call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdCoord","Diameter not coordinate")
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
         
        GPFix_id(nfix)     = pid  + DEM_opt%numPrtcl         ! NOTE HERE
        GPFix_pType(nfix)  = pType
        GPFix_PosR(nfix)   = real3t
        GPFix_PosR(nfix)%w = radius
        nP_in_bin(pType)   = nP_in_bin(pType)+1
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    close(nUnit,IOSTAT=ierror)

    call MPI_ALLREDUCE(nP_in_bin, nP_in_bin_reduce,DEM_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
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
    call MPI_REDUCE(nfix,nfix_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nfix_sum/= numPrtclFix .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdCoord"," nfix_sum/= numPrtclFix " )
    endif
  end subroutine PIO_ReadFixEdCoord

  !**********************************************************************
  ! PIO_RestartCL
  !**********************************************************************
  subroutine PIO_RestartCL(this)
    implicit none
    class(Prtcl_IO_Visu)::this

    ! locals
    character(24)::ch
    character(128)::chFile
    integer,parameter::NumRead=2000
    type(real4)::real4t,TanDel_Un(20)
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd
    type(real3),allocatable,dimension(:)::PosVec
    integer,allocatable,dimension(:)::ncvVec,iCLvec
    integer::itime,nUnit,ierror,nlocal,np,nreal3,pbyte
    integer::CntctVec(20),intvec(2),i,j,ncv,nCntct,nCntctTotal,nLeft,nRead
    integer(kind=8)::disp,disp_pos,disp_ncv,disp_CL,disp_TanStart,disp_Tan

    itime = DEM_Opt%ifirst - 1
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd
#if defined(CFDDEM) || defined(CFDACM)
    write(ch,'(I10.10)')itime/icouple
#else
    write(ch,'(I10.10)')itime
#endif

    nreal3 = 2*(1+GPrtcl_list%tsize+GPrtcl_list%rsize)
#ifdef CFDDEM
    nreal3=nreal3 +3      !(GPrtcl_FpForce, GPrtcl_linVelOld, GPrtcl_VFluid(1,:))
    if(Is_clc_Basset) nreal3= nreal3+ GPrtcl_BassetSeq%nDataLen
    disp=disp+2*int_byte; !skip Is_clc_Basset and HistoryStage
#endif
#ifdef CFDACM
    nreal3=nreal3+2
#endif
    pbyte=int_byte*3 + nreal3*real3_byte  ! corresponding to the subroutine 'PIO_Write_Restart'
    
    ! Begin to read Restart_Contact_List
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_RestartCL","Cannot open file: "//trim(chFile))

    np= DEM_Opt%np_InDomain; disp=1_8 + int_byte  ! firstly skip the np_InDomain
    read(nUnit,pos=disp,IOSTAT=ierror)nCntctTotal; disp=disp+int_byte;

    allocate(PosVec(NumRead),ncvVec(NumRead))
    nlocal=0; nLeft=np
    disp_pos = disp
    disp_ncv = disp + pbyte*np
    disp_CL  = disp_ncv+ int_byte*np
    disp_TanStart= disp_CL + 2*nCntctTotal*int_byte
    DO
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp_pos,IOSTAT=ierror)PosVec(1:nRead)
      read(nUnit,pos=disp_ncv,IOSTAT=ierror)ncvVec(1:nRead)
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
            read(nUnit,pos=disp_CL,IOSTAT=ierror)intvec(1:2); disp_CL=disp_CL+int_byte*2
            disp_Tan = disp_TanStart + real4_byte*(intvec(2)-1)
            read(nUnit,pos=disp_Tan,IOSTAT=ierror)real4t
            CntctVec(j) = intvec(1)
            TanDel_Un(j)= real4t
          enddo
          if(ncv>0) call GPPW_CntctList%Add_RestartCntctlink(nlocal,ncv,CntctVec,TanDel_Un)
        endif
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    deallocate(PosVec,ncvVec)
    close(nUnit,IOSTAT=ierror)
  end subroutine PIO_RestartCL

#if defined(CFDDEM) || defined(CFDACM)
  !**********************************************************************
  ! PIO_ReadInitialCoord
  !**********************************************************************
  subroutine PIO_ReadInitialCoord(this)
    implicit none
    class(Prtcl_IO_visu)::this
 
    ! locals
    type(real3)::real3t
    character(128)::chFile
    integer,parameter::NumRead=100
    integer(kind=8)::byte_total,disp
    integer,allocatable,dimension(:)::nP_in_bin
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd,radius,diam
    integer::i,pid,ierror,nLeft,nRead,numPrtcl,nlocal,nlocal_sum,pType,nUnit,pbyte
    real(RK),dimension(5,NumRead)::InitPrtclIn

    ! InitPrtcl_type
    pbyte=real_byte*5

    numPrtcl = DEM_opt%numPrtcl
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd

    ! The data storage sequence in file "SpheresCoord.dat" is as follow:
    ! Particle1: Position(real3 type), Diameter(real type), Prtcl_Type(real type)
    ! Particle2: Position(real3 type), Diameter(real type), Prtcl_Type(real type) ...
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"SpheresCoord.dat"
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',position='append',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadInitialCoord","Cannot open file: "//trim(chFile))
    inquire(unit=nUnit,Pos=disp); disp=disp-1_8; 
    rewind(unit=nUnit,IOSTAT=ierror)
    byte_total=int(pbyte,8)*int(numPrtcl,8)
    if(disp/=byte_total .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadInitialCoord","file byte wrong")

    allocate(nP_in_bin(DEM_opt%numPrtcl_Type))
    nlocal=0; nLeft=numPrtcl; pid=0; disp= 1_8; nP_in_bin=0
    DO
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp,IOSTAT=ierror)InitPrtclIn(:,1:nRead)
      disp=disp+int(pbyte,8)*int(nRead,8)
      do i=1,nRead
        pid=pid+1
        real3t%x =InitPrtclIn(1,i)
        real3t%y =InitPrtclIn(2,i)
        real3t%z =InitPrtclIn(3,i)
        diam     =InitPrtclIn(4,i)
        pType=int(InitPrtclIn(5,i)+0.2)
        radius= DEMProperty%Prtcl_PureProp(pType)%Radius
        if(abs(two*radius/diam -one)>1.0E-6 .and. nrank==0 )then
          call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadInitialCoord","Diameter not coordinate")
        endif
        if(real3t%x< xSt .or. real3t%y< ySt .or. real3t%z< zSt .or. &
           real3t%x>=xEd .or. real3t%y>=yEd .or. real3t%z>=zEd) cycle
        if(nlocal>=GPrtcl_list%mlocal)  call GPrtcl_list%ReallocatePrtclVar(nlocal)
        nlocal=nlocal+1
       
        GPrtcl_id(nlocal)     = pid
        GPrtcl_pType(nlocal)  = pType
        GPrtcl_PosR(nlocal)   = real3t
        GPrtcl_PosR(nlocal)%w = radius
        nP_in_bin(pType)= nP_in_bin(pType)+1      
      enddo
      nLeft=nLeft-nRead
      if(nLeft==0)exit
    ENDDO
    close(nUnit,IOSTAT=ierror)
    GPrtcl_list%nlocal = nlocal

    call MPI_ALLREDUCE(nP_in_bin, DEMProperty%nPrtcl_in_Bin,DEM_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    deallocate(nP_in_bin)

    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nlocal_sum/= numPrtcl .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadInitialCoord","nlocal_sum/= numPrtcl " )
    endif
  end subroutine PIO_ReadInitialCoord
#endif
#ifdef CFDDEM
  !**********************************************************************
  ! PIO_WriteFixEdRestart
  !**********************************************************************
  subroutine PIO_WriteFixEdRestart(this,itime)
    implicit none 
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime

    ! locals
    character(128)::chFile
    logical::IsWriteBasset
    type(part_io_size_vec)::pvsize
    type(part_io_size_mat)::pmsize
    integer,parameter::NumWrite=100
    real(RK),dimension(5,NumWrite)::FixEdPrtclOut
    integer(kind=MPI_OFFSET_KIND)::disp,disp_pos,FileSize
    integer::bgn_ind,color,key,ierror,PrtclFix_WORLD,fh,nLeft,nWrite,prank,i,pbyte,pid,intVec(2)

    if(DEM_opt%numPrtclFix<1) return
    IsWriteBasset=.false.
    if(is_clc_Basset .and. is_clc_Basset_fixEd) IsWriteBasset=.true.
    if((.not.Is_clc_FluidAcc) .and. (.not.IsWriteBasset)) return

    ! Create the Prtcl_GROUP
    bgn_ind= clc_bgn_ind(GPrtcl_list%mlocalFix)
    color = 1; key=nrank
    if(GPrtcl_list%mlocalFix<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,PrtclFix_WORLD,ierror)
    if(color==2) return

    write(chFile,"(A,I10.10)") trim(DEM_opt%RestartDir)//"FixEdSpheresRestart",itime/icouple
    call MPI_FILE_OPEN(PrtclFix_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_8,ierror)  ! guarantee overwriting
    call MPI_BARRIER(PrtclFix_WORLD,ierror)
    disp = 0_MPI_OFFSET_KIND

    call MPI_COMM_RANK(PrtclFix_WORLD, prank, ierror)
    if(IsWriteBasset) then
      intVec=(/1,GPrtcl_BassetSeq%HistStageFix/)
    else
      intVec=(/0,0/)
    endif
    if(prank==0) then
      call MPI_FILE_WRITE_AT(fh,disp,intVec,2,int_type,MPI_STATUS_IGNORE, ierror)
    endif
    disp=disp+int_byte*2
    call MPI_BARRIER(PrtclFix_WORLD,ierror)

    ! PosD, id, pType
    pbyte=real_byte*5
    
    call MPI_FILE_GET_SIZE(fh,FileSize,ierror)
    FileSize=FileSize+int(pbyte,8)*int(DEM_opt%numPrtclFix,8)
    call MPI_FILE_PREALLOCATE(fh,FileSize,ierror)
    call MPI_BARRIER(PrtclFix_WORLD,ierror)
      
    nLeft=GPrtcl_list%mlocalFix;  pid=0
    disp_pos=disp+int(pbyte,8)*int(bgn_ind,8)
    DO 
      nWrite=min(nLeft,NumWrite)
      do i=1,nWrite
        pid=pid+1
        FixEdPrtclOut(1,i)=GPFix_PosR(pid)%x
        FixEdPrtclOut(2,i)=GPFix_PosR(pid)%y
        FixEdPrtclOut(3,i)=GPFix_PosR(pid)%z
        FixEdPrtclOut(4,i)=real(GPFix_id(pid),RK)
        FixEdPrtclOut(5,i)=real(GPFix_pType(pid),RK)
      enddo
      call MPI_FILE_WRITE_AT(fh,disp_pos,FixEdPrtclOut,5*nWrite,real_type,MPI_STATUS_IGNORE, ierror)
      disp_pos=disp_pos+int(pbyte,8)*int(nWrite,8)
      nLeft=nLeft-nWrite
      if(nLeft==0)exit
    ENDDO
    disp=disp+int(pbyte,8)*int(DEM_opt%numPrtclFix,8)
    call MPI_BARRIER(PrtclFix_WORLD,ierror)

    ! Write GPFix_Vfluid and GPFix_BassetData
    if(Is_clc_FluidAcc) then
      pvsize%sizes(1)   = DEM_opt%numPrtclFix
      pvsize%subsizes(1)= GPrtcl_list%mlocalFix
      pvsize%starts(1)  = bgn_ind
      call Prtcl_dump(fh, disp, GPFix_Vfluid(1,:), pvsize)
    endif
    if(IsWriteBasset) then
      pmsize%sizes(1)   = GPrtcl_BassetSeq%nDataLen;  pmsize%sizes(2)   = DEM_opt%numPrtclFix
      pmsize%subsizes(1)= GPrtcl_BassetSeq%nDataLen;  pmsize%subsizes(2)= GPrtcl_list%mlocalFix
      pmsize%starts(1)  = 0 ;                         pmsize%starts(2)  = bgn_ind
      call Prtcl_dump(fh, disp, GPFix_BassetData, pmsize)
    endif
    call MPI_FILE_CLOSE(fh, ierror) 
    call MPI_COMM_FREE(PrtclFix_WORLD,ierror)
  end subroutine PIO_WriteFixEdRestart

  !**********************************************************************
  ! PIO_ReadFixEdRestart
  !**********************************************************************
  subroutine PIO_ReadFixEdRestart(this)
    implicit none
    class(Prtcl_IO_Visu)::this

    ! locals
    type(real3)::real3t
    logical::IsReadBasset
    character(128)::chFile
    integer,parameter::NumRead=500
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd
    integer(kind=8)::byte_total,disp,disp_Fluid
    type(real3),allocatable,dimension(:)::Real3Vec
    integer,dimension(:),allocatable::nP_in_bin,nP_in_bin_reduce,idVec,IntVec
    integer::i,pbyte,pid,iPos,pType,nUnit,ierror,numPrtclFix,nfix,nfixNew,nRead,nLeft,nfix_sum,int_t(2)
    real(RK)::FixEdPrtclIn(5,NumRead),FixEdPrtclOne(5)

    IsReadBasset=.false.
    if(is_clc_Basset .and. is_clc_Basset_fixEd) IsReadBasset=.true.
    if((.not.Is_clc_FluidAcc) .and. (.not.IsReadBasset)) then
      call this%ReadFixEdCoord()
      return
    endif

    numPrtclFix = DEM_opt%numPrtclFix
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd

    ! PosD, id, pType
    pbyte=real_byte*5

    write(chFile,'(I10.10)')(DEM_Opt%ifirst - 1)/icouple
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"FixEdSpheresRestart"//trim(adjustl(chFile))
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',position='append',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdRestart","Cannot open file:"//trim(chFile))

    byte_total= int_byte*2+pbyte*numPrtclFix ! 'int_byte*2' corresponds to HistStageFix in  PIO_WriteFixEdRestart
    if(is_clc_FluidAcc) byte_total= byte_total+ real3_byte*numPrtclFix
    if(IsReadBasset) byte_total= byte_total+ real3_byte*numPrtclFix*GPrtcl_BassetSeq%nDataLen
    inquire(unit=nUnit,Pos=disp); disp=disp-1_8
    rewind(unit=nUnit,IOSTAT=ierror)
    if(disp/=byte_total .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdCoordRestart","file byte wrong")

    disp=1_8
    read(nUnit,pos=disp,IOSTAT=ierror)int_t(1:2); disp=disp+int_byte*2
    if(IsReadBasset) then
      if(int_t(1)/= 1) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdRestart"," Is_clc_Basset_fixEd Wrong 1" )
      GPrtcl_BassetSeq%HistStageFix= int_t(2)
    else
      if(int_t(1)/= 0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdRestart"," Is_clc_Basset_fixEd Wrong 2" )
    endif   

    nfix=0; nLeft=numPrtclFix; pid=0
    allocate(idVec(GPrtcl_list%mlocalFix))           
    allocate(nP_in_bin(DEM_opt%numPrtcl_Type));        nP_in_bin=0
    allocate(nP_in_bin_reduce(DEM_opt%numPrtcl_Type)); nP_in_bin_reduce=0    
    DO 
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp,IOSTAT=ierror)FixEdPrtclIn(:,1:nRead)
      disp=disp+int(pbyte,8)*int(nRead,8)
      do i=1,nRead
        pid=pid+1
        real3t%x=FixEdPrtclIn(1,i)
        real3t%y=FixEdPrtclIn(2,i)
        real3t%z=FixEdPrtclIn(3,i)
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
    call MPI_ALLREDUCE(nP_in_bin, nP_in_bin_reduce,DEM_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    DEMProperty%nPrtcl_in_Bin= DEMProperty%nPrtcl_in_Bin+ nP_in_bin_reduce
    deallocate(nP_in_bin,nP_in_bin_reduce)

    GPrtcl_list%mlocalFix = nfix
    call MPI_REDUCE(nfix,nfix_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nfix_sum/= numPrtclFix .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_ReadFixEdRestart"," nfix_sum/= numPrtclFix " )
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
        allocate(GPFix_BassetData(GPrtcl_BassetSeq%nDataLen, nfix), Stat= ierror)
        if(ierror/=0) call DEMLogInfo%checkForError(ErrT_Abort,"PIO_ReadFixEdRestart","allocate wrong2")
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
      read(nUnit,pos=disp,IOSTAT=ierror)FixEdPrtclOne
      GPFix_id(pid)   = int(FixEdPrtclOne(4)+0.2)
      pType = int(FixEdPrtclOne(5)+0.2)
      GPFix_pType(pid)= pType
      GPFix_PosR(pid)%x= FixEdPrtclOne(1)
      GPFix_PosR(pid)%y= FixEdPrtclOne(2)
      GPFix_PosR(pid)%z= FixEdPrtclOne(3)
      GPFix_PosR(pid)%w= DEMProperty%Prtcl_PureProp(pType)%Radius
      if(Is_clc_FluidAcc) then
        disp_Fluid= 1_8+int_byte*2+int(numPrtclFix,8)*int(pbyte,8) +int(iPos,8)*int(real3_byte,8)
        read(nUnit,pos=disp_Fluid,IOSTAT=ierror)real3t
        GPFix_Vfluid(1,pid)=real3t
      endif
      if(IsReadBasset) then
        disp_Fluid= 1_8+int_byte*2+int(numPrtclFix,8)*int(pbyte+real3_byte,8) +int(iPos,8)*int(real3_byte*GPrtcl_BassetSeq%nDataLen,8)
        read(nUnit,pos=disp_Fluid,IOSTAT=ierror)Real3Vec(1:GPrtcl_BassetSeq%nDataLen)
        GPFix_BassetData(:,pid)= Real3Vec
      endif
    enddo
    deallocate(idVec,Real3Vec)
    close(nUnit,IOSTAT=ierror)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  end subroutine PIO_ReadFixEdRestart
#endif

  !**********************************************************************
  ! PIO_Read_Restart
  !**********************************************************************
  subroutine PIO_Read_Restart(this)
    implicit none
    class(Prtcl_IO_Visu)::this

    ! locals
    character(24)::ch
    character(128)::chFile
    integer,parameter::NumRead=2000
    real(RK)::xSt,xEd,ySt,yEd,zSt,zEd
    integer,allocatable,dimension(:):: nP_in_bin
    type(real3),allocatable,dimension(:)::real3Vec,PosVec
    integer::nUnit,tsize,rsize,nlocal_sum,nreal3,itype
    integer::itime,ierror,nlocal,i,k,np,nLeft,nRead,int_t(3)
    integer(kind=MPI_OFFSET_KIND)::disp,disp_pos,disp_int,disp_real3
    
    itime = DEM_Opt%ifirst - 1
    xSt=DEM_decomp%xSt; xEd=DEM_decomp%xEd
    ySt=DEM_decomp%ySt; yEd=DEM_decomp%yEd
    zSt=DEM_decomp%zSt; zEd=DEM_decomp%zEd
#if defined(CFDDEM) || defined(CFDACM)
    write(ch,'(I10.10)')itime/icouple
#else
    write(ch,'(I10.10)')itime
#endif

    ! Begin to read Restart_file
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    open(newunit=nUnit,file=trim(chFile),status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart","Cannot open file: "//trim(chFile))
    disp =1_8; read(nUnit,pos=disp,IOSTAT=ierror)np; disp=disp+int_byte
    if(np>DEM_Opt%numPrtcl .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart: "," np_InDomain > numPrtcl " )
    endif
    DEM_Opt%np_InDomain = np
    disp=disp+int_byte   ! skip nCntctTotal
#ifdef CFDDEM
    read(nUnit,pos=disp,IOSTAT=ierror)int_t(1:2); disp=disp+2*int_byte;
    if(Is_clc_Basset) then
      if(int_t(1)/= 1) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart"," Is_clc_Basset Wrong 1" )
      GPrtcl_BassetSeq%HistoryStage= int_t(2)
    else
      if(int_t(1)/= 0) call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart"," Is_clc_Basset Wrong 2" )
    endif
#endif

    tsize=GPrtcl_list%tsize; rsize=GPrtcl_list%rsize;
    nreal3 = 2*(1+tsize+rsize)
#ifdef CFDDEM
    nreal3=nreal3+3
    if(Is_clc_Basset) nreal3=nreal3+GPrtcl_BassetSeq%nDataLen
#endif
#ifdef CFDACM
    nreal3=nreal3+2
#endif
    nreal3=nreal3-1   ! "-1" corresponds to GPrtcl_Pos
    allocate(real3Vec(nreal3),PosVec(NumRead))
    allocate(nP_in_bin(DEM_opt%numPrtcl_Type)); nP_in_bin=0

    nlocal=0; nLeft=np
    disp_pos  = disp
    disp_int  = disp_pos+ real3_byte*np
    disp_real3= disp_int+ int_byte*np*3
    DO
      nRead=min(nLeft,NumRead)
      read(nUnit,pos=disp_pos,IOSTAT=ierror)PosVec(1:nRead)
      disp_pos=disp_pos+int(real3_byte,8)*int(nRead,8)
      do i=1,nRead
        if(PosVec(i)%x>=xSt .and. PosVec(i)%x< xEd .and. PosVec(i)%y>=ySt .and. &
           PosVec(i)%y< yEd .and. PosVec(i)%z>=zSt .and. PosVec(i)%z<zEd) then
          if(nlocal>=GPrtcl_list%mlocal)  call GPrtcl_list%ReallocatePrtclVar(nlocal)
          nlocal=nlocal+1

          read(nUnit,pos=disp_int,IOSTAT=ierror)int_t(1:3)
          GPrtcl_id(nlocal)=int_t(1)      ! id
          itype=int_t(2)
          GPrtcl_pType(nlocal)=itype      ! pType
          nP_in_bin(itype)= nP_in_bin(itype)+1
          GPrtcl_UsrMark(nlocal)=int_t(3) ! Usr_Mark

          GPrtcl_PosR(nlocal)= PosVec(i)  ! PosR
          GPrtcl_PosR(nlocal)%w=DEMProperty%Prtcl_PureProp(itype)%Radius
          k=0;
          read(nUnit,pos=disp_real3,IOSTAT=ierror)real3Vec(1:nreal3)
          GPrtcl_LinVel(1:tsize,nlocal)  =real3Vec(k+1:k+tsize); k=k+tsize ! LinVec
          GPrtcl_LinAcc(1:tsize,nlocal)  =real3Vec(k+1:k+tsize); k=k+tsize ! LinAcc
          GPrtcl_theta(nlocal)           =real3Vec(k+1);         k=k+1     ! Theta
          GPrtcl_RotVel(1:rsize,nlocal)  =real3Vec(k+1:k+rsize); k=k+rsize ! RotVel
          GPrtcl_RotAcc(1:rsize,nlocal)  =real3Vec(k+1:k+rsize); k=k+rsize ! RotAcc
#ifdef CFDACM        
          GPrtcl_FluidIntOld(1:2,nlocal)= real3Vec(k+1:k+2);     k=k+2     ! FluidIntOld
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
    call MPI_ALLREDUCE(nP_in_bin, DEMProperty%nPrtcl_in_Bin,DEM_opt%numPrtcl_Type,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    deallocate(nP_in_bin)

    GPrtcl_list%nlocal = nlocal
    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nlocal_sum/= np .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Read_Restart: "," nlocal_sum/= np_InDomain " )
    endif
    close(nUnit,IOSTAT=ierror)
  end subroutine PIO_Read_Restart

  !**********************************************************************
  ! PIO_Write_Restart
  !**********************************************************************
  subroutine PIO_Write_Restart(this,itime)
    implicit none 
    class(Prtcl_IO_Visu)::this
    integer,intent(in)::itime
   
    ! locals 
    character(24)::ch
    character(128)::chFile
    type(real4)::TanDel_Un
    type(part_io_size_vec)::pvsize
    type(part_io_size_mat)::pmsize
    integer,parameter::NumRestart=100
    integer,allocatable,dimension(:)::IntVec
    integer,allocatable,dimension(:,:)::IntMat
    type(real3),allocatable,dimension(:)::real3Vec
    type(real4),allocatable,dimension(:)::real4Vec
    integer(kind=MPI_OFFSET_KIND)::disp,bgn_byte,FileSize
    integer::pid,i,j,k,nlocal,bgn_ind,prank,nreal3,CntctVec(40),rsize
    integer::nCntct,nCntctTotal,nCnt_ind,nTanDel,nTan_ind,nRestart,nLeft
    integer::ierror,fh,color,key,Prtcl_WORLD1,Prtcl_WORLD2,ncv,prev,now,iCntct,tsize

    if(DEM_Opt%np_InDomain<1) return
    nlocal    = GPrtcl_list%nlocal
    call GPPW_CntctList%Get_numCntcts(nCntct,nTanDel)    
#if defined(CFDDEM) || defined(CFDACM)
    write(ch,'(I10.10)')itime/icouple
#else
    write(ch,'(I10.10)')itime
#endif
  
    ! Calculate the bgn_ind and nlink_ind
    bgn_ind = clc_bgn_ind(nlocal)
    nCnt_ind= clc_bgn_ind(nCntct)
    nTan_ind= clc_bgn_ind(nTanDel)
    call MPI_ALLREDUCE(nCntct,nCntctTotal,1,int_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    call GPPW_CntctList%Prepare_Restart(nTan_ind)

    ! Create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal<=0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD1,ierror)
    if(color==2) return

    ! Begin to write Restart file
    write(chFile,"(A)") trim(DEM_opt%RestartDir)//"RestartFor"//trim(DEM_opt%RunName)//trim(adjustl(ch))
    call MPI_FILE_OPEN(Prtcl_WORLD1, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_8,ierror)  ! Guarantee overwriting
    call MPI_BARRIER(Prtcl_WORLD1,ierror)
    disp = 0_MPI_OFFSET_KIND

    ! Write DEM_Opt%np_InDomain,nCntctTotal, in the begining of the Restart file
    call MPI_COMM_RANK(Prtcl_WORLD1, prank, ierror)
    if(prank==0) then
      call MPI_FILE_WRITE_AT(fh,disp,[DEM_Opt%np_InDomain,nCntctTotal],2,int_type,MPI_STATUS_IGNORE,ierror)       
    endif
    disp = disp + int_byte*2
#ifdef CFDDEM
    if(Is_clc_Basset) then
      i=1; j=GPrtcl_BassetSeq%HistoryStage
    else
      i=0; j=0
    endif
    if(prank==0) then
      call MPI_FILE_WRITE_AT(fh,disp,[i,j],2,int_type,MPI_STATUS_IGNORE,ierror)       
    endif
    disp = disp + int_byte*2
#endif
    call MPI_BARRIER(Prtcl_WORLD1,ierror) 
      
    ! Begin to write
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
    call MPI_FILE_CLOSE(fh,ierror)

    tsize=GPrtcl_list%tsize
    rsize=GPrtcl_list%rsize
    nreal3 = 2*(1+GPrtcl_list%tsize+GPrtcl_list%rsize)
#ifdef CFDDEM
    nreal3=nreal3+3
    if(Is_clc_Basset) nreal3=nreal3+GPrtcl_BassetSeq%nDataLen
#endif
#ifdef CFDACM
    nreal3=nreal3+2
#endif
    nreal3=nreal3-1   ! "-1" corresponds to GPrtcl_Pos
    call MPI_FILE_OPEN(Prtcl_WORLD1, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_GET_SIZE(fh,FileSize,ierror)
    FileSize=FileSize+int(nreal3*real3_byte,8)*int(DEM_Opt%np_InDomain,8)
    call MPI_FILE_PREALLOCATE(fh,FileSize,ierror)
    call MPI_BARRIER(Prtcl_WORLD1,ierror)      
      
    allocate(real3Vec(NumRestart*nreal3))
    nLeft=nlocal; pid=0
    bgn_byte=disp+int(nreal3*real3_byte,8)*int(bgn_ind,8)
    DO
      nRestart=min(nLeft,NumRestart)
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
          real3Vec(k+1:k+j) =GPrtcl_BassetData(1:j,pid);  k=k+j
        endif
#endif
#ifdef CFDACM
        real3Vec(k+1:k+2)=GPrtcl_FluidIntOld(1:2,pid);    k=k+2
#endif
      enddo
      call MPI_FILE_WRITE_AT(fh,bgn_byte,real3Vec,k,real3_type,MPI_STATUS_IGNORE,ierror)
      bgn_byte=bgn_byte+int(nreal3*real3_byte,8)*int(nRestart,8)
      nLeft=nLeft-nRestart
      if(nLeft==0)exit
    ENDDO
    deallocate(real3Vec)
    disp=disp+int(nreal3*real3_byte,8)*int(DEM_Opt%np_InDomain,8)
    call MPI_BARRIER(Prtcl_WORLD1,ierror)
      
    ! ncv: number of particles/walls which have overlap with this particle 
    allocate(IntVec(nlocal)) 
    DO pid=1,nlocal
      call GPPW_CntctList%Count_Cntctlink(pid,ncv)
      IntVec(pid)=ncv
    ENDDO
    call Prtcl_dump(fh,disp,IntVec(1:nlocal),pvsize)
    deallocate(IntVec)
    call MPI_FILE_CLOSE(fh,ierror)

    ! Begin to write Contact List file
    color = 1; key=nrank
    if(nTanDel<=0) color=2
    call MPI_COMM_SPLIT(Prtcl_WORLD1,color,key,Prtcl_WORLD2,ierror)
    call MPI_COMM_FREE( Prtcl_WORLD1, ierror)
    if(color==2) return
    call MPI_FILE_OPEN(Prtcl_WORLD2, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    pvsize%sizes(1)   = 2*nCntctTotal
    pvsize%subsizes(1)= 2*nCntct     
    pvsize%starts(1)  = 2*nCnt_ind   
    allocate(IntVec(2*nCntct))
    iCntct=0
    DO pid=1,nlocal
      call GPPW_CntctList%Gather_Cntctlink_Restart(pid,ncv,CntctVec)
      if(ncv>0) then
        IntVec(iCntct+1:iCntct+2*ncv)=CntctVec(1:2*ncv)
        iCntct=iCntct+2*ncv
      endif
    ENDDO
    call Prtcl_dump(fh,disp,IntVec,pvsize)
    deallocate(IntVec)
    call MPI_FILE_CLOSE(fh,ierror)

    call MPI_ALLREDUCE(nTanDel,nCntctTotal,1,int_type,MPI_SUM,Prtcl_WORLD2,ierror)
    call MPI_FILE_OPEN(Prtcl_WORLD2, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_GET_SIZE(fh,FileSize,ierror)
    FileSize=FileSize+int(real4_byte,8)*int(nCntctTotal,8)
    call MPI_COMM_RANK(Prtcl_WORLD2, prank, ierror)
    call MPI_FILE_PREALLOCATE(fh,FileSize,ierror)
    call MPI_BARRIER(Prtcl_WORLD2,ierror)
      
    prev=1; j=0; k=0
    disp = disp + nTan_ind * real4_byte
    allocate(real4Vec(nlocal))
    DO
      call GPPW_CntctList%GetNextTanDel_Un(TanDel_Un,prev,now);  prev=now+1
      j=j+1; k=k+1; real4Vec(j)=TanDel_Un
      if(k==nTanDel) then
        call MPI_FILE_WRITE_AT(fh,disp, real4Vec, j, real4_type, MPI_STATUS_IGNORE, ierror)
        disp = disp +j*real4_byte; exit
      endif
      if(j==nlocal) then
        call MPI_FILE_WRITE_AT(fh,disp, real4Vec, j, real4_type, MPI_STATUS_IGNORE, ierror)
        disp = disp +j*real4_byte; j=0
      endif
    ENDDO
    deallocate(real4Vec)
    call MPI_FILE_CLOSE(fh, ierror)
    call MPI_COMM_FREE(Prtcl_WORLD2, ierror)
  end subroutine PIO_Write_Restart

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

    ! locals
    character(128)::chFile
    integer(kind=MPI_OFFSET_KIND)::disp
    integer:: indent,nUnitFile,ierror,np,dims,iprec

    if(nrank/=0) return 
    np=DEM_Opt%np_InDomain
    write(chFile,"(A)") trim(DEM_opt%ResultsDir)//"PartVisuFor"//trim(DEM_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=chFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    if(ierror/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Write_XDMF","Cannot open file: "//trim(chFile))
    endif

    indent = 8; disp = 0_MPI_OFFSET_KIND
    write(chFile,"(A)") "PartVisuFor"//trim(DEM_opt%RunName)
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
    write(nUnitFile,'(A,I10.10)')repeat(' ',indent)//trim(chFile),itime
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</DataItem>'
    indent = indent - 4
    write(nUnitFile,'(A)')repeat(' ',indent)//'</Geometry>'

    IF(save_ID) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"ID","Scalar","Int",disp)
    ENDIF
    IF(save_Diameter) THEN
      dims=1; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Diameter","Scalar","Float",disp)
    ENDIF
    IF(save_Type) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Type","Scalar","Int",disp)
    ENDIF
    IF(save_UsrMark) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"UsrMark","Scalar","Int",disp)
    ENDIF
    IF(save_LinVel) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"LinVel","Vector","Float",disp)
    ENDIF
    IF(save_LinAcc) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"LinAcc","Vector","Float",disp)
    ENDIF
    IF(save_Theta) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Theta","Vector","Float",disp)
    ENDIF
    IF(save_RotVel) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"RotVel","Vector","Float",disp)
    ENDIF
    IF(save_RotAcc) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"RotAcc","Vector","Float",disp)
    ENDIF
    IF(save_CntctForce) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"CntctForce","Vector","Float",disp)
    ENDIF
    IF(save_Torque) THEN
      dims=3; iprec=RK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"Torque","Vector","Float",disp)
    ENDIF
#ifdef CFDACM
    IF(save_HighSt) THEN
      dims=1; iprec=IK
      call Write_XDMF_One(nUnitFile,dims,iprec,np,itime,chFile,"IsHighSt","Scalar","Int",disp)
    ENDIF
#endif    
    write(nUnitFile,'(A)')'        </Grid>'
    close(nUnitFile,IOSTAT=ierror)
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
  ! PIO_Final_visu
  !**********************************************************************
  subroutine PIO_Final_visu(this)
    implicit none 
    class(Prtcl_IO_Visu)::this

    ! locals
    integer::nUnitFile,ierror
    character(128)::XdmfFile

    if(nrank/=0 .or. saveXDMFOnce) return
    write(xdmfFile,"(A)") trim(DEM_opt%ResultsDir)//"PartVisuFor"//trim(DEM_opt%RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierror)
    if(ierror /= 0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"PIO_Final_visu","Cannot open file: "//trim(XdmfFile))
    endif
    ! XDMF/XMF Tail
    write(nUnitFile,'(A)') '    </Grid>'
    write(nUnitFile,'(A)') '</Domain>'
    write(nUnitFile,'(A)') '</Xdmf>'
    close(nUnitFile)

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
    integer::ierror,fh,i,color,key,Prtcl_WORLD,nlocal,bgn_ind

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
    write(chFile,"(A,I10.10)") trim(DEM_opt%ResultsDir)//"PartVisuFor"//trim(DEM_opt%RunName),itime
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_8,ierror)  ! guarantee overwriting
    call MPI_BARRIER(Prtcl_WORLD,ierror)
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
        realVec(i)= two*GPrtcl_PosR(i)%w
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
    endif
#endif     
    call MPI_FILE_CLOSE(fh, ierror) 
    call MPI_COMM_FREE( Prtcl_WORLD, ierror)
  end subroutine PIO_Dump_visu

  !**********************************************************************
  ! clc_bgn_ind
  !**********************************************************************
  function clc_bgn_ind(nlocal) result(bgn_ind)
    implicit none
    integer,intent(in)::nlocal
    integer::bgn_ind

    ! locals
    integer::ierror
    integer,dimension(MPI_STATUS_SIZE) :: status1
  
    bgn_ind=0
    if(nproc<=1) return
    IF(nrank==0)THEN
      call MPI_SEND(bgn_ind+nlocal, 1, int_type, nrank+1, 0, MPI_COMM_WORLD, ierror)
    ELSEIF(nrank /= nproc-1) THEN
      call MPI_RECV(bgn_ind,        1, int_type, nrank-1, 0, MPI_COMM_WORLD, status1,  ierror)
      call MPI_SEND(bgn_ind+nlocal, 1, int_type, nrank+1, 0, MPI_COMM_WORLD, ierror)
    ELSE
      call MPI_RECV(bgn_ind,        1, int_type, nrank-1, 0, MPI_COMM_WORLD, status1,  ierror)
    ENDIF 
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

  !**********************************************************************
  ! Prtcl_dump_real3_matrix
  !**********************************************************************
  subroutine Prtcl_dump_real3_matrix(fh,disp,var,pmsize)
    implicit none
    integer,intent(in)::fh
    type(part_io_size_mat),intent(in)::pmsize    
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    type(real3),dimension(1:pmsize%subsizes(1),1:pmsize%subsizes(2)),intent(in)::var

    ! locals
    integer :: ierror,newtype
    integer, dimension(2) :: sizes, subsizes, starts

    ! calculate sizes, subsizes and starts
    sizes     = pmsize%sizes
    subsizes  = pmsize%subsizes
    starts    = pmsize%starts

    ! write the particle relevant real matrix
    call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real3_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,real3_type, newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1)*subsizes(2),real3_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp = disp + sizes(1) * sizes(2) * real3_byte
  end subroutine Prtcl_dump_real3_matrix
end module Prtcl_IOAndVisu
