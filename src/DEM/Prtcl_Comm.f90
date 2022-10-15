module Prtcl_Comm
  use MPI
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Property
  use Prtcl_Decomp_2d
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_Parameters
  use Prtcl_ContactSearchPW
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
#endif
  implicit none
  private

#define xm_axis 1
#define xp_axis 2
#define ym_axis 3
#define yp_axis 4
#define zm_axis 5
#define zp_axis 6

  type(real3)::simLen
  logical,dimension(3)::pbc
  real(RK),dimension(6)::dx_pbc
  real(RK),dimension(6)::dy_pbc
  real(RK),dimension(6)::dz_pbc
  real(RK):: xst0_cs,xst1_cs,xst2_cs
  real(RK):: xed0_cs,xed1_cs,xed2_cs
  real(RK):: yst0_cs,yst1_cs,yst2_cs
  real(RK):: yed0_cs,yed1_cs,yed2_cs
  real(RK):: zst0_cs,zst1_cs,zst2_cs
  real(RK):: zed0_cs,zed1_cs,zed2_cs

  type Prtcl_Comm_info
    integer :: msend
    integer :: GhostCS_size
    integer :: Prtcl_Exchange_size
    real(RK):: LenForCS
  contains
    procedure:: InitComm            => PC_InitComm
    procedure:: Comm_For_Cntct      => PC_Comm_For_Cntct
    procedure:: pack_cntct          => PC_pack_cntct
    procedure:: unpack_cntct        => PC_unpack_cntct
    procedure:: Comm_For_Exchange   => PC_Comm_For_Exchange
    procedure:: pack_Exchange       => PC_pack_Exchange
    procedure:: unpack_Exchange     => PC_unpack_Exchange
    procedure:: ISInThisProc        => PC_ISInThisProc
    procedure:: reallocate_sendlist => PC_reallocate_sendlist
    procedure:: reallocate_ghost_for_Cntct => PC_reallocate_ghost_for_Cntct

    procedure,private:: Comm_For_Cntct_fixed  => PC_Comm_For_Cntct_fixed
  end type Prtcl_Comm_info
  type(Prtcl_Comm_info),public::DEM_Comm
#ifdef CFDACM
  integer,allocatable,dimension(:),public::sendlist
#else
  integer,allocatable,dimension(:)::sendlist
#endif
  integer,allocatable,dimension(:),public:: GhostPFix_id
  integer,allocatable,dimension(:),public:: GhostPFix_pType
  type(real4),allocatable,dimension(:),public:: GhostPFix_PosR

contains

  !**********************************************************************
  ! PC_InitComm
  !**********************************************************************
  subroutine PC_InitComm(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    real(RK)::lx,ly,lz,le_cs,vol1,vol2,vol3
    real(RK)::maxD,vol_tot,vol_ghost_cs,vol_sendlist
    integer::numPrtcl,mCS,msend,mFixedCS,ierrTmp,ierror=0
        
    pbc=DEM_Opt%IsPeriodic
    simLen = DEM_Opt%SimDomain_max-DEM_Opt%SimDomain_min

    maxD = two*maxval(DEMProperty%Prtcl_PureProp%Radius)
#ifdef CFDACM
    this%LenForCS= DEM_Opt%Prtcl_cs_ratio*maxD+maxval(dlub_pp)
#else
    this%LenForCS= DEM_Opt%Prtcl_cs_ratio*maxD
#endif
    if(1.02_RK*this%LenForCS>min(min(simLen%x,simLen%y),simLen%z) ) call DEMLogInfo%CheckForError(ErrT_Abort,"PC_InitComm","so big Diameter")
    le_cs = two*this%LenForCS

    ! (id 1) +(ptype 1) +(PosR 3) +(linvel 3) +(rotvel 3)=11
    this%GhostCS_size   = 11

    ! (id 1) +(ptype 1) +(Mark 1) +(PosR   3) +(linvel 3*tsize) +(linAcc 3*tsize) + &
    !                              (theta  3) +(rotVel 3*rsize) +(rotAcc 3*rsize) = 9+6*(tsize+rsize)
    this%Prtcl_Exchange_size = 9 + 6*(GPrtcl_list%tsize + GPrtcl_list%rsize)
#ifdef CFDDEM
    ! (GPrtcl_FpForce 3)  +(GPrtcl_FpForce_old 3) +(GPrtcl_Vfluid 3+3) +(GPrtcl_linVelOld 3)=15
    this%Prtcl_Exchange_size = this%Prtcl_Exchange_size + 15
    if(is_clc_Basset) this%Prtcl_Exchange_size= this%Prtcl_Exchange_size+ 3*GPrtcl_BassetSeq%nDataLen
#endif
#ifdef CFDACM
    ! FpForce(3)+ FpTorque(3)+ FluidIntOld( 2*3)= 6+ 6 =12
    this%Prtcl_Exchange_size= this%Prtcl_Exchange_size +12
    ! PosOld(3)
    if(IBM_Scheme==2) this%Prtcl_Exchange_size= this%Prtcl_Exchange_size +3
#endif

    lx=DEM_decomp%xEd-DEM_decomp%xSt
    ly=DEM_decomp%yEd-DEM_decomp%ySt
    lz=DEM_decomp%zEd-DEM_decomp%zSt
    xst0_cs = DEM_decomp%xSt-this%LenForCS
    xst1_cs = DEM_decomp%xSt
    xst2_cs = DEM_decomp%xSt+this%LenForCS    
    xed0_cs = DEM_decomp%xEd-this%LenForCS
    xed1_cs = DEM_decomp%xEd
    xed2_cs = DEM_decomp%xEd+this%LenForCS

    yst0_cs = DEM_decomp%ySt-this%LenForCS
    yst1_cs = DEM_decomp%ySt
    yst2_cs = DEM_decomp%ySt+this%LenForCS    
    yed0_cs = DEM_decomp%yEd-this%LenForCS
    yed1_cs = DEM_decomp%yEd
    yed2_cs = DEM_decomp%yEd+this%LenForCS 

    zst0_cs = DEM_decomp%zSt-this%LenForCS
    zst1_cs = DEM_decomp%zSt
    zst2_cs = DEM_decomp%zSt+this%LenForCS    
    zed0_cs = DEM_decomp%zEd-this%LenForCS
    zed1_cs = DEM_decomp%zEd
    zed2_cs = DEM_decomp%zEd+this%LenForCS

    dx_pbc=zero; dy_pbc=zero; dz_pbc=zero
    IF(DEM_decomp%Prtcl_Pencil==x_axis)THEN
      if(pbc(1)) then
        dx_pbc(xm_axis)= simLen%x
        dx_pbc(xp_axis)=-simLen%x     
      endif
      if(pbc(2)) then
        if(DEM_decomp%coord1==0)                 dy_pbc(ym_axis)= simLen%y
        if(DEM_decomp%coord1==DEM_decomp%prow-1) dy_pbc(yp_axis)=-simLen%y
      endif
      if(pbc(3)) then
        if(DEM_decomp%coord2==0)                 dz_pbc(zm_axis)= simLen%z
        if(DEM_decomp%coord2==DEM_decomp%pcol-1) dz_pbc(zp_axis)=-simLen%z
      endif
    ELSEIF(DEM_decomp%Prtcl_Pencil==y_axis)THEN
      if(pbc(1)) then
        if(DEM_decomp%coord1==0)                 dx_pbc(xm_axis)= simLen%x
        if(DEM_decomp%coord1==DEM_decomp%prow-1) dx_pbc(xp_axis)=-simLen%x 
      endif
      if(pbc(2)) then
        dy_pbc(ym_axis)= simLen%y
        dy_pbc(yp_axis)=-simLen%y
      endif
      if(pbc(3)) then
        if(DEM_decomp%coord2==0)                 dz_pbc(zm_axis)= simLen%z
        if(DEM_decomp%coord2==DEM_decomp%pcol-1) dz_pbc(zp_axis)=-simLen%z
      endif
    ELSEIF(DEM_decomp%Prtcl_Pencil==z_axis)THEN
      if(pbc(1)) then
        if(DEM_decomp%coord1==0)                 dx_pbc(xm_axis)= simLen%x
        if(DEM_decomp%coord1==DEM_decomp%prow-1) dx_pbc(xp_axis)=-simLen%x 
      endif
      if(pbc(2)) then
        if(DEM_decomp%coord2==0)                 dy_pbc(ym_axis)= simLen%y
        if(DEM_decomp%coord2==DEM_decomp%pcol-1) dy_pbc(yp_axis)=-simLen%y
      endif
      if(pbc(3)) then
        dz_pbc(zm_axis)= simLen%z
        dz_pbc(zp_axis)=-simLen%z
      endif
    ENDIF

    vol_tot= SimLen%x * SimLen%y * SimLen%z
    vol_ghost_cs=(lx+le_cs)*(ly+le_cs)*(lz+le_cs)-lx*ly*lz

    numPrtcl = DEM_opt%numPrtcl
    mCS= int(numPrtcl*real(vol_ghost_cs/vol_tot,RK))
    mCS= 2*mCS
    mCS= min(mCS, numPrtcl)
    mCS= max(mCS, 10)
    GPrtcl_list%mGhost_CS = mCS

    mFixedCS= int(DEM_opt%numPrtclFix*real(vol_ghost_cs/vol_tot,RK))
    mFixedCS= 2*mFixedCS
    mFixedCS= min(mFixedCS, DEM_opt%numPrtclFix)
    mFixedCS= max(mFixedCS, 10)
    GPrtcl_list%nGhostFix_CS= mFixedCS

    vol1=(lx+le_cs)*(ly+le_cs)*this%LenForCS
    vol2=(lx+le_cs)*(lz+le_cs)*this%LenForCS
    vol3=(ly+le_cs)*(lz+le_cs)*this%LenForCS
    vol_sendlist=max(max(vol1,vol2),vol3)
    msend=int(DEM_opt%numPrtclFix*real(vol_sendlist/vol_tot,RK))
    msend=2*msend
    msend=min(msend,DEM_opt%numPrtclFix)
    msend=max(msend,10)
    this%msend=msend
    if(allocated(sendlist)) deallocate(sendlist);
    allocate(GhostPFix_id(mFixedCS),   Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostPFix_pType(mFixedCS),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostPFix_PosR(mFixedCS), Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(sendlist(msend),          Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"PC_InitComm","Allocation failed-1")
    call this%Comm_For_Cntct_fixed()

    msend=int(numPrtcl*real(vol_sendlist/vol_tot,RK))
    msend=2*msend
    msend=min(msend,numPrtcl)
    msend=max(msend,10)
    this%msend=msend
    if(allocated(sendlist)) deallocate(sendlist);
    allocate(GhostP_id(mCs),    Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostP_pType(mCs), Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostP_PosR(mCs),  Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostP_linVel(mCs),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(GhostP_rotVel(mCS),Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    allocate(sendlist(msend),   Stat=ierrTmp); ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"PC_InitComm","Allocation failed-2")
  end subroutine PC_InitComm

  !**********************************************************************
  ! Prtcl_Comm_For_Cntct
  !**********************************************************************
  subroutine PC_Comm_For_Cntct(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    real(RK)::px,py,pz
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer::nsend,nsend2,nrecv,nrecv2,nsendg,ng,ngp,ngpp
    integer::i,ierror,request(4),SRstatus(MPI_STATUS_SIZE)
    
    ng=0
    SELECT CASE(DEM_decomp%Prtcl_Pencil)
    CASE(x_axis)     ! ccccccccccccccccccccccccc  x-axis  ccccccccccccccccccccccccccc

      ! step1: Handle x-dir
      IF(pbc(1)) THEN
        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px+simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px-simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo
      ENDIF

      ! step2: send to yp_axis, and receive from ym_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%y = py-simLen%y
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py-simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py >=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py >=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,yp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to ym_axis, and receive from yp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%y = py+simLen%y
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_pType(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py+simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py <=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py <=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,ym_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz>=zed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%z = pz-simLen%z
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz-simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,zp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%z = pz+simLen%z
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz+simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,zm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(y_axis)     ! ccccccccccccccccccccccccc  y-axis  ccccccccccccccccccccccccccc

      ! step1: Handle y-dir
      IF(pbc(2)) THEN
        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py+simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py-simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px>=xed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%x = px-simLen%x
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px-simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,xp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%x = px+simLen%x
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_pType(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px+simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,xm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz>=zed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%z = pz-simLen%z
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz-simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,zp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%z = pz+simLen%z
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz+simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostP_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,zm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(z_axis)     ! ccccccccccccccccccccccccc  z-axis  ccccccccccccccccccccccccccc

      ! step1: Handle z-dir
      IF(pbc(3)) THEN
        do i=1,GPrtcl_list%nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz+simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%z = pz-simLen%z
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px>=xed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%x = px-simLen%x
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px-simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,xp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_pType(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%x = px+simLen%x
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_pType(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%x = px+simLen%x
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostP_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%nlocal
          px=GPrtcl_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,xm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to yp_axis, and receive from ym_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py>=yed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%y = py-simLen%y
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py-simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py>=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py>=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,yp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>=GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to ym_axis, and receive from yp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GhostP_id(i)
            GhostP_pType(ng)  = GhostP_ptype(i)
            GhostP_PosR(ng)   = GhostP_PosR(i)
            GhostP_PosR(ng)%y = py+simLen%y
            GhostP_linVel(ng) = GhostP_linVel(i)
            GhostP_rotVel(ng) = GhostP_rotVel(i)
          endif
        enddo

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
            GhostP_id(ng)     = GPrtcl_id(i)
            GhostP_pType(ng)  = GPrtcl_ptype(i)
            GhostP_PosR(ng)   = GPrtcl_PosR(i)
            GhostP_PosR(ng)%y = py+simLen%y
            GhostP_linVel(ng) = GPrtcl_linVel(1,i)
            GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostP_PosR(i)%y
          if(py<=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%nlocal
          py=GPrtcl_PosR(i)%y
          if(py<=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*this%GhostCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*this%GhostCS_size
        allocate(buf_send(nsend2))
        call this%pack_cntct(buf_send,nsendg,nsend,ym_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%mGhost_CS) call this%reallocate_ghost_for_Cntct(ng)
        call this%unpack_cntct(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)
    END SELECT
    GPrtcl_list%nGhost_CS = ng
  end subroutine PC_Comm_For_Cntct

  !**********************************************************************
  ! PC_pack_cntct
  !**********************************************************************
  subroutine PC_pack_cntct(this,buf_send,nsendg,nsend,dir)
    implicit none
    class(Prtcl_Comm_info)::this
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsendg,nsend,dir

    ! locals
    integer::i,id,m
    real(RK)::dx,dy,dz

    m=1
    dx=dx_pbc(dir)
    dy=dy_pbc(dir)
    dz=dz_pbc(dir)
    do i=1,nsendg
      id=sendlist(i)
      buf_send(m)=real(GhostP_id(id));    m=m+1 ! 01
      buf_send(m)=real(GhostP_pType(id)); m=m+1 ! 02
      buf_send(m)=GhostP_PosR(id)%x+dx;   m=m+1 ! 03
      buf_send(m)=GhostP_PosR(id)%y+dy;   m=m+1 ! 04
      buf_send(m)=GhostP_PosR(id)%z+dz;   m=m+1 ! 05
      buf_send(m)=GhostP_linVel(id)%x;    m=m+1 ! 06
      buf_send(m)=GhostP_linVel(id)%y;    m=m+1 ! 07
      buf_send(m)=GhostP_linVel(id)%z;    m=m+1 ! 08
      buf_send(m)=GhostP_rotVel(id)%x;    m=m+1 ! 09
      buf_send(m)=GhostP_rotVel(id)%y;    m=m+1 ! 10
      buf_send(m)=GhostP_rotVel(id)%z;    m=m+1 ! 11 
    enddo
    do i=nsendg+1,nsend
      id=sendlist(i)
      buf_send(m)=real(GPrtcl_id(id));    m=m+1 ! 01
      buf_send(m)=real(GPrtcl_pType(id)); m=m+1 ! 02
      buf_send(m)=GPrtcl_PosR(id)%x+dx;   m=m+1 ! 03
      buf_send(m)=GPrtcl_PosR(id)%y+dy;   m=m+1 ! 04
      buf_send(m)=GPrtcl_PosR(id)%z+dz;   m=m+1 ! 05
      buf_send(m)=GPrtcl_linVel(1,id)%x;  m=m+1 ! 06
      buf_send(m)=GPrtcl_linVel(1,id)%y;  m=m+1 ! 07
      buf_send(m)=GPrtcl_linVel(1,id)%z;  m=m+1 ! 08
      buf_send(m)=GPrtcl_rotVel(1,id)%x;  m=m+1 ! 09
      buf_send(m)=GPrtcl_rotVel(1,id)%y;  m=m+1 ! 10
      buf_send(m)=GPrtcl_rotVel(1,id)%z;  m=m+1 ! 11 
    enddo
  end subroutine PC_pack_cntct

  !**********************************************************************
  ! PC_unpack_cntct
  !**********************************************************************
  subroutine PC_unpack_cntct(this,buf_recv,n1,n2)
    implicit none
    class(Prtcl_Comm_info)::this
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::n1,n2

    ! locals
    integer::i,m,itype
   
    m=1
    do i=n1,n2
      GhostP_id(i)       =int(buf_recv(m)+0.2); m=m+1 ! 01 
      itype              =int(buf_recv(m)+0.2); m=m+1 ! 02
      GhostP_PosR(i)%x   =buf_recv(m);          m=m+1 ! 03
      GhostP_PosR(i)%y   =buf_recv(m);          m=m+1 ! 04
      GhostP_PosR(i)%z   =buf_recv(m);          m=m+1 ! 05
      GhostP_linVel(i)%x =buf_recv(m);          m=m+1 ! 06 
      GhostP_linVel(i)%y =buf_recv(m);          m=m+1 ! 07 
      GhostP_linVel(i)%z =buf_recv(m);          m=m+1 ! 08 
      GhostP_rotVel(i)%x =buf_recv(m);          m=m+1 ! 09 
      GhostP_rotVel(i)%y =buf_recv(m);          m=m+1 ! 10 
      GhostP_rotVel(i)%z =buf_recv(m);          m=m+1 ! 11
      GhostP_pType(i)    =itype;                
      GhostP_PosR(i)%w   =DEMProperty%Prtcl_PureProp(itype)%Radius 
    enddo
  end subroutine PC_unpack_cntct

  !**********************************************************************
  ! PC_Comm_For_Exchange
  !**********************************************************************
  subroutine PC_Comm_For_Exchange(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    integer::i,ierror,request(4),nlocal,nlocalp
    integer::nsend(2),nrecv(2),nlink
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer,dimension(MPI_STATUS_SIZE) :: SRstatus
    real(RK)::px,py,pz

    nlocal=GPrtcl_list%nlocal
    SELECT CASE(DEM_decomp%Prtcl_Pencil)
    CASE(x_axis)     ! ccccccccccccccccccccccccc  x-axis  ccccccccccccccccccccccccccc

      ! step1: Handle x-dir
      IF(pbc(1)) THEN
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            GPrtcl_PosR(i)%x=px - simLen%x
          elseif(px < xst1_cs) then
            GPrtcl_PosR(i)%x=px + simLen%x
          endif
        enddo
      ELSE
        i=1
        do while(i<=nlocal)
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs .or. px < xst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo
      ENDIF

      ! step2: send to yp_axis, and receive from ym_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(3)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            GPrtcl_PosR(i)%y=py - simLen%y
          endif
        enddo

      ELSE
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(3), 9, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(4), 9, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(4),10,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,yp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(3),10,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,yp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to ym_axis, and receive from yp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(4)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            GPrtcl_PosR(i)%y=py + simLen%y
          endif
        enddo

      ELSE
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(4),11, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(3),11, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(3),12,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,ym_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(4),12,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,ym_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(1)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            GPrtcl_PosR(i)%z=pz - simLen%z
          endif
        enddo

      ELSE
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(1),13, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(2),13, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(2),14,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,zp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(1),14,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(2)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            GPrtcl_PosR(i)%z=pz + simLen%z
          endif
        enddo

      ELSE
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(2),15, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(1),15, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(1),16,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,zm_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(2),16,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zm_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(y_axis)     ! ccccccccccccccccccccccccc  y-axis  ccccccccccccccccccccccccccc

      ! step1: Handle y-dir
      IF(pbc(2)) THEN
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            GPrtcl_PosR(i)%y=py - simLen%y
          elseif(py < yst1_cs) then
            GPrtcl_PosR(i)%y=py + simLen%y
          endif
        enddo
      ELSE
        i=1
        do while(i<=nlocal)
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs .or. py < yst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(3)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            GPrtcl_PosR(i)%x=px - simLen%x
          endif
        enddo

      ELSE
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(3), 9, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(4), 9, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(4),10,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,xp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(3),10,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(4)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            GPrtcl_PosR(i)%x=px + simLen%x
          endif
        enddo

      ELSE
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(4),11, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(3),11, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(3),12,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,xm_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(4),12,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xm_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(1)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            GPrtcl_PosR(i)%z=pz - simLen%z
          endif
        enddo

      ELSE
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(1),13, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(2),13, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(2),14,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,zp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(1),14,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(2)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            GPrtcl_PosR(i)%z=pz + simLen%z
          endif
        enddo

      ELSE
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz < zst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(2),15, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(1),15, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(1),16,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,zm_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(2),16,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zm_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(z_axis)     ! ccccccccccccccccccccccccc  z-axis  ccccccccccccccccccccccccccc

      ! step1: Handle z-dir
      IF(pbc(3)) THEN
        do i=1,nlocal
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs) then
            GPrtcl_PosR(i)%z=pz - simLen%z
          elseif(pz < zst1_cs) then
            GPrtcl_PosR(i)%z=pz + simLen%z
          endif
        enddo
      ELSE
        i=1
        do while(i<=nlocal)
          pz=GPrtcl_PosR(i)%z
          if(pz >= zed1_cs .or. pz < zst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(3)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            GPrtcl_PosR(i)%x=px - simLen%x
          endif
        enddo

      ELSE
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px >= xed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(3), 9, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(4), 9, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(4),10,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,xp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(3),10,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(4)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            GPrtcl_PosR(i)%x=px + simLen%x
          endif
        enddo

      ELSE
        do i=1,nlocal
          px=GPrtcl_PosR(i)%x
          if(px < xst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(4),11, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(3),11, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(3),12,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,xm_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(4),12,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xm_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to yp_axis, and receive from ym_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(1)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            GPrtcl_PosR(i)%y=py - simLen%y
          endif
        enddo

      ELSE
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py >= yed1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(1),13, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(2),13, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(2),14,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,yp_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(1),14,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,yp_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to ym_axis, and receive from yp_dir
      nsend =0; nrecv =0
      IF(DEM_decomp%ProcNgh(2)==MPI_PROC_NULL) THEN
        i=1
        do while(i<=nlocal)
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            call GPrtcl_list%copy(i,nlocal)
            call GPPW_CntctList%copy(i,nlocal)
            call DEMContactSearchPW%copy(i,nlocal)
            nlocal = nlocal -1
          else
            i = i + 1
          endif
        enddo

      ELSEIF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            GPrtcl_PosR(i)%y=py + simLen%y
          endif
        enddo

      ELSE
        do i=1,nlocal
          py=GPrtcl_PosR(i)%y
          if(py < yst1_cs) then
            nsend(1)= nsend(1)+1
            if(nsend(1) > this%msend) call this%reallocate_sendlist(nsend(1))
            sendlist(nsend(1)) = i
            nlink = GPPW_CntctList%getPrtcl_nlink(i)
            nsend(2)= nsend(2) +nlink*5+ this%Prtcl_Exchange_size+1 ! The final "1" stands for END_OF_PRTCL
          endif
        enddo

      ENDIF
      call MPI_SENDRECV(nsend, 2, int_type, DEM_decomp%ProcNgh(2),15, &
                        nrecv, 2, int_type, DEM_decomp%ProcNgh(1),15, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv(2)>0) then
        allocate(buf_recv(nrecv(2)))
        call MPI_IRECV(buf_recv,nrecv(2),real_type,DEM_decomp%ProcNgh(1),16,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend(2)>0) then
        allocate(buf_send(nsend(2)))
        nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
        call this%pack_Exchange(buf_send,nsend(1),nlocalp,ym_axis)
        call MPI_SEND(buf_send,nsend(2),real_type,DEM_decomp%ProcNgh(2),16,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv(2)>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
          call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
          call GPPW_CntctList%reallocateCL(nlocal + nrecv(1))
          call DEMContactSearchPW%Reallocate_Bucket(nlocal + nrecv(1))
        endif
        call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,ym_axis)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)
    END SELECT
    GPrtcl_list%nlocal = nlocal
  end subroutine PC_Comm_For_Exchange

  !**********************************************************************
  ! PC_pack_Exchange
  !**********************************************************************
  subroutine PC_pack_Exchange(this,buf_send,nsend,nlocalp,dir)
    implicit none
    class(Prtcl_Comm_info)::this
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsend,dir,nlocalp

    ! locals
    real(RK)::dx,dy,dz
    integer::i,j,id,m,nlocal

    m=1
    nlocal = nlocalp
    dx=dx_pbc(dir)
    dy=dy_pbc(dir)
    dz=dz_pbc(dir)

#ifdef CFDDEM
    IF(is_clc_Basset)  THEN
      DO i=1,nsend
        id = sendlist(i)
        buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
        buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
        buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
        buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
        buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
        buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
        buf_send(m)=GPrtcl_theta(id)%x;       m=m+1 ! 07
        buf_send(m)=GPrtcl_theta(id)%y;       m=m+1 ! 08
        buf_send(m)=GPrtcl_theta(id)%z;       m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
          buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
        enddo
        do j=1,GPrtcl_list%rsize
          buf_send(m)=GPrtcl_rotVel(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%z;  m=m+1 ! 6* rsize 
        enddo
        buf_send(m)=GPrtcl_FpForce(id)%x;     m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%y;     m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%z;     m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%x; m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%y; m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%z; m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%x;   m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%y;   m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%z;   m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%x;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%y;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%z;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%x;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%y;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%z;    m=m+1
        do j=1, GPrtcl_BassetSeq%nDataLen
          buf_send(m)=GPrtcl_BassetData(j,id)%x; m=m+1
          buf_send(m)=GPrtcl_BassetData(j,id)%y; m=m+1
          buf_send(m)=GPrtcl_BassetData(j,id)%z; m=m+1
        enddo
        call GPPW_CntctList%Gather_Cntctlink(id,buf_send,m) ! contact list part
        buf_send(m) = END_OF_PRTCL;           m=m+1         ! END_OF_PRTCL
      ENDDO
    ELSE
      DO i=1,nsend
        id = sendlist(i)
        buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
        buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
        buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
        buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
        buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
        buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
        buf_send(m)=GPrtcl_theta(id)%x;       m=m+1 ! 07
        buf_send(m)=GPrtcl_theta(id)%y;       m=m+1 ! 08
        buf_send(m)=GPrtcl_theta(id)%z;       m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
          buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
        enddo
        do j=1,GPrtcl_list%rsize
          buf_send(m)=GPrtcl_rotVel(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%z;  m=m+1 ! 6* rsize 
        enddo
        buf_send(m)=GPrtcl_FpForce(id)%x;     m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%y;     m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%z;     m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%x; m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%y; m=m+1
        buf_send(m)=GPrtcl_FpForce_old(id)%z; m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%x;   m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%y;   m=m+1
        buf_send(m)=GPrtcl_linVelOld(id)%z;   m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%x;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%y;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(1,id)%z;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%x;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%y;    m=m+1
        buf_send(m)=GPrtcl_Vfluid(2,id)%z;    m=m+1
        call GPPW_CntctList%Gather_Cntctlink(id,buf_send,m) ! contact list part
        buf_send(m) = END_OF_PRTCL;           m=m+1         ! END_OF_PRTCL
      ENDDO
    ENDIF
#elif CFDACM
    if(IBM_Scheme==2) then
      DO i=1,nsend
        id = sendlist(i)
        buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
        buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
        buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
        buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
        buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
        buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
        buf_send(m)=GPrtcl_theta(id)%x;       m=m+1 ! 07
        buf_send(m)=GPrtcl_theta(id)%y;       m=m+1 ! 08
        buf_send(m)=GPrtcl_theta(id)%z;       m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
          buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
        enddo
        do j=1,GPrtcl_list%rsize
          buf_send(m)=GPrtcl_rotVel(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%z;  m=m+1 ! 6* rsize 
        enddo
        buf_send(m)=GPrtcl_FpForce(id)%x;       m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%y;       m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%z;       m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%x;      m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%y;      m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%z;      m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%x; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%y; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%z; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%x; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%y; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%z; m=m+1
        buf_send(m)=GPrtcl_PosOld(id)%x+dx;     m=m+1
        buf_send(m)=GPrtcl_PosOld(id)%y+dy;     m=m+1
        buf_send(m)=GPrtcl_PosOld(id)%z+dz;     m=m+1
        call GPPW_CntctList%Gather_Cntctlink(id,buf_send,m) ! contact list part
        buf_send(m) = END_OF_PRTCL;           m=m+1         ! END_OF_PRTCL
      ENDDO
    else
      DO i=1,nsend
        id = sendlist(i)
        buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
        buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
        buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
        buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
        buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
        buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
        buf_send(m)=GPrtcl_theta(id)%x;       m=m+1 ! 07
        buf_send(m)=GPrtcl_theta(id)%y;       m=m+1 ! 08
        buf_send(m)=GPrtcl_theta(id)%z;       m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
          buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
        enddo
        do j=1,GPrtcl_list%rsize
          buf_send(m)=GPrtcl_rotVel(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotVel(j,id)%z;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%x;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%y;  m=m+1 ! 
          buf_send(m)=GPrtcl_rotAcc(j,id)%z;  m=m+1 ! 6* rsize 
        enddo
        buf_send(m)=GPrtcl_FpForce(id)%x;       m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%y;       m=m+1
        buf_send(m)=GPrtcl_FpForce(id)%z;       m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%x;      m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%y;      m=m+1
        buf_send(m)=GPrtcl_FpTorque(id)%z;      m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%x; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%y; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(1,id)%z; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%x; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%y; m=m+1
        buf_send(m)=GPrtcl_FluidIntOld(2,id)%z; m=m+1
        call GPPW_CntctList%Gather_Cntctlink(id,buf_send,m) ! contact list part
        buf_send(m) = END_OF_PRTCL;           m=m+1         ! END_OF_PRTCL
      ENDDO
    endif
#else
   DO i=1,nsend
      id = sendlist(i)
      buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
      buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
      buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
      buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
      buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
      buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
      buf_send(m)=GPrtcl_theta(id)%x;       m=m+1 ! 07
      buf_send(m)=GPrtcl_theta(id)%y;       m=m+1 ! 08
      buf_send(m)=GPrtcl_theta(id)%z;       m=m+1 ! 09
      do j=1,GPrtcl_list%tsize
        buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
        buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
      enddo
      do j=1,GPrtcl_list%rsize
        buf_send(m)=GPrtcl_rotVel(j,id)%x;  m=m+1 ! 
        buf_send(m)=GPrtcl_rotVel(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_rotVel(j,id)%z;  m=m+1 ! 
        buf_send(m)=GPrtcl_rotAcc(j,id)%x;  m=m+1 ! 
        buf_send(m)=GPrtcl_rotAcc(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_rotAcc(j,id)%z;  m=m+1 ! 6* rsize 
      enddo
      call GPPW_CntctList%Gather_Cntctlink(id,buf_send,m) ! contact list part
      buf_send(m) = END_OF_PRTCL;           m=m+1         ! END_OF_PRTCL
    ENDDO
#endif

    DO i=nsend,1,-1
      id = sendlist(i)
      call GPrtcl_list%copy(id,nlocal)
      call GPPW_CntctList%copy(id,nlocal)
      call DEMContactSearchPW%copy(id,nlocal)
      nlocal = nlocal -1
    ENDDO
  end subroutine PC_pack_Exchange

  !**********************************************************************
  ! PC_unpack_Exchange
  !**********************************************************************
  subroutine PC_unpack_Exchange(this,buf_recv,nrecv,nlocal,dir)
    implicit none
    class(Prtcl_Comm_info)::this
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::nrecv,dir
    integer,intent(inout)::nlocal

    ! locals
    integer::i,j,id,m,itype

    m=1
    id=nlocal+1
#ifdef CFDDEM
    If (is_clc_Basset) then
      DO i=1,nrecv
        if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
        GPrtcl_id(id)     = int(buf_recv(m)+0.2); m=m+1 ! 01
        itype             = int(buf_recv(m)+0.2); m=m+1 ! 02
        GPrtcl_usrMark(id)= int(buf_recv(m)+0.2); m=m+1 ! 03
        GPrtcl_PosR(id)%x = buf_recv(m);          m=m+1 ! 04
        GPrtcl_PosR(id)%y = buf_recv(m);          m=m+1 ! 05
        GPrtcl_PosR(id)%z = buf_recv(m);          m=m+1 ! 06
        GPrtcl_theta(id)%x= buf_recv(m);          m=m+1 ! 07
        GPrtcl_theta(id)%y= buf_recv(m);          m=m+1 ! 08
        GPrtcl_theta(id)%z= buf_recv(m);          m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          GPrtcl_linVel(j,id)%x = buf_recv(m);    m=m+1 ! 6* tsize
          GPrtcl_linVel(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_linVel(j,id)%z = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%z = buf_recv(m);    m=m+1
        enddo
        do j=1,GPrtcl_list%rsize
          GPrtcl_rotVel(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_rotVel(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_rotVel(j,id)%z = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%z = buf_recv(m);    m=m+1 ! 6* rsize
        enddo
        GPrtcl_FpForce(id)%x    =buf_recv(m);     m=m+1
        GPrtcl_FpForce(id)%y    =buf_recv(m);     m=m+1
        GPrtcl_FpForce(id)%z    =buf_recv(m);     m=m+1
        GPrtcl_FpForce_old(id)%x=buf_recv(m);     m=m+1
        GPrtcl_FpForce_old(id)%y=buf_recv(m);     m=m+1
        GPrtcl_FpForce_old(id)%z=buf_recv(m);     m=m+1
        GPrtcl_linVelOld(id)%x  =buf_recv(m);     m=m+1
        GPrtcl_linVelOld(id)%y  =buf_recv(m);     m=m+1
        GPrtcl_linVelOld(id)%z  =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(1,id)%x   =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(1,id)%y   =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(1,id)%z   =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(2,id)%x   =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(2,id)%y   =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(2,id)%z   =buf_recv(m);     m=m+1
        do j=1,GPrtcl_BassetSeq%nDataLen
          GPrtcl_BassetData(j,id)%x= buf_recv(m); m=m+1
          GPrtcl_BassetData(j,id)%y= buf_recv(m); m=m+1
          GPrtcl_BassetData(j,id)%z= buf_recv(m); m=m+1
        enddo
        GPrtcl_pType(id)=itype
        GPrtcl_PosR(id)%w =DEMProperty%Prtcl_PureProp(itype)%Radius
        call GPPW_CntctList%Add_Cntctlink(id,buf_recv,m) ! contact list part
        call DEMContactSearchPW%InsertNearPW(id)
        id = id + 1
        nlocal =nlocal + 1
      ENDDO
    ELSE

      DO i=1,nrecv
        if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
        GPrtcl_id(id)     = int(buf_recv(m)+0.2); m=m+1 ! 01
        itype             = int(buf_recv(m)+0.2); m=m+1 ! 02
        GPrtcl_usrMark(id)= int(buf_recv(m)+0.2); m=m+1 ! 03
        GPrtcl_PosR(id)%x = buf_recv(m);          m=m+1 ! 04
        GPrtcl_PosR(id)%y = buf_recv(m);          m=m+1 ! 05
        GPrtcl_PosR(id)%z = buf_recv(m);          m=m+1 ! 06
        GPrtcl_theta(id)%x= buf_recv(m);          m=m+1 ! 07
        GPrtcl_theta(id)%y= buf_recv(m);          m=m+1 ! 08
        GPrtcl_theta(id)%z= buf_recv(m);          m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          GPrtcl_linVel(j,id)%x = buf_recv(m);    m=m+1 ! 6* tsize
          GPrtcl_linVel(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_linVel(j,id)%z = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%z = buf_recv(m);    m=m+1
        enddo
        do j=1,GPrtcl_list%rsize
          GPrtcl_rotVel(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_rotVel(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_rotVel(j,id)%z = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%z = buf_recv(m);    m=m+1 ! 6* rsize
        enddo
        GPrtcl_FpForce(id)%x    =buf_recv(m);     m=m+1
        GPrtcl_FpForce(id)%y    =buf_recv(m);     m=m+1
        GPrtcl_FpForce(id)%z    =buf_recv(m);     m=m+1
        GPrtcl_FpForce_old(id)%x=buf_recv(m);     m=m+1
        GPrtcl_FpForce_old(id)%y=buf_recv(m);     m=m+1
        GPrtcl_FpForce_old(id)%z=buf_recv(m);     m=m+1
        GPrtcl_linVelOld(id)%x  =buf_recv(m);     m=m+1
        GPrtcl_linVelOld(id)%y  =buf_recv(m);     m=m+1
        GPrtcl_linVelOld(id)%z  =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(1,id)%x   =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(1,id)%y   =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(1,id)%z   =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(2,id)%x   =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(2,id)%y   =buf_recv(m);     m=m+1
        GPrtcl_Vfluid(2,id)%z   =buf_recv(m);     m=m+1
        GPrtcl_pType(id)=itype
        GPrtcl_PosR(id)%w =DEMProperty%Prtcl_PureProp(itype)%Radius

        call GPPW_CntctList%Add_Cntctlink(id,buf_recv,m) ! contact list part
        call DEMContactSearchPW%InsertNearPW(id)
 
        id = id + 1
        nlocal =nlocal + 1
      ENDDO
    ENDIF
#elif CFDACM
    if(IBM_Scheme==2) then
      DO i=1,nrecv
        if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
        GPrtcl_id(id)     = int(buf_recv(m)+0.2); m=m+1 ! 01
        itype             = int(buf_recv(m)+0.2); m=m+1 ! 02
        GPrtcl_usrMark(id)= int(buf_recv(m)+0.2); m=m+1 ! 03
        GPrtcl_PosR(id)%x = buf_recv(m);          m=m+1 ! 04
        GPrtcl_PosR(id)%y = buf_recv(m);          m=m+1 ! 05
        GPrtcl_PosR(id)%z = buf_recv(m);          m=m+1 ! 06
        GPrtcl_theta(id)%x= buf_recv(m);          m=m+1 ! 07
        GPrtcl_theta(id)%y= buf_recv(m);          m=m+1 ! 08
        GPrtcl_theta(id)%z= buf_recv(m);          m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          GPrtcl_linVel(j,id)%x = buf_recv(m);    m=m+1 ! 6* tsize
          GPrtcl_linVel(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_linVel(j,id)%z = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%z = buf_recv(m);    m=m+1
        enddo
        do j=1,GPrtcl_list%rsize
          GPrtcl_rotVel(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_rotVel(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_rotVel(j,id)%z = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%z = buf_recv(m);    m=m+1 ! 6* rsize
        enddo
        GPrtcl_FpForce(id)%x = buf_recv(m);       m=m+1
        GPrtcl_FpForce(id)%y = buf_recv(m);       m=m+1
        GPrtcl_FpForce(id)%z = buf_recv(m);       m=m+1
        GPrtcl_FpTorque(id)%x= buf_recv(m);       m=m+1
        GPrtcl_FpTorque(id)%y= buf_recv(m);       m=m+1
        GPrtcl_FpTorque(id)%z= buf_recv(m);       m=m+1
        GPrtcl_FluidIntold(1,id)%x=buf_recv(m);   m=m+1
        GPrtcl_FluidIntold(1,id)%y=buf_recv(m);   m=m+1
        GPrtcl_FluidIntold(1,id)%z=buf_recv(m);   m=m+1
        GPrtcl_FluidIntold(2,id)%x=buf_recv(m);   m=m+1
        GPrtcl_FluidIntold(2,id)%y=buf_recv(m);   m=m+1
        GPrtcl_FluidIntold(2,id)%z=buf_recv(m);   m=m+1
        GPrtcl_PosOld(id)%x=buf_recv(m);          m=m+1
        GPrtcl_PosOld(id)%y=buf_recv(m);          m=m+1
        GPrtcl_PosOld(id)%z=buf_recv(m);          m=m+1
        GPrtcl_pType(id)=itype
        GPrtcl_PosR(id)%w =DEMProperty%Prtcl_PureProp(itype)%Radius

        call GPPW_CntctList%Add_Cntctlink(id,buf_recv,m) ! contact list part
        call DEMContactSearchPW%InsertNearPW(id)
 
        id = id + 1
        nlocal =nlocal + 1
      ENDDO
    else
      DO i=1,nrecv
        if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
        GPrtcl_id(id)     = int(buf_recv(m)+0.2); m=m+1 ! 01
        itype             = int(buf_recv(m)+0.2); m=m+1 ! 02
        GPrtcl_usrMark(id)= int(buf_recv(m)+0.2); m=m+1 ! 03
        GPrtcl_PosR(id)%x = buf_recv(m);          m=m+1 ! 04
        GPrtcl_PosR(id)%y = buf_recv(m);          m=m+1 ! 05
        GPrtcl_PosR(id)%z = buf_recv(m);          m=m+1 ! 06
        GPrtcl_theta(id)%x= buf_recv(m);          m=m+1 ! 07
        GPrtcl_theta(id)%y= buf_recv(m);          m=m+1 ! 08
        GPrtcl_theta(id)%z= buf_recv(m);          m=m+1 ! 09
        do j=1,GPrtcl_list%tsize
          GPrtcl_linVel(j,id)%x = buf_recv(m);    m=m+1 ! 6* tsize
          GPrtcl_linVel(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_linVel(j,id)%z = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_linAcc(j,id)%z = buf_recv(m);    m=m+1
        enddo
        do j=1,GPrtcl_list%rsize
          GPrtcl_rotVel(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_rotVel(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_rotVel(j,id)%z = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%x = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%y = buf_recv(m);    m=m+1
          GPrtcl_rotAcc(j,id)%z = buf_recv(m);    m=m+1 ! 6* rsize
        enddo
        GPrtcl_FpForce(id)%x = buf_recv(m);       m=m+1
        GPrtcl_FpForce(id)%y = buf_recv(m);       m=m+1
        GPrtcl_FpForce(id)%z = buf_recv(m);       m=m+1
        GPrtcl_FpTorque(id)%x= buf_recv(m);       m=m+1
        GPrtcl_FpTorque(id)%y= buf_recv(m);       m=m+1
        GPrtcl_FpTorque(id)%z= buf_recv(m);       m=m+1
        GPrtcl_FluidIntold(1,id)%x=buf_recv(m);   m=m+1
        GPrtcl_FluidIntold(1,id)%y=buf_recv(m);   m=m+1
        GPrtcl_FluidIntold(1,id)%z=buf_recv(m);   m=m+1
        GPrtcl_FluidIntold(2,id)%x=buf_recv(m);   m=m+1
        GPrtcl_FluidIntold(2,id)%y=buf_recv(m);   m=m+1
        GPrtcl_FluidIntold(2,id)%z=buf_recv(m);   m=m+1
        GPrtcl_pType(id)=itype
        GPrtcl_PosR(id)%w =DEMProperty%Prtcl_PureProp(itype)%Radius

        call GPPW_CntctList%Add_Cntctlink(id,buf_recv,m) ! contact list part
        call DEMContactSearchPW%InsertNearPW(id)
 
        id = id + 1
        nlocal =nlocal + 1
      ENDDO
    endif
#else
    DO i=1,nrecv
      if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
      GPrtcl_id(id)     = int(buf_recv(m)+0.2); m=m+1 ! 01
      itype             = int(buf_recv(m)+0.2); m=m+1 ! 02
      GPrtcl_usrMark(id)= int(buf_recv(m)+0.2); m=m+1 ! 03
      GPrtcl_PosR(id)%x = buf_recv(m);          m=m+1 ! 04
      GPrtcl_PosR(id)%y = buf_recv(m);          m=m+1 ! 05
      GPrtcl_PosR(id)%z = buf_recv(m);          m=m+1 ! 06
      GPrtcl_theta(id)%x= buf_recv(m);          m=m+1 ! 07
      GPrtcl_theta(id)%y= buf_recv(m);          m=m+1 ! 08
      GPrtcl_theta(id)%z= buf_recv(m);          m=m+1 ! 09
      do j=1,GPrtcl_list%tsize
        GPrtcl_linVel(j,id)%x = buf_recv(m);    m=m+1 ! 6* tsize
        GPrtcl_linVel(j,id)%y = buf_recv(m);    m=m+1
        GPrtcl_linVel(j,id)%z = buf_recv(m);    m=m+1
        GPrtcl_linAcc(j,id)%x = buf_recv(m);    m=m+1
        GPrtcl_linAcc(j,id)%y = buf_recv(m);    m=m+1
        GPrtcl_linAcc(j,id)%z = buf_recv(m);    m=m+1
      enddo
      do j=1,GPrtcl_list%rsize
        GPrtcl_rotVel(j,id)%x = buf_recv(m);    m=m+1
        GPrtcl_rotVel(j,id)%y = buf_recv(m);    m=m+1
        GPrtcl_rotVel(j,id)%z = buf_recv(m);    m=m+1
        GPrtcl_rotAcc(j,id)%x = buf_recv(m);    m=m+1
        GPrtcl_rotAcc(j,id)%y = buf_recv(m);    m=m+1
        GPrtcl_rotAcc(j,id)%z = buf_recv(m);    m=m+1 ! 6* rsize
      enddo
      GPrtcl_pType(id)=itype
      GPrtcl_PosR(id)%w =DEMProperty%Prtcl_PureProp(itype)%Radius

      call GPPW_CntctList%Add_Cntctlink(id,buf_recv,m) ! contact list part
      call DEMContactSearchPW%InsertNearPW(id)
 
      id = id + 1
      nlocal =nlocal + 1
    ENDDO
#endif
  end subroutine PC_unpack_Exchange

  !**********************************************************************
  ! ISInThisProc
  !**********************************************************************
  function PC_ISInThisProc(this,buf_recv,m,dir) result(res)
    implicit none
    class(Prtcl_Comm_info)::this
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(inout)::m
    integer,intent(in)::dir
    logical:: res

    !local
    integer:: mp

    res=.true.
    SELECT CASE(dir)
    CASE(xp_axis)
      mp = m+3
      if(buf_recv(mp)>=xed1_cs) res = .false.
    CASE(xm_axis)
      mp = m+3
      if(buf_recv(mp)< xst1_cs) res = .false.
    CASE(yp_axis)
      mp = m+4
      if(buf_recv(mp)>=yed1_cs) res = .false. 
    CASE(ym_axis)
      mp = m+4
      if(buf_recv(mp)< yst1_cs) res = .false.
    CASE(zp_axis)
      mp = m+5
      if(buf_recv(mp)>=zed1_cs) res = .false.
    CASE(zm_axis)
      mp = m+5
      if(buf_recv(mp)< zst1_cs) res = .false.
    END SELECT
    if(res) return
    
    call DEMLogInfo%CheckForError(ErrT_Pass," PC_ISInThisProc"," The following particle is deleted: ")
    call DEMLogInfo%OutInfo(" Exchange direction  is :"//trim(num2str(dir)  ),3)
    call DEMLogInfo%OutInfo("   The particle id   is :"//trim(num2str( int(buf_recv(m)  +0.2))),3)
    call DEMLogInfo%OutInfo("       particle type is :"//trim(num2str( int(buf_recv(m+1)+0.2))),3)
    call DEMLogInfo%OutInfo("       x-coordinate  is :"//trim(num2str( buf_recv(m+3) )),   3)
    call DEMLogInfo%OutInfo("       y-coordinate  is :"//trim(num2str( buf_recv(m+4) )),   3)
    call DEMLogInfo%OutInfo("       z-coordinate  is :"//trim(num2str( buf_recv(m+5) )),   3)
    call DEMLogInfo%OutInfo(" Present processor   is :"//trim(num2str(nrank)), 3)

    m = m + this%Prtcl_Exchange_size
    do while(abs(buf_recv(m)-END_OF_PRTCL)>=1.00E-10_RK)
      m=m+1
    enddo
    m=m+1
  end function PC_ISInThisProc

  !**********************************************************************
  ! PC_reallocate_ghost_for_Cntct
  !**********************************************************************  
  subroutine PC_reallocate_ghost_for_Cntct(this,ng)
    implicit none
    class(Prtcl_Comm_info)::this
    integer,intent(in)::ng

    ! locals
    integer::sizep,sizen,ierrTmp,ierror=0
    integer,dimension(:),allocatable:: IntVec
    type(real3),dimension(:),allocatable::Real3Vec
    type(real4),dimension(:),allocatable::Real4Vec

    sizep= GPrtcl_list%mGhost_CS
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= min(sizen,DEM_Opt%numPrtcl)
    sizen= max(sizen,ng+1)
    GPrtcl_list%mGhost_CS=sizen  ! NOTE HERE, sometimes GPrtcl_list%mGhost_CS CAN bigger than DEM_Opt%numPrtcl

    ! ======= integer vector part =======
    call move_alloc(GhostP_id,IntVec)
    allocate(GhostP_id(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GhostP_id(1:sizep)=IntVec

    call move_alloc(GhostP_pType,IntVec)
    allocate(GhostP_pType(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GhostP_pType(1:sizep)=IntVec
    deallocate(IntVec)  

    ! ======= real3 vercor part =======
    call move_alloc(GhostP_linVel,Real3Vec)
    allocate(GhostP_linVel(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GhostP_linVel(1:sizep)=Real3Vec
 
    call move_alloc(GhostP_rotVel,Real3Vec)
    allocate(GhostP_rotVel(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GhostP_rotVel(1:sizep)=Real3Vec
    deallocate(Real3Vec)

    ! ======= real4 vercor part =======
    call move_alloc(GhostP_PosR,Real4Vec)
    allocate(GhostP_PosR(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GhostP_PosR(1:sizep)=Real4Vec
    deallocate(Real4Vec)

    if(ierror/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort," PC_reallocate_ghost_for_Cntct"," Reallocate wrong!")
      call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif   
    !call DEMLogInfo%CheckForError(ErrT_Pass," reallocate_ghost_for_Cntct"," Need to reallocate Ghost variables")
    !call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    !call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    !call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine PC_reallocate_ghost_for_Cntct

  !**********************************************************************
  ! PC_reallocate_sendlist
  !**********************************************************************  
  subroutine PC_reallocate_sendlist(this,ns)
    implicit none
    class(Prtcl_Comm_info)::this
    integer,intent(in)::ns

    ! locals
    integer::sizep,sizen,ierror
    integer,dimension(:),allocatable:: IntVec

    sizep=this%msend
    sizen=int(1.2_RK*real(sizep,kind=RK))
    sizen=min(sizen,DEM_Opt%numPrtcl)
    sizen=max(sizen,ns+1)
    this%msend=sizen
   
    call move_alloc(sendlist, IntVec)
    allocate(sendlist(sizen),stat=ierror)
    sendlist(1:sizep)=IntVec
    deallocate(IntVec)

    if(ierror/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort," PC_reallocate_sendlist"," Reallocate wrong!")
      call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif   
    !call DEMLogInfo%CheckForError(ErrT_Pass," reallocate_sendlist"," Need to reallocate sendlist")
    !call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    !call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    !call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine PC_reallocate_sendlist

  !**********************************************************************
  ! Prtcl_Comm_For_Cntct_fixed
  !**********************************************************************
  subroutine PC_Comm_For_Cntct_fixed(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    real(RK)::px,py,pz
    integer,dimension(:),allocatable:: IntVec
    type(real4),dimension(:),allocatable::Real4Vec
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer::i,ierror,request(4),SRstatus(MPI_STATUS_SIZE)
    integer::nsend,nsend2,nrecv,nrecv2,nsendg,ng,ngp,ngpp,GhostFixedCS_size
    
    ! (id 1) +(ptype 1) +(PosR 3) =5
    GhostFixedCS_size   = 5

    ng=0
    SELECT CASE(DEM_decomp%Prtcl_Pencil)
    CASE(x_axis)     ! ccccccccccccccccccccccccc  x-axis  ccccccccccccccccccccccccccc

      ! step1: Handle x-dir
      IF(pbc(1)) THEN
        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px+simLen%x
          endif
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px-simLen%x
          endif
        enddo
      ENDIF

      ! step2: send to yp_axis, and receive from ym_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py>=yed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py-simLen%y
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py-simLen%y
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py >=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py >=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,yp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to ym_axis, and receive from yp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py+simLen%y
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py+simLen%y
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py <=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py <=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,ym_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz>=zed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz-simLen%z
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz-simLen%z
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,zp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz+simLen%z
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz+simLen%z
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,zm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(y_axis)     ! ccccccccccccccccccccccccc  y-axis  ccccccccccccccccccccccccccc

      ! step1: Handle y-dir
      IF(pbc(2)) THEN
        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py+simLen%y
          endif
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py-simLen%y
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px>=xed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px-simLen%x
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px-simLen%x
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,xp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px+simLen%x
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px+simLen%x
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,xm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to zp_axis, and receive from zm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz>=zed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz-simLen%z
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz-simLen%z
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz>=zed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,zp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to zm_axis, and receive from zp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz+simLen%z
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz+simLen%z
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          pz=GhostPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,zm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

    CASE(z_axis)     ! ccccccccccccccccccccccccc  z-axis  ccccccccccccccccccccccccccc

      ! step1: Handle z-dir
      IF(pbc(3)) THEN
        do i=1,GPrtcl_list%mlocalFix
          pz=GPFix_PosR(i)%z
          if(pz<=zst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz+simLen%z
          endif
          if(pz>=zed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%z = pz-simLen%z
          endif
        enddo
      ENDIF

      ! step2: send to xp_axis, and receive from xm_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px>=xed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px-simLen%x
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px>=xed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px-simLen%x
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px >=xed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,xp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(1),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step3: send to xm_axis, and receive from xp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px+simLen%x
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px<=xst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%x = px+simLen%x
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          px=GhostPFix_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend
        do i=1,GPrtcl_list%mlocalFix
          px=GPFix_PosR(i)%x
          if(px <=xst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,xm_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(2),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step4: send to yp_axis, and receive from ym_dir
      nsend=0; nrecv=0; ngp=ng; ngpp=ng
      IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py>=yed0_cs ) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py-simLen%y
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py>=yed0_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py-simLen%y
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py>=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py>=yed0_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,yp_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(3),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>=GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)

      ! step5: send to ym_axis, and receive from yp_dir
      nsend=0; nrecv=0; ngp=ng; !ngpp=ng
      IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GhostPFix_id(i)
            GhostPFix_pType(ng)  = GhostPFix_pType(i)
            GhostPFix_PosR(ng)   = GhostPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py+simLen%y
          endif
        enddo

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py<=yst2_cs) then
            ng=ng+1
            if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
            GhostPFix_id(ng)     = GPFix_id(i)
            GhostPFix_pType(ng)  = GPFix_ptype(i)
            GhostPFix_PosR(ng)   = GPFix_PosR(i)
            GhostPFix_PosR(ng)%y = py+simLen%y
          endif
        enddo  
      ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
        
        do i=1,ngpp    ! consider the previous ghost particle firstly
          py=GhostPFix_PosR(i)%y
          if(py<=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
        nsendg=nsend

        do i=1,GPrtcl_list%mlocalFix
          py=GPFix_PosR(i)%y
          if(py<=yst2_cs) then
            nsend=nsend+1
            if(nsend > this%msend) call this%reallocate_sendlist(nsend)
            sendlist(nsend)=i
          endif
        enddo
      ENDIF
      call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                        nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
      if(nrecv>0) then
        nrecv2=nrecv*GhostFixedCS_size
        allocate(buf_recv(nrecv2))
        call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
      endif
      if(nsend>0) then
        nsend2=nsend*GhostFixedCS_size
        allocate(buf_send(nsend2))
        call pack_cntct_fixed(buf_send,nsendg,nsend,ym_axis)
        call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
      endif
      if(nrecv>0) then
        call MPI_WAIT(request(4),SRstatus,ierror)
        ng=ng+nrecv
        if(ng>GPrtcl_list%nGhostFix_CS) call reallocate_ghostFix_for_Cntct(ng)
        call unpack_cntct_fixed(buf_recv,ngp+1,ng)
      endif
      if(allocated(buf_send))deallocate(buf_send) 
      if(allocated(buf_recv))deallocate(buf_recv)
    END SELECT

    if(ng>0) then  ! 2021-02-11, Gong Zheng
      if(GPrtcl_list%mlocalFix>0) then
        call move_alloc(GPFix_id,IntVec)
        allocate(GPFix_id(GPrtcl_list%mlocalFix +ng))
        GPFix_id(1:GPrtcl_list%mlocalFix)= IntVec
        GPFix_id(GPrtcl_list%mlocalFix+1:GPrtcl_list%mlocalFix +ng)   = GhostPFix_id(1:ng)

        call move_alloc(GPFix_pType,IntVec)
        allocate(GPFix_pType(GPrtcl_list%mlocalFix +ng))
        GPFix_pType(1:GPrtcl_list%mlocalFix)= IntVec
        GPFix_pType(GPrtcl_list%mlocalFix+1:GPrtcl_list%mlocalFix +ng)= GhostPFix_pType(1:ng)
        deallocate(IntVec)

        call move_alloc(GPFix_PosR,Real4Vec)
        allocate(GPFix_PosR(GPrtcl_list%mlocalFix +ng))
        GPFix_PosR(1:GPrtcl_list%mlocalFix)= Real4Vec
        GPFix_PosR(GPrtcl_list%mlocalFix+1:GPrtcl_list%mlocalFix +ng) = GhostPFix_PosR(1:ng)
        deallocate(Real4Vec)
      else
        if(allocated(GPFix_id))   deallocate(GPFix_id);   allocate(GPFix_id(ng));   GPFix_id= GhostPFix_id(1:ng)
        if(allocated(GPFix_PosR)) deallocate(GPFix_PosR); allocate(GPFix_PosR(ng)); GPFix_PosR=GhostPFix_PosR(1:ng)
        if(allocated(GPFix_pType))deallocate(GPFix_pType);allocate(GPFix_pType(ng));GPFix_pType=GhostPFix_pType(1:ng)      
      endif
    endif
    deallocate(GhostPFix_id, GhostPFix_pType, GhostPFix_PosR)
    GPrtcl_list%nGhostFix_CS = ng
  end subroutine PC_Comm_For_Cntct_fixed

  !**********************************************************************
  ! pack_cntct_fixed
  !**********************************************************************
  subroutine pack_cntct_fixed(buf_send,nsendg,nsend,dir)
    implicit none
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsendg,nsend,dir

    ! locals
    integer::i,id,m
    real(RK)::dx,dy,dz

    m=1
    dx=dx_pbc(dir)
    dy=dy_pbc(dir)
    dz=dz_pbc(dir)
    do i=1,nsendg
      id=sendlist(i)
      buf_send(m)=real(GhostPFix_id(id));    m=m+1 ! 01
      buf_send(m)=real(GhostPFix_pType(id)); m=m+1 ! 02
      buf_send(m)=GhostPFix_PosR(id)%x+dx;   m=m+1 ! 03
      buf_send(m)=GhostPFix_PosR(id)%y+dy;   m=m+1 ! 04
      buf_send(m)=GhostPFix_PosR(id)%z+dz;   m=m+1 ! 05
    enddo
    do i=nsendg+1,nsend
      id=sendlist(i)
      buf_send(m)=real(GPFix_id(id));    m=m+1 ! 01
      buf_send(m)=real(GPFix_pType(id)); m=m+1 ! 02
      buf_send(m)=GPFix_PosR(id)%x+dx;   m=m+1 ! 03
      buf_send(m)=GPFix_PosR(id)%y+dy;   m=m+1 ! 04
      buf_send(m)=GPFix_PosR(id)%z+dz;   m=m+1 ! 05
    enddo
  end subroutine pack_cntct_fixed

  !**********************************************************************
  ! unpack_cntct_fixed
  !**********************************************************************
  subroutine unpack_cntct_fixed(buf_recv,n1,n2)
    implicit none
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::n1,n2

    ! locals
    integer::i,m,itype
   
    m=1
    do i=n1,n2
      GhostPFix_id(i)       =int(buf_recv(m)+0.2); m=m+1 ! 01 
      itype                 =int(buf_recv(m)+0.2); m=m+1 ! 02
      GhostPFix_PosR(i)%x   =buf_recv(m);          m=m+1 ! 03
      GhostPFix_PosR(i)%y   =buf_recv(m);          m=m+1 ! 04
      GhostPFix_PosR(i)%z   =buf_recv(m);          m=m+1 ! 05
      GhostPFix_pType(i)    =itype;                
      GhostPFix_PosR(i)%w   =DEMProperty%Prtcl_PureProp(itype)%Radius 
    enddo
  end subroutine unpack_cntct_fixed

  !**********************************************************************
  ! reallocate_ghostFix_for_Cntct
  !**********************************************************************  
  subroutine reallocate_ghostFix_for_Cntct(ng)
    implicit none
    integer,intent(in)::ng

    ! locals
    integer:: sizep,sizen
    integer,dimension(:),allocatable:: IntVec
    type(real4),dimension(:),allocatable::Real4Vec

    sizep= GPrtcl_list%nGhostFix_CS
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= min(sizen,DEM_Opt%numPrtclFix)
    sizen= max(sizen,ng+1)
    GPrtcl_list%nGhostFix_CS=sizen

    ! ======= integer vector part =======
    call move_alloc(GhostPFix_id,IntVec)
    allocate(GhostPFix_id(sizen))
    GhostPFix_id(1:sizep)=IntVec

    call move_alloc(GhostPFix_pType,IntVec)
    allocate(GhostPFix_pType(sizen))
    GhostPFix_pType(1:sizep)=IntVec
    deallocate(IntVec)  

    ! ======= real4 verctor part =======
    call move_alloc(GhostPFix_PosR,Real4Vec)
    allocate(GhostPFix_PosR(sizen))
    GhostPFix_PosR(1:sizep)=Real4Vec
    deallocate(Real4Vec)
  end subroutine reallocate_ghostFix_for_Cntct

#undef xm_axis
#undef xp_axis
#undef ym_axis
#undef yp_axis
#undef zm_axis
#undef zp_axis
end module Prtcl_Comm
