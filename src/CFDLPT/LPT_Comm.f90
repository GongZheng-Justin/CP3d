module LPT_Comm
  use MPI
  use m_TypeDef
  use m_LogInfo
  use LPT_Property
  use LPT_Decomp_2d
  use LPT_Variables
  use LPT_Parameters
  use LPT_ContactSearchPW
  use m_Decomp2d,only: nrank
  implicit none
  private

  logical,dimension(3)::pbc
  type(real3)::simLen
  integer,parameter::xm_axis=1
  integer,parameter::xp_axis=2
  integer,parameter::ym_axis=3
  integer,parameter::yp_axis=4
  integer,parameter::zm_axis=5
  integer,parameter::zp_axis=6
  real(RK),dimension(6)::dx_pbc
  real(RK),dimension(6)::dy_pbc
  real(RK),dimension(6)::dz_pbc  

  real(RK):: xst1_cs
  real(RK):: xed1_cs
  real(RK):: yst1_cs
  real(RK):: yed1_cs
  real(RK):: zst1_cs
  real(RK):: zed1_cs
  integer,allocatable,dimension(:)::sendlist

  type Prtcl_Comm_info
    integer :: msend
    integer :: Prtcl_Exchange_size
  contains
    procedure:: InitComm            => PC_InitComm
    procedure:: Comm_For_Exchange   => PC_Comm_For_Exchange
    procedure:: pack_Exchange       => PC_pack_Exchange
    procedure:: unpack_Exchange     => PC_unpack_Exchange
    procedure:: ISInThisProc        => PC_ISInThisProc
    procedure:: reallocate_sendlist => PC_reallocate_sendlist
  end type Prtcl_Comm_info
  type(Prtcl_Comm_info),public::LPTComm

contains

  !**********************************************************************
  ! PC_InitComm
  !**********************************************************************
  subroutine PC_InitComm(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    integer::iErr01
        
    pbc=LPT_Opt%IsPeriodic
    simLen = LPT_Opt%SimDomain_max-LPT_Opt%SimDomain_min

    xst1_cs = LPT_decomp%xSt   
    xed1_cs = LPT_decomp%xEd
    yst1_cs = LPT_decomp%ySt
    yed1_cs = LPT_decomp%yEd 
    zst1_cs = LPT_decomp%zSt
    zed1_cs = LPT_decomp%zEd

    dx_pbc=0.0_RK; dy_pbc=0.0_RK; dz_pbc=0.0_RK
    if(pbc(1)) then
      if(LPT_decomp%coord1==0)                 dx_pbc(xm_axis)= simLen%x
      if(LPT_decomp%coord1==LPT_decomp%prow-1) dx_pbc(xp_axis)=-simLen%x 
    endif
    if(pbc(2)) then
      dy_pbc(ym_axis)= simLen%y
      dy_pbc(yp_axis)=-simLen%y
    endif
    if(pbc(3)) then
      if(LPT_decomp%coord2==0)                 dz_pbc(zm_axis)= simLen%z
      if(LPT_decomp%coord2==LPT_decomp%pcol-1) dz_pbc(zp_axis)=-simLen%z
    endif

    ! (id 1) +(ptype 1) +(Mark 1) +(PosR   3) +(linvel 3*tsize) +(linAcc 3*tsize) +(GPrtcl_FpForce 3) 
    ! = 9+6*tsize
    this%Prtcl_Exchange_size = 9 + 6*GPrtcl_list%tsize
    this%msend=GPrtcl_list%mlocal
    allocate(sendlist(this%msend),   Stat=iErr01)
    if(iErr01/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"PC_InitComm","Allocation failed2")

  end subroutine PC_InitComm

  !**********************************************************************
  ! PC_Comm_For_Exchange
  !**********************************************************************
  subroutine PC_Comm_For_Exchange(this)
    implicit none
    class(Prtcl_Comm_info)::this

    ! locals
    integer::i,ierror,request(4),nlocal,nlocalp
    integer::nsend(2),nrecv(2)
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer,dimension(MPI_STATUS_SIZE) :: SRstatus
    real(RK)::px,py,pz

    nlocal=GPrtcl_list%nlocal

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
          nlocal = nlocal -1
        else
          i = i + 1
        endif
      enddo
    ENDIF

    ! step2: send to xp_axis, and receive from xm_dir
    nsend =0; nrecv =0
    IF(LPT_decomp%ProcNgh(3)==MPI_PROC_NULL) THEN
      i=1
      do while(i<=nlocal)
        px=GPrtcl_PosR(i)%x
        if(px >= xed1_cs) then
          call GPrtcl_list%copy(i,nlocal)
          nlocal = nlocal -1
        else
          i = i + 1
        endif
      enddo

    ELSEIF(LPT_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
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
          nsend(2)= nsend(2) + this%Prtcl_Exchange_size
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 2, int_type, LPT_decomp%ProcNgh(3), 9, &
                      nrecv, 2, int_type, LPT_decomp%ProcNgh(4), 9, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv(2)>0) then
      allocate(buf_recv(nrecv(2)))
      call MPI_IRECV(buf_recv,nrecv(2),real_type,LPT_decomp%ProcNgh(4),10,MPI_COMM_WORLD,request(1),ierror)
    endif
    if(nsend(2)>0) then
      allocate(buf_send(nsend(2)))
      nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
      call this%pack_Exchange(buf_send,nsend(1),nlocalp,xp_axis)
      call MPI_SEND(buf_send,nsend(2),real_type,LPT_decomp%ProcNgh(3),10,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv(2)>0) then
      call MPI_WAIT(request(1),SRstatus,ierror)
      if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
        call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
      endif
      call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xp_axis)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! step3: send to xm_axis, and receive from xp_dir
    nsend =0; nrecv =0
    IF(LPT_decomp%ProcNgh(4)==MPI_PROC_NULL) THEN
      i=1
      do while(i<=nlocal)
        px=GPrtcl_PosR(i)%x
        if(px < xst1_cs) then
          call GPrtcl_list%copy(i,nlocal)
          nlocal = nlocal -1
        else
          i = i + 1
        endif
      enddo

    ELSEIF(LPT_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
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
          nsend(2)= nsend(2)+ this%Prtcl_Exchange_size
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 2, int_type, LPT_decomp%ProcNgh(4),11, &
                      nrecv, 2, int_type, LPT_decomp%ProcNgh(3),11, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv(2)>0) then
      allocate(buf_recv(nrecv(2)))
      call MPI_IRECV(buf_recv,nrecv(2),real_type,LPT_decomp%ProcNgh(3),12,MPI_COMM_WORLD,request(2),ierror)
    endif
    if(nsend(2)>0) then
      allocate(buf_send(nsend(2)))
      nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
      call this%pack_Exchange(buf_send,nsend(1),nlocalp,xm_axis)
      call MPI_SEND(buf_send,nsend(2),real_type,LPT_decomp%ProcNgh(4),12,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv(2)>0) then
      call MPI_WAIT(request(2),SRstatus,ierror)
      if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
        call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
      endif
      call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,xm_axis)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! step4: send to zp_axis, and receive from zm_dir
    nsend =0; nrecv =0
    IF(LPT_decomp%ProcNgh(1)==MPI_PROC_NULL) THEN
      i=1
      do while(i<=nlocal)
        pz=GPrtcl_PosR(i)%z
        if(pz >= zed1_cs) then
          call GPrtcl_list%copy(i,nlocal)
          nlocal = nlocal -1
        else
          i = i + 1
        endif
      enddo

    ELSEIF(LPT_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
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
          nsend(2)= nsend(2)+ this%Prtcl_Exchange_size
        endif
      enddo

    ENDIF
    call MPI_SENDRECV(nsend, 2, int_type, LPT_decomp%ProcNgh(1),13, &
                      nrecv, 2, int_type, LPT_decomp%ProcNgh(2),13, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv(2)>0) then
      allocate(buf_recv(nrecv(2)))
      call MPI_IRECV(buf_recv,nrecv(2),real_type,LPT_decomp%ProcNgh(2),14,MPI_COMM_WORLD,request(3),ierror)
    endif
    if(nsend(2)>0) then
      allocate(buf_send(nsend(2)))
      nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
      call this%pack_Exchange(buf_send,nsend(1),nlocalp,zp_axis)
      call MPI_SEND(buf_send,nsend(2),real_type,LPT_decomp%ProcNgh(1),14,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv(2)>0) then
      call MPI_WAIT(request(3),SRstatus,ierror)
      if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
        call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
      endif
      call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zp_axis)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! step5: send to zm_axis, and receive from zp_dir
    nsend =0; nrecv =0
    IF(LPT_decomp%ProcNgh(2)==MPI_PROC_NULL) THEN
      i=1
      do while(i<=nlocal)
        pz=GPrtcl_PosR(i)%z
        if(pz < zst1_cs) then
          call GPrtcl_list%copy(i,nlocal)
          nlocal = nlocal -1
        else
          i = i + 1
        endif
      enddo

    ELSEIF(LPT_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
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
          nsend(2)= nsend(2)+ this%Prtcl_Exchange_size
        endif
      enddo

    ENDIF
    call MPI_SENDRECV(nsend, 2, int_type, LPT_decomp%ProcNgh(2),15, &
                      nrecv, 2, int_type, LPT_decomp%ProcNgh(1),15, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv(2)>0) then
      allocate(buf_recv(nrecv(2)))
      call MPI_IRECV(buf_recv,nrecv(2),real_type,LPT_decomp%ProcNgh(1),16,MPI_COMM_WORLD,request(4),ierror)
    endif
    if(nsend(2)>0) then
      allocate(buf_send(nsend(2)))
      nlocalp= nlocal;  nlocal = nlocalp - nsend(1)
      call this%pack_Exchange(buf_send,nsend(1),nlocalp,zm_axis)
      call MPI_SEND(buf_send,nsend(2),real_type,LPT_decomp%ProcNgh(2),16,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv(2)>0) then
      call MPI_WAIT(request(4),SRstatus,ierror)
      if(nlocal + nrecv(1)>=GPrtcl_list%mlocal) then
        call GPrtcl_list%ReallocatePrtclVar(nlocal + nrecv(1))
      endif

      call this%unpack_Exchange(buf_recv,nrecv(1),nlocal,zm_axis)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)
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
    integer::i,j,id,m,nlocal
    real(RK)::dx,dy,dz

    m=1
    nlocal = nlocalp
    dx=dx_pbc(dir)
    dy=dy_pbc(dir)
    dz=dz_pbc(dir)
   
    DO i=1,nsend
      id = sendlist(i)
      buf_send(m)=real(GPrtcl_id(id));      m=m+1 ! 01
      buf_send(m)=real(GPrtcl_pType(id));   m=m+1 ! 02
      buf_send(m)=real(GPrtcl_usrMark(id)); m=m+1 ! 03
      buf_send(m)=GPrtcl_PosR(id)%x+dx;     m=m+1 ! 04
      buf_send(m)=GPrtcl_PosR(id)%y+dy;     m=m+1 ! 05
      buf_send(m)=GPrtcl_PosR(id)%z+dz;     m=m+1 ! 06
      do j=1,GPrtcl_list%tsize
        buf_send(m)=GPrtcl_linVel(j,id)%x;  m=m+1 ! 6* tsize
        buf_send(m)=GPrtcl_linVel(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_linVel(j,id)%z;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%x;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%y;  m=m+1 ! 
        buf_send(m)=GPrtcl_linAcc(j,id)%z;  m=m+1 ! 
      enddo
      buf_send(m)=GPrtcl_FpForce(id)%x;     m=m+1
      buf_send(m)=GPrtcl_FpForce(id)%y;     m=m+1
      buf_send(m)=GPrtcl_FpForce(id)%z;     m=m+1
    ENDDO
    DO i=nsend,1,-1
      id = sendlist(i)
      call GPrtcl_list%copy(id,nlocal)
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
    DO i=1,nrecv
      if( .not.(this%ISInThisProc(buf_recv,m,dir)) ) cycle
      GPrtcl_id(id)     = nint(buf_recv(m)); m=m+1 ! 01
      itype             = nint(buf_recv(m)); m=m+1 ! 02
      GPrtcl_usrMark(id)= nint(buf_recv(m)); m=m+1 ! 03
      GPrtcl_PosR(id)%x = buf_recv(m);       m=m+1 ! 04
      GPrtcl_PosR(id)%y = buf_recv(m);       m=m+1 ! 05
      GPrtcl_PosR(id)%z = buf_recv(m);       m=m+1 ! 06
      do j=1,GPrtcl_list%tsize
        GPrtcl_linVel(j,id)%x = buf_recv(m); m=m+1 ! 6* tsize
        GPrtcl_linVel(j,id)%y = buf_recv(m); m=m+1
        GPrtcl_linVel(j,id)%z = buf_recv(m); m=m+1
        GPrtcl_linAcc(j,id)%x = buf_recv(m); m=m+1
        GPrtcl_linAcc(j,id)%y = buf_recv(m); m=m+1
        GPrtcl_linAcc(j,id)%z = buf_recv(m); m=m+1
      enddo
      GPrtcl_FpForce(id)%x    = buf_recv(m); m=m+1
      GPrtcl_FpForce(id)%y    = buf_recv(m); m=m+1
      GPrtcl_FpForce(id)%z    = buf_recv(m); m=m+1
      GPrtcl_pType(id)=itype
      GPrtcl_PosR(id)%w =LPTProperty%Prtcl_PureProp(itype)%Radius
 
      id = id + 1
      nlocal =nlocal + 1
    ENDDO
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
    
    call LPTLogInfo%CheckForError(ErrT_Pass," PC_ISInThisProc"," The following particle is deleted: ")
    call LPTLogInfo%OutInfo(" Exchange direction  is :"//trim(num2str(dir)  ),3)
    call LPTLogInfo%OutInfo("   The particle id   is :"//trim(num2str( nint(buf_recv(m)  ))),3)
    call LPTLogInfo%OutInfo("       particle type is :"//trim(num2str( nint(buf_recv(m+1)))),3)
    call LPTLogInfo%OutInfo("       x-coordinate  is :"//trim(num2str( buf_recv(m+3) )), 3)
    call LPTLogInfo%OutInfo("       y-coordinate  is :"//trim(num2str( buf_recv(m+4) )), 3)
    call LPTLogInfo%OutInfo("       z-coordinate  is :"//trim(num2str( buf_recv(m+5) )), 3)
    call LPTLogInfo%OutInfo(" Present processor   is :"//trim(num2str(nrank)), 3)

    m = m + this%Prtcl_Exchange_size
    m=m+1
  end function PC_ISInThisProc

  !**********************************************************************
  ! PC_reallocate_sendlist
  !**********************************************************************  
  subroutine PC_reallocate_sendlist(this,ns)
    implicit none
    class(Prtcl_Comm_info)::this
    integer,intent(in)::ns

    ! locals
    integer:: sizep,sizen
    integer,dimension(:),allocatable:: IntVec

    sizep= this%msend
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= min(sizen,LPT_Opt%numPrtcl)
    sizen=max(sizen,ns+1)
    this%msend=sizen
   
    call move_alloc(sendlist, IntVec)
    allocate(sendlist(sizen))
    sendlist(1:sizep)=IntVec
    deallocate(IntVec)

    call LPTLogInfo%CheckForError(ErrT_Pass," reallocate_sendlist"," Need to reallocate sendlist")
    call LPTLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    call LPTLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    call LPTLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine PC_reallocate_sendlist

end module LPT_Comm
