#ifndef IBMDistribute2
  !******************************************************************
  ! delta_fun
  !******************************************************************
  real(RK) function delta_fun(ratio_in)
    implicit none
    real(RK),intent(in)::ratio_in

    ! locals
    real(RK)::ratio

    ! Three-point type discrete Delta Function
    ! [1] A.M. Roma, Ch. S. Peskin, M.J. Berger,  J. Comput. Phys. 153 (1999) . 
    !       An adaptive version of the immersed boundary method.
    ! [2] T. Kempe, J. Fröhlich, J. Comput. Phys. 231 (2012) 
    !       An improved immersed boundary method with direct forcing for the simulation of particle laden flows.
    ratio= abs(ratio_in)
    if(ratio> 1.50_RK) then
      delta_fun= 0.0_RK
    elseif(ratio> 0.50_RK) then
      delta_fun= 0.166666666666666666667_RK*(5.0_RK-3.0_RK*ratio-sqrt(1.0_RK-3.0_RK*(1.0_RK-ratio)**2))
    else
      delta_fun= 0.333333333333333333333_RK*(1.0_RK+sqrt(1.0_RK-3.0_RK*ratio*ratio))
    endif
  end function delta_fun
#endif

  !******************************************************************
  ! PComm_FluidIndicator
  !******************************************************************
  subroutine PComm_FluidIndicator()
    implicit none

    ! locals
    real(RK)::px,pxt,pz,pzt
    integer,dimension(MPI_STATUS_SIZE)::SRstatus
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer::i,ierror,request(4),GFluidIndicator_size_Comm
    integer::nlocalNow,nsend,nsend2,nrecv,nrecv2,nsendg,ng,ngp,ngpp    

    ! This part is similar to what is done in subroutine PC_Comm_For_Cntct of file "ACM_comm.f90"(only y_pencil)
    ng=0
    GFluidIndicator_size_Comm=4 ! ptype + pos
    nlocalNow= GPrtcl_list%nlocal

    ! Step1: send to xp_axis, and receive from xm_dir
    nsend=0; nrecv=0; ngp=ng; ngpp=ng
    IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,nlocalNow
        px = GPrtcl_PosR(i)%x
        pxt= px+ GPrtcl_PosR(i)%w
        if(pxt +SMALL> xedCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_pType(ng)    = GPrtcl_ptype(i)
          GhostP_PosR(ng)     = GPrtcl_PosR(i)
          GhostP_PosR(ng)%x   = px-xlx
        endif
      enddo  
    ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) THEN
      nsendg=nsend
      do i=1,nlocalNow
        pxt= GPrtcl_PosR(i)%x +GPrtcl_PosR(i)%w
        if(pxt+SMALL >xedCoord) then
          nsend=nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 1, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GFluidIndicator_size_comm
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GFluidIndicator_size_comm
      allocate(buf_send(nsend2))
      call pack_Comm_FluidIndicator(buf_send,nsendg,nsend,xp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(1),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_Comm_FluidIndicator(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step2: send to xm_axis, and receive from xp_dir
    nsend=0; nrecv=0; ngp=ng; !ngpp=ng  
    IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,nlocalNow
        px = GPrtcl_PosR(i)%x
        pxt= px- GPrtcl_PosR(i)%w
        if(pxt-SMALL<xstCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_pType(ng)  = GPrtcl_ptype(i)
          GhostP_PosR(ng)   = GPrtcl_PosR(i)
          GhostP_PosR(ng)%x = px+xlx
        endif
      enddo 
    ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) THEN
      nsendg=nsend
      do i=1,nlocalNow
        pxt= GPrtcl_PosR(i)%x- GPrtcl_PosR(i)%w
        if(pxt-SMALL<xstCoord) then
          nsend=nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 3, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GFluidIndicator_size_comm
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GFluidIndicator_size_comm
      allocate(buf_send(nsend2))
      call pack_Comm_FluidIndicator(buf_send,nsendg,nsend,xm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(2),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_Comm_FluidIndicator(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step3: send to zp_axis, and receive from zm_dir
    nsend=0; nrecv=0; ngp=ng; ngpp=ng
    IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,ngpp    ! consider the previous ghost particle firstly
        pz = GhostP_PosR(i)%z 
        pzt= pz+ GhostP_PosR(i)%w
        if(pzt+SMALL>zedCoord ) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_pType(ng)  = GPrtcl_ptype(i)
          GhostP_PosR(ng)   = GhostP_PosR(i)
          GhostP_PosR(ng)%z = pz-zlz
        endif
      enddo

      do i=1,nlocalNow
        pz = GPrtcl_PosR(i)%z 
        pzt= pz+ GPrtcl_PosR(i)%w
        if(pzt+SMALL>zedCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_pType(ng)  = GPrtcl_ptype(i)
          GhostP_PosR(ng)   = GPrtcl_PosR(i)
          GhostP_PosR(ng)%z = pz-zlz
        endif
      enddo  
    ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) THEN
      do i=1,ngpp    ! consider the previous ghost particle firstly
        pzt= GhostP_PosR(i)%z+ GhostP_PosR(i)%w
        if(pzt+SMALL>zedCoord) then
          nsend=nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
      nsendg=nsend

      do i=1,nlocalNow
        pzt=GPrtcl_PosR(i)%z+ GPrtcl_PosR(i)%w
        if(pzt+SMALL>zedCoord) then
          nsend=nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 5, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GFluidIndicator_size_comm
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GFluidIndicator_size_comm
      allocate(buf_send(nsend2))
      call pack_Comm_FluidIndicator(buf_send,nsendg,nsend,zp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(3),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_Comm_FluidIndicator(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step4: send to zm_axis, and receive from zp_dir
    nsend=0; nrecv=0; ngp=ng; !ngpp=ng
    IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,ngpp    ! consider the previous ghost particle firstly
        pz = GhostP_PosR(i)%z
        pzt= pz- GhostP_PosR(i)%w
        if(pzt-SMALL<zstCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_pType(ng)  = GPrtcl_ptype(i)
          GhostP_PosR(ng)   = GhostP_PosR(i)
          GhostP_PosR(ng)%z = pz+ zlz
        endif
      enddo

      do i=1,nlocalNow
        pz = GPrtcl_PosR(i)%z 
        pzt= pz- GPrtcl_PosR(i)%w
        if(pzt-SMALL<zstCoord ) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_pType(ng)  = GPrtcl_ptype(i)
          GhostP_PosR(ng)   = GPrtcl_PosR(i)
          GhostP_PosR(ng)%z = pz+ zlz
        endif
      enddo  
    ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) then
      do i=1,ngpp    ! consider the previous ghost particle firstly
        pz = GhostP_PosR(i)%z
        pzt= pz- GhostP_PosR(i)%w
        if(pzt-SMALL<zstCoord) then
          nsend=nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
      nsendg=nsend

      do i=1,nlocalNow
        pzt= GPrtcl_PosR(i)%z- GPrtcl_PosR(i)%w
        if(pzt-SMALL<zstCoord) then
          nsend=nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 7, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GFluidIndicator_size_comm
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GFluidIndicator_size_comm
      allocate(buf_send(nsend2))
      call pack_Comm_FluidIndicator(buf_send,nsendg,nsend,zm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(4),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_Comm_FluidIndicator(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    nGhostIBM= ng
  end subroutine PComm_FluidIndicator

  !**********************************************************************
  ! pack_Comm_FluidIndicator
  !**********************************************************************
  subroutine pack_Comm_FluidIndicator(buf_send,nsendg,nsend,dir)
    implicit none
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsendg,nsend,dir

    ! locals
    integer::i,id,m
    real(RK)::pdx,pdz

    m=1
    pdx=dx_pbc(dir)
    pdz=dz_pbc(dir)
    do i=1,nsendg
      id=sendlist(i)
      buf_send(m)=real(GhostP_pType(id)); m=m+1 ! 01
      buf_send(m)=GhostP_PosR(id)%x+pdx;  m=m+1 ! 02
      buf_send(m)=GhostP_PosR(id)%y;      m=m+1 ! 03
      buf_send(m)=GhostP_PosR(id)%z+pdz;  m=m+1 ! 04
    enddo
    do i=nsendg+1,nsend
      id=sendlist(i)
      buf_send(m)=real(GPrtcl_pType(id)); m=m+1 ! 01
      buf_send(m)=GPrtcl_PosR(id)%x+pdx;  m=m+1 ! 02
      buf_send(m)=GPrtcl_PosR(id)%y;      m=m+1 ! 03
      buf_send(m)=GPrtcl_PosR(id)%z+pdz;  m=m+1 ! 04
    enddo
  end subroutine pack_Comm_FluidIndicator

  !**********************************************************************
  ! unpack_Comm_FluidIndicator
  !**********************************************************************
  subroutine unpack_Comm_FluidIndicator(buf_recv,n1,n2)
    implicit none
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::n1,n2

    ! locals
    integer::i,m,itype
   
    m=1
    do i=n1,n2
      itype              =nint(buf_recv(m)); m=m+1 ! 01
      GhostP_PosR(i)%x   =buf_recv(m);       m=m+1 ! 02
      GhostP_PosR(i)%y   =buf_recv(m);       m=m+1 ! 03
      GhostP_PosR(i)%z   =buf_recv(m);       m=m+1 ! 04
      GhostP_pType(i)    =itype;                
      GhostP_PosR(i)%w   =DEMProperty%Prtcl_PureProp(itype)%Radius 
    enddo
  end subroutine unpack_Comm_FluidIndicator

  !******************************************************************
  ! Init_IBPFix
  !******************************************************************
  subroutine Init_IBPFix()
    implicit none

    ! locals
    type(real3):: Cposition,LagrangeP,PosDiff
    integer:: i,j,itype,nPartition,iErr01,iErr02,iErr03,iErr04,iErrSum
    integer:: ic,jc,kc,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    ! [1] Here ONLY the IBP points within the Processor physical Domains will be considered.
    !       The IBP points outsides will be skipped. 
    nIBPFix = 0
    DO i=1,GPrtcl_list%mlocalFix+ GPrtcl_list%nGhostFix_CS
      itype = GPFix_pType(i)
      Cposition = GPFix_PosR(i)
      nPartition= PrtclIBMProp(itype)%nPartition
      do j=1,nPartition
        PosDiff= SpherePartitionCoord(j,itype)
        LagrangeP = Cposition+ PosDiff
        ic= floor(LagrangeP%x*rdx)+1
        jc= floor(LagrangeP%y*rdyUniform)+1
        kc= floor(LagrangeP%z*rdz)+1
        if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 &
        .or. (jc-y1start(2))*(jc-y1end(2))>0) cycle
        nIBPFix= nIBPFix+1
      enddo
    ENDDO
    if(nIBPFix==0) return

    allocate(IBPFix_indxyz(6,nIBPFix), Stat=iErr01)
    allocate(IBPFix_VolRatio(nIBPFix), Stat=iErr02)
    allocate(IBPFix_Pos(nIBPFix),      Stat=iErr03)
    allocate(IBPFix_Vel(nIBPFix),      Stat=iErr04)
    iErrSum=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)
    if(iErrSum/=0)   call MainLog%CheckForError(ErrT_Abort,"Init_IBPFix","Allocation failed")

    nIBPFix = 0
    DO i=1,GPrtcl_list%mlocalFix+ GPrtcl_list%nGhostFix_CS
      itype = GPFix_pType(i)
      Cposition = GPFix_PosR(i)
      nPartition= PrtclIBMProp(itype)%nPartition
      do j=1,nPartition
        PosDiff= SpherePartitionCoord(j,itype)
        LagrangeP = Cposition+ PosDiff
        ic= floor(LagrangeP%x*rdx)+1
        jc= floor(LagrangeP%y*rdyUniform)+1
        kc= floor(LagrangeP%z*rdz)+1
        if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 &
        .or. (jc-y1start(2))*(jc-y1end(2))>0) cycle
        nIBPFix= nIBPFix+1
        IBPFix_Pos(nIBPFix)= LagrangeP
        IBPFix_VolRatio(nIBPFix)= IbpVolRatio(itype)
#define   clc_Point_indxyz_IBPFix
#include "clc_Point_indxyz_inc.f90"
#undef    clc_Point_indxyz_IBPFix
      enddo
    ENDDO
  end subroutine Init_IBPFix

  !**********************************************************************
  ! clc_VoidFraction_IBM
  !**********************************************************************
  function clc_VoidFraction_IBM(CenterDiff,Radius) result(VoidFraction)
    implicit none
    type(real3),intent(in)::CenterDiff   
    real(RK),intent(in)::Radius
   
    ! locals
    real(RK)::VoidFraction,SumPhi,vertexDX,vertexDY,vertexDZ,CellPhi
    
    ! T. Kempe, J. Fröhlich, J. Comput. Phys. 231 (2012)
    VoidFraction= 0.0_RK; SumPhi=0.0_RK

    ! vertex (0,0,0)
    vertexDX= CenterDiff%x- dxhalf
    vertexDY= CenterDiff%y- dyhalf
    vertexDZ= CenterDiff%z- dzhalf
    CellPhi= sqrt(vertexDX*vertexDX +vertexDY*vertexDY +vertexDZ*vertexDZ)- Radius
    if(CellPhi<0.0_RK) then
      SumPhi= SumPhi- CellPhi;  VoidFraction= VoidFraction- CellPhi
    else
      SumPhi= SumPhi+ CellPhi
    endif

    ! vertex (0,0,1)   
    vertexDX= CenterDiff%x- dxhalf
    vertexDY= CenterDiff%y- dyhalf
    vertexDZ= CenterDiff%z+ dzhalf
    CellPhi= sqrt(vertexDX*vertexDX +vertexDY*vertexDY +vertexDZ*vertexDZ)- Radius
    if(CellPhi<0.0_RK) then
      SumPhi= SumPhi- CellPhi;  VoidFraction= VoidFraction- CellPhi
    else
      SumPhi= SumPhi+ CellPhi
    endif

    ! vertex (0,1,0)
    vertexDX= CenterDiff%x- dxhalf
    vertexDY= CenterDiff%y+ dyhalf
    vertexDZ= CenterDiff%z- dzhalf
    CellPhi= sqrt(vertexDX*vertexDX +vertexDY*vertexDY +vertexDZ*vertexDZ)- Radius
    if(CellPhi<0.0_RK) then
      SumPhi= SumPhi- CellPhi;  VoidFraction= VoidFraction- CellPhi
    else
      SumPhi= SumPhi+ CellPhi
    endif

    ! vertex (1,0,0)
    vertexDX= CenterDiff%x+ dxhalf
    vertexDY= CenterDiff%y- dyhalf
    vertexDZ= CenterDiff%z- dzhalf
    CellPhi= sqrt(vertexDX*vertexDX +vertexDY*vertexDY +vertexDZ*vertexDZ)- Radius
    if(CellPhi<0.0_RK) then
      SumPhi= SumPhi- CellPhi;  VoidFraction= VoidFraction- CellPhi
    else
      SumPhi= SumPhi+ CellPhi
    endif

    ! vertex (1,1,0)
    vertexDX= CenterDiff%x+ dxhalf
    vertexDY= CenterDiff%y+ dyhalf
    vertexDZ= CenterDiff%z- dzhalf
    CellPhi= sqrt(vertexDX*vertexDX +vertexDY*vertexDY +vertexDZ*vertexDZ)- Radius
    if(CellPhi<0.0_RK) then
      SumPhi= SumPhi- CellPhi;  VoidFraction= VoidFraction- CellPhi
    else
      SumPhi= SumPhi+ CellPhi
    endif

    ! vertex (1,0,1)
    vertexDX= CenterDiff%x+ dxhalf
    vertexDY= CenterDiff%y- dyhalf
    vertexDZ= CenterDiff%z+ dzhalf
    CellPhi= sqrt(vertexDX*vertexDX +vertexDY*vertexDY +vertexDZ*vertexDZ)- Radius
    if(CellPhi<0.0_RK) then
      SumPhi= SumPhi- CellPhi;  VoidFraction= VoidFraction- CellPhi
    else
      SumPhi= SumPhi+ CellPhi
    endif

    ! vertex (0,1,1)
    vertexDX= CenterDiff%x- dxhalf
    vertexDY= CenterDiff%y+ dyhalf
    vertexDZ= CenterDiff%z+ dzhalf
    CellPhi= sqrt(vertexDX*vertexDX +vertexDY*vertexDY +vertexDZ*vertexDZ)- Radius
    if(CellPhi<0.0_RK) then
      SumPhi= SumPhi- CellPhi;  VoidFraction= VoidFraction- CellPhi
    else
      SumPhi= SumPhi+ CellPhi
    endif

    ! vertex (1,1,1)
    vertexDX= CenterDiff%x+ dxhalf
    vertexDY= CenterDiff%y+ dyhalf
    vertexDZ= CenterDiff%z+ dzhalf
    CellPhi= sqrt(vertexDX*vertexDX +vertexDY*vertexDY +vertexDZ*vertexDZ)- Radius
    if(CellPhi<0.0_RK) then
      SumPhi= SumPhi- CellPhi;  VoidFraction= VoidFraction- CellPhi
    else
      SumPhi= SumPhi+ CellPhi
    endif

    ! normalize
    VoidFraction=VoidFraction/SumPhi
  end function clc_VoidFraction_IBM
