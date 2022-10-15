module ca_IBM
  use MPI
  use m_LogInfo
  use m_Parameters
  use Prtcl_TypeDef
  use Prtcl_Property
  use Prtcl_Decomp_2d
  use Prtcl_Variables
  use Prtcl_Parameters
  use m_Variables,only:mb1,hi1
  use Prtcl_Comm,only: DEM_Comm,sendlist
  use m_Decomp2d,only:nrank,y1start,y1end
  use ca_BC_and_Halo,only: Gather_Halo_IBMForce
  use m_MeshAndMetries,only: dx,dyUniform,dz,rdx,rdyUniform,rdz,xc,yc,zc
  implicit none
  private

#define THIS_PROC_XP 7
#define THIS_PROC_XM 8
#define THIS_PROC_YP 9
#define THIS_PROC_YM 10
#define THIS_PROC_ZP 11
#define THIS_PROC_ZM 12

  real(RK):: CellRadius,CellVolumeIBM,dxhalf,dyhalf,dzhalf,SMALL

  ! Ghost particle variables for IBM
  integer:: nlocal
  integer:: nGhostIBM
  integer:: mGhostIBM
  logical,dimension(3)::pbc
  real(RK),dimension(6)::dx_pbc
  real(RK),dimension(6)::dz_pbc
  integer:: GIBM_size_forward,GIBM_size_backward
  real(RK)::xstCoord,xedCoord,ystCoord,yedCoord,zstCoord,zedCoord

  integer,dimension(:),allocatable:: GhostP_Direction
  type(real3),dimension(:,:),allocatable:: GhostP_FluidIntegrate

  ! Immersed boundary points parts
  integer:: nIBP         ! number of Immersed Boundary Point in local processor
  integer:: mIBP         ! the possible maxium # of Immersed Boundary Point in local processor
  integer:: numIBP       ! total number of IBPs in all processors
  real(RK),dimension(:),allocatable::    IbpVolRatio

  integer,dimension(:),allocatable::     IBP_idlocal
  integer(kind=2),dimension(:,:),allocatable::IBP_indxyz
  type(real3),dimension(:),allocatable:: IBP_Pos
  type(real3),dimension(:),allocatable:: IBP_Vel
  type(real3),dimension(:),allocatable:: IBP_Force

  ! fixed particle variables
  integer:: nIBPFix                               ! number of Immersed Boundary Point in local processor
  integer(kind=2),dimension(:,:),allocatable::IBPFix_indxyz
  real(RK),dimension(:),allocatable::    IBPFix_VolRatio
  type(real3),dimension(:),allocatable:: IBPFix_Pos
  type(real3),dimension(:),allocatable:: IBPFix_Vel

  ABSTRACT INTERFACE
    subroutine pack_IBM_FpForce_x(buf_send,nsend)
      use m_TypeDef,only:RK
      implicit none
      real(RK),dimension(:),intent(out)::buf_send
      integer,intent(in)::nsend
    end subroutine pack_IBM_FpForce_x

    subroutine unpack_IBM_FpForce_x(buf_recv,nrecv)
      use m_TypeDef,only:RK
      implicit none
      real(RK),dimension(:),intent(in):: buf_recv
      integer,intent(in)::nrecv
    end subroutine unpack_IBM_FpForce_x
  END INTERFACE
  procedure(pack_IBM_FpForce_x  ),pointer:: pack_IBM_FpForce
  procedure(unpack_IBM_FpForce_x),pointer:: unpack_IBM_FpForce
  procedure(),pointer::IntegrateFluidPrtclForce

  ! public variables and functions/subroutines
  public:: Init_IBM
  public:: prepareIBM_interp,  InterpolationAndForcingIBM, SpreadIbpForce, updateRHSIBM
  public:: AdditionalForceIBM, IntegrateFluidPrtclForce,   Update_FluidIndicator

contains
#include "ca_IBM_common_inc.f90"

  !******************************************************************
  ! Init_IBM
  !******************************************************************
  subroutine Init_IBM(chFile)
    implicit none
    character(*),intent(in)::chFile
    
    ! locals
    real(RK)::vol_tot,vol_local,maxR,dxyz,min_xz,max_dxyz
    integer:: i,iErr01,iErr02,iErr03,iErr04,iErr05,ierror

    ! check integer(kind=2) is enough or not.
    if(nrank==0 .and. (nxc>huge(0_2)-20 .or. nyc>huge(0_2)-20 .or. nzc>huge(0_2)-20)) then
      call MainLog%CheckForError(ErrT_Abort,"Init_IBM","kind=2 is not enough for indxyz and IBPFix_indxyz")
    endif
        
    select case(IBM_Scheme)
    case(0)   ! 0: Explicit,Uhlmann(2005,JCP)

      ! (idlocal 1) +(direction 1) +(ptype 1) +(Pos 3) +(linvel 3) +(rotvel 3)=12
      ! (idlocal 1) +(FpForce 3) +(FpTorque 3) = 7
      GIBM_size_forward = 12
      GIBM_size_backward= 7

      pack_IBM_FpForce   => pack_IBM_FpForce_0
      unpack_IBM_FpForce => unpack_IBM_FpForce_0
      IntegrateFluidPrtclForce => IntegrateFluidPrtclForce_0
       
    case(1)   ! 1: Explicit,Kempe(2012,JCP)

      ! (idlocal 1) +(direction 1) +(ptype 1) +(Pos 3) +(linvel 3) +(rotvel 3)=12
      ! (idlocal 1) +(FpForce 3) +(FpTorque 3) +(FVelInt 3) +(FTorInt 3)= 13
      GIBM_size_forward = 12
      GIBM_size_backward= 13

      pack_IBM_FpForce   => pack_IBM_FpForce_1
      unpack_IBM_FpForce => unpack_IBM_FpForce_1
      IntegrateFluidPrtclForce => IntegrateFluidPrtclForce_1
    end select
    
    xstCoord= DEM_decomp%xSt
    xedCoord= DEM_decomp%xEd
    ystCoord= DEM_decomp%ySt
    yedCoord= DEM_decomp%yEd
    zstCoord= DEM_decomp%zSt
    zedCoord= DEM_decomp%zEd
    vol_tot   = xlx *yly *zlz
    vol_local = (xedCoord-xstCoord)*(yedCoord-ystCoord)*(zedCoord-zstCoord)

    maxR = maxval(DEMProperty%Prtcl_PureProp%Radius)
    SMALL=max(1.0E-14*(max(xlx,max(yly,zlz))), 1.0E-9*maxR)
    max_dxyz=max(max(dx,dyUniform),dz)
    min_xz=min(xedCoord-xstCoord,zedCoord-zstCoord)
    if(min_xz<=maxR+ two*max_dxyz) then
      call MainLog%CheckForError(ErrT_Abort,"Init_IBM","so big Diameter")
    endif

    pbc= DEM_Opt%IsPeriodic
    dx_pbc=zero; dz_pbc=zero
    if(DEM_decomp%Prtcl_Pencil/=y_pencil .and. nrank==0)then
      call MainLog%CheckForError(ErrT_Abort,"Init_IBM","only y_pencil is allowed in IBM")
    endif
    if(pbc(1)) then
      if(DEM_decomp%coord1==0)                 dx_pbc(xm_dir)= xlx
      if(DEM_decomp%coord1==DEM_decomp%prow-1) dx_pbc(xp_dir)=-xlx 
    endif
    if(pbc(3)) then
      if(DEM_decomp%coord2==0)                 dz_pbc(zm_dir)= zlz
      if(DEM_decomp%coord2==DEM_decomp%pcol-1) dz_pbc(zp_dir)=-zlz
    endif

    mGhostIBM= GPrtcl_list%mGhost_CS
    allocate(GhostP_Direction(mGhostIBM),Stat=iErr01)
    allocate(GhostP_FluidIntegrate(2,mGhostIBM),  Stat=iErr02)
    ierror=abs(iErr01)+abs(iErr02)
    if(ierror/=0)   call MainLog%CheckForError(ErrT_Abort,"Init_IBM","Allocation failed 1")
    
    !=========================== Immersed boundary points parts ===========================!
    numIBP = 0
    CellRadius= half*sqrt(dx*dx+dyUniform*dyUniform+dz*dz)
    CellVolumeIBM= dx* dyUniform* dz
    dxyz  = CellVolumeIBM**(one/three)
    dxhalf=dx*half; dyhalf=dyUniform*half; dzhalf=dz*half
    allocate(IbpVolRatio(DEM_opt%numPrtcl_Type))
    do i=1,DEM_opt%numPrtcl_Type
      IbpVolRatio(i)= PrtclIBMProp(i)%IBPVolume  *(rdx*rdyUniform*rdz)
      numIBP= numIBP+ PrtclIBMProp(i)%nPartition *DEMProperty%nPrtcl_in_Bin(i)
    enddo
    mIBP = int(real(numIBP,RK)*vol_local/vol_tot)
    mIBP = int(1.25_RK*mIBP)
    mIBP = min(mIBP, numIBP)
    mIBP = max(mIBP, 10)

    allocate(IBP_idlocal(mIBP), Stat=iErr01)
    allocate(IBP_indxyz(6,mIBP),Stat=iErr02)
    allocate(IBP_Pos(mIBP),     Stat=iErr03)
    allocate(IBP_Vel(mIBP),     Stat=iErr04)
    allocate(IBP_Force(mIBP),   Stat=iErr05)
    ierror=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)+abs(iErr05)
    if(ierror/=0)   call MainLog%CheckForError(ErrT_Abort,"Init_IBM","Allocation failed 2")
    IBP_idlocal=0;   IBP_indxyz=0
    IBP_Pos=zero_r3; IBP_Vel=zero_r3; IBP_Force=zero_r3

    ! fixed particle variables
    call Init_IBPFix()
  end subroutine Init_IBM

  !******************************************************************
  ! Update_FluidIndicator
  !******************************************************************
  subroutine Update_FluidIndicator(FluidIndicator)
    implicit none
    character,dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::FluidIndicator

    ! locals
    type(real3)::Cposition,CenterDiff
    real(RK)::Radius,CenterDistY,CenterDistZ,CenterDist
    integer::i,j,k,pid,DomainIx1,DomainIx2,DomainIy1,DomainIy2,DomainIz1,DomainIz2
    
    FluidIndicator= 'F'
    call PComm_FluidIndicator()

    DO pid=1,GPrtcl_list%nlocal
      Cposition= GPrtcl_PosR(pid)
      Radius   = GPrtcl_PosR(pid)%w
      DomainIx1= floor((Cposition%x-Radius)*rdx )+1; DomainIx1=max(y1start(1),DomainIx1)
      DomainIx2= ceiling((Cposition%x+Radius)*rdx ); DomainIx2=min(y1end(1),  DomainIx2)
      DomainIy1= floor((Cposition%y-Radius)*rdyUniform )+1; DomainIy1=max(y1start(2),DomainIy1)
      DomainIy2= ceiling((Cposition%y+Radius)*rdyUniform ); DomainIy2=min(y1end(2),  DomainIy2)
      DomainIz1= floor((Cposition%z-Radius)*rdz )+1; DomainIz1=max(y1start(3),DomainIz1)
      DomainIz2= ceiling((Cposition%z+Radius)*rdz ); DomainIz2=min(y1end(3),  DomainIz2)
      do k=DomainIz1,DomainIz2
        CenterDiff%z= real(k,RK)*dz-dzhalf- Cposition%z
        CenterDistZ=CenterDiff%z*CenterDiff%z
        do j=DomainIy1,DomainIy2
          CenterDiff%y= real(j,RK)*dyUniform-dyhalf- Cposition%y
          CenterDistY=CenterDistZ+CenterDiff%y*CenterDiff%y
          do i=DomainIx1,DomainIx2
            CenterDiff%x= real(i,RK)*dx-dxhalf- Cposition%x
            CenterDist=Radius-sqrt(CenterDiff%x*CenterDiff%x +CenterDistY)
            if(CenterDist>zero) FluidIndicator(i,j,k)='P'
          enddo
        enddo
      enddo
    ENDDO
    DO pid=1,nGhostIBM
      Cposition= GhostP_PosR(pid)
      Radius   = GhostP_PosR(pid)%w
      DomainIx1= floor((Cposition%x-Radius)*rdx )+1; DomainIx1=max(y1start(1),DomainIx1)
      DomainIx2= ceiling((Cposition%x+Radius)*rdx ); DomainIx2=min(y1end(1),  DomainIx2)
      DomainIy1= floor((Cposition%y-Radius)*rdyUniform )+1; DomainIy1=max(y1start(2),DomainIy1)
      DomainIy2= ceiling((Cposition%y+Radius)*rdyUniform ); DomainIy2=min(y1end(2),  DomainIy2)
      DomainIz1= floor((Cposition%z-Radius)*rdz )+1; DomainIz1=max(y1start(3),DomainIz1)
      DomainIz2= ceiling((Cposition%z+Radius)*rdz ); DomainIz2=min(y1end(3),  DomainIz2)
      do k=DomainIz1,DomainIz2
        CenterDiff%z= real(k,RK)*dz-dzhalf- Cposition%z
        CenterDistZ=CenterDiff%z*CenterDiff%z
        do j=DomainIy1,DomainIy2
          CenterDiff%y= real(j,RK)*dyUniform-dyhalf- Cposition%y
          CenterDistY=CenterDistZ+CenterDiff%y*CenterDiff%y
          do i=DomainIx1,DomainIx2
            CenterDiff%x= real(i,RK)*dx-dxhalf- Cposition%x
            CenterDist=Radius-sqrt(CenterDiff%x*CenterDiff%x +CenterDistY)
            if(CenterDist>zero) FluidIndicator(i,j,k)='P'
          enddo
        enddo
      enddo
    ENDDO
    DO pid=1,GPrtcl_list%mlocalFix+ GPrtcl_list%nGhostFix_CS
      Cposition= GPFix_PosR(pid)
      Radius   = GPFix_PosR(pid)%w
      DomainIx1= floor((Cposition%x-Radius)*rdx )+1; DomainIx1=max(y1start(1),DomainIx1)
      DomainIx2= ceiling((Cposition%x+Radius)*rdx ); DomainIx2=min(y1end(1),  DomainIx2)
      DomainIy1= floor((Cposition%y-Radius)*rdyUniform )+1; DomainIy1=max(y1start(2),DomainIy1)
      DomainIy2= ceiling((Cposition%y+Radius)*rdyUniform ); DomainIy2=min(y1end(2),  DomainIy2)
      DomainIz1= floor((Cposition%z-Radius)*rdz )+1; DomainIz1=max(y1start(3),DomainIz1)
      DomainIz2= ceiling((Cposition%z+Radius)*rdz ); DomainIz2=min(y1end(3),  DomainIz2)
      do k=DomainIz1,DomainIz2
        CenterDiff%z= real(k,RK)*dz-dzhalf- Cposition%z
        CenterDistZ=CenterDiff%z*CenterDiff%z
        do j=DomainIy1,DomainIy2
          CenterDiff%y= real(j,RK)*dyUniform-dyhalf- Cposition%y
          CenterDistY=CenterDistZ+CenterDiff%y*CenterDiff%y
          do i=DomainIx1,DomainIx2
            CenterDiff%x= real(i,RK)*dx-dxhalf- Cposition%x
            CenterDist=Radius-sqrt(CenterDiff%x*CenterDiff%x +CenterDistY)
            if(CenterDist>zero) FluidIndicator(i,j,k)='P'
          enddo
        enddo
      enddo
    ENDDO
  end subroutine Update_FluidIndicator

  !******************************************************************
  ! PrepareIBM_interp
  !******************************************************************
  subroutine PrepareIBM_interp()
    implicit none

    ! locals
    integer:: i,j,itype,nPartition,nIbpSum,ierror
    integer:: ic,jc,kc,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    type(real3):: Cposition,IbpPos,PosDiff,LinvelP,RotVelP

    nlocal= GPrtcl_list%nlocal
    call PComm_IBM_forward()
   
    !NOTE (Gong Zheng, 2020/07/01) : 
    ! [1] Here ONLY the IBP points within the Processor physical Domains will be considered.
    !       The IBP points outsides will be skipped. 
    ! [2] ONLY the uniform meshes near y-dir are used.
    nIBP= 0
    DO i=1,nlocal
      itype = GPrtcl_pType(i)
      Cposition = GPrtcl_PosR(i)
      LinvelP= GPrtcl_linVel(1,i)
      RotVelP= GPrtcl_rotVel(1,i)
      nPartition= PrtclIBMProp(itype)%nPartition
      do j=1,nPartition
        PosDiff= SpherePartitionCoord(j,itype)
        IbpPos = Cposition+ PosDiff
        ic= floor(IbpPos%x*rdx)+1
        jc= floor(IbpPos%y*rdyUniform)+1
        kc= floor(IbpPos%z*rdz)+1
        if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 &
        .or. (jc-y1start(2))*(jc-y1end(2))>0) cycle
        nIBP= nIBP+1
        if(nIBP> mIBP) call Reallocate_IbpVar()
        IBP_idlocal(nIBP)= i   !!! note here
        IBP_Pos(nIBP)= IbpPos
        IBP_Vel(nIBP)= LinvelP+ (RotVelP .cross. PosDiff)

        ! IBP_indxyz
        idxc_interp=ic-1
        idyc_interp=jc-1
        idzc_interp=kc-1
        if(IbpPos%x>xc(ic)) then
          idxp_interp= ic
        else
          idxp_interp= ic-1
        endif
        if(IbpPos%y>yc(jc)) then
          idyp_interp= jc
        else
          idyp_interp= jc-1
        endif
        if(IbpPos%z>zc(kc)) then
          idzp_interp= kc
        else
          idzp_interp= kc-1
        endif
        IBP_indxyz(1,nIBP)=idxc_interp
        IBP_indxyz(2,nIBP)=idxp_interp
        IBP_indxyz(3,nIBP)=idyc_interp
        IBP_indxyz(4,nIBP)=idyp_interp
        IBP_indxyz(5,nIBP)=idzc_interp
        IBP_indxyz(6,nIBP)=idzp_interp
      enddo
    ENDDO
    
    DO i=1,nGhostIBM
      itype = GhostP_pType(i)
      Cposition = GhostP_PosR(i)
      LinvelP= GhostP_linVel(i)
      RotVelP= GhostP_rotVel(i)
      nPartition= PrtclIBMProp(itype)%nPartition
      do j=1,nPartition
        PosDiff= SpherePartitionCoord(j,itype)
        IbpPos = Cposition+ PosDiff
        ic= floor(IbpPos%x*rdx)+1
        jc= floor(IbpPos%y*rdyUniform)+1
        kc= floor(IbpPos%z*rdz)+1        
        if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 &
        .or. (jc-y1start(2))*(jc-y1end(2))>0) cycle
        nIBP= nIBP+1
        if(nIBP> mIBP) call Reallocate_IbpVar()
        IBP_idlocal(nIBP)= i+nlocal   !!! note here
        IBP_Pos(nIBP)= IbpPos
        IBP_Vel(nIBP)= LinvelP+ (RotVelP .cross. PosDiff)

        ! IBP_indxyz
        idxc_interp=ic-1
        idyc_interp=jc-1
        idzc_interp=kc-1
        if(IbpPos%x>xc(ic)) then
          idxp_interp= ic
        else
          idxp_interp= ic-1
        endif
        if(IbpPos%y>yc(jc)) then
          idyp_interp= jc
        else
          idyp_interp= jc-1
        endif
        if(IbpPos%z>zc(kc)) then
          idzp_interp= kc
        else
          idzp_interp= kc-1
        endif
        IBP_indxyz(1,nIBP)=idxc_interp
        IBP_indxyz(2,nIBP)=idxp_interp
        IBP_indxyz(3,nIBP)=idyc_interp
        IBP_indxyz(4,nIBP)=idyp_interp
        IBP_indxyz(5,nIBP)=idzc_interp
        IBP_indxyz(6,nIBP)=idzp_interp
      enddo
    ENDDO
    call MPI_REDUCE(nIBP,nIbpSum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      if(nIbpSum>numIBP) call MainLog%CheckForError(ErrT_Abort," PrepareIBM_interp","sum(nIbp)>numIBP, WRONG!!!")
    endif
  end subroutine PrepareIBM_interp

  !******************************************************************
  ! InterpolationAndForcingIBM
  !******************************************************************
  subroutine InterpolationAndForcingIBM(uxStar,uyStar,uzStar)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in) ::uxStar,uyStar,uzStar

    ! locals
    type(real3):: IbpPos
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir
    real(RK),dimension(0:2)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    DO pid=1,nIBP
      IbpPos     = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)

      prxp= (IbpPos%x-xc(idxp_interp))*rdx+half
      RatioXp(0)= Delta_fun(prxp)
      RatioXp(1)= Delta_fun(prxp-one)
      RatioXp(2)= one- RatioXp(0)- RatioXp(1)
      prxc= (IbpPos%x-xc(idxc_interp))*rdx
      RatioXc(0)= Delta_fun(prxc)
      RatioXc(1)= Delta_fun(prxc-one)
      RatioXc(2)= one- RatioXc(0)- RatioXc(1)

      pryp= (IbpPos%y-yc(idyp_interp))*rdyUniform+half
      RatioYp(0)= Delta_fun(pryp)
      RatioYp(1)= Delta_fun(pryp-one)
      RatioYp(2)= one- RatioYp(0)- RatioYp(1)
      pryc= (IbpPos%y-yc(idyc_interp))*rdyUniform
      RatioYc(0)= Delta_fun(pryc)
      RatioYc(1)= Delta_fun(pryc-one)
      RatioYc(2)= one- RatioYc(0)- RatioYc(1)

      przp= (IbpPos%z-zc(idzp_interp))*rdz+half
      RatioZp(0)= Delta_fun(przp)
      RatioZp(1)= Delta_fun(przp-one)
      RatioZp(2)= one- RatioZp(0)- RatioZp(1)
      przc= (IbpPos%z-zc(idzc_interp))*rdz
      RatioZc(0)= Delta_fun(przc)
      RatioZc(1)= Delta_fun(przc-one)
      RatioZc(2)= one- RatioZc(0)- RatioZc(1)

      ! ux grid
      SumXDir=zero
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            SumXDir= SumXDir + uxStar(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uy gird
      SumYDir=zero
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumYDir= SumYDir + uyStar(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uz grid
      SumZDir=zero
      do k=0,2
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumZDir= SumZDir + uzStar(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo
      IBP_Force(pid)= IBP_Vel(pid) -real3(SumXDir,SumYDir,SumZDir)
    ENDDO

    ! fixed IBP part
    DO pid=1,nIBPFix
      IbpPos     = IBPFix_Pos(pid)
      idxc_interp= IBPFix_indxyz(1,pid)
      idxp_interp= IBPFix_indxyz(2,pid)
      idyc_interp= IBPFix_indxyz(3,pid)
      idyp_interp= IBPFix_indxyz(4,pid)
      idzc_interp= IBPFix_indxyz(5,pid)
      idzp_interp= IBPFix_indxyz(6,pid)

      prxp= (IbpPos%x-xc(idxp_interp))*rdx+half
      RatioXp(0)= Delta_fun(prxp)
      RatioXp(1)= Delta_fun(prxp-one)
      RatioXp(2)= one- RatioXp(0)- RatioXp(1)
      prxc= (IbpPos%x-xc(idxc_interp))*rdx
      RatioXc(0)= Delta_fun(prxc)
      RatioXc(1)= Delta_fun(prxc-one)
      RatioXc(2)= one- RatioXc(0)- RatioXc(1)

      pryp= (IbpPos%y-yc(idyp_interp))*rdyUniform+half
      RatioYp(0)= Delta_fun(pryp)
      RatioYp(1)= Delta_fun(pryp-one)
      RatioYp(2)= one- RatioYp(0)- RatioYp(1)
      pryc= (IbpPos%y-yc(idyc_interp))*rdyUniform
      RatioYc(0)= Delta_fun(pryc)
      RatioYc(1)= Delta_fun(pryc-one)
      RatioYc(2)= one- RatioYc(0)- RatioYc(1)

      przp= (IbpPos%z-zc(idzp_interp))*rdz+half
      RatioZp(0)= Delta_fun(przp)
      RatioZp(1)= Delta_fun(przp-one)
      RatioZp(2)= one- RatioZp(0)- RatioZp(1)
      przc= (IbpPos%z-zc(idzc_interp))*rdz
      RatioZc(0)= Delta_fun(przc)
      RatioZc(1)= Delta_fun(przc-one)
      RatioZc(2)= one- RatioZc(0)- RatioZc(1)

      ! ux grid
      SumXDir=zero
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            SumXDir= SumXDir + uxStar(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uy gird
      SumYDir=zero
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumYDir= SumYDir + uyStar(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uz grid
      SumZDir=zero
      do k=0,2
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumZDir= SumZDir + uzStar(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo
      IBPFix_Vel(pid)=  real3(-SumXDir,-SumYDir,-SumZDir)
    ENDDO

  end subroutine InterpolationAndForcingIBM

  !******************************************************************
  ! SpreadIbpForce
  !****************************************************************** 
  subroutine SpreadIbpForce(VolForce_x,VolForce_y,VolForce_z)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout) ::VolForce_x,VolForce_y,VolForce_z

    ! locals
    type(real3):: IbpPos,IbpForce
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,VolRatio
    real(RK),dimension(0:2)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    integer::i,j,k,id,jd,kd,pid,itype,idlocal
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    VolForce_x=zero; VolForce_y=zero; VolForce_z=zero
    DO pid=1,nIBP
      IbpPos     = IBP_Pos(pid)
      IbpForce   = IBP_Force(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
      
      idlocal    = IBP_idlocal(pid)
      if(idlocal> nlocal) then
        itype= GhostP_pType(idlocal- nlocal)
      else
        itype= GPrtcl_Ptype(idlocal)
      endif
      VolRatio= IbpVolRatio(itype)
      
      prxp= (IbpPos%x-xc(idxp_interp))*rdx+half
      RatioXp(0)= Delta_fun(prxp)
      RatioXp(1)= Delta_fun(prxp-one)
      RatioXp(2)= one- RatioXp(0)- RatioXp(1)
      prxc= (IbpPos%x-xc(idxc_interp))*rdx
      RatioXc(0)= Delta_fun(prxc)
      RatioXc(1)= Delta_fun(prxc-one)
      RatioXc(2)= one- RatioXc(0)- RatioXc(1)

      pryp= (IbpPos%y-yc(idyp_interp))*rdyUniform+half
      RatioYp(0)= Delta_fun(pryp)
      RatioYp(1)= Delta_fun(pryp-one)
      RatioYp(2)= one- RatioYp(0)- RatioYp(1)
      pryc= (IbpPos%y-yc(idyc_interp))*rdyUniform
      RatioYc(0)= Delta_fun(pryc)
      RatioYc(1)= Delta_fun(pryc-one)
      RatioYc(2)= one- RatioYc(0)- RatioYc(1)

      przp= (IbpPos%z-zc(idzp_interp))*rdz+half
      RatioZp(0)= Delta_fun(przp)
      RatioZp(1)= Delta_fun(przp-one)
      RatioZp(2)= one- RatioZp(0)- RatioZp(1)
      przc= (IbpPos%z-zc(idzc_interp))*rdz
      RatioZc(0)= Delta_fun(przc)
      RatioZc(1)= Delta_fun(przc-one)
      RatioZc(2)= one- RatioZc(0)- RatioZc(1)

      ! ux grid
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)*VolRatio
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            VolForce_x(id,jd,kd)= VolForce_x(id,jd,kd)+ prx*pry* IbpForce%x
          enddo
        enddo
      enddo

      ! uy gird
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)*VolRatio
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            VolForce_y(id,jd,kd)= VolForce_y(id,jd,kd)+ prx*pry* IbpForce%y
          enddo
        enddo
      enddo

      ! uz grid
      do k=0,2
        kd = k+idzp_interp
        prz= RatioZp(k)*VolRatio
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            VolForce_z(id,jd,kd)= VolForce_z(id,jd,kd)+ prx*pry* IbpForce%z
          enddo
        enddo
      enddo
    ENDDO

    ! fixed IBP part
    DO pid=1,nIBPFix
      IbpPos     = IBPFix_Pos(pid)
      IbpForce   = IBPFix_Vel(pid)
      idxc_interp= IBPFix_indxyz(1,pid)
      idxp_interp= IBPFix_indxyz(2,pid)
      idyc_interp= IBPFix_indxyz(3,pid)
      idyp_interp= IBPFix_indxyz(4,pid)
      idzc_interp= IBPFix_indxyz(5,pid)
      idzp_interp= IBPFix_indxyz(6,pid)
     
      VolRatio= IBPFix_VolRatio(pid)
      
      prxp= (IbpPos%x-xc(idxp_interp))*rdx+half
      RatioXp(0)= Delta_fun(prxp)
      RatioXp(1)= Delta_fun(prxp-one)
      RatioXp(2)= one- RatioXp(0)- RatioXp(1)
      prxc= (IbpPos%x-xc(idxc_interp))*rdx
      RatioXc(0)= Delta_fun(prxc)
      RatioXc(1)= Delta_fun(prxc-one)
      RatioXc(2)= one- RatioXc(0)- RatioXc(1)

      pryp= (IbpPos%y-yc(idyp_interp))*rdyUniform+half
      RatioYp(0)= Delta_fun(pryp)
      RatioYp(1)= Delta_fun(pryp-one)
      RatioYp(2)= one- RatioYp(0)- RatioYp(1)
      pryc= (IbpPos%y-yc(idyc_interp))*rdyUniform
      RatioYc(0)= Delta_fun(pryc)
      RatioYc(1)= Delta_fun(pryc-one)
      RatioYc(2)= one- RatioYc(0)- RatioYc(1)

      przp= (IbpPos%z-zc(idzp_interp))*rdz+half
      RatioZp(0)= Delta_fun(przp)
      RatioZp(1)= Delta_fun(przp-one)
      RatioZp(2)= one- RatioZp(0)- RatioZp(1)
      przc= (IbpPos%z-zc(idzc_interp))*rdz
      RatioZc(0)= Delta_fun(przc)
      RatioZc(1)= Delta_fun(przc-one)
      RatioZc(2)= one- RatioZc(0)- RatioZc(1)

      ! ux grid
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)*VolRatio
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            VolForce_x(id,jd,kd)= VolForce_x(id,jd,kd)+ prx*pry* IbpForce%x
          enddo
        enddo
      enddo

      ! uy gird
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)*VolRatio
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            VolForce_y(id,jd,kd)= VolForce_y(id,jd,kd)+ prx*pry* IbpForce%y
          enddo
        enddo
      enddo

      ! uz grid
      do k=0,2
        kd = k+idzp_interp
        prz= RatioZp(k)*VolRatio
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            VolForce_z(id,jd,kd)= VolForce_z(id,jd,kd)+ prx*pry* IbpForce%z
          enddo
        enddo
      enddo
    ENDDO

    call Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)
  end subroutine SpreadIbpForce

  !******************************************************************
  ! updateRHSIBM
  !******************************************************************
  subroutine updateRHSIBM(RhsX,RhsY,RhsZ,VolForce_x,VolForce_y,VolForce_z)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX,RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ 
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::VolForce_x,VolForce_y,VolForce_z

    ! locals
    integer:: ic,jc,kc

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
           RhsX(ic,jc,kc)= RhsX(ic,jc,kc)+ VolForce_x(ic,jc,kc)
           RhsY(ic,jc,kc)= RhsY(ic,jc,kc)+ VolForce_y(ic,jc,kc)
           RhsZ(ic,jc,kc)= RhsZ(ic,jc,kc)+ VolForce_z(ic,jc,kc)        
        enddo
      enddo
    ENDDO  
  end subroutine updateRHSIBM

  !******************************************************************
  ! AdditionalForceIBM
  !******************************************************************
  subroutine AdditionalForceIBM(ux,uy,uz,VolForce_x,VolForce_y,VolForce_z)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,VolForce_x,VolForce_y,VolForce_z

    ! locals
    type(real3):: IbpPos,IbpForce
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir,VolRatio
    real(RK),dimension(0:2)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp,itype,idlocal

    VolForce_x=zero; VolForce_y=zero; VolForce_z=zero
    DO pid=1,nIBP
      !===== interpolation =====!
      IbpPos     = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)

      prxp= (IbpPos%x-xc(idxp_interp))*rdx+half
      RatioXp(0)= Delta_fun(prxp)
      RatioXp(1)= Delta_fun(prxp-one)
      RatioXp(2)= one- RatioXp(0)- RatioXp(1)
      prxc= (IbpPos%x-xc(idxc_interp))*rdx
      RatioXc(0)= Delta_fun(prxc)
      RatioXc(1)= Delta_fun(prxc-one)
      RatioXc(2)= one- RatioXc(0)- RatioXc(1)

      pryp= (IbpPos%y-yc(idyp_interp))*rdyUniform+half
      RatioYp(0)= Delta_fun(pryp)
      RatioYp(1)= Delta_fun(pryp-one)
      RatioYp(2)= one- RatioYp(0)- RatioYp(1)
      pryc= (IbpPos%y-yc(idyc_interp))*rdyUniform
      RatioYc(0)= Delta_fun(pryc)
      RatioYc(1)= Delta_fun(pryc-one)
      RatioYc(2)= one- RatioYc(0)- RatioYc(1)

      przp= (IbpPos%z-zc(idzp_interp))*rdz+half
      RatioZp(0)= Delta_fun(przp)
      RatioZp(1)= Delta_fun(przp-one)
      RatioZp(2)= one- RatioZp(0)- RatioZp(1)
      przc= (IbpPos%z-zc(idzc_interp))*rdz
      RatioZc(0)= Delta_fun(przc)
      RatioZc(1)= Delta_fun(przc-one)
      RatioZc(2)= one- RatioZc(0)- RatioZc(1)

      SumXDir=zero
      do k=0,2             ! ux grid
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            SumXDir= SumXDir + ux(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      SumYDir=zero
      do k=0,2             ! uy gird
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumYDir= SumYDir + uy(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      SumZDir=zero
      do k=0,2             ! uz grid
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumZDir= SumZDir + uz(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      !=====    forcing    =====!
      IbpForce= IBP_Vel(pid) - real3(SumXDir,SumYDir,SumZDir)
      IBP_Force(pid)= IBP_Force(pid)+IbpForce

      !=====   spreading   =====!
      idlocal    = IBP_idlocal(pid)
      if(idlocal> nlocal) then
        itype= GhostP_pType(idlocal- nlocal)
      else
        itype= GPrtcl_Ptype(idlocal)
      endif
      VolRatio= IbpVolRatio(itype)

      do k=0,2             ! ux grid
        kd = k+idzc_interp
        prz= RatioZc(k)*VolRatio
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            VolForce_x(id,jd,kd)= VolForce_x(id,jd,kd)+ prx*pry* IbpForce%x
          enddo
        enddo
      enddo

      do k=0,2             ! uy gird
        kd = k+idzc_interp
        prz= RatioZc(k)*VolRatio
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            VolForce_y(id,jd,kd)= VolForce_y(id,jd,kd)+ prx*pry* IbpForce%y
          enddo
        enddo
      enddo
      
      do k=0,2              ! uz grid
        kd = k+idzp_interp
        prz= RatioZp(k)*VolRatio
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            VolForce_z(id,jd,kd)= VolForce_z(id,jd,kd)+ prx*pry* IbpForce%z
          enddo
        enddo
      enddo
    ENDDO

    ! fixed IBP part
    DO pid=1,nIBPFix
      !===== interpolation =====!
      IbpPos     = IBPFix_Pos(pid)
      idxc_interp= IBPFix_indxyz(1,pid)
      idxp_interp= IBPFix_indxyz(2,pid)
      idyc_interp= IBPFix_indxyz(3,pid)
      idyp_interp= IBPFix_indxyz(4,pid)
      idzc_interp= IBPFix_indxyz(5,pid)
      idzp_interp= IBPFix_indxyz(6,pid)

      prxp= (IbpPos%x-xc(idxp_interp))*rdx+half
      RatioXp(0)= Delta_fun(prxp)
      RatioXp(1)= Delta_fun(prxp-one)
      RatioXp(2)= one- RatioXp(0)- RatioXp(1)
      prxc= (IbpPos%x-xc(idxc_interp))*rdx
      RatioXc(0)= Delta_fun(prxc)
      RatioXc(1)= Delta_fun(prxc-one)
      RatioXc(2)= one- RatioXc(0)- RatioXc(1)

      pryp= (IbpPos%y-yc(idyp_interp))*rdyUniform+half
      RatioYp(0)= Delta_fun(pryp)
      RatioYp(1)= Delta_fun(pryp-one)
      RatioYp(2)= one- RatioYp(0)- RatioYp(1)
      pryc= (IbpPos%y-yc(idyc_interp))*rdyUniform
      RatioYc(0)= Delta_fun(pryc)
      RatioYc(1)= Delta_fun(pryc-one)
      RatioYc(2)= one- RatioYc(0)- RatioYc(1)

      przp= (IbpPos%z-zc(idzp_interp))*rdz+half
      RatioZp(0)= Delta_fun(przp)
      RatioZp(1)= Delta_fun(przp-one)
      RatioZp(2)= one- RatioZp(0)- RatioZp(1)
      przc= (IbpPos%z-zc(idzc_interp))*rdz
      RatioZc(0)= Delta_fun(przc)
      RatioZc(1)= Delta_fun(przc-one)
      RatioZc(2)= one- RatioZc(0)- RatioZc(1)

      SumXDir=zero
      do k=0,2             ! ux grid
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            SumXDir= SumXDir + ux(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      SumYDir=zero
      do k=0,2             ! uy gird
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumYDir= SumYDir + uy(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      SumZDir=zero
      do k=0,2             ! uz grid
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumZDir= SumZDir + uz(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      !=====    forcing    =====!
      IbpForce= real3(-SumXDir,-SumYDir,-SumZDir)

      !=====   spreading   =====!
      VolRatio= IBPFix_VolRatio(pid)

      do k=0,2             ! ux grid
        kd = k+idzc_interp
        prz= RatioZc(k)*VolRatio
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            VolForce_x(id,jd,kd)= VolForce_x(id,jd,kd)+ prx*pry* IbpForce%x
          enddo
        enddo
      enddo

      do k=0,2             ! uy gird
        kd = k+idzc_interp
        prz= RatioZc(k)*VolRatio
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            VolForce_y(id,jd,kd)= VolForce_y(id,jd,kd)+ prx*pry* IbpForce%y
          enddo
        enddo
      enddo
      
      do k=0,2              ! uz grid
        kd = k+idzp_interp
        prz= RatioZp(k)*VolRatio
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            VolForce_z(id,jd,kd)= VolForce_z(id,jd,kd)+ prx*pry* IbpForce%z
          enddo
        enddo
      enddo
    ENDDO
    call Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)

    ! Velocity corrcetion
    DO k=y1start(3),y1end(3)
      do j=y1start(2),y1end(2)
        do i=y1start(1),y1end(1)
           ux(i,j,k)= ux(i,j,k)+ VolForce_x(i,j,k)
           uy(i,j,k)= uy(i,j,k)+ VolForce_y(i,j,k)
           uz(i,j,k)= uz(i,j,k)+ VolForce_z(i,j,k)        
        enddo
      enddo
    ENDDO  
  end subroutine AdditionalForceIBM

  !******************************************************************
  ! PComm_IBM_forward
  !******************************************************************
  subroutine PComm_IBM_forward()
    implicit none

    ! locals
    real(RK)::px,pxt,pz,pzt
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer::nsend,nsend2,nrecv,nrecv2,nsendg,ng,ngp,ngpp
    integer::i,ierror,request(4),SRstatus(MPI_STATUS_SIZE)
    
    ! This part is similar to what is done in subroutine PC_Comm_For_Cntct of file "Prtcl_comm.f90"(only y_pencil)
    ng=0

    ! Step1: send to xp_axis, and receive from xm_dir
    nsend=0; nrecv=0; ngp=ng; ngpp=ng
    IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,nlocal
        px = GPrtcl_PosR(i)%x
        pxt= px+ GPrtcl_PosR(i)%w
        if(pxt +SMALL> xedCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_id(ng)       = i
          GhostP_pType(ng)    = GPrtcl_ptype(i)
          GhostP_Direction(ng)= THIS_PROC_XP
          GhostP_PosR(ng)     = GPrtcl_PosR(i)
          GhostP_PosR(ng)%x   = px-xlx
          GhostP_linVel(ng)   = GPrtcl_linVel(1,i)
          GhostP_rotVel(ng)   = GPrtcl_rotVel(1,i)
        endif
      enddo  
    ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) THEN
      nsendg=nsend
      do i=1,nlocal
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
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward(buf_send,nsendg,nsend,xp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(1),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step2: send to xm_axis, and receive from xp_dir
    nsend=0; nrecv=0; ngp=ng; !ngpp=ng   
    IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,nlocal
        px = GPrtcl_PosR(i)%x
        pxt= px- GPrtcl_PosR(i)%w
        if(pxt-SMALL<xstCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_id(ng)     = i
          GhostP_pType(ng)  = GPrtcl_pType(i)
          GhostP_Direction(ng)= THIS_PROC_XM
          GhostP_PosR(ng)   = GPrtcl_PosR(i)
          GhostP_PosR(ng)%x = px+xlx
          GhostP_linVel(ng) = GPrtcl_linVel(1,i)
          GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
        endif
      enddo 
    ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) THEN
      nsendg=nsend
      do i=1,nlocal
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
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward(buf_send,nsendg,nsend,xm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(2),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward(buf_recv,ngp+1,ng)
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
          GhostP_id(ng)     = nlocal+ i
          GhostP_pType(ng)  = GhostP_pType(i)
          GhostP_Direction(ng)= THIS_PROC_ZP
          GhostP_PosR(ng)   = GhostP_PosR(i)
          GhostP_PosR(ng)%z = pz-zlz
          GhostP_linVel(ng) = GhostP_linVel(i)
          GhostP_rotVel(ng) = GhostP_rotVel(i)
        endif
      enddo

      do i=1,nlocal
        pz = GPrtcl_PosR(i)%z 
        pzt= pz+ GPrtcl_PosR(i)%w
        if(pzt+SMALL>zedCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_id(ng)     = i
          GhostP_pType(ng)  = GPrtcl_ptype(i)
          GhostP_Direction(ng)= THIS_PROC_ZP
          GhostP_PosR(ng)   = GPrtcl_PosR(i)
          GhostP_PosR(ng)%z = pz-zlz
          GhostP_linVel(ng) = GPrtcl_linVel(1,i)
          GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
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

      do i=1,nlocal
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
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward(buf_send,nsendg,nsend,zp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(3),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward(buf_recv,ngp+1,ng)
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
          GhostP_id(ng)     = nlocal+ i
          GhostP_pType(ng)  = GhostP_pType(i)
          GhostP_Direction(ng)= THIS_PROC_ZM
          GhostP_PosR(ng)   = GhostP_PosR(i)
          GhostP_PosR(ng)%z = pz+ zlz
          GhostP_linVel(ng) = GhostP_linVel(i)
          GhostP_rotVel(ng) = GhostP_rotVel(i)
        endif
      enddo

      do i=1,nlocal
        pz = GPrtcl_PosR(i)%z 
        pzt= pz- GPrtcl_PosR(i)%w
        if(pzt-SMALL<zstCoord ) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_id(ng)     = i
          GhostP_pType(ng)  = GPrtcl_ptype(i)
          GhostP_Direction(ng)= THIS_PROC_ZM
          GhostP_PosR(ng)   = GPrtcl_PosR(i)
          GhostP_PosR(ng)%z = pz+ zlz
          GhostP_linVel(ng) = GPrtcl_linVel(1,i)
          GhostP_rotVel(ng) = GPrtcl_rotVel(1,i)
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

      do i=1,nlocal
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
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward(buf_send,nsendg,nsend,zm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(4),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    nGhostIBM= ng
  end subroutine PComm_IBM_forward

  !******************************************************************
  ! IntegrateFluidPrtclForce_0
  !******************************************************************
  subroutine IntegrateFluidPrtclForce_0(iCountACM,ux,uy,uz)
    implicit none
    integer,intent(in)::iCountACM
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz

    ! locals
    real(RK)::DenAlpha
    type(real3)::FpVolume,Cposition,PosDiff
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer::i,pid,idlocal,itype,idGhost,nsend,nrecv,nsend2,nrecv2,ierror,request(4),SRstatus(MPI_STATUS_SIZE)
    
    ! In order to save memory, GhostP_linVel is used to story FluidPrtclForce for ghost particle temporarily.
    ! Similarily, GhostP_rotVel is used to story FluidPrtclTorque for ghost particle temporarily.
    GPrtcl_FpForce = zero_r3
    GPrtcl_FpTorque= zero_r3
    GhostP_linVel  = zero_r3
    GhostP_rotVel  = zero_r3
    DO pid=1,nIBP
      idlocal    = IBP_idlocal(pid)
      if(idlocal> nlocal) then
        idGhost= idlocal- nlocal
        itype= GhostP_pType(idGhost)
        FpVolume= IBP_Force(pid)*PrtclIBMProp(itype)%IBPVolume
        GhostP_linVel(idGhost)=  GhostP_linVel(idGhost)+  FpVolume

        Cposition = GhostP_PosR(idGhost)
        PosDiff= IBP_Pos(pid)- Cposition
        GhostP_rotVel(idGhost)=  GhostP_rotVel(idGhost)+  (PosDiff  .cross. FpVolume)
      else
        itype= GPrtcl_Ptype(idlocal)
        FpVolume= IBP_Force(pid)*PrtclIBMProp(itype)%IBPVolume
        GPrtcl_FpForce(idlocal) = GPrtcl_FpForce(idlocal)+ FpVolume

        Cposition = GPrtcl_PosR(idlocal)
        PosDiff= IBP_Pos(pid)- Cposition
        GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+ (PosDiff  .cross. FpVolume)
      endif
    ENDDO

    ! Step1: send to zp_axis, and receive from zm_dir
    nsend=0; nrecv=0
    IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.    
      DO i=1, nGhostIBM
        if(GhostP_Direction(i)/= THIS_PROC_ZM) cycle
        idlocal= GhostP_id(i)
        if(idlocal> nlocal) then
          idGhost= idlocal- nlocal
          GhostP_linVel(idGhost)  = GhostP_linVel(idGhost)  +GhostP_linVel(i)
          GhostP_rotVel(idGhost)  = GhostP_rotVel(idGhost)  +GhostP_rotVel(i)
        else
          GPrtcl_FpForce(idlocal) = GPrtcl_FpForce(idlocal) +GhostP_linVel(i)
          GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+GhostP_rotVel(i)
        endif
      ENDDO
    ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) THEN
      DO i=1, nGhostIBM
        if(GhostP_Direction(i) /= zm_dir) cycle
        nsend=nsend+1     
        if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 1, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 1, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),2,MPI_COMM_WORLD,request(1),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),2,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(1),SRstatus,ierror)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step2: send to zm_axis, and receive from zp_dir
    nsend=0; nrecv=0
    IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself. 
      DO i=1, nGhostIBM
        if(GhostP_Direction(i)/= THIS_PROC_ZP) cycle
        idlocal= GhostP_id(i)
        if(idlocal> nlocal) then
          idGhost= idlocal- nlocal
          GhostP_linVel(idGhost)  = GhostP_linVel(idGhost)  +GhostP_linVel(i)
          GhostP_rotVel(idGhost)  = GhostP_rotVel(idGhost)  +GhostP_rotVel(i)
        else
          GPrtcl_FpForce(idlocal) = GPrtcl_FpForce(idlocal) +GhostP_linVel(i)
          GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+GhostP_rotVel(i)
        endif
      ENDDO
    ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) THEN
      DO i=1, nGhostIBM
        if(GhostP_Direction(i) /= zp_dir) cycle
        nsend=nsend+1     
        if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 3, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 3, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),4,MPI_COMM_WORLD,request(2),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),4,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(2),SRstatus,ierror)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step3: send to xp_axis, and receive from xm_dir
    nsend=0; nrecv=0
    IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.    
      DO i=1, nGhostIBM
        if(GhostP_Direction(i)/= THIS_PROC_XM) cycle
        idlocal= GhostP_id(i)
        GPrtcl_FpForce(idlocal)=  GPrtcl_FpForce(idlocal) +GhostP_linVel(i)
        GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+GhostP_rotVel(i)
      ENDDO
    ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) THEN
      DO i=1, nGhostIBM
        if(GhostP_Direction(i) /= xm_dir) cycle
        nsend=nsend+1     
        if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 5, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 5, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),6,MPI_COMM_WORLD,request(3),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),6,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(3),SRstatus,ierror)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step4: send to xm_axis, and receive from xp_dir
    nsend=0; nrecv=0
    IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.    
      DO i=1, nGhostIBM
        if(GhostP_Direction(i)/= THIS_PROC_XP) cycle
        idlocal= GhostP_id(i)
        GPrtcl_FpForce(idlocal)=  GPrtcl_FpForce(idlocal) +GhostP_linVel(i)
        GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+GhostP_rotVel(i)
      ENDDO
    ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) THEN
      DO i=1, nGhostIBM
        if(GhostP_Direction(i) /= xp_dir) cycle
        nsend=nsend+1     
        if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 7, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 7, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),8,MPI_COMM_WORLD,request(4),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),8,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(4),SRstatus,ierror)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Finally, get the final GPrtcl_FpForce
    DenAlpha=-FluidDensity/pmAlpha
    if(iCountACM==1) then
    endif
    do i=1, nlocal
      GPrtcl_FpForce(i) = DenAlpha*GPrtcl_FpForce(i)
      GPrtcl_FpTorque(i)= DenAlpha*GPrtcl_FpTorque(i)
    enddo
  end subroutine IntegrateFluidPrtclForce_0

  !******************************************************************
  ! IntegrateFluidPrtclForce_1
  !******************************************************************
  subroutine IntegrateFluidPrtclForce_1(iCountACM,ux,uy,uz)
    implicit none
    integer,intent(in)::iCountACM
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz

    ! locals
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    type(real3)::FpVolume,Cposition,PosDiff,CenterDiff,CellVel
    real(RK)::CenterDistY,CenterDistZ,CenterDist,Radius,VoidFraction,DenAlpha
    integer::i,j,k,pid,idlocal,idGhost,itype,nsend,nrecv,nsend2,nrecv2,ierror,request(4)
    integer::DomainIx1,DomainIx2,DomainIy1,DomainIy2,DomainIz1,DomainIz2,SRstatus(MPI_STATUS_SIZE)
    
    ! In order to save memory, GhostP_linVel is used to story FluidPrtclForce for ghost particle temporarily.
    ! Similarily, GhostP_rotVel is used to story FluidPrtclTorque for ghost particle temporarily.
    GPrtcl_FpForce = zero_r3
    GPrtcl_FpTorque= zero_r3
    GhostP_linVel  = zero_r3
    GhostP_rotVel  = zero_r3
    DO pid=1,nIBP
      idlocal    = IBP_idlocal(pid)
      if(idlocal> nlocal) then
        idGhost= idlocal- nlocal
        itype= GhostP_pType(idGhost)
        FpVolume= IBP_Force(pid)*PrtclIBMProp(itype)%IBPVolume
        GhostP_linVel(idGhost)=  GhostP_linVel(idGhost)+  FpVolume

        Cposition = GhostP_PosR(idGhost)
        PosDiff= IBP_Pos(pid)- Cposition
        GhostP_rotVel(idGhost)=  GhostP_rotVel(idGhost)+  (PosDiff  .cross. FpVolume)
      else
        itype= GPrtcl_Ptype(idlocal)
        FpVolume= IBP_Force(pid)*PrtclIBMProp(itype)%IBPVolume
        GPrtcl_FpForce(idlocal) = GPrtcl_FpForce(idlocal)+ FpVolume

        Cposition = GPrtcl_PosR(idlocal)
        PosDiff= IBP_Pos(pid)- Cposition
        GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+ (PosDiff  .cross. FpVolume)
      endif
    ENDDO

    GPrtcl_FluidIntegrate = zero_r3
    GhostP_FluidIntegrate = zero_r3
    ! Here I assume the velocities outside the domain are always zero.
    DO pid=1,nlocal
      Cposition= GPrtcl_PosR(pid)
      Radius   = GPrtcl_PosR(pid)%w
      DomainIx1= floor((Cposition%x-Radius)*rdx )+1; DomainIx1=max(y1start(1),DomainIx1)
      DomainIx2= ceiling((Cposition%x+Radius)*rdx ); DomainIx2=min(y1end(1),  DomainIx2)
      DomainIy1= floor((Cposition%y-Radius)*rdyUniform )+1; DomainIy1=max(y1start(2),DomainIy1)
      DomainIy2= ceiling((Cposition%y+Radius)*rdyUniform ); DomainIy2=min(y1end(2),  DomainIy2)
      DomainIz1= floor((Cposition%z-Radius)*rdz )+1; DomainIz1=max(y1start(3),DomainIz1)
      DomainIz2= ceiling((Cposition%z+Radius)*rdz ); DomainIz2=min(y1end(3),  DomainIz2)
      do k=DomainIz1,DomainIz2
        CenterDiff%z= real(k,RK)*dz-dzhalf- Cposition%z
        CenterDistZ=CenterDiff%z*CenterDiff%z
        do j=DomainIy1,DomainIy2
          CenterDiff%y= real(j,RK)*dyUniform-dyhalf- Cposition%y
          CenterDistY=CenterDistZ+CenterDiff%y*CenterDiff%y
          do i=DomainIx1,DomainIx2
            CenterDiff%x= real(i,RK)*dx-dxhalf- Cposition%x
            CenterDist=Radius-sqrt(CenterDiff%x*CenterDiff%x +CenterDistY)
            if(CenterDist >= CellRadius ) then
              CellVel%x= half*(ux(i+1,j,k)+ ux(i,j,k))
              CellVel%y= half*(uy(i,j+1,k)+ uy(i,j,k))
              CellVel%z= half*(uz(i,j,k+1)+ uz(i,j,k))
              GPrtcl_FluidIntegrate(1,pid)=GPrtcl_FluidIntegrate(1,pid)+CellVel
              GPrtcl_FluidIntegrate(2,pid)=GPrtcl_FluidIntegrate(2,pid)+(CenterDiff .cross. CellVel)
            elseif(CenterDist+CellRadius > zero) then
              VoidFraction= clc_VoidFraction_IBM(CenterDiff,Radius)
              CellVel%x= half*(ux(i+1,j,k)+ ux(i,j,k))
              CellVel%y= half*(uy(i,j+1,k)+ uy(i,j,k))
              CellVel%z= half*(uz(i,j,k+1)+ uz(i,j,k))
              GPrtcl_FluidIntegrate(1,pid)=GPrtcl_FluidIntegrate(1,pid)+VoidFraction*CellVel
              GPrtcl_FluidIntegrate(2,pid)=GPrtcl_FluidIntegrate(2,pid)+VoidFraction*(CenterDiff .cross. CellVel)
            endif
          enddo
        enddo
      enddo
      GPrtcl_FluidIntegrate(1,pid)= CellVolumeIBM* GPrtcl_FluidIntegrate(1,pid)
      GPrtcl_FluidIntegrate(2,pid)= CellVolumeIBM* GPrtcl_FluidIntegrate(2,pid)
    ENDDO
    DO pid=1,nGhostIBM
      Cposition = GhostP_PosR(pid)
      Radius    = GhostP_PosR(pid)%w
      DomainIx1= floor((Cposition%x-Radius)*rdx )+1; DomainIx1=max(y1start(1),DomainIx1)
      DomainIx2= ceiling((Cposition%x+Radius)*rdx ); DomainIx2=min(y1end(1),  DomainIx2)
      DomainIy1= floor((Cposition%y-Radius)*rdyUniform )+1; DomainIy1=max(y1start(2),DomainIy1)
      DomainIy2= ceiling((Cposition%y+Radius)*rdyUniform ); DomainIy2=min(y1end(2),  DomainIy2)
      DomainIz1= floor((Cposition%z-Radius)*rdz )+1; DomainIz1=max(y1start(3),DomainIz1)
      DomainIz2= ceiling((Cposition%z+Radius)*rdz ); DomainIz2=min(y1end(3),  DomainIz2)
      do k=DomainIz1,DomainIz2
        CenterDiff%z= real(k,RK)*dz-dzhalf- Cposition%z
        CenterDistZ=CenterDiff%z*CenterDiff%z
        do j=DomainIy1,DomainIy2
          CenterDiff%y= real(j,RK)*dyUniform-dyhalf- Cposition%y
          CenterDistY=CenterDistZ+CenterDiff%y*CenterDiff%y
          do i=DomainIx1,DomainIx2
            CenterDiff%x= real(i,RK)*dx-dxhalf- Cposition%x
            CenterDist=Radius-sqrt(CenterDiff%x*CenterDiff%x +CenterDistY)
            if(CenterDist >= CellRadius ) then
              CellVel%x= half*(ux(i+1,j,k)+ ux(i,j,k))
              CellVel%y= half*(uy(i,j+1,k)+ uy(i,j,k))
              CellVel%z= half*(uz(i,j,k+1)+ uz(i,j,k))
              GhostP_FluidIntegrate(1,pid)=GhostP_FluidIntegrate(1,pid)+CellVel
              GhostP_FluidIntegrate(2,pid)=GhostP_FluidIntegrate(2,pid)+(CenterDiff .cross. CellVel)
            elseif(CenterDist+CellRadius > zero) then
              VoidFraction= clc_VoidFraction_IBM(CenterDiff,Radius)
              CellVel%x= half*(ux(i+1,j,k)+ ux(i,j,k))
              CellVel%y= half*(uy(i,j+1,k)+ uy(i,j,k))
              CellVel%z= half*(uz(i,j,k+1)+ uz(i,j,k))
              GhostP_FluidIntegrate(1,pid)=GhostP_FluidIntegrate(1,pid)+VoidFraction*CellVel
              GhostP_FluidIntegrate(2,pid)=GhostP_FluidIntegrate(2,pid)+VoidFraction*(CenterDiff .cross. CellVel)
            endif
          enddo
        enddo
      enddo
      GhostP_FluidIntegrate(1,pid)= CellVolumeIBM* GhostP_FluidIntegrate(1,pid)
      GhostP_FluidIntegrate(2,pid)= CellVolumeIBM* GhostP_FluidIntegrate(2,pid)
    ENDDO

    ! Step1: send to zp_axis, and receive from zm_dir
    nsend=0; nrecv=0
    IF(DEM_decomp%ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.    
      DO i=1, nGhostIBM
        if(GhostP_Direction(i)/= THIS_PROC_ZM) cycle
        idlocal= GhostP_id(i)
        if(idlocal> nlocal) then
          idGhost= idlocal- nlocal
          GhostP_linVel(idGhost)  = GhostP_linVel(idGhost)  +GhostP_linVel(i)
          GhostP_rotVel(idGhost)  = GhostP_rotVel(idGhost)  +GhostP_rotVel(i)
          GhostP_FluidIntegrate(1,idGhost)= GhostP_FluidIntegrate(1,idGhost)+ GhostP_FluidIntegrate(1,i)
          GhostP_FluidIntegrate(2,idGhost)= GhostP_FluidIntegrate(2,idGhost)+ GhostP_FluidIntegrate(2,i)
        else
          GPrtcl_FpForce(idlocal) = GPrtcl_FpForce(idlocal) +GhostP_linVel(i)
          GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+GhostP_rotVel(i)
          GPrtcl_FluidIntegrate(1,idlocal)= GPrtcl_FluidIntegrate(1,idlocal)+ GhostP_FluidIntegrate(1,i)
          GPrtcl_FluidIntegrate(2,idlocal)= GPrtcl_FluidIntegrate(2,idlocal)+ GhostP_FluidIntegrate(2,i)
        endif
      ENDDO
    ELSEIF(DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) THEN
      DO i=1, nGhostIBM
        if(GhostP_Direction(i) /= zm_dir) cycle
        nsend=nsend+1     
        if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 1, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 1, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),2,MPI_COMM_WORLD,request(1),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),2,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(1),SRstatus,ierror)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step2: send to zm_axis, and receive from zp_dir
    nsend=0; nrecv=0
    IF(DEM_decomp%ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself. 
      DO i=1, nGhostIBM
        if(GhostP_Direction(i)/= THIS_PROC_ZP) cycle
        idlocal= GhostP_id(i)
        if(idlocal> nlocal) then
          idGhost= idlocal- nlocal
          GhostP_linVel(idGhost)  = GhostP_linVel(idGhost)  +GhostP_linVel(i)
          GhostP_rotVel(idGhost)  = GhostP_rotVel(idGhost)  +GhostP_rotVel(i)
          GhostP_FluidIntegrate(1,idGhost)= GhostP_FluidIntegrate(1,idGhost)+ GhostP_FluidIntegrate(1,i)
          GhostP_FluidIntegrate(2,idGhost)= GhostP_FluidIntegrate(2,idGhost)+ GhostP_FluidIntegrate(2,i)
        else
          GPrtcl_FpForce(idlocal) = GPrtcl_FpForce(idlocal) +GhostP_linVel(i)
          GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+GhostP_rotVel(i)
          GPrtcl_FluidIntegrate(1,idlocal)= GPrtcl_FluidIntegrate(1,idlocal)+ GhostP_FluidIntegrate(1,i)
          GPrtcl_FluidIntegrate(2,idlocal)= GPrtcl_FluidIntegrate(2,idlocal)+ GhostP_FluidIntegrate(2,i)
        endif
      ENDDO
    ELSEIF(DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) THEN
      DO i=1, nGhostIBM
        if(GhostP_Direction(i) /= zp_dir) cycle
        nsend=nsend+1     
        if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 3, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 3, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),4,MPI_COMM_WORLD,request(2),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),4,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(2),SRstatus,ierror)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step3: send to xp_axis, and receive from xm_dir
    nsend=0; nrecv=0
    IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.    
      DO i=1, nGhostIBM
        if(GhostP_Direction(i)/= THIS_PROC_XM) cycle
        idlocal= GhostP_id(i)
        GPrtcl_FpForce(idlocal)=  GPrtcl_FpForce(idlocal) +GhostP_linVel(i)
        GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+GhostP_rotVel(i)
        GPrtcl_FluidIntegrate(1,idlocal)= GPrtcl_FluidIntegrate(1,idlocal)+ GhostP_FluidIntegrate(1,i)
        GPrtcl_FluidIntegrate(2,idlocal)= GPrtcl_FluidIntegrate(2,idlocal)+ GhostP_FluidIntegrate(2,i)
      ENDDO
    ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) THEN
      DO i=1, nGhostIBM
        if(GhostP_Direction(i) /= xm_dir) cycle
        nsend=nsend+1     
        if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 5, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 5, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),6,MPI_COMM_WORLD,request(3),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),6,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(3),SRstatus,ierror)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step4: send to xm_axis, and receive from xp_dir
    nsend=0; nrecv=0
    IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.    
      DO i=1, nGhostIBM
        if(GhostP_Direction(i)/= THIS_PROC_XP) cycle
        idlocal= GhostP_id(i)
        GPrtcl_FpForce(idlocal)=  GPrtcl_FpForce(idlocal) +GhostP_linVel(i)
        GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+GhostP_rotVel(i)
        GPrtcl_FluidIntegrate(1,idlocal)= GPrtcl_FluidIntegrate(1,idlocal)+ GhostP_FluidIntegrate(1,i)
        GPrtcl_FluidIntegrate(2,idlocal)= GPrtcl_FluidIntegrate(2,idlocal)+ GhostP_FluidIntegrate(2,i)
      ENDDO
    ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) THEN
      DO i=1, nGhostIBM
        if(GhostP_Direction(i) /= xp_dir) cycle
        nsend=nsend+1     
        if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 7, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 7, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),8,MPI_COMM_WORLD,request(4),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),8,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(4),SRstatus,ierror)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Finally, get the final GPrtcl_FpForce
    DenAlpha=FluidDensity/pmAlpha
    if(iCountACM==1) then
      do i=1,nlocal
        GPrtcl_FpForce(i) = DenAlpha*(zero_r3-GPrtcl_FpForce(i))
        GPrtcl_FpTorque(i)= DenAlpha*(zero_r3-GPrtcl_FpTorque(i))
        GPrtcl_FluidIntOld(1,i)=GPrtcl_FluidIntegrate(1,i)
        GPrtcl_FluidIntOld(2,i)=GPrtcl_FluidIntegrate(2,i)
      enddo
    else
      do i=1, nlocal
        GPrtcl_FpForce(i) = DenAlpha*(GPrtcl_FluidIntegrate(1,i)-GPrtcl_FluidIntOld(1,i)-GPrtcl_FpForce(i))
        GPrtcl_FpTorque(i)= DenAlpha*(GPrtcl_FluidIntegrate(2,i)-GPrtcl_FluidIntOld(2,i)-GPrtcl_FpTorque(i))
        GPrtcl_FluidIntOld(1,i)=GPrtcl_FluidIntegrate(1,i)
        GPrtcl_FluidIntOld(2,i)=GPrtcl_FluidIntegrate(2,i)
      enddo
    endif
  end subroutine IntegrateFluidPrtclForce_1

  !**********************************************************************
  ! pack_IBM_FpForce_0
  !**********************************************************************
  subroutine pack_IBM_FpForce_0(buf_send,nsend)
    implicit none
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsend

    ! locals
    integer::i,pid,m

    m=1
    do i=1,nsend
      pid= sendlist(i)
      buf_send(m)= real(GhostP_id(pid)); m=m+1
      buf_send(m)= GhostP_linVel(pid)%x; m=m+1
      buf_send(m)= GhostP_linVel(pid)%y; m=m+1
      buf_send(m)= GhostP_linVel(pid)%z; m=m+1
      buf_send(m)= GhostP_rotVel(pid)%x; m=m+1
      buf_send(m)= GhostP_rotVel(pid)%y; m=m+1
      buf_send(m)= GhostP_rotVel(pid)%z; m=m+1
    enddo
  end subroutine pack_IBM_FpForce_0

  !**********************************************************************
  ! pack_IBM_FpForce_1
  !**********************************************************************
  subroutine pack_IBM_FpForce_1(buf_send,nsend)
    implicit none
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsend

    ! locals
    integer::i,pid,m

    m=1
    do i=1,nsend
      pid= sendlist(i)
      buf_send(m)= real(GhostP_id(pid)); m=m+1
      buf_send(m)= GhostP_linVel(pid)%x; m=m+1
      buf_send(m)= GhostP_linVel(pid)%y; m=m+1
      buf_send(m)= GhostP_linVel(pid)%z; m=m+1
      buf_send(m)= GhostP_rotVel(pid)%x; m=m+1
      buf_send(m)= GhostP_rotVel(pid)%y; m=m+1
      buf_send(m)= GhostP_rotVel(pid)%z; m=m+1
      buf_send(m)= GhostP_FluidIntegrate(1,pid)%x; m=m+1
      buf_send(m)= GhostP_FluidIntegrate(1,pid)%y; m=m+1
      buf_send(m)= GhostP_FluidIntegrate(1,pid)%z; m=m+1
      buf_send(m)= GhostP_FluidIntegrate(2,pid)%x; m=m+1
      buf_send(m)= GhostP_FluidIntegrate(2,pid)%y; m=m+1
      buf_send(m)= GhostP_FluidIntegrate(2,pid)%z; m=m+1
    enddo
  end subroutine pack_IBM_FpForce_1

  !**********************************************************************
  ! unpack_IBM_FpForce_0
  !**********************************************************************
  subroutine unpack_IBM_FpForce_0(buf_recv,nrecv)
    implicit none
    real(RK),dimension(:),intent(in):: buf_recv
    integer,intent(in)::nrecv

    !locals
    type(real3)::FpForce,FpTorque
    integer::i,pid,m,gid

    m=1
    do i=1,nrecv
      pid= int(buf_recv(m)+0.2); m=m+1
      FpForce%x =buf_recv(m);    m=m+1
      FpForce%y =buf_recv(m);    m=m+1
      FpForce%z =buf_recv(m);    m=m+1
      FpTorque%x=buf_recv(m);    m=m+1
      FpTorque%y=buf_recv(m);    m=m+1
      FpTorque%z=buf_recv(m);    m=m+1    
      if(pid>nlocal) then
        gid= pid-nlocal
        GhostP_linVel(gid)  = GhostP_linVel(gid)  +FpForce
        GhostP_rotVel(gid)  = GhostP_rotVel(gid)  +FpTorque
      else
        GPrtcl_FpForce(pid) = GPrtcl_FpForce(pid) +FpForce
        GPrtcl_FpTorque(pid)= GPrtcl_FpTorque(pid)+FpTorque
      endif   
    enddo
  end subroutine unpack_IBM_FpForce_0

  !**********************************************************************
  ! unpack_IBM_FpForce_1
  !**********************************************************************
  subroutine unpack_IBM_FpForce_1(buf_recv,nrecv)
    implicit none
    real(RK),dimension(:),intent(in):: buf_recv
    integer,intent(in)::nrecv

    ! locals
    integer::i,pid,m,gid
    type(real3)::FpForce,FpTorque,FVelInt,FTorInt

    m=1
    do i=1,nrecv
      pid= int(buf_recv(m)+0.2); m=m+1
      FpForce%x =buf_recv(m);    m=m+1
      FpForce%y =buf_recv(m);    m=m+1
      FpForce%z =buf_recv(m);    m=m+1
      FpTorque%x=buf_recv(m);    m=m+1
      FpTorque%y=buf_recv(m);    m=m+1
      FpTorque%z=buf_recv(m);    m=m+1 
      FVelInt%x =buf_recv(m);    m=m+1
      FVelInt%y =buf_recv(m);    m=m+1
      FVelInt%z =buf_recv(m);    m=m+1
      FTorInt%x =buf_recv(m);    m=m+1
      FTorInt%y =buf_recv(m);    m=m+1
      FTorInt%z =buf_recv(m);    m=m+1   
      if(pid>nlocal) then
        gid= pid-nlocal
        GhostP_linVel(gid)  = GhostP_linVel(gid)  +FpForce
        GhostP_rotVel(gid)  = GhostP_rotVel(gid)  +FpTorque
        GhostP_Fluidintegrate(1,gid) = GhostP_Fluidintegrate(1,gid) +FVelInt
        GhostP_Fluidintegrate(2,gid) = GhostP_Fluidintegrate(2,gid) +FTorInt
      else
        GPrtcl_FpForce(pid) = GPrtcl_FpForce(pid) +FpForce
        GPrtcl_FpTorque(pid)= GPrtcl_FpTorque(pid)+FpTorque
        GPrtcl_FluidIntegrate(1,pid) = GPrtcl_FluidIntegrate(1,pid) +FVelInt
        GPrtcl_FluidIntegrate(2,pid) = GPrtcl_FluidIntegrate(2,pid) +FTorInt
      endif   
    enddo
  end subroutine unpack_IBM_FpForce_1

  !**********************************************************************
  ! pack_IBM_forward
  !**********************************************************************
  subroutine pack_IBM_forward(buf_send,nsendg,nsend,dir)
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
      buf_send(m)=real(nlocal+ id);       m=m+1 ! 01
      buf_send(m)=real(GhostP_pType(id)); m=m+1 ! 02
      buf_send(m)=real(dir);              m=m+1 ! 03
      buf_send(m)=GhostP_PosR(id)%x+pdx;  m=m+1 ! 04
      buf_send(m)=GhostP_PosR(id)%y;      m=m+1 ! 05
      buf_send(m)=GhostP_PosR(id)%z+pdz;  m=m+1 ! 06
      buf_send(m)=GhostP_linVel(id)%x;    m=m+1 ! 07
      buf_send(m)=GhostP_linVel(id)%y;    m=m+1 ! 08
      buf_send(m)=GhostP_linVel(id)%z;    m=m+1 ! 09
      buf_send(m)=GhostP_rotVel(id)%x;    m=m+1 ! 10
      buf_send(m)=GhostP_rotVel(id)%y;    m=m+1 ! 11
      buf_send(m)=GhostP_rotVel(id)%z;    m=m+1 ! 12 
    enddo
    do i=nsendg+1,nsend
      id=sendlist(i)
      buf_send(m)=real(id);               m=m+1 ! 01
      buf_send(m)=real(GPrtcl_pType(id)); m=m+1 ! 02
      buf_send(m)=real(dir);              m=m+1 ! 03
      buf_send(m)=GPrtcl_PosR(id)%x+pdx;  m=m+1 ! 04
      buf_send(m)=GPrtcl_PosR(id)%y;      m=m+1 ! 05
      buf_send(m)=GPrtcl_PosR(id)%z+pdz;  m=m+1 ! 06
      buf_send(m)=GPrtcl_linVel(1,id)%x;  m=m+1 ! 07
      buf_send(m)=GPrtcl_linVel(1,id)%y;  m=m+1 ! 08
      buf_send(m)=GPrtcl_linVel(1,id)%z;  m=m+1 ! 09
      buf_send(m)=GPrtcl_rotVel(1,id)%x;  m=m+1 ! 10
      buf_send(m)=GPrtcl_rotVel(1,id)%y;  m=m+1 ! 11
      buf_send(m)=GPrtcl_rotVel(1,id)%z;  m=m+1 ! 12 
    enddo
  end subroutine pack_IBM_forward

  !**********************************************************************
  ! unpack_IBM_forward
  !**********************************************************************
  subroutine unpack_IBM_forward(buf_recv,n1,n2)
    implicit none
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::n1,n2

    ! locals
    integer::i,m,itype
   
    m=1
    do i=n1,n2
      GhostP_id(i)       =int(buf_recv(m)+0.2); m=m+1 ! 01 
      itype              =int(buf_recv(m)+0.2); m=m+1 ! 02
      GhostP_Direction(i)=int(buf_recv(m)+0.2); m=m+1 ! 03
      GhostP_PosR(i)%x   =buf_recv(m);          m=m+1 ! 04
      GhostP_PosR(i)%y   =buf_recv(m);          m=m+1 ! 05
      GhostP_PosR(i)%z   =buf_recv(m);          m=m+1 ! 06
      GhostP_linVel(i)%x =buf_recv(m);          m=m+1 ! 07 
      GhostP_linVel(i)%y =buf_recv(m);          m=m+1 ! 08 
      GhostP_linVel(i)%z =buf_recv(m);          m=m+1 ! 09 
      GhostP_rotVel(i)%x =buf_recv(m);          m=m+1 ! 10 
      GhostP_rotVel(i)%y =buf_recv(m);          m=m+1 ! 11 
      GhostP_rotVel(i)%z =buf_recv(m);          m=m+1 ! 12
      GhostP_pType(i)    =itype;                
      GhostP_PosR(i)%w   =DEMProperty%Prtcl_PureProp(itype)%Radius 
    enddo
  end subroutine unpack_IBM_forward

  !******************************************************************
  ! Reallocate_ghost_for_IBM
  !******************************************************************
  subroutine Reallocate_ghost_for_IBM(ng)
    implicit none
    integer,intent(in)::ng

    ! locals
    integer:: sizep,ierrTemp,ierror=0
    integer,dimension(:),allocatable::IntVec
    
    call DEM_Comm%reallocate_ghost_for_Cntct(ng)
    sizep= mGhostIBM; mGhostIBM= GPrtcl_list%mGhost_CS

    call move_alloc(GhostP_direction,IntVec)
    allocate(GhostP_direction(mGhostIBM),stat=ierrTemp);ierror=ierror+abs(ierrTemp)
    GhostP_direction(1:sizep)= IntVec
    deallocate(IntVec)

    deallocate(GhostP_FluidIntegrate)
    allocate(GhostP_FluidIntegrate(2,mGhostIBM),stat=ierrTemp);ierror=ierror+abs(ierrTemp)

    if(ierror/=0) then
      call MainLog%CheckForError(ErrT_Abort," Reallocate_ghost_for_IBM"," Reallocate wrong!")
      call MainLog%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif
    !call MainLog%CheckForError(ErrT_Pass," Reallocate_ghost_for_IBM"," Need to reallocate Ghost variables for IBM") 
  end subroutine  Reallocate_ghost_for_IBM

  !******************************************************************
  ! Reallocate_IbpVar
  !******************************************************************
  subroutine Reallocate_IbpVar()
    implicit none

    ! locals
    integer:: sizep,sizen,ierrTemp,ierror=0
    integer,dimension(:),allocatable:: IntVec
    type(real3),dimension(:),allocatable:: Real3Vec
    integer(kind=2),dimension(:,:),allocatable::IntArr   

    sizep= mIBP
    sizen= int(1.1_RK*real(sizep,kind=RK))
    sizen= max(sizen, nIBP+1)
    sizen= min(sizen, numIBP)
    mIBP = sizen    

    ! ======= integer vector part =======
    call move_alloc(IBP_idlocal,IntVec)
    allocate(IBP_idlocal(sizen),stat=ierrTemp);ierror=ierror+abs(ierrTemp)
    IBP_idlocal(1:sizep)=IntVec
    deallocate(IntVec)

    ! ======= integer matrix part =======
    call move_alloc(IBP_indxyz,IntArr)
    allocate(IBP_indxyz(6,sizen),stat=ierrTemp);ierror=ierror+abs(ierrTemp)
    IBP_indxyz(1:6,1:sizep)=IntArr
    deallocate(IntArr)

    ! ======= real3 vercor part =======
    call move_alloc(IBP_Pos,Real3Vec)
    allocate(IBP_Pos(sizen),stat=ierrTemp);ierror=ierror+abs(ierrTemp)
    IBP_Pos(1:sizep)=Real3Vec

    call move_alloc(IBP_Vel,Real3Vec)
    allocate(IBP_Vel(sizen),stat=ierrTemp);ierror=ierror+abs(ierrTemp)
    IBP_Vel(1:sizep)=Real3Vec
    deallocate(Real3Vec)

    deallocate(IBP_Force)
    allocate(IBP_Force(sizen),stat=ierrTemp);ierror=ierror+abs(ierrTemp)

    if(ierror/=0) then
      call MainLog%CheckForError(ErrT_Abort," Reallocate_IbpVar"," Reallocate wrong!")
      call MainLog%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif    
    !call MainLog%CheckForError(ErrT_Pass,"Reallocate_IbpVar","Need to reallocate IBP variables")
    !call MainLog%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    !call MainLog%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    !call MainLog%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine Reallocate_IbpVar

#undef THIS_PROC_XP
#undef THIS_PROC_XM
#undef THIS_PROC_YP
#undef THIS_PROC_YM
#undef THIS_PROC_ZP
#undef THIS_PROC_ZM
end module ca_IBM
