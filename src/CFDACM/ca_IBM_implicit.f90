module ca_IBM_implicit
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Parameters
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
#ifdef IBMDistribute2
#define nDistribute 1
#else
#define nDistribute 2
#endif
#define THIS_PROC_XP 7
#define THIS_PROC_XM 8
#define THIS_PROC_YP 9
#define THIS_PROC_YM 10
#define THIS_PROC_ZP 11
#define THIS_PROC_ZM 12

  real(RK):: CellRadius,CellVolumeIBM,dxhalf,dyhalf,dzhalf,SMALL

  ! Ghost particle variables for IBM
  integer:: nlocalOld
  integer:: mlocalNEW
  integer:: nGhostIBM
  integer:: mGhostIBM
  logical,dimension(3)::pbc
  real(RK),dimension(6)::dx_pbc
  real(RK),dimension(6)::dz_pbc
  integer:: GIBM_size_forward,GIBM_size_backward
  real(RK)::xstCoord,xedCoord,ystCoord,yedCoord,zstCoord,zedCoord

  integer, dimension(:),allocatable:: GhostP_Direction
  type(real3),dimension(:,:),allocatable:: GhostP_FluidIntegrate

  ! Immersed boundary points parts
  integer:: nIBP         ! number of Immersed Boundary Point in local processor
  integer:: mIBP         ! the possible maxium # of Immersed Boundary Point in local processor
  integer:: numIBP       ! total number of IBPs in all processors
  real(RK),dimension(:),allocatable::    IbpVolRatio

  integer,dimension(:),allocatable::     IBP_idlocal
  integer(kind=2),dimension(:,:),allocatable:: IBP_indxyz
  type(real3),dimension(:),allocatable:: IBP_Pos
  type(real3),dimension(:),allocatable:: IBP_Vel

  ! fixed particle variables
  integer:: nIBPFix                               ! number of Immersed Boundary Point in local processor
  integer(kind=2),dimension(:,:),allocatable:: IBPFix_indxyz
  real(RK),dimension(:),allocatable::    IBPFix_VolRatio
  type(real3),dimension(:),allocatable:: IBPFix_Pos
  type(real3),dimension(:),allocatable:: IBPFix_Vel

  ! Variabls for semi-implicit coupling only ==========================================
  integer:: nIBPOld         ! number of Immersed Boundary Point in local processor
  integer,dimension(:),allocatable::     PIBM_pType
  type(real3),dimension(:),allocatable:: PIBM_linVel
  type(real3),dimension(:),allocatable:: PIBM_rotVel
  integer:: GIBM_size_forward_imp,GIBM_size_recover_imp

  ! public variables and functions/subroutines
  public:: PrepareIBM_interp_imp1, PrepareIBM_interp_imp2
  public:: Init_IBM_imp,UpdateRhsIBM_imp1,UpdateRhsIBM_imp2,   Update_FluidIndicator_imp
  public:: InterpolationIBM_imp, IntegrateFluidPrtclForce_imp, SpreadIbpForce_imp
contains
#include "ca_IBM_common_inc.f90"

  !******************************************************************
  ! Init_IBM_imp
  !******************************************************************
  subroutine Init_IBM_imp(chFile)
    implicit none
    character(*),intent(in)::chFile
    
    ! locals
    real(RK)::vol_tot,vol_local,maxR,dxyz,min_xz,max_dxyz
    integer:: i,iErr01,iErr02,iErr03,iErr04,ierror

    ! check integer(kind=2) is enough or not.
    if(nrank==0 .and. (nxc>huge(0_2)-20 .or. nyc>huge(0_2)-20 .or. nzc>huge(0_2)-20)) then
      call MainLog%CheckForError(ErrT_Abort,"Init_IBM_imp","kind=2 is not enough for indxyz and IBPFix_indxyz")
    endif
    
    ! (idlocal 1) +(direction 1) +(ptype 1) +(Pos 3) =6
    ! (idlocal 1) +(direction 1) +(ptype 1) +(Pos 3) +(linvel 3) +(rotvel 3)=12
    ! (idlocal 1) +(FpForce 3) +(FpTorque 3) +(FVelInt 3) +(FTorInt 3)= 13
    ! (pType 1)   +(PosOld  3) +(linVel   3) +(rotVel  3) =10
    GIBM_size_forward_imp=  6
    GIBM_size_forward= 12
    GIBM_size_backward= 13
    GIBM_size_recover_imp= 10
    
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
    if(min_xz<=maxR+ 2.0_RK*max_dxyz) then
      call MainLog%CheckForError(ErrT_Abort,"Init_IBM","so big Diameter")
    endif

    pbc= DEM_Opt%IsPeriodic
    dx_pbc=0.0_RK; dz_pbc=0.0_RK
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
    CellRadius= 0.5_RK*sqrt(dx*dx+dyUniform*dyUniform+dz*dz)
    CellVolumeIBM= dx* dyUniform* dz
    dxyz  = CellVolumeIBM**(0.33333333333333333_RK)
    dxhalf=dx*0.5_RK; dyhalf=dyUniform*0.5_RK; dzhalf=dz*0.5_RK
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
    ierror=abs(iErr01)+abs(iErr02)+abs(iErr03)+abs(iErr04)
    if(ierror/=0)   call MainLog%CheckForError(ErrT_Abort,"Init_IBM","Allocation failed 2")
    IBP_idlocal=0;   IBP_indxyz=0
    IBP_Pos=zero_r3; IBP_Vel=zero_r3

    ! fixed particle variables
    call Init_IBPFix()

  end subroutine Init_IBM_imp

  !******************************************************************
  ! Update_FluidIndicator_imp
  !******************************************************************
  subroutine Update_FluidIndicator_imp(FluidIndicator)
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
            if(CenterDist>0.0_RK) FluidIndicator(i,j,k)='P'
          enddo
        enddo
      enddo
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
            if(CenterDist>0.0_RK) FluidIndicator(i,j,k)='P'
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
            if(CenterDist>0.0_RK) FluidIndicator(i,j,k)='P'
          enddo
        enddo
      enddo
    ENDDO
  end subroutine Update_FluidIndicator_imp

  !******************************************************************
  ! PrepareIBM_interp_imp1
  !******************************************************************
  subroutine PrepareIBM_interp_imp1()
    implicit none

    ! locals
    integer:: i,j,itype,nPartition,nIbpSum,ierror
    integer:: ic,jc,kc,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    type(real3):: Cposition,LagrangeP,PosDiff

    nlocalOld= GPrtcl_list%nlocal
    call PComm_IBM_forward_imp1()
   
    !NOTE (Gong Zheng, 2020/07/01) : 
    ! [1] Here ONLY the IBP points within the Processor physical Domains will be considered.
    !       The IBP points outsides will be skipped. 
    ! [2] ONLY the uniform meshes near y- dir are used.
    nIBP= 0
    DO i=1,nlocalOld
      itype = GPrtcl_pType(i)
      Cposition = GPrtcl_PosR(i)
      GPrtcl_PosOld(i)= GPrtcl_PosR(i)
      nPartition= PrtclIBMProp(itype)%nPartition
      do j=1,nPartition
        PosDiff= SpherePartitionCoord(j,itype)
        LagrangeP = Cposition+ PosDiff
        ic= floor(LagrangeP%x*rdx)+1
        jc= floor(LagrangeP%y*rdyUniform)+1
        kc= floor(LagrangeP%z*rdz)+1        
        if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 &
        .or. (jc-y1start(2))*(jc-y1end(2))>0) cycle
        nIBP= nIBP+1
        if(nIBP> mIBP) call Reallocate_IbpVar()
        IBP_idlocal(nIBP)= i   !!! note here
        IBP_Pos(nIBP)= LagrangeP
#define   clc_Point_indxyz_IBPMove
#include "clc_Point_indxyz_inc.f90"
#undef    clc_Point_indxyz_IBPMove
      enddo
    ENDDO
    
    DO i=1,nGhostIBM
      itype = GhostP_pType(i)
      Cposition = GhostP_PosR(i)
      nPartition= PrtclIBMProp(itype)%nPartition
      do j=1,nPartition
        PosDiff= SpherePartitionCoord(j,itype)
        LagrangeP = Cposition+ PosDiff
        ic= floor(LagrangeP%x*rdx)+1
        jc= floor(LagrangeP%y*rdyUniform)+1
        kc= floor(LagrangeP%z*rdz)+1        
        if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 &
        .or. (jc-y1start(2))*(jc-y1end(2))>0) cycle
        nIBP= nIBP+1
        if(nIBP> mIBP) call Reallocate_IbpVar()
        IBP_idlocal(nIBP)= i+nlocalOld   !!! note here
        IBP_Pos(nIBP)= LagrangeP
#define   clc_Point_indxyz_IBPMove
#include "clc_Point_indxyz_inc.f90"
#undef    clc_Point_indxyz_IBPMove
      enddo
    ENDDO
    call MPI_REDUCE(nIBP,nIbpSum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      if(nIbpSum>numIBP) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp_imp1","sum(nIbp)>numIBP, WRONG!!!")
    endif
    nIBPOld= nIBP
  end subroutine PrepareIBM_interp_imp1

  !******************************************************************
  ! PrepareIBM_interp_imp2
  !******************************************************************
  subroutine PrepareIBM_interp_imp2()
    implicit none

    ! locals 
    integer:: i,j,itype,nPartition
    integer:: ic,jc,kc,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    type(real3):: Cposition,LagrangeP,PosDiff,LinvelP,RotVelP

    mlocalNEW= GPrtcl_list%mlocal
    allocate(PIBM_linVel(mlocalNEW), PIBM_rotVel(mlocalNEW), PIBM_pType(mlocalNEW))

    call PComm_IBM_recover_imp()
    call PComm_IBM_forward_imp2()
   
    !NOTE (Gong Zheng, 2020/07/01) : 
    ! [1] Here ONLY the IBP points within the Processor physical Domains will be considered.
    !       The IBP points outsides will be skipped. 
    ! [2] ONLY the uniform meshes near y- dir are used.
    nIBP= 0
    DO i=1,nlocalOld
      itype  = PIBM_pType(i)
      LinVelP= PIBM_linVel(i)
      RotVelP= PIBM_rotVel(i)
      Cposition = GPrtcl_PosOld(i)
      nPartition= PrtclIBMProp(itype)%nPartition
      do j=1,nPartition
        PosDiff= SpherePartitionCoord(j,itype)
        LagrangeP = Cposition+ PosDiff
        ic= floor(LagrangeP%x*rdx)+1
        jc= floor(LagrangeP%y*rdyUniform)+1
        kc= floor(LagrangeP%z*rdz)+1
        if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 &
        .or. (jc-y1start(2))*(jc-y1end(2))>0) cycle
        nIBP= nIBP+1
        if(nIBP> mIBP) call MainLog%CheckForError(ErrT_Abort," PrepareIBM2","nIBP WRONG 1 !!!")
        IBP_idlocal(nIBP)= i   !!! note here
        IBP_Pos(nIBP)= LagrangeP
        IBP_Vel(nIBP)= LinvelP+ (RotVelP .cross. PosDiff)
#define   clc_Point_indxyz_IBPMove
#include "clc_Point_indxyz_inc.f90"
#undef    clc_Point_indxyz_IBPMove
      enddo
    ENDDO
    
    DO i=1,nGhostIBM
      itype  = GhostP_pType(i)
      LinvelP= GhostP_linVel(i)
      RotVelP= GhostP_rotVel(i)
      Cposition = GhostP_PosR(i)
      nPartition= PrtclIBMProp(itype)%nPartition
      do j=1,nPartition
        PosDiff= SpherePartitionCoord(j,itype)
        LagrangeP = Cposition+ PosDiff
        ic= floor(LagrangeP%x*rdx)+1
        jc= floor(LagrangeP%y*rdyUniform)+1
        kc= floor(LagrangeP%z*rdz)+1
        if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 &
        .or. (jc-y1start(2))*(jc-y1end(2))>0) cycle
        nIBP= nIBP+1
        if(nIBP> mIBP) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp_imp2","nIBP WRONG 2 !!!")
        IBP_idlocal(nIBP)= i+nlocalOld   !!! note here
        IBP_Pos(nIBP)= LagrangeP
        IBP_Vel(nIBP)= LinvelP+ (RotVelP .cross. PosDiff)
#define   clc_Point_indxyz_IBPMove
#include "clc_Point_indxyz_inc.f90"
#undef    clc_Point_indxyz_IBPMove
      enddo
    ENDDO

    if(nIBP /= nIBPOld) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp_imp2","nIBP /= nIBPOld, WRONG!!!")
    deallocate(PIBM_linVel,PIBM_rotVel,PIBM_pType,GPrtcl_PosOld); allocate(GPrtcl_PosOld(GPrtcl_list%mlocal))
  end subroutine PrepareIBM_interp_imp2

  !******************************************************************
  ! InterpolationIBM_imp
  !******************************************************************
  subroutine InterpolationIBM_imp(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in) ::ux,uy,uz

    ! locals
    type(real3):: LagrangeP
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir
    real(RK),dimension(0:nDistribute)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    DO pid=1,nIBP
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "clc_InterpolateCoe_inc.f90"
#include "clc_Point_Interpolation_inc.f90"
      IBP_Vel(pid)=  real3(SumXDir,SumYDir,SumZDir)
    ENDDO

    ! Fixed IBP part
    DO pid=1,nIBPFix
      LagrangeP  = IBPFix_Pos(pid)      
      idxc_interp= IBPFix_indxyz(1,pid)
      idxp_interp= IBPFix_indxyz(2,pid)
      idyc_interp= IBPFix_indxyz(3,pid)
      idyp_interp= IBPFix_indxyz(4,pid)
      idzc_interp= IBPFix_indxyz(5,pid)
      idzp_interp= IBPFix_indxyz(6,pid)
#include "clc_InterpolateCoe_inc.f90"
#include "clc_Point_Interpolation_inc.f90"
      IBPFix_Vel(pid)=  real3(-SumXDir,-SumYDir,-SumZDir)
    ENDDO
  end subroutine InterpolationIBM_imp

  !******************************************************************
  ! SpreadIbpForce_imp
  !****************************************************************** 
  subroutine SpreadIbpForce_imp(VolForce_x,VolForce_y,VolForce_z,Spreadflag)
    implicit none
    integer,intent(in)::Spreadflag
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout) ::VolForce_x,VolForce_y,VolForce_z

    ! locals
    type(real3):: LagrangeP,IbpForce
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,VolRatio
    real(RK),dimension(0:nDistribute)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    integer::i,j,k,id,jd,kd,pid,itype,idlocal
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    VolForce_x=0.0_RK; VolForce_y=0.0_RK; VolForce_z=0.0_RK
    DO pid=1,nIBP
      idlocal    = IBP_idlocal(pid)
      if(idlocal> nlocalOld) then
        itype= GhostP_Ptype(idlocal- nlocalOld)
      else
        itype= GPrtcl_Ptype(idlocal)
      endif
      VolRatio   = IbpVolRatio(itype)
      LagrangeP  = IBP_Pos(pid)
      IbpForce   = IBP_Vel(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "clc_InterpolateCoe_inc.f90"
#include "spreadIBM_force_inc.f90"
    ENDDO

    IF(Spreadflag==2) THEN
      ! Fixed IBP Part
      DO pid=1,nIBPFix
        VolRatio   = IBPFix_VolRatio(pid)
        LagrangeP  = IBPFix_Pos(pid)
        IbpForce   = IBPFix_Vel(pid)
        idxc_interp= IBPFix_indxyz(1,pid)
        idxp_interp= IBPFix_indxyz(2,pid)
        idyc_interp= IBPFix_indxyz(3,pid)
        idyp_interp= IBPFix_indxyz(4,pid)
        idzc_interp= IBPFix_indxyz(5,pid)
        idzp_interp= IBPFix_indxyz(6,pid)
#include "clc_InterpolateCoe_inc.f90"
#include "spreadIBM_force_inc.f90"
      ENDDO
    ENDIF
    call Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)
  end subroutine SpreadIbpForce_imp

  !******************************************************************
  ! UpdateRhsIBM_imp1
  !******************************************************************
  subroutine UpdateRhsIBM_imp1(RhsX,RhsY,RhsZ,VolForce_x,VolForce_y,VolForce_z)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX,RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ 
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::VolForce_x,VolForce_y,VolForce_z

    ! locals
    integer:: ic,jc,kc

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
           RhsX(ic,jc,kc)= RhsX(ic,jc,kc)- VolForce_x(ic,jc,kc)
           RhsY(ic,jc,kc)= RhsY(ic,jc,kc)- VolForce_y(ic,jc,kc)
           RhsZ(ic,jc,kc)= RhsZ(ic,jc,kc)- VolForce_z(ic,jc,kc)        
        enddo
      enddo
    ENDDO  
  end subroutine UpdateRhsIBM_imp1

  !******************************************************************
  ! UpdateRhsIBM_imp2
  !******************************************************************
  subroutine UpdateRhsIBM_imp2(RhsX,RhsY,RhsZ,VolForce_x,VolForce_y,VolForce_z)
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
  end subroutine UpdateRhsIBM_imp2

  !******************************************************************
  ! PComm_IBM_forward_imp1
  !******************************************************************
  subroutine PComm_IBM_forward_imp1()
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
      do i=1,nlocalOld
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
        endif
      enddo  
    ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) THEN
      nsendg=nsend
      do i=1,nlocalOld
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
      nrecv2=nrecv*GIBM_size_forward_imp
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward_imp
      allocate(buf_send(nsend2))
      call pack_IBM_forward_imp1(buf_send,nsendg,nsend,xp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),2,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(1),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward_imp1(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step2: send to xm_axis, and receive from xp_dir
    nsend=0; nrecv=0; ngp=ng; !ngpp=ng  
    IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,nlocalOld
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
        endif
      enddo 
    ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) THEN
      nsendg=nsend
      do i=1,nlocalOld
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
      nrecv2=nrecv*GIBM_size_forward_imp
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward_imp
      allocate(buf_send(nsend2))
      call pack_IBM_forward_imp1(buf_send,nsendg,nsend,xm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),4,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(2),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward_imp1(buf_recv,ngp+1,ng)
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
          GhostP_id(ng)     = nlocalOld+ i
          GhostP_pType(ng)  = GhostP_ptype(i)
          GhostP_Direction(ng)= THIS_PROC_ZP
          GhostP_PosR(ng)   = GhostP_PosR(i)
          GhostP_PosR(ng)%z = pz-zlz
        endif
      enddo

      do i=1,nlocalOld
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

      do i=1,nlocalOld
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
      nrecv2=nrecv*GIBM_size_forward_imp
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward_imp
      allocate(buf_send(nsend2))
      call pack_IBM_forward_imp1(buf_send,nsendg,nsend,zp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),6,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(3),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward_imp1(buf_recv,ngp+1,ng)
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
          GhostP_id(ng)     = nlocalOld+ i
          GhostP_pType(ng)  = GhostP_ptype(i)
          GhostP_Direction(ng)= THIS_PROC_ZM
          GhostP_PosR(ng)   = GhostP_PosR(i)
          GhostP_PosR(ng)%z = pz+ zlz
        endif
      enddo

      do i=1,nlocalOld
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

      do i=1,nlocalOld
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
      nrecv2=nrecv*GIBM_size_forward_imp
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward_imp
      allocate(buf_send(nsend2))
      call pack_IBM_forward_imp1(buf_send,nsendg,nsend,zm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),8,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(4),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward_imp1(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    nGhostIBM= ng
  end subroutine PComm_IBM_forward_imp1

  !**********************************************************************
  ! pack_IBM_forward_imp1
  !**********************************************************************
  subroutine pack_IBM_forward_imp1(buf_send,nsendg,nsend,dir)
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
      buf_send(m)=real(nlocalOld+ id);    m=m+1 ! 01
      buf_send(m)=real(GhostP_pType(id)); m=m+1 ! 02
      buf_send(m)=real(dir);              m=m+1 ! 03
      buf_send(m)=GhostP_PosR(id)%x+pdx;  m=m+1 ! 04
      buf_send(m)=GhostP_PosR(id)%y;      m=m+1 ! 05
      buf_send(m)=GhostP_PosR(id)%z+pdz;  m=m+1 ! 06
    enddo
    do i=nsendg+1,nsend
      id=sendlist(i)
      buf_send(m)=real(id);               m=m+1 ! 01
      buf_send(m)=real(GPrtcl_pType(id)); m=m+1 ! 02
      buf_send(m)=real(dir);              m=m+1 ! 03
      buf_send(m)=GPrtcl_PosR(id)%x+pdx;  m=m+1 ! 04
      buf_send(m)=GPrtcl_PosR(id)%y;      m=m+1 ! 05
      buf_send(m)=GPrtcl_PosR(id)%z+pdz;  m=m+1 ! 06
    enddo
  end subroutine pack_IBM_forward_imp1

  !**********************************************************************
  ! unpack_IBM_forward_imp1
  !**********************************************************************
  subroutine unpack_IBM_forward_imp1(buf_recv,n1,n2)
    implicit none
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::n1,n2

    ! locals
    integer::i,m,itype
   
    m=1
    do i=n1,n2
      GhostP_id(i)       =nint(buf_recv(m)); m=m+1 ! 01 
      itype              =nint(buf_recv(m)); m=m+1 ! 02
      GhostP_Direction(i)=nint(buf_recv(m)); m=m+1 ! 03
      GhostP_PosR(i)%x   =buf_recv(m);       m=m+1 ! 04
      GhostP_PosR(i)%y   =buf_recv(m);       m=m+1 ! 05
      GhostP_PosR(i)%z   =buf_recv(m);       m=m+1 ! 06
      GhostP_pType(i)    =itype;                
      GhostP_PosR(i)%w   =DEMProperty%Prtcl_PureProp(itype)%Radius 
    enddo
  end subroutine unpack_IBM_forward_imp1

  !**********************************************************************
  ! pack_IBM_forward_imp2
  !**********************************************************************
  subroutine pack_IBM_forward_imp2(buf_send,nsendg,nsend,dir)
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
      buf_send(m)=real(nlocalOld+ id);     m=m+1 ! 01
      buf_send(m)=real(GhostP_pType(id));  m=m+1 ! 02
      buf_send(m)=real(dir);               m=m+1 ! 03
      buf_send(m)=GhostP_PosR(id)%x+pdx;   m=m+1 ! 04
      buf_send(m)=GhostP_PosR(id)%y;       m=m+1 ! 05
      buf_send(m)=GhostP_PosR(id)%z+pdz;   m=m+1 ! 06
      buf_send(m)=GhostP_linVel(id)%x;     m=m+1 ! 07
      buf_send(m)=GhostP_linVel(id)%y;     m=m+1 ! 08
      buf_send(m)=GhostP_linVel(id)%z;     m=m+1 ! 09
      buf_send(m)=GhostP_rotVel(id)%x;     m=m+1 ! 10
      buf_send(m)=GhostP_rotVel(id)%y;     m=m+1 ! 11
      buf_send(m)=GhostP_rotVel(id)%z;     m=m+1 ! 12 
    enddo
    do i=nsendg+1,nsend
      id=sendlist(i)
      buf_send(m)=real(id);                m=m+1 ! 01
      buf_send(m)=real(PIBM_pType(id));    m=m+1 ! 02
      buf_send(m)=real(dir);               m=m+1 ! 03
      buf_send(m)=GPrtcl_PosOld(id)%x+pdx; m=m+1 ! 04
      buf_send(m)=GPrtcl_PosOld(id)%y;     m=m+1 ! 05
      buf_send(m)=GPrtcl_PosOld(id)%z+pdz; m=m+1 ! 06
      buf_send(m)=PIBM_linVel(id)%x;       m=m+1 ! 07
      buf_send(m)=PIBM_linVel(id)%y;       m=m+1 ! 08
      buf_send(m)=PIBM_linVel(id)%z;       m=m+1 ! 09
      buf_send(m)=PIBM_rotVel(id)%x;       m=m+1 ! 10
      buf_send(m)=PIBM_rotVel(id)%y;       m=m+1 ! 11
      buf_send(m)=PIBM_rotVel(id)%z;       m=m+1 ! 12 
    enddo
  end subroutine pack_IBM_forward_imp2

  !**********************************************************************
  ! unpack_IBM_forward_imp2
  !**********************************************************************
  subroutine unpack_IBM_forward_imp2(buf_recv,n1,n2)
    implicit none
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::n1,n2

    ! locals
    integer::i,m,itype
   
    m=1
    do i=n1,n2
      GhostP_id(i)       =nint(buf_recv(m)); m=m+1 ! 01 
      itype              =nint(buf_recv(m)); m=m+1 ! 02
      GhostP_Direction(i)=nint(buf_recv(m)); m=m+1 ! 03
      GhostP_PosR(i)%x   =buf_recv(m);       m=m+1 ! 04
      GhostP_PosR(i)%y   =buf_recv(m);       m=m+1 ! 05
      GhostP_PosR(i)%z   =buf_recv(m);       m=m+1 ! 06
      GhostP_linVel(i)%x =buf_recv(m);       m=m+1 ! 07 
      GhostP_linVel(i)%y =buf_recv(m);       m=m+1 ! 08 
      GhostP_linVel(i)%z =buf_recv(m);       m=m+1 ! 09 
      GhostP_rotVel(i)%x =buf_recv(m);       m=m+1 ! 10 
      GhostP_rotVel(i)%y =buf_recv(m);       m=m+1 ! 11 
      GhostP_rotVel(i)%z =buf_recv(m);       m=m+1 ! 12
      GhostP_pType(i)    =itype;                
      GhostP_PosR(i)%w   =DEMProperty%Prtcl_PureProp(itype)%Radius 
    enddo
  end subroutine unpack_IBM_forward_imp2

  !******************************************************************
  ! IntegrateFluidPrtclForce_imp
  !******************************************************************
  subroutine IntegrateFluidPrtclForce_imp(iCountACM,ux,uy,uz)
    implicit none
    integer,intent(in)::iCountACM
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz

    ! locals
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    type(real3)::FpVolume,Cposition,PosDiff,CenterDiff,CellVel
    real(RK)::CenterDistY,CenterDistZ,CenterDist,Radius,VoidFraction,MassL,InertiaL
    integer::i,j,k,pid,idlocal,itype,idGhost,nsend,nrecv,nsend2,nrecv2,ierror,request(4)
    integer::DomainIx1,DomainIx2,DomainIy1,DomainIy2,DomainIz1,DomainIz2,SRstatus(MPI_STATUS_SIZE)
    
    ! In order to save memory, GhostP_linVel is used to story FluidPrtclForce for ghost particle temporarily.
    ! Similarily, GhostP_rotVel is used to story FluidPrtclTorque for ghost particle temporarily.
    GPrtcl_FpForce = zero_r3
    GPrtcl_FpTorque= zero_r3
    GhostP_linVel  = zero_r3
    GhostP_rotVel  = zero_r3
    DO pid=1,nIBP
      idlocal    = IBP_idlocal(pid)
      if(idlocal> nlocalOld) then
        idGhost= idlocal- nlocalOld
        itype= GhostP_Ptype(idGhost)
        FpVolume= IBP_Vel(pid)*PrtclIBMProp(itype)%IBPVolume
        GhostP_linVel(idGhost)=  GhostP_linVel(idGhost)+  FpVolume

        Cposition = GhostP_PosR(idGhost)
        PosDiff= IBP_Pos(pid)- Cposition
        GhostP_rotVel(idGhost)=  GhostP_rotVel(idGhost)+  (PosDiff  .cross. FpVolume)
      else
        itype= GPrtcl_Ptype(idlocal)
        FpVolume= IBP_Vel(pid)*PrtclIBMProp(itype)%IBPVolume
        GPrtcl_FpForce(idlocal) = GPrtcl_FpForce(idlocal)+ FpVolume

        Cposition = GPrtcl_PosR(idlocal)
        PosDiff= IBP_Pos(pid)- Cposition
        GPrtcl_FpTorque(idlocal)= GPrtcl_FpTorque(idlocal)+ (PosDiff  .cross. FpVolume)
      endif
    ENDDO

    ! Here I assume the velocities outside the domain are always 0.0_RK.
    GPrtcl_FluidIntegrate = zero_r3
    GhostP_FluidIntegrate = zero_r3
    DO pid=1,nlocalOld
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
              CellVel%x= 0.5_RK*(ux(i+1,j,k)+ ux(i,j,k))
              CellVel%y= 0.5_RK*(uy(i,j+1,k)+ uy(i,j,k))
              CellVel%z= 0.5_RK*(uz(i,j,k+1)+ uz(i,j,k))
              GPrtcl_FluidIntegrate(1,pid)=GPrtcl_FluidIntegrate(1,pid)+CellVel
              GPrtcl_FluidIntegrate(2,pid)=GPrtcl_FluidIntegrate(2,pid)+(CenterDiff .cross. CellVel)
            elseif(CenterDist+CellRadius > 0.0_RK) then
              VoidFraction= clc_VoidFraction_IBM(CenterDiff,Radius)
              CellVel%x= 0.5_RK*(ux(i+1,j,k)+ ux(i,j,k))
              CellVel%y= 0.5_RK*(uy(i,j+1,k)+ uy(i,j,k))
              CellVel%z= 0.5_RK*(uz(i,j,k+1)+ uz(i,j,k))
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
              CellVel%x= 0.5_RK*(ux(i+1,j,k)+ ux(i,j,k))
              CellVel%y= 0.5_RK*(uy(i,j+1,k)+ uy(i,j,k))
              CellVel%z= 0.5_RK*(uz(i,j,k+1)+ uz(i,j,k))
              GhostP_FluidIntegrate(1,pid)=GhostP_FluidIntegrate(1,pid)+CellVel
              GhostP_FluidIntegrate(2,pid)=GhostP_FluidIntegrate(2,pid)+(CenterDiff .cross. CellVel)
            elseif(CenterDist+CellRadius > 0.0_RK) then
              VoidFraction= clc_VoidFraction_IBM(CenterDiff,Radius)
              CellVel%x= 0.5_RK*(ux(i+1,j,k)+ ux(i,j,k))
              CellVel%y= 0.5_RK*(uy(i,j+1,k)+ uy(i,j,k))
              CellVel%z= 0.5_RK*(uz(i,j,k+1)+ uz(i,j,k))
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
        if(idlocal> nlocalOld) then
          idGhost= idlocal- nlocalOld
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
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 9, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 9, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),10,MPI_COMM_WORLD,request(1),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),10,MPI_COMM_WORLD,ierror)
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
        if(idlocal> nlocalOld) then
          idGhost= idlocal- nlocalOld
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
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 11, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 11, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),12,MPI_COMM_WORLD,request(2),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),12,MPI_COMM_WORLD,ierror)
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
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 13, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 13, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),14,MPI_COMM_WORLD,request(3),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),14,MPI_COMM_WORLD,ierror)
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
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 15, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 15, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),16,MPI_COMM_WORLD,request(4),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),16,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(4),SRstatus,ierror)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Finally, get the final GPrtcl_FpForce
    if(iCountACM==1) then
      do i=1, nlocalOld
        itype   = GPrtcl_pType(i)
        MassL   = PrtclIBMProp(itype)%MassIBL
        InertiaL= PrtclIBMProp(itype)%InertiaIBL
        GPrtcl_FpForce(i) =(FluidDensity*(zero_r3 +GPrtcl_FpForce(i)) -MassL   *GPrtcl_linVel(1,i))/pmAlpha
        GPrtcl_FpTorque(i)=(FluidDensity*(zero_r3 +GPrtcl_FpTorque(i))-InertiaL*GPrtcl_rotVel(1,i))/pmAlpha
        GPrtcl_FluidIntOld(1,i)=GPrtcl_FluidIntegrate(1,i)
        GPrtcl_FluidIntOld(2,i)=GPrtcl_FluidIntegrate(2,i)
      enddo
     else
      do i=1, nlocalOld
        itype   = GPrtcl_pType(i)
        MassL   = PrtclIBMProp(itype)%MassIBL
        InertiaL= PrtclIBMProp(itype)%InertiaIBL
        GPrtcl_FpForce(i) =(FluidDensity*(GPrtcl_FluidIntegrate(1,i)-GPrtcl_FluidIntOld(1,i)+GPrtcl_FpForce(i)) -MassL   *GPrtcl_linVel(1,i))/pmAlpha
        GPrtcl_FpTorque(i)=(FluidDensity*(GPrtcl_FluidIntegrate(2,i)-GPrtcl_FluidIntOld(2,i)+GPrtcl_FpTorque(i))-InertiaL*GPrtcl_rotVel(1,i))/pmAlpha
        GPrtcl_FluidIntOld(1,i)=GPrtcl_FluidIntegrate(1,i)
        GPrtcl_FluidIntOld(2,i)=GPrtcl_FluidIntegrate(2,i)
      enddo
    endif
  end subroutine IntegrateFluidPrtclForce_imp

  !**********************************************************************
  ! pack_IBM_FpForce
  !**********************************************************************
  subroutine pack_IBM_FpForce(buf_send,nsend)
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
  end subroutine pack_IBM_FpForce

  !**********************************************************************
  ! unpack_IBM_FpForce
  !**********************************************************************
  subroutine unpack_IBM_FpForce(buf_recv,nrecv)
    implicit none
    real(RK),dimension(:),intent(in):: buf_recv
    integer,intent(in)::nrecv

    !locals
    type(real3)::FpForce,FpTorque,FVelInt,FTorInt
    integer::i,pid,m,gid

    m=1
    do i=1,nrecv
      pid= nint(buf_recv(m)); m=m+1
      FpForce%x =buf_recv(m); m=m+1
      FpForce%y =buf_recv(m); m=m+1
      FpForce%z =buf_recv(m); m=m+1
      FpTorque%x=buf_recv(m); m=m+1
      FpTorque%y=buf_recv(m); m=m+1
      FpTorque%z=buf_recv(m); m=m+1 
      FVelInt%x =buf_recv(m); m=m+1
      FVelInt%y =buf_recv(m); m=m+1
      FVelInt%z =buf_recv(m); m=m+1
      FTorInt%x =buf_recv(m); m=m+1
      FTorInt%y =buf_recv(m); m=m+1
      FTorInt%z =buf_recv(m); m=m+1   
      if(pid>nlocalOld) then
        gid= pid-nlocalOld
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
  end subroutine unpack_IBM_FpForce

  !******************************************************************
  ! PComm_IBM_recover_imp
  !******************************************************************
  subroutine PComm_IBM_recover_imp()
    implicit none

    ! locals
    real(RK)::px,pz
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer:: i,ierror,nlocalNEW,nsend,nrecv,nsend2,nrecv2,request(4),nlocalp,itype,SRstatus(MPI_STATUS_SIZE)
    
    ! GPrtcl_PosOld, PIBM_pType, PIBM_linVel, PIBM_rotVel
    nlocalNEW= GPrtcl_list%nlocal
    do i=1,nlocalNEW
      itype= GPrtcl_pType(i)
      PIBM_pType(i) = itype
      PIBM_linVel(i)= GPrtcl_linVel(1,i)
      PIBM_rotVel(i)= GPrtcl_rotVel(1,i)
      GPrtcl_PosOld(i)%w= DEMProperty%Prtcl_PureProp(itype)%Radius
    enddo

    ! Step1: send to xp_axis, and receive from xm_dir
    nsend =0; nrecv =0
    IF(DEM_decomp%ProcNgh(3) /= nrank .and. DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) THEN
      do i=1,nlocalNEW
        px=GPrtcl_PosOld(i)%x
        if(px >= xedCoord) then
          nsend= nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 17, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 17, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_recover_imp
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),18,MPI_COMM_WORLD,request(1),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_recover_imp
      allocate(buf_send(nsend2))
      nlocalp= nlocalNEW;  nlocalNEW = nlocalp- nsend
      call pack_IBM_recover_imp(buf_send,nsend,nlocalp,xp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),18,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(1),SRstatus,ierror)
      nlocalp= nlocalNEW;  nlocalNEW= nlocalp+ nrecv
      if(nlocalNEW>mlocalNEW) call Reallocate_PIBMOld(nlocalNEW)
      call unpack_IBM_recover_imp(buf_recv,nrecv,nlocalp,xp_dir)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! step2: send to xm_axis, and receive from xp_dir
    nsend =0; nrecv =0
    IF(DEM_decomp%ProcNgh(4) /= nrank .and. DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) THEN
      do i=1,nlocalNEW
        px=GPrtcl_PosOld(i)%x
        if(px < xstCoord) then
          nsend= nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 19, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 19, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_recover_imp
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),20,MPI_COMM_WORLD,request(2),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_recover_imp
      allocate(buf_send(nsend2))
      nlocalp= nlocalNEW;  nlocalNEW = nlocalp- nsend
      call pack_IBM_recover_imp(buf_send,nsend,nlocalp,xm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),20,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(2),SRstatus,ierror)
      nlocalp= nlocalNEW;  nlocalNEW= nlocalp+ nrecv
      if(nlocalNEW>mlocalNEW) call Reallocate_PIBMOld(nlocalNEW)
      call unpack_IBM_recover_imp(buf_recv,nrecv,nlocalp,xm_dir)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)  

    ! Step3: send to zp_axis, and receive from zm_dir
    nsend =0; nrecv =0
    IF(DEM_decomp%ProcNgh(1) /= nrank .and. DEM_decomp%ProcNgh(1) /= MPI_PROC_NULL) THEN
      do i=1,nlocalNEW
        pz=GPrtcl_PosOld(i)%z
        if(pz >= zedCoord) then
          nsend= nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 21, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 21, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_recover_imp
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),22,MPI_COMM_WORLD,request(3),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_recover_imp
      allocate(buf_send(nsend2))
      nlocalp= nlocalNEW;  nlocalNEW = nlocalp- nsend
      call pack_IBM_recover_imp(buf_send,nsend,nlocalp,zp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),22,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(3),SRstatus,ierror)
      nlocalp= nlocalNEW;  nlocalNEW= nlocalp+ nrecv
      if(nlocalNEW>mlocalNEW) call Reallocate_PIBMOld(nlocalNEW)
      call unpack_IBM_recover_imp(buf_recv,nrecv,nlocalp,zp_dir)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv) 

    ! Step4: send to zm_axis, and receive from zp_dir
    nsend =0; nrecv =0
    IF(DEM_decomp%ProcNgh(2) /= nrank .and. DEM_decomp%ProcNgh(2) /= MPI_PROC_NULL) THEN
      do i=1,nlocalNEW
        pz=GPrtcl_PosOld(i)%z
        if(pz < zstCoord) then
          nsend= nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 23, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 23, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_recover_imp
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),24,MPI_COMM_WORLD,request(4),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_recover_imp
      allocate(buf_send(nsend2))
      nlocalp= nlocalNEW;  nlocalNEW = nlocalp- nsend
      call pack_IBM_recover_imp(buf_send,nsend,nlocalp,zm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),24,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(4),SRstatus,ierror)
      nlocalp= nlocalNEW;  nlocalNEW= nlocalp+ nrecv
      if(nlocalNEW>mlocalNEW) call Reallocate_PIBMOld(nlocalNEW)
      call unpack_IBM_recover_imp(buf_recv,nrecv,nlocalp,zm_dir)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    if(nlocalNEW /= nlocalOld) call MainLog%CheckForError(ErrT_Abort," PComm_IBM_recover_imp","nlocalNEW /= nlocalOld, WRONG!!!")
  end subroutine PComm_IBM_recover_imp

  !******************************************************************
  ! PComm_IBM_forward_imp2
  !******************************************************************
  subroutine PComm_IBM_forward_imp2()
    implicit none

    ! locals
    real(RK)::px,pxt,pz,pzt
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer::nsend,nsend2,nrecv,nrecv2,nsendg,ng,ngp,ngpp
    integer::i,ierror,request(4),SRstatus(MPI_STATUS_SIZE)
    
    ! This part is similar to what is done in subroutine PC_Comm_For_Cntct of file "ACM_comm.f90"(only y_pencil)
    ng=0

    ! Step1: send to xp_axis, and receive from xm_dir
    nsend=0; nrecv=0; ngp=ng; ngpp=ng
    IF(DEM_decomp%ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,nlocalOld
        px = GPrtcl_PosOld(i)%x
        pxt= px+ GPrtcl_PosOld(i)%w
        if(pxt +SMALL> xedCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_id(ng)       = i
          GhostP_pType(ng)    = PIBM_ptype(i)
          GhostP_Direction(ng)= THIS_PROC_XP
          GhostP_PosR(ng)     = GPrtcl_PosOld(i)
          GhostP_PosR(ng)%x   = px-xlx
          GhostP_linVel(ng)   = PIBM_linVel(i)
          GhostP_rotVel(ng)   = PIBM_rotVel(i)
        endif
      enddo  
    ELSEIF(DEM_decomp%ProcNgh(3) /= MPI_PROC_NULL) THEN
      nsendg=nsend
      do i=1,nlocalOld
        pxt= GPrtcl_PosOld(i)%x +GPrtcl_PosOld(i)%w
        if(pxt+SMALL >xedCoord) then
          nsend=nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(3), 25, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(4), 25, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(4),26,MPI_COMM_WORLD,request(1),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward_imp2(buf_send,nsendg,nsend,xp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(3),26,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(1),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward_imp2(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step2: send to xm_axis, and receive from xp_dir
    nsend=0; nrecv=0; ngp=ng; !ngpp=ng   
    IF(DEM_decomp%ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,nlocalOld
        px = GPrtcl_PosOld(i)%x
        pxt= px- GPrtcl_PosOld(i)%w
        if(pxt-SMALL<xstCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_id(ng)     = i
          GhostP_pType(ng)  = PIBM_pType(i)
          GhostP_Direction(ng)= THIS_PROC_XM
          GhostP_PosR(ng)   = GPrtcl_PosOld(i)
          GhostP_PosR(ng)%x = px+xlx
          GhostP_linVel(ng) = PIBM_linVel(i)
          GhostP_rotVel(ng) = PIBM_rotVel(i)
        endif
      enddo 
    ELSEIF(DEM_decomp%ProcNgh(4) /= MPI_PROC_NULL) THEN
      nsendg=nsend
      do i=1,nlocalOld
        pxt= GPrtcl_PosOld(i)%x- GPrtcl_PosOld(i)%w
        if(pxt-SMALL<xstCoord) then
          nsend=nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(4), 27, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(3), 27, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(3),28,MPI_COMM_WORLD,request(2),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward_imp2(buf_send,nsendg,nsend,xm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(4),28,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(2),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward_imp2(buf_recv,ngp+1,ng)
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
          GhostP_id(ng)     = nlocalOld+ i
          GhostP_pType(ng)  = GhostP_ptype(i)
          GhostP_Direction(ng)= THIS_PROC_ZP
          GhostP_PosR(ng)   = GhostP_PosR(i)
          GhostP_PosR(ng)%z = pz-zlz
          GhostP_linVel(ng) = GhostP_linVel(i)
          GhostP_rotVel(ng) = GhostP_rotVel(i)
        endif
      enddo

      do i=1,nlocalOld
        pz = GPrtcl_PosOld(i)%z 
        pzt= pz+ GPrtcl_PosOld(i)%w
        if(pzt+SMALL>zedCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_id(ng)     = i
          GhostP_pType(ng)  = PIBM_pType(i)
          GhostP_Direction(ng)= THIS_PROC_ZP
          GhostP_PosR(ng)   = GPrtcl_PosOld(i)
          GhostP_PosR(ng)%z = pz-zlz
          GhostP_linVel(ng) = PIBM_linVel(i)
          GhostP_rotVel(ng) = PIBM_rotVel(i)
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

      do i=1,nlocalOld
        pzt=GPrtcl_PosOld(i)%z+ GPrtcl_PosOld(i)%w
        if(pzt+SMALL>zedCoord) then
          nsend=nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(1), 29, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(2), 29, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(2),30,MPI_COMM_WORLD,request(3),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward_imp2(buf_send,nsendg,nsend,zp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(1),30,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(3),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward_imp2(buf_recv,ngp+1,ng)
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
          GhostP_id(ng)     = nlocalOld+ i
          GhostP_pType(ng)  = GhostP_ptype(i)
          GhostP_Direction(ng)= THIS_PROC_ZM
          GhostP_PosR(ng)   = GhostP_PosR(i)
          GhostP_PosR(ng)%z = pz+ zlz
          GhostP_linVel(ng) = GhostP_linVel(i)
          GhostP_rotVel(ng) = GhostP_rotVel(i)
        endif
      enddo

      do i=1,nlocalOld
        pz = GPrtcl_PosOld(i)%z 
        pzt= pz- GPrtcl_PosOld(i)%w
        if(pzt-SMALL<zstCoord ) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
          GhostP_id(ng)     = i
          GhostP_pType(ng)  = PIBM_ptype(i)
          GhostP_Direction(ng)= THIS_PROC_ZM
          GhostP_PosR(ng)   = GPrtcl_PosOld(i)
          GhostP_PosR(ng)%z = pz+ zlz
          GhostP_linVel(ng) = PIBM_linVel(i)
          GhostP_rotVel(ng) = PIBM_rotVel(i)
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

      do i=1,nlocalOld
        pzt= GPrtcl_PosOld(i)%z- GPrtcl_PosOld(i)%w
        if(pzt-SMALL<zstCoord) then
          nsend=nsend+1
          if(nsend > DEM_Comm%msend) call DEM_Comm%reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, int_type, DEM_decomp%ProcNgh(2), 31, &
                      nrecv, 1, int_type, DEM_decomp%ProcNgh(1), 31, MPI_COMM_WORLD,SRstatus,ierror)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,DEM_decomp%ProcNgh(1),32,MPI_COMM_WORLD,request(4),ierror)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward_imp2(buf_send,nsendg,nsend,zm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,DEM_decomp%ProcNgh(2),32,MPI_COMM_WORLD,ierror)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(4),SRstatus,ierror)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng)
      call unpack_IBM_forward_imp2(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    nGhostIBM= ng
  end subroutine PComm_IBM_forward_imp2


  !******************************************************************
  ! pack_IBM_recover_imp
  !****************************************************************** 
  subroutine pack_IBM_recover_imp(buf_send,nsend,nlocalp,dir)
    implicit none
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsend,dir
    integer,intent(inout)::nlocalp

    ! locals
    integer::i,id,m
    real(RK)::pdx,pdz

    m=1
    pdx=dx_pbc(dir)
    pdz=dz_pbc(dir)
    DO i=1,nsend
      id = sendlist(i)
      buf_send(m)=GPrtcl_PosOld(id)%x+pdx; m=m+1 ! 01
      buf_send(m)=GPrtcl_PosOld(id)%y;     m=m+1 ! 02
      buf_send(m)=GPrtcl_PosOld(id)%z+pdz; m=m+1 ! 03
      buf_send(m)=real(PIBM_pType(id));    m=m+1 ! 04
      buf_send(m)=PIBM_linVel(id)%x;       m=m+1 ! 05
      buf_send(m)=PIBM_linVel(id)%y;       m=m+1 ! 06
      buf_send(m)=PIBM_linVel(id)%z;       m=m+1 ! 07
      buf_send(m)=PIBM_rotVel(id)%x;       m=m+1 ! 08
      buf_send(m)=PIBM_rotVel(id)%y;       m=m+1 ! 09
      buf_send(m)=PIBM_rotVel(id)%z;       m=m+1 ! 10
    ENDDO
    DO i=nsend,1,-1
      id = sendlist(i)
      GPrtcl_PosOld(id)= GPrtcl_PosOld(nlocalp)
      PIBM_pType(id)   = PIBM_pType(nlocalp)
      PIBM_linVel(id)  = PIBM_linVel(nlocalp)
      PIBM_rotVel(id)  = PIBM_rotVel(nlocalp)
      nlocalp = nlocalp -1
    ENDDO
  end subroutine pack_IBM_recover_imp

  !******************************************************************
  ! unpack_IBM_recover_imp
  !******************************************************************
  subroutine unpack_IBM_recover_imp(buf_recv,nrecv,nlocalp,dir)
    implicit none
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in):: nrecv,nlocalp,dir

    ! locals
    logical:: IsInDomain
    integer::i,id,m,mp,itype

    m=1
    id=nlocalp+1
    DO i=1,nrecv
      IsInDomain= .true.
      SELECT CASE(dir)    
      CASE(xp_dir)
        mp = m
        if(buf_recv(mp)>=xedCoord) IsInDomain = .false.
      CASE(xm_dir)
        mp = m
        if(buf_recv(mp)< xstCoord) IsInDomain = .false.
      CASE(zp_dir)
        mp = m+2
        if(buf_recv(mp)>=zedCoord) IsInDomain = .false.
      CASE(zm_dir)
        mp = m+2
        if(buf_recv(mp)< zstCoord) IsInDomain = .false.
      END SELECT
      IF( .not. IsInDomain) THEN
        call MainLog%CheckForError(ErrT_Abort,"unpack_IBM_recover_imp","particle coordinate WRONG !!!")
      ENDIF

      GPrtcl_PosOld(id)%x= buf_recv(m); m=m+1   ! 01
      GPrtcl_PosOld(id)%y= buf_recv(m); m=m+1   ! 02
      GPrtcl_PosOld(id)%z= buf_recv(m); m=m+1   ! 03
      itype= nint(buf_recv(m));         m=m+1   ! 04
      PIBM_linVel(id)%x  = buf_recv(m); m=m+1   ! 05
      PIBM_linVel(id)%y  = buf_recv(m); m=m+1   ! 06
      PIBM_linVel(id)%z  = buf_recv(m); m=m+1   ! 07
      PIBM_rotVel(id)%x  = buf_recv(m); m=m+1   ! 08
      PIBM_rotVel(id)%y  = buf_recv(m); m=m+1   ! 09
      PIBM_rotVel(id)%z  = buf_recv(m); m=m+1   ! 10
      PIBM_pType(id)= itype
      GPrtcl_PosOld(id)%w= DEMProperty%Prtcl_PureProp(itype)%Radius
      id= id+1
    ENDDO
  end subroutine unpack_IBM_recover_imp

  !******************************************************************
  ! Reallocate_ghost_for_IBM
  !******************************************************************
  subroutine Reallocate_ghost_for_IBM(ng)
    implicit none
    integer,intent(in)::ng

    ! locals
    integer::sizep,ierrTemp,ierror=0
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
    !call MainLog%CheckForError(ErrT_Pass,"Reallocate_ghost_for_IBM"," Need to reallocate Ghost variables for IBM")
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

    if(ierror/=0) then
      call MainLog%CheckForError(ErrT_Abort," Reallocate_IbpVar"," Reallocate wrong!")
      call MainLog%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif    
    !call MainLog%CheckForError(ErrT_Pass,"Reallocate_IbpVar"," Need to reallocate IBP variables")
    !call MainLog%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    !call MainLog%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    !call MainLog%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine Reallocate_IbpVar

  !******************************************************************
  ! Reallocate_PIBMOld
  !******************************************************************
  subroutine Reallocate_PIBMOld(nlocalNEW)
    implicit none
    integer,intent(in)::nlocalNEW
    
    ! locals
    integer:: sizep,sizen,ierrTemp,ierror=0
    integer,dimension(:),allocatable:: IntVec
    type(real3),dimension(:),allocatable:: Real3Vec  
    type(real4),dimension(:),allocatable:: Real4Vec 

    sizep= mlocalNEW
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen, nlocalNEW)
    sizen= min(sizen, DEM_Opt%numPrtcl+10)
    mlocalNEW = sizen
    
    ! ======= integer vector part =======
    call move_alloc(PIBM_pType,IntVec)
    allocate(PIBM_pType(sizen),stat=ierrTemp);ierror=ierror+abs(ierrTemp)
    PIBM_pType(1:sizep)=IntVec
    deallocate(IntVec)  

    ! ======= real3 vercor part =======
    call move_alloc(PIBM_linVel,Real3Vec)
    allocate(PIBM_linVel(sizen),stat=ierrTemp);ierror=ierror+abs(ierrTemp)
    PIBM_linVel(1:sizep)=Real3Vec

    call move_alloc(PIBM_rotVel,Real3Vec)
    allocate(PIBM_rotVel(sizen),stat=ierrTemp);ierror=ierror+abs(ierrTemp)
    PIBM_rotVel(1:sizep)=Real3Vec
    deallocate(Real3Vec)

    ! ======= real4 vercor part =======  
    call move_alloc(GPrtcl_PosOld,Real4Vec)
    allocate(GPrtcl_PosOld(sizen),stat=ierrTemp);ierror=ierror+abs(ierrTemp)
    GPrtcl_PosOld(1:sizep)=Real4Vec
    deallocate(Real4Vec)

    if(ierror/=0) then
      call MainLog%CheckForError(ErrT_Abort," Reallocate_PIBMOld"," Reallocate wrong!")
      call MainLog%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif
    !call MainLog%CheckForError(ErrT_Pass,"Reallocate_PIBMOld"," Need to reallocate PIBPOld variables")
    !call MainLog%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    !call MainLog%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    !call MainLog%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine Reallocate_PIBMOld
#undef nDistribute
#undef THIS_PROC_XP
#undef THIS_PROC_XM
#undef THIS_PROC_YP
#undef THIS_PROC_YM
#undef THIS_PROC_ZP
#undef THIS_PROC_ZM
end module ca_IBM_implicit
