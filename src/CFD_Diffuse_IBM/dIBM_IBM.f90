#include "definitions_inc.f90"
module dIBM_IBM
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d,only:nrank,y1start,y1end,real_type
  use f2_Parameters
  use f2_Variables,only:mb1
  use f2_MeshAndMetries,only: dx,dyUniform,dz,rdx,rdyUniform,rdz,xc,yc,zc,dyp
  use dIBM_BC_and_Halo,only: Gather_Halo_IBMForce, SetBC_and_UpdateHalo_VolForce
  implicit none
  private
#ifdef IBMDistributeLinear
#define nDistribute 1
#else
#define nDistribute 2
#endif

  character(len=:),allocatable :: IBP_info_file_
  real(RK):: CellRadius,CellVolumeIBM,dxhalf,dyhalf,dzhalf,SMALL
  real(RK):: xstCoord,xedCoord,ystCoord,yedCoord,zstCoord,zedCoord

  ! Immersed boundary points parts
  integer:: nIBP         ! number of Immersed Boundary Point in local processor
  integer:: mIBP         ! the possible maxium # of Immersed Boundary Point in local processor
  real(RK),dimension(:),allocatable::    IBP_VolRatio
  integer(kind=2),dimension(:,:),allocatable::IBP_indxyz
  type(real3),dimension(:),allocatable:: IBP_Pos
  type(real3),dimension(:),allocatable:: IBP_Vel
  type(real3),dimension(:),allocatable:: IBP_Force
  type(real3),dimension(:),allocatable:: IBP_ForceAmplify

  ABSTRACT INTERFACE
    function deltaFunction_(ratio_in) result(delta)
      use mc_TypeDef,only:RK; implicit none
      real(RK),intent(in)::ratio_in
      real(RK):: delta
    end function deltaFunction_
  END INTERFACE
  procedure(deltaFunction_),pointer::deltaFunction
  procedure(),pointer::AdditionalForceIBM

  integer:: IBMForceScheme = 0
  integer:: nForcingExtra  = 2
  real(RK)::ForceAmplifyCoe= 1.0_RK
  
  ! public variables and functions/subroutines
  public:: nForcingExtra
  public:: Init_IBM, updateRhsIBM, prepareIBM_interp, Clc_NoSlipErr, AdditionalForceIBM
contains
#include "dIBM_IBM_common_inc.f90"

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Init_IBM
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Init_IBM(chFile); implicit none
    character(len=*),intent(in)::chFile
    
    ! locals
    character(512):: IBP_info_file
    integer:: ierr,ierrTmp,iUnit, deltaFunction_Scheme
    Namelist/IBMForce_Options/ IBMForceScheme, ForceAmplifyCoe, nforcingextra, deltaFunction_Scheme, IBP_info_file

    open(newunit=iUnit, file=chFile, status='old', form='formatted', IOSTAT=ierr)
    if(ierr /= 0) call MainLog%CheckForError(ErrT_Abort,"Init_IBM","Cannot open file:"//strip(chFile))
    read(iUnit, nml=IBMForce_Options)
    if(nrank==0) write(MainLog%iUnit, nml=IBMForce_Options)
    close(unit=iUnit, IOSTAT=ierr)
    IBP_info_file_ = strip(IBP_info_file)
    
    ! check integer(kind=2) is enough or not.
    if(nrank==0 .and. (nxc>huge(int(0,2))-2 .or. nyc>huge(int(0,2))-2 .or. nzc>huge(int(0,2))-2)) then
      call MainLog%CheckForError(ErrT_Abort,"Init_IBM","kind=2 is not enough for indxyz and IBPFix_indxyz")
    endif

#ifndef IBMDistributeLinear    
    if(deltaFunction_Scheme ==0) then
      deltaFunction => deltaFunction_Roma
    else
      deltaFunction => deltaFunction_Yang
    endif     
#endif

    select case(IBMForceScheme)
    case(0)  ! Kempe(2012,JCP), Gsell(2021,JCP)
      AdditionalForceIBM => AdditionalForceIBM_0
    case(1)  ! Zhao(2021,JCP), Cheylan(2023,JCP)
      AdditionalForceIBM => AdditionalForceIBM_1
    end select
        
    xstCoord= real(y1start(1)-1, kind=RK)*dx
    xedCoord= real(y1end(1),     kind=RK)*dx
    ystCoord= 0.0_RK
    yedCoord= yly
    zstCoord= real(y1start(3)-1, kind=RK)*dz
    zedCoord= real(y1end(3),     kind=RK)*dz

    SMALL= 1.0E-9_RK*dx
    
    !=========================== Immersed boundary points parts ===========================!
    CellRadius= 0.5_RK*sqrt(dx*dx+dyUniform*dyUniform+dz*dz)
    CellVolumeIBM= dx* dyUniform* dz
    dxhalf=dx*0.5_RK; dyhalf=dyUniform*0.5_RK; dzhalf=dz*0.5_RK
    mIBP = 1000

    allocate(IBP_VolRatio(mIBP),Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(IBP_indxyz(6,mIBP),Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(IBP_Pos(mIBP),     Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(IBP_Vel(mIBP),     Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(IBP_Force(mIBP),   Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"Init_IBM","Allocation failed 2")
    IBP_indxyz=0
    IBP_Pos=zero_r3; IBP_Vel=zero_r3; IBP_Force=zero_r3
  end subroutine Init_IBM

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PrepareIBM_interp
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PrepareIBM_interp(VolForce_x,VolForce_y,VolForce_z); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::VolForce_x,VolForce_y,VolForce_z
    
    ! locals
    type :: stl_ele
      sequence
      real(4) :: normal(3)
      real(4) :: p1(3)
      real(4) :: p2(3)
      real(4) :: p3(3)
      integer(2) :: marker = int(0, 2)
    endtype stl_ele
    type(stl_ele) :: element

    type(real3):: LagrangeP   
    real(RK) :: area, dxyz, p1(3), p2(3), p3(3)
    integer:: iUnit,ierr, ntri, itri, i,ic,jc,kc
    integer:: idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    
    dxyz  = CellVolumeIBM**(0.33333333333333333333_RK)
    
    !NOTE (Gong Zheng, 2020/07/01):
    ! [1] Here ONLY the IBP points within the Processor physical Domains will be considered.
    !       The IBP points outsides will be skipped.
    ! [2] ONLY the uniform meshes near y-dir are used.
    nIBP= 0
    open(newunit=iUnit, file= IBP_info_file_,form='unformatted',access='stream',status='old', iostat=ierr)
    read (unit = iUnit, pos = 81) ntri
    
    DO itri = 1, ntri
      read(unit=iUnit, pos = 35 + 50 * itri, iostat=ierr) element
      if (ierr /= 0) then
         call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp"," STL file wrong!")
      endif
      p1 = element%p1;  p2 = element%p2;  p3 = element%p3
      LagrangeP = (p1 + p2 + p3) / 3.0_RK
      call clc_area(p1, p2, p3, area)

      ic= floor(LagrangeP%x*rdx)+1
      jc= floor(LagrangeP%y*rdyUniform)+1
      kc= floor(LagrangeP%z*rdz)+1
      if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 &
        .or. (jc-y1start(2))*(jc-y1end(2))>0) cycle        
      nIBP= nIBP+1
      if(nIBP> mIBP) call Reallocate_IbpVar()
      IBP_Pos(nIBP)= LagrangeP
      IBP_Vel(nIBP)= zero_r3
      IBP_VolRatio(nIBP) = area/(dxyz *dxyz)
#define   clc_Point_indxyz_IBPMove
#include "dIBM_clc_Point_indxyz_inc.f90"
#undef    clc_Point_indxyz_IBPMove
    ENDDO
    close(unit=iUnit, iostat=ierr)
    do i=1,nIBP
      IBP_Force(i) =zero_r3
    enddo
    write(*, '(A, I0, A, I0, A)') '   rank ', nrank, ' has ', nIBP, ' IBP points.'
    
    IF(IBMForceScheme==1) call PrepareIBM_ForceAmplify(VolForce_x,VolForce_y,VolForce_z)
  end subroutine PrepareIBM_interp

subroutine clc_area(tp1, tp2, tp3, area)
   implicit none
   real(RK), intent(in) :: tp1(3), tp2(3), tp3(3)
   real(RK), intent(out) :: area
     
   ! locals
   real(RK) :: pt1(3), pt2(3), pt3(3)

   pt1 = tp2 - tp1;  pt2 = tp3 - tp1
   pt3(1) = pt1(2) * pt2(3) - pt1(3) * pt2(2)
   pt3(2) = pt1(3) * pt2(1) - pt1(1) * pt2(3)
   pt3(3) = pt1(1) * pt2(2) - pt1(2) * pt2(1)
   area = 0.5_RK * sqrt(pt3(1) * pt3(1) + pt3(2) * pt3(2) + pt3(3) * pt3(3))  
endsubroutine clc_area

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PrepareIBM_ForceAmplify
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PrepareIBM_ForceAmplify(VolForce_x,VolForce_y,VolForce_z); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::VolForce_x,VolForce_y,VolForce_z

    ! locals
    type(real3):: LagrangeP
    real(RK),dimension(0:nDistribute)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir,VolRatio
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    if(allocated(IBP_ForceAmplify)) deallocate(IBP_ForceAmplify)
    if(nIBP >0) allocate(IBP_ForceAmplify(nIBP))
    VolForce_x=0.0_RK; VolForce_y=0.0_RK; VolForce_z=0.0_RK
    
    ! step1: spread unit Lagrangian force
    DO pid=1,nIBP
      VolRatio   = IBP_VolRatio(pid)
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "dIBM_clc_InterpolateCoe_inc.f90"

      !=====   spreading   =====!
#define spreadIBM_force_inc_unity
#include "dIBM_spreadIBM_force_inc.f90"
#undef  spreadIBM_force_inc_unity
    ENDDO
    call Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)
    call SetBC_and_UpdateHalo_VolForce(VolForce_x,VolForce_y,VolForce_z)
    
    ! step 2: Interpolate unify force to Lagrangian point
    DO pid=1,nIBP
      VolRatio   = IBP_VolRatio(pid)
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "dIBM_clc_InterpolateCoe_inc.f90"

      !===== interpolation =====!
#define spreadIBM_force_inc_unity
#include "dIBM_clc_Point_Interpolation_inc.f90"
#undef  spreadIBM_force_inc_unity

      !=====    forcing    =====!
      IBP_ForceAmplify(pid)= real3(1.0_RK/SumXDir,1.0_RK/SumYDir,1.0_RK/SumZDir)
    ENDDO
  end subroutine PrepareIBM_ForceAmplify
        
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Clc_NoSlipErr
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Clc_NoSlipErr(ux,uy,uz); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz

    ! locals
    type(real3):: LagrangeP
    real(RK),dimension(0:nDistribute)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir,VolRatio
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    real(RK):: NoSlipErr(4), ErrTmp(3), ErrSum(4)
    
    ErrTmp=-1.0_RK; NoSlipErr=-1.0_RK; ErrSum=0.0_RK
    DO pid=1,nIBP
      VolRatio   = IBP_VolRatio(pid)
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "dIBM_clc_InterpolateCoe_inc.f90"

      !===== interpolation =====!
#include "dIBM_clc_Point_Interpolation_inc.f90"

      !=====    forcing    =====!
      ErrSum(1) = ErrSum(1) +abs(IBP_Vel(pid)%x -SumXDir)
      ErrSum(2) = ErrSum(2) +abs(IBP_Vel(pid)%y -SumYDir)
      ErrSum(3) = ErrSum(3) +abs(IBP_Vel(pid)%z -SumZDir)
      ErrSum(4) = ErrSum(4) +1.0_RK
      ErrTmp(1) = max(ErrTmp(1), abs(IBP_Vel(pid)%x -SumXDir))
      ErrTmp(2) = max(ErrTmp(2), abs(IBP_Vel(pid)%y -SumYDir))
      ErrTmp(3) = max(ErrTmp(3), abs(IBP_Vel(pid)%z -SumZDir))
    ENDDO
    call MPI_REDUCE(ErrTmp, NoSlipErr, 3, real_type, MPI_MAX, 0, MPI_COMM_WORLD, i)
    if(nrank==0) print*, itime,'^^^^^^^^^^^^ Max No-Slip Error = ', NoSlipErr(1:3)
    call MPI_REDUCE(ErrSum, NoSlipErr, 4, real_type, MPI_SUM, 0, MPI_COMM_WORLD, i)
    if(nrank==0) print*, itime,'^^^^^^^^^^^^ Ave No-Slip Error = ', NoSlipErr(1:3)/NoSlipErr(4)
  end subroutine Clc_NoSlipErr

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! updateRhsIBM
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine updateRhsIBM(uxStar,uyStar,uzStar,ux,uy,uz,RhsX,RhsY,RhsZ); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::uxStar,uyStar,uzStar
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in) ::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsX,RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    
    ! locals
    integer:: ic,jc,kc

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
           RhsX(ic,jc,kc)= uxStar(ic,jc,kc) -ux(ic,jc,kc)
           RhsY(ic,jc,kc)= uyStar(ic,jc,kc) -uy(ic,jc,kc)
           RhsZ(ic,jc,kc)= uzStar(ic,jc,kc) -uz(ic,jc,kc)        
        enddo
      enddo
    ENDDO  
  end subroutine updateRhsIBM

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! NormVolForce_x
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine NormVolForce_x(VolForce_x); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::VolForce_x
    
    ! locals
    integer::ic,jc,kc,ierr
    real(RK)::SumVolForceX,SumVolForceXTot,ForcedCoe
    
    SumVolForceX=0.0_RK
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        ForcedCoe=dyp(jc)
        do ic=y1start(1),y1end(1)
          SumVolForceX=SumVolForceX+ForcedCoe*VolForce_x(ic,jc,kc)
        enddo
      enddo
    ENDDO
    call MPI_ALLREDUCE(SumVolForceX,SumVolForceXTot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    SumVolForceX = -SumVolForceXTot/(real(nxc,RK)*real(nzc,RK))/yly
    VolForce_x =VolForce_x+SumVolForceX
    
    PrGradData(3)=PrGradData(3) +SumVolForceX/pmAlpha
    PrGradData(1)=PrGradData(1) +SumVolForceX/dt
  end subroutine NormVolForce_x
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! AdditionalForceIBM_0
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine AdditionalForceIBM_0(ux,uy,uz,VolForce_x,VolForce_y,VolForce_z); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,VolForce_x,VolForce_y,VolForce_z

    ! locals
    type(real3):: LagrangeP,IbpForce
    real(RK),dimension(0:nDistribute)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir,VolRatio
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    VolForce_x=0.0_RK; VolForce_y=0.0_RK; VolForce_z=0.0_RK
    DO pid=1,nIBP
      VolRatio   = IBP_VolRatio(pid)
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "dIBM_clc_InterpolateCoe_inc.f90"

      !===== interpolation =====!
#include "dIBM_clc_Point_Interpolation_inc.f90"

      !=====    forcing    =====!
      IbpForce= (IBP_Vel(pid) -real3(SumXDir,SumYDir,SumZDir))*ForceAmplifyCoe
      IBP_Force(pid)= IBP_Force(pid)+IbpForce

      !=====   spreading   =====!
#include "dIBM_spreadIBM_force_inc.f90"
    ENDDO
    call Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)
    IF(IsUxConst) call NormVolForce_x(VolForce_x)
    
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
  end subroutine AdditionalForceIBM_0

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! AdditionalForceIBM_1
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine AdditionalForceIBM_1(ux,uy,uz,VolForce_x,VolForce_y,VolForce_z); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,VolForce_x,VolForce_y,VolForce_z

    ! locals
    type(real3):: LagrangeP,IbpForce
    real(RK),dimension(0:nDistribute)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir,VolRatio
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    VolForce_x=0.0_RK; VolForce_y=0.0_RK; VolForce_z=0.0_RK
    DO pid=1,nIBP
      VolRatio   = IBP_VolRatio(pid)
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "dIBM_clc_InterpolateCoe_inc.f90"

      !===== interpolation =====!
#include "dIBM_clc_Point_Interpolation_inc.f90"

      !=====    forcing    =====!
      IbpForce= (IBP_Vel(pid) -real3(SumXDir,SumYDir,SumZDir))*IBP_ForceAmplify(pid)
      IBP_Force(pid)= IBP_Force(pid)+IbpForce

      !=====   spreading   =====!
#include "dIBM_spreadIBM_force_inc.f90"
    ENDDO
    call Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)
    IF(IsUxConst) call NormVolForce_x(VolForce_x)
    
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
  end subroutine AdditionalForceIBM_1

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Reallocate_IbpVar
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Reallocate_IbpVar(); implicit none

    ! locals
    integer:: sizep,sizen,ierrTemp, ierr
    real(RK),dimension(:),allocatable::RealVec
    type(real3),dimension(:),allocatable:: Real3Vec
    integer(kind=2),dimension(:,:),allocatable::IntArr   
 
    ierr=0
    sizep= mIBP
    sizen= int(1.1_RK*real(sizep,kind=RK))
    sizen= max(sizen, nIBP+1)
    mIBP = sizen    

    ! ======= real vector part =======
    call move_alloc(IBP_VolRatio,RealVec)
    allocate(IBP_VolRatio(sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)
    IBP_VolRatio(1:sizep)=RealVec
    deallocate(RealVec)

    ! ======= integer matrix part =======
    call move_alloc(IBP_indxyz,IntArr)
    allocate(IBP_indxyz(6,sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)
    IBP_indxyz(1:6,1:sizep)=IntArr
    deallocate(IntArr)

    ! ======= real3 vercor part =======
    call move_alloc(IBP_Pos,Real3Vec)
    allocate(IBP_Pos(sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)
    IBP_Pos(1:sizep)=Real3Vec

    call move_alloc(IBP_Vel,Real3Vec)
    allocate(IBP_Vel(sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)
    IBP_Vel(1:sizep)=Real3Vec
    deallocate(Real3Vec)

    deallocate(IBP_Force)
    allocate(IBP_Force(sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)

    if(ierr/=0) then
      call MainLog%CheckForError(ErrT_Abort," Reallocate_IbpVar"," Reallocate wrong!")
      call MainLog%OutInfo("The present processor  is :"//strip(num2str(nrank)),3)
    endif    
    !call MainLog%CheckForError(ErrT_Pass,"Reallocate_IbpVar","Need to reallocate IBP variables")
    !call MainLog%OutInfo("The present processor  is :"//strip(num2str(nrank)),3)
    !call MainLog%OutInfo("Previous matirx length is :"//strip(num2str(sizep)),3)
    !call MainLog%OutInfo("Updated  matirx length is :"//strip(num2str(sizen)),3)
  end subroutine Reallocate_IbpVar
#undef nDistribute

end module dIBM_IBM
