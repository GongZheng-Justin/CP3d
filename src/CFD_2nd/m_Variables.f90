module m_Variables
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  implicit none
  private
  
  ! define all major arrays here 
  real(RK), public,allocatable, dimension(:,:,:) :: ux
  real(RK), public,allocatable, dimension(:,:,:) :: uy
  real(RK), public,allocatable, dimension(:,:,:) :: uz
  real(RK), public,allocatable, dimension(:,:,:) :: HistxOld
  real(RK), public,allocatable, dimension(:,:,:) :: HistyOld
  real(RK), public,allocatable, dimension(:,:,:) :: HistzOld
  real(RK), public,allocatable, dimension(:,:,:) :: pressure
  
  real(RK), public,allocatable, dimension(:,:,:) :: RealArr1
  real(RK), public,allocatable, dimension(:,:,:) :: RealArr2
  real(RK), public,allocatable, dimension(:,:,:) :: Realhalo

  real(RK), public,allocatable, dimension(:,:,:) :: nut
  real(RK), public,allocatable, dimension(:,:) :: uy_ym, duy_ym
  real(RK), public,allocatable, dimension(:,:,:)::OutFlowInfoX,OutFlowInfoY

  type(MatBound),public:: mb1  ! matrix bound type 1
  type(HaloInfo),public:: hi1  ! halo info type 1

#ifdef CFDDEM
  type(MatBound),public:: mb_dist  ! matrix bound for distribution
  type(HaloInfo),public:: hi_dist  ! halo info for distribution
  real(RK), public,allocatable, dimension(:,:,:)::FpForce_x,FpForce_y,FpForce_z
#endif

#ifdef CFDACM
  type(MatBound),public:: mb_dist  ! matrix bound for distribution
  type(HaloInfo),public:: hi_dist  ! halo info for distribution
  real(RK), public,allocatable, dimension(:,:,:):: IBMArr1,IBMArr2,IBMArr3
  character,public,allocatable, dimension(:,:,:):: FluidIndicator
#endif

#ifdef CFDLPT_TwoWay
  real(RK),public,allocatable, dimension(:,:,:)::FpForce_x,FpForce_y,FpForce_z
#endif
    
  public:: AllocateVariables
contains  

  !******************************************************************
  ! InverseTridiagonal
  !****************************************************************** 
  subroutine AllocateVariables()
    implicit none
      
    ! locals
    integer::ierrTmp,ierror=0

    mb1%pencil = y_pencil  
    mb1%xme=1;  mb1%xpe=2
    mb1%yme=1;  mb1%ype=2
    mb1%zme=1;  mb1%zpe=2
    
    !-------------------------------------------------
    ! Arrays with ghost cells
    !-------------------------------------------------
    call myallocate(ux, mb1, opt_global=.true.)
    call myallocate(uy, mb1, opt_global=.true.)
    call myallocate(uz, mb1, opt_global=.true.)
    call myallocate(pressure, mb1, opt_global=.true.)
    call myallocate(RealHalo, mb1, opt_global=.true.)
    if(LES_type>0) then 
      call myallocate(nut, mb1, opt_global=.true.)
      nut=zero
    endif
#ifdef CFDACM
    call myallocate(IBMArr1, mb1, opt_global=.true.)
    call myallocate(IBMArr2, mb1, opt_global=.true.)
    call myallocate(IBMArr3, mb1, opt_global=.true.)
    IBMArr1=zero; IBMArr2=zero; IBMArr3=zero
#endif
#ifdef CFDLPT_TwoWay
    call myallocate(FpForce_x, mb1, opt_global=.true.)
    call myallocate(FpForce_y, mb1, opt_global=.true.)
    call myallocate(FpForce_z, mb1, opt_global=.true.)
    FpForce_x=zero; FpForce_y=zero; FpForce_z=zero
#endif
   
    !-------------------------------------------------
    ! Arrays without ghost cells
    !-------------------------------------------------
    allocate(HistxOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(HistyOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(HistzOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(RealArr1(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(RealArr2(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
#ifdef CFDACM
    allocate(FluidIndicator(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
#endif
    allocate(uy_ym(y1start(1):y1end(1),  y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(duy_ym(y1start(1):y1end(1), y1start(3):y1end(3)),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0)call MainLog%CheckForError(ErrT_Abort,"AllocateVariables","Allocation failed 1")
    ux=zero;         uy=zero;        uz=zero
    HistxOld=zero;   HistyOld=zero;  HistzOld=zero
    pressure=zero;   RealArr1=zero;  RealArr2=zero;  RealHalo=zero;
    uy_ym=zero;      duy_ym=zero
     
    ! xp - outflow
    if(myProcNghBC(y_pencil,3)==BC_OutFlow) then
      allocate(OutFlowInfoX(6,y1start(2):y1end(2),y1start(3):y1end(3)),stat=ierror); OutFlowInfoX=zero
    endif
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"AllocateVariables","Allocation failed 2")
    
    ! yp - outflow
    if(BcOption(yp_dir)==BC_OutFlow) then
      allocate(OutFlowInfoY(6,y1start(1):y1end(1),y1start(3):y1end(3)),stat=ierror); OutFlowInfoY=zero    
    endif
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"AllocateVariables","Allocation failed 3")
  end subroutine AllocateVariables
    
end module m_Variables
