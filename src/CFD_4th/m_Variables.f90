module m_Variables
  use m_TypeDef
  use m_LogInfo
  use m_decomp2d
  use m_Parameters
  implicit none
  private
  
  ! define all major arrays here 
  real(RK),public,allocatable, dimension(:,:,:) :: ux
  real(RK),public,allocatable, dimension(:,:,:) :: uy
  real(RK),public,allocatable, dimension(:,:,:) :: uz
  real(RK),public,allocatable, dimension(:,:,:) :: HistxOld
  real(RK),public,allocatable, dimension(:,:,:) :: HistyOld
  real(RK),public,allocatable, dimension(:,:,:) :: HistzOld
  real(RK),public,allocatable, dimension(:,:,:) :: pressure
  
  real(RK),public,allocatable, dimension(:,:,:) :: RealArr1
  real(RK),public,allocatable, dimension(:,:,:) :: RealArr2
  real(RK),public,allocatable, dimension(:,:,:) :: Realhalo

  real(RK),public,allocatable, dimension(:,:) :: uy_ym, duy_ym

  type(MatBound),public:: mb1        ! matrix bound type 1

#ifdef CFDLPT_TwoWay
  real(RK),public,allocatable,dimension(:,:,:)::FpForce_x,FpForce_y,FpForce_z
#endif
#ifdef ScalarFlow
  real(RK),public,allocatable,dimension(:,:,:):: scalar
  real(RK),public,allocatable,dimension(:,:,:):: HistCOld
  real(RK),public,allocatable,dimension(:,:,:):: RealArrC
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
    mb1%xme=2;  mb1%xpe=2
    mb1%yme=1;  mb1%ype=1
    mb1%zme=2;  mb1%zpe=2
       
    !-------------------------------------------------
    ! Arrays with ghost cells
    !-------------------------------------------------
    call myallocate(ux, mb1, opt_global=.true.)
    call myallocate(uy, mb1, opt_global=.true.)
    call myallocate(uz, mb1, opt_global=.true.)
    call myallocate(pressure, mb1, opt_global=.true.)
    call myallocate(RealHalo, mb1, opt_global=.true.)
#ifdef CFDLPT_TwoWay
    call myallocate(FpForce_x, mb1, opt_global=.true.)
    call myallocate(FpForce_y, mb1, opt_global=.true.)
    call myallocate(FpForce_z, mb1, opt_global=.true.)
    FpForce_x=zero; FpForce_y=zero; FpForce_z=zero
#endif
#ifdef ScalarFlow
    call myallocate(scalar,   mb1, opt_global=.true.)
    call myallocate(RealArrC, mb1, opt_global=.true.)
    allocate(HistCOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)))
    scalar=Scalar_InitValue; RealArrC=zero; HistCOld=zero
#endif
    !-------------------------------------------------
    ! Arrays without ghost cells
    !-------------------------------------------------
    allocate(HistxOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(HistyOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(HistzOld(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(RealArr1(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(RealArr2(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(uy_ym(y1start(1):y1end(1),  y1start(3):y1end(3)), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(duy_ym(y1start(1):y1end(1), y1start(3):y1end(3)), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"AllocateVariables","Allocation failed")    
    
    ux=zero;         uy=zero;        uz=zero
    HistxOld=zero;   HistyOld=zero;  HistzOld=zero
    pressure=zero;   RealArr1=zero;  RealArr2=zero;  RealHalo=zero;
    uy_ym=zero;      duy_ym=zero
      
  end subroutine AllocateVariables
    
end module m_Variables
