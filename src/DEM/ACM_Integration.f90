module Prtcl_Integration
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_Parameters
  use m_TScheme,only:IsUxConst
  use m_Parameters,only:gravity,PrGradNow
  implicit none
  private

  real(RK),parameter,dimension(2):: AB2C = (/ three/two,     -one/two /)
  real(RK),parameter,dimension(3):: AB3C = (/ twentythree/twelve, -sixteen/twelve, five/twelve /)
  public::Prtcl_Integrate
#ifdef ObliqueWallTest
  logical::IsRotate=.false.
#endif
contains

  !******************************************************************
  ! Prtcl_Integrate
  !******************************************************************
  subroutine Prtcl_Integrate()
    implicit none
    integer::pid,itype
    real(RK)::dt,dth,MassEff,InertiaEff,MassInFluid
    type(real3)::linVel1,rotVel1,FpLinAcc,FpRotAcc,GravityToTal
#ifdef ObliqueWallTest
    real(RK)::rMagnitude,rxDir,ryDir
#endif
    
    GravityToTal=DEM_Opt%Gravity
    if(IsAddFluidPressureGradient) then
      if(IsUxConst) then
        GravityToTal%x=GravityToTal%x + PrGradNow
      else
        GravityToTal%x=GravityToTal%x + gravity(1)
      endif
    endif

    dt=DEM_opt%dt;   dth=dt*half
    DO pid =1,GPrtcl_list%nlocal
      itype  = GPrtcl_pType(pid)
      MassEff    = PrtclIBMProp(itype)%MassEff
      InertiaEff = PrtclIBMProp(itype)%InertiaEff
      MassInFluid= PrtclIBMProp(itype)%MassinFluid
    
      ! linear velocity and position
      linVel1= GPrtcl_linVel(1,pid)
      GPrtcl_linAcc(1,pid) = (GPrtcl_cntctForce(pid)/MassEff)
      if(GPrtcl_HighSt(pid)=="Y") then   ! Turn off fluid forces for large St collsions
        FpLinAcc= (MassInFluid*GravityToTal)/MassEff
      else
        FpLinAcc= (GPrtcl_FpForce(pid)+ MassInFluid*GravityToTal)/MassEff
      endif
      GPrtcl_linVel(1,pid)=linVel1 + dt*(FpLinAcc+ AB2C(1)*GPrtcl_linAcc(1,pid)+ AB2C(2)*GPrtcl_linAcc(2,pid))
      GPrtcl_linVel(2,pid)=linVel1
      GPrtcl_linAcc(2,pid)=GPrtcl_linAcc(1,pid)

      ! rotate position
      rotVel1= GPrtcl_rotVel(1,pid)
      GPrtcl_rotAcc(1,pid) = (GPrtcl_torque(pid)/InertiaEff)
      FpRotAcc = (GPrtcl_FpTorque(pid)/InertiaEff)
      GPrtcl_rotVel(1,pid)=rotVel1 + dt*(FpRotAcc+ AB2C(1)*GPrtcl_rotAcc(1,pid)+ AB2C(2)*GPrtcl_rotAcc(2,pid))
      GPrtcl_rotVel(2,pid)=rotVel1
      GPrtcl_rotAcc(2,pid)=GPrtcl_rotAcc(1,pid)

#ifdef RotateOnly
      GPrtcl_linVel(1,pid)=zero_r3
#endif
#ifdef ObliqueWallTest
      GPrtcl_linVel(1,pid)%z=zero
      if(.not.IsRotate) then
        GPrtcl_rotVel(1,pid)=zero_r3
        rMagnitude=sqrt(DEM_Opt%Gravity%x*DEM_Opt%Gravity%x+DEM_Opt%Gravity%y*DEM_Opt%Gravity%y)
        rxDir=DEM_Opt%Gravity%x/rMagnitude
        ryDir=DEM_Opt%Gravity%y/rMagnitude
        rMagnitude=sqrt(GPrtcl_linVel(1,pid)%x*GPrtcl_linVel(1,pid)%x+GPrtcl_linVel(1,pid)%y*GPrtcl_linVel(1,pid)%y)
        GPrtcl_linVel(1,pid)%x=rxDir*rMagnitude
        GPrtcl_linVel(1,pid)%y=ryDir*rMagnitude
      endif
#endif
      GPrtcl_PosR(pid)=GPrtcl_PosR(pid)+dth*(linVel1+ GPrtcl_linVel(1,pid))
#ifdef ObliqueWallTest    
      if(GPrtcl_PosR(pid)%y<two*GPrtcl_PosR(pid)%w) IsRotate=.true.
#endif
    ENDDO

  end subroutine Prtcl_Integrate

end module Prtcl_Integration
