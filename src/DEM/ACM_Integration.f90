module Prtcl_Integration
  use m_TypeDef
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_Parameters
  use m_Parameters,only:gravity,PrGradData,IsUxConst
  implicit none
  private
  real(RK),parameter,dimension(2):: AB2C = [1.5_RK,-0.5_RK]
  real(RK),parameter,dimension(3):: AB3C = [23.0_RK,-16.0_RK,5.0_RK]/12.0_RK
  
  public::Prtcl_Integrate
#ifdef ObliqueWallTest
  logical::IsRotate=.false.
#endif
contains

  !******************************************************************
  ! Prtcl_Integrate
  !******************************************************************
  subroutine Prtcl_Integrate(iCountACM)
    implicit none
    integer,intent(in)::iCountACM
    
    ! locals
    integer::pid,itype
    type(real3)::linVel1,rotVel1,FpLinAcc,FpRotAcc,GravityToTal
    real(RK)::dt,dth,MassEff,InertiaEff,MassInFluid,TimeIntCoe(2)
#ifdef ObliqueWallTest
    real(RK)::rMagnitude,rxDir,ryDir
#endif
    
    if(iCountACM==1) then
      TimeIntCoe(1)=1.0_RK
      TimeIntCoe(2)=0.0_RK
    else
      TimeIntCoe(1)=AB2C(1)
      TimeIntCoe(2)=AB2C(2)
    endif
    GravityToTal=DEM_Opt%Gravity
    if(IsAddFluidPressureGradient .and. IsUxConst) GravityToTal%x= PrGradData(2)
    
    dt=DEM_opt%dt;   dth=dt*0.5_RK
    DO pid =1,GPrtcl_list%nlocal
      itype  = GPrtcl_pType(pid)
      InertiaEff = PrtclIBMProp(itype)%InertiaEff
      MassInFluid= PrtclIBMProp(itype)%MassinFluid
    
      ! linear velocity and position
      linVel1= GPrtcl_linVel(1,pid)
      if(GPrtcl_HighSt(pid)=="Y") then   ! Turn off fluid forces for large St collisions
        MassEff = DEMProperty%Prtcl_PureProp(itype)%Mass
        FpLinAcc= (MassInFluid*GravityToTal)/MassEff
      else
        MassEff = PrtclIBMProp(itype)%MassEff
        FpLinAcc= (GPrtcl_FpForce(pid)+ MassInFluid*GravityToTal)/MassEff
      endif
      GPrtcl_linAcc(1,pid) = (GPrtcl_cntctForce(pid)/MassEff)
      GPrtcl_linVel(1,pid)=linVel1 + dt*(FpLinAcc+ TimeIntCoe(1)*GPrtcl_linAcc(1,pid)+ TimeIntCoe(2)*GPrtcl_linAcc(2,pid))
      GPrtcl_linVel(2,pid)=linVel1
      GPrtcl_linAcc(2,pid)=GPrtcl_linAcc(1,pid)

      ! rotate position
      rotVel1= GPrtcl_rotVel(1,pid)
      GPrtcl_rotAcc(1,pid) = (GPrtcl_torque(pid)/InertiaEff)
      FpRotAcc = (GPrtcl_FpTorque(pid)/InertiaEff)
      GPrtcl_rotVel(1,pid)=rotVel1 + dt*(FpRotAcc+ TimeIntCoe(1)*GPrtcl_rotAcc(1,pid)+ TimeIntCoe(2)*GPrtcl_rotAcc(2,pid))
      GPrtcl_rotVel(2,pid)=rotVel1
      GPrtcl_rotAcc(2,pid)=GPrtcl_rotAcc(1,pid)

#ifdef RotateOnly
      GPrtcl_linVel(1,pid)=zero_r3
#endif
#ifdef ObliqueWallTest
      GPrtcl_linVel(1,pid)%z=0.0_RK
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
      if(GPrtcl_PosR(pid)%y<2.0_RK*GPrtcl_PosR(pid)%w) IsRotate=.true.
#endif
    ENDDO

  end subroutine Prtcl_Integrate

end module Prtcl_Integration
