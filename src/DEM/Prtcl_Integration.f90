module Prtcl_Integration
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_Parameters
#ifdef CFDDEM
  use m_TScheme,only:IsUxConst
  use m_Parameters,only:gravity,PrGradNow
#endif
  implicit none
  private
 
  real(RK),parameter,dimension(2):: AB2C = (/ three/two,     -one/two /)
  real(RK),parameter,dimension(3):: AB3C = (/ twentythree/twelve, -sixteen/twelve, five/twelve /)
    
  public::Prtcl_Integrate
contains
  
  !**********************************************************************
  ! calculating acceleration of particles (linear and angular) 
  !**********************************************************************
#ifdef CFDDEM
  subroutine clc_Acceleration()
    implicit none

    ! locals
    integer:: i,itype
    type(real3):: Fpforce,GravityToTal
    real(RK)::Mass,Inertia,MassInFluid,MassTot,PrtclDensity

    GravityToTal=DEM_Opt%Gravity
    if(IsAddFluidPressureGradient) then
      if(IsUxConst) then
        GravityToTal%x=GravityToTal%x + PrGradNow
      else
        GravityToTal%x=GravityToTal%x + gravity(1)
      endif
    endif
    
    do i=1,GPrtcl_list%nlocal
      itype   = GPrtcl_pType(i)
      Mass    = DEMProperty%Prtcl_PureProp(itype)%Mass
      
      MassTot    = DEMProperty%Prtcl_PureProp(itype)%MassOfFluid*half +Mass
      MassInFluid= DEMProperty%Prtcl_PureProp(itype)%MassInFluid
      Inertia    = DEMProperty%Prtcl_PureProp(itype)%Inertia
      Fpforce = 1.50_RK*GPrtcl_FpForce(i)-0.50_RK*GPrtcl_FpForce_old(i)
      
      GPrtcl_linAcc(1,i) = ((GPrtcl_cntctForce(i)+Fpforce)/MassTot)+ MassInFluid/MassTot*GravityToTal
      GPrtcl_rotAcc(1,i) = ((GPrtcl_torque(i)/Inertia))
    enddo
  end subroutine clc_Acceleration

#else
  subroutine clc_Acceleration()
    implicit none

    ! locals
    integer:: i,itype
    real(RK)::Mass, Inertia
    
    do i=1,GPrtcl_list%nlocal
      itype   = GPrtcl_pType(i)
      Mass    = DEMProperty%Prtcl_PureProp(itype)%Mass
      Inertia = DEMProperty%Prtcl_PureProp(itype)%Inertia
      GPrtcl_linAcc(1,i) = ((GPrtcl_cntctForce(i))/Mass)+DEM_opt%gravity
      GPrtcl_rotAcc(1,i) = ((GPrtcl_torque(i)/Inertia))
    enddo
  end subroutine clc_Acceleration
#endif

  !**********************************************************************
  ! calculating acceleration of particles (linear and angular) 
  !**********************************************************************
  subroutine Prtcl_Integrate()
    implicit none
    
    ! locals
    real(RK)::dt
    integer::i,nlocal
    type(real3)::linVel1,linVel2,rotVel1,rotVel2
    
    call clc_Acceleration()
    
    dt=DEM_opt%dt
    nlocal = GPrtcl_list%nlocal

    ! linear position
    if(DEM_Opt%PI_Method==PIM_FE) then
      DO i = 1,nlocal 
        GPrtcl_PosR(i) = GPrtcl_PosR(i) + GPrtcl_linVel(1,i) *dt
        GPrtcl_linVel(1,i) = GPrtcl_linVel(1,i) + GPrtcl_linAcc(1,i) *dt
      ENDDO

    elseif(DEM_Opt%PI_Method==PIM_AB2) then
      DO i = 1,nlocal 
        linVel1=GPrtcl_linVel(1,i)
        GPrtcl_PosR(i)=GPrtcl_PosR(i)+(AB2C(1)*linVel1 + AB2C(2)*GPrtcl_linVel(2,i))*dt
        GPrtcl_linVel(1,i)=linVel1+(AB2C(1)*GPrtcl_linAcc(1,i)+AB2C(2)*GPrtcl_linAcc(2,i))*dt
        GPrtcl_linVel(2,i)=linVel1
        GPrtcl_linAcc(2,i)=GPrtcl_linAcc(1,i)
      ENDDO

    elseif(DEM_Opt%PI_Method==PIM_AB3 ) then
      DO i=1,nlocal
        linVel1=GPrtcl_linVel(1,i)
        linVel2=GPrtcl_linVel(2,i)
                
        GPrtcl_PosR(i)=GPrtcl_PosR(i)+(AB3C(1)*linVel1+AB3C(2)*linVel2+AB3C(3)*GPrtcl_linVel(3,i))*dt
        GPrtcl_linVel(1,i) =linVel1+(AB3C(1)*GPrtcl_linAcc(1,i)+AB3C(2)*GPrtcl_linAcc(2,i)+ AB3C(3)*GPrtcl_linAcc(3,i))*dt

        GPrtcl_linVel(3,i) = linVel2
        GPrtcl_linVel(2,i) = linVel1
        GPrtcl_linAcc(3,i) = GPrtcl_linAcc(2,i)
        GPrtcl_linAcc(2,i) = GPrtcl_linAcc(1,i)
      ENDDO
    endif
    
    ! rotate position
    if(DEM_Opt%PRI_Method==PIM_FE) then
      DO i=1,nlocal
        GPrtcl_theta(i)= GPrtcl_theta(i)+ GPrtcl_rotVel(1,i)*dt
        GPrtcl_rotVel(1,i) = GPrtcl_rotVel(1,i) + GPrtcl_rotAcc(1,i) *dt
      ENDDO        

    elseif(DEM_Opt%PRI_Method==PIM_AB2) then
      DO i=1,nlocal 
        rotVel1=GPrtcl_linVel(1,i)
                
        GPrtcl_theta(i)=GPrtcl_theta(i)+(AB2C(1)*rotVel1+AB2C(2)*GPrtcl_rotVel(2,i))*dt
        GPrtcl_rotVel(1,i)=rotVel1+(AB2C(1)*GPrtcl_rotAcc(1,i)+AB2C(2)*GPrtcl_rotAcc(2,i))*dt
                
        GPrtcl_rotVel(2,i)=rotVel1
        GPrtcl_rotAcc(2,i)=GPrtcl_rotAcc(1,i)
      ENDDO        

    elseif(DEM_Opt%PRI_Method==PIM_AB3 ) then
      DO i=1,nlocal 
        rotVel1=GPrtcl_rotVel(1,i)
        rotVel2=GPrtcl_rotVel(2,i)
        GPrtcl_theta(i)=GPrtcl_theta(i)+(AB3C(1)*rotVel1+AB3C(2)*rotVel2+AB3C(3)*GPrtcl_rotVel(3,i))*dt
        GPrtcl_rotVel(1,i)=rotVel1+(AB3C(1)*GPrtcl_rotAcc(1,i)+AB3C(2)*GPrtcl_rotAcc(2,i)+ &
                                    AB3C(3)*GPrtcl_rotAcc(3,i))*dt
                
        GPrtcl_rotVel(3,i) = rotVel2
        GPrtcl_rotVel(2,i) = rotVel1
        GPrtcl_rotAcc(3,i) = GPrtcl_rotAcc(2,i)
        GPrtcl_rotAcc(2,i) = GPrtcl_rotAcc(1,i)
      ENDDO
    endif
    
  end subroutine Prtcl_Integrate

end module Prtcl_Integration
