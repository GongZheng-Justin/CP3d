module ATP_Integration
  use m_TypeDef
  use ATP_Property
  use ATP_Variables
  use ATP_Parameters
  implicit none
  private
  logical::AddRandInOrientation=.true.
  real(RK)::twopi=6.283185307179586477_RK
  real(RK),parameter,dimension(2)::AB2C=[1.5_RK, -0.5_RK]
  real(RK),parameter,dimension(3)::AB3C=[23.0_RK,-16.0_RK,5.0_RK]/12.0_RK
    
  public::Prtcl_Integrate
contains

  !**********************************************************************
  ! Prtcl_Integrate
  !**********************************************************************  
  subroutine Prtcl_Integrate(iCountATP)
    implicit none
    integer,intent(in)::iCountATP
    
    ! locals
    integer::pid,nlocal,itype
    real(RK)::dt,Theta,Phi,normValue,DiffR,TimeIntCoe(3)
    real(RK),allocatable,dimension(:)::randTmp1,randTmp2
    type(real3)::PosOld,SwimOld,real3T,exPrime,eyPrime,ezPrime,MoveIncreasement
    
    dt=ATP_opt%dt
    nlocal = GPrtcl_list%nlocal
    
    ! linear position
    if(ATP_Opt%PI_Method==PIM_AB2) then
      if(iCountATP==1) then
        TimeIntCoe(1)=1.0_RK
        TimeIntCoe(2)=0.0_RK
      else
        TimeIntCoe(1)=AB2C(1)
        TimeIntCoe(2)=AB2C(2)
      endif
      DO pid=1,nlocal
        PosOld =GPrtcl_PosR(pid)
        GPrtcl_PosOld(pid)=PosOld
        MoveIncreasement=(TimeIntCoe(1)*GPrtcl_linVel(1,pid)+TimeIntCoe(2)*GPrtcl_linVel(2,pid))*dt
        GPrtcl_PosR(pid)=PosOld+MoveIncreasement
        GPrtcl_MoveDistance(pid)=GPrtcl_MoveDistance(pid)+MoveIncreasement
        
        SwimOld = GPrtcl_SwimDir(pid)+(TimeIntCoe(1)*GPrtcl_SwimAcc(1,pid)+TimeIntCoe(2)*GPrtcl_SwimAcc(2,pid))*dt
        normValue=norm(SwimOld)
        GPrtcl_SwimDir(pid)=SwimOld/normValue
      ENDDO

    elseif(ATP_Opt%PI_Method==PIM_AB3) then
      if(iCountATP==1) then
        TimeIntCoe(1)=1.0_RK
        TimeIntCoe(2)=0.0_RK
        TimeIntCoe(3)=0.0_RK
      elseif(iCountATP==2) then
        TimeIntCoe(1)=AB2C(1)
        TimeIntCoe(2)=AB2C(2)
        TimeIntCoe(3)=0.0_RK 
      else
        TimeIntCoe(1)=AB3C(1)
        TimeIntCoe(2)=AB3C(2)
        TimeIntCoe(3)=AB3C(3)    
      endif
      DO pid=1,nlocal
        PosOld = GPrtcl_PosR(pid)
        GPrtcl_PosOld(pid)= PosOld
        MoveIncreasement=(TimeIntCoe(1)*GPrtcl_linVel(1,pid)+TimeIntCoe(2)*GPrtcl_linVel(2,pid)+TimeIntCoe(3)*GPrtcl_linVel(3,pid))*dt
        GPrtcl_PosR(pid)=PosOld+MoveIncreasement
        GPrtcl_MoveDistance(pid)=GPrtcl_MoveDistance(pid)+MoveIncreasement
        
        SwimOld = GPrtcl_SwimDir(pid)+(TimeIntCoe(1)*GPrtcl_SwimAcc(1,pid)+TimeIntCoe(2)*GPrtcl_SwimAcc(2,pid)+TimeIntCoe(3)*GPrtcl_SwimAcc(3,pid))*dt
        normValue=norm(SwimOld)
        GPrtcl_SwimDir(pid)=SwimOld/normValue
      ENDDO
    endif
    if(AddRandInOrientation) then
      if(nlocal<1) return
      allocate(randTmp1(nlocal),randTmp2(nlocal))
      call normrnd(randTmp1)
      call random_number(randTmp2)
      DO pid=1,nlocal
        itype=GPrtcl_pType(pid)
        Phi=twopi*randTmp2(pid)
        DiffR=ATPProperty%Prtcl_PureProp(itype)%DiffuseR
        Theta=sqrt(4.0_RK*dt*DiffR)*randTmp1(pid)
        real3T%x=sin(Theta)*cos(Phi)
        real3T%y=sin(Theta)*sin(Phi)
        real3T%z=cos(Theta)
        ezPrime=GPrtcl_SwimDir(pid)
        call generateCoord(exPrime,eyPrime,ezPrime)
        SwimOld%x= exPrime%x*real3T%x +eyPrime%x*real3T%y +ezPrime%x*real3T%z
        SwimOld%y= exPrime%y*real3T%x +eyPrime%y*real3T%y +ezPrime%y*real3T%z
        SwimOld%z= exPrime%z*real3T%x +eyPrime%z*real3T%y +ezPrime%z*real3T%z
        normValue=norm(SwimOld)
        GPrtcl_SwimDir(pid)=SwimOld/normValue    
      ENDDO
      deallocate(randTmp1,randTmp2)
    endif
  end subroutine Prtcl_Integrate

  !**********************************************************************
  ! generateCoord
  !**********************************************************************
  subroutine generateCoord(exPrime,eyPrime,ezPrime)
    implicit none
    type(real3),intent(in) ::ezPrime
    type(real3),intent(out)::exPrime,eyPrime 
    
    ! locals
    integer::iMin
    real(RK)::MinValue,rTmp
    
    iMin=1;
    MinValue=abs(ezPrime%x)
    rTmp=abs(ezPrime%y)
    if(rTmp<MinValue) then
      iMin=2; MinValue=rTmp
    endif     
    rTmp=abs(ezPrime%z)
    if(rTmp<MinValue) then
      iMin=3; MinValue=rTmp
    endif
    
    if(iMin==1) then
      rTmp=sqrt(1.0_RK-ezPrime%x*ezPrime%x)
      exPrime%x= 0.0_RK
      exPrime%y= ezPrime%z/rTmp
      exPrime%z=-ezPrime%y/rTmp
    endif
    if(iMin==2) then
      rTmp=sqrt(1.0_RK-ezPrime%y*ezPrime%y)
      exPrime%y= 0.0_RK
      exPrime%x= ezPrime%z/rTmp
      exPrime%z=-ezPrime%x/rTmp
    endif
    if(iMin==3) then
      rTmp=sqrt(1.0_RK-ezPrime%z*ezPrime%z)
      exPrime%z= 0.0_RK
      exPrime%x= ezPrime%y/rTmp
      exPrime%y=-ezPrime%x/rTmp
    endif
    eyPrime=(ezPrime .cross. exPrime)
  end subroutine generateCoord
  
  !**********************************************************************
  ! normrnd
  !**********************************************************************
  subroutine normrnd(x)
    implicit none
    real(RK),dimension(:),intent(out)::x

    ! locals
    integer::n,m
    real(RK),dimension(:),allocatable::r
     
    n=size(x,1)
    allocate(r(n+1))
    if(n==1) then
      call random_number(harvest=r(1:2))
      x(1)= sqrt(-2.0_RK*log(r(1)))*cos(twopi*r(2))
    elseif(mod(n,2)==0) then
      m=n/2
      call random_number(r(1:n))  
      x(1:n-1:2)=&
        sqrt(-2.0_RK*log(r(1:2*m-1:2))) &
        *cos(twopi*r(2:2*m:2) )
      x(2:n:2)=&
        sqrt(-2.0_RK*log(r(1:2*m-1:2))) &
        *sin(twopi*r(2:2*m:2) )
    else
      m = (n+1)/2
      call random_number(r(1:2*m))
      x(1:n-2:2) = &
        sqrt ( -2.0_RK * log ( r(1:2*m-3:2) ) ) &
        *cos(twopi*r(2:2*m-2:2) )
      x(2:n-1:2) = &
        sqrt( -2.0_RK * log ( r(1:2*m-3:2) ) ) &
        *sin(twopi*r(2:2*m-2:2) )
      x(n) = sqrt ( -2.0_RK * log ( r(2*m-1) ) ) &
        *cos(twopi*r(2*m))
    endif  
  end subroutine normrnd  
end module ATP_Integration
