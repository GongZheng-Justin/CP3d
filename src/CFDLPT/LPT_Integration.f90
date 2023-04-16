module LPT_Integration
  use m_TypeDef
  use LPT_Property
  use LPT_Variables
  use LPT_Parameters
  implicit none
  private    
  real(RK),parameter,dimension(2):: AB2C = [1.5_RK,-0.5_RK]
  real(RK),parameter,dimension(3):: AB3C = [23.0_RK,-16.0_RK,5.0_RK]/12.0_RK
    
  public::Prtcl_Integrate
contains

  !******************************************************************
  ! Prtcl_Integrate
  !******************************************************************  
  subroutine Prtcl_Integrate(iCountLPT)
    implicit none
    integer,intent(in)::iCountLPT
     
    ! locals
    integer::pid,nlocal,itype
    real(RK)::dt,Mass,TimeIntCoe(3)
    type(real3)::linVel1,linVel2,PosOld,Gravity
    
    dt=LPT_opt%dt
    nlocal = GPrtcl_list%nlocal
    Gravity= LPT_opt%gravity
        
    ! linear position
    if(LPT_Opt%PI_Method==PIM_AB2) then
      if(iCountLPT==1) then
        TimeIntCoe(1)=1.0_RK
        TimeIntCoe(2)=0.0_RK
      else
        TimeIntCoe(1)=AB2C(1)
        TimeIntCoe(2)=AB2C(2)
      endif
      DO pid = 1,nlocal 
        itype= GPrtcl_pType(pid)
        Mass = LPTProperty%Prtcl_PureProp(itype)%Mass      
        GPrtcl_linAcc(1,pid)= GPrtcl_FpForce(pid)/Mass + Gravity
        
        PosOld = GPrtcl_PosR(pid)
        GPrtcl_PosOld(pid) = PosOld

        linVel1=GPrtcl_linVel(1,pid)
        GPrtcl_PosR(pid)=PosOld+(TimeIntCoe(1)*linVel1 + TimeIntCoe(2)*GPrtcl_linVel(2,pid))*dt
        GPrtcl_linVel(1,pid)=linVel1+(TimeIntCoe(1)*GPrtcl_linAcc(1,pid)+TimeIntCoe(2)*GPrtcl_linAcc(2,pid))*dt
        GPrtcl_linVel(2,pid)=linVel1
        GPrtcl_linAcc(2,pid)=GPrtcl_linAcc(1,pid)
      ENDDO

    elseif(LPT_Opt%PI_Method==PIM_AB3 ) then
      if(iCountLPT==1) then
        TimeIntCoe(1)=1.0_RK
        TimeIntCoe(2)=0.0_RK
        TimeIntCoe(3)=0.0_RK
      elseif(iCountLPT==2) then
        TimeIntCoe(1)=AB2C(1)
        TimeIntCoe(2)=AB2C(2)
        TimeIntCoe(3)=0.0_RK 
      else
        TimeIntCoe(1)=AB3C(1)
        TimeIntCoe(2)=AB3C(2)
        TimeIntCoe(3)=AB3C(3)    
      endif
      DO pid=1,nlocal
        itype= GPrtcl_pType(pid)
        Mass = LPTProperty%Prtcl_PureProp(itype)%Mass      
        GPrtcl_linAcc(1,pid)= GPrtcl_FpForce(pid)/Mass + Gravity
        
        PosOld = GPrtcl_PosR(pid)
        GPrtcl_PosOld(pid) = PosOld

        linVel1=GPrtcl_linVel(1,pid)
        linVel2=GPrtcl_linVel(2,pid)
                
        GPrtcl_PosR(pid)=PosOld+(TimeIntCoe(1)*linVel1+TimeIntCoe(2)*linVel2+TimeIntCoe(3)*GPrtcl_linVel(3,pid))*dt
        GPrtcl_linVel(1,pid) =linVel1+(TimeIntCoe(1)*GPrtcl_linAcc(1,pid)+TimeIntCoe(2)*GPrtcl_linAcc(2,pid)+ &
                                     TimeIntCoe(3)*GPrtcl_linAcc(3,pid))*dt
        GPrtcl_linVel(3,pid) = linVel2
        GPrtcl_linVel(2,pid) = linVel1
        GPrtcl_linAcc(3,pid) = GPrtcl_linAcc(2,pid)
        GPrtcl_linAcc(2,pid) = GPrtcl_linAcc(1,pid)
      ENDDO
    endif
    
  end subroutine Prtcl_Integrate

end module LPT_Integration
