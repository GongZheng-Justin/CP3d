module LPT_Integration
  use LPT_TypeDef
  use LPT_LogInfo
  use LPT_Property
  use LPT_Variables
  use LPT_Parameters
  implicit none
  private
    
  real(RK),parameter,dimension(2):: AB2C = (/ three/two,     -one/two /)
  real(RK),parameter,dimension(3):: AB3C = (/ twentythree/twelve, -sixteen/twelve, five/twelve /)
    
  public::Prtcl_Integrate
contains
  
  subroutine Prtcl_Integrate()
    implicit none
    real(RK)::dt,Mass
    integer::pid,nlocal,itype
    type(real3):: linVel1,linVel2,PosOld,Gravity
    
    dt=LPT_opt%dt
    nlocal = GPrtcl_list%nlocal
    Gravity= LPT_opt%gravity
    
    ! linear position
    if(LPT_Opt%PI_Method==PIM_AB2) then
      DO pid = 1,nlocal 
        itype= GPrtcl_pType(pid)
        Mass = LPTProperty%Prtcl_PureProp(itype)%Mass      
        GPrtcl_linAcc(1,pid)= GPrtcl_FpForce(pid)/Mass + Gravity
        
        PosOld = GPrtcl_PosR(pid)
        GPrtcl_PosOld(pid) = PosOld

        linVel1=GPrtcl_linVel(1,pid)
        GPrtcl_PosR(pid)=PosOld+(AB2C(1)*linVel1 + AB2C(2)*GPrtcl_linVel(2,pid))*dt
        GPrtcl_linVel(1,pid)=linVel1+(AB2C(1)*GPrtcl_linAcc(1,pid)+AB2C(2)*GPrtcl_linAcc(2,pid))*dt
        GPrtcl_linVel(2,pid)=linVel1
        GPrtcl_linAcc(2,pid)=GPrtcl_linAcc(1,pid)
      ENDDO

    elseif(LPT_Opt%PI_Method==PIM_AB3 ) then
      DO pid=1,nlocal
        itype= GPrtcl_pType(pid)
        Mass = LPTProperty%Prtcl_PureProp(itype)%Mass      
        GPrtcl_linAcc(1,pid)= GPrtcl_FpForce(pid)/Mass + Gravity
        
        PosOld = GPrtcl_PosR(pid)
        GPrtcl_PosOld(pid) = PosOld

        linVel1=GPrtcl_linVel(1,pid)
        linVel2=GPrtcl_linVel(2,pid)
                
        GPrtcl_PosR(pid)=PosOld+(AB3C(1)*linVel1+AB3C(2)*linVel2+AB3C(3)*GPrtcl_linVel(3,pid))*dt
        GPrtcl_linVel(1,pid) =linVel1+(AB3C(1)*GPrtcl_linAcc(1,pid)+AB3C(2)*GPrtcl_linAcc(2,pid)+ &
                                     AB3C(3)*GPrtcl_linAcc(3,pid))*dt
        GPrtcl_linVel(3,pid) = linVel2
        GPrtcl_linVel(2,pid) = linVel1
        GPrtcl_linAcc(3,pid) = GPrtcl_linAcc(2,pid)
        GPrtcl_linAcc(2,pid) = GPrtcl_linAcc(1,pid)
      ENDDO
    endif
    
  end subroutine Prtcl_Integrate

end module LPT_Integration
