module ATP_ContactSearchPW
  use MPI
  use m_TypeDef
  use ATP_Property
  use ATP_Variables
  use ATP_Parameters
  use m_Decomp2d,only: nrank
  use m_Parameters,only:yly,xlx,zlz  
  use ATP_Decomp_2d,only:int_type,real_type
  implicit none
  private

  type::ContactSearchPW
  contains
    procedure:: FindContactsPW
  end type ContactSearchPW
  type(ContactSearchPW),public::ATPContactSearchPW
    
contains

  !**********************************************************************
  ! Performing contact search to determine particle-wall contacts 
  !**********************************************************************
  subroutine FindContactsPW(this)
    implicit none
    class(ContactSearchPW):: this

    ! locals
    integer::pid
    real(RK)::PosY

#ifdef UseDEMWallContact
    integer::wid
    DO pid=1,GPrtcl_list%nlocal
      do wid=1,ATPGeometry%nPW_local
         if(ATPGeometry%pWall(wid)%isInContact(GPrtcl_PosR(pid),ovrlp,nv))then
           GPrtcl_MoveDistance(pid)%y=GPrtcl_MoveDistance(pid)%y+GPrtcl_PosOld(pid)%y-GPrtcl_PosR(pid)%y
           GPrtcl_PosR(pid)%y=  GPrtcl_PosOld(pid)%y
           GPrtcl_linVel(1,pid)%y= -GPrtcl_linVel(2,pid)%y
           GPrtcl_SwimDir(pid)%y = -GPrtcl_SwimDir(pid)%y           
         endif 
      enddo
    ENDDO
#else
    DO pid=1,GPrtcl_list%nlocal
      PosY=GPrtcl_PosR(pid)%y
      if(PosY<=0.0_RK .or. PosY>=yly) then
        GPrtcl_MoveDistance(pid)%y=GPrtcl_MoveDistance(pid)%y+GPrtcl_PosOld(pid)%y-PosY
        GPrtcl_PosR(pid)%y    =  GPrtcl_PosOld(pid)%y
        GPrtcl_linVel(1,pid)%y= -GPrtcl_linVel(2,pid)%y
        GPrtcl_SwimDir(pid)%y = -GPrtcl_SwimDir(pid)%y
      endif
    ENDDO
#endif
  end subroutine FindContactsPW
    
end module ATP_ContactSearchPW
