module LPT_ContactSearchPW
  use MPI
  use LPT_TypeDef
  use LPT_logInfo
  use LPT_Property
  use LPT_Geometry
  use LPT_Variables
  use LPT_Parameters
  use m_Decomp2d,only: nrank
  use m_Parameters,only:yly,xlx,zlz  
  use LPT_Decomp_2d,only:int_type,real_type
  implicit none
  private

  type::ContactSearchPW
  contains
    procedure:: FindContactsPW
  end type ContactSearchPW
  type(ContactSearchPW),public::LPTContactSearchPW
    
contains

  !**********************************************************************
  ! Performing contact search to determine particle-wall contacts 
  !**********************************************************************
  subroutine FindContactsPW(this)
    implicit none
    class(ContactSearchPW):: this

    ! locals
    integer:: pid
    real(RK)::Posy,radius

#ifdef UseDEMWallContact
    integer::wid
    DO pid=1,GPrtcl_list%nlocal
      do wid=1,LPTGeometry%nPW_local
         if(LPTGeometry%pWall(wid)%isInContact(GPrtcl_PosR(pid),ovrlp,nv))then
           GPrtcl_PosR(pid)%y=  GPrtcl_PosOld(pid)%y
           GPrtcl_linVel(1,pid)%y= - GPrtcl_linVel(2,pid)%y               
         endif 
      enddo
    ENDDO
#else
    DO pid=1,GPrtcl_list%nlocal
      Posy=GPrtcl_PosR(pid)%y
      radius= GPrtcl_PosR(pid)%w
      if(Posy<=radius .or. Posy+radius>=yly) then
        GPrtcl_PosR(pid)%y    =   GPrtcl_PosOld(pid)%y
        GPrtcl_linVel(1,pid)%y= - GPrtcl_linVel(2,pid)%y
      endif
    ENDDO
#endif
  end subroutine FindContactsPW
    
end module LPT_ContactSearchPW
