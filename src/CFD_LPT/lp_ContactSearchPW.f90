#include "definitions_inc.f90"
module lp_ContactSearchPW
  use MPI
  use mc_TypeDef
  use mc_Decomp2d,only: nrank
#ifdef CFDSecondOrder
  use f2_Parameters,only:yly,xlx,zlz  
#else  
  use f4_Parameters,only:yly,xlx,zlz  
#endif    
  use lp_Property
  use lp_Geometry
  use lp_Variables
  use lp_Parameters
  use lp_Decomp_2d,only:int_type,real_type
  implicit none
  private

  type::ContactSearchPW
  contains
    procedure,nopass:: FindContactsPW
  end type ContactSearchPW
  type(ContactSearchPW),public::LPTContactSearchPW
    
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Performing contact search to determine particle-wall contacts 
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine FindContactsPW(); implicit none

    ! locals
    integer:: pid
    real(RK)::Posx, Posy, Posz, radius

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

    if(.not. LPT_opt%IsPeriodic(1)) then
       DO pid=1,GPrtcl_list%nlocal
         Posx=GPrtcl_PosR(pid)%x
         radius= GPrtcl_PosR(pid)%w
         if(Posx<=radius .or. Posx+radius>=xlx) then
           GPrtcl_PosR(pid)%x    =   GPrtcl_PosOld(pid)%x
           GPrtcl_linVel(1,pid)%x= - GPrtcl_linVel(2,pid)%x
         endif
       ENDDO
    endif

    if(.not. LPT_opt%IsPeriodic(2)) then
       DO pid=1,GPrtcl_list%nlocal
         Posy=GPrtcl_PosR(pid)%y
         radius= GPrtcl_PosR(pid)%w
         if(Posy<=radius .or. Posy+radius>=yly) then
           GPrtcl_PosR(pid)%y    =   GPrtcl_PosOld(pid)%y
           GPrtcl_linVel(1,pid)%y= - GPrtcl_linVel(2,pid)%y
         endif
       ENDDO
    endif

    if(.not. LPT_opt%IsPeriodic(3)) then
       DO pid=1,GPrtcl_list%nlocal
         Posz=GPrtcl_PosR(pid)%z
         radius= GPrtcl_PosR(pid)%w
         if(Posz<=radius .or. Posz+radius>=zlz) then
           GPrtcl_PosR(pid)%z    =   GPrtcl_PosOld(pid)%z
           GPrtcl_linVel(1,pid)%z= - GPrtcl_linVel(2,pid)%z
         endif
       ENDDO
    endif

#endif
  end subroutine FindContactsPW
    
end module lp_ContactSearchPW
