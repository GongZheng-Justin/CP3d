#include "definitions_inc.f90"
module ap_ContactSearchPW
  use MPI
  use mc_TypeDef
  use mc_Decomp2d,only: nrank
#ifdef CFDSecondOrder
  use f2_Parameters,only:yly,xlx,zlz
#else
  use f4_Parameters,only:yly,xlx,zlz
#endif
  use ap_Property
  use ap_Variables
  use ap_Parameters
  use ap_Decomp_2d,only:int_type,real_type
  implicit none
  private

  type::ContactSearchPW
  contains
    procedure,nopass:: FindContactsPW
  end type ContactSearchPW
  type(ContactSearchPW),public::ATPContactSearchPW
    
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Performing contact search to determine particle-wall contacts 
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine FindContactsPW(); implicit none

    ! locals
    integer::pid
    real(RK)::PosX,PosY,PosZ

    if(.not. ATP_opt%IsPeriodic(1)) then
      DO pid=1,GPrtcl_list%nlocal
        PosX = GPrtcl_PosR(pid)%x
        if(PosX<=0.0_RK .or. PosX>=xlx) then
          GPrtcl_MoveDistance(pid)%x=GPrtcl_MoveDistance(pid)%x+GPrtcl_PosOld(pid)%x-PosX
          GPrtcl_PosR(pid)%x    =  GPrtcl_PosOld(pid)%x
          GPrtcl_linVel(1,pid)%x= -GPrtcl_linVel(2,pid)%x
          GPrtcl_SwimDir(pid)%x = -GPrtcl_SwimDir(pid)%x
        endif
      ENDDO
    endif
    
    if(.not. ATP_opt%IsPeriodic(2)) then
      DO pid=1,GPrtcl_list%nlocal
        PosY=GPrtcl_PosR(pid)%y
        if(PosY<=0.0_RK .or. PosY>=yly) then
          GPrtcl_MoveDistance(pid)%y=GPrtcl_MoveDistance(pid)%y+GPrtcl_PosOld(pid)%y-PosY
          GPrtcl_PosR(pid)%y    =  GPrtcl_PosOld(pid)%y
          GPrtcl_linVel(1,pid)%y= -GPrtcl_linVel(2,pid)%y
          GPrtcl_SwimDir(pid)%y = -GPrtcl_SwimDir(pid)%y
        endif
      ENDDO    
    endif
    
    if(.not. ATP_opt%IsPeriodic(3)) then
      DO pid=1,GPrtcl_list%nlocal        
        PosZ = GPrtcl_PosR(pid)%z
        if(PosZ<=0.0_RK .or. PosZ>=zlz) then
          GPrtcl_MoveDistance(pid)%z=GPrtcl_MoveDistance(pid)%z+GPrtcl_PosOld(pid)%z-PosZ
          GPrtcl_PosR(pid)%z    =  GPrtcl_PosOld(pid)%z
          GPrtcl_linVel(1,pid)%z= -GPrtcl_linVel(2,pid)%z
          GPrtcl_SwimDir(pid)%z = -GPrtcl_SwimDir(pid)%z
        endif
      ENDDO
    endif
  end subroutine FindContactsPW
    
end module ap_ContactSearchPW
