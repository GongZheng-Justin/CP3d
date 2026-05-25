#include "definitions_inc.f90"
module f2_FlowType_AddedNew
  use MPI
  use mc_TypeDef
  use mc_Decomp2d
  use f2_Parameters
  use f2_MeshAndMetries
  use f2_Variables,only: mb1
  use f2_Tools,only:CalcUxAver
  implicit none
  private    

  public:: InitVelocity_AN, InitStatVar_AN,  clcStat_AN
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitVelocity_AN
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitVelocity_AN(ux,uy,uz); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
  
    ! locals
    integer :: ic,jc,kc
    real(RK):: VelRef,Ratiot,yct
    
    VelRef= uxBcValue(ym_dir)
    Ratiot=(uxBcValue(yp_dir)-uxBcValue(ym_dir))/yly
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        yct=yc(jc)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)=  VelRef+Ratiot*yct
        enddo
      enddo
    enddo
    uy=0.0_RK
    uz=0.0_RK
  end subroutine InitVelocity_AN

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitStatVar_AN
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitStatVar_AN(); implicit none


  end subroutine InitStatVar_AN

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clcStat_CH
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine clcStat_AN(); implicit none
   
  end subroutine clcStat_AN
    
end module f2_FlowType_AddedNew
