#include "definitions_inc.f90"
#ifndef IBMDistributeLinear
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! deltaFunction_Roma
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  function deltaFunction_Roma(ratio_in) result(delta); implicit none
    real(RK),intent(in)::ratio_in
    real(RK)::delta

    ! locals
    real(RK)::ratio

    ! Three-point type discrete Delta Function
    ! 
    ! [1] A.M. Roma, Ch. S. Peskin, M.J. Berger,  J. Comput. Phys. 153 (1999) . 
    !       An adaptive version of the immersed boundary method.
    ! [2] T. Kempe, J. Fröhlich, J. Comput. Phys. 231 (2012) 
    !       An improved immersed boundary method with direct forcing for the simulation of particle laden flows.
    ratio= abs(ratio_in)
    if(ratio> 1.50_RK) then
      delta= 0.0_RK
    elseif(ratio> 0.50_RK) then
      delta= 0.166666666666666666667_RK*(5.0_RK-3.0_RK*ratio-sqrt(1.0_RK-3.0_RK*(1.0_RK-ratio)**2))
    else
      delta= 0.333333333333333333333_RK*(1.0_RK+sqrt(1.0_RK-3.0_RK*ratio*ratio))
    endif
  end function deltaFunction_Roma

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! deltaFunction_Yang
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  function deltaFunction_Yang(ratio_in) result(delta); implicit none
    real(RK),intent(in)::ratio_in
    real(RK)::delta

    ! locals
    real(RK)::ratio
    ! 
    ! Three-point type discrete Delta Function 
    ! [1] Yang Xiaolei, Zhang Xing, Li Zhilin, He Guowei, JCP, 2009.
    ! 
    ratio= abs(ratio_in)
    if(ratio > 1.50_RK) then
      delta = 0.0_RK
    elseif(ratio > 0.50_RK) then
      delta = (ratio -3.0_RK)*ratio*0.5_RK +1.125_RK
    else
      delta = 0.75_RK -ratio*ratio
    endif
  end function deltaFunction_Yang
#endif
