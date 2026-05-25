#include "definitions_inc.f90"
#ifdef IBMDistributeLinear
  prxp= (LagrangeP%x-xc(idxp_interp))*rdx+0.5_RK
  RatioXp(0)= 1.0_RK-prxp
  RatioXp(1)= prxp
  prxc= (LagrangeP%x-xc(idxc_interp))*rdx
  RatioXc(0)= 1.0_RK-prxc
  RatioXc(1)= prxc

  pryp= (LagrangeP%y-yc(idyp_interp))*rdyUniform+0.5_RK
  RatioYp(0)= 1.0_RK-pryp
  RatioYp(1)= pryp
  pryc= (LagrangeP%y-yc(idyc_interp))*rdyUniform
  RatioYc(0)= 1.0_RK-pryc
  RatioYc(1)= pryc

  przp= (LagrangeP%z-zc(idzp_interp))*rdz+0.5_RK
  RatioZp(0)= 1.0_RK-przp
  RatioZp(1)= przp
  przc= (LagrangeP%z-zc(idzc_interp))*rdz
  RatioZc(0)= 1.0_RK-przc
  RatioZc(1)= przc
#else
  prxp= (LagrangeP%x-xc(idxp_interp))*rdx+0.5_RK
  RatioXp(0)= deltaFunction(prxp)
  RatioXp(1)= deltaFunction(prxp-1.0_RK)
  RatioXp(2)= 1.0_RK- RatioXp(0)- RatioXp(1)
  prxc= (LagrangeP%x-xc(idxc_interp))*rdx
  RatioXc(0)= deltaFunction(prxc)
  RatioXc(1)= deltaFunction(prxc-1.0_RK)
  RatioXc(2)= 1.0_RK- RatioXc(0)- RatioXc(1)

  pryp= (LagrangeP%y-yc(idyp_interp))*rdyUniform+0.5_RK
  RatioYp(0)= deltaFunction(pryp)
  RatioYp(1)= deltaFunction(pryp-1.0_RK)
  RatioYp(2)= 1.0_RK- RatioYp(0)- RatioYp(1)
  pryc= (LagrangeP%y-yc(idyc_interp))*rdyUniform
  RatioYc(0)= deltaFunction(pryc)
  RatioYc(1)= deltaFunction(pryc-1.0_RK)
  RatioYc(2)= 1.0_RK- RatioYc(0)- RatioYc(1)

  przp= (LagrangeP%z-zc(idzp_interp))*rdz+0.5_RK
  RatioZp(0)= deltaFunction(przp)
  RatioZp(1)= deltaFunction(przp-1.0_RK)
  RatioZp(2)= 1.0_RK- RatioZp(0)- RatioZp(1)
  przc= (LagrangeP%z-zc(idzc_interp))*rdz
  RatioZc(0)= deltaFunction(przc)
  RatioZc(1)= deltaFunction(przc-1.0_RK)
  RatioZc(2)= 1.0_RK- RatioZc(0)- RatioZc(1)
#endif
