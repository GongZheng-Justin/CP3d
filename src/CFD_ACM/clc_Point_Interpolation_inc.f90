#include "definitions_inc.f90"
  SumXDir=0.0_RK  ! ux grid
  do k=0,nDistribute
    kd = k+idzc_interp
    prz= RatioZc(k)
    do j=0,nDistribute
      jd  = j+idyc_interp
      pry = RatioYc(j)*prz
      do i=0,nDistribute
        id = i+idxp_interp
        prx= RatioXp(i)
#ifdef spreadIBM_force_inc_unity
        SumXDir= SumXDir + VolForce_x(id,jd,kd)*prx*pry
#else
        SumXDir= SumXDir + ux(id,jd,kd)*prx*pry
#endif
      enddo
    enddo
  enddo
      
  SumYDir=0.0_RK  ! uy gird
  do k=0,nDistribute
    kd = k+idzc_interp
    prz= RatioZc(k)
    do j=0,nDistribute
      jd = j+idyp_interp
      pry= RatioYp(j)*prz
      do i=0,nDistribute
        id = i+idxc_interp
        prx= RatioXc(i)
#ifdef spreadIBM_force_inc_unity
        SumYDir= SumYDir + VolForce_y(id,jd,kd)*prx*pry
#else
        SumYDir= SumYDir + uy(id,jd,kd)*prx*pry
#endif
      enddo
    enddo
  enddo
      
  SumZDir=0.0_RK  ! uz grid
  do k=0,nDistribute
    kd = k+idzp_interp
    prz= RatioZp(k)
    do j=0,nDistribute
      jd = j+idyc_interp
      pry= RatioYc(j)*prz
      do i=0,nDistribute
        id = i+idxc_interp
        prx= RatioXc(i)
#ifdef spreadIBM_force_inc_unity
        SumZDir= SumZDir + VolForce_z(id,jd,kd)*prx*pry
#else
        SumZDir= SumZDir + uz(id,jd,kd)*prx*pry
#endif
      enddo
    enddo
  enddo
