#include "definitions_inc.f90"
  do k=0,nDistribute  ! ux grid
    kd = k+idzc_interp
    prz= RatioZc(k)*VolRatio
    do j=0,nDistribute
      jd  = j+idyc_interp
      pry = RatioYc(j)*prz
      do i=0,nDistribute
        id = i+idxp_interp
        prx= RatioXp(i)
#ifdef spreadIBM_force_inc_unity
        VolForce_x(id,jd,kd)= VolForce_x(id,jd,kd)+ prx*pry
#else
        VolForce_x(id,jd,kd)= VolForce_x(id,jd,kd)+ prx*pry* IbpForce%x
#endif
      enddo
    enddo
  enddo
      
  do k=0,nDistribute  ! uy gird
    kd = k+idzc_interp
    prz= RatioZc(k)*VolRatio
    do j=0,nDistribute
      jd = j+idyp_interp
      pry= RatioYp(j)*prz
      do i=0,nDistribute
        id = i+idxc_interp
        prx= RatioXc(i)
#ifdef spreadIBM_force_inc_unity
        VolForce_y(id,jd,kd)= VolForce_y(id,jd,kd)+ prx*pry
#else
        VolForce_y(id,jd,kd)= VolForce_y(id,jd,kd)+ prx*pry* IbpForce%y
#endif
      enddo
    enddo
  enddo

  do k=0,nDistribute  ! uz grid
    kd = k+idzp_interp
    prz= RatioZp(k)*VolRatio
    do j=0,nDistribute
      jd = j+idyc_interp
      pry= RatioYc(j)*prz
      do i=0,nDistribute
        id = i+idxc_interp
        prx= RatioXc(i)
#ifdef spreadIBM_force_inc_unity
        VolForce_z(id,jd,kd)= VolForce_z(id,jd,kd)+ prx*pry
#else
        VolForce_z(id,jd,kd)= VolForce_z(id,jd,kd)+ prx*pry* IbpForce%z
#endif
      enddo
    enddo
  enddo
