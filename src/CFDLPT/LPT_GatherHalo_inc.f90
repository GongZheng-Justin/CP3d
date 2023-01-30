  !******************************************************************
  ! Gather_Halo_dist_1
  !******************************************************************
  subroutine Gather_Halo_dist_1()
    implicit none
    type(HaloInfo)::hi_Force
    
    ! Force_x
    hi_Force%pencil = y_pencil
    hi_Force%xmh=0;  hi_Force%xph=1
    hi_Force%ymh=0;  hi_Force%yph=0
    hi_Force%zmh=1;  hi_Force%zph=1
    call sum_halo(FpForce_x,mb1,hi_Force)   
    
    ! Force_y
    hi_Force%pencil = y_pencil
    hi_Force%xmh=1;  hi_Force%xph=1
    hi_Force%ymh=0;  hi_Force%yph=0
    hi_Force%zmh=1;  hi_Force%zph=1
    call sum_halo(FpForce_y,mb1,hi_Force)

    ! Force_z
    hi_Force%pencil = y_pencil
    hi_Force%xmh=1;  hi_Force%xph=1
    hi_Force%ymh=0;  hi_Force%yph=0
    hi_Force%zmh=0;  hi_Force%zph=1
    call sum_halo(FpForce_z,mb1,hi_Force)
  end subroutine Gather_Halo_dist_1
  
  !******************************************************************
  ! Gather_Halo_dist_2
  !******************************************************************
  subroutine Gather_Halo_dist_2()
    implicit none
    type(HaloInfo)::hi_Force
    
    ! Force_x
    hi_Force%pencil = y_pencil
    hi_Force%xmh=1;  hi_Force%xph=2
    hi_Force%ymh=0;  hi_Force%yph=0
    hi_Force%zmh=1;  hi_Force%zph=1
    call sum_halo(FpForce_x,mb1,hi_Force)   
    
    ! Force_y
    hi_Force%pencil = y_pencil
    hi_Force%xmh=1;  hi_Force%xph=1
    hi_Force%ymh=0;  hi_Force%yph=0
    hi_Force%zmh=1;  hi_Force%zph=1
    call sum_halo(FpForce_y,mb1,hi_Force)

    ! Force_z
    hi_Force%pencil = y_pencil
    hi_Force%xmh=1;  hi_Force%xph=1
    hi_Force%ymh=0;  hi_Force%yph=0
    hi_Force%zmh=1;  hi_Force%zph=2
    call sum_halo(FpForce_z,mb1,hi_Force)
  end subroutine Gather_Halo_dist_2
