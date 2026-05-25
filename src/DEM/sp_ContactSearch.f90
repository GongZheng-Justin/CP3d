#include "definitions_inc.f90"
module sp_ContactSearch
  use sp_Parameters
  use sp_NBS_Munjiza
  use sp_Hrchl_Munjiza
  implicit none
  private
    
  public:: get_num_Cnsv_cntct_PP, InitContactSearchPP, FindContactsPP  
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Initializing particle-particle contact search object
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitContactSearchPP(); implicit none
    
    SELECT CASE(DEM_Opt%CS_Method)
    CASE( CSM_NBS_Munjiza )
      allocate( m_NBS_Munjiza )
      call m_NBS_Munjiza%Init_NBSM()
    CASE( CSM_NBS_Munjiza_Hrchl )
      allocate( m_NBS_Munjiza_Hrchl)
      call m_NBS_Munjiza_Hrchl%Init_Munjiza_Hrchl()
    END SELECT
  end subroutine InitContactSearchPP                           

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! finding contact pairs of particles
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine FindContactsPP(); implicit none

    SELECT CASE( DEM_Opt%CS_Method)
    CASE( CSM_NBS_Munjiza )
      call m_NBS_Munjiza%ContactSearch()
    CASE( CSM_NBS_Munjiza_Hrchl )
      call m_NBS_Munjiza_Hrchl%ContactSearch()
    END SELECT
  end subroutine FindContactsPP

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! get_num_Cnsv_cntct_PP
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine get_num_Cnsv_cntct_PP(res); implicit none
    integer,dimension(2),intent(out):: res
    
    res = 0
    SELECT CASE( DEM_Opt%CS_Method)
    CASE( CSM_NBS_Munjiza )
      res(1) = m_NBS_Munjiza%num_Cnsv_cntct
    CASE( CSM_NBS_Munjiza_Hrchl )
      res(1) = m_NBS_Munjiza_Hrchl%num_Cnsv_cntct
      res(2) = m_NBS_Munjiza_Hrchl%lvl_num_cnsv_cntct
    END SELECT
  end subroutine get_num_Cnsv_cntct_PP

end module sp_ContactSearch
