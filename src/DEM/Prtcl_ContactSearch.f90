module Prtcl_ContactSearch
  use Prtcl_Parameters
  use Prtcl_NBS_Munjiza
  use Prtcl_Hrchl_Munjiza
  implicit none
  private
    
  type::ContactSearch
  contains
    procedure:: InitContactSearch => CS_InitContactSearch
    procedure:: FindContacts      => CS_FindContacts
    procedure:: get_numContact    => CS_get_numContact
  end type ContactSearch
  type(ContactSearch),public:: DEMContactSearch
    
contains

  !**********************************************************************
  ! Initializing particle-particle contact search object
  !**********************************************************************
  subroutine CS_InitContactSearch(this )
    implicit none
    class(ContactSearch)::this
    
    SELECT CASE(DEM_opt%CS_Method)
    CASE( CSM_NBS_Munjiza )
      allocate( m_NBS_Munjiza )
      call m_NBS_Munjiza%Init_NBSM()
    CASE( CSM_NBS_Munjiza_Hrchl )
      allocate( m_NBS_Munjiza_Hrchl)
      call m_NBS_Munjiza_Hrchl%Init_Munjiza_Hrchl()
    END SELECT
  end subroutine CS_InitContactSearch                            

  !**********************************************************************
  ! finding contact pairs of particles
  !**********************************************************************
  subroutine CS_FindContacts( this )
    implicit none
    class(ContactSearch)::this

    SELECT CASE( DEM_opt%CS_Method)
    CASE( CSM_NBS_Munjiza )
      call m_NBS_Munjiza%ContactSearch()
    CASE( CSM_NBS_Munjiza_Hrchl )
      call m_NBS_Munjiza_Hrchl%ContactSearch()
    END SELECT
  end subroutine CS_FindContacts

  !**********************************************************************
  ! CS_get_numContact
  !**********************************************************************
  function CS_get_numContact( this ) result (res)
    implicit none
    class(ContactSearch)::this
    integer,dimension(2):: res
    
    res = 0
    SELECT CASE( DEM_opt%CS_Method)
    CASE( CSM_NBS_Munjiza )
      res(1) = m_NBS_Munjiza%num_Cnsv_cntct
    CASE( CSM_NBS_Munjiza_Hrchl )
      res(1) = m_NBS_Munjiza_Hrchl%num_Cnsv_cntct
      res(2) = m_NBS_Munjiza_Hrchl%lvl_num_cnsv_cntct
    END SELECT
  end function CS_get_numContact

end module Prtcl_ContactSearch
