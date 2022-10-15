module m_BC_and_Halo
  use m_TypeDef
  use m_LogInfo
  use m_decomp2d
  use m_Parameters  
#ifdef ScalarFlow
  use m_MeshAndMetries,only: dyc,rdyc
#endif
  use m_Variables,only: mb1
  implicit none
  private
  type(HaloInfo)::hi1,hi_uxPrSrc,hi_uzPrSrc,hi_pr
  
  public:: Init_Halo, SetBC_and_UpdateHalo_pr
  public:: SetBC_and_UpdateHalo,SetBC_and_UpdateHaloForPrSrc
#ifdef ScalarFlow
  public::SetBC_and_UpdateHalo_scalar
#endif
contains

  !******************************************************************
  ! Init_Halo
  !******************************************************************
  subroutine Init_Halo()
    implicit none

    hi1%pencil = y_pencil
    hi1%xmh=2;  hi1%xph=2
    hi1%ymh=0;  hi1%yph=0
    hi1%zmh=2;  hi1%zph=2
        
    hi_uxPrSrc%pencil = y_pencil
    hi_uxPrSrc%xmh=1;  hi_uxPrSrc%xph=2
    hi_uxPrSrc%ymh=0;  hi_uxPrSrc%yph=0
    hi_uxPrSrc%zmh=0;  hi_uxPrSrc%zph=0

    hi_uzPrSrc%pencil = y_pencil
    hi_uzPrSrc%xmh=0;  hi_uzPrSrc%xph=0
    hi_uzPrSrc%ymh=0;  hi_uzPrSrc%yph=0
    hi_uzPrSrc%zmh=1;  hi_uzPrSrc%zph=2
    
    hi_pr%pencil = y_pencil
    hi_pr%xmh=2;  hi_pr%xph=1
    hi_pr%ymh=0;  hi_pr%yph=0
    hi_pr%zmh=2;  hi_pr%zph=1
  end subroutine Init_Halo
      
  !******************************************************************
  ! SetBC_and_UpdateHalo
  !******************************************************************
  subroutine SetBC_and_UpdateHalo(ux,uy,uz,uy_ym)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in)::uy_ym
    
    ! locals
    integer::ic,kc

    ! yp-dir
    SELECT CASE(FlowType)
    CASE(FT_CH)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,nyp,  kc) = uxBcValue(yp_dir)*two -ux(ic, nyc, kc)
          uy(ic,nyp,  kc) = uyBcValue(yp_dir)
          uz(ic,nyp,  kc) = uzBcValue(yp_dir)*two -uz(ic, nyc, kc)
        enddo
      enddo 
    CASE(FT_HC)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic, nyp,  kc) = ux(ic, nyc, kc)
          uy(ic, nyp,  kc) = zero
          uz(ic, nyp,  kc) = uz(ic, nyc, kc)
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    do kc=y1start(3),y1end(3)
      do ic=y1start(1),y1end(1)
        ux(ic, 0, kc) = uxBcValue(ym_dir)*two-ux(ic, 1, kc)
        uy(ic, 1, kc) = uy_ym(ic,kc)
        uz(ic, 0, kc) = uzBcValue(ym_dir)*two-uz(ic, 1, kc)
      enddo
    enddo

    ! update halo
    call myupdate_halo(ux, mb1, hi1)
    call myupdate_halo(uy, mb1, hi1)
    call myupdate_halo(uz, mb1, hi1) 
  end subroutine SetBC_and_UpdateHalo

  !******************************************************************
  ! SetBC_and_UpdateHaloForPrSrc
  !******************************************************************   
  subroutine SetBC_and_UpdateHaloForPrSrc(ux,uy,uz,uy_ym)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in)::uy_ym
    
    ! locals
    integer::ic,kc

    ! yp-dir
    SELECT CASE(FlowType)
    CASE(FT_CH)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          uy(ic,nyp,  kc) = uyBcValue(yp_dir)
        enddo
      enddo 
    CASE(FT_HC)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          uy(ic, nyp,  kc) = zero
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    do kc=y1start(3),y1end(3)
      do ic=y1start(1),y1end(1)
        uy(ic, 1, kc) = uy_ym(ic,kc)
      enddo
    enddo

    ! update halo
    call myupdate_halo(ux, mb1, hi_uxPrSrc)
    call myupdate_halo(uz, mb1, hi_uzPrSrc)
  end subroutine SetBC_and_UpdateHaloForPrSrc
  
  !******************************************************************
  ! SetBC_and_UpdateHalo_pr
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo_pr(pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure

    ! locals
    integer::ic,kc    
        
    do kc=y1start(3),y1end(3)
      do ic=y1start(1),y1end(1)
        pressure(ic, 0,   kc) = pressure(ic, 1,   kc)
        pressure(ic, nyp, kc) = pressure(ic, nyc, kc)
      enddo
    enddo
    call myupdate_halo(pressure, mb1, hi_pr)
  end subroutine SetBC_and_UpdateHalo_pr
  
#ifdef ScalarFlow
  !******************************************************************
  ! SetBC_and_UpdateHalo_scalar
  !****************************************************************** 
  subroutine SetBC_and_UpdateHalo_scalar(scalar)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::scalar
    
    ! locals
    integer::ic,kc
    real(RK)::rcoe1,rcoe2
    
    ! ym-dir
    SELECT CASE(ScalarBcOption(1))
    CASE(-1) ! Dirichlet Bc. C     = F
      rcoe1= -one
      rcoe2=  two*ScalarBcValues(1)  
    CASE(-2) ! Neumann Bc.   dC/dy = S
      rcoe1=  one
      rcoe2= -dyc(1)*ScalarBcValues(3)
    CASE(-3) ! Robin Bc.     dC/dy + F*C = S
      rcoe1= (two+ScalarBcValues(1)*dyc(1))/(two-ScalarBcValues(1)*dyc(1))
      rcoe2= -two*ScalarBcValues(3)*dyc(1) /(two-ScalarBcValues(1)*dyc(1))
    CASE DEFAULT
      rcoe1=zero; rcoe2=zero
      call MainLog%CheckForError(ErrT_Abort,"SetBC_and_UpdateHalo_scalar","ScalarBcOption(1) wrong !!!")
    END SELECT
    do kc=y1start(3),y1end(3)
      do ic=y1start(1),y1end(1)
        scalar(ic, 0, kc) = rcoe1*scalar(ic, 1, kc) +rcoe2
      enddo
    enddo

    ! yp-dir
    SELECT CASE(ScalarBcOption(2))
    CASE(-1)
      rcoe1= -one
      rcoe2=  two*ScalarBcValues(2)
    CASE(-2)
      rcoe1= one
      rcoe2= dyc(nyp)*ScalarBcValues(4)
    CASE(-3)
      rcoe1= (two-ScalarBcValues(2)*dyc(nyp))/(two+ScalarBcValues(2)*dyc(nyp))
      rcoe2=  two*ScalarBcValues(4)*dyc(nyp) /(two+ScalarBcValues(2)*dyc(nyp))
    CASE DEFAULT
      rcoe1=zero; rcoe2=zero
      call MainLog%CheckForError(ErrT_Abort,"SetBC_and_UpdateHalo_scalar","ScalarBcOption(2) wrong !!!")
    END SELECT
    do kc=y1start(3),y1end(3)
      do ic=y1start(1),y1end(1)
        scalar(ic, nyp, kc) = rcoe1*scalar(ic, nyc, kc) +rcoe2
      enddo
    enddo
    call myupdate_halo(scalar, mb1, hi1)
  end subroutine SetBC_and_UpdateHalo_scalar
#endif 
end module m_BC_and_Halo
