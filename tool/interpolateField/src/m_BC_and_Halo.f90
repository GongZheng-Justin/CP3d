module m_BC_and_Halo
  use m_TypeDef
  use m_decomp2d
  use m_Parameters
#ifdef ScalarFlow
  use m_MeshAndMetries,only: ypOld
#endif
  implicit none
  private
  
  public:: SetBC_and_UpdateHalo_ux,SetBC_and_UpdateHalo_uy
  public:: SetBC_and_UpdateHalo_uz,SetBC_and_UpdateHalo_pr
#ifdef ScalarFlow
  public::SetBC_and_UpdateHalo_scalar
#endif
contains
    
  !******************************************************************
  ! SetBC_and_UpdateHalo_ux
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo_ux(ux,mbOld)
    implicit none
    type(MatBound),intent(in)::mbOld
    real(RK),dimension(mbOld%xmm:mbOld%xpm,mbOld%ymm:mbOld%ypm,mbOld%zmm:mbOld%zpm),intent(inout)::ux
    
    ! locals
    type(HaloInfo)::hiOld
    integer::ic,kc,nxp,nyp,nzp,nxc,nyc,nzc

    nxp=nxpOld;nyp=nypOld;nzp=nzpOld
    nxc=nxcOld;nyc=nycOld;nzc=nzcOld

    ! yp-dir
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          ux(ic,nyp,  kc) = uxBcValue(yp_dir)*two -ux(ic, nyc, kc)
        enddo
      enddo 
    CASE(BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          ux(ic, nyp,  kc) = ux(ic, nyc, kc)
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          ux(ic, 0, kc) = uxBcValue(ym_dir)*two-ux(ic, 1, kc)
        enddo
      enddo
    CASE(BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          ux(ic, 0, kc) = ux(ic, 1, kc)
        enddo
      enddo    
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC( 3))
    CASE(BC_NSLIP)
      ux(nxp,  0:nyp,ystart(3):yend(3))= uxBcValue(xp_dir)
    CASE(BC_FSLIP)
      ux(nxp,  0:nyp,ystart(3):yend(3))= zero      
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC( 4))
    CASE(BC_NSLIP)
      ux(1,0:nyp,ystart(3):yend(3)) = uxBcValue(xm_dir)      
    CASE(BC_FSLIP)
      ux(1,0:nyp,ystart(3):yend(3)) = zero       
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC( 1))
    CASE(BC_NSLIP) 
      ux(ystart(1):yend(1),0:nyp,nzp) = uxBcValue(zp_dir)*two -ux(ystart(1):yend(1), 0:nyp, nzc) 
    CASE(BC_FSLIP)
      ux(ystart(1):yend(1),0:nyp,nzp) = ux(ystart(1):yend(1), 0:nyp, nzc) 
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC( 2))
    CASE(BC_NSLIP)
      ux(ystart(1):yend(1),0:nyp,0) = uxBcValue(zm_dir)*two -ux(ystart(1):yend(1), 0:nyp, 1)
    CASE(BC_FSLIP)
      ux(ystart(1):yend(1),0:nyp,0) = ux(ystart(1):yend(1), 0:nyp, 1)    
    END SELECT     

    ! update halo
    hiOld%xmh=0;  hiOld%xph=1
    hiOld%zmh=1;  hiOld%zph=1
    if(BcOption(ym_dir)==BC_PERIOD) then
      hiOld%ymh=1;  hiOld%yph=1
    else
      hiOld%ymh=0;  hiOld%yph=0
    endif
    call decomp2d_updateHalo(ux, mbOld, hiOld)
  end subroutine SetBC_and_UpdateHalo_ux

  !******************************************************************
  ! SetBC_and_UpdateHalo_uy
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo_uy(uy,mbOld)
    implicit none
    type(MatBound),intent(in)::mbOld
    real(RK),dimension(mbOld%xmm:mbOld%xpm,mbOld%ymm:mbOld%ypm,mbOld%zmm:mbOld%zpm),intent(inout)::uy

    ! locals
    type(HaloInfo)::hiOld
    integer::ic,kc,nxp,nyp,nzp,nxc,nyc,nzc

    nxp=nxpOld;nyp=nypOld;nzp=nzpOld
    nxc=nxcOld;nyc=nycOld;nzc=nzcOld

    ! yp-dir
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uy(ic,nyp,  kc) = uyBcValue(yp_dir)
        enddo
      enddo 
    CASE(BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uy(ic, nyp,  kc) = zero
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uy(ic, 1, kc) = uyBcValue(ym_dir)
        enddo
      enddo
    CASE(BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uy(ic, 1, kc) = zero
        enddo
      enddo    
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(3))
    CASE(BC_NSLIP)
      uy(nxp,  0:nyp,ystart(3):yend(3))= uyBcValue(xp_dir)*two -uy(nxc, 0:nyp,ystart(3):yend(3))
    CASE(BC_FSLIP)
      uy(nxp,  0:nyp,ystart(3):yend(3))= uy(nxc, 0:nyp,ystart(3):yend(3))      
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(4))
    CASE(BC_NSLIP)
      uy(0,0:nyp,ystart(3):yend(3)) = uyBcValue(xm_dir)*two -uy(1,0:nyp,ystart(3):yend(3))     
    CASE(BC_FSLIP)
      uy(0,0:nyp,ystart(3):yend(3)) = uy(1,0:nyp,ystart(3):yend(3))       
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(1))
    CASE(BC_NSLIP)  
      uy(ystart(1):yend(1), 0:nyp, nzp  ) = uyBcValue(zp_dir)*two -uy(ystart(1):yend(1), 0:nyp, nzc)
    CASE(BC_FSLIP)
      uy(ystart(1):yend(1), 0:nyp, nzp  ) = uy(ystart(1):yend(1), 0:nyp, nzc)
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(2))
    CASE(BC_NSLIP)
      uy(ystart(1):yend(1), 0:nyp, 0) = uyBcValue(zm_dir)*two -uy(ystart(1):yend(1), 0:nyp, 1)
    CASE(BC_FSLIP)
      uy(ystart(1):yend(1), 0:nyp, 0) = uy(ystart(1):yend(1), 0:nyp, 1)      
    END SELECT

    ! update halo
    hiOld%xmh=1;  hiOld%xph=1
    hiOld%zmh=1;  hiOld%zph=1
    if(BcOption(ym_dir)==BC_PERIOD) then
      hiOld%ymh=0;  hiOld%yph=1
    else
      hiOld%ymh=0;  hiOld%yph=0
    endif
    call decomp2d_updateHalo(uy, mbOld, hiOld)
  end subroutine SetBC_and_UpdateHalo_uy

  !******************************************************************
  ! SetBC_and_UpdateHalo_uz
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo_uz(uz,mbOld)
    implicit none
    type(MatBound),intent(in)::mbOld
    real(RK),dimension(mbOld%xmm:mbOld%xpm,mbOld%ymm:mbOld%ypm,mbOld%zmm:mbOld%zpm),intent(inout)::uz

    ! locals
    type(HaloInfo)::hiOld
    integer::ic,kc,nxp,nyp,nzp,nxc,nyc,nzc

    nxp=nxpOld;nyp=nypOld;nzp=nzpOld
    nxc=nxcOld;nyc=nycOld;nzc=nzcOld

    ! yp-dir
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uz(ic,nyp,  kc) = uzBcValue(yp_dir)*two -uz(ic, nyc, kc)
        enddo
      enddo 
    CASE(BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uz(ic, nyp,  kc) = uz(ic, nyc, kc)
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uz(ic, 0, kc) = uzBcValue(ym_dir)*two-uz(ic, 1, kc)
        enddo
      enddo
    CASE(BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          uz(ic, 0, kc) = uz(ic, 1, kc)
        enddo
      enddo    
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(3))
    CASE(BC_NSLIP)
      uz(nxp,  0:nyp,ystart(3):yend(3))=  uzBcValue(xp_dir)*two -uz(nxc, 0:nyp,ystart(3):yend(3))
    CASE(BC_FSLIP)
      uz(nxp,  0:nyp,ystart(3):yend(3))= uz(nxc, 0:nyp,ystart(3):yend(3))        
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(4))
    CASE(BC_NSLIP)
      uz(0,0:nyp,ystart(3):yend(3)) = uzBcValue(xm_dir)*two -uz(1,0:nyp,ystart(3):yend(3))       
    CASE(BC_FSLIP)
      uz(0,0:nyp,ystart(3):yend(3)) = uz(1,0:nyp,ystart(3):yend(3))         
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(1))
    CASE(BC_NSLIP) 
      uz(ystart(1):yend(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
    CASE(BC_FSLIP)
      uz(ystart(1):yend(1), 0:nyp, nzp  ) = zero
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(2))
    CASE(BC_NSLIP)
      uz(ystart(1):yend(1), 0:nyp, 1) = uzBcValue(zm_dir)
    CASE(BC_FSLIP)
      uz(ystart(1):yend(1), 0:nyp, 1) = zero     
    END SELECT     

    ! update halo
    hiOld%xmh=1;  hiOld%xph=1
    hiOld%zmh=0;  hiOld%zph=1
    if(BcOption(ym_dir)==BC_PERIOD) then
      hiOld%ymh=1;  hiOld%yph=1
    else
      hiOld%ymh=0;  hiOld%yph=0
    endif
    call decomp2d_updateHalo(uz, mbOld, hiOld)
  end subroutine SetBC_and_UpdateHalo_uz

  !******************************************************************
  ! SetBC_and_UpdateHalo_pr
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo_pr(pressure,mbOld)
    implicit none
    type(MatBound),intent(in)::mbOld
    real(RK),dimension(mbOld%xmm:mbOld%xpm,mbOld%ymm:mbOld%ypm,mbOld%zmm:mbOld%zpm),intent(inout)::pressure

    ! locals
    type(HaloInfo)::hiOld
    integer::ic,kc,nxp,nyp,nzp,nxc,nyc,nzc

    nxp=nxpOld;nyp=nypOld;nzp=nzpOld
    nxc=nxcOld;nyc=nycOld;nzc=nzcOld

    ! yp-dir
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NSLIP, BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          pressure(ic, nyp, kc) = pressure(ic, nyc, kc)
        enddo
      enddo
    END SELECT         
        
    ! ym-dir
    SELECT CASE(BcOption(ym_dir))        
    CASE(BC_NSLIP, BC_FSLIP)
      do kc=ystart(3),yend(3)
        do ic=ystart(1),yend(1)
          pressure(ic, 0, kc) = pressure(ic, 1, kc)
        enddo
      enddo     
    END SELECT

    ! xp-dir
    SELECT CASE(myProcNghBC(3))
    CASE(BC_NSLIP, BC_FSLIP)
      pressure(nxp, 0:nyp,ystart(3):yend(3)) = pressure(nxc, 0:nyp,ystart(3):yend(3))        
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(4))
    CASE(BC_NSLIP, BC_FSLIP)
      pressure(0,0:nyp,ystart(3):yend(3)) = pressure(1,0:nyp,ystart(3):yend(3))       
    END SELECT      
        
    ! zp-dir
    SELECT CASE(myProcNghBC(1))
    CASE(BC_NSLIP, BC_FSLIP)
      pressure(ystart(1):yend(1), 0:nyp, nzp) = pressure(ystart(1):yend(1), 0:nyp, nzc)         
    END SELECT         
    
    ! zm-dir
    SELECT CASE(myProcNghBC(2))
    CASE(BC_NSLIP, BC_FSLIP)
      pressure(ystart(1):yend(1), 0:nyp, 0) = pressure(ystart(1):yend(1), 0:nyp, 1)         
    END SELECT

    ! update halo
    hiOld%xmh=1;  hiOld%xph=1
    hiOld%zmh=1;  hiOld%zph=1
    if(BcOption(ym_dir)==BC_PERIOD) then
      hiOld%ymh=1;  hiOld%yph=1
    else
      hiOld%ymh=0;  hiOld%yph=0
    endif
    call decomp2d_updateHalo(pressure, mbOld, hiOld)
  end subroutine SetBC_and_UpdateHalo_pr

#ifdef ScalarFlow
  !******************************************************************
  ! SetBC_and_UpdateHalo_scalar
  !****************************************************************** 
  subroutine SetBC_and_UpdateHalo_scalar(scalar,mbOld)
    implicit none
    type(MatBound),intent(in)::mbOld
    real(RK),dimension(mbOld%xmm:mbOld%xpm,mbOld%ymm:mbOld%ypm,mbOld%zmm:mbOld%zpm),intent(inout)::scalar
    
    ! locals
    type(HaloInfo)::hiOld
    integer::ic,kc,nyp,nyc
    real(RK)::rcoe1,rcoe2,dyc1,dycn
    
    nyp=nypOld
    nyc=nycOld
    dyc1=ypOld(2)-ypOld(1)
    dycn=ypOld(nypOld)-ypOld(nycOld)
    
    ! ym-dir
    rcoe1=one; rcoe2=one
    SELECT CASE(ScalarBcOption(1))
    CASE(-1) ! Dirichlet Bc. C     = F
      rcoe1= -one
      rcoe2=  two*ScalarBcValues(1)  
    CASE(-2) ! Neumann Bc.   dC/dy = S
      rcoe1=  one
      rcoe2= -dyc1*ScalarBcValues(3)
    CASE(-3) ! Robin Bc.     dC/dy + F*C = S
      rcoe1= (two+ScalarBcValues(1)*dyc1)/(two-ScalarBcValues(1)*dyc1)
      rcoe2= -two*ScalarBcValues(3)*dyc1 /(two-ScalarBcValues(1)*dyc1)
    END SELECT
    do kc=ystart(3),yend(3)
      do ic=ystart(1),yend(1)
        scalar(ic, 0, kc) = rcoe1*scalar(ic, 1, kc) +rcoe2
      enddo
    enddo

    ! yp-dir
    SELECT CASE(ScalarBcOption(1))
    CASE(-1)
      rcoe1= -one
      rcoe2=  two*ScalarBcValues(2)
    CASE(-2)
      rcoe1= one
      rcoe2= dycn*ScalarBcValues(4)
    CASE(-3)
      rcoe1= (two-ScalarBcValues(2)*dycn)/(two+ScalarBcValues(2)*dycn)
      rcoe2=  two*ScalarBcValues(4)*dycn /(two+ScalarBcValues(2)*dycn)
    END SELECT
    do kc=ystart(3),yend(3)
      do ic=ystart(1),yend(1)
        scalar(ic, nyp, kc) = rcoe1*scalar(ic, nyc, kc) +rcoe2
      enddo
    enddo

    ! update halo
    hiOld%xmh=1;  hiOld%xph=1
    hiOld%zmh=1;  hiOld%zph=1
    if(BcOption(ym_dir)==BC_PERIOD) then
      hiOld%ymh=1;  hiOld%yph=1
    else
      hiOld%ymh=0;  hiOld%yph=0
    endif
    call decomp2d_updateHalo(scalar, mbOld, hiOld)
  end subroutine SetBC_and_UpdateHalo_scalar
#endif 
end module m_BC_and_Halo
