module ca_BC_and_Halo
  use MPI
  use m_TypeDef
  use m_Parameters
  use Prtcl_Decomp_2d
  use m_MeshAndMetries
  use m_Variables,only:mb1
  use m_Decomp2d,only:nrank,y1start,y1end,HaloInfo,myProcNghBC,update_halo,sum_halo
  implicit none
  private

  public:: SetBC_and_UpdateHalo_VelIBM,Gather_Halo_IBMForce
contains
  !******************************************************************
  ! SetBC_and_UpdateHalo_VelIBM
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo_VelIBM(ux,uy,uz,uy_ym)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in)::uy_ym
    
    ! locals
    integer::ic,kc
    type(HaloInfo)::hi_interp

   ! yp-dir
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NoSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,nyp,  kc)= uxBcValue(yp_dir)*2.0_RK -ux(ic, nyc, kc)
          uy(ic,nyp,  kc)= uyBcValue(yp_dir)
          uy(ic,nyp+1,kc)= uyBcValue(yp_dir)*2.0_RK -uy(ic, nyc, kc)
          uz(ic,nyp,  kc)= uzBcValue(yp_dir)*2.0_RK -uz(ic, nyc, kc)
        enddo
      enddo 
    CASE(BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,nyp,  kc)= uxBcValue(yp_dir)*dyp(nyp) +ux(ic, nyc, kc)
          uy(ic,nyp,  kc)= uyBcValue(yp_dir)
          uy(ic,nyp+1,kc)= uyBcValue(yp_dir)*2.0_RK   -uy(ic, nyc, kc)
          uz(ic,nyp,  kc)= uzBcValue(yp_dir)*dyp(nyp) +uz(ic, nyc, kc)
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NoSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic, 0, kc) = uxBcValue(ym_dir)*2.0_RK-ux(ic, 1, kc)
          uy(ic, 1, kc) = uy_ym(ic,kc)
          uy(ic, 0, kc) = uy_ym(ic,kc)*2.0_RK -uy(ic,2,kc)
          uz(ic, 0, kc) = uzBcValue(ym_dir)*2.0_RK-uz(ic, 1, kc)
        enddo
      enddo
    CASE(BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic, 0, kc) = ux(ic, 1, kc)-uxBcValue(ym_dir)*dyp(1)
          uy(ic, 1, kc) = uy_ym(ic,kc)
          uy(ic, 0, kc) = uy_ym(ic,kc)*2.0_RK -uy(ic,2,kc)
          uz(ic, 0, kc) = uz(ic, 1, kc)-uzBcValue(ym_dir)*dyp(1)
        enddo
      enddo    
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(y_pencil, 3))
    CASE(BC_NoSlip)
      ux(nxp,  0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)
      ux(nxp+1,0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)*2.0_RK -ux(nxc, 0:nyp,y1start(3):y1end(3))
      uy(nxp,  0:nyp,y1start(3):y1end(3))= uyBcValue(xp_dir)*2.0_RK -uy(nxc, 0:nyp,y1start(3):y1end(3))
      uz(nxp,  0:nyp,y1start(3):y1end(3))= uzBcValue(xp_dir)*2.0_RK -uz(nxc, 0:nyp,y1start(3):y1end(3))
    CASE(BC_FreeSlip)
      ux(nxp,  0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)
      ux(nxp+1,0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)*2.0_RK -ux(nxc, 0:nyp,y1start(3):y1end(3))
      uy(nxp,  0:nyp,y1start(3):y1end(3))= uy(nxc, 0:nyp,y1start(3):y1end(3)) +uyBcValue(xp_dir)*dx
      uz(nxp,  0:nyp,y1start(3):y1end(3))= uz(nxc, 0:nyp,y1start(3):y1end(3)) +uzBcValue(xp_dir)*dx        
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(y_pencil, 4))
    CASE(BC_NoSlip)
      ux(1,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)
      ux(0,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)*2.0_RK -ux(2,0:nyp,y1start(3):y1end(3)) 
      uy(0,0:nyp,y1start(3):y1end(3)) = uyBcValue(xm_dir)*2.0_RK -uy(1,0:nyp,y1start(3):y1end(3))
      uz(0,0:nyp,y1start(3):y1end(3)) = uzBcValue(xm_dir)*2.0_RK -uz(1,0:nyp,y1start(3):y1end(3))       
    CASE(BC_FreeSlip)
      ux(1,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)
      ux(0,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)*2.0_RK -ux(2,0:nyp,y1start(3):y1end(3)) 
      uy(0,0:nyp,y1start(3):y1end(3)) = uy(1,0:nyp,y1start(3):y1end(3)) -uyBcValue(xm_dir)*dx
      uz(0,0:nyp,y1start(3):y1end(3)) = uz(1,0:nyp,y1start(3):y1end(3)) -uzBcValue(xm_dir)*dx         
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(y_pencil, 1))
    CASE(BC_NoSlip) 
      ux(y1start(1):y1end(1), 0:nyp, nzp  ) = uxBcValue(zp_dir)*2.0_RK -ux(y1start(1):y1end(1), 0:nyp, nzc) 
      uy(y1start(1):y1end(1), 0:nyp, nzp  ) = uyBcValue(zp_dir)*2.0_RK -uy(y1start(1):y1end(1), 0:nyp, nzc)
      uz(y1start(1):y1end(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
      uz(y1start(1):y1end(1), 0:nyp, nzp+1) = uzBcValue(zp_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, nzc)
    CASE(BC_FreeSlip)
      ux(y1start(1):y1end(1), 0:nyp, nzp  ) = ux(y1start(1):y1end(1), 0:nyp, nzc) +uxBcValue(zp_dir)*dz 
      uy(y1start(1):y1end(1), 0:nyp, nzp  ) = uy(y1start(1):y1end(1), 0:nyp, nzc) +uyBcValue(zp_dir)*dz 
      uz(y1start(1):y1end(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
      uz(y1start(1):y1end(1), 0:nyp, nzp+1) = uzBcValue(zp_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, nzc)
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(y_pencil, 2))
    CASE(BC_NoSlip)
      ux(y1start(1):y1end(1), 0:nyp, 0) = uxBcValue(zm_dir)*2.0_RK -ux(y1start(1):y1end(1), 0:nyp, 1)
      uy(y1start(1):y1end(1), 0:nyp, 0) = uyBcValue(zm_dir)*2.0_RK -uy(y1start(1):y1end(1), 0:nyp, 1)
      uz(y1start(1):y1end(1), 0:nyp, 1) = uzBcValue(zm_dir)
      uz(y1start(1):y1end(1), 0:nyp, 0) = uzBcValue(zm_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, 2)
    CASE(BC_FreeSlip)
      ux(y1start(1):y1end(1), 0:nyp, 0) = ux(y1start(1):y1end(1), 0:nyp, 1) -uxBcValue(zm_dir)*dz 
      uy(y1start(1):y1end(1), 0:nyp, 0) = uy(y1start(1):y1end(1), 0:nyp, 1) -uyBcValue(zm_dir)*dz 
      uz(y1start(1):y1end(1), 0:nyp, 1) = uzBcValue(zm_dir)
      uz(y1start(1):y1end(1), 0:nyp, 0) = uzBcValue(zm_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, 2)     
    END SELECT     

    ! update halo
    hi_interp%pencil = y_pencil
    hi_interp%ymh=0;  hi_interp%yph=0
#ifdef IBMDistribute2
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_interp%ymh=1;  hi_interp%yph=1
    endif
    hi_interp%xmh=0;  hi_interp%xph=1
    hi_interp%zmh=1;  hi_interp%zph=1
    call update_halo(ux, mb1, hi_interp)

    hi_interp%xmh=1;  hi_interp%xph=1
    hi_interp%zmh=0;  hi_interp%zph=1
    call update_halo(uz, mb1, hi_interp)

    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_interp%ymh=0;  hi_interp%yph=1    
    endif    
    hi_interp%xmh=1;  hi_interp%xph=1
    hi_interp%zmh=1;  hi_interp%zph=1    
    call update_halo(uy, mb1, hi_interp)
#else
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_interp%ymh=1;  hi_interp%yph=1    
    endif
    hi_interp%xmh=1;  hi_interp%xph=2
    hi_interp%zmh=1;  hi_interp%zph=1
    call update_halo(ux, mb1, hi_interp)

    hi_interp%xmh=1;  hi_interp%xph=1
    hi_interp%zmh=1;  hi_interp%zph=2
    call update_halo(uz, mb1, hi_interp)

    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_interp%ymh=1;  hi_interp%yph=2
    endif    
    hi_interp%xmh=1;  hi_interp%xph=1
    hi_interp%zmh=1;  hi_interp%zph=1    
    call update_halo(uy, mb1, hi_interp)
#endif
  end subroutine SetBC_and_UpdateHalo_VelIBM

  !******************************************************************
  ! Gather_Halo_IBMForce
  !******************************************************************
  subroutine Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::VolForce_x,VolForce_y,VolForce_z

    ! locals
    integer::ic,kc
    type(HaloInfo)::hi_Force

    hi_Force%pencil = y_pencil    
    hi_Force%ymh=0;  hi_Force%yph=0
#ifdef IBMDistribute2
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_Force%ymh=1;  hi_Force%yph=1    
    endif
    hi_Force%xmh=0;  hi_Force%xph=1
    hi_Force%zmh=1;  hi_Force%zph=1
    call sum_halo(VolForce_x,mb1,hi_Force)   

    hi_Force%xmh=1;  hi_Force%xph=1
    hi_Force%zmh=0;  hi_Force%zph=1
    call sum_halo(VolForce_z,mb1,hi_Force)
        
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_Force%ymh=0;  hi_Force%yph=1    
    endif
    hi_Force%xmh=1;  hi_Force%xph=1
    hi_Force%zmh=1;  hi_Force%zph=1
    call sum_halo(VolForce_y,mb1,hi_Force) 
#else
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_Force%ymh=1;  hi_Force%yph=1    
    endif
    hi_Force%xmh=1;  hi_Force%xph=2
    hi_Force%zmh=1;  hi_Force%zph=1
    call sum_halo(VolForce_x,mb1,hi_Force)   

    hi_Force%xmh=1;  hi_Force%xph=1
    hi_Force%zmh=1;  hi_Force%zph=2
    call sum_halo(VolForce_z,mb1,hi_Force)
        
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_Force%ymh=1;  hi_Force%yph=2    
    endif
    hi_Force%xmh=1;  hi_Force%xph=1
    hi_Force%zmh=1;  hi_Force%zph=1
    call sum_halo(VolForce_y,mb1,hi_Force)
#endif    
    
#ifdef ReserveIBMForce
    if(BcOption(ym_dir)/=BC_PERIOD) then
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          VolForce_x(ic,  1,kc)=VolForce_x(ic,  1,kc)+VolForce_x(ic,    0,k);  VolForce_x(ic,  0,  kc)=0.0_RK
          VolForce_x(ic,nyc,kc)=VolForce_x(ic,nyc,kc)+VolForce_x(ic,  nyp,k);  VolForce_x(ic,nyp,  kc)=0.0_RK
          
          VolForce_z(ic,  1,kc)=VolForce_z(ic,  1,kc)+VolForce_z(ic,    0,k);  VolForce_z(ic,  0,  kc)=0.0_RK        
          VolForce_z(ic,nyc,kc)=VolForce_z(ic,nyc,kc)+VolForce_z(ic,  nyp,k);  VolForce_z(ic,nyp,  kc)=0.0_RK
          
          VolForce_y(ic,  1,kc)=VolForce_y(ic,  1,kc)+VolForce_y(ic,    0,k);  VolForce_y(ic,    0,kc)=0.0_RK
          VolForce_y(ic,nyp,kc)=VolForce_y(ic,nyp,kc)+VolForce_y(ic,nyp+1,k);  VolForce_y(ic,nyp+1,kc)=0.0_RK          
        enddo
      enddo
    endif
#endif
  end subroutine Gather_Halo_IBMForce
end module ca_BC_and_Halo
