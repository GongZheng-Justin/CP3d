module m_BC_and_Halo
  use MPI
  use m_TypeDef
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1,hi1,OutFlowInfoX,OutFlowInfoY
  implicit none
  private
  type(HaloInfo)::hi_uxPrSrc,hi_uyPrSrc,hi_uzPrSrc  ! halo info type 1
    
  public:: Init_Halo, clcOutFlowVelocity, correctOutFlowFaceVelocity
  public:: SetBC_and_UpdateHalo, SetBC_and_UpdateHaloForPrSrc, SetBC_and_UpdateHalo_pr
contains

  !******************************************************************
  ! Init_Halo
  !******************************************************************
  subroutine Init_Halo()
    implicit none

    hi1%pencil = y_pencil
    hi1%xmh=1;  hi1%xph=1
    hi1%zmh=1;  hi1%zph=1
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi1%ymh=1;  hi1%yph=1
    else
      hi1%ymh=0;  hi1%yph=0
    endif
                
    hi_uxPrSrc%pencil = y_pencil
    hi_uxPrSrc%xmh=0;  hi_uxPrSrc%xph=1
    hi_uxPrSrc%ymh=0;  hi_uxPrSrc%yph=0
    hi_uxPrSrc%zmh=0;  hi_uxPrSrc%zph=0

    hi_uyPrSrc%pencil = y_pencil
    hi_uyPrSrc%xmh=0;  hi_uyPrSrc%xph=0;
    hi_uyPrSrc%zmh=0;  hi_uyPrSrc%zph=0;
    if(BcOption(ym_dir)==BC_PERIOD) then
      hi_uyPrSrc%ymh=0;  hi_uyPrSrc%yph=1
    else
      hi_uyPrSrc%ymh=0;  hi_uyPrSrc%yph=0
    endif

    hi_uzPrSrc%pencil = y_pencil
    hi_uzPrSrc%xmh=0;  hi_uzPrSrc%xph=0
    hi_uzPrSrc%ymh=0;  hi_uzPrSrc%yph=0
    hi_uzPrSrc%zmh=0;  hi_uzPrSrc%zph=1
  end subroutine Init_Halo
  
  !******************************************************************
  ! clcOutFlowVelocity
  !******************************************************************   
  subroutine clcOutFlowVelocity(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    
    ! locals
    integer::ic,jc,kc,im,jm,ierror
    real(RK)::ConVectVel,rTemp,convEd,rCoeC,rCoeP
    
    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyp; jm=jc-1
      ConVectVel=zero
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ConVectVel=ConVectVel+uy(ic,jc,kc)
        enddo
      enddo
      call MPI_ALLREDUCE(ConVectVel,rTemp,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      ConVectVel=rTemp/(real(nxc,RK)*real(nzc,RK))
      !print*, 'ConVectVel_y=',ConVectVel
      rCoeC=rdyc(jc); rCoeP=rdyp(jm)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          convEd=-ConVectVel*(ux(ic,jc,kc)-ux(ic,jm,kc))*rCoeC
          OutFlowInfoY(4,ic,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoY(1,ic,kc)
          ux(ic,jc,kc)=ux(ic,jc,kc)+OutFlowInfoY(4,ic,kc)
          OutFlowInfoY(1,ic,kc)=convEd

          convEd=-ConVectVel*(uy(ic,jc,kc)-uy(ic,jm,kc))*rCoeP
          OutFlowInfoY(5,ic,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoY(2,ic,kc)
          uy(ic,jc,kc)=uy(ic,jc,kc)+OutFlowInfoY(5,ic,kc)
          OutFlowInfoY(2,ic,kc)=convEd

          convEd=-ConVectVel*(uz(ic,jc,kc)-uz(ic,jm,kc))*rCoeC
          OutFlowInfoY(6,ic,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoY(3,ic,kc)
          uz(ic,jc,kc)=uz(ic,jc,kc)+OutFlowInfoY(6,ic,kc)
          OutFlowInfoY(3,ic,kc)=convEd
        enddo
      enddo    
    endif
    if(myProcNghBC(y_pencil,3)==BC_outFlow) then  
      ic=nxp; im=ic-1
      ConVectVel=zero
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          ConVectVel=ConVectVel+ux(ic,jc,kc)*dyp(jc)
        enddo
      enddo
      call MPI_ALLREDUCE(ConVectVel,rTemp,1,real_type,MPI_SUM,DECOMP_2D_COMM_ROW,ierror)
      ConVectVel=rTemp/(yly*real(nzc,RK))
      !print*, 'ConVectVel_x=',ConVectVel
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          convEd=-ConVectVel*(ux(ic,jc,kc)-ux(im,jc,kc))*rdx
          OutFlowInfoX(4,jc,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoX(1,jc,kc)
          ux(ic,jc,kc)=ux(ic,jc,kc)+OutFlowInfoX(4,jc,kc)
          OutFlowInfoX(1,jc,kc)=convEd        
        
          convEd=-ConVectVel*(uy(ic,jc,kc)-uy(im,jc,kc))*rdx
          OutFlowInfoX(5,jc,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoX(2,jc,kc)
          uy(ic,jc,kc)=uy(ic,jc,kc)+OutFlowInfoX(5,jc,kc)
          OutFlowInfoX(2,jc,kc)=convEd
          
          convEd=-ConVectVel*(uz(ic,jc,kc)-uz(im,jc,kc))*rdx
          OutFlowInfoX(6,jc,kc)=pmGamma*convEd+ pmTheta*OutFlowInfoX(3,jc,kc)
          uz(ic,jc,kc)=uz(ic,jc,kc)+OutFlowInfoX(6,jc,kc)
          OutFlowInfoX(3,jc,kc)=convEd           
        enddo
      enddo
    endif
  end subroutine clcOutFlowVelocity

  !******************************************************************
  ! correctOutFlowFaceVelocity
  !******************************************************************   
  subroutine correctOutFlowFaceVelocity(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    
    ! locals
    integer::ic,jc,kc,ierror
    real(RK)::CorrectVel,rTemp,OutFlowArea
    
    if(BcOption(xp_dir)/=BC_OutFlow .and. BcOption(yp_dir)/=BC_OutFlow) return
    CorrectVel=zero
    
    ! yp-dir
    if(BcOption(yp_dir)<0) then
      jc=nyp; rTemp=zero
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          rTemp=rTemp+uy(ic,jc,kc)
        enddo
      enddo
      CorrectVel=CorrectVel+rTemp*dx*dz
    endif
    
    ! ym-dir
    if(BcOption(ym_dir)<0) then
      jc=1; rTemp=zero
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          rTemp=rTemp+uy(ic,jc,kc)
        enddo
      enddo
      CorrectVel=CorrectVel-rTemp*dx*dz
    endif
    
    ! xp-dir
    if(myProcNghBC(y_pencil,3)<0) then
      ic=nxp; rTemp=zero
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          rTemp=rTemp+ux(ic,jc,kc)*dyp(jc)    
        enddo
      enddo
      CorrectVel=CorrectVel+rTemp*dz      
    endif
    
    ! xm-dir
    if(myProcNghBC(y_pencil,4)<0) then
      ic=1; rTemp=zero
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          rTemp=rTemp+ux(ic,jc,kc)*dyp(jc)    
        enddo
      enddo
      CorrectVel=CorrectVel-rTemp*dz       
    endif
    
    ! zp-dir
    if(myProcNghBC(y_pencil,1)<0) then
      kc=nzp; rtemp=zero
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          rTemp=rTemp+uz(ic,jc,kc)*dyp(jc)    
        enddo
      enddo
      CorrectVel=CorrectVel+rTemp*dx      
    endif
    
    ! zm-dir
    if(myProcNghBC(y_pencil,2)<0) then
      kc=1; rtemp=zero
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          rTemp=rTemp+uz(ic,jc,kc)*dyp(jc)    
        enddo
      enddo
      CorrectVel=CorrectVel-rTemp*dx      
    endif
    
    OutFlowArea=zero
    if(BcOption(xp_dir)==BC_outFlow) OutFlowArea=OutFlowArea+yly*zlz
    if(BcOption(yp_dir)==BC_outFlow) OutFlowArea=OutFlowArea+xlx*zlz
    call MPI_ALLREDUCE(CorrectVel,rTemp,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    CorrectVel=-rTemp/OutFlowArea
    
    if(BcOption(yp_dir)==BC_outFlow) then
      jc=nyp
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          uy(ic,jc,kc)=uy(ic,jc,kc)+CorrectVel
        enddo
      enddo    
    endif
    if(myProcNghBC(y_pencil,3)==BC_outFlow) then
      ic=nxp
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          ux(ic,jc,kc)=ux(ic,jc,kc)+CorrectVel
        enddo
      enddo       
    endif
  end subroutine correctOutFlowFaceVelocity
        
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
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NoSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,nyp,  kc)= uxBcValue(yp_dir)*two -ux(ic, nyc, kc)
          uy(ic,nyp,  kc)= uyBcValue(yp_dir)
          uy(ic,nyp+1,kc)= uyBcValue(yp_dir)*two -uy(ic, nyc, kc)
          uz(ic,nyp,  kc)= uzBcValue(yp_dir)*two -uz(ic, nyc, kc)
        enddo
      enddo 
    CASE(BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,nyp,  kc)= uxBcValue(yp_dir)*dyp(nyp) +ux(ic, nyc, kc)
          uy(ic,nyp,  kc)= uyBcValue(yp_dir)
          uy(ic,nyp+1,kc)= uyBcValue(yp_dir)*two      -uy(ic, nyc, kc)
          uz(ic,nyp,  kc)= uzBcValue(yp_dir)*dyp(nyp) +uz(ic, nyc, kc)
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NoSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic, 0, kc) = uxBcValue(ym_dir)*two-ux(ic, 1, kc)
          uy(ic, 1, kc) = uy_ym(ic,kc)
          uy(ic, 0, kc) = uy_ym(ic,kc)*two -uy(ic,2,kc)
          uz(ic, 0, kc) = uzBcValue(ym_dir)*two-uz(ic, 1, kc)
        enddo
      enddo
    CASE(BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic, 0, kc) = ux(ic, 1, kc)-uxBcValue(ym_dir)*dyp(1)
          uy(ic, 1, kc) = uy_ym(ic,kc)
          uy(ic, 0, kc) = uy_ym(ic,kc)*two -uy(ic,2,kc)
          uz(ic, 0, kc) = uz(ic, 1, kc)-uzBcValue(ym_dir)*dyp(1)
        enddo
      enddo    
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(y_pencil, 3))
    CASE(BC_NoSlip)
      ux(nxp,  0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)
      ux(nxp+1,0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)*two -ux(nxc, 0:nyp,y1start(3):y1end(3))
      uy(nxp,  0:nyp,y1start(3):y1end(3))= uyBcValue(xp_dir)*two -uy(nxc, 0:nyp,y1start(3):y1end(3))
      uz(nxp,  0:nyp,y1start(3):y1end(3))= uzBcValue(xp_dir)*two -uz(nxc, 0:nyp,y1start(3):y1end(3))
    CASE(BC_FreeSlip)
      ux(nxp,  0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)
      ux(nxp+1,0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)*two -ux(nxc, 0:nyp,y1start(3):y1end(3))
      uy(nxp,  0:nyp,y1start(3):y1end(3))= uy(nxc, 0:nyp,y1start(3):y1end(3)) +uyBcValue(xp_dir)*dx
      uz(nxp,  0:nyp,y1start(3):y1end(3))= uz(nxc, 0:nyp,y1start(3):y1end(3)) +uzBcValue(xp_dir)*dx        
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(y_pencil, 4))
    CASE(BC_NoSlip)
      ux(1,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)
      ux(0,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)*two -ux(2,0:nyp,y1start(3):y1end(3)) 
      uy(0,0:nyp,y1start(3):y1end(3)) = uyBcValue(xm_dir)*two -uy(1,0:nyp,y1start(3):y1end(3))
      uz(0,0:nyp,y1start(3):y1end(3)) = uzBcValue(xm_dir)*two -uz(1,0:nyp,y1start(3):y1end(3))       
    CASE(BC_FreeSlip)
      ux(1,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)
      ux(0,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)*two -ux(2,0:nyp,y1start(3):y1end(3)) 
      uy(0,0:nyp,y1start(3):y1end(3)) = uy(1,0:nyp,y1start(3):y1end(3)) -uyBcValue(xm_dir)*dx
      uz(0,0:nyp,y1start(3):y1end(3)) = uz(1,0:nyp,y1start(3):y1end(3)) -uzBcValue(xm_dir)*dx         
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(y_pencil, 1))
    CASE(BC_NoSlip) 
      ux(y1start(1):y1end(1), 0:nyp, nzp  ) = uxBcValue(zp_dir)*two -ux(y1start(1):y1end(1), 0:nyp, nzc) 
      uy(y1start(1):y1end(1), 0:nyp, nzp  ) = uyBcValue(zp_dir)*two -uy(y1start(1):y1end(1), 0:nyp, nzc)
      uz(y1start(1):y1end(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
      uz(y1start(1):y1end(1), 0:nyp, nzp+1) = uzBcValue(zp_dir)*two -uz(y1start(1):y1end(1), 0:nyp, nzc)
    CASE(BC_FreeSlip)
      ux(y1start(1):y1end(1), 0:nyp, nzp  ) = ux(y1start(1):y1end(1), 0:nyp, nzc) +uxBcValue(zp_dir)*dz 
      uy(y1start(1):y1end(1), 0:nyp, nzp  ) = uy(y1start(1):y1end(1), 0:nyp, nzc) +uyBcValue(zp_dir)*dz 
      uz(y1start(1):y1end(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
      uz(y1start(1):y1end(1), 0:nyp, nzp+1) = uzBcValue(zp_dir)*two -uz(y1start(1):y1end(1), 0:nyp, nzc)
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(y_pencil, 2))
    CASE(BC_NoSlip)
      ux(y1start(1):y1end(1), 0:nyp, 0) = uxBcValue(zm_dir)*two -ux(y1start(1):y1end(1), 0:nyp, 1)
      uy(y1start(1):y1end(1), 0:nyp, 0) = uyBcValue(zm_dir)*two -uy(y1start(1):y1end(1), 0:nyp, 1)
      uz(y1start(1):y1end(1), 0:nyp, 1) = uzBcValue(zm_dir)
      uz(y1start(1):y1end(1), 0:nyp, 0) = uzBcValue(zm_dir)*two -uz(y1start(1):y1end(1), 0:nyp, 2)
    CASE(BC_FreeSlip)
      ux(y1start(1):y1end(1), 0:nyp, 0) = ux(y1start(1):y1end(1), 0:nyp, 1) -uxBcValue(zm_dir)*dz 
      uy(y1start(1):y1end(1), 0:nyp, 0) = uy(y1start(1):y1end(1), 0:nyp, 1) -uyBcValue(zm_dir)*dz 
      uz(y1start(1):y1end(1), 0:nyp, 1) = uzBcValue(zm_dir)
      uz(y1start(1):y1end(1), 0:nyp, 0) = uzBcValue(zm_dir)*two -uz(y1start(1):y1end(1), 0:nyp, 2)     
    END SELECT     

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
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NoSlip,BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          uy(ic,nyp,kc) = uyBcValue(yp_dir)
        enddo
      enddo
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NoSlip,BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          uy(ic, 1, kc) = uy_ym(ic,kc)
        enddo
      enddo
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(y_pencil,3))
    CASE(BC_NoSlip,BC_FreeSlip)
      ux(nxp,0:nyp,y1start(3):y1end(3))=  uxBcValue(xp_dir)    
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(y_pencil,4))
    CASE(BC_NoSlip,BC_FreeSlip)
      ux(1,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)        
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(y_pencil,1))
    CASE(BC_NoSlip,BC_FreeSlip)
      uz(y1start(1):y1end(1),0:nyp,nzp) = uzBcValue(zp_dir)
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(y_pencil,2))
    CASE(BC_NoSlip,BC_FreeSlip)
      uz(y1start(1):y1end(1),0:nyp,1) = uzBcValue(zm_dir)  
    END SELECT

    ! update halo
    call myupdate_halo(ux, mb1, hi_uxPrSrc)
    call myupdate_halo(uy, mb1, hi_uyPrSrc)
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

    ! yp-dir
    if(BcOption(yp_dir)<0) then
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          pressure(ic, nyp, kc) = pressure(ic, nyc, kc)
        enddo
      enddo
    endif         
        
    ! ym-dir
    if(BcOption(ym_dir)<0) then
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          pressure(ic, 0, kc) = pressure(ic, 1, kc)
        enddo
      enddo     
    endif

    ! xp-dir
    if(myProcNghBC(y_pencil,3)<0) then
      pressure(nxp,0:nyp,y1start(3):y1end(3))= pressure(nxc,0:nyp,y1start(3):y1end(3))        
    endif
    
    ! xm-dir
    if(myProcNghBC(y_pencil,4)<0) then
      pressure(0,0:nyp,y1start(3):y1end(3))  = pressure(1,0:nyp,y1start(3):y1end(3))       
    endif      
        
    ! zp-dir
    if(myProcNghBC(y_pencil,1)<0) then
      pressure(y1start(1):y1end(1),0:nyp,nzp)= pressure(y1start(1):y1end(1),0:nyp,nzc)         
    endif        
    
    ! zm-dir
    if(myProcNghBC(y_pencil,2)<0) then
      pressure(y1start(1):y1end(1),0:nyp,0)  = pressure(y1start(1):y1end(1),0:nyp,1)       
    endif
    
    call myupdate_halo(pressure, mb1, hi1)
  end subroutine SetBC_and_UpdateHalo_pr
end module m_BC_and_Halo
