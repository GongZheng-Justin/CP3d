
  ! This file is included in the module m_TScheme

  !******************************************************************
  ! clcRhsX_FEXP
  !******************************************************************
  subroutine  clcRhsX_FEXP(ux,uy,uz,RhsX,HistXold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsX
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistXold
    
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp,ierror
    real(RK)::qdx1,qdx3,h11,h12,h13,sucaj,dp1ns,Forced,Forced_tot,ForcedCoe,ForcedXnu
    real(RK)::d11q1,d22q1,d33q1,dcq123,convEd1,gradp1,InterpY1,InterpY2,InterpY3,InterpY4
   
    Forced=0.0_RK 
    qdx1=0.25_RK*rdx
    qdx3=0.25_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        sucaj= 0.5_RK*rdyp(jc)
        ForcedCoe=dyp(jc); ForcedXnu=ForcedCoe*xnu 
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=1.0_RK-InterpY3
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1

          h11=( (ux(ip,jc,kc)+ux(ic,jc,kc))* (ux(ip,jc,kc)+ux(ic,jc,kc))  &
               -(ux(ic,jc,kc)+ux(im,jc,kc))* (ux(ic,jc,kc)+ux(im,jc,kc)) )*qdx1 
          h12=( (uy(ic,jp,kc)+uy(im,jp,kc))* (InterpY3*ux(ic,jc,kc) +InterpY4*ux(ic,jp,kc))  &
               -(uy(ic,jc,kc)+uy(im,jc,kc))* (InterpY1*ux(ic,jm,kc) +InterpY2*ux(ic,jc,kc)) )*sucaj            
          h13=( (uz(ic,jc,kp)+uz(im,jc,kp))* (ux(ic,jc,kp)+ux(ic,jc,kc))  &
               -(uz(ic,jc,kc)+uz(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jc,km)) )*qdx3
                
          d11q1= (ux(ip,jc,kc)-2.0_RK*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          d22q1= ap2c(jc)*ux(ic,jp,kc) +ac2c(jc)*ux(ic,jc,kc) +am2c(jc)*ux(ic,jm,kc)                
          d33q1= ap3c(kc)*ux(ic,jc,kp) +ac3c(kc)*ux(ic,jc,kc) +am3c(kc)*ux(ic,jc,km)
          dcq123= d11q1+d22q1+d33q1
             
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
#ifdef CFDDEM
          Forced_tot= 0.5_RK*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
          Forced=Forced +ForcedXnu*dcq123 +ForcedCoe*Forced_tot
          convEd1= -h11-h12-h13 + xnu*dcq123 + gravity(1) +Forced_tot
#elif CFDLPT_TwoWay
          Forced_tot= FpForce_x(ic,jc,kc)
          Forced=Forced +ForcedXnu*dcq123 +ForcedCoe*Forced_tot
          convEd1= -h11-h12-h13 + xnu*dcq123 + gravity(1) +Forced_tot
#else
          Forced=Forced +ForcedXnu*dcq123
          convEd1= -h11-h12-h13 + xnu*dcq123 + gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1 +pmTheta*HistXold(ic,jc,kc) -pmAlpha*gradp1
          HistXold(ic,jc,kc)=convEd1
        enddo
      enddo
    ENDDO
  
    ! in dp1ns there is the mean pressure gradient to keep constant mass
    IF(IsUxConst) THEN
      call MPI_ALLREDUCE(Forced,Forced_tot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      Forced= -Forced_tot/(real(nxc*nzc,kind=RK))/yly
      dp1ns = pmGamma*(Forced-gravity(1)) +pmTheta*PrGradData(4)
      PrGradData(4) =Forced-gravity(1)
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
            RhsX(ic,jc,kc)=RhsX(ic,jc,kc) +dp1ns
          enddo
        enddo
      enddo
      PrGradData(3)= dp1ns/pmAlpha
      PrGradData(1)= PrGradData(1) +dp1ns/dt
    ENDIF
  end subroutine clcRhsX_FEXP

  !******************************************************************
  ! clcRhsY_FEXP
  !******************************************************************    
  subroutine  clcRhsY_FEXP(ux,uy,uz,RhsY,HistYold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsY
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistYold
   
    ! locals 
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::hdx1,hdx3,h21,h22,h23,sucac,qsucac
    real(RK)::d11q2,d22q2,d33q2,dcq123,convEd2,gradp2,InterpY1,InterpY2    
    
    hdx1=0.5_RK*rdx
    hdx3=0.5_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        sucac = rdyc(jc)
        qsucac= 0.25_RK*sucac
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1

          h21=( (InterpY1*ux(ip,jm,kc)+InterpY2*ux(ip,jc,kc))* (uy(ip,jc,kc)+uy(ic,jc,kc)) &
               -(InterpY1*ux(ic,jm,kc)+InterpY2*ux(ic,jc,kc))* (uy(ic,jc,kc)+uy(im,jc,kc)) )*hdx1
          h22=( (uy(ic,jp,kc)+uy(ic,jc,kc))* (uy(ic,jp,kc)+uy(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jm,kc)) )*qsucac
          h23=( (InterpY1*uz(ic,jm,kp)+InterpY2*uz(ic,jc,kp))* (uy(ic,jc,kp)+uy(ic,jc,kc)) &
               -(InterpY1*uz(ic,jm,kc)+InterpY2*uz(ic,jc,kc))* (uy(ic,jc,kc)+uy(ic,jc,km)) )*hdx3
                
          d11q2= ap1c(ic)*uy(ip,jc,kc)+ac1c(ic)*uy(ic,jc,kc)+am1c(ic)*uy(im,jc,kc)
          d22q2= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          d33q2= ap3c(kc)*uy(ic,jc,kp)+ac3c(kc)*uy(ic,jc,kc)+am3c(kc)*uy(ic,jc,km)
          dcq123= d11q2+d22q2+d33q2
       
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
#ifdef CFDDEM
          convEd2= -h21-h22-h23+xnu*dcq123+ gravity(2) +InterpY1*FpForce_y(ic,jm,kc)+InterpY2*FpForce_y(ic,jc,kc)
#elif CFDLPT_TwoWay
          convEd2= -h21-h22-h23+xnu*dcq123+ gravity(2) +FpForce_y(ic,jc,kc)
#else
          convEd2= -h21-h22-h23+xnu*dcq123+ gravity(2)
#endif
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2
          HistYold(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO
  end subroutine clcRhsY_FEXP

  !******************************************************************
  ! clcRhsZ_FEXP
  !****************************************************************** 
  subroutine  clcRhsZ_FEXP(ux,uy,uz,RhsZ,HistZold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistZold
   
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::d11q3,d22q3,d33q3,dcq123,convEd3,gradp3
    real(RK)::qdx1,qdx3,h31,h32,h33,sucaj,InterpY1,InterpY2,InterpY3,InterpY4
    
    qdx1=0.25_RK*rdx
    qdx3=0.25_RK*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        sucaj =0.5_RK*rdyp(jc)
        InterpY1= YinterpCoe(jc); InterpY2=1.0_RK-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=1.0_RK-InterpY3
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1

          h31=( (ux(ip,jc,kc)+ux(ip,jc,km))* (uz(ip,jc,kc)+uz(ic,jc,kc)) &
               -(ux(ic,jc,kc)+ux(ic,jc,km))* (uz(ic,jc,kc)+uz(im,jc,kc)) )*qdx1
          h32=( (uy(ic,jp,kc)+uy(ic,jp,km))* (InterpY3*uz(ic,jc,kc) +InterpY4*uz(ic,jp,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jc,km))* (InterpY1*uz(ic,jm,kc) +InterpY2*uz(ic,jc,kc)) )*sucaj                
          h33=( (uz(ic,jc,kp)+uz(ic,jc,kc))* (uz(ic,jc,kp)+uz(ic,jc,kc)) &
               -(uz(ic,jc,kc)+uz(ic,jc,km))* (uz(ic,jc,kc)+uz(ic,jc,km)) )*qdx3

          d11q3= ap1c(ic)*uz(ip,jc,kc)+ac1c(ic)*uz(ic,jc,kc)+am1c(ic)*uz(im,jc,kc)
          d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
          d33q3= (uz(ic,jc,kp)-2.0_RK*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          dcq123= d11q3+d22q3+d33q3
                
          gradp3= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz                
#ifdef CFDDEM
          convEd3= -h31-h32-h33+xnu*dcq123+ gravity(3)+0.5_RK*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#elif CFDLPT_TwoWay
          convEd3= -h31-h32-h33+xnu*dcq123+ gravity(3)+FpForce_z(ic,jc,kc)
#else
          convEd3= -h31-h32-h33+xnu*dcq123+ gravity(3)
#endif
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3
          HistZold(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO 
  end subroutine clcRhsZ_FEXP

  !******************************************************************
  ! clcU1Hat_FEXP
  !******************************************************************    
  subroutine clcU1Hat_FEXP(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1) 
          ux(ic,jc,kc)= ux(ic,jc,kc)+ RhsX(ic,jc,kc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FEXP
  
  !******************************************************************
  ! clcU2Hat_FEXP
  !******************************************************************  
  subroutine clcU2Hat_FEXP(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in):: duy_ym
    
    ! locals
    integer::ic,jc,kc

    DO kc=y1start(3),y1end(3) 
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1) 
          uy(ic,jc,kc)= uy(ic,jc,kc)+ RhsY(ic,jc,kc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_FEXP  
  
  !******************************************************************
  ! clcU3Hat_FEXP
  !******************************************************************  
  subroutine clcU3Hat_FEXP(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz
  
    ! locals
    integer::ic,jc,kc
    
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+ RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FEXP

  !******************************************************************
  ! PressureUpdate_FEXP
  !******************************************************************
  subroutine PressureUpdate_FEXP(pressure, prphiHalo)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout):: pressure
    
    ! locals
    integer::ic,jc,kc

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          pressure(ic,jc,kc)= pressure(ic,jc,kc)+ prphiHalo(ic,jc,kc)
        enddo
      enddo
    ENDDO
  end subroutine PressureUpdate_FEXP
