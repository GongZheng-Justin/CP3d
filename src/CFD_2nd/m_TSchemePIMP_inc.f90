
  ! This file is included in the module m_TScheme

  !******************************************************************
  ! clcRhsX_PIMP
  !******************************************************************    
  subroutine  clcRhsX_PIMP(ux,uy,uz,RhsX,HistXold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsX
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistXold
    
    ! locals    
    integer::im,ic,ip,jc,jm,jp,km,kc,kp,ierror
    real(RK)::qdx1,qdx3,h11,h12,h13,sucaj,s1tot,s1tot1,s1totCoe,dpmdxns
    real(RK)::d11q1,d22q1,d33q1,dcq13,convEd1,gradp1,InterpY1,InterpY2,InterpY3,InterpY4
    
    s1tot=zero
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        s1totCoe=dyp(jc)
        sucaj=half*rdyp(jc)
        if(yc(jc)<ybulk1 .or. yc(jc)>=ybulk2) s1totCoe=zero
        InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=one-InterpY3            
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1

          h11=( (ux(ip,jc,kc)+ux(ic,jc,kc))* (ux(ip,jc,kc)+ux(ic,jc,kc))  &
               -(ux(ic,jc,kc)+ux(im,jc,kc))* (ux(ic,jc,kc)+ux(im,jc,kc)) )*qdx1
          h12=( (uy(ic,jp,kc)+uy(im,jp,kc))* (InterpY3*ux(ic,jc,kc) +InterpY4*ux(ic,jp,kc))  &
               -(uy(ic,jc,kc)+uy(im,jc,kc))* (InterpY1*ux(ic,jm,kc) +InterpY2*ux(ic,jc,kc)) )*sucaj           
          h13=( (uz(ic,jc,kp)+uz(im,jc,kp))* (ux(ic,jc,kp)+ux(ic,jc,kc))  &
               -(uz(ic,jc,kc)+uz(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jc,km)) )*qdx3
                
          d11q1= (ux(ip,jc,kc)-two*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          d22q1= ap2c(jc)*ux(ic,jp,kc) + ac2c(jc)*ux(ic,jc,kc)+ am2c(jc)*ux(ic,jm,kc)                
          d33q1= ap3c(kc)*ux(ic,jc,kp) + ac3c(kc)*ux(ic,jc,kc)+ am3c(kc)*ux(ic,jc,km)
          dcq13= d11q1+d33q1
             
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
          s1tot=s1tot+ux(ic,jc,kc)*s1totCoe
#ifdef CFDDEM
          convEd1= -h11-h12-h13 + xnu*dcq13 + gravity(1)+half*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
#elif CFDLPT_TwoWay
          convEd1= -h11-h12-h13 + xnu*dcq13 + gravity(1)+FpForce_x(ic,jc,kc)
#else
          convEd1= -h11-h12-h13 + xnu*dcq13 + gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1+ pmTheta*HistXold(ic,jc,kc)- pmAlpha*gradp1+ two*pmBeta*d22q1
          HistXold(ic,jc,kc)=convEd1
        enddo
      enddo
    ENDDO
  
    ! in dp1ns there is the mean pressure gradient to keep constant mass
    IF(IsUxConst) THEN
      call MPI_ALLREDUCE(s1tot,s1tot1,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      s1tot1=s1tot1/(real(nxc*nzc,kind=RK))/(ybulk2-ybulk1)
      dpmdxns= pmAlphaC*(ubulk - s1tot1)
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
            RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+ dpmdxns
          enddo
        enddo
      enddo
      PrGradNow=(ubulk - s1tot1)/dt
      PrGradAver = PrGradAver+ pmAlphaC*PrGradNow
    ENDIF
  end subroutine clcRhsX_PIMP

  !******************************************************************
  ! clcRhsX_PIMP_LES
  !******************************************************************    
  subroutine  clcRhsX_PIMP_LES(ux,uy,uz,RhsX,HistXold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsX
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistXold
    
    ! locals    
    integer::im,ic,ip,jc,jm,jp,km,kc,kp,ierror
    real(RK)::visa,visb,visc,visd,sgs1,sgs2,sgs3
    real(RK)::d11q1,d22q1,d33q1,dcq13,convEd1,gradp1
    real(RK)::qdx1,qdx3,h11,h12,h13,sucaj,s1tot,s1tot1,s1totCoe,dpmdxns
    
    s1tot=zero
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1; kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1; jp=jc+1
        s1totCoe=dyp(jc)
        sucaj=quarter*rdyp(jc)
        if(yc(jc)<ybulk1 .or. yc(jc)>=ybulk2) s1totCoe=zero 
        do ic=y1start(1),y1end(1)
          im=ic-1
          ip=ic+1
          h11=( (ux(ip,jc,kc)+ux(ic,jc,kc))* (ux(ip,jc,kc)+ux(ic,jc,kc))  &
               -(ux(ic,jc,kc)+ux(im,jc,kc))* (ux(ic,jc,kc)+ux(im,jc,kc)) )*qdx1

          h12=( (uy(ic,jp,kc)+uy(im,jp,kc))* (ux(ic,jp,kc)+ux(ic,jc,kc))  &
               -(uy(ic,jc,kc)+uy(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jm,kc)) )*sucaj            
          h13=( (uz(ic,jc,kp)+uz(im,jc,kp))* (ux(ic,jc,kp)+ux(ic,jc,kc))  &
               -(uz(ic,jc,kc)+uz(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jc,km)) )*qdx3

          visa=quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc))
          visb=quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jp,kc)+nut(im,jp,kc))
          visc=quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,km)+nut(im,jc,km))
          visd=quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jc,kp)+nut(im,jc,kp))
          sgs1=rdx2*(nut(ic,jc,kc)*(ux(ip,jc,kc)-ux(ic,jc,kc))-nut(im,jc,kc)*(ux(ic,jc,kc)-ux(im,jc,kc)))
          sgs2=rdx*rdyp(jc)*(visb*(uy(ic,jp,kc)-uy(im,jp,kc))-visa*(uy(ic,jc,kc)-uy(im,jc,kc)))
          sgs3=rdx*rdz     *(visd*(uz(ic,jc,kp)-uz(im,jc,kp))-visc*(uz(ic,jc,kc)-uz(im,jc,kc)))
                
          d11q1= sgs1
          d22q1= ap2c(jc)*visb*(ux(ic,jp,kc)-ux(ic,jc,kc))-am2c(jc)*visa*(ux(ic,jc,kc)-ux(ic,jm,kc))
          d33q1= ap3c(kc)*visd*(ux(ic,jc,kp)-ux(ic,jc,kc))-am3c(kc)*visc*(ux(ic,jc,kc)-ux(ic,jc,km))
          dcq13= d11q1+d33q1
                
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
          s1tot=s1tot+ux(ic,jc,kc)*s1totCoe        
#ifdef CFDDEM
          convEd1= -h11-h12-h13 +sgs1+sgs2+sgs3+dcq13+gravity(1)+half*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
#elif CFDLPT_TwoWay
          convEd1= -h11-h12-h13 +sgs1+sgs2+sgs3+dcq13+gravity(1)+FpForce_x(ic,jc,kc)
#else
          convEd1= -h11-h12-h13 +sgs1+sgs2+sgs3+dcq13+gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1+ pmTheta*HistXold(ic,jc,kc)- pmAlpha*gradp1+ two*pmBetaT*d22q1
          HistXold(ic,jc,kc)=convEd1
        enddo
      enddo
    ENDDO
  
    ! in dp1ns there is the mean pressure gradient to keep constant mass
    IF(IsUxConst) THEN
      call MPI_ALLREDUCE(s1tot,s1tot1,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      s1tot1=s1tot1/(real(nxc*nzc,kind=RK))/(ybulk2-ybulk1)
      dpmdxns= pmAlphaC*(ubulk - s1tot1)
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
            RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+ dpmdxns
          enddo
        enddo
      enddo
      PrGradNow=(ubulk - s1tot1)/dt
      PrGradAver = PrGradAver+ pmAlphaC*PrGradNow
    ENDIF
  end subroutine clcRhsX_PIMP_LES

  !******************************************************************
  ! clcRhsY_PIMP
  !******************************************************************    
  subroutine  clcRhsY_PIMP(ux,uy,uz,RhsY,HistYold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure 
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsY
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistYold
   
    ! locals 
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::hdx1,hdx3,h21,h22,h23,sucac,qsucac
    real(RK)::d11q2,d22q2,d33q2,dcq13,convEd2,gradp2,InterpY1,InterpY2     
    
    hdx1=half*rdx
    hdx3=half*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1; kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1; jp=jc+1
        sucac = rdyc(jc)
        qsucac= quarter*sucac
        InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1
        do ic=y1start(1),y1end(1)
          im=ic-1; ip=ic+1

          h21=( (InterpY1*ux(ip,jm,kc)+InterpY2*ux(ip,jc,kc))* (uy(ip,jc,kc)+uy(ic,jc,kc)) &
               -(InterpY1*ux(ic,jm,kc)+InterpY2*ux(ic,jc,kc))* (uy(ic,jc,kc)+uy(im,jc,kc)) )*hdx1
          h22=( (uy(ic,jp,kc)+uy(ic,jc,kc))* (uy(ic,jp,kc)+uy(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jm,kc)) )*qsucac
          h23=( (InterpY1*uz(ic,jm,kp)+InterpY2*uz(ic,jc,kp))* (uy(ic,jc,kp)+uy(ic,jc,kc)) &
               -(InterpY1*uz(ic,jm,kc)+InterpY2*uz(ic,jc,kc))* (uy(ic,jc,kc)+uy(ic,jc,km)) )*hdx3
                
          d11q2= ap1c(ic)*uy(ip,jc,kc)+ac1c(ic)*uy(ic,jc,kc)+am1c(ic)*uy(im,jc,kc)
          d22q2= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          d33q2= ap3c(kc)*uy(ic,jc,kp)+ac3c(kc)*uy(ic,jc,kc)+am3c(kc)*uy(ic,jc,km)
          dcq13= d11q2+d33q2
       
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
#ifdef CFDDEM
          convEd2= -h21-h22-h23+xnu*dcq13+ gravity(2) +InterpY1*FpForce_y(ic,jm,kc)+InterpY2*FpForce_y(ic,jc,kc)
#elif CFDLPT_TwoWay
          convEd2= -h21-h22-h23+xnu*dcq13+ gravity(2) +FpForce_y(ic,jc,kc)
#else
          convEd2= -h21-h22-h23+xnu*dcq13+ gravity(2)
#endif
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2+ two*pmBeta*d22q2
          HistYold(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO
  end subroutine clcRhsY_PIMP

  !******************************************************************
  ! clcRhsY_PIMP_LES
  !******************************************************************    
  subroutine  clcRhsY_PIMP_LES(ux,uy,uz,RhsY,HistYold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsY
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistYold
   
    ! locals 
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::qdx1,qdx3,h21,h22,h23,sucac,qsucac
    real(RK)::d11q2,d22q2,d33q2,dcq13,convEd2,gradp2  
    real(RK)::visa,visb,visc,visd,sgs1,sgs2,sgs3  
    
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1
      kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1
        jp=jc+1
        sucac= rdyc(jc)
        qsucac=quarter*sucac
        do ic=y1start(1),y1end(1)
          im=ic-1                
          ip=ic+1
          h21=( (ux(ip,jc,kc)+ux(ip,jm,kc))* (uy(ip,jc,kc)+uy(ic,jc,kc)) &
               -(ux(ic,jc,kc)+ux(ic,jm,kc))* (uy(ic,jc,kc)+uy(im,jc,kc)) )*qdx1
          h22=( (uy(ic,jp,kc)+uy(ic,jc,kc))* (uy(ic,jp,kc)+uy(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jm,kc)) )*qsucac
          h23=( (uz(ic,jc,kp)+uz(ic,jm,kp))* (uy(ic,jc,kp)+uy(ic,jc,kc)) &
               -(uz(ic,jc,kc)+uz(ic,jm,kc))* (uy(ic,jc,kc)+uy(ic,jc,km)) )*qdx3

          visa=quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc))
          visb=quarter*(nut(ic,jc,kc)+nut(ip,jc,kc)+nut(ic,jm,kc)+nut(ip,jm,kc))
          visc=quarter*(nut(ic,jc,kc)+nut(ic,jm,kc)+nut(ic,jc,km)+nut(ic,jm,km))
          visd=quarter*(nut(ic,jc,kc)+nut(ic,jm,kc)+nut(ic,jc,kp)+nut(ic,jm,kp))
          sgs1=rdx*sucac*(visb*(ux(ip,jc,kc)-ux(ip,jm,kc))-visa*(ux(ic,jc,kc)-ux(ic,jm,kc)))
          sgs2=ap2p(jc)*nut(ic,jc,kc)*(uy(ic,jp,kc)-uy(ic,jc,kc))-am2p(jc)*nut(ic,jm,kc)*(uy(ic,jc,kc)-uy(ic,jm,kc))
          sgs3=rdz*sucac*(visd*(uz(ic,jc,kp)-uz(ic,jm,kp))-visc*(uz(ic,jc,kc)-uz(ic,jm,kc)))

          d11q2= ap1c(ic)*visb*(uy(ip,jc,kc)-uy(ic,jc,kc) ) -am1c(ic)*visa*(uy(ic,jc,kc)- uy(im,jc,kc)) 
          d22q2= sgs2
          d33q2= ap3c(kc)*visd*(uy(ic,jc,kp)-uy(ic,jc,kc) ) -am3c(kc)*visc*(uy(ic,jc,kc)-uy(ic,jc,km))
          dcq13= d11q2+d33q2
                
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
#ifdef CFDDEM
          convEd2= -h21-h22-h23+sgs1+sgs2+sgs3+dcq13+ gravity(2)+half*(FpForce_y(ic,jc,kc)+FpForce_y(ic,jm,kc))
#elif CFDLPT_TwoWay
          convEd2= -h21-h22-h23+sgs1+sgs2+sgs3+dcq13+ gravity(2)+FpForce_y(ic,jc,kc)
#else
          convEd2= -h21-h22-h23+sgs1+sgs2+sgs3+dcq13+ gravity(2)
#endif
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2+ two*pmBetaT*d22q2
          HistYold(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO
  end subroutine clcRhsY_PIMP_LES

  !******************************************************************
  ! clcRhsZ_PIMP
  !****************************************************************** 
  subroutine  clcRhsZ_PIMP(ux,uy,uz,RhsZ,HistZold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistZold
   
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::d11q3,d22q3,d33q3,dcq13,convEd3,gradp3
    real(RK)::qdx1,qdx3,h31,h32,h33,sucaj,InterpY1,InterpY2,InterpY3,InterpY4
    
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1; kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1; jp=jc+1
        sucaj=half*rdyp(jc)
        InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=one-InterpY3
        do ic=y1start(1),y1end(1)
          im=ic-1; ip=ic+1

          h31=( (ux(ip,jc,kc)+ux(ip,jc,km))* (uz(ip,jc,kc)+uz(ic,jc,kc)) &
               -(ux(ic,jc,kc)+ux(ic,jc,km))* (uz(ic,jc,kc)+uz(im,jc,kc)) )*qdx1
          h32=( (uy(ic,jp,kc)+uy(ic,jp,km))* (InterpY3*uz(ic,jc,kc) +InterpY4*uz(ic,jp,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jc,km))* (InterpY1*uz(ic,jm,kc) +InterpY2*uz(ic,jc,kc)) )*sucaj                
          h33=( (uz(ic,jc,kp)+uz(ic,jc,kc))* (uz(ic,jc,kp)+uz(ic,jc,kc)) &
               -(uz(ic,jc,kc)+uz(ic,jc,km))* (uz(ic,jc,kc)+uz(ic,jc,km)) )*qdx3

          d11q3= ap1c(ic)*uz(ip,jc,kc)+ac1c(ic)*uz(ic,jc,kc)+am1c(ic)*uz(im,jc,kc)
          d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
          d33q3= (uz(ic,jc,kp)-two*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          dcq13= d11q3+d33q3
                
          gradp3= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz                
#ifdef CFDDEM
          convEd3= -h31-h32-h33+xnu*dcq13+ gravity(3)+half*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#elif CFDLPT_TwoWay
          convEd3= -h31-h32-h33+xnu*dcq13+ gravity(3)+FpForce_z(ic,jc,kc)
#else
          convEd3= -h31-h32-h33+xnu*dcq13+ gravity(3)
#endif
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3+ two*pmBeta*d22q3
          HistZold(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO 
  end subroutine clcRhsZ_PIMP

  !******************************************************************
  ! clcRhsZ_PIMP_LES
  !******************************************************************
  subroutine  clcRhsZ_PIMP_LES(ux,uy,uz,RhsZ,HistZold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistZold
   
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::qdx1,qdx3,h31,h32,h33,sucaj
    real(RK)::d11q3,d22q3,d33q3,dcq13,convEd3,gradp3
    real(RK)::visa,visb,visc,visd,sgs1,sgs2,sgs3
    
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1
      kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1
        jp=jc+1
        sucaj=quarter*rdyp(jc)
        do ic=y1start(1),y1end(1)
          im=ic-1
          ip=ic+1
          h31=( (ux(ip,jc,kc)+ux(ip,jc,km))* (uz(ip,jc,kc)+uz(ic,jc,kc)) &
               -(ux(ic,jc,kc)+ux(ic,jc,km))* (uz(ic,jc,kc)+uz(im,jc,kc)) )*qdx1
          h32=( (uy(ic,jp,kc)+uy(ic,jp,km))* (uz(ic,jp,kc)+uz(ic,jc,kc)) &
               -(uy(ic,jc,kc)+uy(ic,jc,km))* (uz(ic,jc,kc)+uz(ic,jm,kc)) )*sucaj                
          h33=( (uz(ic,jc,kp)+uz(ic,jc,kc))* (uz(ic,jc,kp)+uz(ic,jc,kc)) &
               -(uz(ic,jc,kc)+uz(ic,jc,km))* (uz(ic,jc,kc)+uz(ic,jc,km)) )*qdx3

          visa=quarter*(nut(ic,jc,km)+nut(im,jc,km)+nut(ic,jc,kc)+nut(im,jc,kc))
          visb=quarter*(nut(ic,jc,km)+nut(ip,jc,km)+nut(ic,jc,kc)+nut(ip,jc,kc))
          visc=quarter*(nut(ic,jc,km)+nut(ic,jm,km)+nut(ic,jm,kc)+nut(ic,jc,kc))
          visd=quarter*(nut(ic,jc,km)+nut(ic,jp,km)+nut(ic,jp,kc)+nut(ic,jc,kc))
          sgs1=rdx*rdz     *(visb*(ux(ip,jc,kc)-ux(ip,jc,km))-visa*(ux(ic,jc,kc)-ux(ic,jc,km)))
          sgs2=rdz*rdyp(jc)*(visd*(uy(ic,jp,kc)-uy(ic,jp,km))-visc*(uy(ic,jc,kc)-uy(ic,jc,km)))
          sgs3=rdz2*(nut(ic,jc,kc)*(uz(ic,jc,kp)-uz(ic,jc,kc))-nut(ic,jc,km)*(uz(ic,jc,kc)-uz(ic,jc,km)))

          d11q3= ap1c(ic)*visb*(uz(ip,jc,kc)-uz(ic,jc,kc)) -am1c(ic)*visa*(uz(ic,jc,kc)-uz(im,jc,kc)) 
          d22q3= ap2c(jc)*visd*(uz(ic,jp,kc)-uz(ic,jc,kc)) -am2c(jc)*visc*(uz(ic,jc,kc)-uz(ic,jm,kc))        
          d33q3= sgs3
          dcq13= d11q3+d33q3
                
          gradp3= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz 
#ifdef CFDDEM
          convEd3=-h31-h32-h33+sgs1+sgs2+sgs3+dcq13+gravity(3)+half*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#elif CFDLPT_TwoWay
          convEd3=-h31-h32-h33+sgs1+sgs2+sgs3+dcq13+gravity(3)+FpForce_z(ic,jc,kc)
#else
          convEd3=-h31-h32-h33+sgs1+sgs2+sgs3+dcq13+gravity(3)
#endif   
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3+ two*pmBetaT*d22q3
          HistZold(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO 
  end subroutine clcRhsZ_PIMP_LES

  !******************************************************************
  ! clcU1Hat_PIMP
  !******************************************************************    
  subroutine clcU1Hat_PIMP(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
 
    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+rTemp*OutFlowInfoY(4,ic,kc)
        enddo
      ENDDO
    endif    
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc) 
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc) 
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+one
      enddo
    enddo   
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) =  RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_PIMP

  !******************************************************************
  ! clcU1Hat_PIMP_0
  !******************************************************************    
  subroutine clcU1Hat_PIMP_0(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK)::mjc,cjc,pjc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridfj
    
    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*two +one
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) =  RhsX(ic,jc,kc)
        enddo
      enddo  
      call InversePTriFixedCoe(mjc,cjc,pjc, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_PIMP_0

  !******************************************************************
  ! clcU1Hat_LES
  !******************************************************************    
  subroutine clcU1Hat_PIMP_LES(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX  
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer:: ic,jc,kc,im,jm,jp
    real(RK)::visc,visd,rt1,rt2
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        jm=jc-1
        jp=jc+1
        do ic=y1start(1),y1end(1)
          im=ic-1
          visc= quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc))
          visd= quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jp,kc)+nut(im,jp,kc))
          rt1 = pmBetaT*visc*am2cForCN(jc)
          rt2 = pmBetaT*visd*ap2cForCN(jc)
          tridmj(ic,jc) = -rt1 
          tridcj(ic,jc) =  rt1+rt2+one
          tridpj(ic,jc) = -rt2
          tridfj(ic,jc) =  RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_PIMP_LES

  !******************************************************************
  ! clcU1Hat_PIMP_LES_0
  !******************************************************************    
  subroutine clcU1Hat_PIMP_LES_0(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX  
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer:: ic,jc,kc,im,jm,jp
    real(RK)::visc,visd,rt1,rt2
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        jm=jc-1
        jp=jc+1
        do ic=y1start(1),y1end(1)
          im=ic-1
          visc= quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jm,kc)+nut(im,jm,kc))
          visd= quarter*(nut(ic,jc,kc)+nut(im,jc,kc)+nut(ic,jp,kc)+nut(im,jp,kc))
          rt1 = pmBetaT*visc*am2cForCN(jc)
          rt2 = pmBetaT*visd*ap2cForCN(jc)
          tridmj(ic,jc) = -rt1 
          tridcj(ic,jc) =  rt1+rt2+one
          tridpj(ic,jc) = -rt2
          tridfj(ic,jc) =  RhsX(ic,jc,kc)
        enddo
      enddo  
      call InversePeriodicTridiagonal(tridmj, tridcj, tridpj, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_PIMP_LES_0
  
  !******************************************************************
  ! clcU2Hat_PIMP
  !******************************************************************  
  subroutine clcU2Hat_PIMP(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in):: duy_ym
    
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2p(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsY(ic,jc,kc)=RhsY(ic,jc,kc)+rTemp*OutFlowInfoY(5,ic,kc)
        enddo
      ENDDO
    endif
    do ic=y1start(1),y1end(1) 
      tridpj(ic,1)=zero
      tridcj(ic,1)=one
      tridmj(ic,1)=zero
    enddo
    do jc=2,nyc
      do ic=y1start(1),y1end(1) 
        tridpj(ic,jc) = -pmBeta*ap2p(jc)
        tridcj(ic,jc) = -pmBeta*ac2p(jc)+one
        tridmj(ic,jc) = -pmBeta*am2p(jc)
      enddo
    enddo
    DO kc=y1start(3),y1end(3) 
      do ic=y1start(1),y1end(1)
        tridfj(ic,1)=duy_ym(ic,kc)
      enddo
      do jc=2,nyc
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsY(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=1,nyc
        do ic=y1start(1),y1end(1) 
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo  
    ENDDO
  end subroutine clcU2Hat_PIMP 

  !******************************************************************
  ! clcU2Hat_PIMP_0
  !******************************************************************  
  subroutine clcU2Hat_PIMP_0(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in):: duy_ym
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mjc,cjc,pjc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridfj

    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*two +one
    DO kc=y1start(3),y1end(3) 
      do jc=1,nyc
        do ic=y1start(1),y1end(1) 
          tridfj(ic,jc) = RhsY(ic,jc,kc)
        enddo
       enddo
       call InversePTriFixedCoe(mjc,cjc,pjc,tridfj,y1size(1),nyc)
       do jc=1,nyc
         do ic=y1start(1),y1end(1) 
           uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
         enddo
       enddo  
     ENDDO
  end subroutine clcU2Hat_PIMP_0 

  !******************************************************************
  ! clcU2Hat_PIMP_LES
  !******************************************************************  
  subroutine clcU2Hat_PIMP_LES(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in):: duy_ym
    
    ! locals
    integer::ic,jc,kc,jm
    real(RK)::rt1,rt2
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj

    DO kc=y1start(3),y1end(3) 
      do ic=y1start(1),y1end(1) 
        tridpj(ic,1)=zero
        tridcj(ic,1)=one
        tridmj(ic,1)=zero
        tridfj(ic,1)=duy_ym(ic,kc)
      enddo
      do jc=2,nyc
        jm=jc-1
        do ic=y1start(1),y1end(1) 
          rt1= pmBetaT*nut(ic,jm,kc)*am2p(jc)
          rt2= pmBetaT*nut(ic,jc,kc)*ap2p(jc)
          tridmj(ic,jc) = -rt1
          tridcj(ic,jc) =  rt1+rt2+one
          tridpj(ic,jc) = -rt2
          tridfj(ic,jc) =  RhsY(ic,jc,kc)
        enddo
       enddo
       call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
       do jc=1,nyc
         do ic=y1start(1),y1end(1) 
           uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
         enddo
       enddo  
     ENDDO
  end subroutine clcU2Hat_PIMP_LES

  !******************************************************************
  ! clcU2Hat_PIMP_LES_0
  !******************************************************************  
  subroutine clcU2Hat_PIMP_LES_0(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in):: duy_ym
    
    ! locals
    integer::ic,jc,kc,jm
    real(RK)::rt1,rt2
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj

    DO kc=y1start(3),y1end(3)
      do jc=1,nyc
        jm=jc-1
        do ic=y1start(1),y1end(1) 
          rt1= pmBetaT*nut(ic,jm,kc)*rdy2
          rt2= pmBetaT*nut(ic,jc,kc)*rdy2
          tridmj(ic,jc) = -rt1
          tridcj(ic,jc) =  rt1+rt2+one
          tridpj(ic,jc) = -rt2
          tridfj(ic,jc) =  RhsY(ic,jc,kc)
        enddo
       enddo
       call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
       do jc=1,nyc
         do ic=y1start(1),y1end(1) 
           uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
         enddo
       enddo  
     ENDDO
  end subroutine clcU2Hat_PIMP_LES_0
  
  !******************************************************************
  ! clcU3Hat_PIMP
  !******************************************************************  
  subroutine clcU3Hat_PIMP(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz
  
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsZ(ic,jc,kc)=RhsZ(ic,jc,kc)+rTemp*OutFlowInfoY(6,ic,kc)
        enddo
      ENDDO
    endif
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        tridpj(ic,jc) = -pmBeta*ap2cForCN(jc)
        tridmj(ic,jc) = -pmBeta*am2cForCN(jc)
        tridcj(ic,jc) = -tridpj(ic,jc)-tridmj(ic,jc)+one
      enddo
    enddo  
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsZ(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_PIMP

  !******************************************************************
  ! clcU3Hat_PIMP_0
  !******************************************************************  
  subroutine clcU3Hat_PIMP_0(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz
  
    ! locals
    integer::ic,jc,kc
    real(RK):: mjc,cjc,pjc
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridfj

    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*two +one
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) = RhsZ(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mjc,cjc,pjc,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_PIMP_0 

  !******************************************************************
  ! clcU3Hat_LES
  !******************************************************************  
  subroutine clcU3Hat_PIMP_LES(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz
  
    ! locals
    integer::ic,jc,kc,jm,jp,km
    real(RK)::visc,visd,rt1,rt2
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    
    DO kc=y1start(3),y1end(3)
      km=kc-1
      do jc=y1start(2),y1end(2)
        jm=jc-1
        jp=jc+1
        do ic=y1start(1),y1end(1)
          visc=quarter*(nut(ic,jc,km)+nut(ic,jm,km)+nut(ic,jc,kc)+nut(ic,jm,kc))
          visd=quarter*(nut(ic,jc,km)+nut(ic,jp,km)+nut(ic,jc,kc)+nut(ic,jp,kc))
          rt1= pmBetaT*visc*am2cForCN(jc)
          rt2= pmBetaT*visd*ap2cForCN(jc)
          tridmj(ic,jc) = -rt1
          tridcj(ic,jc) =  rt1+rt2+one
          tridpj(ic,jc) = -rt2
          tridfj(ic,jc) =  RhsZ(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_PIMP_LES

  !******************************************************************
  ! clcU3Hat_PIMP_LES_0
  !******************************************************************  
  subroutine clcU3Hat_PIMP_LES_0(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz
  
    ! locals
    integer::ic,jc,kc,jm,jp,km
    real(RK)::visc,visd,rt1,rt2
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    
    DO kc=y1start(3),y1end(3)
      km=kc-1
      do jc=y1start(2),y1end(2)
        jm=jc-1
        jp=jc+1
        do ic=y1start(1),y1end(1)
          visc=quarter*(nut(ic,jc,km)+nut(ic,jm,km)+nut(ic,jc,kc)+nut(ic,jm,kc))
          visd=quarter*(nut(ic,jc,km)+nut(ic,jp,km)+nut(ic,jc,kc)+nut(ic,jp,kc))
          rt1= pmBetaT*visc*am2cForCN(jc)
          rt2= pmBetaT*visd*ap2cForCN(jc)
          tridmj(ic,jc) = -rt1
          tridcj(ic,jc) =  rt1+rt2+one
          tridpj(ic,jc) = -rt2
          tridfj(ic,jc) =  RhsZ(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_PIMP_LES_0

  !******************************************************************
  ! PressureUpdate_PIMP
  !******************************************************************
  subroutine PressureUpdate_PIMP(pressure, prphiHalo)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout):: pressure
    
    integer::ic,jc,kc,jp,jm
    real(RK)::pmBetap,pmBetac,pmBetam

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        jp=jc+1
        jm=jc-1
        pmBetap= -pmBeta*ap2Pr(jc)
        pmBetac= -pmBeta*ac2Pr(jc) + one
        pmBetam= -pmBeta*am2Pr(jc)
        do ic=y1start(1),y1end(1)
          pressure(ic,jc,kc)= pressure(ic,jc,kc)+ pmBetap*prphiHalo(ic,jp,kc)+ pmBetac*prphiHalo(ic,jc,kc)+ &
                                                  pmBetam*prphiHalo(ic,jm,kc)
        enddo
      enddo
    ENDDO
  end subroutine PressureUpdate_PIMP
  
  !******************************************************************
  ! PressureUpdate_PIMP_LES
  !******************************************************************
  subroutine PressureUpdate_PIMP_LES(pressure, prphiHalo)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout):: pressure
    
    !locals
    integer::ic,jc,kc,jp,jm
    real(RK)::visd,visc,metr1,metr2,pmBetap,pmBetac,pmBetam

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        jp=jc+1
        jm=jc-1
        metr1=rdyp(jc)*rdyc(jc)
        metr2=rdyp(jc)*rdyc(jp)
        do ic=y1start(1),y1end(1)
          visd= half*(nut(ic,jp,kc)+nut(ic,jc,kc))
          visc= half*(nut(ic,jc,kc)+nut(ic,jm,kc))
          pmBetap= -pmBetaT*metr2*visd
          pmBetam= -pmBetaT*metr1*visc
          pmBetac= -(pmBetap + pmBetam) +one
          pressure(ic,jc,kc)= pressure(ic,jc,kc)+ pmBetap*prphiHalo(ic,jp,kc)+ pmBetac*prphiHalo(ic,jc,kc)+ &
                                                  pmBetam*prphiHalo(ic,jm,kc)
        enddo
      enddo
    ENDDO
  end subroutine PressureUpdate_PIMP_LES
