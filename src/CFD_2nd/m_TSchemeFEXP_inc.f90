
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
    real(RK)::qdx1,qdx3,h11,h12,h13,sucaj,s1tot,s1tot1,s1totCoe,dpmdxns
    real(RK)::d11q1,d22q1,d33q1,dcq123,convEd1,gradp1,InterpY1,InterpY2,InterpY3,InterpY4
   
    s1tot=zero 
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        s1totCoe=dyp(jc)
        sucaj= half*rdyp(jc)
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
          d22q1= ap2c(jc)*ux(ic,jp,kc) +ac2c(jc)*ux(ic,jc,kc) +am2c(jc)*ux(ic,jm,kc)                
          d33q1= ap3c(kc)*ux(ic,jc,kp) +ac3c(kc)*ux(ic,jc,kc) +am3c(kc)*ux(ic,jc,km)
          dcq123= d11q1+d22q1+d33q1
             
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
          s1tot=s1tot+ux(ic,jc,kc)*s1totCoe
#ifdef CFDDEM
          convEd1= -h11-h12-h13 + xnu*dcq123 + gravity(1)+half*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
#elif CFDLPT_TwoWay 
          convEd1= -h11-h12-h13 + xnu*dcq123 + gravity(1)+FpForce_x(ic,jc,kc)
#else
          convEd1= -h11-h12-h13 + xnu*dcq123 + gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1+ pmTheta*HistXold(ic,jc,kc)- pmAlpha*gradp1
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
  end subroutine clcRhsX_FEXP

  !******************************************************************
  ! clcRhsX_FEXP_LES
  !******************************************************************    
  subroutine  clcRhsX_FEXP_LES(ux,uy,uz,RhsX,HistXold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsX
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistXold
    
    ! locals    
    integer::im,ic,ip,jc,jm,jp,km,kc,kp,ierror
    real(RK)::visa,visb,visc,visd,sgs1,sgs2,sgs3
    real(RK)::d11q1,d22q1,d33q1,dcq123,convEd1,gradp1
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
          im=ic-1; ip=ic+1
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
          dcq123= d11q1+d22q1+d33q1
                
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
          s1tot=s1tot+ux(ic,jc,kc)*s1totCoe         
#ifdef CFDDEM
          convEd1= -h11-h12-h13 +sgs1+sgs2+sgs3+dcq123+gravity(1)+half*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
#elif CFDLPT_TwoWay
          convEd1= -h11-h12-h13 +sgs1+sgs2+sgs3+dcq123+gravity(1)+FpForce_x(ic,jc,kc)
#else
          convEd1= -h11-h12-h13 +sgs1+sgs2+sgs3+dcq123+gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1+ pmTheta*HistXold(ic,jc,kc)- pmAlpha*gradp1
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
  end subroutine clcRhsX_FEXP_LES

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
    
    hdx1=half*rdx
    hdx3=half*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        sucac = rdyc(jc)
        qsucac= quarter*sucac
        InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1
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
  ! clcRhsY_FEXP_LES
  !******************************************************************    
  subroutine  clcRhsY_FEXP_LES(ux,uy,uz,RhsY,HistYold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsY
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistYold
   
    ! locals 
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::qdx1,qdx3,h21,h22,h23,sucac,qsucac
    real(RK)::d11q2,d22q2,d33q2,dcq123,convEd2,gradp2  
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
          dcq123= d11q2+d22q2+d33q2
                
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
#ifdef CFDDEM
          convEd2= -h21-h22-h23+ sgs1+sgs2+sgs3 +dcq123+ gravity(2)+half*(FpForce_y(ic,jc,kc)+FpForce_y(ic,jm,kc))
#elif CFDLPT_TwoWay         
          convEd2= -h21-h22-h23+ sgs1+sgs2+sgs3 +dcq123+ gravity(2)+FpForce_y(ic,jc,kc)
#else
          convEd2= -h21-h22-h23+ sgs1+sgs2+sgs3 +dcq123+ gravity(2)
#endif
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2
          HistYold(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO
  end subroutine clcRhsY_FEXP_LES

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
    
    qdx1=quarter*rdx
    qdx3=quarter*rdz
    DO kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;jp=jc+1
        sucaj =half*rdyp(jc)
        InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=one-InterpY3
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
          d33q3= (uz(ic,jc,kp)-two*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          dcq123= d11q3+d22q3+d33q3
                
          gradp3= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz                
#ifdef CFDDEM
          convEd3= -h31-h32-h33+xnu*dcq123+ gravity(3)+half*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
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
  ! clcRhsZ_FEXP_LES
  !****************************************************************** 
  subroutine  clcRhsZ_FEXP_LES(ux,uy,uz,RhsZ,HistZold,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::HistZold
   
    ! locals
    integer::im,ic,ip,jc,jm,jp,km,kc,kp
    real(RK)::qdx1,qdx3,h31,h32,h33,sucaj
    real(RK)::d11q3,d22q3,d33q3,dcq123,convEd3,gradp3
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
          dcq123= d11q3+d22q3+d33q3
                
          gradp3= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz 
#ifdef CFDDEM
          convEd3=-h31-h32-h33+sgs1+sgs2+sgs3+dcq123+gravity(3)+half*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#elif CFDLPT_TwoWay
          convEd3=-h31-h32-h33+sgs1+sgs2+sgs3+dcq123+gravity(3)+FpForce_z(ic,jc,kc)
#else
          convEd3=-h31-h32-h33+sgs1+sgs2+sgs3+dcq123+gravity(3)
#endif
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3
          HistZold(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO 
  end subroutine clcRhsZ_FEXP_LES

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
  ! clcU1Hat_FEXP_LES
  !******************************************************************    
  subroutine clcU1Hat_FEXP_LES(ux,RhsX)
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
  end subroutine clcU1Hat_FEXP_LES
  
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
  ! clcU2Hat_FEXP_LES
  !******************************************************************  
  subroutine clcU2Hat_FEXP_LES(uy,RhsY,duy_ym)
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
  end subroutine clcU2Hat_FEXP_LES
  
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
  ! clcU3Hat_FEXP_LES
  !******************************************************************  
  subroutine clcU3Hat_FEXP_LES(uz,RhsZ)
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
  end subroutine clcU3Hat_FEXP_LES  

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
  
  !******************************************************************
  ! PressureUpdate_FEXP_LES
  !******************************************************************
  subroutine PressureUpdate_FEXP_LES(pressure, prphiHalo)
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
  end subroutine PressureUpdate_FEXP_LES
