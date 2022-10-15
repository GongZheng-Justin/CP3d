
  ! This file is included in the module m_TScheme

  !******************************************************************
  ! clcRhsX_FIMP
  !******************************************************************    
  subroutine  clcRhsX_FIMP(ux,uy,uz,RhsX,HistXold,pressure)
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
      km=kc-1; kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1; jp=jc+1
        s1totCoe=dyp(jc)
        sucaj=half*rdyp(jc)
        if(yc(jc)<ybulk1 .or. yc(jc)>=ybulk2) s1totCoe=zero
        InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1  
        InterpY3= YinterpCoe(jp); InterpY4=one-InterpY3            
        do ic=y1start(1),y1end(1)
          im=ic-1; ip=ic+1

          h11=( (ux(ip,jc,kc)+ux(ic,jc,kc))* (ux(ip,jc,kc)+ux(ic,jc,kc))  &
               -(ux(ic,jc,kc)+ux(im,jc,kc))* (ux(ic,jc,kc)+ux(im,jc,kc)) )*qdx1 
          h12=( (uy(ic,jp,kc)+uy(im,jp,kc))* (InterpY3*ux(ic,jc,kc) +InterpY4*ux(ic,jp,kc))  &
               -(uy(ic,jc,kc)+uy(im,jc,kc))* (InterpY1*ux(ic,jm,kc) +InterpY2*ux(ic,jc,kc)) )*sucaj            
          h13=( (uz(ic,jc,kp)+uz(im,jc,kp))* (ux(ic,jc,kp)+ux(ic,jc,kc))  &
               -(uz(ic,jc,kc)+uz(im,jc,kc))* (ux(ic,jc,kc)+ux(ic,jc,km)) )*qdx3
                
          d11q1= (ux(ip,jc,kc)-two*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          d22q1= ap2c(jc)*ux(ic,jp,kc) + ac2c(jc)*ux(ic,jc,kc)+ am2c(jc)*ux(ic,jm,kc)                
          d33q1= ap3c(kc)*ux(ic,jc,kp) + ac3c(kc)*ux(ic,jc,kc)+ am3c(kc)*ux(ic,jc,km)
          dcq123= d11q1+d22q1+d33q1
             
          gradp1= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
          s1tot=s1tot+ux(ic,jc,kc)*s1totCoe
#ifdef CFDDEM
          convEd1= -h11-h12-h13 + gravity(1)+half*(FpForce_x(ic,jc,kc)+FpForce_x(im,jc,kc))
#elif CFDLPT_TwoWay
          convEd1= -h11-h12-h13 + gravity(1)+FpForce_x(ic,jc,kc)
#else
          convEd1= -h11-h12-h13 + gravity(1)
#endif
          RhsX(ic,jc,kc)=pmGamma*convEd1+ pmTheta*HistXold(ic,jc,kc)- pmAlpha*gradp1+ two*pmBeta*dcq123
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
  end subroutine clcRhsX_FIMP

  !******************************************************************
  ! clcRhsY_FIMP
  !******************************************************************    
  subroutine  clcRhsY_FIMP(ux,uy,uz,RhsY,HistYold,pressure)
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
          dcq123= d11q2+d22q2+d33q2
       
          gradp2= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
#ifdef CFDDEM
          convEd2= -h21-h22-h23+ gravity(2) +InterpY1*FpForce_y(ic,jm,kc)+InterpY2*FpForce_y(ic,jc,kc)
#elif CFDLPT_TwoWay
          convEd2= -h21-h22-h23+ gravity(2) +FpForce_y(ic,jc,kc)
#else
          convEd2= -h21-h22-h23+ gravity(2)
#endif
          RhsY(ic,jc,kc)=pmGamma*convEd2+ pmTheta*HistYold(ic,jc,kc)- pmAlpha*gradp2+ two*pmBeta*dcq123
          HistYold(ic,jc,kc)=convEd2   
        enddo
      enddo
    ENDDO
  end subroutine clcRhsY_FIMP

  !******************************************************************
  ! clcRhsZ_FIMP
  !****************************************************************** 
  subroutine  clcRhsZ_FIMP(ux,uy,uz,RhsZ,HistZold,pressure)
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
          convEd3= -h31-h32-h33+ gravity(3)+half*(FpForce_z(ic,jc,kc)+FpForce_z(ic,jc,km))
#elif CFDLPT_TwoWay
          convEd3= -h31-h32-h33+ gravity(3)+FpForce_z(ic,jc,kc)
#else
          convEd3= -h31-h32-h33+ gravity(3)
#endif
          RhsZ(ic,jc,kc)= pmGamma*convEd3+ pmTheta*HistZold(ic,jc,kc)- pmAlpha*gradp3+ two*pmBeta*dcq123
          HistZold(ic,jc,kc)=convEd3
        enddo
      enddo
    ENDDO 
  end subroutine clcRhsZ_FIMP

  !******************************************************************
  ! clcU1Hat_FIMP_000
  !******************************************************************    
  subroutine clcU1Hat_FIMP_000(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mjc,cjc,pjc,mkc,ckc,pkc
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y1_to_x1(RhsX,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsX) 

    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_y1_to_z1(RhsX,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsX) 

    ! compute dq2* sweeping in the y direction
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
  end subroutine clcU1Hat_FIMP_000

  !******************************************************************
  ! clcU1Hat_FIMP_010
  !******************************************************************    
  subroutine clcU1Hat_FIMP_010(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+rTemp*OutFlowInfoY(4,ic,kc)
        enddo
      ENDDO
    endif
        
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y1_to_x1(RhsX,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsX)

    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_y1_to_z1(RhsX,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,Rhsx)

    ! compute dq2* sweeping in the y direction
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
          tridfj(ic,jc) = RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_010
  
  !******************************************************************
  ! clcU1Hat_FIMP_011
  !******************************************************************    
  subroutine clcU1Hat_FIMP_011(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+rTemp*OutFlowInfoY(4,ic,kc)
        enddo
      ENDDO
    endif
        
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y1_to_x1(RhsX,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsX) 

    ! compute dq3* sweeping in the z direction
    do kc=1,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+one
      enddo
    enddo
    call transpose_y1_to_z1(RhsX,arrz1)
    Do jc=1,z1size(2)
      do kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsX)

    ! compute dq2* sweeping in the y direction
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
          tridfj(ic,jc) = RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_011

  !******************************************************************
  ! clcU1Hat_FIMP_111
  !******************************************************************    
  subroutine clcU1Hat_FIMP_111(ux,RhsX)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    
    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+rTemp*OutFlowInfoY(4,ic,kc)
        enddo
      ENDDO
    endif
    if(myProcNghBC(y_pencil,3)==BC_outFlow) then
      ic=nxc; rTemp=pmBeta*rdx2
      DO kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          RhsX(ic,jc,kc)=RhsX(ic,jc,kc)+rTemp*OutFlowInfoX(4,jc,kc)
        enddo
      ENDDO      
    endif
        
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    do jc=1,x1size(2)
      tridmi(jc,1) = zero
      tridci(jc,1) = one
      tridpi(jc,1) = zero
    enddo
    do ic=2,x1size(1)
      do jc=1,x1size(2)  
        tridmi(jc,ic) =  mic
        tridci(jc,ic) =  cic
        tridpi(jc,ic) =  pic
      enddo
    enddo
    call transpose_y1_to_x1(RhsX,arrx1)
    DO kc=1,x1size(3)
      do jc=1,x1size(2)
        tridfi(jc,1) = zero
        do ic=2,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmi,tridci,tridpi,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsX)

    ! compute dq3* sweeping in the z direction
    do kc=1,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+one
      enddo
    enddo
    call transpose_y1_to_z1(RhsX,arrz1)
    Do jc=1,z1size(2)
      do kc=1,z1size(3)  
        do ic=1,z1size(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsX)

    ! compute dq2* sweeping in the y direction
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
          tridfj(ic,jc) = RhsX(ic,jc,kc)
        enddo
      enddo  
      call InverseTridiagonal(tridmj, tridcj, tridpj, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= ux(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU1Hat_FIMP_111

  !******************************************************************
  ! clcU2Hat_FIMP_000
  !******************************************************************    
  subroutine clcU2Hat_FIMP_000(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in):: duy_ym    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mjc,cjc,pjc,mkc,ckc,pkc
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y1_to_x1(RhsY,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsY) 

    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_y1_to_z1(RhsY,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsY)

    ! compute dq2* sweeping in the y direction
    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*two +one
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) =  RhsY(ic,jc,kc)
        enddo
      enddo  
      call InversePTriFixedCoe(mjc,cjc,pjc, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uy(ic,jc,kc)= uy(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU2Hat_FIMP_000

  !******************************************************************
  ! clcU2Hat_FIMP_010
  !******************************************************************    
  subroutine clcU2Hat_FIMP_010(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in):: duy_ym    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2p(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsY(ic,jc,kc)=RhsY(ic,jc,kc)+rTemp*OutFlowInfoY(5,ic,kc)
        enddo
      ENDDO
    endif
        
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y1_to_x1(RhsY,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 
    call transpose_x1_to_y1(arrx1,RhsY)
    
    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_y1_to_z1(RhsY,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsY)
    
    ! compute dq2* sweeping in the y direction
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
  end subroutine clcU2Hat_FIMP_010

  !******************************************************************
  ! clcU2Hat_FIMP_011
  !******************************************************************    
  subroutine clcU2Hat_FIMP_011(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in):: duy_ym    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2p(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsY(ic,jc,kc)=RhsY(ic,jc,kc)+rTemp*OutFlowInfoY(5,ic,kc)
        enddo
      ENDDO
    endif
        
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y1_to_x1(RhsY,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsY)
    
    ! compute dq3* sweeping in the z direction
    do kc=1,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+one
      enddo
    enddo
    call transpose_y1_to_z1(RhsY,arrz1)
    Do jc=1,z1size(2)
      do kc=1,z1size(3)  
        do ic=1,z1size(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsY)

    ! compute dq2* sweeping in the y direction
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
  end subroutine clcU2Hat_FIMP_011

  !******************************************************************
  ! clcU2Hat_FIMP_111
  !******************************************************************    
  subroutine clcU2Hat_FIMP_111(uy,RhsY,duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uy
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in):: duy_ym    
    
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc
    real(RK),dimension(x1size(2),x1size(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2p(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsY(ic,jc,kc)=RhsY(ic,jc,kc)+rTemp*OutFlowInfoY(5,ic,kc)
        enddo
      ENDDO
    endif
    if(myProcNghBC(y_pencil,3)==BC_outFlow) then
      ic=nxc; rTemp=pmBeta*rdx2
      DO kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          RhsY(ic,jc,kc)=RhsY(ic,jc,kc)+rTemp*OutFlowInfoX(5,jc,kc)
        enddo
      ENDDO      
    endif
            
    ! compute dq1* sweeping in the x direction
    do ic=1,x1size(1)
      do jc=1,x1size(2)  
        tridpi(jc,ic)= -pmBeta*ap1cForCN(ic)
        tridmi(jc,ic)= -pmBeta*am1cForCN(ic)
        tridci(jc,ic)= -tridpi(jc,ic)-tridmi(jc,ic)+one
      enddo
    enddo
    call transpose_y1_to_x1(RhsY,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic)=  arrx1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmi,tridci,tridpi,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsY)
    
    ! compute dq3* sweeping in the z direction
    do kc=1,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) = -pmBeta*ap3cForCN(kc)
        tridmk(ic,kc) = -pmBeta*am3cForCN(kc)
        tridck(ic,kc) = -tridpk(ic,kc)-tridmk(ic,kc)+one
      enddo
    enddo
    call transpose_y1_to_z1(RhsY,arrz1)
    Do jc=1,z1size(2)
      do kc=1,z1size(3)  
        do ic=1,z1size(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk, tridck, tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsY)

    ! compute dq2* sweeping in the y direction
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
  end subroutine clcU2Hat_FIMP_111

  !******************************************************************
  ! clcU3Hat_FIMP_000
  !******************************************************************    
  subroutine clcU3Hat_FIMP_000(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mjc,cjc,pjc,mkc,ckc,pkc
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::RhsZ_temp

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y1_to_x1(RhsZ_temp,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO
    call transpose_x1_to_y1(arrx1,RhsZ_temp) 

    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_y1_to_z1(RhsZ_temp,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO 
    call transpose_z1_to_y1(arrz1,RhsZ_temp)
    
    ! compute dq2* sweeping in the y direction
    mjc= -pmBeta*rdy2
    pjc= -pmBeta*rdy2
    cjc=  pmBeta*rdy2*two +one
    DO kc=y1start(3),y1end(3)  
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          tridfj(ic,jc) =  RhsZ_temp(ic,jc,kc)
        enddo
      enddo  
      call InversePTriFixedCoe(mjc,cjc,pjc, tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_000

  !******************************************************************
  ! clcU3Hat_FIMP_010
  !******************************************************************    
  subroutine clcU3Hat_FIMP_010(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::RhsZ_temp

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsZ(ic,jc,kc)=RhsZ(ic,jc,kc)+rTemp*OutFlowInfoY(6,ic,kc)
        enddo
      ENDDO
    endif
    
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
      
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y1_to_x1(RhsZ_temp,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 
    call transpose_x1_to_y1(arrx1,RhsZ_temp)
    
    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    call transpose_y1_to_z1(RhsZ_temp,arrz1)
    Do jc=1,z1size(2)
      dO kc=1,z1size(3)  
        do ic=1,z1size(1) 
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mkc,ckc,pkc,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsZ_temp)

    ! compute dq2* sweeping in the y direction
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
          tridfj(ic,jc) = RhsZ_temp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_010

  !******************************************************************
  ! clcU3Hat_FIMP_011
  !******************************************************************    
  subroutine clcU3Hat_FIMP_011(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mic,cic,pic,mkc,ckc,pkc,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::RhsZ_temp

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsZ(ic,jc,kc)=RhsZ(ic,jc,kc)+rTemp*OutFlowInfoY(6,ic,kc)
        enddo
      ENDDO
    endif
    
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* sweeping in the x direction
    mic= -pmBeta*rdx2
    pic= -pmBeta*rdx2
    cic=  pmBeta*rdx2*two +one
    call transpose_y1_to_x1(RhsZ_temp,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InversePTriFixedCoe(mic,cic,pic,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 
    call transpose_x1_to_y1(arrx1,RhsZ_temp)
    
    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    do ic=1,z1size(1)
      tridpk(ic,1) =  zero
      tridck(ic,1) =  one
      tridmk(ic,1) =  zero
    enddo
    dO kc=2,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) =  pkc
        tridck(ic,kc) =  ckc
        tridmk(ic,kc) =  mkc
      enddo
    enddo
    call transpose_y1_to_z1(RhsZ_temp,arrz1)
    Do jc=1,z1size(2)
      do ic=1,z1size(1)
        tridfk(ic,1) =  zero
      enddo
      dO kc=2,z1size(3)  
        do ic=1,z1size(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk,tridck,tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsZ_temp)

    ! compute dq2* sweeping in the y direction
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
          tridfj(ic,jc) = RhsZ_temp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_011

  !******************************************************************
  ! clcU3Hat_FIMP_111 (unfinished !!!)
  !******************************************************************    
  subroutine clcU3Hat_FIMP_111(uz,RhsZ)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uz    
    
    ! locals
    integer::ic,jc,kc
    real(RK):: mkc,ckc,pkc,rTemp
    real(RK),dimension(x1size(2),x1size(1))::tridmi,tridci,tridpi,tridfi
    real(RK),dimension(z1size(1),z1size(3))::tridmk,tridck,tridpk,tridfk
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::RhsZ_temp

    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyc; rTemp=pmBeta*ap2cForCN(jc)
      DO kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          RhsZ(ic,jc,kc)=RhsZ(ic,jc,kc)+rTemp*OutFlowInfoY(6,ic,kc)
        enddo
      ENDDO
    endif
    if(myProcNghBC(y_pencil,3)==BC_outFlow) then
      ic=nxc; rTemp=pmBeta*rdx2
      DO kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          RhsZ(ic,jc,kc)=RhsZ(ic,jc,kc)+rTemp*OutFlowInfoX(6,jc,kc)
        enddo
      ENDDO      
    endif
        
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          RhsZ_temp(ic,jc,kc)= RhsZ(ic,jc,kc)
        enddo
      enddo
    ENDDO
    
    ! compute dq1* sweeping in the x direction
    do jc=1,x1size(2)
      do ic=1,x1size(1)
        tridpi(jc,ic)= -pmBeta*ap1cForCN(ic)
        tridmi(jc,ic)= -pmBeta*am1cForCN(ic)
        tridci(jc,ic)= -tridpi(jc,ic)-tridmi(jc,ic)+one
      enddo
    enddo
    call transpose_y1_to_x1(RhsZ_temp,arrx1)
    DO kc=1,x1size(3)  
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          tridfi(jc,ic) =  arrx1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmi,tridci,tridpi,tridfi,x1size(2),nxc)
      do jc=1,x1size(2)
        do ic=1,x1size(1)
          arrx1(ic,jc,kc)= tridfi(jc,ic)
        enddo
      enddo
    ENDDO 
    call transpose_x1_to_y1(arrx1,RhsZ_temp)
    
    ! compute dq3* sweeping in the z direction
    mkc= -pmBeta*rdz2
    pkc= -pmBeta*rdz2
    ckc=  pmBeta*rdz2*two +one
    do ic=1,z1size(1)
      tridpk(ic,1) =  zero
      tridck(ic,1) =  one
      tridmk(ic,1) =  zero
    enddo
    dO kc=2,z1size(3)  
      do ic=1,z1size(1)
        tridpk(ic,kc) =  pkc
        tridck(ic,kc) =  ckc
        tridmk(ic,kc) =  mkc
      enddo
    enddo
    call transpose_y1_to_z1(RhsZ_temp,arrz1)
    Do jc=1,z1size(2)
      do ic=1,z1size(1)
        tridfk(ic,1) =  zero
      enddo
      dO kc=2,z1size(3)  
        do ic=1,z1size(1)
          tridfk(ic,kc) =  arrz1(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmk,tridck,tridpk,tridfk,z1size(1),nzc)
      do kc=1,z1size(3)
        do ic=1,z1size(1)
          arrz1(ic,jc,kc)= tridfk(ic,kc)
        enddo
      enddo
    ENDDO
    call transpose_z1_to_y1(arrz1,RhsZ_temp)

    ! compute dq2* sweeping in the y direction
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
          tridfj(ic,jc) = RhsZ_temp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y1size(1),nyc)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uz(ic,jc,kc)= uz(ic,jc,kc)+tridfj(ic,jc)
        enddo
      enddo
    ENDDO
  end subroutine clcU3Hat_FIMP_111

  !******************************************************************
  ! PressureUpdate_FIMP
  !******************************************************************
  subroutine PressureUpdate_FIMP(pressure, prphiHalo)
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
  end subroutine PressureUpdate_FIMP
