#ifdef ScalarFlow
  !******************************************************************
  ! clcCosImagX
  !******************************************************************
  subroutine clcCosImagX(arrx1,arrx2,m1,m2)
    implicit none
    integer,intent(in)::m1,m2
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx1,arrx2

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp,2)

    normEgy= 1.0_RK/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,x1size(2)
      EgyX=0.0_RK
      do kc=1,x1size(3)
        EgyX(1,1)   = EgyX(1,1)   + arrx1(1,jc,kc)   *arrx2(1,jc,kc)
        EgyX(nxhp,1)= EgyX(nxhp,1)+ arrx1(nxhp,jc,kc)*arrx2(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic,1)= EgyX(ic,1)+ (arrx1(ic,jc,kc)*arrx2(ic,jc,kc)+arrx1(icCnter,jc,kc)*arrx2(icCnter,jc,kc))*2.0_RK
          EgyX(ic,2)= EgyX(ic,2)+ (arrx1(ic,jc,kc)*arrx2(icCnter,jc,kc)-arrx1(icCnter,jc,kc)*arrx2(ic,jc,kc))*2.0_RK
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m1)= EnergySpecX(ic,jc,m1) +EgyX(ic,1)*normEgy
        EnergySpecX(ic,jc,m2)= EnergySpecX(ic,jc,m2) +EgyX(ic,2)*normEgy
      enddo
    enddo
  end subroutine clcCosImagX
      
  !******************************************************************
  ! clcCosImagZ
  !******************************************************************
  subroutine clcCosImagZ(arrz1,arrz2,m1,m2)
    implicit none
    integer,intent(in)::m1,m2
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in)::arrz1,arrz2

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp,2)

    normEgy= 1.0_RK/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,z1size(2)
      EgyZ=0.0_RK
      do ic=1,z1size(1)
        EgyZ(1,1)   = EgyZ(1,1)   + arrz1(ic,jc,1)   *arrz2(ic,jc,1)
        EgyZ(nzhp,1)= EgyZ(nzhp,1)+ arrz1(ic,jc,nzhp)*arrz2(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc,1)= EgyZ(kc,1)+ (arrz1(ic,jc,kc)*arrz2(ic,jc,kc)+arrz1(ic,jc,kcCnter)*arrz2(ic,jc,kcCnter))*2.0_RK
          EgyZ(kc,2)= EgyZ(kc,2)+ (arrz1(ic,jc,kc)*arrz2(ic,jc,kcCnter)-arrz1(ic,jc,kcCnter)*arrz2(ic,jc,kc))*2.0_RK
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m1)= EnergySpecZ(kc,jc,m1) +EgyZ(kc,1)*normEgy
        EnergySpecZ(kc,jc,m2)= EnergySpecZ(kc,jc,m2) +EgyZ(kc,2)*normEgy
      enddo
    enddo
  end subroutine clcCosImagZ
#endif

#if defined(EnergySpectra_staggered_2nd) || defined(EnergySpectra_staggered_4th)
  !******************************************************************
  ! transform_uy
  !******************************************************************
  subroutine transform_uy(uyIn,uyOut)
    implicit none
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in) ::uyIn
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::uyOut
    
    ! locals
    integer::ic,jc,kc,jp
    
    if(FlowType==FT_CH) then
      DO kc=1,y1size(3)
        do jc=1,nyc
          jp=jc+1
          if(jc>nyc/2)jp=jc
          do ic=1,y1size(1)
            uyOut(ic,jc,kc)=uyIn(ic,jp,kc)
          enddo
        enddo
      ENDDO    
    else
      DO kc=1,y1size(3)
        do jc=1,nyc-1
          jp=jc+1
          do ic=1,y1size(1)
            uyOut(ic,jc,kc)=uyIn(ic,jp,kc)
          enddo
        enddo
      ENDDO
      uyOut(:,nyc,:)=0.0_RK    
    endif
  end subroutine transform_uy
#endif

  !******************************************************************
  ! clcEnergySpectraX
  !******************************************************************
  subroutine clcEnergySpectraX(arrx,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= 1.0_RK/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,x1size(2)
      EgyX=0.0_RK
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx(1,jc,kc)   *arrx(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx(nxhp,jc,kc)*arrx(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx(ic,jc,kc)*arrx(ic,jc,kc)+arrx(icCnter,jc,kc)*arrx(icCnter,jc,kc))*2.0_RK
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcEnergySpectraX
   
  !******************************************************************
  ! clcCospectraX
  !******************************************************************
  subroutine clcCospectraX(arrx1,arrx2,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx1,arrx2

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= 1.0_RK/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,x1size(2)
      EgyX=0.0_RK
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx1(1,jc,kc)   *arrx2(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx1(nxhp,jc,kc)*arrx2(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic)+ (arrx1(ic,jc,kc)*arrx2(ic,jc,kc)+arrx1(icCnter,jc,kc)*arrx2(icCnter,jc,kc))*2.0_RK
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcCospectraX
  
  !******************************************************************
  ! clcEnergySpectraZ
  !******************************************************************
  subroutine clcEnergySpectraZ(arrz,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in)::arrz

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= 1.0_RK/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,z1size(2)
      EgyZ=0.0_RK
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz(ic,jc,1)   *arrz(ic,jc,1)
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz(ic,jc,nzhp)*arrz(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc)+ (arrz(ic,jc,kc)*arrz(ic,jc,kc)+arrz(ic,jc,kcCnter)*arrz(ic,jc,kcCnter))*2.0_RK
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcEnergySpectraZ

  !******************************************************************
  ! clcCospectraZ
  !******************************************************************
  subroutine clcCospectraZ(arrz1,arrz2,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in)::arrz1,arrz2

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= 1.0_RK/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,z1size(2)
      EgyZ=0.0_RK
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz1(ic,jc,1)   *arrz2(ic,jc,1)
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz1(ic,jc,nzhp)*arrz2(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz1(ic,jc,kc)*arrz2(ic,jc,kc)+arrz1(ic,jc,kcCnter)*arrz2(ic,jc,kcCnter))*2.0_RK
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcCospectraZ
  
  !******************************************************************
  ! clcLTS_x
  !******************************************************************  
  subroutine clcLTS_x(arrx1,arrx2,ArrTemp,n1,n2,n3)
    implicit none
    integer,intent(in)::n1,n2,n3
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in) ::arrx1
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(out)::arrx2
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::ArrTemp

    ! locals
    integer::ic,jc,kc,jt,icCnter
    real(RK)::normEgy,EgyX(nxhp)
    real(RK),allocatable,dimension(:,:,:)::arrYplane
    
    if(FlowType==FT_HC) then
      allocate(arrYplane(y1size(1),y1size(3),1))
    else
      allocate(arrYplane(y1size(1),y1size(3),2))
    endif
    normEgy= 1.0_RK/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
 
    ! LCS -Real -Ref1
    call transpose_x1_to_y1(arrx1,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(1),kc)
        enddo
      enddo
    ELSE
      jt=nyc+1-jForLCS(1)
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=ArrTemp(ic,jt,kc)
        enddo
      enddo
    ENDIF
    
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_x1(ArrTemp,arrx2)
    do jc=1,x1size(2)
      EgyX=0.0_RK
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx2(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx2(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx2(ic,jc,kc)+arrx2(icCnter,jc,kc)) ! Note, there is NO "*2.0_RK" here !
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,n1)= EnergySpecX(ic,jc,n1) +EgyX(ic)*normEgy
      enddo
    enddo
    
    ! LCS -Imag -Ref1
    arrx2(1,:,:)=arrx1(1,:,:)
    do kc=1,x1size(3)
      do jc=1,x1size(2)
        do ic=2,x1size(1)
          icCnter=nxc+2-ic
          arrx2(ic,jc,kc)=arrx1(icCnter,jc,kc)
        enddo
      enddo
    enddo
    call transpose_x1_to_y1(arrx2,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_x1(ArrTemp,arrx2)
    do jc=1,x1size(2)
      EgyX=0.0_RK
      ! EgyX(1)   = EgyX(nxhp)= 0.0_RK
      do kc=1,x1size(3)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx2(ic,jc,kc)-arrx2(icCnter,jc,kc)) ! Note here !
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,n2)= EnergySpecX(ic,jc,n2) +EgyX(ic)*normEgy
      enddo
    enddo   

    ! LCS -Real -Ref2
    call transpose_x1_to_y1(arrx1,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(2),kc)
        enddo
      enddo
    ELSE
      jt=nyc+1-jForLCS(2)
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(2),kc)
          arrYplane(ic,kc,2)=ArrTemp(ic,jt,kc)
        enddo
      enddo
    ENDIF
    
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_x1(ArrTemp,arrx2)
    do jc=1,x1size(2)
      EgyX=0.0_RK
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx2(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx2(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx2(ic,jc,kc)+arrx2(icCnter,jc,kc)) ! Note, there is NO "*2.0_RK" here !
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,n3)= EnergySpecX(ic,jc,n3) +EgyX(ic)*normEgy
      enddo
    enddo

    deallocate(arrYplane)
  end subroutine clcLTS_x

  !******************************************************************
  ! clcLTS_z
  !******************************************************************
  subroutine clcLTS_z(arrz1,arrz2,ArrTemp,n1,n2,n3)
    implicit none
    integer,intent(in)::n1,n2,n3
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in) ::arrz1
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(out)::arrz2
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::ArrTemp

    ! locals
    integer::ic,jc,kc,jt,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)
    real(RK),allocatable,dimension(:,:,:)::arrYplane
    
    if(FlowType==FT_HC) then
      allocate(arrYplane(y1size(1),y1size(3),1))
    else
      allocate(arrYplane(y1size(1),y1size(3),2))        
    endif
    normEgy= 1.0_RK/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
        
    ! LCS -Real -Ref1
    call transpose_z1_to_y1(arrz1,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(1),kc)
        enddo
      enddo
    ELSE
      jt=nyc+1-jForLCS(1)
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=ArrTemp(ic,jt,kc)
        enddo
      enddo        
    ENDIF
    
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_z1(ArrTemp,arrz2)
    do jc=1,z1size(2)
      EgyZ=0.0_RK
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz2(ic,jc,1)   
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz2(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz2(ic,jc,kc)+arrz2(ic,jc,kcCnter)) ! Note, there is NO "*2.0_RK" here !
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,n1)= EnergySpecZ(kc,jc,n1) +EgyZ(kc)*normEgy
      enddo
    enddo
    
    ! LCS -Imag -Ref1
    arrz2(:,:,1)=arrz1(:,:,1)
    do kc=2,z1size(3)
      kcCnter=nzc+2-kc
      do jc=1,z1size(2)
        do ic=1,z1size(1)
          arrz2(ic,jc,kc)=arrz1(ic,jc,kcCnter)
        enddo
      enddo
    enddo
    call transpose_z1_to_y1(arrz2,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_z1(ArrTemp,arrz2)
    do jc=1,z1size(2)
      EgyZ=0.0_RK
      ! EgyZ(1)   = EgyZ(nzhp)= 0.0_RK
      do ic=1,z1size(1)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz2(ic,jc,kc)-arrz2(ic,jc,kcCnter)) ! Note here !
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,n2)= EnergySpecZ(kc,jc,n2) +EgyZ(kc)*normEgy
      enddo
    enddo

    ! LCS -Real -Ref2
    call transpose_z1_to_y1(arrz1,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(2),kc)
        enddo
      enddo
    ELSE
      jt=nyc+1-jForLCS(2)
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(2),kc)
          arrYplane(ic,kc,2)=ArrTemp(ic,jt,kc)
        enddo
      enddo        
    ENDIF
    
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_z1(ArrTemp,arrz2)
    do jc=1,z1size(2)
      EgyZ=0.0_RK
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz2(ic,jc,1)   
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz2(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz2(ic,jc,kc)+arrz2(ic,jc,kcCnter)) ! Note, there is NO "*2.0_RK" here !
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,n3)= EnergySpecZ(kc,jc,n3) +EgyZ(kc)*normEgy
      enddo
    enddo
    deallocate(arrYplane)
  end subroutine clcLTS_z

  !******************************************************************
  ! clcSpectraAndLCSR2D
  !******************************************************************
  subroutine clcSpectraAndLCSR2D(ArrIN,m1,m2)
    implicit none
    integer,intent(in)::m1,m2
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in)::ArrIn
    
    ! locals
    real(RK)::normEgy
    integer::ic,jc,kc,ict,jct,kct,jcs
    real(RK),allocatable,dimension(:,:,:)::arrYplane
    real(RK),allocatable,dimension(:,:,:)::arrx,arrx_xhzf
    real(RK),allocatable,dimension(:,:,:)::arry_xhzh,arry_xhzf
    real(RK),allocatable,dimension(:,:,:)::arrz_xhzf,arrz_xhzh,arrz_xhzh_LCS
    
    normEgy= 1.0_RK/(real(nxc,RK)*real(nzc,RK))
    normEgy= normEgy*normEgy
    IF(FlowType==FT_CH)normEgy=normEgy*0.5_RK
    allocate(arrx(nxc,x1size(2),x1size(3)))
    call transpose_y1_to_x1(ArrIN,arrx)
    call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    allocate(arrx_xhzf(nxhp,x1size(2),x1size(3)),   &
             arrz_xhzf(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    
    ! Real part
    arrx_xhzf=arrx(1:nxhp,:,:)
    call transpose_x1_to_z2(arrx_xhzf,arrz_xhzf,decomp_xhzf)
    call my_execute_FFTW_r2r_z(fft_plan_z,decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3),arrz_xhzf)
    arrx_xhzf(1,:,:)=0.0_RK
    arrx_xhzf(nxhp,:,:)=0.0_RK
    do ic=2,nxh
      ict=nxc+2-ic
      arrx_xhzf(ic,:,:)=arrx(ict,:,:)
    enddo
    deallocate(arrx)
    allocate(arrz_xhzh(decomp_xhzh%z2sz(1),decomp_xhzh%z2sz(2),decomp_xhzh%z2sz(3)))
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc )*arrz_xhzf(ic,jc,kc )+ &
                               arrz_xhzf(ic,jc,kct)*arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc) 
      enddo
    enddo
    
    ! LCS
    allocate(arry_xhzf(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(2),decomp_xhzf%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzf,arry_xhzf,decomp_xhzf)
    IF(FlowType==FT_HC) THEN
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),1))
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,decomp_xhzf%y2sz(2)
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
      enddo
    else
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),2))
      jct=nyc+1-jForLCS(1)
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=arry_xhzf(ic,jct,kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,nyc/2
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo        
      enddo
    ENDIF
    call transpose_y2_to_z2(arry_xhzf,arrz_xhzf,decomp_xhzf)
    deallocate(arry_xhzf,arrYplane)
    allocate(arrz_xhzh_LCS(decomp_xhzh%z2sz(1),decomp_xhzh%z2sz(2),decomp_xhzh%z2sz(3)))
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc )+ arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc)
      enddo
    enddo
    
    ! Imag part
    call transpose_x1_to_z2(arrx_xhzf,arrz_xhzf,decomp_xhzf)
    call my_execute_FFTW_r2r_z(fft_plan_z,decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3),arrz_xhzf)
    deallocate(arrx_xhzf)
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc )*arrz_xhzf(ic,jc,kc )+ &
                               arrz_xhzf(ic,jc,kct)*arrz_xhzf(ic,jc,kct)+ arrz_xhzh(ic,jc,kc)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
      enddo
    enddo
    
    ! LCS
    allocate(arry_xhzf(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(2),decomp_xhzf%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzf,arry_xhzf,decomp_xhzf)
    IF(FlowType==FT_HC) THEN
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),1))
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,decomp_xhzf%y2sz(2)
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
      enddo
    else
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),2))
      jct=nyc+1-jForLCS(1)
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=arry_xhzf(ic,jct,kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,nyc/2
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo        
      enddo
    ENDIF
    call transpose_y2_to_z2(arry_xhzf,arrz_xhzf,decomp_xhzf)
    deallocate(arry_xhzf,arrYplane)
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzh_LCS(ic,jc,kc)+arrz_xhzf(ic,jc,kc )+ arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh_LCS(ic,jc,kc)=arrz_xhzh_LCS(ic,jc,kc)+ arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh_LCS(ic,jc,kc)=arrz_xhzh_LCS(ic,jc,kc)+ arrz_xhzf(ic,jc,kc)
      enddo
    enddo
    
    deallocate(arrz_xhzf)
    allocate(arry_xhzh(decomp_xhzh%y2sz(1),decomp_xhzh%y2sz(2),decomp_xhzh%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzh,arry_xhzh,decomp_xhzh)
    deallocate(arrz_xhzh)

    IF(FlowType==FT_HC) THEN
      do jc=1,nyc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m1)=EnergySpec2D(ic,jc,kc,m1)+normEgy*arry_xhzh(ic,jc,kc)
          enddo
        enddo
      enddo
    ELSE 
      do jc=1,nyc/2
        jcs=nyc+1-jc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m1)=EnergySpec2D(ic,jc,kc,m1)+ &
                                      normEgy*(arry_xhzh(ic,jc,kc)+arry_xhzh(ic,jcs,kc))
          enddo
        enddo
      enddo
    ENDIF
    
    ! LCS
    call transpose_z2_to_y2(arrz_xhzh_LCS,arry_xhzh,decomp_xhzh)
    deallocate(arrz_xhzh_LCS)
    IF(FlowType==FT_HC) THEN
      do jc=1,nyc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m2)=EnergySpec2D(ic,jc,kc,m2)+normEgy*arry_xhzh(ic,jc,kc)
          enddo
        enddo
      enddo
    ELSE 
      do jc=1,nyc/2
        jcs=nyc+1-jc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m2)=EnergySpec2D(ic,jc,kc,m2)+ &
                                      normEgy*(arry_xhzh(ic,jc,kc)+arry_xhzh(ic,jcs,kc))
          enddo
        enddo
      enddo
    ENDIF
    deallocate(arry_xhzh)
  end subroutine clcSpectraAndLCSR2D

  !******************************************************************
  ! clcSpecAndLCSR2DFrom1D
  !******************************************************************
  subroutine clcSpecAndLCSR2DFrom1D(arrx,m1,m2)
    implicit none
    integer,intent(in)::m1,m2
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx
    
    ! locals
    real(RK)::normEgy
    integer::ic,jc,kc,ict,jct,kct,jcs
    real(RK),allocatable,dimension(:,:,:)::arrYplane
    real(RK),allocatable,dimension(:,:,:)::arrx_xhzf,arry_xhzh,arry_xhzf
    real(RK),allocatable,dimension(:,:,:)::arrz_xhzf,arrz_xhzh,arrz_xhzh_LCS
    
    normEgy= 1.0_RK/(real(nxc,RK)*real(nzc,RK))
    normEgy= normEgy*normEgy
    IF(FlowType==FT_CH)normEgy=normEgy*0.5_RK
    allocate(arrx_xhzf(nxhp,x1size(2),x1size(3)),   &
             arrz_xhzf(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    
    ! Real part
    arrx_xhzf=arrx(1:nxhp,:,:)
    call transpose_x1_to_z2(arrx_xhzf,arrz_xhzf,decomp_xhzf)
    call my_execute_FFTW_r2r_z(fft_plan_z,decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3),arrz_xhzf)
    arrx_xhzf(1,:,:)=0.0_RK
    arrx_xhzf(nxhp,:,:)=0.0_RK
    do ic=2,nxh
      ict=nxc+2-ic
      arrx_xhzf(ic,:,:)=arrx(ict,:,:)
    enddo
    allocate(arrz_xhzh(decomp_xhzh%z2sz(1),decomp_xhzh%z2sz(2),decomp_xhzh%z2sz(3)))
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc )*arrz_xhzf(ic,jc,kc )+ &
                               arrz_xhzf(ic,jc,kct)*arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc) 
      enddo
    enddo
    
    ! LCS
    allocate(arry_xhzf(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(2),decomp_xhzf%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzf,arry_xhzf,decomp_xhzf)
    IF(FlowType==FT_HC) THEN
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),1))
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,decomp_xhzf%y2sz(2)
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
      enddo
    else
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),2))
      jct=nyc+1-jForLCS(1)
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=arry_xhzf(ic,jct,kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,nyc/2
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo        
      enddo
    ENDIF
    call transpose_y2_to_z2(arry_xhzf,arrz_xhzf,decomp_xhzf)
    deallocate(arry_xhzf,arrYplane)
    allocate(arrz_xhzh_LCS(decomp_xhzh%z2sz(1),decomp_xhzh%z2sz(2),decomp_xhzh%z2sz(3)))
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc )+ arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc)
      enddo
    enddo
    
    ! Imag part
    call transpose_x1_to_z2(arrx_xhzf,arrz_xhzf,decomp_xhzf)
    call my_execute_FFTW_r2r_z(fft_plan_z,decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3),arrz_xhzf)
    deallocate(arrx_xhzf)
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc )*arrz_xhzf(ic,jc,kc )+ &
                               arrz_xhzf(ic,jc,kct)*arrz_xhzf(ic,jc,kct)+ arrz_xhzh(ic,jc,kc)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
      enddo
    enddo
    
    ! LCS
    allocate(arry_xhzf(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(2),decomp_xhzf%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzf,arry_xhzf,decomp_xhzf)
    IF(FlowType==FT_HC) THEN
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),1))
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,decomp_xhzf%y2sz(2)
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
      enddo
    else
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),2))
      jct=nyc+1-jForLCS(1)
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=arry_xhzf(ic,jct,kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,nyc/2
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo        
      enddo
    ENDIF
    call transpose_y2_to_z2(arry_xhzf,arrz_xhzf,decomp_xhzf)
    deallocate(arry_xhzf,arrYplane)
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzh_LCS(ic,jc,kc)+arrz_xhzf(ic,jc,kc )+ arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh_LCS(ic,jc,kc)=arrz_xhzh_LCS(ic,jc,kc)+ arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh_LCS(ic,jc,kc)=arrz_xhzh_LCS(ic,jc,kc)+ arrz_xhzf(ic,jc,kc)
      enddo
    enddo
    
    deallocate(arrz_xhzf)
    allocate(arry_xhzh(decomp_xhzh%y2sz(1),decomp_xhzh%y2sz(2),decomp_xhzh%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzh,arry_xhzh,decomp_xhzh)
    deallocate(arrz_xhzh)

    IF(FlowType==FT_HC) THEN
      do jc=1,nyc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m1)=EnergySpec2D(ic,jc,kc,m1)+normEgy*arry_xhzh(ic,jc,kc)
          enddo
        enddo
      enddo
    ELSE 
      do jc=1,nyc/2
        jcs=nyc+1-jc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m1)=EnergySpec2D(ic,jc,kc,m1)+ &
                                      normEgy*(arry_xhzh(ic,jc,kc)+arry_xhzh(ic,jcs,kc))
          enddo
        enddo
      enddo
    ENDIF
    
    ! LCS
    call transpose_z2_to_y2(arrz_xhzh_LCS,arry_xhzh,decomp_xhzh)
    deallocate(arrz_xhzh_LCS)
    IF(FlowType==FT_HC) THEN
      do jc=1,nyc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m2)=EnergySpec2D(ic,jc,kc,m2)+normEgy*arry_xhzh(ic,jc,kc)
          enddo
        enddo
      enddo
    ELSE 
      do jc=1,nyc/2
        jcs=nyc+1-jc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m2)=EnergySpec2D(ic,jc,kc,m2)+ &
                                      normEgy*(arry_xhzh(ic,jc,kc)+arry_xhzh(ic,jcs,kc))
          enddo
        enddo
      enddo
    ENDIF 
    deallocate(arry_xhzh)
  end subroutine clcSpecAndLCSR2DFrom1D
  
  !******************************************************************
  ! clcCosSpectra2D
  !******************************************************************
  subroutine clcCosSpectra2D(ArrIN1,ArrIN2,m,IsSymmetry)
    implicit none
    integer,intent(in)::m
    logical,intent(in)::IsSymmetry
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in)::ArrIn1,ArrIN2
    
    ! locals
    real(RK)::normEgy,rcoe
    integer::ic,jc,kc,ict,kct,jcs
    real(RK),dimension(:,:,:),allocatable::arry_xhzh
    real(RK),dimension(:,:,:),allocatable::arrx,arrx_xhzf1,arrx_xhzf2
    real(RK),dimension(:,:,:),allocatable::arrz_xhzf1,arrz_xhzf2,arrz_xhzh
    
    normEgy= 1.0_RK/(real(nxc,RK)*real(nzc,RK))
    normEgy= normEgy*normEgy
    IF(FlowType==FT_CH)normEgy=normEgy*0.5_RK
    allocate(arrx(nxc,x1size(2),x1size(3)))
    call transpose_y1_to_x1(ArrIN1,arrx)
    call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    allocate(arrx_xhzf1(nxhp,x1size(2),x1size(3)),   &
             arrz_xhzf1(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
                          
    ! Real part
    arrx_xhzf1=arrx(1:nxhp,:,:)
    call transpose_x1_to_z2(arrx_xhzf1,arrz_xhzf1,decomp_xhzf)
    call my_execute_FFTW_r2r_z(fft_plan_z,decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3),arrz_xhzf1)
    arrx_xhzf1(1,:,:)=0.0_RK
    arrx_xhzf1(nxhp,:,:)=0.0_RK
    do ic=2,nxh
      ict=nxc+2-ic
      arrx_xhzf1(ic,:,:)=arrx(ict,:,:)
    enddo
    call transpose_y1_to_x1(ArrIN2,arrx)
    call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    allocate(arrx_xhzf2(nxhp,x1size(2),x1size(3)),   &
             arrz_xhzf2(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    arrx_xhzf2=arrx(1:nxhp,:,:)
    call transpose_x1_to_z2(arrx_xhzf2,arrz_xhzf2,decomp_xhzf)
    call my_execute_FFTW_r2r_z(fft_plan_z,decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3),arrz_xhzf2)
    arrx_xhzf2(1,:,:)=0.0_RK
    arrx_xhzf2(nxhp,:,:)=0.0_RK
    do ic=2,nxh
      ict=nxc+2-ic
      arrx_xhzf2(ic,:,:)=arrx(ict,:,:)
    enddo
    deallocate(arrx)
    allocate(arrz_xhzh(decomp_xhzh%z2sz(1),decomp_xhzh%z2sz(2),decomp_xhzh%z2sz(3)))
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc )*arrz_xhzf2(ic,jc,kc )+ &
                               arrz_xhzf1(ic,jc,kct)*arrz_xhzf2(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc)*arrz_xhzf2(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc)*arrz_xhzf2(ic,jc,kc) 
      enddo
    enddo    
 
    ! Imag part
    call transpose_x1_to_z2(arrx_xhzf1,arrz_xhzf1,decomp_xhzf)
    call my_execute_FFTW_r2r_z(fft_plan_z,decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3),arrz_xhzf1)
    deallocate(arrx_xhzf1)
    call transpose_x1_to_z2(arrx_xhzf2,arrz_xhzf2,decomp_xhzf)
    call my_execute_FFTW_r2r_z(fft_plan_z,decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3),arrz_xhzf2)
    deallocate(arrx_xhzf2)   
    
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc )*arrz_xhzf2(ic,jc,kc )+ &
                               arrz_xhzf1(ic,jc,kct)*arrz_xhzf2(ic,jc,kct)+ arrz_xhzh(ic,jc,kc)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc)*arrz_xhzf2(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc)*arrz_xhzf2(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
      enddo
    enddo
    deallocate(arrz_xhzf1,arrz_xhzf2)    
    allocate(arry_xhzh(decomp_xhzh%y2sz(1),decomp_xhzh%y2sz(2),decomp_xhzh%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzh,arry_xhzh,decomp_xhzh)
    deallocate(arrz_xhzh)
            
    IF(FlowType==FT_HC) THEN
      do jc=1,nyc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m)=EnergySpec2D(ic,jc,kc,m)+normEgy*arry_xhzh(ic,jc,kc)
          enddo
        enddo
      enddo
    ELSE
      if(IsSymmetry) then
        rCoe= 1.0_RK      
      else
        rCoe=-1.0_RK
      endif 
      do jc=1,nyc/2
        jcs=nyc+1-jc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m)=EnergySpec2D(ic,jc,kc,m)+ &
                                     normEgy*(arry_xhzh(ic,jc,kc)+rCoe*arry_xhzh(ic,jcs,kc))
          enddo
        enddo
      enddo
    ENDIF
    deallocate(arry_xhzh)     
  end subroutine clcCosSpectra2D
