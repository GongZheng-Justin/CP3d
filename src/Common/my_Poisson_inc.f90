  !******************************************************************
  ! Initialize_ReduceMatrix
  !******************************************************************
  subroutine Initialize_ReduceMatrix(startIn,endIn,x_or_z)
    implicit none
    character(len=1),intent(in)::x_or_z
    integer,dimension(3),intent(in)::startIn,endIn

    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc
    logical,dimension(6)::initializeIn
    real(RK),dimension(:,:,:),allocatable::ArrM,ArrC,ArrP
    
    initializeIn=.false.
    if(x_or_z=='x') then
      initializeIn(4)=.true.
      call decomp_info_init(nxc,2*p_col,nzc,decomp_PPE,initialize=initializeIn)
    elseif(x_or_z=='z') then
      initializeIn(5)=.true.
      call decomp_info_init(nxc,2*p_row,nzc,decomp_PPE,initialize=initializeIn)
    endif
    allocate(ArrM(startIn(1):endIn(1),startIn(2):endIn(2),startIn(3):endIn(3)))
    allocate(ArrC(startIn(1):endIn(1),startIn(2):endIn(2),startIn(3):endIn(3)))
    allocate(ArrP(startIn(1):endIn(1),startIn(2):endIn(2),startIn(3):endIn(3)))
    do kc=startIn(3),endIn(3)
      do jc=startIn(2),endIn(2)
        do ic=startIn(1),endIn(1)
          ArrM(ic,jc,kc)=am2Pr(jc)
          ArrC(ic,jc,kc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          ArrP(ic,jc,kc)=ap2Pr(jc)
        enddo
      enddo
    enddo
    IF(nrank==0) THEN
      ArrM(1,1,1)=0.0_RK
      ArrC(1,1,1)=1.0_RK
      ArrP(1,1,1)=0.0_RK
    ENDIF
    do kc=startIn(3),endIn(3)
      do ic=startIn(1),endIn(1)
        jc=startIn(2)
        ArrM(ic,jc,kc)=ArrM(ic,jc,kc)/ArrC(ic,jc,kc)
        ArrP(ic,jc,kc)=ArrP(ic,jc,kc)/ArrC(ic,jc,kc)
        jc=startIn(2)+1
        ArrM(ic,jc,kc)=ArrM(ic,jc,kc)/ArrC(ic,jc,kc)
        ArrP(ic,jc,kc)=ArrP(ic,jc,kc)/ArrC(ic,jc,kc)
        do jc=startIn(2)+2,endIn(2)
          rTemp=1.0_RK/(ArrC(ic,jc,kc)-ArrM(ic,jc,kc)*ArrP(ic,jc-1,kc))
          ArrP(ic,jc,kc)= rTemp*ArrP(ic,jc,kc)
          ArrM(ic,jc,kc)=-rTemp*ArrM(ic,jc,kc)*ArrM(ic,jc-1,kc)
        enddo
        do jc=endIn(2)-2,startIn(2)+1,-1
          ArrM(ic,jc,kc)= ArrM(ic,jc,kc)-ArrP(ic,jc,kc)*ArrM(ic,jc+1,kc)
          ArrP(ic,jc,kc)=-ArrP(ic,jc,kc)*ArrP(ic,jc+1,kc)
        enddo
        jc=startIn(2)
        rTemp=1.0_RK/(1.0_RK-ArrM(ic,jc+1,kc)*ArrP(ic,jc,kc))
        ArrM(ic,jc,kc)= rTemp*ArrM(ic,jc,kc)
        ArrP(ic,jc,kc)=-rTemp*ArrP(ic,jc,kc)*ArrP(ic,jc+1,kc)        
      enddo
    enddo    
    deallocate(ArrC)
    allocate(ArrC(startIn(1):endIn(1),2,startIn(3):endIn(3)))
    
    do kc=startIn(3),endIn(3)
      do ic=startIn(1),endIn(1)
        ArrC(ic,1,kc)=ArrM(ic,startIn(2),kc)
        ArrC(ic,2,kc)=ArrM(ic,endIn(2),  kc)
      enddo
    enddo
    if(x_or_z=='x') then
      call transpose_x2_to_y2(ArrC,a_reduce,decomp_PPE)
    elseif(x_or_z=='z') then
      call transpose_z2_to_y2(ArrC,a_reduce,decomp_PPE)
    endif
    do kc=startIn(3),endIn(3)
      do ic=startIn(1),endIn(1)
        ArrC(ic,1,kc)=ArrP(ic,startIn(2),kc)
        ArrC(ic,2,kc)=ArrP(ic,endIn(2),  kc)
      enddo
    enddo
    if(x_or_z=='x') then
      call transpose_x2_to_y2(ArrC,c_reduce,decomp_PPE)
    elseif(x_or_z=='z') then
      call transpose_z2_to_y2(ArrC,c_reduce,decomp_PPE)
    endif  
    deallocate(ArrM,ArrC,ArrP)         
  end subroutine Initialize_ReduceMatrix

  !******************************************************************
  ! Destory_Poisson_FFT_Plan
  !******************************************************************
  subroutine Destory_Poisson_FFT_Plan()
    implicit none
    nullify(execute_FFTW_r2r_z)
    call dfftw_destroy_plan(fwd_plan_x,bwd_plan_x) 
    call dfftw_destroy_plan(fwd_plan_z,bwd_plan_z)
    if(allocated(decomp_PPE)) then
      call decomp_info_finalize(decomp_PPE)
      deallocate(decomp_PPE)
    endif
  end subroutine Destory_Poisson_FFT_Plan

  !******************************************************************
  ! clcPPE_x1
  !******************************************************************
  subroutine clcPPE_x1(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx
    real(RK),dimension(z2size(1),z2size(2),z2size(3))::arrz
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2),y2start(3):y2end(3))::prTemp
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2))::tridmj,tridcj,tridpj,tridfj
    
    call transpose_y1_to_x1(prsrc,arrx)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_z2(arrx,arrz)
    call execute_FFTW_r2r_z(fwd_plan_z,z2size(1),z2size(2),z2size(3),arrz)
    call transpose_z2_to_y2(arrz,prTemp)    

    do jc=y2start(2),y2end(2)
      do ic=y2start(1),y2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridpj(ic,jc)=ap2Pr(jc)
      enddo
    enddo
    IF(nrank==0) THEN
      do kc=2,y2end(3)
        do jc=1,nyc
          do ic=1,y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InverseTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo

      ! for nrank=0, y1start(1)=y1start(3)=1
      ic=1; jc=1; kc=1;
      tridpj(ic,jc)=0.0_RK
      tridcj(ic,jc)=1.0_RK
      tridmj(ic,jc)=0.0_RK
      tridfj(ic,jc)=0.0_RK
      do ic=2,y2end(1)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridfj(ic,jc)=prTemp(ic,jc,kc)
      enddo
      do jc=2,nyc
        do ic=1,y2end(1)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridfj(ic,jc)=prTemp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y2size(1),nyc)
      do jc=1,nyc
        do ic=1,y2end(1)
          prTemp(ic,jc,kc)= tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=y2start(3),y2end(3)
        do jc=y2start(2),y2end(2)
          do ic=y2start(1),y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InverseTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo       
    ENDIF
    call transpose_y2_to_z2(prTemp,arrz)
    call execute_FFTW_r2r_z(bwd_plan_z,z2size(1),z2size(2),z2size(3),arrz)
    call transpose_z2_to_x1(arrz,arrx)
    call my_execute_FFTW_r2r_x(bwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_y1(arrx,prsrc)
    
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_x1

  !******************************************************************
  ! clcPPE_x1_reduce
  !******************************************************************
  subroutine clcPPE_x1_reduce(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc,jp,jm
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx
    real(RK),dimension(z2start(1):z2end(1),2,z2start(3):z2end(3))::arrz2_reduce
    real(RK),dimension(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3))::arry2_reduce    
    real(RK),dimension(z2start(1):z2end(1),z2start(2):z2end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(z2start(1):z2end(1),z2start(2):z2end(2),z2start(3):z2end(3))::arrz1,arrz2,arrz3  

    call transpose_y1_to_x1(prsrc,arrx)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_z2(arrx,arrz2)
    call execute_FFTW_r2r_z(fwd_plan_z,z2size(1),z2size(2),z2size(3),arrz2)

    IF(nrank==0) THEN
      do kc=2,z2end(3)
        do jc=z2start(2),z2end(2)
          do ic=z2start(1),z2end(1)
            tridmj(ic,jc)=am2Pr(jc)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
            tridpj(ic,jc)=ap2Pr(jc)
            tridfj(ic,jc)=arrz2(ic,jc,kc)
          enddo
        enddo      
        do ic=z2start(1),z2end(1)
          jc=z2start(2)
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
          jc=z2start(2)+1
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        enddo
        do jc=z2start(2)+2,z2end(2)
          jm=jc-1
          do ic=z2start(1),z2end(1)
            rTemp=1.0_RK/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
            tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
            tridpj(ic,jc)= rTemp*tridpj(ic,jc)
            tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
          enddo
        enddo
        do jc=z2end(2)-2,z2start(2)+1,-1
          jp=jc+1
          do ic=z2start(1),z2end(1)
            tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
            tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
            tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
          enddo
        enddo
        do ic=z2start(1),z2end(1)
          jc=z2start(2); jp=jc+1
          rTemp=1.0_RK/(1.0_RK-tridmj(ic,jp)*tridpj(ic,jc))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
          tridmj(ic,jc)= rTemp*tridmj(ic,jc)
          tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
          arrz2_reduce(ic,1,kc)=tridfj(ic,z2start(2))
          arrz2_reduce(ic,2,kc)=tridfj(ic,z2end(2))      
        enddo
        do jc=z2start(2),z2end(2)
          do ic=z2start(1),z2end(1)
            arrz1(ic,jc,kc)=tridmj(ic,jc)
            arrz3(ic,jc,kc)=tridpj(ic,jc)
            arrz2(ic,jc,kc)=tridfj(ic,jc)
          enddo
        enddo
      enddo
      
      ! for nrank=0, z2start(1)=z2start(2)=1
      ic=1; jc=1; kc=1 
      tridmj(ic,jc)=0.0_RK
      tridcj(ic,jc)=1.0_RK
      tridpj(ic,jc)=0.0_RK
      tridfj(ic,jc)=0.0_RK
      do ic=2,z2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridpj(ic,jc)=ap2Pr(jc)
        tridfj(ic,jc)=arrz2(ic,jc,kc)
      enddo
      do jc=2,z2end(2)
        do ic=z2start(1),z2end(1)
          tridmj(ic,jc)=am2Pr(jc)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridpj(ic,jc)=ap2Pr(jc)
          tridfj(ic,jc)=arrz2(ic,jc,kc)
        enddo
      enddo            
      do ic=z2start(1),z2end(1)
        jc=z2start(2)
        tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
        tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
        tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        jc=z2start(2)+1
        tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
        tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
        tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
      enddo
      do jc=z2start(2)+2,z2end(2)
        jm=jc-1
        do ic=z2start(1),z2end(1)
          rTemp=1.0_RK/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
          tridpj(ic,jc)= rTemp*tridpj(ic,jc)
          tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
        enddo
      enddo
      do jc=z2end(2)-2,z2start(2)+1,-1
        jp=jc+1
        do ic=z2start(1),z2end(1)
          tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
          tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
          tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
        enddo
      enddo
      do ic=z2start(1),z2end(1)
        jc=z2start(2); jp=jc+1
        rTemp=1.0_RK/(1.0_RK-tridmj(ic,jp)*tridpj(ic,jc))
        tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
        tridmj(ic,jc)= rTemp*tridmj(ic,jc)
        tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
        arrz2_reduce(ic,1,kc)=tridfj(ic,z2start(2))
        arrz2_reduce(ic,2,kc)=tridfj(ic,z2end(2))      
      enddo
      do jc=z2start(2),z2end(2)
        do ic=z2start(1),z2end(1)
          arrz1(ic,jc,kc)=tridmj(ic,jc)
          arrz3(ic,jc,kc)=tridpj(ic,jc)
          arrz2(ic,jc,kc)=tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=z2start(3),z2end(3)
        do jc=z2start(2),z2end(2)
          do ic=z2start(1),z2end(1)
            tridmj(ic,jc)=am2Pr(jc)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
            tridpj(ic,jc)=ap2Pr(jc)
            tridfj(ic,jc)=arrz2(ic,jc,kc)
          enddo
        enddo      
        do ic=z2start(1),z2end(1)
          jc=z2start(2)
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
          jc=z2start(2)+1
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        enddo
        do jc=z2start(2)+2,z2end(2)
          jm=jc-1
          do ic=z2start(1),z2end(1)
            rTemp=1.0_RK/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
            tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
            tridpj(ic,jc)= rTemp*tridpj(ic,jc)
            tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
          enddo
        enddo
        do jc=z2end(2)-2,z2start(2)+1,-1
          jp=jc+1
          do ic=z2start(1),z2end(1)
            tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
            tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
            tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
          enddo
        enddo
        do ic=z2start(1),z2end(1)
          jc=z2start(2); jp=jc+1
          rTemp=1.0_RK/(1.0_RK-tridmj(ic,jp)*tridpj(ic,jc))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
          tridmj(ic,jc)= rTemp*tridmj(ic,jc)
          tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
          arrz2_reduce(ic,1,kc)=tridfj(ic,z2start(2))
          arrz2_reduce(ic,2,kc)=tridfj(ic,z2end(2))      
        enddo
        do jc=z2start(2),z2end(2)
          do ic=z2start(1),z2end(1)
            arrz1(ic,jc,kc)=tridmj(ic,jc)
            arrz3(ic,jc,kc)=tridpj(ic,jc)
            arrz2(ic,jc,kc)=tridfj(ic,jc)
          enddo
        enddo
      enddo
    ENDIF
    call transpose_z2_to_y2(arrz2_reduce,arry2_reduce,decomp_PPE)
    do kc=y2start(3),y2end(3)
      call InverseTridiagonal_reduce(a_reduce(:,:,kc),c_reduce(:,:,kc),arry2_reduce(:,:,kc),y2size(1),2*p_row)
    enddo
    call transpose_y2_to_z2(arry2_reduce,arrz2_reduce,decomp_PPE)
    do kc=z2start(3),z2end(3)
      do ic=z2start(1),z2end(1)
        arrz2(ic,z2start(2),kc)=arrz2_reduce(ic,1,kc)
        arrz2(ic,z2end(2),  kc)=arrz2_reduce(ic,2,kc)
      enddo
      do jc=z2start(2)+1,z2end(2)-1
        do ic=z2start(1),z2end(1)    
          arrz2(ic,jc,kc)=arrz2(ic,jc,kc)-arrz1(ic,jc,kc)*arrz2_reduce(ic,1,kc)-arrz3(ic,jc,kc)*arrz2_reduce(ic,2,kc)
        enddo
      enddo
    enddo

    call execute_FFTW_r2r_z(bwd_plan_z,z2size(1),z2size(2),z2size(3),arrz2)
    call transpose_z2_to_x1(arrz2,arrx)
    call my_execute_FFTW_r2r_x(bwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_y1(arrx,prsrc)
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_x1_reduce

  !******************************************************************
  ! clcPPE_z1
  !******************************************************************
  subroutine clcPPE_z1(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(x2size(1),x2size(2),x2size(3))::arrx
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2),y2start(3):y2end(3))::prTemp
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2))::tridmj,tridcj,tridpj,tridfj

    call transpose_y1_to_z1(prsrc,arrz)
    call execute_FFTW_r2r_z(fwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_x2(arrz,arrx)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x2size(1),x2size(2),x2size(3),arrx)
    call transpose_x2_to_y2(arrx,prTemp)
    do jc=y2start(2),y2end(2)
      do ic=y2start(1),y2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridpj(ic,jc)=ap2Pr(jc)
      enddo
    enddo
    IF(nrank==0) THEN
      do kc=2,y2end(3)
        do jc=1,nyc
          do ic=1,y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InverseTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo

      ! for nrank=0, y1start(1)=y1start(3)=1
      ic=1; jc=1; kc=1;
      tridpj(ic,jc)=0.0_RK
      tridcj(ic,jc)=1.0_RK
      tridmj(ic,jc)=0.0_RK
      tridfj(ic,jc)=0.0_RK
      do ic=2,y2end(1)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridfj(ic,jc)=prTemp(ic,jc,kc)
      enddo
      do jc=2,nyc
        do ic=1,y2end(1)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridfj(ic,jc)=prTemp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y2size(1),nyc)
      do jc=1,nyc
        do ic=1,y2end(1)
          prTemp(ic,jc,kc)= tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=y2start(3),y2end(3)
        do jc=y2start(2),y2end(2)
          do ic=y2start(1),y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InverseTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo       
    ENDIF
    call transpose_y2_to_x2(prTemp,arrx)
    call my_execute_FFTW_r2r_x(bwd_plan_x,x2size(1),x2size(2),x2size(3),arrx)
    call transpose_x2_to_z1(arrx,arrz)
    call execute_FFTW_r2r_z(bwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_y1(arrz,prsrc)
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_z1

  !******************************************************************
  ! clcPPE_z1_reduce
  !******************************************************************
  subroutine clcPPE_z1_reduce(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc,jp,jm
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz
    real(RK),dimension(x2start(1):x2end(1),2,x2start(3):x2end(3))::arrx2_reduce
    real(RK),dimension(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3))::arry2_reduce    
    real(RK),dimension(x2start(1):x2end(1),x2start(2):x2end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x2start(1):x2end(1),x2start(2):x2end(2),x2start(3):x2end(3))::arrx1,arrx2,arrx3

    call transpose_y1_to_z1(prsrc,arrz)  
    call execute_FFTW_r2r_z(fwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_x2(arrz,arrx2)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x2size(1),x2size(2),x2size(3),arrx2)
    IF(nrank==0) THEN      
      do kc=2,x2end(3)
        do jc=x2start(2),x2end(2)
          do ic=x2start(1),x2end(1)
            tridmj(ic,jc)=am2Pr(jc)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
            tridpj(ic,jc)=ap2Pr(jc)
            tridfj(ic,jc)=arrx2(ic,jc,kc)
          enddo
        enddo
        do ic=x2start(1),x2end(1)
          jc=x2start(2)
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
          jc=x2start(2)+1
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        enddo
        do jc=x2start(2)+2,x2end(2)
          jm=jc-1
          do ic=x2start(1),x2end(1)
            rTemp=1.0_RK/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
            tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
            tridpj(ic,jc)= rTemp*tridpj(ic,jc)
            tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
          enddo
        enddo
        do jc=x2end(2)-2,x2start(2)+1,-1
          jp=jc+1
          do ic=x2start(1),x2end(1)
            tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
            tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
            tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
          enddo
        enddo
        do ic=x2start(1),x2end(1)
          jc=x2start(2); jp=jc+1
          rTemp=1.0_RK/(1.0_RK-tridmj(ic,jp)*tridpj(ic,jc))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
          tridmj(ic,jc)= rTemp*tridmj(ic,jc)
          tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
          arrx2_reduce(ic,1,kc)=tridfj(ic,x2start(2))
          arrx2_reduce(ic,2,kc)=tridfj(ic,x2end(2))       
        enddo
        do jc=x2start(2),x2end(2)
          do ic=x2start(1),x2end(1)
            arrx1(ic,jc,kc)=tridmj(ic,jc)
            arrx3(ic,jc,kc)=tridpj(ic,jc)
            arrx2(ic,jc,kc)=tridfj(ic,jc)
          enddo
        enddo
      enddo

      ! for nrank=0, x2start(2)=x2start(3)=1
      ic=1; jc=1; kc=1
      tridmj(ic,jc)=0.0_RK
      tridcj(ic,jc)=1.0_RK
      tridpj(ic,jc)=0.0_RK
      tridfj(ic,jc)=0.0_RK
      do ic=2,x2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridpj(ic,jc)=ap2Pr(jc)
        tridfj(ic,jc)=arrx2(ic,jc,kc)
      enddo
      do jc=2,x2end(2)
        do ic=1,x2end(1)
          tridmj(ic,jc)=am2Pr(jc)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridpj(ic,jc)=ap2Pr(jc)
          tridfj(ic,jc)=arrx2(ic,jc,kc)
        enddo
      enddo
      do ic=x2start(1),x2end(1)
        jc=x2start(2)
        tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
        tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
        tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        jc=x2start(2)+1
        tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
        tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
        tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
      enddo
      do jc=x2start(2)+2,x2end(2)
        jm=jc-1
        do ic=x2start(1),x2end(1)
          rTemp=1.0_RK/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
          tridpj(ic,jc)= rTemp*tridpj(ic,jc)
          tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
        enddo
      enddo
      do jc=x2end(2)-2,x2start(2)+1,-1
        jp=jc+1
        do ic=x2start(1),x2end(1)
          tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
          tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
          tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
        enddo
      enddo
      do ic=x2start(1),x2end(1)
        jc=x2start(2); jp=jc+1
        rTemp=1.0_RK/(1.0_RK-tridmj(ic,jp)*tridpj(ic,jc))
        tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
        tridmj(ic,jc)= rTemp*tridmj(ic,jc)
        tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
        arrx2_reduce(ic,1,kc)=tridfj(ic,x2start(2))
        arrx2_reduce(ic,2,kc)=tridfj(ic,x2end(2))       
      enddo
      do jc=x2start(2),x2end(2)
        do ic=x2start(1),x2end(1)
          arrx1(ic,jc,kc)=tridmj(ic,jc)
          arrx3(ic,jc,kc)=tridpj(ic,jc)
          arrx2(ic,jc,kc)=tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=x2start(3),x2end(3)
        do jc=x2start(2),x2end(2)
          do ic=x2start(1),x2end(1)
            tridmj(ic,jc)=am2Pr(jc)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
            tridpj(ic,jc)=ap2Pr(jc)
            tridfj(ic,jc)=arrx2(ic,jc,kc)
          enddo
        enddo
        do ic=x2start(1),x2end(1)
          jc=x2start(2)
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
          jc=x2start(2)+1
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        enddo
        do jc=x2start(2)+2,x2end(2)
          jm=jc-1
          do ic=x2start(1),x2end(1)
            rTemp=1.0_RK/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
            tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
            tridpj(ic,jc)= rTemp*tridpj(ic,jc)
            tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
          enddo
        enddo
        do jc=x2end(2)-2,x2start(2)+1,-1
          jp=jc+1
          do ic=x2start(1),x2end(1)
            tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
            tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
            tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
          enddo
        enddo
        do ic=x2start(1),x2end(1)
          jc=x2start(2); jp=jc+1
          rTemp=1.0_RK/(1.0_RK-tridmj(ic,jp)*tridpj(ic,jc))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
          tridmj(ic,jc)= rTemp*tridmj(ic,jc)
          tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
          arrx2_reduce(ic,1,kc)=tridfj(ic,x2start(2))
          arrx2_reduce(ic,2,kc)=tridfj(ic,x2end(2))       
        enddo
        do jc=x2start(2),x2end(2)
          do ic=x2start(1),x2end(1)
            arrx1(ic,jc,kc)=tridmj(ic,jc)
            arrx3(ic,jc,kc)=tridpj(ic,jc)
            arrx2(ic,jc,kc)=tridfj(ic,jc)
          enddo
        enddo
      enddo
    ENDIF
    call transpose_x2_to_y2(arrx2_reduce,arry2_reduce,decomp_PPE)
    do kc=y2start(3),y2end(3)
      call InverseTridiagonal_reduce(a_reduce(:,:,kc),c_reduce(:,:,kc),arry2_reduce(:,:,kc),y2size(1),2*p_col)
    enddo
    call transpose_y2_to_x2(arry2_reduce,arrx2_reduce,decomp_PPE)
    do kc=x2start(3),x2end(3)
      do ic=x2start(1),x2end(1)
        arrx2(ic,x2start(2),kc)=arrx2_reduce(ic,1,kc)
        arrx2(ic,x2end(2),  kc)=arrx2_reduce(ic,2,kc)
      enddo
      do jc=x2start(2)+1,x2end(2)-1
        do ic=x2start(1),x2end(1)    
          arrx2(ic,jc,kc)=arrx2(ic,jc,kc)-arrx1(ic,jc,kc)*arrx2_reduce(ic,1,kc)-arrx3(ic,jc,kc)*arrx2_reduce(ic,2,kc)
        enddo
      enddo
    enddo

    call my_execute_FFTW_r2r_x(bwd_plan_x,x2size(1),x2size(2),x2size(3),arrx2)
    call transpose_x2_to_z1(arrx2,arrz)
    call execute_FFTW_r2r_z(bwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_y1(arrz,prsrc)
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_z1_reduce

  !******************************************************************
  ! InverseTridiagonal_reduce
  !******************************************************************  
  subroutine InverseTridiagonal_reduce(aj,cj,fj,m,n)
    implicit none
    integer,intent(in)::m,n
    real(RK),dimension(m,n),intent(in):: aj,cj
    real(RK),dimension(m,n),intent(inout)::fj
    
    ! locals
    integer:: i,j
    real(RK),dimension(m):: vecm
    real(RK),dimension(m,n)::arrmn
      
    do i=1,m
      vecm(i)=1.0_RK
    enddo
    do j=2,n
      do i=1,m
        arrmn(i,j)=cj(i,j-1)/vecm(i)
        vecm(i)=1.0_RK-aj(i,j)*arrmn(i,j)
        fj(i,j)=(fj(i,j)-aj(i,j)*fj(i,j-1))/vecm(i)
      enddo
    enddo
    do j=n-1,1,-1
      do i=1,m
        fj(i,j)=fj(i,j)-arrmn(i,j+1)*fj(i,j+1)
       enddo
    enddo    
  end subroutine InverseTridiagonal_reduce
  
#ifdef my_Poisson_inc_add_Periodic_2d
  !******************************************************************
  ! clcPPE_x1_periodic
  !******************************************************************
  subroutine clcPPE_x1_periodic(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx
    real(RK),dimension(z2size(1),z2size(2),z2size(3))::arrz
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2),y2start(3):y2end(3))::prTemp
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2))::tridmj,tridcj,tridpj,tridfj

    call transpose_y1_to_x1(prsrc,arrx)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_z2(arrx,arrz)
    call execute_FFTW_r2r_z(fwd_plan_z,z2size(1),z2size(2),z2size(3),arrz)
    call transpose_z2_to_y2(arrz,prTemp)    

    do jc=y2start(2),y2end(2)
      do ic=y2start(1),y2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridpj(ic,jc)=ap2Pr(jc)
      enddo
    enddo
    IF(nrank==0) THEN
      do kc=2,y2end(3)
        do jc=1,nyc
          do ic=1,y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo

      ! for nrank=0, y1start(1)=y1start(3)=1
      ic=1; jc=1; kc=1;
      tridpj(ic,jc)=0.0_RK
      tridcj(ic,jc)=1.0_RK
      tridmj(ic,jc)=0.0_RK
      tridfj(ic,jc)=0.0_RK
      do ic=2,y2end(1)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridfj(ic,jc)=prTemp(ic,jc,kc)
      enddo
      do jc=2,nyc
        do ic=1,y2end(1)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridfj(ic,jc)=prTemp(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,tridfj,y2size(1),nyc)
      do jc=1,nyc
        do ic=1,y2end(1)
          prTemp(ic,jc,kc)= tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=y2start(3),y2end(3)
        do jc=y2start(2),y2end(2)
          do ic=y2start(1),y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo       
    ENDIF
    call transpose_y2_to_z2(prTemp,arrz)
    call execute_FFTW_r2r_z(bwd_plan_z,z2size(1),z2size(2),z2size(3),arrz)
    call transpose_z2_to_x1(arrz,arrx)
    call my_execute_FFTW_r2r_x(bwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_y1(arrx,prsrc)
    
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_x1_periodic

  !******************************************************************
  ! clcPPE_z1_periodic
  !******************************************************************
  subroutine clcPPE_z1_periodic(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(x2size(1),x2size(2),x2size(3))::arrx
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2),y2start(3):y2end(3))::prTemp
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2))::tridmj,tridcj,tridpj,tridfj

    call transpose_y1_to_z1(prsrc,arrz)
    call execute_FFTW_r2r_z(fwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_x2(arrz,arrx)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x2size(1),x2size(2),x2size(3),arrx)
    call transpose_x2_to_y2(arrx,prTemp)
    do jc=y2start(2),y2end(2)
      do ic=y2start(1),y2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridpj(ic,jc)=ap2Pr(jc)
      enddo
    enddo
    IF(nrank==0) THEN
      do kc=2,y2end(3)
        do jc=1,nyc
          do ic=1,y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo

      ! for nrank=0, y1start(1)=y1start(3)=1
      ic=1; jc=1; kc=1;
      tridpj(ic,jc)=0.0_RK
      tridcj(ic,jc)=1.0_RK
      tridmj(ic,jc)=0.0_RK
      tridfj(ic,jc)=0.0_RK
      do ic=2,y2end(1)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridfj(ic,jc)=prTemp(ic,jc,kc)
      enddo
      do jc=2,nyc
        do ic=1,y2end(1)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridfj(ic,jc)=prTemp(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,tridfj,y2size(1),nyc)
      do jc=1,nyc
        do ic=1,y2end(1)
          prTemp(ic,jc,kc)= tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=y2start(3),y2end(3)
        do jc=y2start(2),y2end(2)
          do ic=y2start(1),y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo       
    ENDIF
    call transpose_y2_to_x2(prTemp,arrx)
    call my_execute_FFTW_r2r_x(bwd_plan_x,x2size(1),x2size(2),x2size(3),arrx)
    call transpose_x2_to_z1(arrx,arrz)
    call execute_FFTW_r2r_z(bwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_y1(arrz,prsrc)
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_z1_periodic
! end of my_Poisson_inc_add_Periodic_2d
#endif
