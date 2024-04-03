    !=========================== Calculate Energy Spectra ===========================!
    IF(mod(itime,ivSpec)==0) THEN
    
    ! Determine ux in cell center
#ifdef EnergySpectra_staggered_4th
      DO kc=y1start(3),y1end(3)
        kt=kc-y1start(3)+1
        do jc=y1start(2),y1end(2)
          jt=jc
          do ic=y1start(1),y1end(1)
            it=ic-y1start(1)+1
            im=ic-1;ip=ic+1;iu=ic+2
            ArrTemp1(it,jt,kt)=InterpCoe1*ux(im,jc,kc) +InterpCoe2*ux(ic,jc,kc) + &
                               InterpCoe3*ux(ip,jc,kc) +InterpCoe4*ux(iu,jc,kc) + uCRF
          enddo
        enddo
      ENDDO
#else
      DO kc=y1start(3),y1end(3)
        kt=kc-y1start(3)+1
        do jc=y1start(2),y1end(2)
          jt=jc
          do ic=y1start(1),y1end(1)
            it=ic-y1start(1)+1
            ip=ic+1;
            ArrTemp1(it,jt,kt)=0.5_RK*ux(ic,jc,kc) +0.5_RK*ux(ip,jc,kc)
          enddo
        enddo
      ENDDO
#endif
      
      !=========================== Calculate 1-D Energy Spectra
      IF(clcSpectra1D) THEN
        !=============== Spectra in x-dir ===============
        allocate(arrx1(x1size(1),x1size(2),x1size(3)))
        allocate(arrx2(x1size(1),x1size(2),x1size(3)))
        call transpose_y1_to_x1(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1)
#ifdef EnergySpectra_staggered_4th
        arrx1=arrx1+uCRF
#endif
        call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DUU)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DUU,iLCSI1DUU,iLCSR1DUU2) ! LCS -Real -Imag -ux
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DUU,iLCSR2DUU)
        
        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc; jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=0.5_RK*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO
        call transpose_y1_to_x1(ArrTemp2,arrx2)
        call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx2)
        call clcCospectraX(arrx1,arrx2,iSpec1DUV2)
        call transpose_y1_to_x1(ArrTemp1,arrx1)
        call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx1)
        call clcCospectraX(arrx1,arrx2,iSpec1DUV)    
#ifdef ScalarFlow
        call transpose_y1_to_x1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1)
        call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx1)
        call clcCosImagX(arrx1,arrx2,iSpec1DVC,iImag1DVC)
#endif
        call transform_uy(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),ArrTemp2)
        call transpose_y1_to_x1(ArrTemp2,arrx1)
        call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DVV)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DVV,iLCSI1DVV,iLCSR1DVV2) ! LCS -Real -Imag -uy
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DVV,iLCSR2DVV)
        
        call transpose_y1_to_x1(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1)
        call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DWW)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DWW,iLCSI1DWW,iLCSR1DWW2) ! LCS -Real -Imag -uz
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DWW,iLCSR2DWW)
        
        call transpose_y1_to_x1(pr(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1)
        call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DPP)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DPP,iLCSI1DPP,iLCSR1DPP2) ! LCS -Real -Imag -pp
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DPP,iLCSR2DPP)
#ifdef ScalarFlow
        call transpose_y1_to_x1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1)
        call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DCC)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DCC,iLCSI1DCC,iLCSR1DCC2) ! LCS -Real -Imag -cc
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DCC,iLCSR2DCC)
        
        call transpose_y1_to_x1(ArrTemp1,arrx2)
        call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arrx2)
        call clcCosImagX(arrx1,arrx2,iSpec1DUC,iImag1DUC)
#endif
        deallocate(arrx1,arrx2)
        !
        !=============== Spectra in z-dir ===============
        allocate(arrz1(z1size(1),z1size(2),z1size(3)))
        allocate(arrz2(z1size(1),z1size(2),z1size(3)))
        call transpose_y1_to_z1(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1)
#ifdef EnergySpectra_staggered_4th
        arrz1=arrz1+uCRF
#endif
        call my_execute_FFTW_r2r_z(fft_plan_z,z1size(1),z1size(2),z1size(3),arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DUU)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DUU,iLCSI1DUU,iLCSR1DUU2) ! LCS -Real -Imag -ux

        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc; jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=0.5_RK*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO
        call transpose_y1_to_z1(ArrTemp2,arrz2)
        call my_execute_FFTW_r2r_z(fft_plan_z,z1size(1),z1size(2),z1size(3),arrz2)
        call clcCospectraZ(arrz1,arrz2,iSpec1DUV2)
        call transpose_y1_to_z1(ArrTemp1,arrz1)
        call my_execute_FFTW_r2r_z(fft_plan_z,z1size(1),z1size(2),z1size(3),arrz1)
        call clcCospectraZ(arrz1,arrz2,iSpec1DUV)  
#ifdef ScalarFlow
        call transpose_y1_to_z1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1)
        call my_execute_FFTW_r2r_z(fft_plan_z,z1size(1),z1size(2),z1size(3),arrz1)
        call clcCosImagZ(arrz1,arrz2,iSpec1DVC,iImag1DVC)
#endif

        call transform_uy(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),ArrTemp2)
        call transpose_y1_to_z1(ArrTemp2,arrz1)
        call my_execute_FFTW_r2r_z(fft_plan_z,z1size(1),z1size(2),z1size(3),arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DVV)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DVV,iLCSI1DVV,iLCSR1DVV2) ! LCS -Real -Imag -uy
        
        call transpose_y1_to_z1(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1)
        call my_execute_FFTW_r2r_z(fft_plan_z,z1size(1),z1size(2),z1size(3),arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DWW)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DWW,iLCSI1DWW,iLCSR1DWW2) ! LCS -Real -Imag -uz

        call transpose_y1_to_z1(pr(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1)
        call my_execute_FFTW_r2r_z(fft_plan_z,z1size(1),z1size(2),z1size(3),arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DPP)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DPP,iLCSI1DPP,iLCSR1DPP2) ! LCS -Real -Imag -pp
#ifdef ScalarFlow
        call transpose_y1_to_z1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1)
        call my_execute_FFTW_r2r_z(fft_plan_z,z1size(1),z1size(2),z1size(3),arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DCC)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DCC,iLCSI1DCC,iLCSR1DCC2) ! LCS -Real -Imag -cc
        
        call transpose_y1_to_z1(ArrTemp1,arrz2)
        call my_execute_FFTW_r2r_z(fft_plan_z,z1size(1),z1size(2),z1size(3),arrz2)
        call clcCosImagZ(arrz1,arrz2,iSpec1DUC,iImag1DUC)
#endif
        deallocate(arrz1,arrz2)
      ENDIF
      !=========================== Calculate 1-D Energy Spectra

      !=========================== Calculate 2-D Energy Spectra
      IF(clcSpectra2D) THEN
        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc
            jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=0.5_RK*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO
        call clcCosSpectra2D(ArrTemp1,ArrTemp2,iSpec2DUV,.false.)
        if(.not. clcSpectra1D) then
          call transform_uy(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),ArrTemp2)
#ifdef EnergySpectra_staggered_4th
          call clcSpectraAndLCSR2D(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))+uCRF,  &
                                   iSpec2DUU,iLCSR2DUU)
#else
          call clcSpectraAndLCSR2D(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),       &
                                   iSpec2DUU,iLCSR2DUU)
#endif
          call clcSpectraAndLCSR2D(ArrTemp2,                                                              &
                                   iSpec2DVV,iLCSR2DVV)
          call clcSpectraAndLCSR2D(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),       &
                                   iSpec2DWW,iLCSR2DWW)
          call clcSpectraAndLCSR2D(pr(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)), &
                                   iSpec2DPP,iLCSR2DPP)
#ifdef ScalarFlow
          call clcSpectraAndLCSR2D(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),   &
                                   iSpec2DCC,iLCSR2DCC)
#endif
        endif
      ENDIF
      !=========================== Calculate 2-D Energy Spectra
      nSpectime= nSpectime+1
    ENDIF
