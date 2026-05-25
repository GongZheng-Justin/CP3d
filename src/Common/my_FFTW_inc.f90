#include "definitions_inc.f90"
#ifdef my_FFTW_inc_add_y
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! my_execute_FFTW_r2r_y
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine my_execute_FFTW_r2r_y(plan_flag,nx_FFT,ny_FFT,nz_FFT,arrFFT); implicit none
    type(C_PTR),intent(in)::plan_flag
    integer,intent(in)::nx_FFT,ny_FFT,nz_FFT
    real(RK),dimension(nx_FFT,ny_FFT,nz_FFT),intent(inout)::arrFFT
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(ny_FFT)::VecIn,VecOut
    
    do kc=1,nz_FFT
      do ic=1,nx_FFT
        do jc=1,ny_FFT
          VecIn(jc)=arrFFT(ic,jc,kc)
        enddo
        call dfftw_execute_r2r(plan_flag,VecIn,VecOut)
        do jc=1,ny_FFT
          arrFFT(ic,jc,kc)=VecOut(jc)
        enddo
      enddo
    enddo
  end subroutine my_execute_FFTW_r2r_y
#endif

#ifdef my_FFTW_inc_add_z2
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! my_execute_FFTW_r2r_z_2
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine my_execute_FFTW_r2r_z_2(plan_flag,nx_FFT,ny_FFT,nz_FFT,arrFFT); implicit none
    type(C_PTR),intent(in)::plan_flag
    integer,intent(in)::nx_FFT,ny_FFT,nz_FFT
    real(RK),dimension(nx_FFT,ny_FFT,nz_FFT),intent(inout)::arrFFT
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(nz_FFT,nx_FFT)::VecIn,VecOut
    
    do jc=1,ny_FFT
      do kc=1,nz_FFT
        do ic=1,nx_FFT
          VecIn(kc,ic)=arrFFT(ic,jc,kc)
        enddo
      enddo
      do ic=1,nx_FFT
        call dfftw_execute_r2r(plan_flag,VecIn(:,ic),VecOut(:,ic))
      enddo
      do kc=1,nz_FFT
        do ic=1,nx_FFT
          arrFFT(ic,jc,kc)=VecOut(kc,ic)
        enddo
      enddo
    enddo
  end subroutine my_execute_FFTW_r2r_z_2
#endif
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! my_execute_FFTW_r2r_x
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine my_execute_FFTW_r2r_x(plan_flag,nx_FFT,ny_FFT,nz_FFT,arrFFT); implicit none
    type(C_PTR),intent(in)::plan_flag
    integer,intent(in)::nx_FFT,ny_FFT,nz_FFT
    real(RK),dimension(nx_FFT,ny_FFT,nz_FFT),intent(inout)::arrFFT
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(nx_FFT)::VecIn,VecOut
    
    do kc=1,nz_FFT
      do jc=1,ny_FFT
        do ic=1,nx_FFT
          VecIn(ic)=arrFFT(ic,jc,kc)
        enddo
        call dfftw_execute_r2r(plan_flag,VecIn,VecOut)
        do ic=1,nx_FFT
          arrFFT(ic,jc,kc)=VecOut(ic)
        enddo
      enddo
    enddo
  end subroutine my_execute_FFTW_r2r_x
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! my_execute_FFTW_r2r_z
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine my_execute_FFTW_r2r_z(plan_flag,nx_FFT,ny_FFT,nz_FFT,arrFFT); implicit none
    type(C_PTR),intent(in)::plan_flag
    integer,intent(in)::nx_FFT,ny_FFT,nz_FFT
    real(RK),dimension(nx_FFT,ny_FFT,nz_FFT),intent(inout)::arrFFT
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(nz_FFT)::VecIn,VecOut
    
    do jc=1,ny_FFT
      do ic=1,nx_FFT
        do kc=1,nz_FFT
          VecIn(kc)=arrFFT(ic,jc,kc)
        enddo
        call dfftw_execute_r2r(plan_flag,VecIn,VecOut)
        do kc=1,nz_FFT
          arrFFT(ic,jc,kc)=VecOut(kc)
        enddo
      enddo
    enddo
  end subroutine my_execute_FFTW_r2r_z
