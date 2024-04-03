!******************************************************************
! Zheng Gong, 2023-07-18
!******************************************************************  

! Constant number
integer, parameter :: C_FFTW_R2R_KIND = C_INT32_T
integer(C_INT), parameter :: FFTW_R2HC = 0
integer(C_INT), parameter :: FFTW_HC2R = 1
integer(C_INT), parameter :: FFTW_REDFT01 = 4
integer(C_INT), parameter :: FFTW_REDFT10 = 5
integer(C_INT), parameter :: FFTW_FORWARD = -1
integer(C_INT), parameter :: FFTW_BACKWARD = +1
integer(C_INT), parameter :: FFTW_MEASURE = 0
integer(C_INT), parameter :: FFTW_EXHAUSTIVE = 8
integer(C_INT), parameter :: FFTW_PATIENT = 32
integer(C_INT), parameter :: FFTW_ESTIMATE = 64    


INTERFACE
!========== C-Fortran interface for creating FFTW_PLAN
type(C_PTR) function fftw_plan_r2r_1d(n,VecIn,VecOut,FFTW_kind,flags) bind(C, name='fftw_plan_r2r_1d')
  import
  integer(C_INT), value :: n
  real(C_DOUBLE), dimension(*), intent(out) :: VecIn
  real(C_DOUBLE), dimension(*), intent(out) :: VecOut
  integer(C_FFTW_R2R_KIND), value :: FFTW_kind
  integer(C_INT), value :: flags
end function fftw_plan_r2r_1d

type(C_PTR) function fftw_plan_dft_r2c_1d(n,VecIn,VecOut,flags) bind(C, name='fftw_plan_dft_r2c_1d')
  import
  integer(C_INT), value :: n
  real(C_DOUBLE), dimension(*), intent(out) :: VecIn
  real(C_DOUBLE), dimension(*), intent(out) :: VecOut
  !complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: VecOut
  integer(C_INT), value :: flags
end function fftw_plan_dft_r2c_1d
    
type(C_PTR) function fftw_plan_dft_c2r_1d(n,VecIn,VecOut,flags) bind(C, name='fftw_plan_dft_c2r_1d')
  import
  integer(C_INT), value :: n
  !complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: VecIn
  real(C_DOUBLE), dimension(*), intent(out) :: VecIn
  real(C_DOUBLE), dimension(*), intent(out) :: VecOut
  integer(C_INT), value :: flags
end function fftw_plan_dft_c2r_1d

!========== C-Fortran interface for executing   
subroutine fftw_execute_r2r(p,VecIn,VecOut) bind(C, name='fftw_execute_r2r')
  import
  type(C_PTR), value :: p
  real(C_DOUBLE), dimension(*), intent(inout) :: VecIn
  real(C_DOUBLE), dimension(*), intent(out) :: VecOut
end subroutine fftw_execute_r2r

subroutine fftw_execute_dft_r2c(p,VecIn,VecOut) bind(C, name='fftw_execute_dft_r2c')
  import
  type(C_PTR), value :: p
  real(C_DOUBLE), dimension(*), intent(inout) :: VecIn
  real(C_DOUBLE), dimension(*), intent(out) :: VecOut
  !complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: VecOut
end subroutine fftw_execute_dft_r2c
    
subroutine fftw_execute_dft_c2r(p,VecIn,VecOut) bind(C, name='fftw_execute_dft_c2r')
  import
  type(C_PTR), value :: p
  !complex(C_DOUBLE_COMPLEX), dimension(*), intent(inout) :: VecIn
  real(C_DOUBLE), dimension(*), intent(inout) :: VecIn
  real(C_DOUBLE), dimension(*), intent(out) :: VecOut
end subroutine fftw_execute_dft_c2r    

!========== C-Fortran interface for nullifing FFTW_PLAN
subroutine fftw_destroy_plan(p) bind(C, name='fftw_destroy_plan')
  import
  type(C_PTR), value :: p
end subroutine fftw_destroy_plan

END INTERFACE
