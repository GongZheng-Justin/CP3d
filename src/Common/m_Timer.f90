module m_Timer
#ifdef Compiled_With_MPI
  use MPI
#endif
  use m_TypeDef,only:RK
  implicit none
  private
  
  type,public:: Timer
    integer:: numCycle  = 0
    real(RK)::tot_time  = 0.0_RK
    real(RK)::last_time = 0.0_RK
    real(RK)::time_start= 0.0_RK
    real(RK)::time_end  = 0.0_RK
  contains
    procedure:: reset  => t_reset    ! resetting the timer
    procedure:: start  => t_start    ! start of execution
    procedure:: finish => t_finish   ! end of the execution, end of event
    procedure:: average=> t_average  ! average execution time per event
  end type Timer
    
  public:: time2str ! a character(len=16), e.g., 2022-07-29 14:36
contains

  !*****************************************************************
  ! t_reset
  !*****************************************************************
  subroutine t_reset(this)
    implicit none
    class(timer) this
    this%numCycle = 0
    this%tot_time = 0.0_RK
    this%last_time = 0.0_RK
  end subroutine t_reset

  !*****************************************************************
  ! t_start
  !*****************************************************************
  subroutine t_start(this)
    implicit none
    class(timer) this
#ifdef Compiled_With_MPI
    this%time_start = MPI_WTIME()
#else
    call cpu_time(this%time_start)
#endif
  end subroutine t_start

  !*****************************************************************
  ! t_finish
  !!*****************************************************************
  subroutine t_finish(this)
    implicit none
    class(timer) this
#ifdef Compiled_With_MPI
    this%time_end = MPI_WTIME()
#else
    call cpu_time(this%time_end)
#endif
    this%last_time= this%time_end - this%time_start
    this%tot_time = this%tot_time+ this%last_time
    this%numCycle = this%numCycle + 1
  end subroutine t_finish

  !*****************************************************************
  ! t_average
  !*****************************************************************
  real(RK) function t_average(this)
    implicit none
    class(timer) this
    if(this%numCycle==0)then
      t_average= 0.0_RK
    else
      t_average= this%tot_time/real(this%numCycle,kind=RK)
    endif
  end function t_average

  !*****************************************************************
  ! time2str
  !*****************************************************************  
  function time2str() result(res)
    implicit none
    character(len=16)::res
    integer::TimeValue(8)
    call date_and_time(values=TimeValue)
    write(res,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')TimeValue(1),'-', &
      TimeValue(2),'-',TimeValue(3),' ',TimeValue(5),':',TimeValue(6)
  end function time2str
end module m_Timer
