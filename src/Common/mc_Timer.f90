#include "definitions_inc.f90"
module mc_Timer
#ifdef Compiled_With_MPI
  use MPI
#endif
  implicit none
  private

#define RKTimer 8
  type,public:: Timer
    integer:: numCycle  = 0
    real(RKTimer)::tot_time  = real(0.0, RKTimer)
    real(RKTimer)::last_time = real(0.0, RKTimer)
    real(RKTimer)::time_start= real(0.0, RKTimer)
    real(RKTimer)::time_end  = real(0.0, RKTimer)
  contains
    procedure:: reset  => t_reset    ! resetting the timer
    procedure:: start  => t_start    ! start of execution
    procedure:: finish => t_finish   ! end of the execution, end of event
    procedure:: average=> t_average  ! average execution time per event
  end type Timer
    
  public:: time2str ! a character(len=16), e.g., 2022-07-29 14:36
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! t_reset
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine t_reset(this); implicit none
    class(timer):: this
    this%numCycle = 0
    this%tot_time = real(0.0, RKTimer)
    this%last_time= real(0.0, RKTimer)
  end subroutine t_reset

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! t_start
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine t_start(this); implicit none
    class(timer):: this
#ifdef Compiled_With_MPI
    this%time_start = MPI_WTIME()
#else
    call cpu_time(this%time_start)
#endif
  end subroutine t_start

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! t_finish
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine t_finish(this); implicit none
    class(timer):: this
#ifdef Compiled_With_MPI
    this%time_end = MPI_WTIME()
#else
    call cpu_time(this%time_end)
#endif
    this%last_time= this%time_end - this%time_start
    this%tot_time = this%tot_time + this%last_time
    this%numCycle = this%numCycle + 1
  end subroutine t_finish

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! t_average
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  real(RKTimer) function t_average(this); implicit none
    class(timer):: this
    if(this%numCycle==0) then
      t_average= real(0.0, RKTimer)
    else
      t_average= this%tot_time/real(this%numCycle, RKTimer)
    endif
  end function t_average

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! time2str
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90  
  function time2str() result(res); implicit none
    character(len=16)::res
    integer::TimeValue(8)
    call date_and_time(values=TimeValue)
    write(res,'(I4.4,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')TimeValue(1),'-', &
      TimeValue(2),'-',TimeValue(3),' ',TimeValue(5),':',TimeValue(6)
  end function time2str
#undef RKTimer
end module mc_Timer
