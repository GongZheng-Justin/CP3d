module m_TypeDef 
  implicit none
  public

  integer,parameter::  RK = KIND(0.0D0)
  real(RK),parameter:: Pi   = 3.141592653589793238462643383279502884_RK    ! Pi constant
  real(RK),parameter:: twoPi= 6.283185307179586476925286766559005768_RK

  real(RK),parameter:: zero     =  0.000_RK
  real(RK),parameter:: quarter  =  0.250_RK
  real(RK),parameter:: oneEighth=  0.125_RK
  real(RK),parameter:: half     =  0.500_RK
  real(RK),parameter:: one      =  1.000_RK
  real(RK),parameter:: two      =  2.000_RK
  real(RK),parameter:: three    =  3.000_RK
  real(RK),parameter:: four     =  4.000_RK
  real(RK),parameter:: five     =  5.000_RK
  real(RK),parameter:: eight    =  8.000_RK
  real(RK),parameter:: ten      = 10.000_RK    
  real(RK),parameter:: twelve   = 12.000_RK
  real(RK),parameter:: fifteen  = 15.000_RK
  real(RK),parameter:: sixteen  = 16.000_RK
  real(RK),parameter:: seventeen= 17.000_RK
  real(RK),parameter:: sixty    = 60.000_RK

  interface num2str
    module procedure num2strI, num2strR4, num2strR8
  end interface

contains

  ! integer number to string 
  character(64) function num2strI( i )
    implicit none
    integer, intent(in) :: i
    character(64):: ch
    
    write(ch,"(I20)") i
    num2strI = trim(adjustl(ch))
  end function

  ! float number to string 
  character(64) function num2strR4( i )
    implicit none
    real(4), intent(in) :: i
    character(64):: ch
    
    write(ch,"(ES15.6)") i
    num2strR4 = trim(adjustl(ch))
  end function

  ! double number to string  
  character(64) function num2strR8( i )
    implicit none
    real(8), intent(in) :: i
    character(64):: ch
    
    write(ch,"(ES20.10)") i
    num2strR8 = trim(adjustl(ch))
  end function

end module m_TypeDef 
