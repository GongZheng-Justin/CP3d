!==================== TODO ====================
module m_FileOperation
  implicit none
  private


  public:: is_directory_exist
contains

  !******************************************************************
  ! is_directory_exist
  !******************************************************************
  function is_directory_exist(dirStr) result(flagout)
    implicit none
    character(*),intent(in)::dirStr
    
    ! locals
    logical::flagout
    
    
  
  end function is_directory_exist
  
end module m_FileOperation
