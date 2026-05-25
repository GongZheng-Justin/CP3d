#include "definitions_inc.f90"
module mc_FileOperator
  use iso_c_binding, only:c_null_char, c_int16_t
  implicit none
  private

  interface 
    function mkdir_c(path,mode) result(mkdir_stat) bind(c,name='mkdir') 
      use iso_c_binding, only:c_char, c_int16_t, c_int
      character(kind=c_char),intent(in):: path(*) 
      integer(c_int16_t), intent(in),value:: mode
      integer(c_int):: mkdir_stat
    end function mkdir_c
    
    function access_c(path, mode) result(access_stat) bind(c,name='access')
      use iso_c_binding, only:c_char, c_int16_t, c_int
      character(kind=c_char),intent(in):: path(*) 
      integer(c_int16_t), intent(in),value:: mode      
      integer(c_int):: access_stat
    end function access_c
  end interface
    
  public:: mkdir, access_
contains
#define F_OK 0

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! access_
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  logical function access_(path); implicit none
    character(len=*),intent(in)::path
    if(access_c(trim(adjustl(path)) // c_null_char, int(F_OK, c_int16_t) ) ==0) then
      access_ = .true.
    else
      access_ = .false.
    endif
  end function access_
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! mkdir
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine mkdir(StrIn, stat); implicit none
    character(len=*),intent(in)::StrIn
    integer,intent(out)::stat
    
    ! locals
    integer :: i, nLen, ierr
    character(len=:),allocatable::path,path_use

    path = trim(adjustl(StrIn))
    nLen = len(path)
    if(nLen <1) then
      stat = -1; return
    endif
    do i=1, nLen
      if(path(i:i)==achar(92)) path(i:i)='/'
    enddo
    if(nLen ==1 .and. path(1:1)=='/') then
      stat = -2; return
    endif
    if(path(nLen:nLen) /= '/') then
      nLen = nLen +1
      path_use = path // '/'
    else
      path_use = path
    endif
    path = path_use
    if(access_(path)) then
      stat=1; return
    endif
    
    do i=1, nLen
      if(path_use(i:i) =='/') then
        path = path_use(1:i)
        if(.not. access_(path)) then
          ierr = mkdir_c(path // c_null_char, int(o'755', c_int16_t))
          if(ierr /=0) then
            stat =-3; return
          endif
        endif
      endif
    enddo
    
    if(.not. access_(path_use)) then 
      stat = -4
    else
      stat =  0
    endif
  end subroutine mkdir
  
#undef F_OK
end module mc_FileOperator
