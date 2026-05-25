#include "definitions_inc.f90"
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
! Prtcl_dump_int_vector
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
subroutine Prtcl_dump_int_vector(fh,disp,var,pvsize); implicit none
  integer,intent(in)::fh
  type(part_io_size_vec),intent(in)::pvsize  
  integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
  integer,dimension(1:pvsize%subsizes(1)),intent(in)::var

  ! locals
  integer:: ierr,newtype
  integer,dimension(1) :: sizes, subsizes, starts

  ! calculate sizes, subsizes and starts
  sizes    = pvsize%sizes
  subsizes = pvsize%subsizes
  starts   = pvsize%starts

  ! write the particle revelant integer vector
  call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, int_type, newtype, ierr)
  call MPI_TYPE_COMMIT(newtype,ierr)
  call MPI_FILE_SET_VIEW(fh,disp,int_type, newtype,'native',MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),int_type, MPI_STATUS_IGNORE, ierr)
  call MPI_TYPE_FREE(newtype,ierr)
  disp = disp + int(sizes(1), MPI_OFFSET_KIND) *int(int_byte,MPI_OFFSET_KIND)
end subroutine Prtcl_dump_int_vector

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
! Prtcl_dump_int_matrix
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
subroutine Prtcl_dump_int_matrix(fh,disp,var,pmsize); implicit none
  integer,intent(in)::fh
  type(part_io_size_mat),intent(in)::pmsize    
  integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
  integer,dimension(1:pmsize%subsizes(1),1:pmsize%subsizes(2)),intent(in)::var

  ! locals
  integer:: ierr,newtype
  integer, dimension(2) :: sizes, subsizes, starts

  ! calculate sizes, subsizes and starts
  sizes     = pmsize%sizes
  subsizes  = pmsize%subsizes
  starts    = pmsize%starts

  ! write the particle relevant real matrix
  call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, int_type, newtype, ierr)
  call MPI_TYPE_COMMIT(newtype,ierr)
  call MPI_FILE_SET_VIEW(fh,disp,int_type, newtype,'native',MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(fh, var, subsizes(1)*subsizes(2),int_type, MPI_STATUS_IGNORE, ierr)
  call MPI_TYPE_FREE(newtype,ierr)
  disp = disp + int(sizes(1), MPI_OFFSET_KIND) *int(sizes(2), MPI_OFFSET_KIND) *int(int_byte,MPI_OFFSET_KIND) 
end subroutine Prtcl_dump_int_matrix

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
! Prtcl_dump_real_vector
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
subroutine Prtcl_dump_real_vector(fh,disp,var,pvsize); implicit none
  integer,intent(in)::fh
  type(part_io_size_vec),intent(in)::pvsize  
  integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
  real(RK),dimension(1:pvsize%subsizes(1)),intent(in)::var

  ! locals
  integer::ierr,newtype
  integer,dimension(1):: sizes,subsizes,starts

  ! calculate sizes, subsizes and starts
  sizes    = pvsize%sizes
  subsizes = pvsize%subsizes
  starts   = pvsize%starts

  ! write the particle revelant real vector
  call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real_type, newtype, ierr)
  call MPI_TYPE_COMMIT(newtype,ierr)
  call MPI_FILE_SET_VIEW(fh,disp,real_type, newtype,'native',MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),real_type, MPI_STATUS_IGNORE, ierr)
  call MPI_TYPE_FREE(newtype,ierr)
  disp = disp + int(sizes(1), MPI_OFFSET_KIND) *int(real_byte,MPI_OFFSET_KIND) 
end subroutine Prtcl_dump_real_vector

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
! Prtcl_dump_real3_vector
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
subroutine Prtcl_dump_real3_vector(fh,disp,var,pvsize); implicit none
  integer,intent(in)::fh
  type(part_io_size_vec),intent(in)::pvsize  
  integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
  type(real3),dimension(1:pvsize%subsizes(1)),intent(in)::var

  ! locals
  integer::ierr,newtype
  integer,dimension(1)::sizes,subsizes,starts

  ! calculate sizes, subsizes and starts
  sizes    = pvsize%sizes
  subsizes = pvsize%subsizes
  starts   = pvsize%starts

  ! write the particle revelant real3 vector
  call MPI_TYPE_CREATE_SUBARRAY(1, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real3_type, newtype, ierr)
  call MPI_TYPE_COMMIT(newtype,ierr)
  call MPI_FILE_SET_VIEW(fh,disp,real3_type, newtype,'native',MPI_INFO_NULL,ierr)
  call MPI_FILE_WRITE_ALL(fh, var, subsizes(1),real3_type, MPI_STATUS_IGNORE, ierr)
  call MPI_TYPE_FREE(newtype,ierr)
  disp = disp +int(sizes(1), MPI_OFFSET_KIND) *int(real3_byte,MPI_OFFSET_KIND) 
end subroutine Prtcl_dump_real3_vector

#ifdef CFDDEM
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_dump_real3_matrix
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Prtcl_dump_real3_matrix(fh,disp,var,pmsize); implicit none
    integer,intent(in)::fh
    type(part_io_size_mat),intent(in)::pmsize    
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp
    type(real3),dimension(1:pmsize%subsizes(1),1:pmsize%subsizes(2)),intent(in)::var

    ! locals
    integer::ierr,newtype
    integer,dimension(2)::sizes,subsizes,starts

    ! calculate sizes, subsizes and starts
    sizes     = pmsize%sizes
    subsizes  = pmsize%subsizes
    starts    = pmsize%starts

    ! write the particle relevant real matrix
    call MPI_TYPE_CREATE_SUBARRAY(2, sizes, subsizes, starts, MPI_ORDER_FORTRAN, real3_type, newtype, ierr)
    call MPI_TYPE_COMMIT(newtype,ierr)
    call MPI_FILE_SET_VIEW(fh,disp,real3_type, newtype,'native',MPI_INFO_NULL,ierr)
    call MPI_FILE_WRITE_ALL(fh, var, subsizes(1)*subsizes(2),real3_type, MPI_STATUS_IGNORE, ierr)
    call MPI_TYPE_FREE(newtype,ierr)
    disp = disp + int(sizes(1), MPI_OFFSET_KIND) *int(sizes(2), MPI_OFFSET_KIND)* int(real3_byte, MPI_OFFSET_KIND)
  end subroutine Prtcl_dump_real3_matrix
#endif
