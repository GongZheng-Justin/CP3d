module m_TypeDef
  implicit none
  public
  integer,parameter::  RK = KIND(0.0D0)
  real(RK),parameter:: Pi = 3.141592653589793238462643383279502884_RK    ! Pi constant
  interface num2str
    module procedure num2strI, num2strR4, num2strR8
#ifdef VectorOperator
    module procedure num2strInt3
#endif
  end interface num2str

#ifdef VectorOperator
  type integer3
    integer::x
    integer::y
    integer::z
  end type integer3
  type real3
    real(RK)::x
    real(RK)::y
    real(RK)::z
  end type real3  
  type real4
    real(RK)::x
    real(RK)::y
    real(RK)::z
    real(RK)::w
  end type real4
  type(real3),parameter:: zero_r3=real3(0.0_RK,0.0_RK,0.0_RK)
  type(real4),parameter:: zero_r4=real4(0.0_RK,0.0_RK,0.0_RK,0.0_RK)
    
  interface operator(+)
    module procedure real3_add1, real3_add2, real3_add3, integer3_add
    module procedure real3_add_real4, real4_add_real3
  end interface

  interface operator(-)
    module procedure real3_sub
  end interface

  interface operator(*)      
    module procedure real3_prod1, real3_prod2, real3_prod3
  end interface

  interface operator(/)
    module procedure  real3_division2
  end interface
    
  interface assignment(=)
    module procedure real3_real4, real4_real3, real3_array3
  end interface

  interface operator(-)
    module procedure real4_sub
  end interface
    
  interface operator(.dot.)
    module procedure real3_dotProd
  end interface

  interface operator(.cross.)
    module procedure real3_CrossProd
  end interface

  interface operator(.ovlp.)
    module procedure real4_ovlp
  end interface
    
  interface operator(.dist.)
    module procedure real4_dist, real3_dist
  end interface
    
  interface operator(.nv.)
    module procedure real4_nv, real3_nv
  end interface
    
  interface norm
    module procedure real4_norm, real3_norm
  end interface
#endif

contains
  !*****************************************************************
  ! integer number to string
  !*****************************************************************
  character(64) function num2strI( i )
    implicit none
    integer,intent(in) :: i
    character(64):: ch
    
    write(ch,"(I20)") i
    num2strI=trim(adjustl(ch))
  end function num2strI

  !*****************************************************************
  ! float number to string 
  !*****************************************************************
  character(64) function num2strR4( i )
    implicit none
    real(4), intent(in) :: i
    character(64):: ch
    
    write(ch,"(ES15.6)") i
    num2strR4=trim(adjustl(ch))
  end function num2strR4

  !*****************************************************************
  ! double number to string 
  !*****************************************************************
  character(64) function num2strR8(i)
    implicit none
    real(8),intent(in)::i
    character(64)::ch
    
    write(ch,"(ES20.10)") i
    num2strR8=trim(adjustl(ch))
  end function num2strR8
  
#ifdef VectorOperator
  !*****************************************************************
  ! integer3 variable to string 
  !*****************************************************************
  character(64) function num2strInt3(i)
    implicit none
    type(integer3),intent(in)::i
    character(64)::ch    
    write(ch,*)i%x,",",i%y,",",i%z
    num2strInt3="["// trim(adjustl(ch)) // "]"
  end function num2strInt3

  !*****************************************************************
  ! Vector1 + Vector2 
  !*****************************************************************  
  function real3_add1(p1,p2) result(p3)
    implicit none
    type(real3),intent(in)::p1,p2
    type(real3)::p3
    p3= real3(p1%x+p2%x, p1%y+p2%y, p1%z+p2%z) 
  end function real3_add1

  !*****************************************************************
  ! Vector + Scalar
  !*****************************************************************
  function real3_add2(p1,v) result(p3)
    implicit none
    type(real3),intent(in)::p1
    real(RK),intent(in)::v
    type(real3):: p3
    p3= real3(p1%x+ v, p1%y+v, p1%z+v) 
  end function real3_add2

  !*****************************************************************
  ! Scalar + Vector
  !*****************************************************************
  function real3_add3(v,p1) result(p3)
    implicit none
    real(RK),intent(in)::v
    type(real3),intent(in)::p1
    type(real3)::p3
    p3= real3(p1%x+v, p1%y+v , p1%z+v) 
  end function real3_add3

  !*****************************************************************
  ! Vector + LongVector
  !*****************************************************************
  function real3_add_real4(p1,p2) result(p3)
    implicit none
    type(real3),intent(in)::p1
    type(real4),intent(in)::p2
    type(real4)::p3
    p3=real4(p1%x+p2%x, p1%y+p2%y, p1%z+p2%z, p2%w) 
  end function real3_add_real4

  !*****************************************************************
  ! LongVector + Vector
  !*****************************************************************   
  function real4_add_real3(p1,p2) result(p3)
    implicit none
    type(real4),intent(in)::p1
    type(real3),intent(in)::p2
    type(real4)::p3
    p3= real4(p1%x+p2%x, p1%y+p2%y, p1%z+p2%z, p1%w) 
  end function real4_add_real3   

  !*****************************************************************
  ! Vector1 - Vector2
  !*****************************************************************    
  function real3_sub(p1,p2) result(p3)
    implicit none
    type(real3),intent(in)::p1,p2
    type(real3)::p3
    p3 = real3(p1%x-p2%x, p1%y-p2%y, p1%z-p2%z) 
  end function real3_sub

  !*****************************************************************
  ! real3_prod1
  !***************************************************************** 
  function real3_prod1(p1,p2) result(p3)
    implicit none
    type(real3),intent(in)::p1,p2
    type(real3)::p3
    p3 = real3(p1%x*p2%x, p1%y*p2%y, p1%z*p2%z)
  end function real3_prod1

  !*****************************************************************
  ! real3_prod2
  !*****************************************************************
  function real3_prod2(p1,v) result(p3)
    implicit none
    type(real3),intent(in)::p1
    real(RK),intent(in)::v
    type(real3)::p3
    p3 = real3( p1%x*v, p1%y*v, p1%z*v) 
  end function real3_prod2

  !*****************************************************************
  ! real3_prod3
  !*****************************************************************
  function real3_prod3(v,p1) result(p3)
    implicit none
    type(real3),intent(in)::p1
    real(RK),intent(in)::v
    type(real3)::p3
    p3 = real3(p1%x*v, p1%y*v, p1%z*v) 
  end function real3_prod3

  !*****************************************************************
  ! real3_division2
  !*****************************************************************
  function real3_division2(p1,v) result(p3)
    implicit none
    type(real3),intent(in)::p1
    real(RK),intent(in)::v
    type(real3) p3
    p3 = real3( p1%x/v, p1%y/v, p1%z/v)
  end function

  !*****************************************************************
  ! real3_real4
  !*****************************************************************
  subroutine real3_real4(x,y)
    implicit none
    type(real3),intent(out)::x
    type(real4),intent(in) ::y
     x = real3(y%x, y%y, y%z)
  end subroutine real3_real4

  !*****************************************************************
  ! real4_real3
  !*****************************************************************    
  subroutine real4_real3(x,y)
    implicit none
    type(real4),intent(out):: x
    type(real3),intent(in) :: y
    x = real4(y%x, y%y, y%z, x%w)
  end subroutine real4_real3

  !*****************************************************************
  ! real3_array3
  !*****************************************************************
  subroutine real3_array3(x,y)
    implicit none
    type(real3),intent(out)::x
    real(RK),dimension(3),intent(in)::y
    x = real3(y(1), y(2), y(3))
  end subroutine real3_array3    

  !*****************************************************************
  ! real3_dotProd
  !*****************************************************************    
  function real3_dotProd(p1,p2) result(v)
    implicit none
    type(real3),intent(in)::p1,p2
    real(RK)::v
    v = p1%x*p2%x + p1%y*p2%y + p1%z*p2%z
  end function real3_dotProd

  !*****************************************************************
  ! real3_CrossProd
  !*****************************************************************
  function real3_CrossProd(u,v) result(p)
    implicit none
    type(real3),intent(in)::u,v
    type(real3):: p
    p = real3( u%y*v%z-u%z*v%y,  u%z*v%x-u%x*v%z,  u%x*v%y-u%y*v%x )
  end function real3_CrossProd

  !*****************************************************************
  ! real3_dist
  !*****************************************************************    
  function real3_dist(p1,p2) result(dist)
    implicit none
    type(real3),intent(in):: p1,p2
    real(RK)::dist
    dist = norm(p1-p2)
  end function real3_dist

  !*****************************************************************
  ! real3_nv
  !*****************************************************************    
  function real3_nv(p1,p2) result(nv)
    implicit none
    type(real3),intent(in)::p1,p2
    type(real3)::nv
    type(real3)::s
    real(RK)::normv
                
    s = p1-p2
    normv=norm(s)
    if( normv>1.0E-50_RK) then
      nv = s / normv
    else
      nv = zero_r3
    endif
  end function real3_nv

  !*****************************************************************
  ! real3_norm
  !*****************************************************************    
  function real3_norm(p1) result (norm)
    implicit none
    type(real3),intent(in):: p1
    real(RK)::norm
    norm = sqrt(p1%x*p1%x + p1%y*p1%y + p1%z*p1%z)
  end function

  !*****************************************************************
  ! real4_sub
  !*****************************************************************        
  function real4_sub( p1, p2 ) result(p3)
    implicit none
    type(real4),intent(in)::  p1, p2
    type(real4)::p3
    p3 = real4( p1%x-p2%x , p1%y - p2%y , p1%z - p2%z, p1%w - p2%w ) 
  end function    

  !*****************************************************************
  ! real4_ovlp
  !*****************************************************************
  function real4_ovlp(p1, p2 ) result(ovlp)
    implicit none
    type(real4),intent(in):: p1, p2
    real(RK)::ovlp
    ovlp = p1%w+p2%w - norm(p1-p2)
  end function real4_ovlp

  !*****************************************************************
  ! real4_dist
  !*****************************************************************
  function real4_dist(p1, p2 ) result(dist)
    implicit none
    type(real4),intent(in):: p1, p2
    real(RK) dist
    dist = norm(p1-p2)
  end function real4_dist

  !*****************************************************************
  ! real4_nv
  !*****************************************************************    
  function real4_nv(p1, p2) result(nv)
    implicit none
    type(real4),intent(in):: p1, p2
    type(real3)::nv,s
    real(RK)::normv
        
    s = p1-p2
    normv=norm(s)
    if( normv>1.0E-50_RK) then
      nv = s / normv
    else
      nv = zero_r3
    endif
  end function real4_nv

  !*****************************************************************
  ! real4_norm
  !*****************************************************************    
  function real4_norm(p1) result (norm)
    implicit none
    type(real4),intent(in):: p1
    real(RK)::norm
    norm = sqrt(p1%x*p1%x + p1%y*p1%y + p1%z*p1%z)
  end function real4_norm

  !*****************************************************************
  ! integer3_add
  !*****************************************************************    
  function integer3_add( p1, p2 ) result(p3)
    implicit none
    type(integer3),intent(in)::  p1, p2
    type(integer3)::p3
    p3 = integer3( p1%x+p2%x , p1%y + p2%y , p1%z + p2%z ) 
  end function integer3_add
#endif

end module m_TypeDef 
