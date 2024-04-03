module m_math
  implicit none
  private  
  integer,parameter::RK=8
  real(RK),parameter::SmallValue=1.0D-15
  
  public::SymmetricMatrixEigen3
contains

  !**************************************************************************80********90
  ! SymmetricMatrixEigen3, Zheng Gong, 2023-08-03
  !**************************************************************************80********90
  subroutine SymmetricMatrixEigen3(Mat, EigenValue, EigenVector, EigenErr)
    implicit none
    real(RK),dimension(3,3),intent(in)::Mat
    real(RK),dimension(3),intent(out)::EigenValue
    real(RK),dimension(3,3),intent(out)::EigenVector
    real(RK),optional,intent(out)::EigenErr
  
    ! locals
    real(RK)::a11,a22,a33,a12,a13,a23,EigenErr1,EigenErr2,rnorm
  
    EigenErr1=0.0_RK; EigenErr2=0.0_RK; rnorm= 0.0_RK
    a11= Mat(1,1); a22= Mat(2,2); a33= Mat(3,3)
    a12= Mat(1,2); a13= Mat(1,3); a23= Mat(2,3)
    if(abs(a11)>rnorm) rnorm= abs(a11)
    if(abs(a22)>rnorm) rnorm= abs(a22)
    if(abs(a33)>rnorm) rnorm= abs(a33)
    if(abs(a12)>rnorm) rnorm= abs(a12)
    if(abs(a13)>rnorm) rnorm= abs(a13)  
    if(abs(a23)>rnorm) rnorm= abs(a23)   
    if(rnorm < SmallValue) then
      EigenVector=0.0_RK
      EigenValue=[0.0_RK, 0.0_RK, 0.0_RK]
      EigenVector(1,1)=1.0_RK
      EigenVector(2,2)=1.0_RK
      EigenVector(3,3)=1.0_RK
      return
    endif
    a11= a11/rnorm; a22= a22/rnorm; a33= a33/rnorm
    a12= a12/rnorm; a13= a13/rnorm; a23= a23/rnorm
  
    call NISymmetricEigen3()
    call clc_EigenError(EigenErr1)
    call ItSymmetricEigen3()
    call clc_EigenError(EigenErr2)
    if(EigenErr1 < EigenErr2) then
      call NISymmetricEigen3()
      if(present(EigenErr)) EigenErr=EigenErr1
    else
      call ItSymmetricEigen3()
      if(present(EigenErr)) EigenErr=EigenErr2
    endif
    EigenValue= EigenValue*rnorm
    !print*,">>>>>>>>>>>", EigenErr1,EigenErr2
  CONTAINS
    !***************************************
    ! www.geometrictools.com/GTE/Mathematics/SymmetricMatrixEigen3x3.h
    !***************************************
    subroutine NISymmetricEigen3()
      implicit none
  
      ! locals
      real(RK)::b11,b22,b33,c11,c12,c13
      real(RK)::norm,p,q,halfdet,angle,beta1,beta2,beta3
  
      norm= a12*a12 +a13*a13 +a23*a23
      if(norm < SmallValue) then
        EigenVector=0.0_RK
        EigenValue=[a11, a22, a33]
        EigenVector(1,1)=1.0_RK
        EigenVector(2,2)=1.0_RK
        EigenVector(3,3)=1.0_RK    
      else
        q  =(a11 +a22 +a33)/3.0_RK
        b11= a11 -q
        b22= a22 -q
        b33= a33 -q
        p  = sqrt((b11*b11 +b22*b22 +b33*b33 +norm*2.0_RK)/6.0_RK)
        c11= b22*b33 -a23*a23
        c12= a12*b33 -a23*a13
        c13= a12*a23 -b22*a13
        halfdet= 0.5_RK*(b11*c11 -a12*c12 +a13*c13)/(p*p*p)
        halfdet= min(max(halfdet, -1.0_RK), 1.0_RK)
        angle  = acos(halfdet)/3.0_RK
        beta3  = cos(angle)*2.0_RK
        beta1  = cos(angle +2.0943951023931954923_RK)*2.0_RK
        beta2  = -(beta1 +beta3)
        EigenValue(1) = q +p*beta1
        EigenValue(2) = q +p*beta2
        EigenValue(3) = q +p*beta3
        if(halfdet >= 0.0_RK) then
          call ClcEigenVec0(EigenValue(3), EigenVector(:,3))
          call ClcEigenVec1(EigenVector(:,3), EigenValue(2), EigenVector(:,2))
          EigenVector(:,1)= Vector3Cross(EigenVector(:,2), EigenVector(:,3))
        else
          call ClcEigenVec0(EigenValue(1), EigenVector(:,1))
          call ClcEigenVec1(EigenVector(:,1), EigenValue(2), EigenVector(:,2))
          EigenVector(:,3)= Vector3Cross(EigenVector(:,1), EigenVector(:,2))
        endif
      endif
      call NormalizeEigen()
    end subroutine NISymmetricEigen3
    
    !***************************************
    ! people.sc.fsu.edu/~jburkardt/f_src/jacobi_eigenvalue/jacobi_eigenvalue.f90
    !***************************************
    subroutine ItSymmetricEigen3()
      implicit none

      ! locals
      integer::it_num,m,i,j,k
      real(RK)::MatIn(3,3),bw(3),row1(3)
      real(RK)::c,g,gapq,h,s,t,tau,term,termp,termq,theta,thresh

      MatIn(1,1)=a11; MatIn(1,2)=a12; MatIn(1,3)=a13
      MatIn(2,1)=a12; MatIn(2,2)=a22; MatIn(2,3)=a23
      MatIn(3,1)=a13; MatIn(3,2)=a23; MatIn(3,3)=a33
      do j=1,3
        do i=1,3
          EigenVector(i,j)= 0.0_RK
        enddo
        EigenVector(j,j) = 1.0_RK
      enddo
      do i=1,3
        EigenValue(i) = MatIn(i,i)
      enddo

      bw= EigenValue
      row1= 0.0_RK
      it_num = 0
      do while(it_num < 1000)
        it_num =it_num+1
 
        ! Convergence based on element sizes in the upper triangle of the matrix.
        thresh = 0.0_RK
        do j= 1,3
          do i= 1,j-1
            thresh =thresh +MatIn(i,j)*MatIn(i,j)
          enddo
        enddo
        thresh=sqrt(thresh)/real(12,kind=rk)
        if(thresh == 0.0_RK) exit
        do k= 1,3
          do m= k+1,3
            gapq = 10.0_RK *abs(MatIn(k,m))
            termp= gapq +abs(EigenValue(k))
            termq= gapq +abs(EigenValue(m))

            ! Annihilate tiny offdiagonal elements.
            if(it_num >4 .and. termp==abs(EigenValue(k)) .and. &
                               termq==abs(EigenValue(m))) then
              MatIn(k,m)= 0.0_RK
            elseif(thresh <=abs(MatIn(k,m))) then ! Otherwise, apply a rotation.
              h = EigenValue(m)-EigenValue(k)
              term= abs(h)+gapq

              if(term == abs(h)) then
                t = MatIn(k,m) / h
              else
                theta =0.5_RK*h/MatIn(k,m)
                t = 1.0_RK/(abs(theta)+sqrt(1.0_RK + theta*theta))
                if(theta <0.0_RK) t = -t
              endif
              c= 1.0_RK/sqrt(1.0_RK +t*t)
              s= t*c
              tau= s/(1.0_RK +c)
              h = t*MatIn(k,m)

              ! Accumulate corrections to diagonal elements.
              row1(k)= row1(k) -h                  
              row1(m)= row1(m) +h
              EigenValue(k)= EigenValue(k) -h
              EigenValue(m)= EigenValue(m) +h
              MatIn(k,m) = 0.0_RK

              ! Rotate, using information from the upper triangle of A only.
              do j =1,k-1
                g= MatIn(j,k)
                h= MatIn(j,m)
                MatIn(j,k)= g -s*(h +g*tau)
                MatIn(j,m)= h +s*(g -h*tau)
              enddo

              do j =k+1,m-1
                g= MatIn(k,j)
                h= MatIn(j,m)
                MatIn(k,j)= g -s*(h +g*tau)
                MatIn(j,m)= h +s*(g -h*tau)
              enddo

              do j =m+1,3
                g= MatIn(k,j)
                h= MatIn(m,j)
                MatIn(k,j)= g -s*(h +g*tau)
                MatIn(m,j)= h +s*(g -h*tau)
              enddo

              ! Accumulate information in the eigenvector matrix.
              do j =1,3
                g= EigenVector(j,k)
                h= EigenVector(j,m)
                EigenVector(j,k)= g -s*(h +g*tau)
                EigenVector(j,m)= h +s*(g -h*tau)
              enddo
            endif
          enddo
        enddo
        bw = bw +row1
        EigenValue = bw
        row1 = 0.0_RK
      enddo
      call NormalizeEigen()
    end subroutine ItSymmetricEigen3

    !=======================================
    subroutine NormalizeEigen()
      implicit none
    
      ! locals
      integer::i,j,k
      real(RK)::rt,row1(3)
    
      do k=1,2
        j=k
        do i=k+1,3
          if(abs(EigenValue(i)) > abs(EigenValue(j))) then
            j = i
          endif
        enddo
        if(j /= k) then
          rt= EigenValue(j)
          EigenValue(j) = EigenValue(k)
          EigenValue(k) = rt
          row1= EigenVector(:,j)
          EigenVector(:,j) = EigenVector(:,k)
          EigenVector(:,k) = row1
        endif
      enddo
      EigenVector(:,3)= Vector3Cross(EigenVector(:,1),EigenVector(:,2))
      do k=1,3
       rt=EigenVector(1,k)*EigenVector(1,k)+EigenVector(2,k)*EigenVector(2,k) &
         +EigenVector(3,k)*EigenVector(3,k)
       rt= 1.0_RK/sqrt(rt)
       EigenVector(:,k) =rt*EigenVector(:,k)
      enddo
    end subroutine NormalizeEigen
  
    !=======================================
    subroutine clc_EigenError(EigenError)
      implicit none
      real(RK),intent(out)::EigenError
  
      ! locals
      integer::j
      real(RK),dimension(3,3)::MatIn,MatError

      MatIn(1,1)=a11; MatIn(1,2)=a12; MatIn(1,3)=a13
      MatIn(2,1)=a12; MatIn(2,2)=a22; MatIn(2,3)=a23
      MatIn(3,1)=a13; MatIn(3,2)=a23; MatIn(3,3)=a33
      MatError = matmul(MatIn, EigenVector)
      do j=1,3
        MatError(:,j) = MatError(:,j) -EigenValue(j)*EigenVector(:,j)
      enddo
      EigenError =sum(abs(MatError))/9.0_RK               
    end subroutine clc_EigenError

    !=======================================
    subroutine ClcEigenVec0(eval1, evec1)
      implicit none
      real(RK),intent(in)::eval1
      real(RK),dimension(3),intent(out)::evec1
    
      ! locals
      integer::i
      real(RK)::d1,d2,d3,dmax
      real(RK),dimension(3)::row1,row2,row3, r1xr2,r1xr3,r2xr3
    
      row1=[a11-eval1, a12, a13]
      row2=[a12, a22-eval1, a23]
      row3=[a13, a23, a33-eval1]
      r1xr2= Vector3Cross(row1, row2)
      r1xr3= Vector3Cross(row1, row3)
      r2xr3= Vector3Cross(row2, row3)
      d1= r1xr2(1)*r1xr2(1) +r1xr2(2)*r1xr2(2) +r1xr2(3)*r1xr2(3)
      d2= r1xr3(1)*r1xr3(1) +r1xr3(2)*r1xr3(2) +r1xr3(3)*r1xr3(3)
      d3= r2xr3(1)*r2xr3(1) +r2xr3(2)*r2xr3(2) +r2xr3(3)*r2xr3(3)
    
      dmax=d1; i=0
      if(d2 >dmax) then
        dmax=d2; i=1
      endif
      if(d3 >dmax) i=2
      if(i==0) then
        evec1= r1xr2/sqrt(d1)
      elseif(i==1) then
        evec1= r1xr3/sqrt(d2)
      else
        evec1= r2xr3/sqrt(d3)
      endif
    end subroutine ClcEigenVec0

    !=======================================
    subroutine ClcEigenVec1(evec1, eval2, evec2)
      implicit none
      real(RK),dimension(3),intent(in)::evec1
      real(RK),intent(in)::eval2
      real(RK),dimension(3),intent(out)::evec2
    
      ! locals
      real(RK),dimension(3)::u,v,au,av
      real(RK)::invLength, m11,m12,m22, absM11,absM12,absM22, maxAbsComp

      if(abs(evec1(1)) > abs(evec1(2))) then
        invLength=1.0_RK/sqrt(evec1(1)*evec1(1) +evec1(3)*evec1(3))
        u=[-evec1(3)*invLength, 0.0_RK, evec1(1)*invLength]
      else
        invLength=1.0_RK/sqrt(evec1(2)*evec1(2) +evec1(3)*evec1(3))
        u=[0.0_RK, evec1(3)*invLength, -evec1(2)*invLength]
      endif
      v= Vector3Cross(evec1,u)
    
      au(1)= a11*u(1) +a12*u(2) +a13*u(3)
      au(2)= a12*u(1) +a22*u(2) +a23*u(3)
      au(3)= a13*u(1) +a23*u(2) +a33*u(3)
      av(1)= a11*v(1) +a12*v(2) +a13*v(3)
      av(2)= a12*v(1) +a22*v(2) +a23*v(3)
      av(3)= a13*v(1) +a23*v(2) +a33*v(3)
      m11= u(1)*au(1) +u(2)*au(2) +u(3)*au(3) -eval2
      m12= u(1)*av(1) +u(2)*av(2) +u(3)*av(3)
      m22= v(1)*av(1) +v(2)*av(2) +v(3)*av(3) -eval2
      absM11=abs(m11); absM12=abs(m12); absM22=abs(m22)
      IF(absM11 >= absM22) THEN
        maxAbsComp =max(absM11, absM12)
        if(maxAbsComp >0.0_RK) then
          if(absM11 >= absM12) then
            m12= m12/m11
            m11= 1.0_RK/sqrt(1.0_RK +m12*m12)
            m12= m12*m11
          else
            m11= m11/m12
            m12= 1.0_RK/sqrt(1.0_RK +m11*m11)
            m11= m11*m12   
          endif
          evec2= m12*u -m11*v
        else
          evec2= u
        endif
      ELSE
        maxAbsComp =max(absM22, absM12)
        if(maxAbsComp >0.0_RK) then
          if(absM22 >= absM12) then
            m12= m12/m22
            m22= 1.0_RK/sqrt(1.0_RK +m12*m12)
            m12= m12*m22
          else
            m22= m22/m12
            m12= 1.0_RK/sqrt(1.0_RK +m22*m22)
            m22= m22*m12
          endif
          evec2= m22*u -m12*v
        else
          evec2= u
        endif
      ENDIF
    end subroutine ClcEigenVec1

    !=======================================
    function Vector3Cross(u,v) result(VecOut)
      implicit none
      real(RK),dimension(3)::u,v,VecOut
      VecOut=[u(2)*v(3)-u(3)*v(2), u(3)*v(1)-u(1)*v(3), u(1)*v(2)-u(2)*v(1)]
    end function Vector3Cross
  end subroutine SymmetricMatrixEigen3

  !**************************************************************************80********90
  !
  !**************************************************************************80********90
end module m_math

#define Test_math
#ifdef Test_math
program main
  use m_math
  implicit none
  integer,parameter::RK=8 
  print*," "
 
  ! test SymmetricMatrixEigen3
  BLOCK    
    real(RK)::Mat(3,3),EigVec(3,3),EigVal(3),EigErr
    
    call random_number(Mat)
    Mat= Mat+transpose(Mat)
    call SymmetricMatrixEigen3(Mat, EigVal, EigVec, EigErr)
    print*,"SymmetricMatrixEigen3: Err=",EigErr
    print*," "
  END BLOCK
  
end program main
#endif
