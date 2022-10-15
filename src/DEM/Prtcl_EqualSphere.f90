!********************************************************************
! 
!  Purpose::
!    1) Provide a much cut down module to performing partition of the   
!       surface of a sphere into regions of equal area.
! 
!  References:
!    1) Leopardi P. A partition of the unit sphere into regions of 
!       equal area and small diameter[J]. Electronic Transactions on 
!       Numerical Analysis Etna, 2006, 25(1):2006.
!    2) http://eqsp.sourceforge.net/
!    3) https://sourceforge.net/projects/eqsp/
! 
!  Author: Zheng Gong
!  Date:   2021-08-13
! 
!********************************************************************
module Prtcl_EqualSphere
  implicit none
  private
  integer,parameter::RK=8
  public::eq_sphere
contains
  !******************************************************************
  ! eq_sphere
  !******************************************************************
  subroutine eq_sphere(Point)
    implicit none
    real(RK),dimension(:,:),intent(out)::Point

    ! locals
    integer,dimension(:),allocatable::n_regions
    integer::k,m,nMarker,nColumn,n_collars,nCount,n_top,n_bot
    real(RK),parameter::PI=3.141592653589793238462643383279502884_RK
    real(RK)::rx,ry,rz,rnorm,r_regions,a_top,a_bot,area_tot,Psi,Phi,aTemp
    real(RK)::area_cap,c_polar,discrepancy,a_fitting,area_top,area_bot,offset

    nMarker=size(Point,1)
    nColumn=size(Point,2)
    if(nMarker<1 .or. nColumn/=3) then
      print*,"Error in eq_sphere, nMarker or nColumn wrong:",nMarker,nColumn; stop
    endif
    Point(1,:)=(/0.0_RK, 0.0_RK, 1.0_RK/);       if(nMarker==1)return
    Point(2,:)=(/1.0_RK, 0.0_RK, 0.0_RK/); 
    Point(nMarker,:)=(/0.0_RK, 0.0_RK,-1.0_RK/); if(nMarker< 4)return

    area_cap=4.0_RK*PI/real(nMarker,RK)
    c_polar =2.0_RK*asin(sqrt(area_cap/(4.0_RK*PI)))
    n_collars= max(1,nint((PI-2.0_RK*c_polar)/sqrt(area_cap)))
    allocate(n_regions(n_collars+2))
    n_regions(1)=1;n_regions(n_collars+2)=1
    discrepancy = 0.0_RK; area_top=area_cap
    a_fitting = (PI-2.0_RK*c_polar)/real(n_collars,RK)
    do k=1,n_collars
      area_bot=c_polar+real(k,RK)*a_fitting
      area_bot=sin(0.5_RK*area_bot)
      area_bot=4*PI*area_bot*area_bot
      r_regions=(area_bot-area_top)/area_cap
      area_top=area_bot
      n_regions(k+1)=nint(r_regions+discrepancy)
      discrepancy=discrepancy+r_regions-n_regions(k+1)
    enddo

    nCount= 2; offset= 0.0_RK
    a_top = c_polar; area_tot=area_cap
    do k=1,n_collars
      n_top=n_regions(k+1)
      n_bot=n_regions(k+2)
      area_tot=area_tot+n_top*area_cap
      a_bot= 2.0_RK*asin(sqrt(area_tot/(4.0_RK*PI)))
      Psi= 0.5_RK*(a_top+a_bot); a_top=a_bot
      do m=1,n_top
        aTemp=real(2*m-1,RK)*PI/real(n_top,RK) +2.0_RK*PI*offset
        Phi = aTemp -2.0*PI*floor(aTemp/(2.0*PI))
        rx=sin(Psi)*cos(Phi)
        ry=sin(Psi)*sin(Phi)
        rz=cos(Psi)
        rnorm=sqrt(rx*rx+ry*ry+rz*rz)
        Point(nCount,:)=(/rx,ry,rz/)/rnorm
        nCount=nCount+1
      enddo
      offset=offset+real(n_top-n_bot+gcd(n_top,n_bot),RK)/real(2*n_top*n_bot,RK)
      offset=offset-floor(offset)
    enddo
    if(nMarker /= nCount) then
      print*,'Error in eq_sphere: nMarker /= nCount',nMarker,nCount; stop
    endif
  end subroutine eq_sphere

  !******************************************************************
  ! gcd: Greatest Common Divisor
  !******************************************************************
  function gcd(m1,n1) result(r)
    implicit none
    integer,intent(in)::m1,n1
    integer::m,n,r
    m=m1; n=n1
    do 
      r=mod(m,n); m=n; n=r
      if(r==0)exit
    enddo
    r=m
  end function gcd
end module Prtcl_EqualSphere

#ifdef test_EqualSphere
!******************************************************************
! Main program
!******************************************************************
program main
  use Prtcl_EqualSphere
  implicit none
  integer::k,nMarker
  real(8),dimension(:,:),allocatable::Point

  print*,'nMarker='
  read*,nMarker
  allocate(Point(nMarker,3))
  call eq_sphere(Point)
  print*,'    nid,      x,      y,      z'
  do k=1,nMarker
    print*,k,Point(k,:)
  enddo
end program main
#endif
