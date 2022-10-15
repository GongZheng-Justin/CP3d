module m_Tools
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1
  implicit none
  private
  
  public:: CalcUxAver, Clc_Q_vor, Clc_lamda2
  public:: CalcMaxCFL, CheckDivergence, CalcVmax, CalcDissipationRate
  public:: InverseTridiagonal, InversePeriodicTridiagonal,InversePTriFixedCoe
#ifdef ScalarFlow
  public:: ClcMeanValue
#endif
contains    
   
  !******************************************************************
  ! CalcMaxCFL
  !******************************************************************     
  !  uddxmax = max{ |u1|/dx + |u2|/dy  + |u3|/dz } at cell center
  subroutine CalcMaxCFL(ux,uy,uz,uddxmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),intent(out):: uddxmax
    
    ! locals
    real(RK)::vcf
    integer::ic,jc,kc,ip,jp,kp,ierror
    
    uddxmax=zero
    do kc=y1start(3),y1end(3)
      kp=kc+1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        do ic=y1start(1),y1end(1)
          ip=ic+1
          vcf= abs(ux(ic,jc,kc)+ux(ip,jc,kc))*rdx+abs(uy(ic,jc,kc)+uy(ic,jp,kc))*rdyp(jc)+ &
               abs(uz(ic,jc,kc)+uz(ic,jc,kp))*rdz
          if(vcf>uddxmax)uddxmax=vcf
        enddo
      enddo
    enddo      
    vcf = uddxmax*half

    call MPI_ALLREDUCE(vcf,uddxmax,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
  end subroutine CalcMaxCFL
    
  !******************************************************************
  ! CheckDivergence
  !******************************************************************     
  subroutine  CheckDivergence(ux,uy,uz, divmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),intent(out)::divmax
    
    ! locals
    real(RK)::sucaj,udiv,divmax1,dudx,dvdy,dwdz
    integer::ic,jc,kc,ip,jp,kp,ierror,km,im,ku,iu

    !  ***** compute the divg(U)
    divmax1=zero
    do kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1;ku=kc+2
      do jc=y1start(2),y1end(2)
        jp=jc+1;sucaj=rdyp(jc)
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1;iu=ic+2
          dudx=  dxCoe1*ux(im,jc,kc) +dxCoe2*ux(ic,jc,kc) +dxCoe3*ux(ip,jc,kc) +dxCoe4*ux(iu,jc,kc)
          dvdy=  (uy(ic,jp,kc)-uy(ic,jc,kc))*sucaj
          dwdz=  dzCoe1*uz(ic,jc,km) +dzCoe2*uz(ic,jc,kc) +dzCoe3*uz(ic,jc,kp) +dzCoe4*uz(ic,jc,ku)

          udiv= abs(dudx+dvdy+dwdz)
          if(udiv>divmax1) divmax1=udiv
        enddo
      enddo
    enddo
    call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)
  end subroutine CheckDivergence
  
  !******************************************************************
  ! CalcVmax
  !******************************************************************  
  function CalcVmax(ux,uy,uz) result(vmaxabs)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(3)::vmaxabs
    
    ! locals
    integer::ic,jc,kc,ierror
    real(RK)::vfm,vmaxabs1(3)

    vmaxabs=zero; vmaxabs1=zero
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          vfm=abs(ux(ic,jc,kc))
          if(vfm>vmaxabs1(1))vmaxabs1(1)=vfm
          vfm=abs(uy(ic,jc,kc))
          if(vfm>vmaxabs1(2))vmaxabs1(2)=vfm
          vfm=abs(uz(ic,jc,kc))
          if(vfm>vmaxabs1(3))vmaxabs1(3)=vfm                
        enddo
      enddo
    enddo
    call MPI_REDUCE(vmaxabs1,vmaxabs,3,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)

  end function CalcVmax

  !******************************************************************
  ! InversePeriodicTridiagonal
  !******************************************************************
  subroutine InversePeriodicTridiagonal(aj,bj,cj,fj,m,n) ! my periodic tridiagonal solvers (variable coefficients )
    implicit none
    integer,intent(in)::m,n
    real(RK),dimension(m,n),intent(in):: aj,bj,cj
    real(RK),dimension(m,n),intent(inout)::fj

    ! locals
    integer::  i,j
    real(RK):: ppj,arrmn1(m,n),arrmn2(m,n),arrmn3(m,n)

    do i=1,m                                                     
      arrmn1(i,1)= -cj(i,1)/bj(i,1)                                        
      arrmn2(i,1)= -aj(i,1)/bj(i,1)                                     
      fj(i,1)= fj(i,1)/bj(i,1)                                         
    enddo                                                          
                                                                      
    ! forward elimination sweep                                                                      
    do j=2,n-1                                                     
      do i=1,m                                                     
        ppj =one/(bj(i,j)+ aj(i,j)*arrmn1(i,j-1))                              
        arrmn1(i,j) = -cj(i,j)*ppj                                            
        arrmn2(i,j) = -aj(i,j)*arrmn2(i,j-1)*ppj                                   
        fj(i,j) = (fj(i,j)-aj(i,j)*fj(i,j-1))*ppj                         
      enddo                                                         
    enddo                                         
                                                                      
    ! backward pass                                                                   
    do i=1,m                                                     
      arrmn2(i,n)= one                                                     
      arrmn3(i,n)= zero                                                    
    enddo                                                         
    do j=n-1,1,-1                                                                                                           
      do i=1,m                                                     
        arrmn2(i,j)= arrmn2(i,j) + arrmn1(i,j)*arrmn2(i,j+1)                                 
        arrmn3(i,j)= fj(i,j) + arrmn1(i,j)*arrmn3(i,j+1)                               
      enddo                                                          
    enddo                                                          
    do i=1,m                                                     
      fj(i,n)=(fj(i,n)-cj(i,n)*arrmn3(i,1)-aj(i,n)*arrmn3(i,n-1))/(cj(i,n)*arrmn2(i,1)+aj(i,n)*arrmn2(i,n-1)+bj(i,n))           
    enddo
                                                                    
    ! backward elimination pass                                                                     
    do j=n-1,1,-1                                                    
      do i=1,m                                                     
        fj(i,j)= fj(i,n)*arrmn2(i,j)+arrmn3(i,j)                                 
      enddo                                                         
    enddo
  end subroutine InversePeriodicTridiagonal

  !******************************************************************
  ! InversePTriFixedCoe
  !******************************************************************
  subroutine InversePTriFixedCoe(ajf,bjf,cjf,fj,m,n) ! my periodic tridiagonal solvers (fixedcoefficients )
    implicit none
    integer,intent(in)::m,n
    real(RK),intent(in):: ajf,bjf,cjf
    real(RK),dimension(m,n),intent(inout)::fj

    ! locals
    integer::  i,j
    real(RK),dimension(n):: ppj,vecn1,vecn2
    real(RK),dimension(m,n)::arrmn

    vecn1(1)= -cjf/bjf                                     
    vecn2(1)= -ajf/bjf
    vecn2(n)= one  
    do j=2,n-1                                                     
      ppj(j)   = one/(bjf+ ajf*vecn1(j-1))                              
      vecn1(j) = -cjf*ppj(j)                                            
      vecn2(j) = -ajf*vecn2(j-1)*ppj(j)                                                       
    enddo 
    do j=n-1,1,-1                                                     
      vecn2(j)= vecn2(j) + vecn1(j)*vecn2(j+1)                                                        
    enddo
         
    ! forward elimination sweep
    do i=1,m                     
      fj(i,1)= fj(i,1)/bjf                                         
    enddo                                                             
    do j=2,n-1                                                     
      do i=1,m                                   
        fj(i,j) = (fj(i,j)-ajf*fj(i,j-1))*ppj(j)                         
      enddo                                                         
    enddo                                         
                                                                      
    ! backward pass                                                 
    do i=1,m                                              
      arrmn(i,n)= zero                                                    
    enddo
    do j=n-1,1,-1                                                                                                           
      do i=1,m                                                                      
        arrmn(i,j)= fj(i,j) + vecn1(j)*arrmn(i,j+1)                               
      enddo                                                          
    enddo                                                             
    do i=1,m                                                     
      fj(i,n)=(fj(i,n)-cjf*arrmn(i,1)-ajf*arrmn(i,n-1))/(cjf*vecn2(1)+ajf*vecn2(n-1)+bjf)
    enddo
                                                                    
    ! backward elimination pass                                                                     
    do j=n-1,1,-1                                                    
      do i=1,m                                                     
        fj(i,j)= fj(i,n)*vecn2(j)+arrmn(i,j)                                
      enddo                                                         
    enddo
  end subroutine InversePTriFixedCoe

  !******************************************************************
  ! InverseTridiagonal
  !******************************************************************  
  subroutine InverseTridiagonal(aj,bj,cj,fj,m,n)   ! my tridiagonal solvers (variable coefficients )
    implicit none
    integer,intent(in)::m,n
    real(RK),dimension(m,n),intent(in):: aj,bj,cj
    real(RK),dimension(m,n),intent(inout)::fj
    
    ! locals
    integer:: i,j
    real(RK),dimension(m):: vecm
    real(RK),dimension(m,n)::arrmn
      
    do i=1,m
      vecm(i)=bj(i,1)
      fj(i,1)=fj(i,1)/vecm(i)
    enddo
    do j=2,n
      do i=1,m
        arrmn(i,j)=cj(i,j-1)/vecm(i)
        vecm(i)=bj(i,j)-aj(i,j)*arrmn(i,j)
        fj(i,j)=(fj(i,j)-aj(i,j)*fj(i,j-1))/vecm(i)
      enddo
    enddo
    do j=n-1,1,-1
      do i=1,m
        fj(i,j)=fj(i,j)-arrmn(i,j+1)*fj(i,j+1)
       enddo
    enddo    
  end subroutine InverseTridiagonal

  !******************************************************************
  ! CalcUxAver
  !******************************************************************
  function CalcUxAver(ux) result(res)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux
    real(RK):: res

    ! locals
    integer:: ic,jc,kc,ierror
    real(RK):: res1

    res1=zero
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          res1=res1+ux(ic,jc,kc)*dyp(jc)
        enddo
      enddo
    enddo
    call MPI_REDUCE(res1,res,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    res=res/(real(nxc*nzc,kind=RK))/yly
  end function CalcUxAver
  
#ifdef ScalarFlow
  !******************************************************************
  ! ClcMeanValue
  !******************************************************************
  function ClcMeanValue(ux,scalar) result(res)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,scalar
    real(RK)::res(3),res1(3),scValue,dypValue
    integer:: ic,jc,kc,ierror

    res1=zero
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        dypValue=dyp(jc)
        do ic=y1start(1),y1end(1)
          scValue=scalar(ic,jc,kc)
          res1(1)=res1(1)+dypValue*ux(ic,jc,kc)
          res1(2)=res1(2)+dypValue*scValue
          res1(3)=res1(3)+dypValue*scValue*scValue
        enddo
      enddo
    enddo
    call MPI_REDUCE(res1,res,3,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    res=res/(real(nxc*nzc,kind=RK))/yly
  end function ClcMeanValue
#endif
  
  !******************************************************************
  ! Clc_Q_vor
  !******************************************************************
  subroutine Clc_Q_vor(ux,uy,uz,Q_vor)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::Q_vor

    ! locals
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::udiv,caj,cac1,cac2,cac12
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    
    ! Q = 0.5*(ui,i *uj,j - ui,j *uj,i)
    ! Similar expression in subroutine "ClcVelStrain"
    do kc=y1start(3),y1end(3)
      kp=kc+1
      km=kc-1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        jm=jc-1
        caj  = rdyp(jc)
        cac1 = rdyc(jc)
        cac2 = rdyc(jp)    
        cac12= cac1 - cac2
        do ic=y1start(1),y1end(1)
          ip=ic+1
          im=ic-1
          dudx=  (ux(ip,jc,kc) -ux(ic,jc,kc))*rdx
          dudy= ((ux(ip,jp,kc) +ux(ic,jp,kc))*cac2    &
                +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12   &
                -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1    )    *quarter
          dudz=  (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz *quarter
      
          dvdx=  (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx *quarter
          dvdy=  (uy(ic,jp,kc) -uy(ic,jc,kc))*rdyp(jc)
          dvdz=  (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz *quarter

          dwdx=  (uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx *quarter
          dwdy= ((uz(ic,jp,kp) +uz(ic,jp,kc))*cac2    &
                +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12   &
                -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1    )    *quarter
          dwdz=  (uz(ic,jc,kp) -uz(ic,jc,kc))*rdz

          udiv= dudx +dvdy +dwdz 
          Q_vor(ic,jc,kc)= half*(udiv*udiv -dudx*dudx -dvdy*dvdy -dwdz*dwdz) -dudy*dvdx -dudz*dwdx -dvdz*dwdy
        enddo
      enddo
    enddo
  end subroutine Clc_Q_vor

  !******************************************************************
  ! Clc_lamda2
  !******************************************************************
  subroutine Clc_lamda2(ux,uy,uz,lamda2_vor)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::lamda2_vor

    ! locals
    logical:: AllRealFlag
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,rb,rc,rd
    real(RK)::Ome12,Ome31,Ome23,Eps11,Eps12,Eps13,Eps22,Eps23,Eps33
    real(RK)::caj,cac1,cac2,cac12,Lam11,Lam12,Lam13,Lam22,Lam23,Lam33,Root(3)
    
    ! Similar expression in subroutine "ClcVelStrain"

    ! Ome(i,j)= 0.5*(uj,i - ui,j)
    ! Eps(i,j)= 0.5*(uj,i + ui,j)
    ! Lam(i,j)= Ome(i,k) *Ome(k,j) +Eps(i,k) *Eps(k,j)
    ! lamda2 is the second eigvalue of tensor Lam.

    ! Lam12= Ome(1,k) *Ome(k,2) +Eps(1,k) *Eps(k,2)= Eps11*Eps12 +Eps12*Eps22 +Eps13*Eps23 +Ome31*Ome23
    ! Lam13= Ome(1,k) *Ome(k,3) +Eps(1,k) *Eps(k,3)= Eps11*Eps13 +Eps12*Eps23 +Eps13*Eps33 +Ome12*Ome23
    ! Lam23= Ome(2,k) *Ome(k,3) +Eps(2,k) *Eps(k,3)= Eps12*Eps13 +Eps22*Eps23 +Eps23*Eps33 +Ome12*Ome31
    ! Lam11= Ome(1,k) *Ome(k,1) +Eps(1,k) *Eps(k,1)= Eps11*Eps11 +Eps12*Eps12 +Eps13*Eps13 -Ome12*Ome12 -Ome31*Ome31
    ! Lam22= Ome(2,k) *Ome(k,2) +Eps(2,k) *Eps(k,2)= Eps12*Eps12 +Eps22*Eps22 +Eps23*Eps23 -Ome12*Ome12 -Ome23*Ome23
    ! Lam33= Ome(3,k) *Ome(k,3) +Eps(3,k) *Eps(k,3)= Eps13*Eps13 +Eps23*Eps23 +Eps33*Eps33 -Ome31*Ome31 -Ome23*Ome23

    do kc=y1start(3),y1end(3)
      kp=kc+1
      km=kc-1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        jm=jc-1
        caj  = rdyp(jc)
        cac1 = rdyc(jc)
        cac2 = rdyc(jp)    
        cac12= cac1 - cac2
        do ic=y1start(1),y1end(1)
          ip=ic+1
          im=ic-1
          dudx=  (ux(ip,jc,kc) -ux(ic,jc,kc))*rdx
          dudy= ((ux(ip,jp,kc) +ux(ic,jp,kc))*cac2    &
                +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12   &
                -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1    )    *quarter
          dudz=  (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz *quarter
      
          dvdx=  (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx *quarter
          dvdy=  (uy(ic,jp,kc) -uy(ic,jc,kc))*rdyp(jc)
          dvdz=  (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz *quarter

          dwdx=  (uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx *quarter
          dwdy= ((uz(ic,jp,kp) +uz(ic,jp,kc))*cac2    &
                +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12   &
                -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1    )    *quarter
          dwdz=  (uz(ic,jc,kp) -uz(ic,jc,kc))*rdz

          Ome12= dvdx - dudy
          Ome23= dwdy - dvdz
          Ome31= dudz - dwdx
          Eps11= dudx + dudx
          Eps12= dudy + dvdx
          Eps13= dudz + dwdx
          Eps22= dvdy + dvdy
          Eps23= dvdz + dwdy
          Eps33= dwdz + dwdz
          Lam12= Eps11*Eps12 +Eps12*Eps22 +Eps13*Eps23 +Ome31*Ome23
          Lam13= Eps11*Eps13 +Eps12*Eps23 +Eps13*Eps33 +Ome12*Ome23
          Lam23= Eps12*Eps13 +Eps22*Eps23 +Eps23*Eps33 +Ome12*Ome31
          Lam11= Eps11*Eps11 +Eps12*Eps12 +Eps13*Eps13 -Ome12*Ome12 -Ome31*Ome31
          Lam22= Eps12*Eps12 +Eps22*Eps22 +Eps23*Eps23 -Ome12*Ome12 -Ome23*Ome23
          Lam33= Eps13*Eps13 +Eps23*Eps23 +Eps33*Eps33 -Ome31*Ome31 -Ome23*Ome23

          rb= -(Lam11+Lam22+Lam33)
          rc= Lam11*Lam22 +Lam22*Lam33 +Lam33*Lam11 -Lam12*Lam12 -Lam13*Lam13 -Lam23*Lam23
          rd= Lam12*Lam12*Lam33 +Lam13*Lam13*Lam22 +Lam23*Lam23*Lam11 -Lam11*Lam22*Lam33 -2.0_RK*Lam12*Lam23*Lam13

          call CubicRoot(rb,rc,rd,Root,AllRealFlag)
          if(AllRealFlag) then
            lamda2_vor(ic,jc,kc)=quarter*Root(2)
          else
            call MainLog%CheckForError(ErrT_Abort,"Clc_lamda2","lamda2 wrong")
          endif
        enddo
      enddo
    enddo
  end subroutine Clc_lamda2

  !******************************************************************
  ! CalcDissipationRate
  !******************************************************************
  subroutine CalcDissipationRate(ux,uy,uz,dissp)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::dissp

    ! locals
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::st1,st2,st3,st4,st5,st6,caj,cac1,cac2,cac12

    ! epsilon_ij= (2*xnu*<S_ij*S_ij>)/xnu, where S_ij = 0.5*(dui/dxj + duj/dxi)
    ! Similar expression in subroutine "ClcVelStrain"
    DO kc=y1start(3),y1end(3)
      kp=kc+1
      km=kc-1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        jm=jc-1
        caj  = rdyp(jc)
        cac1 = rdyc(jc)
        cac2 = rdyc(jp)    
        cac12= cac1 - cac2
        do ic=y1start(1),y1end(1)
          ip=ic+1
          im=ic-1
          st1=(ux(ip,jc,kc)-ux(ic,jc,kc))*rdx
          st2=(uy(ic,jp,kc)-uy(ic,jc,kc))*caj
          st3=(uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
          st4=oneEighth*( (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx  &
                         +(ux(ip,jp,kc) +ux(ic,jp,kc))*cac2                             &
                         +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12                            &
                         -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1                             )
          st5=oneEighth*( (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz  &  
                         +(uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx  )
          st6=oneEighth*( (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz  &
                         +(uz(ic,jp,kp) +uz(ic,jp,kc))*cac2                             &
                         +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12                            &
                         -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1                             )
          dissp(ic,jc,kc)= two*(st1*st1+ st2*st2+ st3*st3)+ four*(st4*st4 + st5*st5 + st6*st6)
        enddo
      enddo
    ENDDO

  end subroutine CalcDissipationRate

  !******************************************************************
  ! CubicRoot
  !******************************************************************
  subroutine CubicRoot(b,c,d,Root,AllRealFlag)
    implicit none
    real(RK),intent(in):: b,c,d
    logical,intent(out):: AllRealFlag
    real(RK),dimension(3),intent(out)::Root

    ! locals
    real(RK):: SMALL=1.0E-12_RK
    real(RK):: p,q,t,rho,eta,root1,root2,root3
    real(RK):: OneThird=0.3333333333333333333333_RK

    p = c-OneThird*b*b
    q = (2.0_RK*b*b*b-9.0_RK*b*c +27.0_RK*d)/27.0_RK
    t = q*q/4.0_RK+p*p*p/27.0_RK

    !  as q^2 / 4 + p^3/27 < 0, p < 0, -p and rho > 0
    if(t<0.0_RK) then
      AllRealFlag=.true.
      t  = -OneThird*b
      rho= ((-p)**1.5_RK)/sqrt(27.0_RK)
      eta= acos(0.5_RK*q/rho)
      rho= -2.0_RK*(rho**OneThird)
      root1= cos(OneThird*eta)*rho +t
      root2= cos(OneThird*(2.0_RK*Pi +eta))*rho +t
      root3= cos(OneThird*(4.0_RK*Pi +eta))*rho +t
      if(root1<root2 .and. root1<root3) then
        if(root2<root3) then
          Root=(/root1,root2,root3/);return
        else
          Root=(/root1,root3,root2/);return
        endif
      elseif(root1>root2 .and. root1>root3) then
        if(root2<root3) then
          Root=(/root2,root3,root1/);return
        else
          Root=(/root3,root2,root1/);return
        endif
      else
        if(root2<root3) then
          Root=(/root2,root1,root3/);return
        else
          Root=(/root3,root1,root2/);return
        endif
      endif

    elseif(t<SMALL) then
      AllRealFlag=.true.  
      rho= 0.5_RK*q
      if(rho<0.0_RK) then
        rho= ((-rho)**OneThird)
      else
        rho= -rho**OneThird
      endif
      Root(1)= 2.0_RK*rho-OneThird*b
      Root(2)= -rho-OneThird*b
      Root(3)= Root(2)

    else
      AllRealFlag=.false.
      t= sqrt(t)
      rho= 0.5_RK*q+t
      eta= 0.5_RK*q-t
      if(rho<0.0_RK) then
        rho= ((-rho)**OneThird)
      else
        rho= -rho**OneThird
      endif
      if(eta<0.0_RK) then
        eta= ((-eta)**OneThird)
      else
        eta= -eta**OneThird
      endif
      Root(1)= rho + eta-OneThird*b
      Root(2)= -0.5_RK*(rho+eta)-OneThird*b
      Root(3)= (0.5_RK*sqrt(3.0))*(eta-rho)
    endif
  end subroutine CubicRoot
    
end module m_Tools
