!================================== Pressure Gradient force ==================================!
subroutine clc_PrGradForce(pressure,RatioYp_interp,RatioYc_interp)
  implicit none
  real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
  type(real3),dimension(GPrtcl_list%mlocal),intent(in)::RatioYp_interp,RatioYc_interp

  ! locals
  type(real4)::pos
  type(real3)::PrGradForce
  type(HaloInfo):: hi_px_interp,hi_pz_interp
  integer::ic,jc,kc,im,jm,km,ip,jp,kp,id,jd,kd,pid,itype,i,j,k,nlocal
  integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
  real(RK)::prx,pry,prz,SumXDir,SumYDir,SumZDir,sucac
  real(RK),dimension(0:2)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
  real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm):: dpdx,dpdy,dpdz

  nlocal=GPrtcl_list%nlocal   
  dpdx=0.0_RK; dpdy=0.0_RK; dpdz=0.0_RK

  ! dp/dx
  DO kc=y1start(3)-1,y1end(3)+1
    do jc=0,nyp
      do ic=y1start(1),y1end(1)+1
        im=ic-1
        dpdx(ic,jc,kc)= (pressure(ic,jc,kc)-pressure(im,jc,kc))*rdx
      enddo
    enddo
  ENDDO

  hi_px_interp%pencil = y_pencil
  hi_px_interp%xmh=1; hi_px_interp%xph=2
  hi_px_interp%ymh=0; hi_px_interp%yph=0
  hi_px_interp%zmh=0; hi_px_interp%zph=0
  call myupdate_halo(dpdx, mb1, hi_px_interp)

  ! dp/dy
  DO kc=y1start(3)-1,y1end(3)+1
    do jc=1,nyp
      jm=jc-1
      sucac= rdyc(jc)
      do ic=y1start(1)-1,y1end(1)+1
        dpdy(ic,jc,kc)= (pressure(ic,jc,kc)-pressure(ic,jm,kc))*sucac
      enddo
    enddo
  ENDDO

  ! dp/dz
  DO kc=y1start(3),y1end(3)+1
    km=kc-1
    do jc=0,nyp
      do ic=y1start(1)-1,y1end(1)+1
        dpdz(ic,jc,kc)= (pressure(ic,jc,kc)-pressure(ic,jc,km))*rdz
      enddo
    enddo
  ENDDO

  hi_pz_interp%pencil = y_pencil
  hi_pz_interp%xmh=0; hi_pz_interp%xph=0
  hi_pz_interp%ymh=0; hi_pz_interp%yph=0
  hi_pz_interp%zmh=1; hi_pz_interp%zph=2
  call myupdate_halo(dpdz, mb1, hi_pz_interp)

  ! interpolation to get Pressure Gradient at particle position
  DO pid=1,nlocal
    pos = GPrtcl_posR(pid) 
    itype= GPrtcl_pType(pid)
    SumXDir=0.0_RK;  SumYDir=0.0_RK;  SumZDir=0.0_RK

    idxc_interp=indxyz(1,pid)
    idxp_interp=indxyz(2,pid)
    idyc_interp=indxyz(3,pid)
    idyp_interp=indxyz(4,pid)
    idzc_interp=indxyz(5,pid)
    idzp_interp=indxyz(6,pid)

    prx=(pos%x-xc(idxp_interp))*rdx+0.5_RK
    RatioXp(0)=0.5_RK*(prx-1.0_RK)*(prx-2.0_RK)
    RatioXp(1)=      prx*(2.0_RK-prx)
    RatioXp(2)=0.5_RK* prx*(prx-1.0_RK)
    prx=(pos%x-xc(idxc_interp))*rdx
    RatioXc(0)=0.5_RK*(prx-1.0_RK)*(prx-2.0_RK)
    RatioXc(1)=      prx*(2.0_RK-prx)
    RatioXc(2)=0.5_RK* prx*(prx-1.0_RK)          

    RatioYp(0) = RatioYp_interp(pid)%x
    RatioYp(1) = RatioYp_interp(pid)%y
    RatioYp(2) = RatioYp_interp(pid)%z
    RatioYc(0) = RatioYc_interp(pid)%x
    RatioYc(1) = RatioYc_interp(pid)%y
    RatioYc(2) = RatioYc_interp(pid)%z

    prz=(pos%z-zc(idzp_interp))*rdz+0.5_RK
    RatioZp(0)=0.5_RK*(prz-1.0_RK)*(prz-2.0_RK)
    RatioZp(1)=      prz*(2.0_RK-prz)
    RatioZp(2)=0.5_RK* prz*(prz-1.0_RK) 
    prz=(pos%z-zc(idzc_interp))*rdz
    RatioZc(0)=0.5_RK*(prz-1.0_RK)*(prz-2.0_RK)
    RatioZc(1)=      prz*(2.0_RK-prz)
    RatioZc(2)=0.5_RK* prz*(prz-1.0_RK)

    ! dp/dx grid, the same to ux grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd  = j+idyc_interp
        pry = RatioYc(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumXDir= SumXDir + dpdx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dp/dy grid, the same to uy grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumYDir= SumYDir + dpdy(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dp/dz grid, the same to uz grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumZDir= SumZDir + dpdz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    PrGradForce= (-DEMProperty%Prtcl_PureProp(itype)%MassOfFluid)* real3(SumXDir,SumYDir,SumZDir)
    GPrtcl_FpForce(pid)=GPrtcl_FpForce(pid)+ PrGradForce
  ENDDO

  ! fixed particle part
  DO pid=1,mlocalFix
    pos  = GPFix_posR(pid) 
    itype= GPFix_pType(pid)
    SumXDir=0.0_RK;  SumYDir=0.0_RK;  SumZDir=0.0_RK

    idxc_interp=indxyzFixed(1,pid)
    idxp_interp=indxyzFixed(2,pid)
    idyc_interp=indxyzFixed(3,pid)
    idyp_interp=indxyzFixed(4,pid)
    idzc_interp=indxyzFixed(5,pid)
    idzp_interp=indxyzFixed(6,pid)

    prx=(pos%x-xc(idxp_interp))*rdx+0.5_RK
    RatioXp(0)=0.5_RK*(prx-1.0_RK)*(prx-2.0_RK)
    RatioXp(1)=      prx*(2.0_RK-prx)
    RatioXp(2)=0.5_RK* prx*(prx-1.0_RK)
    prx=(pos%x-xc(idxc_interp))*rdx
    RatioXc(0)=0.5_RK*(prx-1.0_RK)*(prx-2.0_RK)
    RatioXc(1)=      prx*(2.0_RK-prx)
    RatioXc(2)=0.5_RK* prx*(prx-1.0_RK)          

    RatioYp(0) = RYpFixed_interp(pid)%x
    RatioYp(1) = RYpFixed_interp(pid)%y
    RatioYp(2) = RYpFixed_interp(pid)%z
    RatioYc(0) = RYcFixed_interp(pid)%x
    RatioYc(1) = RYcFixed_interp(pid)%y
    RatioYc(2) = RYcFixed_interp(pid)%z

    prz=(pos%z-zc(idzp_interp))*rdz+0.5_RK
    RatioZp(0)=0.5_RK*(prz-1.0_RK)*(prz-2.0_RK)
    RatioZp(1)=      prz*(2.0_RK-prz)
    RatioZp(2)=0.5_RK* prz*(prz-1.0_RK) 
    prz=(pos%z-zc(idzc_interp))*rdz
    RatioZc(0)=0.5_RK*(prz-1.0_RK)*(prz-2.0_RK)
    RatioZc(1)=      prz*(2.0_RK-prz)
    RatioZc(2)=0.5_RK* prz*(prz-1.0_RK)

    ! dp/dx grid, the same to ux grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd  = j+idyc_interp
        pry = RatioYc(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumXDir= SumXDir + dpdx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dp/dy grid, the same to uy grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumYDir= SumYDir + dpdy(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dp/dz grid, the same to uz grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumZDir= SumZDir + dpdz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    PrGradForce= (-DEMProperty%Prtcl_PureProp(itype)%MassOfFluid)* real3(SumXDir,SumYDir,SumZDir)
    GPFix_FpForce(pid)=GPFix_FpForce(pid)+ PrGradForce
  ENDDO

end subroutine clc_PrGradForce
