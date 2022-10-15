!===================================== Fluid Acceletare force  =====================================!

subroutine clc_FluidAccForce(ux,uy,uz,RatioYp_interp,RatioYc_interp)
  implicit none
  real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
  type(real3),dimension(GPrtcl_list%mlocal),intent(in)::RatioYp_interp,RatioYc_interp

  ! locals
  type(real4)::pos
  type(real3)::FpVel,Vfluid,VfluidOld,FluidAccForce
  integer::ic,jc,kc,im,jm,km,ip,jp,kp,id,jd,kd,pid,itype,i,j,k,nlocal
  integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
  real(RK)::prx,pry,prz,SumXDir,SumYDir,SumZDir,sucac,sucaj,massFeff
  real(RK),dimension(0:2)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
  type(HaloInfo):: hi_dudz_interp, hi_dvdx_interp, hi_dvdz_interp, hi_dwdx_interp
  real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm):: Mat1,Mat2,Mat3

  if(.not.IsHaveVfluidOld) then
   IsHaveVfluidOld= .true.
   return
  endif
  nlocal=GPrtcl_list%nlocal
  DO pid=1,nlocal
    itype    = GPrtcl_pType(pid)
    Vfluid   = GPrtcl_Vfluid(1,pid)
    VfluidOld= GPrtcl_Vfluid(2,pid)
    massFeff = FluidAccCoe*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid
    FluidAccForce      = massFeff*(Vfluid-VfluidOld)/dt
    GPrtcl_FpForce(pid)= GPrtcl_FpForce(pid)+ FluidAccForce
  ENDDO

  ! Fixed particle part
  DO pid=1,mlocalFix
    itype    = GPFix_pType(pid)
    Vfluid   = GPFix_Vfluid(1,pid)
    VfluidOld= GPFix_Vfluid(2,pid)
    massFeff = FluidAccCoe*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid
    FluidAccForce     = massFeff*(Vfluid-VfluidOld)/dt 
    GPFix_FpForce(pid)= GPFix_FpForce(pid)+ FluidAccForce
  ENDDO

  !================================== advection of ux
  asso_FluidAcc1: associate(dudx =>Mat1, dudy =>Mat2, dudz =>Mat3)
  dudx=zero; dudy=zero; dudz=zero

  ! dudx
  DO kc=y1start(3)-1,y1end(3)+1
    do jc=0,nyp
      do ic=y1start(1)-1,y1end(1)+1
        ip=ic+1
        dudx(ic,jc,kc)=(ux(ip,jc,kc)-ux(ic,jc,kc))*rdx
      enddo
    enddo
  ENDDO

  ! dudy
  DO kc=y1start(3)-1,y1end(3)+1
    do jc=1,nyp
      jm=jc-1
      sucac=rdyc(jc)
      do ic=y1start(1)-1,y1end(1)+2
        dudy(ic,jc,kc)=(ux(ic,jc,kc)-ux(ic,jm,kc))*sucac
      enddo
    enddo
  ENDDO

  ! dudz
  DO kc=y1start(3),y1end(3)+1
    km=kc-1
    do jc=0,nyp
      do ic=y1start(1)-1,y1end(1)+2
        dudz(ic,jc,kc)=(ux(ic,jc,kc)-ux(ic,jc,km))*rdz
      enddo
    enddo
  ENDDO
  hi_dudz_interp%pencil = y_pencil
  hi_dudz_interp%xmh=0; hi_dudz_interp%xph=0
  hi_dudz_interp%ymh=0; hi_dudz_interp%yph=0
  hi_dudz_interp%zmh=1; hi_dudz_interp%zph=2
  call myupdate_halo(dudz, mb1, hi_dudz_interp)

  DO pid=1,nlocal
    pos    = GPrtcl_posR(pid)
    itype  = GPrtcl_pType(pid)
    FpVel  = GPrtcl_Vfluid(1,pid)- GPrtcl_linVel(1,pid)
    SumXDir= zero;  SumYDir=zero;  SumZDir=zero
    massFeff = FluidAccCoe*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid

    idxc_interp=indxyz(1,pid)
    idxp_interp=indxyz(2,pid)
    idyc_interp=indxyz(3,pid)
    idyp_interp=indxyz(4,pid)
    idzc_interp=indxyz(5,pid)
    idzp_interp=indxyz(6,pid)

    prx=(pos%x-xc(idxp_interp))*rdx+half
    RatioXp(0)=half*(prx-one)*(prx-two)
    RatioXp(1)=      prx*(two-prx)
    RatioXp(2)=half* prx*(prx-one)
    prx=(pos%x-xc(idxc_interp))*rdx
    RatioXc(0)=half*(prx-one)*(prx-two)
    RatioXc(1)=      prx*(two-prx)
    RatioXc(2)=half* prx*(prx-one)          

    RatioYp(0) = RatioYp_interp(pid)%x
    RatioYp(1) = RatioYp_interp(pid)%y
    RatioYp(2) = RatioYp_interp(pid)%z
    RatioYc(0) = RatioYc_interp(pid)%x
    RatioYc(1) = RatioYc_interp(pid)%y
    RatioYc(2) = RatioYc_interp(pid)%z

    prz=(pos%z-zc(idzp_interp))*rdz+half
    RatioZp(0)=half*(prz-one)*(prz-two)
    RatioZp(1)=      prz*(two-prz)
    RatioZp(2)=half* prz*(prz-one) 
    prz=(pos%z-zc(idzc_interp))*rdz
    RatioZc(0)=half*(prz-one)*(prz-two)
    RatioZc(1)=      prz*(two-prz)
    RatioZc(2)=half* prz*(prz-one)

    ! dudx grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd  = j+idyc_interp
        pry = RatioYc(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumXDir= SumXDir + dudx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dudy gird
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumYDir= SumYDir + dudy(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dudz grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumZDir= SumZDir + dudz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    GPrtcl_FpForce(pid)%x = GPrtcl_FpForce(pid)%x + massFeff*(SumXDir*FpVel%x +SumYDir*FpVel%y +SumZDir*FpVel%z)
  ENDDO

  ! Fixed particle part
  DO pid=1,mlocalFix
    pos    = GPFix_posR(pid)
    itype  = GPFix_pType(pid)
    FpVel  = GPFix_Vfluid(1,pid)
    SumXDir= zero;  SumYDir=zero;  SumZDir=zero
    massFeff = FluidAccCoe*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid

    idxc_interp=indxyzFixed(1,pid)
    idxp_interp=indxyzFixed(2,pid)
    idyc_interp=indxyzFixed(3,pid)
    idyp_interp=indxyzFixed(4,pid)
    idzc_interp=indxyzFixed(5,pid)
    idzp_interp=indxyzFixed(6,pid)

    prx=(pos%x-xc(idxp_interp))*rdx+half
    RatioXp(0)=half*(prx-one)*(prx-two)
    RatioXp(1)=      prx*(two-prx)
    RatioXp(2)=half* prx*(prx-one)
    prx=(pos%x-xc(idxc_interp))*rdx
    RatioXc(0)=half*(prx-one)*(prx-two)
    RatioXc(1)=      prx*(two-prx)
    RatioXc(2)=half* prx*(prx-one)          

    RatioYp(0) = RYpFixed_interp(pid)%x
    RatioYp(1) = RYpFixed_interp(pid)%y
    RatioYp(2) = RYpFixed_interp(pid)%z
    RatioYc(0) = RYcFixed_interp(pid)%x
    RatioYc(1) = RYcFixed_interp(pid)%y
    RatioYc(2) = RYcFixed_interp(pid)%z

    prz=(pos%z-zc(idzp_interp))*rdz+half
    RatioZp(0)=half*(prz-one)*(prz-two)
    RatioZp(1)=      prz*(two-prz)
    RatioZp(2)=half* prz*(prz-one) 
    prz=(pos%z-zc(idzc_interp))*rdz
    RatioZc(0)=half*(prz-one)*(prz-two)
    RatioZc(1)=      prz*(two-prz)
    RatioZc(2)=half* prz*(prz-one)

    ! dudx grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd  = j+idyc_interp
        pry = RatioYc(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumXDir= SumXDir + dudx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dudy gird
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumYDir= SumYDir + dudy(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dudz grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumZDir= SumZDir + dudz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    GPFix_FpForce(pid)%x = GPFix_FpForce(pid)%x + massFeff*(SumXDir*FpVel%x +SumYDir*FpVel%y +SumZDir*FpVel%z)
  ENDDO
  end associate asso_FluidAcc1

  !================================== advection of uy
  asso_FluidAcc2: associate(dvdx =>Mat1, dvdy =>Mat2, dvdz =>Mat3)
  dvdx=zero; dvdy=zero; dvdz=zero

  ! dvdx
  DO kc=y1start(3)-1,y1end(3)+1
    do jc=1,nyp
      do ic=y1start(1),y1end(1)+1
        im=ic-1
        dvdx(ic,jc,kc)=(uy(ic,jc,kc)-uy(im,jc,kc))*rdx
      enddo
    enddo
  ENDDO
  hi_dvdx_interp%pencil = y_pencil
  hi_dvdx_interp%xmh=1; hi_dvdx_interp%xph=2
  hi_dvdx_interp%ymh=0; hi_dvdx_interp%yph=0
  hi_dvdx_interp%zmh=0; hi_dvdx_interp%zph=0
  call myupdate_halo(dvdx, mb1, hi_dvdx_interp)

  ! dvdy
  DO kc=y1start(3)-1,y1end(3)+1
    do jc=0,nyp
      jp=jc+1; sucaj=rdyp(jc)
      do ic=y1start(1)-1,y1end(1)+1
        dvdy(ic,jc,kc)=(uy(ic,jp,kc)-uy(ic,jc,kc))*sucaj
      enddo
    enddo
  ENDDO

  ! dvdz
  DO kc=y1start(3),y1end(3)+1
    km=kc-1
    do jc=1,nyp
      do ic=y1start(1)-1,y1end(1)+1
        dvdz(ic,jc,kc)=(uy(ic,jc,kc)-uy(ic,jc,km))*rdz
      enddo
    enddo
  ENDDO 
  hi_dvdz_interp%pencil = y_pencil
  hi_dvdz_interp%xmh=0; hi_dvdz_interp%xph=0
  hi_dvdz_interp%ymh=0; hi_dvdz_interp%yph=0
  hi_dvdz_interp%zmh=1; hi_dvdz_interp%zph=2
  call myupdate_halo(dvdz, mb1, hi_dvdz_interp)
  
  DO pid=1,nlocal
    pos    = GPrtcl_posR(pid)
    itype  = GPrtcl_pType(pid)
    FpVel  = GPrtcl_Vfluid(1,pid)- GPrtcl_linVel(1,pid)
    SumXDir= zero;  SumYDir=zero;  SumZDir=zero
    massFeff = FluidAccCoe*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid

    idxc_interp=indxyz(1,pid)
    idxp_interp=indxyz(2,pid)
    idyc_interp=indxyz(3,pid)
    idyp_interp=indxyz(4,pid)
    idzc_interp=indxyz(5,pid)
    idzp_interp=indxyz(6,pid)

    prx=(pos%x-xc(idxp_interp))*rdx+half
    RatioXp(0)=half*(prx-one)*(prx-two)
    RatioXp(1)=      prx*(two-prx)
    RatioXp(2)=half* prx*(prx-one)
    prx=(pos%x-xc(idxc_interp))*rdx
    RatioXc(0)=half*(prx-one)*(prx-two)
    RatioXc(1)=      prx*(two-prx)
    RatioXc(2)=half* prx*(prx-one)          

    RatioYp(0) = RatioYp_interp(pid)%x
    RatioYp(1) = RatioYp_interp(pid)%y
    RatioYp(2) = RatioYp_interp(pid)%z
    RatioYc(0) = RatioYc_interp(pid)%x
    RatioYc(1) = RatioYc_interp(pid)%y
    RatioYc(2) = RatioYc_interp(pid)%z

    prz=(pos%z-zc(idzp_interp))*rdz+half
    RatioZp(0)=half*(prz-one)*(prz-two)
    RatioZp(1)=      prz*(two-prz)
    RatioZp(2)=half* prz*(prz-one) 
    prz=(pos%z-zc(idzc_interp))*rdz
    RatioZc(0)=half*(prz-one)*(prz-two)
    RatioZc(1)=      prz*(two-prz)
    RatioZc(2)=half* prz*(prz-one)

    ! dvdx grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd  = j+idyp_interp
        pry = RatioYp(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumXDir= SumXDir + dvdx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dvdy gird
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumYDir= SumYDir + dvdy(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dvdz grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumZDir= SumZDir + dvdz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    GPrtcl_FpForce(pid)%y = GPrtcl_FpForce(pid)%y + massFeff*(SumXDir*FpVel%x +SumYDir*FpVel%y +SumZDir*FpVel%z)
  ENDDO

  ! Fixed particle part
  DO pid=1,mlocalFix
    pos    = GPFix_posR(pid)
    itype  = GPFix_pType(pid)
    FpVel  = GPFix_Vfluid(1,pid)
    SumXDir= zero;  SumYDir=zero;  SumZDir=zero
    massFeff = FluidAccCoe*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid

    idxc_interp=indxyzFixed(1,pid)
    idxp_interp=indxyzFixed(2,pid)
    idyc_interp=indxyzFixed(3,pid)
    idyp_interp=indxyzFixed(4,pid)
    idzc_interp=indxyzFixed(5,pid)
    idzp_interp=indxyzFixed(6,pid)

    prx=(pos%x-xc(idxp_interp))*rdx+half
    RatioXp(0)=half*(prx-one)*(prx-two)
    RatioXp(1)=      prx*(two-prx)
    RatioXp(2)=half* prx*(prx-one)
    prx=(pos%x-xc(idxc_interp))*rdx
    RatioXc(0)=half*(prx-one)*(prx-two)
    RatioXc(1)=      prx*(two-prx)
    RatioXc(2)=half* prx*(prx-one)          

    RatioYp(0) = RYpFixed_interp(pid)%x
    RatioYp(1) = RYpFixed_interp(pid)%y
    RatioYp(2) = RYpFixed_interp(pid)%z
    RatioYc(0) = RYcFixed_interp(pid)%x
    RatioYc(1) = RYcFixed_interp(pid)%y
    RatioYc(2) = RYcFixed_interp(pid)%z

    prz=(pos%z-zc(idzp_interp))*rdz+half
    RatioZp(0)=half*(prz-one)*(prz-two)
    RatioZp(1)=      prz*(two-prz)
    RatioZp(2)=half* prz*(prz-one) 
    prz=(pos%z-zc(idzc_interp))*rdz
    RatioZc(0)=half*(prz-one)*(prz-two)
    RatioZc(1)=      prz*(two-prz)
    RatioZc(2)=half* prz*(prz-one)

    ! dvdx grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd  = j+idyp_interp
        pry = RatioYp(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumXDir= SumXDir + dvdx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dvdy gird
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumYDir= SumYDir + dvdy(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dvdz grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumZDir= SumZDir + dvdz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    GPFix_FpForce(pid)%y = GPFix_FpForce(pid)%y + massFeff*(SumXDir*FpVel%x +SumYDir*FpVel%y +SumZDir*FpVel%z)
  ENDDO
  end associate asso_FluidAcc2

  !================================== advection of uz
  asso_FluidAcc3: associate(dwdx =>Mat1, dwdy =>Mat2, dwdz =>Mat3)
  dwdx=zero; dwdy=zero; dwdz=zero

  ! dwdx
  DO kc=y1start(3)-1,y1end(3)+2
    do jc=0,nyp
      do ic=y1start(1),y1end(1)+1
        im=ic-1
        dwdx(ic,jc,kc)=(uz(ic,jc,kc)-uz(im,jc,kc))*rdx
      enddo
    enddo
  ENDDO
  hi_dwdx_interp%pencil = y_pencil
  hi_dwdx_interp%xmh=1; hi_dwdx_interp%xph=2
  hi_dwdx_interp%ymh=0; hi_dwdx_interp%yph=0
  hi_dwdx_interp%zmh=0; hi_dwdx_interp%zph=0
  call myupdate_halo(dwdx, mb1, hi_dwdx_interp)

  ! dwdy
  DO kc=y1start(3)-1,y1end(3)+2
    do jc=1,nyp
      jm=jc-1
      sucac=rdyc(jc)
      do ic=y1start(1)-1,y1end(1)+1
        dwdy(ic,jc,kc)=(uz(ic,jc,kc)-uz(ic,jm,kc))*sucac
      enddo
    enddo
  ENDDO

  ! dwdz
  DO kc=y1start(3)-1,y1end(3)+1
    kp=kc+1
    do jc=0,nyp
      do ic=y1start(1)-1,y1end(1)+1
        dwdz(ic,jc,kc)=(uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
      enddo
    enddo
  ENDDO

  DO pid=1,nlocal
    pos    = GPrtcl_posR(pid)
    itype  = GPrtcl_pType(pid)
    FpVel  = GPrtcl_Vfluid(1,pid)- GPrtcl_linVel(1,pid)
    SumXDir= zero;  SumYDir=zero;  SumZDir=zero
    massFeff = FluidAccCoe*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid

    idxc_interp=indxyz(1,pid)
    idxp_interp=indxyz(2,pid)
    idyc_interp=indxyz(3,pid)
    idyp_interp=indxyz(4,pid)
    idzc_interp=indxyz(5,pid)
    idzp_interp=indxyz(6,pid)

    prx=(pos%x-xc(idxp_interp))*rdx+half
    RatioXp(0)=half*(prx-one)*(prx-two)
    RatioXp(1)=      prx*(two-prx)
    RatioXp(2)=half* prx*(prx-one)
    prx=(pos%x-xc(idxc_interp))*rdx
    RatioXc(0)=half*(prx-one)*(prx-two)
    RatioXc(1)=      prx*(two-prx)
    RatioXc(2)=half* prx*(prx-one)          

    RatioYp(0) = RatioYp_interp(pid)%x
    RatioYp(1) = RatioYp_interp(pid)%y
    RatioYp(2) = RatioYp_interp(pid)%z
    RatioYc(0) = RatioYc_interp(pid)%x
    RatioYc(1) = RatioYc_interp(pid)%y
    RatioYc(2) = RatioYc_interp(pid)%z

    prz=(pos%z-zc(idzp_interp))*rdz+half
    RatioZp(0)=half*(prz-one)*(prz-two)
    RatioZp(1)=      prz*(two-prz)
    RatioZp(2)=half* prz*(prz-one) 
    prz=(pos%z-zc(idzc_interp))*rdz
    RatioZc(0)=half*(prz-one)*(prz-two)
    RatioZc(1)=      prz*(two-prz)
    RatioZc(2)=half* prz*(prz-one)

    ! dwdx grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd  = j+idyc_interp
        pry = RatioYc(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumXDir= SumXDir + dwdx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dwdy gird
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumYDir= SumYDir + dwdy(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dwdz grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumZDir= SumZDir + dwdz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    GPrtcl_FpForce(pid)%z = GPrtcl_FpForce(pid)%z + massFeff*(SumXDir*FpVel%x +SumYDir*FpVel%y +SumZDir*FpVel%z)
  ENDDO

  ! Fixed particle part
  DO pid=1,mlocalFix
    pos    = GPFix_posR(pid)
    itype  = GPFix_pType(pid)
    FpVel  = GPFix_Vfluid(1,pid)
    SumXDir= zero;  SumYDir=zero;  SumZDir=zero
    massFeff = FluidAccCoe*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid

    idxc_interp=indxyzFixed(1,pid)
    idxp_interp=indxyzFixed(2,pid)
    idyc_interp=indxyzFixed(3,pid)
    idyp_interp=indxyzFixed(4,pid)
    idzc_interp=indxyzFixed(5,pid)
    idzp_interp=indxyzFixed(6,pid)

    prx=(pos%x-xc(idxp_interp))*rdx+half
    RatioXp(0)=half*(prx-one)*(prx-two)
    RatioXp(1)=      prx*(two-prx)
    RatioXp(2)=half* prx*(prx-one)
    prx=(pos%x-xc(idxc_interp))*rdx
    RatioXc(0)=half*(prx-one)*(prx-two)
    RatioXc(1)=      prx*(two-prx)
    RatioXc(2)=half* prx*(prx-one)          

    RatioYp(0) = RYpFixed_interp(pid)%x
    RatioYp(1) = RYpFixed_interp(pid)%y
    RatioYp(2) = RYpFixed_interp(pid)%z
    RatioYc(0) = RYcFixed_interp(pid)%x
    RatioYc(1) = RYcFixed_interp(pid)%y
    RatioYc(2) = RYcFixed_interp(pid)%z

    prz=(pos%z-zc(idzp_interp))*rdz+half
    RatioZp(0)=half*(prz-one)*(prz-two)
    RatioZp(1)=      prz*(two-prz)
    RatioZp(2)=half* prz*(prz-one) 
    prz=(pos%z-zc(idzc_interp))*rdz
    RatioZc(0)=half*(prz-one)*(prz-two)
    RatioZc(1)=      prz*(two-prz)
    RatioZc(2)=half* prz*(prz-one)

    ! dwdx grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd  = j+idyc_interp
        pry = RatioYc(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumXDir= SumXDir + dwdx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dwdy gird
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumYDir= SumYDir + dwdy(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! dwdz grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumZDir= SumZDir + dwdz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    GPFix_FpForce(pid)%z = GPFix_FpForce(pid)%z + massFeff*(SumXDir*FpVel%x +SumYDir*FpVel%y +SumZDir*FpVel%z)
  ENDDO
  end associate asso_FluidAcc3

end subroutine clc_FluidAccForce
