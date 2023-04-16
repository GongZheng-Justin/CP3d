!==================================      viscous force      ==================================!
subroutine clc_VisForce(ux,uy,uz,RatioYp_interp,RatioYc_interp)
  implicit none
  real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
  type(real3),dimension(GPrtcl_list%mlocal),intent(in)::RatioYp_interp,RatioYc_interp

  ! locals
  type(real4)::pos
  type(real3)::VisForce
  type(HaloInfo):: hi_visx_interp, hi_visz_interp
  integer::ic,jc,kc,im,jm,km,ip,jp,kp,id,jd,kd,pid,itype,i,j,k,nlocal
  integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
  real(RK)::d11q1,d22q1,d33q1,d11q2,d22q2,d33q2,d11q3,d22q3,d33q3
  real(RK)::prx,pry,prz,SumXDir,SumYDir,SumZDir,qdx1,qdx3
  real(RK),dimension(0:2)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
  real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm):: visx,visy,visz

  nlocal=GPrtcl_list%nlocal
  visx=0.0_RK; visy=0.0_RK; visz=0.0_RK
  qdx1=0.25_RK*rdx;  qdx3=0.25_RK*rdz

  !================================== visx
  DO kc=y1start(3),y1end(3)
    km=kc-1;  kp=kc+1
    do jc=y1start(2),y1end(2)
      jm=jc-1;  jp=jc+1           
      do ic=y1start(1)+1,y1end(1)
        im=ic-1;  ip=ic+1
        d11q1= (ux(ip,jc,kc)-2.0_RK*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
        d22q1= ap2c(jc)*ux(ic,jp,kc)+ ac2c(jc)*ux(ic,jc,kc)+ am2c(jc)*ux(ic,jm,kc)                
        d33q1= ap3c(kc)*ux(ic,jc,kp)+ ac3c(kc)*ux(ic,jc,kc)+ am3c(kc)*ux(ic,jc,km)
        visx(ic,jc,kc)=d11q1+d22q1+d33q1
      enddo
    enddo
  ENDDO
  IF(myProcNghBC(y_pencil,4)<0) THEN  !xm_dir
    ic=1
    DO kc=y1start(3),y1end(3)
      km=kc-1;  kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;  jp=jc+1  

        d11q1= (ux(3,jc,kc)-2.0_RK*ux(2,jc,kc)+ux(1,jc,kc))*rdx2
        d22q1= ap2c(jc)*ux(ic,jp,kc) + ac2c(jc)*ux(ic,jc,kc)+ am2c(jc)*ux(ic,jm,kc)                
        d33q1= ap3c(kc)*ux(ic,jc,kp) + ac3c(kc)*ux(ic,jc,kc)+ am3c(kc)*ux(ic,jc,km)
        visx(ic,jc,kc)=d11q1+d22q1+d33q1
      enddo
    ENDDO
  ELSE
    ic=y1start(1); im=ic-1;  ip=ic+1
    DO kc=y1start(3),y1end(3)
      km=kc-1;  kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;  jp=jc+1  

        d11q1= (ux(ip,jc,kc)-2.0_RK*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
        d22q1= ap2c(jc)*ux(ic,jp,kc) + ac2c(jc)*ux(ic,jc,kc)+ am2c(jc)*ux(ic,jm,kc)                
        d33q1= ap3c(kc)*ux(ic,jc,kp) + ac3c(kc)*ux(ic,jc,kc)+ am3c(kc)*ux(ic,jc,km)
        visx(ic,jc,kc)=d11q1+d22q1+d33q1
      enddo
    ENDDO
  ENDIF

  IF(myProcNghBC(y_pencil,3)<0) THEN  !xp_dir
    ic=nxp
    DO kc=y1start(3),y1end(3)
      km=kc-1;  kp=kc+1
      do jc=y1start(2),y1end(2)
        jm=jc-1;  jp=jc+1   

        d11q1= (ux(nxp,jc,kc)-2.0_RK*ux(nxc,jc,kc)+ux(nxc-1,jc,kc))*rdx2
        d22q1= ap2c(jc)*ux(ic,jp,kc) + ac2c(jc)*ux(ic,jc,kc)+ am2c(jc)*ux(ic,jm,kc)                
        d33q1= ap3c(kc)*ux(ic,jc,kp) + ac3c(kc)*ux(ic,jc,kc)+ am3c(kc)*ux(ic,jc,km)
        visx(ic,jc,kc)=d11q1+d22q1+d33q1
      enddo
    ENDDO
  ENDIF

  ! Here I set 0.0_RK gradient Bcs for visx TERM for all the non-periodic directions !!!
  ! yp-dir and ym_dir
  do kc=y1start(3),y1end(3)
    do ic=y1start(1),y1end(1)
      visx(ic, nyp, kc) = visx(ic, nyc, kc)
      visx(ic, 0, kc)   = visx(ic, 1,   kc)
    enddo
  enddo

  ! zp-dir and zm_dir
  IF(myProcNghBC(y_pencil, 1)<0) THEN
    visx(y1start(1):y1end(1), 0:nyp, nzp) = visx(y1start(1):y1end(1), 0:nyp, nzc)
  ENDIF
  IF(myProcNghBC(y_pencil, 2)<0) THEN
    visx(y1start(1):y1end(1), 0:nyp, 0) = visx(y1start(1):y1end(1), 0:nyp, 1)
  ENDIF

  !!! update Halo here !!!
  hi_visx_interp%pencil = y_pencil
  hi_visx_interp%xmh=1;  hi_visx_interp%xph=2
  hi_visx_interp%ymh=0;  hi_visx_interp%yph=0
  hi_visx_interp%zmh=1;  hi_visx_interp%zph=1
  call update_halo(visx, mb1, hi_visx_interp) 

  !================================== visy
  DO kc=y1start(3),y1end(3)
    km=kc-1; kp=kc+1
    do jc=2,nyc
      jm=jc-1; jp=jc+1
      do ic=y1start(1),y1end(1)
        im=ic-1; ip=ic+1  
        d11q2= ap1c(ic)*uy(ip,jc,kc)+ac1c(ic)*uy(ic,jc,kc)+am1c(ic)*uy(im,jc,kc)
        d22q2= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
        d33q2= ap3c(kc)*uy(ic,jc,kp)+ac3c(kc)*uy(ic,jc,kc)+am3c(kc)*uy(ic,jc,km)
        visy(ic,jc,kc)=d11q2+d22q2+d33q2
      enddo
    enddo
  ENDDO
  jc=1
  DO kc=y1start(3),y1end(3)
    km=kc-1; kp=kc+1
    do ic=y1start(1),y1end(1)
      im=ic-1; ip=ic+1  
      d11q2= ap1c(ic)*uy(ip,jc,kc)+ac1c(ic)*uy(ic,jc,kc)+am1c(ic)*uy(im,jc,kc)
      d22q2= 2.0_RK*rdyp(2)*( (uy(ic,3,kc)-uy(ic,1,kc))/(yp(3)-yp(1))           &
                          -(uy(ic,2,kc)-uy(ic,1,kc))*rdyp(1) )
      d33q2= ap3c(kc)*uy(ic,jc,kp)+ac3c(kc)*uy(ic,jc,kc)+am3c(kc)*uy(ic,jc,km)
      visy(ic,jc,kc)=d11q2+d22q2+d33q2
    enddo
  ENDDO
  jc=nyp
  DO kc=y1start(3),y1end(3)
    km=kc-1; kp=kc+1
    do ic=y1start(1),y1end(1)
      im=ic-1; ip=ic+1  
      d11q2= ap1c(ic)*uy(ip,jc,kc)+ac1c(ic)*uy(ic,jc,kc)+am1c(ic)*uy(im,jc,kc)
      d22q2= 2.0_RK*rdyp(nyc-1)*( (uy(ic,nyp,kc)-uy(ic,nyc,  kc))*rdyp(nyc)     &
                              -(uy(ic,nyp,kc)-uy(ic,nyc-1,kc))/(yp(nyp)-yp(nyc-1)))
      d33q2= ap3c(kc)*uy(ic,jc,kp)+ac3c(kc)*uy(ic,jc,kc)+am3c(kc)*uy(ic,jc,km)
      visy(ic,jc,kc)=d11q2+d22q2+d33q2
    enddo
  ENDDO

  ! Here I set 0.0_RK gradient Bcs for visy TERM for all the non-periodic directions !!!
  ! xp-dir and xm_dir
  IF(myProcNghBC(y_pencil, 3)<0) THEN
    visy(nxp, 1:nyp,y1start(3):y1end(3))= visy(nxc, 1:nyp,y1start(3):y1end(3)) 
  ENDIF
  IF(myProcNghBC(y_pencil, 4)<0) THEN
    visy(0, 1:nyp,y1start(3):y1end(3))  = visy(1, 1:nyp,y1start(3):y1end(3)) 
  ENDIF 

  ! zp-dir and zm_dir
  IF(myProcNghBC(y_pencil, 1)<0) THEN
    visy(y1start(1):y1end(1), 1:nyp, nzp) = visy(y1start(1):y1end(1), 1:nyp, nzc)
  ENDIF
  IF(myProcNghBC(y_pencil, 2)<0) THEN
    visy(y1start(1):y1end(1), 1:nyp, 0) = visy(y1start(1):y1end(1), 1:nyp, 1)
  ENDIF

  !!! update Halo here !!!
  call update_halo(visy, mb1, hi1)

  !================================== visz
  DO kc=y1start(3)+1,y1end(3)
    km=kc-1;  kp=kc+1
    do jc=y1start(2),y1end(2)
      jm=jc-1;  jp=jc+1
      do ic=y1start(1),y1end(1)
        im=ic-1;  ip=ic+1
        d11q3= ap1c(ic)*uz(ip,jc,kc)+ac1c(ic)*uz(ic,jc,kc)+am1c(ic)*uz(im,jc,kc)
        d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
        d33q3= (uz(ic,jc,kp)-2.0_RK*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
        visz(ic,jc,kc)=d11q3+d22q3+d33q3
      enddo
    enddo
  ENDDO 

  IF(myProcNghBC(y_pencil,2)<0) THEN  !zm_dir
    kc=1
    do jc=y1start(2),y1end(2)
      jm=jc-1;  jp=jc+1
      do ic=y1start(1),y1end(1)
        im=ic-1;  ip=ic+1
        d11q3= ap1c(ic)*uz(ip,jc,kc)+ac1c(ic)*uz(ic,jc,kc)+am1c(ic)*uz(im,jc,kc)
        d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
        d33q3= (uz(ic,jc,3)-2.0_RK*uz(ic,jc,2)+uz(ic,jc,1))*rdz2
        visz(ic,jc,kc)=d11q3+d22q3+d33q3
      enddo
    enddo
  ELSE
    kc=y1start(3); km=kc-1;  kp=kc+1
    do jc=y1start(2),y1end(2)
      jm=jc-1;  jp=jc+1
      do ic=y1start(1),y1end(1)
        im=ic-1;  ip=ic+1
        d11q3= ap1c(ic)*uz(ip,jc,kc)+ac1c(ic)*uz(ic,jc,kc)+am1c(ic)*uz(im,jc,kc)
        d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
        d33q3= (uz(ic,jc,kp)-2.0_RK*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
        visz(ic,jc,kc)=d11q3+d22q3+d33q3
      enddo
    enddo
  ENDIF

  IF(myProcNghBC(y_pencil,1)<0) THEN  !zp_dir
    kc=nzp
    do jc=y1start(2),y1end(2)
      jm=jc-1;  jp=jc+1
      do ic=y1start(1),y1end(1)
        im=ic-1;  ip=ic+1
        d11q3= ap1c(ic)*uz(ip,jc,kc)+ac1c(ic)*uz(ic,jc,kc)+am1c(ic)*uz(im,jc,kc)
        d22q3= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)                
        d33q3= (uz(ic,jc,nzp)-2.0_RK*uz(ic,jc,nzc)+uz(ic,jc,nzc-1))*rdz2
        visz(ic,jc,kc)=d11q3+d22q3+d33q3
      enddo
    enddo
  ENDIF

  ! Here I set 0.0_RK gradient Bcs for visz TERM for all the non-periodic directions !!!
  ! yp-dir and ym_dir
  do kc=y1start(3),y1end(3)
    do ic=y1start(1),y1end(1)
      visz(ic, nyp, kc) = visz(ic, nyc, kc)
      visz(ic, 0, kc)   = visz(ic, 1,   kc)
    enddo
  enddo

  ! zp-dir and zm_dir
  IF(myProcNghBC(y_pencil, 1)<0) THEN
    visz(y1start(1):y1end(1), 0:nyp, nzp) = visz(y1start(1):y1end(1), 0:nyp, nzc)
  ENDIF
  IF(myProcNghBC(y_pencil, 2)<0) THEN
    visz(y1start(1):y1end(1), 0:nyp, 0) = visz(y1start(1):y1end(1), 0:nyp, 1)
  ENDIF

  !!! update Halo here !!!
  hi_visz_interp%pencil = y_pencil
  hi_visz_interp%xmh=1;  hi_visz_interp%xph=1
  hi_visz_interp%ymh=0;  hi_visz_interp%yph=0
  hi_visz_interp%zmh=1;  hi_visz_interp%zph=2
  call update_halo(visz, mb1, hi_visz_interp)

  !================================== interpolation to get viscous force at particle position
  DO pid=1,nlocal
    pos  = GPrtcl_posR(pid) 
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

    ! dux_dxmxm grid, the same to ux grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd  = j+idyc_interp
        pry = RatioYc(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumXDir= SumXDir + visx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! duy_dxmdxm grid, the same to uy grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumYDir= SumYDir + visy(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! duz_dxmdxm grid, the same to uz grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumZDir= SumZDir + visz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    VisForce= DEMProperty%Prtcl_PureProp(itype)%MassOfFluid*xnu*real3(SumXDir,SumYDir,SumZDir)
    GPrtcl_FpForce(pid)=GPrtcl_FpForce(pid)+ VisForce
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

    ! dux_dxmxm grid, the same to ux grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd  = j+idyc_interp
        pry = RatioYc(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumXDir= SumXDir + visx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! duy_dxmdxm grid, the same to uy grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumYDir= SumYDir + visy(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! duz_dxmdxm grid, the same to uz grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumZDir= SumZDir + visz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    VisForce= DEMProperty%Prtcl_PureProp(itype)%MassOfFluid* xnu*real3(SumXDir,SumYDir,SumZDir)
    GPFix_FpForce(pid)=GPFix_FpForce(pid)+ VisForce
  ENDDO

end subroutine clc_VisForce
