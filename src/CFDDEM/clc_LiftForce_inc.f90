  !==================================       Lift force        ==================================!
  ! [1] Crowe, C. T., Schwarzkopf, J. D., Sommerfeld, M., & Tsuji, Y. (2011).  Boca Raton, FL, U.S.: CRC Press.
  !       Multiphase flows with droplets and particles, p96-p99
  ! [2] Norouzi H R , Zarghami R , Sotudeh-Gharebagh R , et al. (2016)
  !       Coupled CFD‐DEM Modeling: Formulation, impleMentation and application to Multiphase Flows, p300-p303
  ! [3] Mei R . International Journal of Multiphase Flow, 1992, 18(1):145-147.
  !       An approximate expression for the shear lift force on a spherical particle at finite reynolds number
  ! [4] E. Loth. AIAA JOURNAL. Vol. 46, No. 4, April 2008
  !       Lift of a Solid Spherical Particle Subject to Vorticity and/or Spin
  ! [5] Spandan V., Rodolfo O., Roberto V., Detelef L. (2016) J. Fluid Mech.
  !       Drag reduction in numerical two-phase Taylor-Couette turbulence using an Eular-Largrange approach

#define SaffManType2
!#define CLCMagnusForce

subroutine clc_LiftForce(ux,uy,uz,RatioYp_interp,RatioYc_interp)
  implicit none
  real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
  type(real3),dimension(GPrtcl_list%mlocal),intent(in)::RatioYp_interp,RatioYc_interp

  ! locals
  type(real4)::pos
  type(real3)::FpVel,SaffmanLift,vorticity_interp,MagnusLift,PRotVel
  type(haloInfo):: hi_vorx_interp, hi_vory_interp, hi_vorz_interp
  integer::ic,jc,kc,im,jm,km,ip,jp,kp,id,jd,kd,pid,itype,i,j,k,nlocal
  integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
  real(RK)::rep,diam,prx,pry,prz,SumXDir,SumYDir,SumZDir,sucac,veldiff
  real(RK)::SaffmanCoe,c_saff,beta_saff,magVorticity,MagnusCoe,c_magnus_veldiff,MagRotvel
  real(RK),dimension(0:2)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
  real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm):: vorx,vory,vorz

  nlocal=GPrtcl_list%nlocal
  vorx= zero;  vory=zero;  vorz=zero
  !================================== vorticity in x-dir
  DO kc=y1start(3),y1end(3)+1
    km=kc-1
    do jc=1,nyp
      jm=jc-1;  sucac=rdyc(jc)
      do ic=y1start(1)-1,y1end(1)+1
        vorx(ic,jc,kc)=(uz(ic,jc,kc)-uz(ic,jm,kc))*sucac -(uy(ic,jc,kc)-uy(ic,jc,km))*rdz
      enddo
    enddo
  ENDDO 
  hi_vorx_interp%pencil = y_pencil
  hi_vorx_interp%xmh=0; hi_vorx_interp%xph=0
  hi_vorx_interp%ymh=0; hi_vorx_interp%yph=0
  hi_vorx_interp%zmh=1; hi_vorx_interp%zph=2
  call myupdate_halo(vorx, mb1, hi_vorx_interp)

  !================================== vorticity in y-dir
  DO kc=y1start(3),y1end(3)+1
    km=kc-1
    do jc=0,nyp
      do ic=y1start(1),y1end(1)+1
        im=ic-1
        vory(ic,jc,kc)=(ux(ic,jc,kc)-ux(ic,jc,km))*rdz -(uz(ic,jc,kc)-uz(im,jc,kc))*rdx
      enddo
    enddo
  ENDDO
  hi_vory_interp%pencil = y_pencil
  hi_vory_interp%xmh=1; hi_vory_interp%xph=2
  hi_vory_interp%ymh=0; hi_vory_interp%yph=0
  hi_vory_interp%zmh=1; hi_vory_interp%zph=2
  call myupdate_halo(vory, mb1, hi_vory_interp)

  !================================== vorticity in z-dir
  DO kc=y1start(3)-1,y1end(3)+1
    do jc=1,nyp
      jm=jc-1; sucac=rdyc(jc)
      do ic=y1start(1),y1end(1)+1
        im=ic-1
        vorz(ic,jc,kc)=(uy(ic,jc,kc)-uy(im,jc,kc))*rdx -(ux(ic,jc,kc)-ux(ic,jm,kc))*sucac 
      enddo
    enddo
  ENDDO
  hi_vorz_interp%pencil = y_pencil
  hi_vorz_interp%xmh=1; hi_vorz_interp%xph=2
  hi_vorz_interp%ymh=0; hi_vorz_interp%yph=0
  hi_vorz_interp%zmh=0; hi_vorz_interp%zph=0
  call myupdate_halo(vorz, mb1, hi_vorz_interp)

  SaffmanCoe= 1.615_RK*FluidDensity*sqrt(xnu)
  MagnusCoe = FluidDensity*PI/8.00_RK
  DO pid=1,nlocal
    pos   = GPrtcl_posR(pid)
    diam  = pos%w*two 
    FpVel = GPrtcl_Vfluid(1,pid)- GPrtcl_linVel(1,pid)
    veldiff=sqrt(FpVel%x*FpVel%x+ FpVel%y*FpVel%y+ FpVel%z*FpVel%z)
    rep= veldiff *diam/xnu
  
    SumXDir= zero;  SumYDir=zero;  SumZDir=zero
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

    ! vorx grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd  = j+idyp_interp
        pry = RatioYp(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumXDir= SumXDir + vorx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! vory gird
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumYDir= SumYDir + vory(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! vorz grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumZDir= SumZDir + vorz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    vorticity_interp =real3(SumXDir, SumYDir, SumZDir)

#ifdef SaffManType1
    magVorticity = sqrt(SumXDir*SumXDir+SumYDir*SumYDir+SumZDir*SumZDir)

    ! Saffman Lift, Eq.(4.146) in Crowe(2011). Or Eq(6.126) in Norouzi(2016)
    IF(magVorticity >1.0E-10_RK) THEN
      beta_saff = diam*sqrt(magVorticity/(two*xnu))
      if(rep<0.5_RK) then
        c_saff= one-0.1_RK *rep + 0.03314_RK*sqrt(rep)*beta_saff
      elseif(rep<=40_RK) then
        c_saff= exp(-0.1_RK*rep) +0.3314_RK*beta_saff/sqrt(rep)*(one-exp(-0.1_RK*rep))
      else
        c_saff= 0.0524_RK*beta_saff
      endif
      c_saff=min(one,c_saff) ! c_saff SMALLER THEN unity (according to [3])

      SaffmanLift = SaffmanCoe*c_saff*diam*diam/sqrt(magVorticity)* ( FpVel .cross. vorticity_interp)
    ELSE
      SaffmanLift = zero_r3 
    ENDIF
#else
    itype = GPrtcl_pType(pid)
    c_saff= SaffmanConst*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid    ! (according to [5])
    SaffmanLift = c_saff* ( FpVel .cross. vorticity_interp)
#endif
    
#ifdef CLCMagnusForce
    ! Magnus Lift, Eq(4.154) in Crowe(2011)
    PRotVel  = GPrtcl_RotVel(1,pid)
    MagRotvel= sqrt(PRotVel%x *PRotVel%x +PRotVel%y *PRotVel%y +PRotVel%z *PRotVel%z) 
    if(MagRotvel >1.0E-10_RK) then
      c_magnus_veldiff=0.45_RK*veldiff + (diam*MagRotvel -0.45_RK*veldiff) &
                       *exp(-0.075_RK*((half*diam*diam*MagRotvel/xnu)**0.4_RK)*(rep**0.3_RK))
      PRotVel = (PRotVel/MagRotvel)  ! normalized rotational vector
      MagnusLift = MagnusCoe*c_magnus_veldiff*diam*diam*(FpVel .cross. PRotVel )
    else
      MagnusLift = zero_r3
    endif
#else
    MagnusLift = zero_r3
#endif

    ! According to Eq.(20) in Loth(2008)
    GPrtcl_FpForce(pid)= GPrtcl_FpForce(pid)+ SaffmanLift + MagnusLift
  ENDDO

  ! fixed particle part
  DO pid=1,mlocalFix
    pos   = GPFix_posR(pid)
    diam  = pos%w*two 
    FpVel = GPFix_Vfluid(1,pid)
    veldiff=sqrt(FpVel%x*FpVel%x+ FpVel%y*FpVel%y+ FpVel%z*FpVel%z)
    rep= veldiff *diam/xnu
  
    SumXDir= zero;  SumYDir=zero;  SumZDir=zero
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

    ! vorx grid
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd  = j+idyp_interp
        pry = RatioYp(j)*prz
        do i=0,2
          id = i+idxc_interp
          prx= RatioXc(i)
          SumXDir= SumXDir + vorx(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! vory gird
    do k=0,2
      kd = k+idzp_interp
      prz= RatioZp(k)
      do j=0,2
        jd = j+idyc_interp
        pry= RatioYc(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumYDir= SumYDir + vory(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo

    ! vorz grid
    do k=0,2
      kd = k+idzc_interp
      prz= RatioZc(k)
      do j=0,2
        jd = j+idyp_interp
        pry= RatioYp(j)*prz
        do i=0,2
          id = i+idxp_interp
          prx= RatioXp(i)
          SumZDir= SumZDir + vorz(id,jd,kd)*prx*pry
        enddo
      enddo
    enddo
    vorticity_interp =real3(SumXDir, SumYDir, SumZDir)

#ifdef SaffManType1
    magVorticity = sqrt(SumXDir*SumXDir+SumYDir*SumYDir+SumZDir*SumZDir)

    ! Saffman Lift, Eq.(4.146) in Crowe(2011). Or Eq(6.126) in Norouzi(2016)
    IF(magVorticity >1.0E-10_RK) THEN
      beta_saff = diam*sqrt(magVorticity/(two*xnu))
      if(rep<0.5_RK) then
        c_saff= one-0.1_RK *rep + 0.03314_RK*sqrt(rep)*beta_saff
      elseif(rep<=40_RK) then
        c_saff= exp(-0.1_RK*rep) +0.3314_RK*beta_saff/sqrt(rep)*(one-exp(-0.1_RK*rep))
      else
        c_saff= 0.0524_RK*beta_saff
      endif
      c_saff=min(one,c_saff) ! c_saff SMALLER THEN unity (according to [3])

      SaffmanLift = SaffmanCoe*c_saff*diam*diam/sqrt(magVorticity)* ( FpVel .cross. vorticity_interp)
    ELSE
      SaffmanLift = zero_r3 
    ENDIF
#else 
    itype = GPFix_pType(pid)
    c_saff= SaffmanConst*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid   ! (according to [5])
    SaffmanLift = c_saff* ( FpVel .cross. vorticity_interp)
#endif
    
    ! Magnus Lift, Eq(4.154) in Crowe(2011)
    MagnusLift = zero_r3

    ! According to Eq.(20) in Loth(2008)
    GPFix_FpForce(pid)= GPFix_FpForce(pid)+ SaffmanLift + MagnusLift
  ENDDO

end subroutine clc_LiftForce
