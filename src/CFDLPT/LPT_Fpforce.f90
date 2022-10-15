module LPT_Fpforce
  use MPI
  use m_LogInfo
  use m_Decomp2d
  use LPT_TypeDef
  use m_Parameters
  use LPT_Property
  use LPT_Variables
  use LPT_parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1
#ifdef CFDLPT_TwoWay
  use m_Variables,only: FpForce_x,FpForce_y,FpForce_z
#endif
  implicit none
  private

#ifdef CFDFourthOrder
  real(RK),public,allocatable,dimension(:):: xc   ! center coordinate in x-dir
  real(RK),public,allocatable,dimension(:):: zc   ! center coordinate in z-dir
#endif
  real(RK),allocatable,dimension(:)::iDistRatioYp,iDistRatioYc

  integer(kind=2),allocatable,dimension(:,:)::indxyz
  type(real3),allocatable,dimension(:)::RatioYp_interp,RatioYc_interp

  type(HaloInfo):: hi_ux_interp,   hi_uz_interp       ! halo info type for interpolation(velocity)
  integer:: xmp_interp, xep_interp, zmp_interp, zep_interp ! index constraints for interpolation in xp_dir,xm_dir,zp_dir,zm_dir
  integer:: xmc_interp, xec_interp, zmc_interp, zec_interp ! index constraints for interpolation in xp_dir,xm_dir,zp_dir,zm_dir
  
  procedure(),pointer::PrepareInterpolation,clc_VelInterpolation,distribute_FpForce
  public:: InitFpForce,PrepareInterpolation,clc_VelInterpolation
  public:: clc_FpForce,distribute_FpForce,FinalFpForce
contains
 
  !******************************************************************
  ! InitFpForce
  !******************************************************************
  subroutine InitFpForce(chFile)
    implicit none
    character(*),intent(in)::chFile
    
    ! locals
    integer::j,nUnitFile,ierror,InterpAccuracy
    namelist/CFDLPT_interpolation/InterpAccuracy
 
    ! check integer(kind=2) is enough or not.
    if(nrank==0 .and. (nxc>huge(0_2)-20 .or. nyc>huge(0_2)-20 .or. nzc>huge(0_2)-20)) then
      call MainLog%CheckForError(ErrT_Abort,"InitFpForce","kind=2 is not enough for indxyz")
    endif
    
#ifdef CFDFourthOrder
    ! xc,zc, center coordinate interval in x-dir and z-dir
    allocate(xc(0:nxp))
    allocate(zc(0:nzp))
    do j=0,nxp
      xc(j)=dx*(real(j,RK)-half)
    enddo
    do j=0,nzp
      zc(j)=dz*(real(j,RK)-half)
    enddo    
#endif
  
    open(newunit=nUnitFile, file=chFile, status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitFpForce","Cannot open file: "//trim(chFile))
    read(nUnitFile, nml=CFDLPT_interpolation)
    close(nUnitFile,IOSTAT=ierror)
    
    if(InterpAccuracy==1) then      ! tri-linear interpolation
      PrepareInterpolation => PrepareInterpolation_1
      clc_VelInterpolation => clc_VelInterpolation_1 
#if defined(CFDLPT_TwoWay)
      distribute_FpForce   => distribute_FpForce_1
#endif      
      xmp_interp=y1start(1)
      xep_interp=y1end(1)     
      zmp_interp=y1start(3)
      zep_interp=y1end(3)
      xmc_interp=y1start(1)-1;   if( myProcNghBC(y_pencil,4)<0 ) xmc_interp=1
      xec_interp=y1end(1);       if( myProcNghBC(y_pencil,3)<0 ) xec_interp=nxc -1
      zmc_interp=y1start(3)-1;   if( myProcNghBC(y_pencil,2)<0 ) zmc_interp=1
      zec_interp=y1end(3);       if( myProcNghBC(y_pencil,1)<0 ) zec_interp=nzc -1
             
    elseif(InterpAccuracy==2) then  ! Quadratic interpolation
      PrepareInterpolation => PrepareInterpolation_2
      clc_VelInterpolation => clc_VelInterpolation_2
#if defined(CFDLPT_TwoWay)
      distribute_FpForce   => distribute_FpForce_2
#endif            
      ! ux
      hi_ux_interp%pencil = y_pencil
      hi_ux_interp%xmh=0;  hi_ux_interp%xph=2
      hi_ux_interp%ymh=0;  hi_ux_interp%yph=0
      hi_ux_interp%zmh=0;  hi_ux_interp%zph=0

      ! uz
      hi_uz_interp%pencil = y_pencil
      hi_uz_interp%xmh=0;  hi_uz_interp%xph=0
      hi_uz_interp%ymh=0;  hi_uz_interp%yph=0
      hi_uz_interp%zmh=0;  hi_uz_interp%zph=2

      ! index for interpolation in xp_dir,xm_dir,zp_dir,zm_dir
      xmp_interp=y1start(1)-1;   if( myProcNghBC(y_pencil,4)<0 ) xmp_interp=1
      xep_interp=y1end(1);       if( myProcNghBC(y_pencil,3)<0 ) xep_interp=nxc -1 ! Modified by Zheng Gong,2021-09-23
      zmp_interp=y1start(3)-1;   if( myProcNghBC(y_pencil,2)<0 ) zmp_interp=1
      zep_interp=y1end(3);       if( myProcNghBC(y_pencil,1)<0 ) zep_interp=nzc -1 ! Modified by Zheng Gong,2021-09-23
      xmc_interp=y1start(1)-1;   if( myProcNghBC(y_pencil,4)<0 ) xmc_interp=1
      xec_interp=y1end(1)  -1;   if( myProcNghBC(y_pencil,3)<0 ) xec_interp=nxc -2
      zmc_interp=y1start(3)-1;   if( myProcNghBC(y_pencil,2)<0 ) zmc_interp=1
      zec_interp=y1end(3)  -1;   if( myProcNghBC(y_pencil,1)<0 ) zec_interp=nzc -2
    else
      call MainLog%CheckForError(ErrT_Abort,"InitFpForce","wrong LPT_opt%InterpAccuracy")
    endif
    
    ! calculate inverse distribution retio
    allocate(iDistRatioYp(0:nyp),iDistRatioYc(0:nyp))
    do j=0,nyp
      iDistRatioYp(j)=rdx*rdyc(j)*rdz/FluidDensity ! Note Here
      iDistRatioYc(j)=rdx*rdyp(j)*rdz/FluidDensity 
    enddo
    
  end subroutine InitFpForce

  !******************************************************************
  ! PrepareInterpolation_1
  !******************************************************************
  subroutine PrepareInterpolation_1()
    implicit none
    
    ! locals
    type(real3)::pos,RatioYp,RatioYc
    integer::pid,ic,jc,kc,nlocal,js,je
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    
    nlocal= GPrtcl_list%nlocal
    if(nlocal > 0) then
      allocate(indxyz(6,nlocal),RatioYp_interp(nlocal),RatioYc_interp(nlocal))
    endif

    DO pid=1,nlocal
      pos = GPrtcl_posR(pid)

      ! if pos%y is within [0,yly), jc will be within [1,nyc]
      js=0
      je=nyp+1
      do
        jc=(js+je)/2
        if(je-js==1) exit
        if(pos%y< yp(jc)) then
          je =jc
        else
          js =jc
        endif
      enddo
      ic= floor(pos%x*rdx)+1; ic=min(ic,y1end(1)); ic=max(ic,y1start(1));
      kc= floor(pos%z*rdz)+1; kc=min(kc,y1end(3)); kc=max(kc,y1start(3));
      if(jc>nyc) call MainLog%CheckForError(ErrT_Abort,"InitDistribute","wrong jc")

      ! index for interpolation as follow, first-order largrange-interpolation:
      idxp_interp=ic
      if(pos%x>xc(ic)) then
        idxc_interp= min(ic,  xec_interp)
      else
        idxc_interp= max(ic-1,xmc_interp)
      endif

      idyp_interp=jc
      if(pos%y>yc(jc)) then
        idyc_interp=min(jc,nyc-1)
      else
        idyc_interp=max(jc-1,1)
      endif

      idzp_interp=kc
      if(pos%z>zc(kc)) then
        idzc_interp= min(kc,  zec_interp)     
      else
        idzc_interp= max(kc-1,zmc_interp)
      endif
      
      indxyz(1,pid)=int(idxc_interp,2)
      indxyz(2,pid)=int(idxp_interp,2)
      indxyz(3,pid)=int(idyc_interp,2)
      indxyz(4,pid)=int(idyp_interp,2)
      indxyz(5,pid)=int(idzc_interp,2)
      indxyz(6,pid)=int(idzp_interp,2)

      RatioYp%x=(yp(idyp_interp+1)-pos%y)*rdyp(idyp_interp) !(y2cord-pos%y)/(y2cord-y1cord)
      RatioYp%y=one-RatioYp%x
      RatioYc%x=(yc(idyc_interp+1)-pos%y)*rdyc(idyc_interp+1)
      RatioYc%y=one-RatioYc%x
      RatioYp_interp(pid)=RatioYp
      RatioYc_interp(pid)=RatioYc
    ENDDO
  end subroutine PrepareInterpolation_1
    
  !******************************************************************
  ! PrepareInterpolation_2
  !******************************************************************
  subroutine PrepareInterpolation_2()
    implicit none
    
    ! locals
    real(RK)::y1cord,y2cord,y3cord
    type(real3)::pos,RatioYp,RatioYc
    integer::pid,ic,jc,kc,nlocal,js,je
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    
    nlocal= GPrtcl_list%nlocal
    if(nlocal > 0) then
      allocate(indxyz(6,nlocal),RatioYp_interp(nlocal),RatioYc_interp(nlocal))
    endif

    DO pid=1,nlocal
      pos = GPrtcl_posR(pid)

      ! if pos%y is within [0,yly), jc will be within [1,nyc]
      js=0
      je=nyp+1
      do
        jc=(js+je)/2
        if(je-js==1) exit
        if(pos%y< yp(jc)) then
          je =jc
        else
          js =jc
        endif
      enddo
      ic= floor(pos%x*rdx)+1; ic=min(ic,y1end(1)); ic=max(ic,y1start(1));
      kc= floor(pos%z*rdz)+1; kc=min(kc,y1end(3)); kc=max(kc,y1start(3));
      if(jc>nyc) call MainLog%CheckForError(ErrT_Abort,"InitDistribute","wrong jc")

      ! index for interpolation as follow, second-order largrange-interpolation:
      idxc_interp=max(xmc_interp,ic-1)          ! Modified by Zheng Gong,2021-09-23
      idxc_interp=min(xec_interp,idxc_interp)   ! Modified by Zheng Gong,2021-09-23
      if(pos%x>xc(ic)) then
        idxp_interp= min(ic,  xep_interp)       ! Modified by Zheng Gong,2021-09-23, "xep_interp-2" => "xep_interp"
      else
        idxp_interp= max(ic-1,xmp_interp)
      endif

      idyc_interp=max(1,    jc-1)               ! Modified by Zheng Gong,2021-09-23
      idyc_interp=min(nyc-2,idyc_interp)        ! Modified by Zheng Gong,2021-09-23
      if(pos%y>yc(jc)) then
        idyp_interp=min(jc,nyc-1)
      else
        idyp_interp=max(jc-1,1)
      endif

      idzc_interp=max(zmc_interp,kc-1)          ! Modified by Zheng Gong,2021-09-23
      idzc_interp=min(zec_interp,idzc_interp)   ! Modified by Zheng Gong,2021-09-23
      if(pos%z>zc(kc)) then
        idzp_interp= min(kc,  zep_interp)       ! Modified by Zheng Gong,2021-09-23, "zep_interp-2" => "zep_interp"
      else
        idzp_interp= max(kc-1,zmp_interp)
      endif

      indxyz(1,pid)=int(idxc_interp,2)
      indxyz(2,pid)=int(idxp_interp,2)
      indxyz(3,pid)=int(idyc_interp,2)
      indxyz(4,pid)=int(idyp_interp,2)
      indxyz(5,pid)=int(idzc_interp,2)
      indxyz(6,pid)=int(idzp_interp,2)

      y1cord= yp(idyp_interp  )
      y2cord= yp(idyp_interp+1)
      y3cord= yp(idyp_interp+2)
      RatioYp%x=((pos%y-y2cord)*(pos%y-y3cord))/((y1cord-y2cord)*(y1cord-y3cord))
      RatioYp%y=((pos%y-y1cord)*(pos%y-y3cord))/((y2cord-y1cord)*(y2cord-y3cord))
      RatioYp%z=one-RatioYp%x-RatioYp%y

      y1cord= yc(idyc_interp  )
      y2cord= yc(idyc_interp+1)
      y3cord= yc(idyc_interp+2)
      RatioYc%x=((pos%y-y2cord)*(pos%y-y3cord))/((y1cord-y2cord)*(y1cord-y3cord))
      RatioYc%y=((pos%y-y1cord)*(pos%y-y3cord))/((y2cord-y1cord)*(y2cord-y3cord))
      RatioYc%z=one-RatioYc%x-RatioYc%y

      RatioYp_interp(pid)=RatioYp
      RatioYc_interp(pid)=RatioYc
    ENDDO
  end subroutine PrepareInterpolation_2
  
  !******************************************************************
  ! clc_VelInterpolation_1
  !******************************************************************
  subroutine clc_VelInterpolation_1(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz

    ! locals
    type(real4)::pos
    integer::id,jd,kd,pid,i,j,k,nlocal
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    real(RK)::prx,pry,prz,SumXDir,SumYDir,SumZDir,RatioXc(0:1),RatioYc(0:1),RatioZc(0:1),RatioXp(0:1),RatioYp(0:1),RatioZp(0:1)

    nlocal=GPrtcl_list%nlocal
    DO pid=1,nlocal
      pos = GPrtcl_posR(pid)
      SumXDir=zero;  SumYDir=zero;  SumZDir=zero

      idxc_interp=indxyz(1,pid)
      idxp_interp=indxyz(2,pid)
      idyc_interp=indxyz(3,pid)
      idyp_interp=indxyz(4,pid)
      idzc_interp=indxyz(5,pid)
      idzp_interp=indxyz(6,pid)

      RatioXp(1)=(pos%x-xc(idxp_interp))*rdx+half
      RatioXp(0)=one-RatioXp(1)
      RatioXc(1)=(pos%x-xc(idxc_interp))*rdx
      RatioXc(0)=one-RatioXc(1)
 
      RatioYp(0)= RatioYp_interp(pid)%x
      RatioYp(1)= RatioYp_interp(pid)%y
      RatioYc(0)= RatioYc_interp(pid)%x
      RatioYc(1)= RatioYc_interp(pid)%y

      RatioZp(1)=(pos%z-zc(idzp_interp))*rdz+half
      RatioZp(0)=one-RatioZp(1)
      RatioZc(1)=(pos%z-zc(idzc_interp))*rdz
      RatioZc(0)=one-RatioZc(1)
      
      ! ux grid
      do k=0,1
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,1
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,1
            id = i+idxp_interp
            prx= RatioXp(i)
            SumXDir= SumXDir + ux(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uy gird
      do k=0,1
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,1
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,1
            id = i+idxc_interp
            prx= RatioXc(i)
            SumYDir= SumYDir + uy(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uz grid
      do k=0,1
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,1
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,1
            id = i+idxc_interp
            prx= RatioXc(i)
            SumZDir= SumZDir + uz(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo
      GPrtcl_Vfluid(pid)  = real3(SumXDir,SumYDir,SumZDir)
    ENDDO
  end subroutine clc_VelInterpolation_1
  
  !******************************************************************
  ! clc_VelInterpolation_2
  !******************************************************************
  subroutine clc_VelInterpolation_2(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz

    ! locals
    type(real4)::pos
    integer::id,jd,kd,pid,i,j,k,nlocal
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    real(RK)::prx,pry,prz,SumXDir,SumYDir,SumZDir,RatioXc(0:2),RatioYc(0:2),RatioZc(0:2),RatioXp(0:2),RatioYp(0:2),RatioZp(0:2)

#ifdef CFDSecondOrder
    call myupdate_halo(ux, mb1, hi_ux_interp)
    call myupdate_halo(uz, mb1, hi_uz_interp)
#endif
    nlocal=GPrtcl_list%nlocal
    DO pid=1,nlocal
      pos = GPrtcl_posR(pid)
      SumXDir=zero;  SumYDir=zero;  SumZDir=zero

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
      
      ! ux grid
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            SumXDir= SumXDir + ux(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uy gird
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumYDir= SumYDir + uy(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uz grid
      do k=0,2
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            SumZDir= SumZDir + uz(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo
      GPrtcl_Vfluid(pid)  = real3(SumXDir,SumYDir,SumZDir)
    ENDDO
  end subroutine clc_VelInterpolation_2

  !******************************************************************
  ! FinalFpForce
  !******************************************************************
  subroutine FinalFpForce()
    implicit none
    
    if(GPrtcl_list%nlocal > 0) then
      deallocate(indxyz,RatioYp_interp,RatioYc_interp)
    endif
  end subroutine FinalFpForce

  !******************************************************************
  ! clc_FpForce
  !******************************************************************
  subroutine clc_FpForce()!clc_FpForce(ux,uy,uz,pressure)
    implicit none
    !real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,pressure

    ! locals
    !type(real4)::pos
    type(real3)::FpVelDiff
    integer::pid,nlocal,itype
    real(RK)::cd,rep,diam,udiff,vdiff,wdiff,veldiff,taup,Mass

    nlocal=GPrtcl_list%nlocal
    DO pid=1,nlocal
      itype= GPrtcl_pType(pid)
      diam = two*GPrtcl_PosR(pid)%w
      taup = LPTProperty%Prtcl_PureProp(itype)%RelaxionTime
      Mass = LPTProperty%Prtcl_PureProp(itype)%Mass

      FpVelDiff = GPrtcl_Vfluid(pid)-GPrtcl_LinVel(1,pid) 
      udiff= FpVelDiff%x; vdiff= FpVelDiff%y; wdiff= FpVelDiff%z
      veldiff=sqrt(udiff*udiff+ vdiff*vdiff+ wdiff*wdiff)

      rep= veldiff *diam/xnu
      cd = one+0.15_RK*(rep**0.687_RK)
      GPrtcl_FpForce(pid)=(Mass*cd/taup)*FpVelDiff
    ENDDO
  end subroutine clc_FpForce

#ifdef CFDLPT_TwoWay
  !******************************************************************
  ! distribute_FpForce_1
  !****************************************************************** 
  subroutine distribute_FpForce_1()
    implicit none

    ! locals
    real(RK)::prx,pry,prz
    type(real3)::Pos,FpForce
    integer::i,j,k,id,jd,kd,pid
    real(RK),dimension(0:1)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    FpForce_x=zero; FpForce_y=zero; FpForce_z=zero
    DO pid=1,GPrtcl_list%nlocal
      Pos     = GPrtcl_PosR(pid)
      FpForce = zero_r3-GPrtcl_FpForce(pid)
      idxc_interp=indxyz(1,pid)
      idxp_interp=indxyz(2,pid)
      idyc_interp=indxyz(3,pid)
      idyp_interp=indxyz(4,pid)
      idzc_interp=indxyz(5,pid)
      idzp_interp=indxyz(6,pid)

      RatioXp(1)=(pos%x-xc(idxp_interp))*rdx+half
      RatioXp(0)=one-RatioXp(1)
      RatioXc(1)=(pos%x-xc(idxc_interp))*rdx
      RatioXc(0)=one-RatioXc(1)
 
      RatioYp(0)= RatioYp_interp(pid)%x
      RatioYp(1)= RatioYp_interp(pid)%y
      RatioYc(0)= RatioYc_interp(pid)%x
      RatioYc(1)= RatioYc_interp(pid)%y

      RatioZp(1)=(pos%z-zc(idzp_interp))*rdz+half
      RatioZp(0)=one-RatioZp(1)
      RatioZc(1)=(pos%z-zc(idzc_interp))*rdz
      RatioZc(0)=one-RatioZc(1)

      ! ux grid
      do k=0,1
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,1
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz*iDistRatioYc(jd)
          do i=0,1
            id = i+idxp_interp
            prx= RatioXp(i)
            FpForce_x(id,jd,kd)= FpForce_x(id,jd,kd)+ prx*pry* FpForce%x
          enddo
        enddo
      enddo

      ! uy gird
      do k=0,1
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,1
          jd = j+idyp_interp
          pry= RatioYp(j)*prz*iDistRatioYp(jd)
          do i=0,1
            id = i+idxc_interp
            prx= RatioXc(i)
            FpForce_y(id,jd,kd)= FpForce_y(id,jd,kd)+ prx*pry* FpForce%y
          enddo
        enddo
      enddo

      ! uz grid
      do k=0,1
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,1
          jd = j+idyc_interp
          pry= RatioYc(j)*prz*iDistRatioYc(jd)
          do i=0,1
            id = i+idxc_interp
            prx= RatioXc(i)
            FpForce_z(id,jd,kd)= FpForce_z(id,jd,kd)+ prx*pry* FpForce%z
          enddo
        enddo
      enddo
    ENDDO

    call Gather_Halo_dist_1()
  end subroutine distribute_FpForce_1
  
  !******************************************************************
  ! distribute_FpForce_2
  !****************************************************************** 
  subroutine distribute_FpForce_2()
    implicit none

    ! locals
    real(RK)::prx,pry,prz
    type(real3)::Pos,FpForce
    integer::i,j,k,id,jd,kd,pid
    real(RK),dimension(0:2)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    FpForce_x=zero; FpForce_y=zero; FpForce_z=zero
    DO pid=1,GPrtcl_list%nlocal
      Pos     = GPrtcl_PosR(pid)
      FpForce = zero_r3-GPrtcl_FpForce(pid)
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

      ! ux grid
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd  = j+idyc_interp
          pry = RatioYc(j)*prz*iDistRatioYc(jd)
          do i=0,2
            id = i+idxp_interp
            prx= RatioXp(i)
            FpForce_x(id,jd,kd)= FpForce_x(id,jd,kd)+ prx*pry* FpForce%x
          enddo
        enddo
      enddo

      ! uy gird
      do k=0,2
        kd = k+idzc_interp
        prz= RatioZc(k)
        do j=0,2
          jd = j+idyp_interp
          pry= RatioYp(j)*prz*iDistRatioYp(jd)
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            FpForce_y(id,jd,kd)= FpForce_y(id,jd,kd)+ prx*pry* FpForce%y
          enddo
        enddo
      enddo

      ! uz grid
      do k=0,2
        kd = k+idzp_interp
        prz= RatioZp(k)
        do j=0,2
          jd = j+idyc_interp
          pry= RatioYc(j)*prz*iDistRatioYc(jd)
          do i=0,2
            id = i+idxc_interp
            prx= RatioXc(i)
            FpForce_z(id,jd,kd)= FpForce_z(id,jd,kd)+ prx*pry* FpForce%z
          enddo
        enddo
      enddo
    ENDDO
    
    call Gather_Halo_dist_2()
  end subroutine distribute_FpForce_2
#include "LPT_GatherHalo_inc.f90"
#endif
end module LPT_Fpforce
