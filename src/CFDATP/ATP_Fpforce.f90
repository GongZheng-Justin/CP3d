module ATP_Fpforce
  use MPI
  use m_LogInfo
  use m_TypeDef
  use m_Decomp2d
  use m_Parameters
  use ATP_Property
  use ATP_Variables
  use ATP_parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1
  implicit none
  private

  ! Active information
  type(real3)::ReorientationVec=real3(0.0_RK,1.0_RK,0.0_RK)
  
  real(RK),allocatable,dimension(:):: xc   ! center coordinate in x-dir
  real(RK),allocatable,dimension(:):: zc   ! center coordinate in z-dir
  real(RK),allocatable,dimension(:)::iDistRatioYp,iDistRatioYc

  integer(kind=2),allocatable,dimension(:,:)::indxyz
  type(real3),allocatable,dimension(:)::RatioYp_interp,RatioYc_interp

  type(HaloInfo):: hi_ux_interp,   hi_uz_interp       ! halo info type for interpolation(velocity)
  integer:: xmp_interp, xep_interp, zmp_interp, zep_interp ! index constraints for interpolation in xp_dir,xm_dir,zp_dir,zm_dir
  integer:: xmc_interp, xec_interp, zmc_interp, zec_interp ! index constraints for interpolation in xp_dir,xm_dir,zp_dir,zm_dir
  
  public:: InitFpForce,PrepareInterpolation,clc_FluidInterpolation,FinalFpForce
contains
 
  !******************************************************************
  ! InitFpForce
  !******************************************************************
  subroutine InitFpForce()
    implicit none
    
    ! locals
    integer::j,nUnitFile,ierror
 
    ! check integer(kind=2) is enough or not.
    if(nrank==0 .and. (nxc>huge(0_2)-20 .or. nyc>huge(0_2)-20 .or. nzc>huge(0_2)-20)) then
      call MainLog%CheckForError(ErrT_Abort,"InitFpForce","kind=2 is not enough for indxyz")
    endif
    
    ! xc,zc, center coordinate interval in x-dir and z-dir
    allocate(xc(0:nxp))
    allocate(zc(0:nzp))
    do j=0,nxp
      xc(j)=dx*(real(j,RK)-0.5_RK)
    enddo
    do j=0,nzp
      zc(j)=dz*(real(j,RK)-0.5_RK)
    enddo
              
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
    
    ! calculate inverse distribution retio
    allocate(iDistRatioYp(0:nyp),iDistRatioYc(0:nyp))
    do j=0,nyp
      iDistRatioYp(j)=rdx*rdyc(j)*rdz/FluidDensity ! Note Here
      iDistRatioYc(j)=rdx*rdyp(j)*rdz/FluidDensity 
    enddo
  end subroutine InitFpForce
    
  !******************************************************************
  ! PrepareInterpolation
  !******************************************************************
  subroutine PrepareInterpolation()
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
      RatioYp%z=1.0_RK-RatioYp%x-RatioYp%y

      y1cord= yc(idyc_interp  )
      y2cord= yc(idyc_interp+1)
      y3cord= yc(idyc_interp+2)
      RatioYc%x=((pos%y-y2cord)*(pos%y-y3cord))/((y1cord-y2cord)*(y1cord-y3cord))
      RatioYc%y=((pos%y-y1cord)*(pos%y-y3cord))/((y2cord-y1cord)*(y2cord-y3cord))
      RatioYc%z=1.0_RK-RatioYc%x-RatioYc%y

      RatioYp_interp(pid)=RatioYp
      RatioYc_interp(pid)=RatioYc
    ENDDO
  end subroutine PrepareInterpolation
  
  !******************************************************************
  ! clc_FluidInterpolation
  !******************************************************************
  subroutine clc_FluidInterpolation(ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz

    ! locals
    type(real4)::pos
    type(haloInfo)::hi_vor
    type(real3)::OmegaAtP,SwimDir
    integer::id,jd,kd,pid,ic,jc,kc,im,jm,km,nlocal,iType
    real(RK)::prx,pry,prz,SumXDir,SumYDir,SumZDir,sucac,OneDtwoB,SwimVelocityMag
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm)::vorx,vory,vorz
    real(RK)::RatioXc(0:2),RatioYc(0:2),RatioZc(0:2),RatioXp(0:2),RatioYp(0:2),RatioZp(0:2)
      
    vorx= 0.0_RK;  vory=0.0_RK;  vorz=0.0_RK
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
    hi_vor%pencil = y_pencil
    hi_vor%xmh=0; hi_vor%xph=0
    hi_vor%ymh=0; hi_vor%yph=0
    hi_vor%zmh=1; hi_vor%zph=2
    call update_halo(vorx, mb1, hi_vor)

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
    hi_vor%pencil = y_pencil
    hi_vor%xmh=1; hi_vor%xph=2
    hi_vor%ymh=0; hi_vor%yph=0
    hi_vor%zmh=1; hi_vor%zph=2
    call update_halo(vory, mb1, hi_vor)

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
    hi_vor%pencil = y_pencil
    hi_vor%xmh=1; hi_vor%xph=2
    hi_vor%ymh=0; hi_vor%yph=0
    hi_vor%zmh=0; hi_vor%zph=0
    call update_halo(vorz, mb1, hi_vor)

    nlocal=GPrtcl_list%nlocal
    DO pid=1,nlocal
      do jc=GPrtcl_list%tsize-1,1,-1
        kc=jc+1
        GPrtcl_linVel(kc,pid) =GPrtcl_linVel(jc,pid)
        GPrtcl_SwimAcc(kc,pid)=GPrtcl_SwimAcc(jc,pid)
      enddo
    enddo
    
    ! perform interpolation
    DO pid=1,nlocal
      pos = GPrtcl_posR(pid)
      itype=GPrtcl_pType(pid)
      SwimDir=GPrtcl_SwimDir(pid)
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
      
      ! ux grid
      SumXDir=0.0_RK
      do kc=0,2
        kd = kc+idzc_interp
        prz= RatioZc(kc)
        do jc=0,2
          jd  = jc+idyc_interp
          pry = RatioYc(jc)*prz
          do ic=0,2
            id = ic+idxp_interp
            prx= RatioXp(ic)
            SumXDir= SumXDir + ux(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uy grid
      SumYDir=0.0_RK
      do kc=0,2
        kd = kc+idzc_interp
        prz= RatioZc(kc)
        do jc=0,2
          jd = jc+idyp_interp
          pry= RatioYp(jc)*prz
          do ic=0,2
            id = ic+idxc_interp
            prx= RatioXc(ic)
            SumYDir= SumYDir + uy(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo

      ! uz grid
      SumZDir=0.0_RK
      do kc=0,2
        kd = kc+idzp_interp
        prz= RatioZp(kc)
        do jc=0,2
          jd = jc+idyc_interp
          pry= RatioYc(jc)*prz
          do ic=0,2
            id = ic+idxc_interp
            prx= RatioXc(ic)
            SumZDir= SumZDir + uz(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo
      SwimVelocityMag=ATPProperty%Prtcl_PureProp(itype)%SwimVelocityMag
      GPrtcl_linVel(1,pid) =real3(SumXDir,SumYDir,SumZDir)+ SwimVelocityMag*SwimDir
      
      ! vor_x grid
      SumXDir=0.0_RK
      do kc=0,2
        kd = kc+idzp_interp
        prz= RatioZp(kc)
        do jc=0,2
          jd  = jc+idyp_interp
          pry = RatioYp(jc)*prz
          do ic=0,2
            id = ic+idxc_interp
            prx= RatioXc(ic)
            SumXDir= SumXDir + vorx(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo     
            
      ! vor_y grid
      SumYDir=0.0_RK
      do kc=0,2
        kd = kc+idzp_interp
        prz= RatioZp(kc)
        do jc=0,2
          jd  = jc+idyc_interp
          pry = RatioYc(jc)*prz
          do ic=0,2
            id = ic+idxp_interp
            prx= RatioXp(ic)
            SumYDir= SumYDir + vory(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo        
           
      ! vor_z grid
      SumZDir=0.0_RK
      do kc=0,2
        kd = kc+idzc_interp
        prz= RatioZc(kc)
        do jc=0,2
          jd  = jc+idyp_interp
          pry = RatioYp(jc)*prz
          do ic=0,2
            id = ic+idxp_interp
            prx= RatioXp(ic)
            SumZDir= SumZDir + vorz(id,jd,kd)*prx*pry
          enddo
        enddo
      enddo
      OneDtwoB=ATPProperty%Prtcl_PureProp(itype)%OneDtwoB
      OmegaAtP=real3(SumXDir,SumYDir,SumZDir)
      GPrtcl_SwimAcc(1,pid)= OneDtwoB*(ReorientationVec-(ReorientationVec .dot. SwimDir)*SwimDir) +0.5_RK*(OmegaAtP .cross. SwimDir)     
    ENDDO
  end subroutine clc_FluidInterpolation
  
  !******************************************************************
  ! FinalFpForce
  !******************************************************************
  subroutine FinalFpForce()
    implicit none
    
    if(GPrtcl_list%nlocal > 0) then
      deallocate(indxyz,RatioYp_interp,RatioYc_interp)
    endif
  end subroutine FinalFpForce
end module ATP_Fpforce
