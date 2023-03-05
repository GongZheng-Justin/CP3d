module cd_FpForce
  use MPI
  use m_LogInfo
  use m_TypeDef
  use m_Decomp2d
  use m_Parameters  
  use Prtcl_Property
  use Prtcl_Variables
  use m_MeshAndMetries
  use Prtcl_parameters
  use Prtcl_Decomp_2d,only:DEM_decomp
  use m_Variables,only:mb1,hi1,mb_dist,hi_dist
  use m_Variables,only:FpForce_x,FpForce_y,FpForce_z
  implicit none
  private

  integer(kind=2),allocatable,dimension(:)::idxc,idyc,idzc
  integer(kind=2),allocatable,dimension(:,:)::indxyz
  real(RK),allocatable,dimension(:,:,:):: xphalo,xmhalo,xphalos,xmhalos,zphalo,zmhalo
  type(HaloInfo):: hi_ux_interp,   hi_uz_interp   ! halo info type for interpolation(velocity)

  integer:: xm_dist, xe_dist, zm_dist, ze_dist  ! index for distribution in xp_dir,xm_dir,zp_dir,zm_dir
  integer:: xmMin,   xeMax,   zmMin,   zeMax
  integer:: distwx,  distwz                     ! distribution width in x_dir and z_dir

  integer:: xm_interp, xe_interp, zm_interp, ze_interp  ! index for interpolation in xp_dir,xm_dir,zp_dir,zm_dir

  ! Basset force
  logical:: IsHaveVfluidOld=.false.
  real(RK),allocatable,dimension(:)  :: aTailBasset,tTailBasset,exp_tTailBasset,WindowsCoeBasset
  real(RK),allocatable,dimension(:,:):: tailDiCoeBasset

  ! fixed particle variables
  integer:: mlocalFix
  real(RK),allocatable,dimension(:)   :: GPFix_VolFluid
  integer(kind=2),allocatable,dimension(:)  :: idxcFixed,idycFixed,idzcFixed,GPFix_nFluidCell
  integer(kind=2),allocatable,dimension(:,:):: indxyzFixed
  type(real3),allocatable,dimension(:):: RYpFixed_interp,RYcFixed_interp,GPFix_FpForce

  public:: InitDistribute, PrepareDistribute, FinalDistribute
  public:: clc_FpForce,    distribute_FpForce
contains
 
#include "clc_VisForce_inc.f90"
#include "clc_LiftForce_inc.f90"
#include "clc_PrGradForce_inc.f90"
#include "clc_BassetForce_inc.f90"
#include "clc_FluidAccForce_inc.f90"

  !******************************************************************
  ! InitDistribute
  !******************************************************************
  subroutine InitDistribute()
    implicit none

    ! locals
    type(real4)::pos
    type(real3)::RatioYp,RatioYc
    real(RK)::maxR,y1cord,y2cord,y3cord
    integer::iErr01,iErr02,iErr03,iErr04,iErr05,iErr06,iErr07,iErr08,iErr09,iErrSum,js,je
    integer::i,j,k,ic,jc,kc,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    integer::pid,itype,nFluidCell,idy2,idxs,idys,idzs,idxe,idye,idze
    real(RK)::dcell2,SearchR,ry,dist,volp,rad
    real(RK)::volt,dist2,distx,disty,distz,disty2,distz2

    ! check integer(kind=2) is enough or not.
    if(nrank==0 .and. (nxc>huge(0_2)-20 .or. nyc>huge(0_2)-20 .or. nzc>huge(0_2)-20)) then
      call MainLog%CheckForError(ErrT_Abort,"InitDistribute","kind=2 is not enough for indxyz and indxyzFixed")
    endif
    
    maxR = maxval( DEMProperty%Prtcl_PureProp%Radius )
    distwx = ceiling(1.02_RK*RatioSR*maxR*rdx -0.5_RK); distwx = max(distwx,1)
    distwz = ceiling(1.02_RK*RatioSR*maxR*rdz -0.5_RK); distwz = max(distwz,1)

    xe_dist= y1end(1)  +distwx;     xeMax=xe_dist
    if( myProcNghBC(y_pencil,3)<0 ) xeMax=y1end(1)
    xm_dist= y1start(1)-distwx;     xmMin=xm_dist
    if( myProcNghBC(y_pencil,4)<0 ) xmMin=y1start(1)

    ze_dist= y1end(3)  +distwz;     zeMax=ze_dist
    if( myProcNghBC(y_pencil,1)<0 ) zeMax=y1end(3)
    zm_dist= y1start(3)-distwz;     zmMin=zm_dist
    if( myProcNghBC(y_pencil,2)<0 ) zmMin=y1start(3)

    mb_dist%pencil = y_pencil  
    mb_dist%xme=distwx;   mb_dist%xpe=distwx
    mb_dist%yme=1;        mb_dist%ype=1
    mb_dist%zme=distwz;   mb_dist%zpe=distwz
    hi_dist%pencil = y_pencil
    hi_dist%xmh=1;        hi_dist%xph=1
    hi_dist%ymh=0;        hi_dist%yph=0
    hi_dist%zmh=1;        hi_dist%zph=1
    call myallocate(FpForce_x, mb_dist, opt_global=.true.);  FpForce_x=0.0_RK
    call myallocate(FpForce_y, mb_dist, opt_global=.true.);  FpForce_y=0.0_RK
    call myallocate(FpForce_z, mb_dist, opt_global=.true.);  FpForce_z=0.0_RK

    allocate(xphalo(y1end(1)-distwx+1:y1end(1),     y1start(2)-1:y1end(2)+1, zm_dist:ze_dist))
    allocate(xphalos(y1end(1)-distwx+1:y1end(1),    y1start(2)-1:y1end(2)+1, zm_dist:ze_dist))
    allocate(xmhalo(y1start(1):y1start(1)+distwx-1, y1start(2)-1:y1end(2)+1, zm_dist:ze_dist))
    allocate(xmhalos(y1start(1):y1start(1)+distwx-1,y1start(2)-1:y1end(2)+1, zm_dist:ze_dist))
    allocate(zphalo(xm_dist:xe_dist,y1start(2)-1:y1end(2)+1,    y1end(3)-distwz+1:y1end(3)  ))
    allocate(zmhalo(xm_dist:xe_dist,y1start(2)-1:y1end(2)+1,  y1start(3):y1start(3)+distwz-1))

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
    xm_interp=y1start(1)-1;   if( myProcNghBC(y_pencil,4)<0 ) xm_interp=1
    xe_interp=y1end(1)  +2;   if( myProcNghBC(y_pencil,3)<0 ) xe_interp=nxp
    zm_interp=y1start(3)-1;   if( myProcNghBC(y_pencil,2)<0 ) zm_interp=1
    ze_interp=y1end(3)  +2;   if( myProcNghBC(y_pencil,1)<0 ) ze_interp=nzp

    ! Basset history force properties
    mlocalFix= GPrtcl_list%mlocalFix
    if(is_clc_Basset) then
      allocate(WindowsCoeBasset(mWinBasset+1 )); WindowsCoeBasset=clc_muCoeBasset(mWinBasset)
      allocate(aTailBasset(mTailBasset),tTailBasset(mTailBasset),exp_tTailBasset(mTailBasset),tailDiCoeBasset(mTailBasset,4))
      call clc_TailCoeBasset()
      if(is_clc_Basset_fixed .and. (.not.DEM_Opt%RestartFlag) .and. mlocalFix>0) then
        allocate (GPFix_BassetData(GPrtcl_BassetSeq%nDataLen, mlocalFix), Stat= iErr01)
        if(iErr01 /= 0 ) call MainLog%CheckForError(ErrT_Abort,"InitDistribute","allocate wrong1")
        GPFix_BassetData=zero_r3
      endif
    endif

    ! fixed particle variables
    if(mlocalFix>0) then
      allocate(idxcFixed(mlocalFix),        Stat= iErr01)
      allocate(idycFixed(mlocalFix),        Stat= iErr02)
      allocate(idzcFixed(mlocalFix),        Stat= iErr03)
      allocate(indxyzFixed(6,mlocalFix),    Stat= iErr04)
      allocate(RYpFixed_interp(mlocalFix),  Stat= iErr05)
      allocate(RYcFixed_interp(mlocalFix),  Stat= iErr06)
      allocate(GPFix_FpForce(mlocalFix),    Stat= iErr07)
      allocate(GPFix_VolFluid(mlocalFix),   Stat= iErr08)
      allocate(GPFix_nFluidCell(mlocalFix), Stat= iErr09)
      iErrSum= abs(iErr01) +abs(iErr02) +abs(iErr03) +abs(iErr04) +abs(iErr05) &
              +abs(iErr06) +abs(iErr07) +abs(iErr08) +abs(iErr09)
      if(iErrSum/=0) call MainLog%CheckForError(ErrT_Abort,"InitDistribute","allocate wrong2")
      idxcFixed=0;     idycFixed=0; idzcFixed=0; indxyzFixed=0; 
      RYpFixed_interp=zero_r3; RYcFixed_interp=zero_r3;
      GPFix_FpForce=zero_r3;   GPFix_VolFluid=0.0_RK;     GPFix_nFluidCell=0
    endif

    DO pid=1,mlocalFix
      pos = GPFix_PosR(pid)

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

      idxcFixed(pid)=ic
      idycFixed(pid)=jc
      idzcFixed(pid)=kc

      ! index for interpolation as follow, second-order largrange-interpolation:
      idxc_interp=ic-1
      if(pos%x>xc(ic)) then
        idxp_interp= min(ic,  xe_interp-2)
      else
        idxp_interp= max(ic-1,xm_interp)
      endif

      idyc_interp=jc-1
      if(pos%y>yc(jc)) then
        idyp_interp=min(jc,nyc-1)
      else
        idyp_interp=max(jc-1,1)
      endif

      idzc_interp=kc-1
      if(pos%z>zc(kc)) then
        idzp_interp= min(kc,  ze_interp-2)
      else
        idzp_interp= max(kc-1,zm_interp)
      endif

      indxyzFixed(1,pid)=idxc_interp
      indxyzFixed(2,pid)=idxp_interp
      indxyzFixed(3,pid)=idyc_interp
      indxyzFixed(4,pid)=idyp_interp
      indxyzFixed(5,pid)=idzc_interp
      indxyzFixed(6,pid)=idzp_interp

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

      RYpFixed_interp(pid)=RatioYp
      RYcFixed_interp(pid)=RatioYc

      ! Following is to calculate the GPFix_VolFluid, and GPFix_nFluidCell
      itype= GPFix_pType(pid);   rad = pos%w
      volp= DEMProperty%Prtcl_PureProp(itype)%Volume

      SearchR  = RatioSR * rad
      idxs= max(ceiling((pos%x-SearchR)*rdx), xmMin)
      idxe= min(ceiling((pos%x+SearchR)*rdx), xeMax)
      idzs= max(ceiling((pos%z-SearchR)*rdz), zmMin)
      idze= min(ceiling((pos%z+SearchR)*rdz), zeMax)

      ry = pos%y-SearchR
      idy2= 1
      do j=jc,1,-1
        if(ry >= yp(j)) then
         idy2 = j; exit
        endif      
      enddo
      idys = idy2

      ry  = pos%y+SearchR
      idy2= nyc
      do j=jc+1,nyp
        if(ry < yp(j)) then
         idy2 = j-1; exit
        endif      
      enddo
      idye = idy2
 
      volt= 0.0_RK; nFluidCell=0
      do k=idzs,idze
        distz = zc(k) - pos%z
        distz2= distz*distz
        do j=idys,idye
          disty = yc(j) - pos%y
          disty2= disty*disty
          do i=idxs,idxe
            distx = xc(i) - pos%x
            dist  = RatioSR-sqrt(distx*distx +disty2 +distz2)/rad
            if(dist>0.0_RK) then
              nFluidCell=nFluidCell+1
              volt=volt+VolCell(j)*exp(-0.125_RK*dist*dist)
            endif
          enddo
        enddo
      enddo

      IF(nFluidCell<27) THEN
        volt=0.0_RK
        dcell2=DeltaCell(jc)*DeltaCell(jc)
        do k= kc-1,kc+1
          distz = zc(k) - pos%z
          distz2= distz*distz
          do j=jc-1,jc+1
            disty = yc(j) - pos%y
            disty2= disty*disty
            do i=ic-1,ic+1
              distx = xc(i) - pos%x
              dist2 = (distx*distx+disty2+distz2)/dcell2
              volt  = volt + VolCell(j)/(16.0_RK**dist2)
            enddo
          enddo
        enddo
      ENDIF
      GPFix_VolFluid(pid)  =volt
      GPFix_nFluidCell(pid)=nFluidCell
    ENDDO
    if(DEM_Opt%RestartFlag) IsHaveVfluidOld=.true.
  end subroutine InitDistribute

  !******************************************************************
  ! PrepareDistribute
  !******************************************************************
  subroutine PrepareDistribute(RatioYp_interp,RatioYc_interp)
    implicit none
    type(real3),dimension(GPrtcl_list%mlocal),intent(out)::RatioYp_interp,RatioYc_interp
    
    ! locals
    type(real3)::pos,RatioYp,RatioYc
    integer::pid,j,ic,jc,kc,nlocal,js,je
    integer::idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    real(RK)::y1cord,y2cord,y3cord

    nlocal= GPrtcl_list%nlocal
    if(nlocal > 0) then
      allocate(idxc(nlocal), idyc(nlocal), idzc(nlocal))
      allocate(indxyz(6,nlocal))
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
      if(jc>nyc) call MainLog%CheckForError(ErrT_Abort,"PrepareDistribute","wrong jc")

      idxc(pid)=ic
      idyc(pid)=jc
      idzc(pid)=kc

      ! index for interpolation as follow, second-order largrange-interpolation:
      idxc_interp=ic-1
      if(pos%x>xc(ic)) then
        idxp_interp= min(ic,  xe_interp-2)
      else
        idxp_interp= max(ic-1,xm_interp)
      endif

      idyc_interp=jc-1
      if(pos%y>yc(jc)) then
        idyp_interp=min(jc,nyc-1)
      else
        idyp_interp=max(jc-1,1)
      endif

      idzc_interp=kc-1
      if(pos%z>zc(kc)) then
        idzp_interp= min(kc,  ze_interp-2)
      else
        idzp_interp= max(kc-1,zm_interp)
      endif

      indxyz(1,pid)=idxc_interp
      indxyz(2,pid)=idxp_interp
      indxyz(3,pid)=idyc_interp
      indxyz(4,pid)=idyp_interp
      indxyz(5,pid)=idzc_interp
      indxyz(6,pid)=idzp_interp

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
  end subroutine PrepareDistribute

  !******************************************************************
  ! FinalDistribute
  !******************************************************************
  subroutine FinalDistribute()
    implicit none
    integer::pid
    
    if(GPrtcl_list%nlocal > 0) then
      deallocate(idxc,idyc,idzc,indxyz)
    endif
    DO pid=1,GPrtcl_list%nlocal
      GPrtcl_linVelOld(pid) =GPrtcl_LinVel(1,pid)
    ENDDO
  end subroutine FinalDistribute

  !******************************************************************
  ! clc_FpForce
  !******************************************************************
  subroutine clc_FpForce(ux,uy,uz,pressure,RatioYp_interp,RatioYc_interp)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,pressure
    type(real3),dimension(GPrtcl_list%mlocal),intent(in)::RatioYp_interp,RatioYc_interp

    ! locals
    type(real4)::pos
    type(real3)::FpVel,Drag
    real(RK)::cd,rep,diam,veldiff
    real(RK)::prx,pry,prz,SumXDir,SumYDir,SumZDir,DragCoe
    integer::idzp_interp,idxc_interp,idyc_interp,idzc_interp
    integer::id,jd,kd,pid,i,j,k,nlocal,idxp_interp,idyp_interp
    real(RK),dimension(0:2)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp

    nlocal=GPrtcl_list%nlocal
    !=====================================    Drag force     =====================================!
    ! [1] F.F. Abraham, Phys. Fluids 13 (1970) 2194–2195
    !       Functional dependence of drag coefficient of a sphere on Reynolds number.
    ! [2] Brown P. P., and Lawler D. F., J. Environ. Eng., 2003, 129(3): 222-231
    !       Sphere drag and settling velocity revisited
#define DragType1
    call update_halo(ux, mb1, hi_ux_interp)
    call update_halo(uz, mb1, hi_uz_interp)
    DragCoe=3.0_RK*pi*FluidDensity*xnu
    DO pid=1,nlocal
      pos = GPrtcl_posR(pid)
      diam= pos%w*2.0_RK 
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
      FpVel  = real3(SumXDir,SumYDir,SumZDir)
   
      GPrtcl_Vfluid(2,pid) =GPrtcl_Vfluid(1,pid)
      GPrtcl_Vfluid(1,pid) =FpVel
      FpVel = FpVel-GPrtcl_LinVel(1,pid)

      veldiff=sqrt(FpVel%x*FpVel%x+ FpVel%y*FpVel%y+ FpVel%z*FpVel%z)
      rep= veldiff *diam/xnu
#ifdef DragType1
      cd =(1.0_RK+0.11037527593819_RK*sqrt(rep))**2                                 ! Ref.[1]
#else
      cd = 1.0_RK+0.150_RK*(rep**0.681_RK)+ 0.0169583333_RK*rep*rep/(rep+8710.0_RK) ! Ref.[2]
#endif
      Drag= cd*(DragCoe*diam)*FpVel

      GPrtcl_FpForce_old(pid)=GPrtcl_FpForce(pid)
      GPrtcl_FpForce(pid) = Drag
    ENDDO
 
    ! Fixed particle part
    DO pid=1,mlocalFix
      pos = GPFix_PosR(pid)
      diam= pos%w*2.0_RK 
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
      FpVel=real3(SumXDir,SumYDir,SumZDir)
      GPFix_Vfluid(2,pid) =GPFix_Vfluid(1,pid)
      GPFix_Vfluid(1,pid) =FpVel
    
      veldiff=sqrt(FpVel%x*FpVel%x+ FpVel%y*FpVel%y+ FpVel%z*FpVel%z)
      rep= veldiff *diam/xnu
#ifdef DragType1
      cd =(1.0_RK+0.11037527593819_RK*sqrt(rep))**2                                 ! Ref.[1]
#else
      cd = 1.0_RK+0.150_RK*(rep**0.681_RK)+ 0.0169583333_RK*rep*rep/(rep+8710.0_RK) ! Ref.[2]
#endif
      Drag= cd*(DragCoe*diam)*FpVel

      GPFix_FpForce(pid) = Drag
    ENDDO 

    if(is_clc_FluidAcc) call clc_FluidAccForce(ux,uy,uz,RatioYp_interp,RatioYc_interp)
    if(is_clc_Lift)     call clc_LiftForce(ux,uy,uz,RatioYp_interp,RatioYc_interp)
    if(is_clc_Basset)   call clc_BassetForce(ux,uy,uz,RatioYp_interp,RatioYc_interp)
    if(is_clc_ViscousForce)     call clc_VisForce(ux,uy,uz,RatioYp_interp,RatioYc_interp)
    if(is_clc_PressureGradient) call clc_PrGradForce(pressure,RatioYp_interp,RatioYc_interp)

  end subroutine clc_FpForce

  !******************************************************************
  ! dist_FpForce 
  !******************************************************************
  subroutine distribute_FpForce()
    implicit none

    ! locals
    type(real4)::pos
    type(real3)::FpForce
    integer::ic,jc,kc,pid,itype,i,j,k,nFluidCell,idy2,idxs,idys,idzs,idxe,idye,idze
    real(RK)::dcell2,SearchR,ry,dist,volp,rad,massFH,fpx,fpy,fpz
    real(RK)::volt,dist2,distx,disty,distz,disty2,distz2,vcell,xi

    ! Refrence: An Euler–Lagrange strategy for simulating particle-laden flows, Eq.(52)
    !           J. Capecelatro, O. Desjardins/Journal of Computational Physics 238 (2013)
    FpForce_x=0.0_RK; FpForce_y=0.0_RK; FpForce_z=0.0_RK
    DO pid=1,GPrtcl_list%nlocal
      ic  = idxc(pid)
      jc  = idyc(pid)
      kc  = idzc(pid)
      itype = GPrtcl_pType(pid)
      massFH = 0.5_RK*DEMProperty%Prtcl_PureProp(itype)%MassOfFluid
      FpForce= massFH*(GPrtcl_linVel(1,pid)-GPrtcl_linVelOld(pid))/dt -GPrtcl_FpForce(pid)

      pos = GPrtcl_posR(pid);    rad = pos%w
      volp= DEMProperty%Prtcl_PureProp(itype)%Volume

      SearchR  = RatioSR * rad
      idxs= max(ceiling((pos%x-SearchR)*rdx), xmMin)
      idxe= min(ceiling((pos%x+SearchR)*rdx), xeMax)
      idzs= max(ceiling((pos%z-SearchR)*rdz), zmMin)
      idze= min(ceiling((pos%z+SearchR)*rdz), zeMax)

      ry = pos%y-SearchR
      idy2= 1
      do j=jc,1,-1
        if(ry >= yp(j)) then
         idy2 = j; exit
        endif      
      enddo
      idys = idy2

      ry  = pos%y+SearchR
      idy2= nyc
      do j=jc+1,nyp
        if(ry < yp(j)) then
         idy2 = j-1; exit
        endif      
      enddo
      idye = idy2
 
      volt= 0.0_RK; nFluidCell=0
      do k=idzs,idze
        distz = zc(k) - pos%z
        distz2= distz*distz
        do j=idys,idye
          vcell = VolCell(j)
          disty = yc(j) - pos%y
          disty2= disty*disty
          do i=idxs,idxe
            distx = xc(i) - pos%x
            dist  = RatioSR-sqrt(distx*distx +disty2 +distz2)/rad
            if(dist>0.0_RK) then
              nFluidCell=nFluidCell+1
              volt=volt+vcell*exp(-0.125_RK*dist*dist)
            endif
          enddo
        enddo
      enddo
      if(nFluidCell<27) then
        volt=0.0_RK
        dcell2=DeltaCell(jc)*DeltaCell(jc)
        do k= kc-1,kc+1
          distz = zc(k) - pos%z
          distz2= distz*distz
          do j=jc-1,jc+1
            vcell = VolCell(j)
            disty = yc(j) - pos%y
            disty2= disty*disty
            do i=ic-1,ic+1
              distx = xc(i) - pos%x
              dist2 = (distx*distx+disty2+distz2)/dcell2
              volt  = volt + vcell/(16.0_RK**dist2)
            enddo
          enddo
        enddo
        do k= kc-1,kc+1
          distz = zc(k) - pos%z
          distz2= distz*distz
          do j=jc-1,jc+1
            disty = yc(j) - pos%y
            disty2= disty*disty
            do i=ic-1,ic+1
              distx = xc(i) - pos%x
              dist2 = (distx*distx+disty2+distz2)/dcell2
              xi =  1.0_RK/((16.0_RK**dist2)*volt*FluidDensity)
              FpForce_x(i,j,k)=FpForce_x(i,j,k)+ xi*FpForce%x
              FpForce_y(i,j,k)=FpForce_y(i,j,k)+ xi*FpForce%y
              FpForce_z(i,j,k)=FpForce_z(i,j,k)+ xi*FpForce%z          
            enddo
          enddo
        enddo
        cycle
      endif

      do k=idzs,idze
        distz = zc(k) - pos%z
        distz2= distz*distz
        do j=idys,idye
          disty = yc(j) - pos%y
          disty2= disty*disty
          do i=idxs,idxe
            distx = xc(i) - pos%x
            dist  = RatioSR-sqrt(distx*distx +disty2 +distz2)/rad
            if(dist>0.0_RK) then
              xi = exp(-0.125_RK*dist*dist)/(volt*FluidDensity)
              FpForce_x(i,j,k)=FpForce_x(i,j,k)+ xi*FpForce%x
              FpForce_y(i,j,k)=FpForce_y(i,j,k)+ xi*FpForce%y
              FpForce_z(i,j,k)=FpForce_z(i,j,k)+ xi*FpForce%z
            endif
          enddo
        enddo
      enddo
    ENDDO

    ! fixed particle part
    DO pid=1,mlocalFix
      ic  = idxcFixed(pid)
      jc  = idycFixed(pid)
      kc  = idzcFixed(pid)
      fpx=  -GPFix_FpForce(pid)%x
      fpy=  -GPFix_FpForce(pid)%y
      fpz=  -GPFix_FpForce(pid)%z

      pos = GPFix_posR(pid);    rad = pos%w

      SearchR  = RatioSR * rad
      idxs= max(ceiling((pos%x-SearchR)*rdx), xmMin)
      idxe= min(ceiling((pos%x+SearchR)*rdx), xeMax)
      idzs= max(ceiling((pos%z-SearchR)*rdz), zmMin)
      idze= min(ceiling((pos%z+SearchR)*rdz), zeMax)

      ry = pos%y-SearchR
      idy2= 1
      do j=jc,1,-1
        if(ry >= yp(j)) then
         idy2 = j; exit
        endif      
      enddo
      idys = idy2

      ry  = pos%y+SearchR
      idy2= nyc
      do j=jc+1,nyp
        if(ry < yp(j)) then
         idy2 = j-1; exit
        endif      
      enddo
      idye = idy2
  
      volt      = GPFix_VolFluid(pid)
      nFluidCell= GPFix_nFluidCell(pid)
      IF(nFluidCell<27) THEN
        do k= kc-1,kc+1
          distz = zc(k) - pos%z
          distz2= distz*distz
          do j=jc-1,jc+1
            disty = yc(j) - pos%y
            disty2= disty*disty
            do i=ic-1,ic+1
              distx = xc(i) - pos%x
              dist2 = (distx*distx+disty2+distz2)/dcell2
              xi =  1.0_RK/((16.0_RK**dist2)*volt*FluidDensity)
              FpForce_x(i,j,k)=FpForce_x(i,j,k)+ xi*fpx
              FpForce_y(i,j,k)=FpForce_y(i,j,k)+ xi*fpy
              FpForce_z(i,j,k)=FpForce_z(i,j,k)+ xi*fpz          
            enddo
          enddo
        enddo
      ELSE
        do k=idzs,idze
          distz = zc(k) - pos%z
          distz2= distz*distz
          do j=idys,idye
            disty = yc(j) - pos%y
            disty2= disty*disty
            do i=idxs,idxe
              distx = xc(i) - pos%x
              dist  = RatioSR-sqrt(distx*distx +disty2 +distz2)/rad
              if(dist>0.0_RK) then
                xi = exp(-0.125_RK*dist*dist)/(volt*FluidDensity)
                FpForce_x(i,j,k)=FpForce_x(i,j,k)+ xi*fpx
                FpForce_y(i,j,k)=FpForce_y(i,j,k)+ xi*fpy
                FpForce_z(i,j,k)=FpForce_z(i,j,k)+ xi*fpz
              endif
            enddo
          enddo
        enddo
      ENDIF
    ENDDO

#ifdef slgjskdljflksdg
    do k=mb_dist%zmm,mb_dist%zpm
      do i=mb_dist%xmm,mb_dist%xpm
        FpForce_x(i,1,k)=FpForce_x(i,1,k)+FpForce_x(i,0,k); FpForce_x(i,0,k)=0.0_RK
        FpForce_y(i,1,k)=FpForce_y(i,1,k)+FpForce_y(i,0,k); FpForce_y(i,0,k)=0.0_RK
        FpForce_z(i,1,k)=FpForce_z(i,1,k)+FpForce_z(i,0,k); FpForce_z(i,0,k)=0.0_RK
        FpForce_x(i,nyc,k)=FpForce_x(i,nyc,k)+FpForce_x(i,nyp,k); FpForce_x(i,nyp,k)=0.0_RK
        FpForce_y(i,nyc,k)=FpForce_y(i,nyc,k)+FpForce_y(i,nyp,k); FpForce_y(i,nyp,k)=0.0_RK
        FpForce_z(i,nyc,k)=FpForce_z(i,nyc,k)+FpForce_z(i,nyp,k); FpForce_z(i,nyp,k)=0.0_RK
      enddo
    enddo
#endif

    call Gather_Halo_dist(FpForce_x)
    call Gather_Halo_dist(FpForce_y)
    call Gather_Halo_dist(FpForce_z)
    call SetBC_and_UpdateHalo_dist(FpForce_x)
    call SetBC_and_UpdateHalo_dist(FpForce_y)
    call SetBC_and_UpdateHalo_dist(FpForce_z)
  end subroutine distribute_FpForce

  !******************************************************************
  ! Gather_Halo_dist
  !******************************************************************
  subroutine Gather_Halo_dist(mat)
    implicit none
    real(RK),dimension(mb_dist%xmm:mb_dist%xpm,mb_dist%ymm:mb_dist%ypm,mb_dist%zmm:mb_dist%zpm),intent(inout)::mat

    ! locals
    integer:: xs1,xs2,ys1,ys2,zs1,zs2                ! index for halo sending
    integer,dimension(MPI_STATUS_SIZE,2) :: SRstatus
    integer:: i,j,k,icount,ierror,requests(2),ProcNgh(4)

    do j=1,4
      if(myProcNghBC(y_pencil,j)<0) then
        ProcNgh(j)=MPI_PROC_NULL
      else
      ProcNgh(j)=myProcNghBC(y_pencil,j)
      endif
    enddo

    ! step1: receive from xm_dir, and send to xp_dir
    xs1=y1end(1)+1;          xs2=y1end(1)+distwx
    ys1=y1start(2)-1;        ys2=y1end(2)+1    
    zs1=zm_dist;             zs2=ze_dist
    if(ProcNgh(3)==nrank) then  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      xmhalo = mat(xs1:xs2,ys1:ys2, zs1:zs2)
    else
      icount = (ys2-ys1+1)*(zs2-zs1+1)*distwx
      xmhalos= mat(xs1:xs2,ys1:ys2, zs1:zs2)
      call MPI_IRECV( xmhalo, icount,real_type,ProcNgh(4),1,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(xmhalos,icount,real_type,ProcNgh(3),1,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    endif
    if(ProcNgh(4)>=0) then
      do k=zm_dist, ze_dist
        do j=y1start(2), y1end(2)
          do i=y1start(1), y1start(1)+distwx-1
            mat(i,j,k)=mat(i,j,k)+xmhalo(i,j,k)
          enddo
        enddo
      enddo
    endif

    ! step2: receive from xp_dir, and send to xm_dir
    xs1=y1start(1)-distwx;   xs2=y1start(1)-1
    ys1=y1start(2)-1;        ys2=y1end(2)+1    
    zs1=zm_dist;            zs2=ze_dist
    if(ProcNgh(4)==nrank) then  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
      xphalo  = mat(xs1:xs2,ys1:ys2, zs1:zs2)
    else
      icount  = (ys2-ys1+1)*(zs2-zs1+1)*distwx
      xphalos = mat(xs1:xs2,ys1:ys2, zs1:zs2)
      call MPI_IRECV( xphalo, icount,real_type,ProcNgh(3),2,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(xphalos,icount,real_type,ProcNgh(4),2,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    endif
    if(ProcNgh(3)>=0) then
      do k=zm_dist, ze_dist
        do j=y1start(2), y1end(2)
          do i=y1end(1)-distwx+1, y1end(1)
            mat(i,j,k)=mat(i,j,k)+xphalo(i,j,k)
          enddo
        enddo
      enddo
    endif
    
    ! step3: receive from zm_dir, and send to zp_dir
    xs1=xm_dist;            xs2=xe_dist
    ys1=y1start(2)-1;        ys2=y1end(2)+1
    zs1=y1end(3)+1;          zs2=y1end(3)+distwz
    if(ProcNgh(1)==nrank) then  ! neighbour is nrank itself, and ProcNgh(2)==nrank also !!!
      zmhalo = mat(xs1:xs2,ys1:ys2, zs1:zs2)
    else
      icount = (xs2-xs1+1)*(ys2-ys1+1)*distwz
      call MPI_IRECV( zmhalo,          icount,real_type,ProcNgh(2),3,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(mat(xs1,ys1,zs1),icount,real_type,ProcNgh(1),3,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    endif
    if(ProcNgh(2)>=0) then
      do k=y1start(3), y1start(3)+distwz-1
        do j=y1start(2), y1end(2)
          do i=xm_dist, xe_dist
            mat(i,j,k)=mat(i,j,k)+zmhalo(i,j,k)
          enddo
        enddo
      enddo
    endif

    ! step4: receive from zp_dir, and send to zm_dir
    xs1=xm_dist;           xs2=xe_dist
    ys1=y1start(2)-1;       ys2=y1end(2)+1
    zs1=y1start(3)-distwz;  zs2=y1start(3)-1
    if(ProcNgh(2)==nrank) then  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      zphalo = mat(xs1:xs2,ys1:ys2, zs1:zs2)
    else
      icount = (xs2-xs1+1)*(ys2-ys1+1)*distwz
      call MPI_IRECV( zphalo,          icount,real_type,ProcNgh(1),4,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(mat(xs1,ys1,zs1),icount,real_type,ProcNgh(2),4,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    endif
    if(ProcNgh(1)>=0) then
      do k=y1end(3)-distwz+1, y1end(3)
        do j=y1start(2), y1end(2)
          do i=xm_dist, xe_dist
            mat(i,j,k)=mat(i,j,k)+zphalo(i,j,k)
          enddo
        enddo
      enddo
    endif

  end subroutine Gather_Halo_dist

  !******************************************************************
  ! SetBC_and_UpdateHalo_dist
  !******************************************************************
  subroutine SetBC_and_UpdateHalo_dist(mat)
    implicit none
    real(RK),dimension(mb_dist%xmm:mb_dist%xpm,mb_dist%ymm:mb_dist%ypm,mb_dist%zmm:mb_dist%zpm),intent(inout)::mat

    ! locals
    integer:: i,j,k,ProcNgh(4)

    do j=1,4
      if(myProcNghBC(y_pencil,j)<0) then
        ProcNgh(j)=MPI_PROC_NULL
      else
      ProcNgh(j)=myProcNghBC(y_pencil,j)
      endif
    enddo

    ! step5: set Bc for mat, Here 0.0_RK-gradient Bcs for all non-periodic directions
    ! yp-dir, and ym-dir
    do k=y1start(3),y1end(3)
      do i=y1start(1),y1end(1)
        mat(i, nyp, k) = mat(i, nyc, k)
        mat(i, 0,   k) = mat(i, 1,   k)
      enddo
    enddo   

    ! xp-dir, and  xm-dir
    if(ProcNgh(3)<0) then
      do k=y1start(3),y1end(3)
        do j=y1start(2),y1end(2)
          mat(nxp,j,k) = mat(nxc,j,k)
        enddo
      enddo
    endif
    if(ProcNgh(4)<0) then
      do k=y1start(3),y1end(3)
        do j=y1start(2),y1end(2)
          mat(0,j,k) = mat(1,j,k)
        enddo
      enddo
    endif

    ! zp-dir, and zm-dir
    if(ProcNgh(1)<0) then
      do j=y1start(2),y1end(2)
        do i=y1start(1),y1end(1)
          mat(i,j,nzp) = mat(i,j,nzc)
        enddo
      enddo
    endif
    if(ProcNgh(2)<0) then
      do j=y1start(2),y1end(2)
        do i=y1start(1),y1end(1)
          mat(i,j,0) = mat(i,j,1)
        enddo
      enddo
    endif
   
    ! step6: update halo
    call update_halo(mat, mb_dist, hi_dist)
  end subroutine SetBC_and_UpdateHalo_dist

end module cd_FpForce
