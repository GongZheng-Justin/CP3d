module m_FlowCase
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_decomp2d
  use m_Parameters
  use iso_c_binding
  use m_MeshAndMetries
  use m_Variables,only:mb1
  use m_Tools,only:CalcUxAver
#ifdef CFDLPT_TwoWay
  use m_Variables,only:FpForce_x,FpForce_y,FpForce_z
#endif
  implicit none
  private
  include "fftw3.f03"
#define HighOrderGradStat
#define SAVE_SINGLE_Spec2D

  ! statistics variabls
  integer:: nfstime
  real(RK):: PrGradsum
  real(RK),allocatable,dimension(:,:)::SumStat
#ifdef HighOrderGradStat 
  real(RK),allocatable,dimension(:,:)::SumGrad
#endif

  ! SpectraOptions
  integer:: ivSpec,jForLCS(2)
  logical:: clcSpectra1D,clcSpectra2D
  
  ! Spectra variables
  integer:: nxh,nxhp,nzh,nzhp,nSpectime
  type(C_PTR)::fft_plan_x,fft_plan_z,fft_plan_z2
  type(decomp_info),allocatable::decomp_xhzf,decomp_xhzh
  real(RK),allocatable,dimension(:,:,:,:)::EnergySpec2D
  real(RK),allocatable,dimension(:,:,:)::EnergySpecX,EnergySpecZ
  
  public:: InitVelocity, Update_uy_ym, InitStatVar, clcStat
contains
#define iSpec1DUU  1
#define iSpec1DVV  2
#define iSpec1DWW  3
#define iSpec1DPP  4
#define iSpec1DUV  5
#define iSpec1DUV2 6
#define iLCSR1DUU  7
#define iLCSI1DUU  8
#define iLCSR1DUU2 9
#define iLCSR1DVV  10
#define iLCSI1DVV  11
#define iLCSR1DVV2 12
#define iLCSR1DWW  13
#define iLCSI1DWW  14
#define iLCSR1DWW2 15
#define iLCSR1DPP  16
#define iLCSI1DPP  17
#define iLCSR1DPP2 18
#define iSpec1DCC  19
#define iSpec1DUC  20
#define iSpec1DVC  21
#define iLCSR1DCC  22
#define iLCSI1DCC  23
#define iLCSR1DCC2 24

#define iSpec2DUU  1
#define iSpec2DVV  2
#define iSpec2DWW  3
#define iSpec2DPP  4
#define iSpec2DUV  5
#define iLCSR2DUU  6
#define iLCSR2DVV  7
#define iLCSR2DWW  8
#define iLCSR2DPP  9
#define iSpec2DCC  10
#define iLCSR2DCC  11

#if defined(CFDLPT_TwoWay)
#define NCHASTAT      45
#define NEnergySpec1D 18
#define NEnergySpec2D 9
#elif defined(ScalarFlow)
#define NCHASTAT      40
#define NEnergySpec1D 24
#define NEnergySpec2D 11
#else
#define NCHASTAT      35
#define NEnergySpec1D 18
#define NEnergySpec2D 9
#endif

#if defined(HighOrderGradStat)
#define nGradStat     98
#endif

  !******************************************************************
  ! InitVelocity
  !******************************************************************
  subroutine InitVelocity(ux,uy,uz,Deviation)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
  
    ! locals
    integer :: ii,code,i,j,k,m1,m2
    real(RK):: retau_guass,utau_guass,height,rem,wx,wz,xlxPlus,zlzPlus
    real(RK):: xplus,yplus,zplus,yct,ybar,xp,zp,ratiot,uzmean(nyc),uzmeanR(nyc)

    height=yly
    if(FlowType==FT_CH)height=half*yly

    ux=zero; uy=zero; uz=zero
    rem = ubulk * height / xnu
    retau_guass = 0.1538_RK*rem**0.887741_RK
    utau_guass   = retau_guass*xnu/height
    if(nrank==0)print*,'************** retau_gauss=',retau_guass
    if(nrank==0)print*,'************** utau_gauss= ',utau_guass

    call system_clock(count=code); !code=0
    call random_seed(size = ii)
    call random_seed(put = code+63946*(/ (i - 1, i = 1, ii) /))
    call random_number(Deviation)
    Deviation= 0.2_RK* Deviation + 0.9_RK ! [0.8, 1.2]

    !modulation of the random noise + initial velocity profile
    uzmean=zero
    wx=twopi/500.0_RK; wz=twopi/200.0_RK
    xlxPlus=xlx*utau_guass/xnu;   zlzPlus=zlz*utau_guass/xnu;
    m1=floor(xlxPlus*wx/twopi)+1; wx=real(m1,RK)*twopi/xlxPlus
    m2=floor(zlzPlus*wz/twopi)+1; wz=real(m2,RK)*twopi/zlzPlus
    do j=y1start(2),y1end(2)
      yct = height-abs(height-yc(j)) 
      ybar= yct/height; yplus=utau_guass*yct/xnu
      do k=y1start(3),y1end(3)
        zp   =real(k-1,kind=RK)*dz+dz*half
        zplus=utau_guass*zp/xnu
        do i=y1start(1),y1end(1)
          xp   =real(i-1,kind=RK)*dx+dx*half
          xplus=utau_guass*xp/xnu
          !ux(i,j,k) = 0.0052_RK*ubulk*yplus*exp(-yplus*yplus/1800.0_RK)*cos(wz*zplus)*Deviation(i,j,k) ! original expression
          !uz(i,j,k) = 0.0050_RK*ubulk*yplus*exp(-yplus*yplus/1800.0_RK)*sin(wx*xplus)*Deviation(i,j,k) ! original expression
          ux(i,j,k) = ubulk*ybar*exp(-4.5_RK*ybar*ybar)*cos(wz*zplus)*Deviation(i,j,k)
          uz(i,j,k) = ubulk*ybar*exp(-4.5_RK*ybar*ybar)*sin(wx*xplus)*Deviation(i,j,k)
          ux(i,j,k) = ux(i,j,k)+ three*ubulk*(ybar-half*ybar*ybar)
          uzmean(j) = uzmean(j)+ uz(i,j,k)
        enddo
      enddo
      uzmean(j)=uzmean(j)/real(nxc*nzc,RK)
    enddo
    ratiot=ubulk/CalcUxAver(ux)
    call MPI_Bcast(ratiot,1,real_type,0,MPI_COMM_WORLD,code)
    call MPI_ALLREDUCE(uzmean,uzmeanR,nyc,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    ux= ux*ratiot -uCRF
    do j=y1start(2),y1end(2)
      uz(:,j,:)=uz(:,j,:)-uzmeanR(j)
    enddo   
  end subroutine InitVelocity

  !******************************************************************
  ! Update_uy_ym
  !******************************************************************   
  subroutine Update_uy_ym(uy_ym, duy_ym)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(out):: uy_ym
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(inout):: duy_ym
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = zero
    duy_ym= uy_ym - duy_ym
  end subroutine Update_uy_ym

  !******************************************************************
  ! InitStatVar
  !******************************************************************
  subroutine InitStatVar(ChannelPrm)
    implicit none 
    character(*),intent(in)::ChannelPrm

    ! locals
    character(len=128)::filename
    integer::ierror,nUnit,plan_type,nySpec2D
#ifdef OverWriteFFT
    real(RK),dimension(:,:,:),allocatable::Arr
#else
    real(RK),dimension(:,:,:),allocatable::Arr1,Arr2
#endif
    type(fftw_iodim),dimension(1)::iodim,iodim_howmany
    NAMELIST/SpectraOptions/clcSpectra1D,clcSpectra2D,ivSpec,jForLCS
    
    open(newunit=nUnit, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar", "Cannot open file: "//trim(ChannelPrm))
    read(nUnit, nml=SpectraOptions)
    close(nUnit,IOSTAT=ierror)
    
    if(nrank==0) then
      write(MainLog%nUnit, nml=SpectraOptions)
      if(mod(saveStat,ivstats)/=0 )  call MainLog%CheckForError(ErrT_Abort,"InitStatVar","ivstats wrong !!!")
      if(clcSpectra1D .and. (mod(saveStat,ivSpec)/=0 .or. mod(ivSpec,ivstats)/=0 )) then
        call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","ivSpec wrong !!!")
      endif
      if(IsUxConst)then
        write(filename,'(A,I10.10)')trim(ResultsDir)//"PrGrad",ilast
        open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
        if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitCAStatistics","Cannot open file: "//trim(filename))
        close(nUnit,IOSTAT=ierror)
      endif
    endif
    allocate(SumStat(NCHASTAT,nyp),Stat=ierror)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar: ","Allocation failed 1")  
    nfstime=0; nSpectime=0; SumStat=zero; PrGradsum=zero
#ifdef HighOrderGradStat
    allocate(SumGrad(nGradStat,nyp),Stat=ierror)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar: ","Allocation failed 2")  
    SumGrad=zero
#endif

    if(FFTW_plan_type == 1) then
      plan_type=FFTW_MEASURE
    else
      plan_type=FFTW_ESTIMATE
    endif    
    nxh=nxc/2; nxhp=nxh+1
    nzh=nzc/2; nzhp=nzh+1
    IF(clcSpectra1D) THEN
      iodim(1)%n  = x1size(1)
      iodim(1)%is = 1
      iodim(1)%os = 1
      iodim_howmany(1)%n  = x1size(2)*x1size(3)
      iodim_howmany(1)%is = x1size(1)
      iodim_howmany(1)%os = x1size(1)
      allocate(EnergySpecX(nxhp,x1size(2),NEnergySpec1D));EnergySpecX=zero
#ifdef OverWriteFFT
      allocate(Arr(x1size(1),x1size(2),x1size(3)))
      fft_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr,Arr,[FFTW_R2HC],plan_type)
      deallocate(Arr)
#else
      allocate(Arr1(x1size(1),x1size(2),x1size(3)),Arr2(x1size(1),x1size(2),x1size(3)))
      fft_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,[FFTW_R2HC],plan_type)
      deallocate(Arr1,Arr2)
#endif

      iodim(1)%n  = z1size(3)
      iodim(1)%is = z1size(1)*z1size(2)
      iodim(1)%os = z1size(1)*z1size(2)
      iodim_howmany(1)%n  = z1size(1)*z1size(2)
      iodim_howmany(1)%is = 1
      iodim_howmany(1)%os = 1    
      allocate(EnergySpecZ(nzhp,z1size(2),NEnergySpec1D));EnergySpecZ=zero
#ifdef OverWriteFFT
      allocate(Arr(z1size(1),z1size(2),z1size(3)))
      fft_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr,Arr,[FFTW_R2HC],plan_type)
      deallocate(Arr)
#else
      allocate(Arr1(z1size(1),z1size(2),z1size(3)),Arr2(z1size(1),z1size(2),z1size(3)))
      fft_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,[FFTW_R2HC],plan_type)
      deallocate(Arr1,Arr2)
#endif
    ENDIF
    IF(clcSpectra2D) THEN
      allocate(decomp_xhzf,decomp_xhzh)
      Block
        logical,dimension(6)::initializeIn
        initializeIn=.false.; initializeIn(5:6)=.true.
        call decomp_info_init(nxhp,nyc,nzc, decomp_xhzf,initialize=initializeIn)
        initializeIn=.false.; initializeIn(5)=.true.
        call decomp_info_init(nxhp,nyc,nzhp,decomp_xhzh,initialize=initializeIn)
      end Block
    
      iodim(1)%n  = decomp_xhzf%z2sz(3)
      iodim(1)%is = decomp_xhzf%z2sz(1)*decomp_xhzf%z2sz(2)
      iodim(1)%os = decomp_xhzf%z2sz(1)*decomp_xhzf%z2sz(2)
      iodim_howmany(1)%n  = decomp_xhzf%z2sz(1)*decomp_xhzf%z2sz(2)
      iodim_howmany(1)%is = 1
      iodim_howmany(1)%os = 1
#ifdef OverWriteFFT
      allocate(Arr(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
      fft_plan_z2 = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr,Arr,[FFTW_R2HC],plan_type)
      deallocate(Arr)
#else
      allocate(Arr1(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)),&
               Arr2(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)) )
      fft_plan_z2 = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,[FFTW_R2HC],plan_type)
      deallocate(Arr1,Arr2)
#endif
      if(FlowType==FT_CH) then
        nySpec2D=nyc/2
      else
        nySpec2D=nyc
      endif
      allocate(EnergySpec2D(decomp_xhzh%y2sz(1),nySpec2D,decomp_xhzh%y2sz(3),NEnergySpec2D),Stat=ierror)
      if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar","Allocation failed For Spectra2D")
      EnergySpec2D=zero
    ENDIF
  end subroutine InitStatVar

  !******************************************************************
  ! clcStat
  !******************************************************************
#ifdef ScalarFlow
  subroutine clcStat(ux,uy,uz,pr,scalar,ArrTemp1,ArrTemp2)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pr,scalar
#else
  subroutine clcStat(ux,uy,uz,pr,ArrTemp1,ArrTemp2)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pr
#endif
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::ArrTemp1,ArrTemp2
   
    ! locals
    character(len=128)::filename
    integer(kind=8)::disp,disp_inc
#ifndef OverWriteFFT
    real(RK),dimension(:,:,:),allocatable::ArrFFT
#endif
    real(RK),dimension(:,:,:),allocatable::arrx1,arrx2,arrz1,arrz2
    real(RK),allocatable,dimension(:,:,:)::EnergySpecXR,EnergySpecZR
    integer::ic,jc,kc,im,jm,km,ip,jp,kp,ierror,is,ks,iu,ku,it,jt,kt,nrankX,nrankZ,nUnit
    real(RK)::uxloc,uyloc,uzloc,prloc,uxCell,uyCell,uzCell,prloc2,inxz,infstime,cac,caj,cacU,dudyU,uy_xc,uy_xp
    real(RK)::dudxx,dudyy,dudzz,dvdxx,dvdyy,dvdzz,dwdxx,dwdyy,dwdzz,SumStatR(NCHASTAT,nyp),SumVec(NCHASTAT),rdxt,rdzt
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vor_x,vor_y,vor_z,InterpY1,InterpY2,dudzC,dvdzC,dudxC,dvdxU
#ifdef HighOrderGradStat
    real(RK)::dpdx,dpdy,dpdz,dpdxx,dpdyy,dpdzz,SumGradR(nGradStat,nyp),SumVec2(nGradStat)
#endif
#ifdef CFDLPT_TwoWay
    real(RK)::uyCells,uyCellm,uyCellp
#endif
#ifdef ScalarFlow
    real(RK)::scloc
#endif

    rdxt=rdx/12.0_RK
    rdzt=rdz/12.0_RK
    inxz=one/(real(nxc,RK)*real(nzc,RK))
    DO jc=1,nyc
      jp=jc+1; jm=jc-1
      InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1
      cac=rdyc(jc); cacU=rdyc(jp); caj=rdyp(jc); SumVec=zero
#ifdef HighOrderGradStat
      SumVec2=zero
#endif
      do kc=y1start(3),y1end(3)
        ks=kc-2;km=kc-1;kp=kc+1;ku=kc+2
        do ic=y1start(1),y1end(1)
          is=ic-2;im=ic-1;ip=ic+1;iu=ic+2

          uxloc= ux(ic,jc,kc)+uCRF
          uyloc= uy(ic,jc,kc)
          uzloc= uz(ic,jc,kc)
          prloc= pr(ic,jc,kc)
          uxCell= InterpCoe1*ux(im,jc,kc) +InterpCoe2*ux(ic,jc,kc) +InterpCoe3*ux(ip,jc,kc) +InterpCoe4*ux(iu,jc,kc)+uCRF
          uyCell= (uyloc+ uy(ic,jp,kc))*half
          uy_xc = InterpCoe1*uy(is,jc,kc) +InterpCoe2*uy(im,jc,kc) +InterpCoe3*uy(ic,jc,kc) +InterpCoe4*uy(ip,jc,kc)
          uy_xp = InterpCoe1*uy(is,jp,kc) +InterpCoe2*uy(im,jp,kc) +InterpCoe3*uy(ic,jp,kc) +InterpCoe4*uy(ip,jp,kc)
          uzCell= InterpCoe1*uz(ic,jc,km) +InterpCoe2*uz(ic,jc,kc) +InterpCoe3*uz(ic,jc,kp) +InterpCoe4*uz(ic,jc,ku)
          prloc2= InterpY1*(InterpCoe1*pr(is,jm,kc) +InterpCoe2*pr(im,jm,kc) + &
                            InterpCoe3*pr(ic,jm,kc) +InterpCoe4*pr(ip,jm,kc))+ &
                  InterpY2*(InterpCoe1*pr(is,jc,kc) +InterpCoe2*pr(im,jc,kc) + &
                            InterpCoe3*pr(ic,jc,kc) +InterpCoe4*pr(ip,jc,kc))  !2021-12-14
 
          dudx= dxCoe1*ux(im,jc,kc) +dxCoe2*ux(ic,jc,kc) +dxCoe3*ux(ip,jc,kc) +dxCoe4*ux(iu,jc,kc)
          dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
          dudz= dzCoe1*ux(ic,jc,ks) +dzCoe2*ux(ic,jc,km) +dzCoe3*ux(ic,jc,kc) +dzCoe4*ux(ic,jc,kp)
          dvdx= dxCoe1*uy(is,jc,kc) +dxCoe2*uy(im,jc,kc) +dxCoe3*uy(ic,jc,kc) +dxCoe4*uy(ip,jc,kc)
          dvdy= (uy(ic,jp,kc)-uyloc)*caj
          dvdz= dzCoe1*uy(ic,jc,ks) +dzCoe2*uy(ic,jc,km) +dzCoe3*uy(ic,jc,kc) +dzCoe4*uy(ic,jc,kp)
          dwdx= dxCoe1*uz(is,jc,kc) +dxCoe2*uz(im,jc,kc) +dxCoe3*uz(ic,jc,kc) +dxCoe4*uz(ip,jc,kc)
          dwdy= (uzloc-uz(ic,jm,kc))*cac
          dwdz= dzCoe1*uz(ic,jc,km) +dzCoe2*uz(ic,jc,kc) +dzCoe3*uz(ic,jc,kp) +dzCoe4*uz(ic,jc,ku)

          dudxx= dxxCoe1*ux(is,jc,kc) +dxxCoe2*ux(im,jc,kc) +dxxCoe3*ux(ic,jc,kc) +dxxCoe4*ux(ip,jc,kc) +dxxCoe5*ux(iu,jc,kc)
          dudyy= ap2c(jc)*ux(ic,jp,kc)+ac2c(jc)*ux(ic,jc,kc)+am2c(jc)*ux(ic,jm,kc)
          dudzz= dzzCoe1*ux(ic,jc,ks) +dzzCoe2*ux(ic,jc,km) +dzzCoe3*ux(ic,jc,kc) +dzzCoe4*ux(ic,jc,kp) +dzzCoe5*ux(ic,jc,ku)
          dvdxx= dxxCoe1*uy(is,jc,kc) +dxxCoe2*uy(im,jc,kc) +dxxCoe3*uy(ic,jc,kc) +dxxCoe4*uy(ip,jc,kc) +dxxCoe5*uy(iu,jc,kc)
          dvdyy= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          dvdzz= dzzCoe1*uy(ic,jc,ks) +dzzCoe2*uy(ic,jc,km) +dzzCoe3*uy(ic,jc,kc) +dzzCoe4*uy(ic,jc,kp) +dzzCoe5*uy(ic,jc,ku)
          dwdxx= dxxCoe1*uz(is,jc,kc) +dxxCoe2*uz(im,jc,kc) +dxxCoe3*uz(ic,jc,kc) +dxxCoe4*uz(ip,jc,kc) +dxxCoe5*uz(iu,jc,kc)
          dwdyy= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)
          dwdzz= dzzCoe1*uz(ic,jc,ks) +dzzCoe2*uz(ic,jc,km) +dzzCoe3*uz(ic,jc,kc) +dzzCoe4*uz(ic,jc,kp) +dzzCoe5*uz(ic,jc,ku)
          vor_x= dwdy-dvdz
          vor_y= dudz-dwdx
          vor_z= dvdx-dudy

          dudxC= (-ux(iu,jc,kc)+8.0_RK*ux(ip,jc,kc)-8.0_RK*ux(im,jc,kc)+ux(is,jc,kc))*rdxt
          dvdxU= dxCoe1*uy(is,jp,kc) +dxCoe2*uy(im,jp,kc) +dxCoe3*uy(ic,jp,kc) +dxCoe4*uy(ip,jp,kc)
          dudyU= (ux(ic,jp,kc)-ux(ic,jc,kc))*cacU
          dudzC= (InterpY1*(-ux(ic,jm,ku)+8.0_RK*ux(ic,jm,kp)-8.0_RK*ux(ic,jm,km)+ux(ic,jm,ks) ) +&
                  InterpY2*(-ux(ic,jc,ku)+8.0_RK*ux(ic,jc,kp)-8.0_RK*ux(ic,jc,km)+ux(ic,jc,ks) ))*rdzt !2021-12-14
          dvdzC=((-uy(is,jc,ku)+8.0_RK*uy(is,jc,kp)-8.0_RK*uy(is,jc,km)+uy(is,jc,ks))*InterpCoe1+ &
                 (-uy(im,jc,ku)+8.0_RK*uy(im,jc,kp)-8.0_RK*uy(im,jc,km)+uy(im,jc,ks))*InterpCoe2+ &
                 (-uy(ic,jc,ku)+8.0_RK*uy(ic,jc,kp)-8.0_RK*uy(ic,jc,km)+uy(ic,jc,ks))*InterpCoe3+ &
                 (-uy(ip,jc,ku)+8.0_RK*uy(ip,jc,kp)-8.0_RK*uy(ip,jc,km)+uy(ip,jc,ks))*InterpCoe4)*rdzt !2021-12-14

          SumVec( 1)=SumVec( 1)+ uxloc
          SumVec( 2)=SumVec( 2)+ uyloc                     ! yp !
          SumVec( 3)=SumVec( 3)+ uzloc
          SumVec( 4)=SumVec( 4)+ prloc
          SumVec( 5)=SumVec( 5)+ uxloc*uxloc
          SumVec( 6)=SumVec( 6)+ uyloc*uyloc               ! yp !  
          SumVec( 7)=SumVec( 7)+ uzloc*uzloc
          SumVec( 8)=SumVec( 8)+ prloc*prloc
          SumVec( 9)=SumVec( 9)+ uxCell*uyCell
          SumVec(10)=SumVec(10)+ uyCell*uzCell
          SumVec(11)=SumVec(11)+ uxCell*uzCell
          SumVec(12)=SumVec(12)+ uxCell*prloc
          SumVec(13)=SumVec(13)+ uyCell*prloc
          SumVec(14)=SumVec(14)+ uzCell*prloc
          SumVec(15)=SumVec(15)+ uxCell*uxCell*uyCell 
          SumVec(16)=SumVec(16)+ uyloc *uyloc *uyloc       ! yp !
          SumVec(17)=SumVec(17)+ uzCell*uzCell*uyCell
          SumVec(18)=SumVec(18)+ uxCell*uyCell*uyCell
          SumVec(19)=SumVec(19)+ uxloc*uxloc*uxloc
          SumVec(20)=SumVec(20)+ uzloc*uzloc*uzloc
          SumVec(21)=SumVec(21)+ uxloc*uxloc*uxloc*uxloc
          SumVec(22)=SumVec(22)+ uyloc*uyloc*uyloc*uyloc   ! yp !
          SumVec(23)=SumVec(23)+ uzloc*uzloc*uzloc*uzloc
          SumVec(24)=SumVec(24)+ prloc*dudx
          SumVec(25)=SumVec(25)+ prloc*dvdy
          SumVec(26)=SumVec(26)+ prloc*dwdz
          SumVec(27)=SumVec(27)+ prloc2*(dudy+dvdx)        ! yp !
          SumVec(28)=SumVec(28)+ uxloc*(dudxx+dudyy+dudzz)
          SumVec(29)=SumVec(29)+ uyloc*(dvdxx+dvdyy+dvdzz) ! yp !
          SumVec(30)=SumVec(30)+ uzloc*(dwdxx+dwdyy+dwdzz)
          SumVec(31)=SumVec(31)+ (dudxC*(dvdxU+dvdx)+ (dudyU+dudy)*(uy_xp-uy_xc)*caj)*half !2021-12-14
          SumVec(32)=SumVec(32)+ dudzC*dvdzC               ! yp !
          SumVec(33)=SumVec(33)+ vor_x*vor_x               ! yp !
          SumVec(34)=SumVec(34)+ vor_y*vor_y
          SumVec(35)=SumVec(35)+ vor_z*vor_z               ! yp !
#ifdef CFDLPT_TwoWay
          SumVec(36)=SumVec(36)+ FpForce_x(ic,jc,kc)
          SumVec(37)=SumVec(37)+ FpForce_y(ic,jc,kc)       ! yp !
          SumVec(38)=SumVec(38)+ FpForce_z(ic,jc,kc)
          SumVec(39)=SumVec(39)+ FpForce_x(ic,jc,kc)*FpForce_x(ic,jc,kc)
          SumVec(40)=SumVec(40)+ FpForce_y(ic,jc,kc)*FpForce_y(ic,jc,kc)  ! yp !
          SumVec(41)=SumVec(41)+ FpForce_z(ic,jc,kc)*FpForce_z(ic,jc,kc)
                    
          SumVec(42)=SumVec(42)+ FpForce_x(ic,jc,kc)*uxloc
          SumVec(43)=SumVec(43)+ FpForce_y(ic,jc,kc)*uyloc ! yp !
          SumVec(44)=SumVec(44)+ FpForce_z(ic,jc,kc)*uzloc
          uyCells=(uy(is,jc,kc)+ uy(is,jp,kc))*half
          uyCellm=(uy(im,jc,kc)+ uy(im,jp,kc))*half
          uyCellp=(uy(ip,jc,kc)+ uy(ip,jp,kc))*half
          SumVec(45)=SumVec(45)+ half*(FpForce_y(ic,jc,kc)+FpForce_y(ic,jp,kc))*uxCell+ &
                     FpForce_x(ic,jc,kc)*(InterpCoe1*uyCells +InterpCoe2*uyCellm +InterpCoe3*uyCell +InterpCoe4*uyCellp)
#endif
#ifdef ScalarFlow
          scloc=scalar(ic,jc,kc)
          SumVec(36)=SumVec(36)+ scloc
          SumVec(37)=SumVec(37)+ scloc*scloc
          SumVec(38)=SumVec(38)+ scloc*uxCell
          SumVec(39)=SumVec(39)+ scloc*uyCell
          SumVec(40)=SumVec(40)+ scloc*uzCell
#endif

#ifdef HighOrderGradStat
          dpdx=dxCoe1*pr(is,jc,kc) +dxCoe2*pr(im,jc,kc) +dxCoe3*pr(ic,jc,kc) +dxCoe4*pr(ip,jc,kc)
          dpdy=(pr(ic,jc,kc)-pr(ic,jm,kc))*cac          ! yp ! 
          dpdz=dzCoe1*pr(ic,jc,ks) +dzCoe2*pr(ic,jc,km) +dzCoe3*pr(ic,jc,kc) +dzCoe4*pr(ic,jc,kp)
          dpdxx= dxxCoe1*pr(is,jc,kc) +dxxCoe2*pr(im,jc,kc) +dxxCoe3*pr(ic,jc,kc) +dxxCoe4*pr(ip,jc,kc) +dxxCoe5*pr(iu,jc,kc)
          dpdyy= ap2c(jc)*pr(ic,jp,kc)+ac2c(jc)*pr(ic,jc,kc)+am2c(jc)*pr(ic,jm,kc)
          dpdzz= dzzCoe1*pr(ic,jc,ks) +dzzCoe2*pr(ic,jc,km) +dzzCoe3*pr(ic,jc,kc) +dzzCoe4*pr(ic,jc,kp) +dzzCoe5*pr(ic,jc,ku)
          
          SumVec2( 1)=SumVec2( 1)+ prloc*prloc*prloc
          SumVec2( 2)=SumVec2( 2)+ prloc*prloc*prloc*prloc
          
          SumVec2( 3)=SumVec2( 3)+ dudx
          SumVec2( 4)=SumVec2( 4)+ dudx*dudx
          SumVec2( 5)=SumVec2( 5)+ dudx*dudx*dudx
          SumVec2( 6)=SumVec2( 6)+ dudx*dudx*dudx*dudx
          SumVec2( 7)=SumVec2( 7)+ dudy                 ! yp !
          SumVec2( 8)=SumVec2( 8)+ dudy*dudy            ! yp !
          SumVec2( 9)=SumVec2( 9)+ dudy*dudy*dudy       ! yp !
          SumVec2(10)=SumVec2(10)+ dudy*dudy*dudy*dudy  ! yp !
          SumVec2(11)=SumVec2(11)+ dudz
          SumVec2(12)=SumVec2(12)+ dudz*dudz
          SumVec2(13)=SumVec2(13)+ dudz*dudz*dudz
          SumVec2(14)=SumVec2(14)+ dudz*dudz*dudz*dudz

          SumVec2(15)=SumVec2(15)+ dvdx                 ! yp !          
          SumVec2(16)=SumVec2(16)+ dvdx*dvdx            ! yp !
          SumVec2(17)=SumVec2(17)+ dvdx*dvdx*dvdx       ! yp ! 
          SumVec2(18)=SumVec2(18)+ dvdx*dvdx*dvdx*dvdx  ! yp !
          SumVec2(19)=SumVec2(19)+ dvdy                 
          SumVec2(20)=SumVec2(20)+ dvdy*dvdy          
          SumVec2(21)=SumVec2(21)+ dvdy*dvdy*dvdy     
          SumVec2(22)=SumVec2(22)+ dvdy*dvdy*dvdy*dvdy  
          SumVec2(23)=SumVec2(23)+ dvdz                 ! yp !
          SumVec2(24)=SumVec2(24)+ dvdz*dvdz            ! yp !
          SumVec2(25)=SumVec2(25)+ dvdz*dvdz*dvdz       ! yp !
          SumVec2(26)=SumVec2(26)+ dvdz*dvdz*dvdz*dvdz  ! yp !
          
          SumVec2(27)=SumVec2(27)+ dwdx
          SumVec2(28)=SumVec2(28)+ dwdx*dwdx
          SumVec2(29)=SumVec2(29)+ dwdx*dwdx*dwdx
          SumVec2(30)=SumVec2(30)+ dwdx*dwdx*dwdx*dwdx
          SumVec2(31)=SumVec2(31)+ dwdy                 ! yp !
          SumVec2(32)=SumVec2(32)+ dwdy*dwdy            ! yp !
          SumVec2(33)=SumVec2(33)+ dwdy*dwdy*dwdy       ! yp !
          SumVec2(34)=SumVec2(34)+ dwdy*dwdy*dwdy*dwdy  ! yp !
          SumVec2(35)=SumVec2(35)+ dwdz
          SumVec2(36)=SumVec2(36)+ dwdz*dwdz
          SumVec2(37)=SumVec2(37)+ dwdz*dwdz*dwdz
          SumVec2(38)=SumVec2(38)+ dwdz*dwdz*dwdz*dwdz
          
          SumVec2(39)=SumVec2(39)+ dpdx
          SumVec2(40)=SumVec2(40)+ dpdx*dpdx
          SumVec2(41)=SumVec2(41)+ dpdx*dpdx*dpdx
          SumVec2(42)=SumVec2(42)+ dpdx*dpdx*dpdx*dpdx
          SumVec2(43)=SumVec2(43)+ dpdy                 ! yp !
          SumVec2(44)=SumVec2(44)+ dpdy*dpdy            ! yp !
          SumVec2(45)=SumVec2(45)+ dpdy*dpdy*dpdy       ! yp !
          SumVec2(46)=SumVec2(46)+ dpdy*dpdy*dpdy*dpdy  ! yp !
          SumVec2(47)=SumVec2(47)+ dpdz
          SumVec2(48)=SumVec2(48)+ dpdz*dpdz
          SumVec2(49)=SumVec2(49)+ dpdz*dpdz*dpdz
          SumVec2(50)=SumVec2(50)+ dpdz*dpdz*dpdz*dpdz
          
          SumVec2(51)=SumVec2(51)+ dudxx
          SumVec2(52)=SumVec2(52)+ dudxx*dudxx
          SumVec2(53)=SumVec2(53)+ dudxx*dudxx*dudxx
          SumVec2(54)=SumVec2(54)+ dudxx*dudxx*dudxx*dudxx
          SumVec2(55)=SumVec2(55)+ dudyy
          SumVec2(56)=SumVec2(56)+ dudyy*dudyy
          SumVec2(57)=SumVec2(57)+ dudyy*dudyy*dudyy
          SumVec2(58)=SumVec2(58)+ dudyy*dudyy*dudyy*dudyy
          SumVec2(59)=SumVec2(59)+ dudzz
          SumVec2(60)=SumVec2(60)+ dudzz*dudzz
          SumVec2(61)=SumVec2(61)+ dudzz*dudzz*dudzz
          SumVec2(62)=SumVec2(62)+ dudzz*dudzz*dudzz*dudzz
          
          SumVec2(63)=SumVec2(63)+ dvdxx                   ! yp !
          SumVec2(64)=SumVec2(64)+ dvdxx*dvdxx             ! yp !
          SumVec2(65)=SumVec2(65)+ dvdxx*dvdxx*dvdxx       ! yp !
          SumVec2(66)=SumVec2(66)+ dvdxx*dvdxx*dvdxx*dvdxx ! yp !
          SumVec2(67)=SumVec2(67)+ dvdyy                   ! yp !
          SumVec2(68)=SumVec2(68)+ dvdyy*dvdyy             ! yp !
          SumVec2(69)=SumVec2(69)+ dvdyy*dvdyy*dvdyy       ! yp !
          SumVec2(70)=SumVec2(70)+ dvdyy*dvdyy*dvdyy*dvdyy ! yp !
          SumVec2(71)=SumVec2(71)+ dvdzz                   ! yp !
          SumVec2(72)=SumVec2(72)+ dvdzz*dvdzz             ! yp !
          SumVec2(73)=SumVec2(73)+ dvdzz*dvdzz*dvdzz       ! yp !
          SumVec2(74)=SumVec2(74)+ dvdzz*dvdzz*dvdzz*dvdzz ! yp !                                   

          SumVec2(75)=SumVec2(75)+ dwdxx
          SumVec2(76)=SumVec2(76)+ dwdxx*dwdxx
          SumVec2(77)=SumVec2(77)+ dwdxx*dwdxx*dwdxx
          SumVec2(78)=SumVec2(78)+ dwdxx*dwdxx*dwdxx*dwdxx
          SumVec2(79)=SumVec2(79)+ dwdyy
          SumVec2(80)=SumVec2(80)+ dwdyy*dwdyy
          SumVec2(81)=SumVec2(81)+ dwdyy*dwdyy*dwdyy
          SumVec2(82)=SumVec2(82)+ dwdyy*dwdyy*dwdyy*dwdyy
          SumVec2(83)=SumVec2(83)+ dwdzz
          SumVec2(84)=SumVec2(84)+ dwdzz*dwdzz
          SumVec2(85)=SumVec2(85)+ dwdzz*dwdzz*dwdzz
          SumVec2(86)=SumVec2(86)+ dwdzz*dwdzz*dwdzz*dwdzz
          
          SumVec2(87)=SumVec2(87)+ dpdxx
          SumVec2(88)=SumVec2(88)+ dpdxx*dpdxx
          SumVec2(89)=SumVec2(89)+ dpdxx*dpdxx*dpdxx
          SumVec2(90)=SumVec2(90)+ dpdxx*dpdxx*dpdxx*dpdxx
          SumVec2(91)=SumVec2(91)+ dpdyy
          SumVec2(92)=SumVec2(92)+ dpdyy*dpdyy
          SumVec2(93)=SumVec2(93)+ dpdyy*dpdyy*dpdyy
          SumVec2(94)=SumVec2(94)+ dpdyy*dpdyy*dpdyy*dpdyy
          SumVec2(95)=SumVec2(95)+ dpdzz
          SumVec2(96)=SumVec2(96)+ dpdzz*dpdzz
          SumVec2(97)=SumVec2(97)+ dpdzz*dpdzz*dpdzz
          SumVec2(98)=SumVec2(98)+ dpdzz*dpdzz*dpdzz*dpdzz
#endif
        enddo
      enddo
      do kc=1,NCHASTAT
        SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
      enddo
#ifdef HighOrderGradStat
      do kc=1,nGradStat
        SumGrad(kc,jc)=SumGrad(kc,jc)+ SumVec2(kc)*inxz
      enddo
#endif
    ENDDO

    ! nyp only
    jc=nyp; jm=jc-1; cac=rdyc(jc); SumVec=zero
#ifdef HighOrderGradStat
    SumVec2=zero
#endif
    InterpY1= YinterpCoe(jc); InterpY2=one-InterpY1
    do kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do ic=y1start(1),y1end(1)
        is=ic-2;im=ic-1;ip=ic+1
        prloc2= InterpY1*(InterpCoe1*pr(is,jm,kc) +InterpCoe2*pr(im,jm,kc) +InterpCoe3*pr(ic,jm,kc) +InterpCoe4*pr(ip,jm,kc))+ &
                InterpY2*(InterpCoe1*pr(is,jc,kc) +InterpCoe2*pr(im,jc,kc) +InterpCoe3*pr(ic,jc,kc) +InterpCoe4*pr(ip,jc,kc)) !2021-12-14
        dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
        dwdy= (uz(ic,jc,kc)-uz(ic,jm,kc))*cac
        vor_x=  dwdy
        vor_z= -dudy
        SumVec(27)=SumVec(27)+ prloc2*(dudy+zero)     ! yp !
        SumVec(33)=SumVec(33)+ vor_x*vor_x            ! yp !
        SumVec(35)=SumVec(35)+ vor_z*vor_z            ! yp !
#if defined(CFDLPT_TwoWay)
        SumVec(37)=SumVec(37)+ FpForce_y(ic,jc,kc)                     ! yp !
        SumVec(40)=SumVec(40)+ FpForce_y(ic,jc,kc)*FpForce_y(ic,jc,kc) ! yp !
        SumVec(43)=SumVec(43)+ FpForce_y(ic,jc,kc)*uy(ic,jc,kc)        ! yp !
#endif
#ifdef HighOrderGradStat
        SumVec2( 7)=SumVec2( 7)+ dudy                 ! yp !
        SumVec2( 8)=SumVec2( 8)+ dudy*dudy            ! yp !
        SumVec2( 9)=SumVec2( 9)+ dudy*dudy*dudy       ! yp !
        SumVec2(10)=SumVec2(10)+ dudy*dudy*dudy*dudy  ! yp !
        
        SumVec2(31)=SumVec2(31)+ dwdy                 ! yp !
        SumVec2(32)=SumVec2(32)+ dwdy*dwdy            ! yp !
        SumVec2(33)=SumVec2(33)+ dwdy*dwdy*dwdy       ! yp !
        SumVec2(34)=SumVec2(34)+ dwdy*dwdy*dwdy*dwdy  ! yp !
#endif
      enddo
    enddo
    do kc=1,NCHASTAT
      SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
    enddo
#ifdef HighOrderGradStat
    do kc=1,nGradStat
      SumGrad(kc,jc)=SumGrad(kc,jc)+ SumVec2(kc)*inxz
    enddo
#endif

    ! shear stress and pr gradient
    PrGradsum   = PrGradsum+ PrGradAver
    if(nrank==0 .and. IsUxConst) then
      write(filename,'(A,I10.10)')trim(ResultsDir)//"PrGrad",ilast
      open(newunit=nUnit,file=filename,status='old',position='append',form='formatted',IOSTAT=ierror)
      if(ierror/=0) then
        call MainLog%CheckForError(ErrT_Pass,"clcStat","Cannot open file: "//trim(filename))
      else
        write(nUnit,'(I7,2ES24.15)')itime,SimTime,PrGradAver
      endif
      close(nUnit,IOSTAT=ierror)
    endif
    nfstime= nfstime + 1

    IF(mod(itime,ivSpec)==0) THEN
      ! Determine ux in cell center
      DO kc=y1start(3),y1end(3)
        kt=kc-y1start(3)+1
        do jc=y1start(2),y1end(2)
          jt=jc
          do ic=y1start(1),y1end(1)
            it=ic-y1start(1)+1
            im=ic-1;ip=ic+1;iu=ic+2
            ArrTemp1(it,jt,kt)=InterpCoe1*ux(im,jc,kc) +InterpCoe2*ux(ic,jc,kc) + &
                               InterpCoe3*ux(ip,jc,kc) +InterpCoe4*ux(iu,jc,kc) + uCRF
          enddo
        enddo
      ENDDO
      
      !=========================== Calculate 1-D Energy Spectra ===========================!
      IF(clcSpectra1D) THEN
        !=============== Spectra in x-dir ===============
        allocate(arrx1(x1size(1),x1size(2),x1size(3)))
        allocate(arrx2(x1size(1),x1size(2),x1size(3)))
#ifdef OverWriteFFT
        call transpose_y1_to_x1(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1)
        arrx1=arrx1+uCRF
        call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DUU)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DUU,iLCSI1DUU,iLCSR1DUU2) ! LCS -Real -Imag -ux
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DUU,iLCSR2DUU)
        
        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc; jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=half*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO
        call transpose_y1_to_x1(ArrTemp2,arrx2)
        call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
        call clcCospectraX(arrx1,arrx2,iSpec1DUV2)
        call transpose_y1_to_x1(ArrTemp1,arrx1)
        call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
        call clcCospectraX(arrx1,arrx2,iSpec1DUV)    
#ifdef ScalarFlow
        call transpose_y1_to_x1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1)
        call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
        call clcCospectraX(arrx1,arrx2,iSpec1DVC)
#endif

        call transform_uy(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),ArrTemp2)
        call transpose_y1_to_x1(ArrTemp2,arrx1)
        call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DVV)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DVV,iLCSI1DVV,iLCSR1DVV2) ! LCS -Real -Imag -uy
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DVV,iLCSR2DVV)
        
        call transpose_y1_to_x1(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1)
        call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DWW)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DWW,iLCSI1DWW,iLCSR1DWW2) ! LCS -Real -Imag -uz
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DWW,iLCSR2DWW)
        
        call transpose_y1_to_x1(pr(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1)
        call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DPP)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DPP,iLCSI1DPP,iLCSR1DPP2) ! LCS -Real -Imag -pp
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DPP,iLCSR2DPP)
#ifdef ScalarFlow
        call transpose_y1_to_x1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrx1)
        call dfftw_execute_r2r(fft_plan_x,arrx1,arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DCC)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DCC,iLCSI1DCC,iLCSR1DCC2) ! LCS -Real -Imag -cc
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DCC,iLCSR2DCC)
        
        call transpose_y1_to_x1(ArrTemp1,arrx2)
        call dfftw_execute_r2r(fft_plan_x,arrx2,arrx2)
        call clcCospectraX(arrx1,arrx2,iSpec1DUC)
#endif
        deallocate(arrx1,arrx2)
#else
        allocate(arrFFT(x1size(1),x1size(2),x1size(3)))
        call transpose_y1_to_x1(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        arrFFT=arrFFT+uCRF
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DUU)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DUU,iLCSI1DUU,iLCSR1DUU2) ! LCS -Real -Imag -ux
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DUU,iLCSR2DUU)
        
        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc; jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=half*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO
        call transpose_y1_to_x1(ArrTemp2,arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx2)
        call clcCospectraX(arrx1,arrx2,iSpec1DUV2)   
        call transpose_y1_to_x1(ArrTemp1,arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx1)
        call clcCospectraX(arrx1,arrx2,iSpec1DUV)
#ifdef ScalarFlow
        call transpose_y1_to_x1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx1)
        call clcCospectraX(arrx1,arrx2,iSpec1DVC)
#endif

        call transform_uy(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),ArrTemp2)
        call transpose_y1_to_x1(ArrTemp2,arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DVV)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DVV,iLCSI1DVV,iLCSR1DVV2) ! LCS -Real -Imag -uy
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DVV,iLCSR2DVV)
        
        call transpose_y1_to_x1(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DWW)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DWW,iLCSI1DWW,iLCSR1DWW2) ! LCS -Real -Imag -uz
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DWW,iLCSR2DWW)
        
        call transpose_y1_to_x1(pr(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DPP)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DPP,iLCSI1DPP,iLCSR1DPP2) ! LCS -Real -Imag -pp
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DPP,iLCSR2DPP)
#ifdef ScalarFlow
        call transpose_y1_to_x1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx1)
        call clcEnergySpectraX(arrx1,iSpec1DCC)
        call clcLTS_x(arrx1,arrx2,ArrTemp2,iLCSR1DCC,iLCSI1DCC,iLCSR1DCC2) ! LCS -Real -Imag -cc
        if(clcSpectra2D) call clcSpecAndLCSR2DFrom1D(arrx1,iSpec2DCC,iLCSR2DCC)
        
        call transpose_y1_to_x1(ArrTemp1,arrFFT)
        call dfftw_execute_r2r(fft_plan_x,arrFFT,arrx2)
        call clcCospectraX(arrx1,arrx2,iSpec1DUC)
#endif       
        deallocate(arrx1,arrx2,arrFFT)
#endif
        
        !=============== Spectra in z-dir ===============
        allocate(arrz1(z1size(1),z1size(2),z1size(3)))
        allocate(arrz2(z1size(1),z1size(2),z1size(3)))
#ifdef OverWriteFFT
        call transpose_y1_to_z1(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1)
        arrz1=arrz1+uCRF
        call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DUU)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DUU,iLCSI1DUU,iLCSR1DUU2) ! LCS -Real -Imag -ux

        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc; jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=half*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO
        call transpose_y1_to_z1(ArrTemp2,arrz2)
        call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
        call clcCospectraZ(arrz1,arrz2,iSpec1DUV2)
        call transpose_y1_to_z1(ArrTemp1,arrz1)
        call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
        call clcCospectraZ(arrz1,arrz2,iSpec1DUV)  
#ifdef ScalarFlow
        call transpose_y1_to_z1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1)
        call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
        call clcCospectraZ(arrz1,arrz2,iSpec1DVC)
#endif

        call transform_uy(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),ArrTemp2)
        call transpose_y1_to_z1(ArrTemp2,arrz1)
        call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DVV)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DVV,iLCSI1DVV,iLCSR1DVV2) ! LCS -Real -Imag -uy
        
        call transpose_y1_to_z1(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1)
        call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DWW)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DWW,iLCSI1DWW,iLCSR1DWW2) ! LCS -Real -Imag -uz

        call transpose_y1_to_z1(pr(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1)
        call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DPP)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DPP,iLCSI1DPP,iLCSR1DPP2) ! LCS -Real -Imag -pp
#ifdef ScalarFlow
        call transpose_y1_to_z1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrz1)
        call dfftw_execute_r2r(fft_plan_z,arrz1,arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DCC)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DCC,iLCSI1DCC,iLCSR1DCC2) ! LCS -Real -Imag -cc
        
        call transpose_y1_to_z1(ArrTemp1,arrz2)
        call dfftw_execute_r2r(fft_plan_z,arrz2,arrz2)
        call clcCospectraZ(arrz1,arrz2,iSpec1DUC)
#endif
        deallocate(arrz1,arrz2)
#else
        allocate(arrFFT(z1size(1),z1size(2),z1size(3)))
        call transpose_y1_to_z1(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        arrFFT=arrFFT+uCRF
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DUU)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DUU,iLCSI1DUU,iLCSR1DUU2) ! LCS -Real -Imag -ux

        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc; jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=half*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO
        call transpose_y1_to_z1(ArrTemp2,arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz2)
        call clcCospectraZ(arrz1,arrz2,iSpec1DUV2)
        call transpose_y1_to_z1(ArrTemp1,arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz1)
        call clcCospectraZ(arrz1,arrz2,iSpec1DUV) 
         
#ifdef ScalarFlow
        call transpose_y1_to_z1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz1)
        call clcCospectraZ(arrz1,arrz2,iSpec1DVC)
#endif

        call transform_uy(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),ArrTemp2)
        call transpose_y1_to_z1(ArrTemp2,arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DVV)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DVV,iLCSI1DVV,iLCSR1DVV2) ! LCS -Real -Imag -uy
        
        call transpose_y1_to_z1(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DWW)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DWW,iLCSI1DWW,iLCSR1DWW2) ! LCS -Real -Imag -uz

        call transpose_y1_to_z1(pr(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DPP)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DPP,iLCSI1DPP,iLCSR1DPP2) ! LCS -Real -Imag -pp
#ifdef ScalarFlow
        call transpose_y1_to_z1(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz1)
        call clcEnergySpectraZ(arrz1,iSpec1DCC)
        call clcLTS_z(arrz1,arrz2,ArrTemp2,iLCSR1DCC,iLCSI1DCC,iLCSR1DCC2) ! LCS -Real -Imag -cc
        
        call transpose_y1_to_z1(ArrTemp1,arrFFT)
        call dfftw_execute_r2r(fft_plan_z,arrFFT,arrz2)
        call clcCospectraZ(arrz1,arrz2,iSpec1DUC)
#endif  
        deallocate(arrz1,arrz2,arrFFT)
#endif
      ENDIF
      !=========================== Calculate 1-D Energy Spectra ===========================!

      !=========================== Calculate 2-D Energy Spectra ===========================!
      IF(clcSpectra2D) THEN
        ! Determine uy in cell center
        DO kc=y1start(3),y1end(3)
          kt=kc-y1start(3)+1
          do jc=y1start(2),y1end(2)
            jt=jc
            jp=jc+1
            do ic=y1start(1),y1end(1)
              it=ic-y1start(1)+1
              ArrTemp2(it,jt,kt)=half*(uy(ic,jc,kc)+uy(ic,jp,kc))
            enddo
          enddo
        ENDDO
        call clcCosSpectra2D(ArrTemp1,ArrTemp2,iSpec2DUV,.false.)
        if(.not. clcSpectra1D) then
          call transform_uy(uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),ArrTemp2)
          call clcSpectraAndLCSR2D(ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))+uCRF,  &
                                   iSpec2DUU,iLCSR2DUU)
          call clcSpectraAndLCSR2D(ArrTemp2,                                                              &
                                   iSpec2DVV,iLCSR2DVV)
          call clcSpectraAndLCSR2D(uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),       &
                                   iSpec2DWW,iLCSR2DWW)
          call clcSpectraAndLCSR2D(pr(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)), &
                                   iSpec2DPP,iLCSR2DPP)
#ifdef ScalarFlow
          call clcSpectraAndLCSR2D(scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),   &
                                   iSpec2DCC,iLCSR2DCC)
#endif
        endif
      ENDIF
      !=========================== Calculate 2-D Energy Spectra ===========================!
      nSpectime= nSpectime+1
    ENDIF
    if(mod(itime,SaveStat)/=0) return
    
    ! Write statistics
    call MPI_REDUCE(SumStat,SumStatR,NCHASTAT*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
#ifdef HighOrderGradStat
    call MPI_REDUCE(SumGrad,SumGradR,nGradStat*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
#endif
    if(nrank==0) then
      infstime = one/real(nfstime,RK)
      write(filename,"(A,I10.10)") trim(ResultsDir)//'stats',itime
      open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
      IF(ierror/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat","Cannot open file: "//trim(filename))
      ELSE
        write(nUnit,'(a,I7,a,I7,a,I7)')'  The time step range for this fluid statistics is ', &
                                    itime-(nfstime-1)*ivstats, ':', ivstats, ':', itime
        write(nUnit,'(A)')'  '
        if(IsUxConst) then
          write(nUnit,'(A)')'  Constant velocity in x-dir by adding a pressure gradient.'
          write(nUnit,'(A, ES24.15)')'    time averaged pressure gradient is: ',abs(PrGradsum)*infstime
        else
          dudy = abs(SumStatR(1,1))*infstime*2.0_RK*rdyc(1)
          if(FlowType==FT_CH) then
            dudyU= abs(SumStatR(1,nyc))*infstime*2.0_RK*rdyc(nyp)
            dudy = dudy+dudyU
          endif
          write(nUnit,'(A, ES24.15)')'  Variable velocity in x-dir while adding a constant body force.',xnu*dudy/yly
        endif
        write(nUnit,'(A)')'  '
        
        Block 
        character(len=128)::FormatStr
        write(FormatStr,'(A,I3,A)')'(',NCHASTAT,'ES24.15)'
        do jc=1,nyp
          write(nUnit,FormatStr)SumStatR(1:NCHASTAT,jc)*infstime
        enddo
        End Block
      ENDIF
      close(nUnit,IOSTAT=ierror)
#ifdef HighOrderGradStat
      write(filename,"(A,I10.10)") trim(ResultsDir)//'vGrad',itime
      open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
      IF(ierror/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat","Cannot open file: "//trim(filename))
      ELSE
        Block 
        character(len=128)::FormatStr
        write(FormatStr,'(A,I3,A)')'(',nGradStat,'ES24.15)'
        do jc=1,nyp
          write(nUnit,FormatStr)SumGradR(1:nGradStat,jc)*infstime
        enddo
        End Block
      ENDIF
      close(nUnit,IOSTAT=ierror)
#endif
    endif

    ! Write 1-D Energy Spectra
    IF(clcSpectra1D) THEN
      infstime= one/real(nSpectime,RK)

      ! Spectra in x-dir
      allocate(EnergySpecXR(nxhp,x1size(2),NEnergySpec1D));EnergySpecXR=zero
      call MPI_REDUCE(EnergySpecX,EnergySpecXR,NEnergySpec1D*x1size(2)*nxhp,real_type,MPI_SUM,0,DECOMP_2D_COMM_ROW,ierror)
      call MPI_COMM_RANK(DECOMP_2D_COMM_ROW,nrankX,ierror)
      if(nrankX==0) then
        EnergySpecXR=EnergySpecXR*infstime
        disp_inc=int(mytype_bytes,8)*int(nyc,8)*int(nxhp,8)
        disp=int(mytype_bytes,8)*int(x1start(2)-1,8)*int(nxhp,8)
        
        write(filename,"(A,A,I10.10)") trim(ResultsDir),'SpecX',itime
        call MPI_FILE_OPEN(DECOMP_2D_COMM_COL,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,nUnit,ierror)
        call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
        do kc=1,NEnergySpec1D
          call MPI_FILE_WRITE_AT_ALL(nUnit,disp,EnergySpecXR(:,:,kc),x1size(2)*nxhp,real_type,MPI_STATUS_IGNORE,ierror)
          disp=disp+disp_inc
        enddo
        call MPI_FILE_CLOSE(nUnit,ierror)
      endif
      deallocate(EnergySpecXR)
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)

      ! Spectra in z-dir
      allocate(EnergySpecZR(nzhp,z1size(2),NEnergySpec1D));EnergySpecZR=zero
      call MPI_REDUCE(EnergySpecZ,EnergySpecZR,NEnergySpec1D*z1size(2)*nzhp,real_type,MPI_SUM,0,DECOMP_2D_COMM_COL,ierror)
      call MPI_COMM_RANK(DECOMP_2D_COMM_COL,nrankZ,ierror)
      if(nrankZ==0) then
        EnergySpecZR=EnergySpecZR*infstime
        disp_inc=int(mytype_bytes,8)*int(nyc,8)*int(nzhp,8)
        disp=int(mytype_bytes,8)*int(z1start(2)-1,8)*int(nzhp,8)
        
        write(filename,"(A,A,I10.10)") trim(ResultsDir),'SpecZ',itime
        call MPI_FILE_OPEN(DECOMP_2D_COMM_ROW,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,nUnit,ierror)
        call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
        do kc=1,NEnergySpec1D
          call MPI_FILE_WRITE_AT_ALL(nUnit,disp,EnergySpecZR(:,:,kc),z1size(2)*nzhp,real_type,MPI_STATUS_IGNORE,ierror)
          disp=disp+disp_inc
        enddo
        call MPI_FILE_CLOSE(nUnit,ierror)
      endif
      deallocate(EnergySpecZR)
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    ENDIF

    ! Write 2-D Energy Spectra
    IF(clcSpectra2D) THEN
      infstime= one/real(nSpectime,RK)
      Block
      real(RK),dimension(:,:),allocatable::Ratio2D
      allocate(Ratio2D(decomp_xhzh%y2sz(1),decomp_xhzh%y2sz(3)))
      Ratio2D=4.0_RK
      if(decomp_xhzh%y2st(1)==1)    Ratio2D(1,:)=2.0_RK
      if(decomp_xhzh%y2st(3)==1)    Ratio2D(:,1)=2.0_RK
      if(decomp_xhzh%y2en(1)==nxhp) Ratio2D(decomp_xhzh%y2sz(1),:)=2.0_RK
      if(decomp_xhzh%y2en(3)==nzhp) Ratio2D(:,decomp_xhzh%y2sz(3))=2.0_RK
      if(decomp_xhzh%y2st(1)==1    .and. decomp_xhzh%y2st(3)==1   )  then
        Ratio2D(1,1)=1.0_RK
      endif
      if(decomp_xhzh%y2st(1)==1    .and. decomp_xhzh%y2en(3)==nzhp)  then
        Ratio2D(1,decomp_xhzh%y2sz(3))=1.0_RK
      endif
      if(decomp_xhzh%y2en(1)==nxhp .and. decomp_xhzh%y2en(3)==nzhp)  then
        Ratio2D(decomp_xhzh%y2sz(1),decomp_xhzh%y2sz(3))=1.0_RK
      endif
      if(decomp_xhzh%y2en(1)==nxhp .and. decomp_xhzh%y2st(3)==1   )  then
        Ratio2D(decomp_xhzh%y2sz(1),1)=1.0_RK
      endif
      Ratio2D=infstime*Ratio2D
      do kc=1,decomp_xhzh%y2sz(3)
        do ic=1,decomp_xhzh%y2sz(1)
          EnergySpec2D(ic,:,kc,:)=Ratio2D(ic,kc)*EnergySpec2D(ic,:,kc,:)
        enddo
      enddo   
      deallocate(Ratio2D)
      end Block
 
      Block
      integer(MPI_OFFSET_KIND)::disp
      integer::data_type,data_byte,newtype,nySpec2D
      integer,dimension(3)::sizes,subsizes,starts
#ifdef SAVE_SINGLE_Spec2D
      real(4),allocatable,dimension(:,:,:,:)::EnergySpecOut
#endif      
      if(FlowType==FT_CH) then
        nySpec2D=nyc/2
      else
        nySpec2D=nyc
      endif
#ifdef SAVE_SINGLE_Spec2D
      data_type= MPI_REAL
      allocate(EnergySpecOut(decomp_xhzh%y2sz(1),nySpec2D,decomp_xhzh%y2sz(3),NEnergySpec2D))
      do kt=1,NEnergySpec2D
        do kc=1,decomp_xhzh%y2sz(3)
          do jc=1,nySpec2D
            do ic=1,decomp_xhzh%y2sz(1)
              EnergySpecOut(ic,jc,kc,kt)=real(EnergySpec2D(ic,jc,kc,kt),4)
            enddo
          enddo
        enddo
      enddo
#else
      data_type= MPI_DOUBLE_PRECISION
#endif
      call MPI_TYPE_SIZE(data_type,data_byte,ierror)
      sizes(1)= nxhp
      sizes(2)= nySpec2D
      sizes(3)= nzhp 
      subsizes(1)=decomp_xhzh%y2sz(1)
      subsizes(2)=nySpec2D
      subsizes(3)=decomp_xhzh%y2sz(3)
      starts(1)=decomp_xhzh%y2st(1)-1
      starts(2)=0
      starts(3)=decomp_xhzh%y2st(3)-1
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      
      disp = 0_MPI_OFFSET_KIND
      write(filename,"(A,A,I10.10)") trim(ResultsDir),'Spec2D',itime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      do jc=1,NEnergySpec2D
        call MPI_FILE_SET_VIEW(nUnit,disp,data_type,newtype,'native',MPI_INFO_NULL,ierror)

#ifdef SAVE_SINGLE_Spec2D
        call MPI_FILE_WRITE_ALL(nUnit,EnergySpecOut(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
#else
        call MPI_FILE_WRITE_ALL(nUnit,EnergySpec2D(:,:,:,jc), subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
#endif
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(data_byte,8) ! Update displacement
      enddo
      call MPI_FILE_CLOSE(nUnit,ierror)
      call MPI_TYPE_FREE(newtype,ierror)
      End Block
    ENDIF
    
    nfstime=0; SumStat=zero; PrGradsum=zero; nSpectime=0; 
    if(clcSpectra1D) then
      EnergySpecX=zero; EnergySpecZ=zero
    endif
    if(clcSpectra2D) EnergySpec2D=zero
#ifdef HighOrderGradStat
    SumGrad=zero
#endif
  end subroutine clcStat

  !******************************************************************
  ! clcEnergySpectraX
  !******************************************************************
  subroutine clcEnergySpectraX(arrx,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,x1size(2)
      EgyX=zero
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx(1,jc,kc)   *arrx(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx(nxhp,jc,kc)*arrx(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx(ic,jc,kc)*arrx(ic,jc,kc)+arrx(icCnter,jc,kc)*arrx(icCnter,jc,kc))*two
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcEnergySpectraX
   
  !******************************************************************
  ! clcCospectraX
  !******************************************************************
  subroutine clcCospectraX(arrx1,arrx2,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx1,arrx2

    ! locals
    integer::ic,jc,kc,icCnter
    real(RK)::normEgy,EgyX(nxhp)

    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
    do jc=1,x1size(2)
      EgyX=zero
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx1(1,jc,kc)   *arrx2(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx1(nxhp,jc,kc)*arrx2(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx1(ic,jc,kc)*arrx2(ic,jc,kc)+arrx1(icCnter,jc,kc)*arrx2(icCnter,jc,kc))*two
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,m)= EnergySpecX(ic,jc,m) +EgyX(ic)*normEgy
      enddo
    enddo
  end subroutine clcCospectraX

  !******************************************************************
  ! clcEnergySpectraZ
  !******************************************************************
  subroutine clcEnergySpectraZ(arrz,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in)::arrz

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,z1size(2)
      EgyZ=zero
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz(ic,jc,1)   *arrz(ic,jc,1)
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz(ic,jc,nzhp)*arrz(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz(ic,jc,kc)*arrz(ic,jc,kc)+arrz(ic,jc,kcCnter)*arrz(ic,jc,kcCnter))*two
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcEnergySpectraZ
    
  !******************************************************************
  ! clcCospectraZ
  !******************************************************************
  subroutine clcCospectraZ(arrz1,arrz2,m)
    implicit none
    integer,intent(in)::m
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in)::arrz1,arrz2

    ! locals
    integer::ic,jc,kc,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)

    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
    do jc=1,z1size(2)
      EgyZ=zero
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz1(ic,jc,1)   *arrz2(ic,jc,1)
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz1(ic,jc,nzhp)*arrz2(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz1(ic,jc,kc)*arrz2(ic,jc,kc)+arrz1(ic,jc,kcCnter)*arrz2(ic,jc,kcCnter))*two
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,m)= EnergySpecZ(kc,jc,m) +EgyZ(kc)*normEgy
      enddo
    enddo
  end subroutine clcCospectraZ

  !******************************************************************
  ! clcLTS_x
  !******************************************************************  
  subroutine clcLTS_x(arrx1,arrx2,ArrTemp,n1,n2,n3)
    implicit none
    integer,intent(in)::n1,n2,n3
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in) ::arrx1
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(out)::arrx2
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::ArrTemp

    ! locals
    integer::ic,jc,kc,jt,icCnter
    real(RK)::normEgy,EgyX(nxhp)
    real(RK),allocatable,dimension(:,:,:)::arrYplane
    
    if(FlowType==FT_HC) then
      allocate(arrYplane(y1size(1),y1size(3),1))
    else
      allocate(arrYplane(y1size(1),y1size(3),2))
    endif
    normEgy= one/(real(nxc,RK)*real(nxc,RK)*real(nzc,RK))
 
    ! LCS -Real -Ref1
    call transpose_x1_to_y1(arrx1,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(1),kc)
        enddo
      enddo
    ELSE
      jt=nyc+1-jForLCS(1)
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=ArrTemp(ic,jt,kc)
        enddo
      enddo
    ENDIF
    
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_x1(ArrTemp,arrx2)
    do jc=1,x1size(2)
      EgyX=zero
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx2(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx2(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx2(ic,jc,kc)+arrx2(icCnter,jc,kc)) ! Note, there is NO "*two" here !
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,n1)= EnergySpecX(ic,jc,n1) +EgyX(ic)*normEgy
      enddo
    enddo
    
    ! LCS -Imag -Ref1
    arrx2(1,:,:)=arrx1(1,:,:)
    do kc=1,x1size(3)
      do jc=1,x1size(2)
        do ic=2,x1size(1)
          icCnter=nxc+2-ic
          arrx2(ic,jc,kc)=arrx1(icCnter,jc,kc)
        enddo
      enddo
    enddo
    call transpose_x1_to_y1(arrx2,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_x1(ArrTemp,arrx2)
    do jc=1,x1size(2)
      EgyX=zero
      ! EgyX(1)   = EgyX(nxhp)= zero
      do kc=1,x1size(3)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx2(ic,jc,kc)-arrx2(icCnter,jc,kc)) ! Note here !
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,n2)= EnergySpecX(ic,jc,n2) +EgyX(ic)*normEgy
      enddo
    enddo   

    ! LCS -Real -Ref2
    call transpose_x1_to_y1(arrx1,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(2),kc)
        enddo
      enddo
    ELSE
      jt=nyc+1-jForLCS(2)
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(2),kc)
          arrYplane(ic,kc,2)=ArrTemp(ic,jt,kc)
        enddo
      enddo
    ENDIF
    
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_x1(ArrTemp,arrx2)
    do jc=1,x1size(2)
      EgyX=zero
      do kc=1,x1size(3)
        EgyX(1)   = EgyX(1)   + arrx2(1,jc,kc)
        EgyX(nxhp)= EgyX(nxhp)+ arrx2(nxhp,jc,kc)
        do ic=2,nxh
          icCnter=nxc+2-ic
          EgyX(ic)= EgyX(ic) + (arrx2(ic,jc,kc)+arrx2(icCnter,jc,kc)) ! Note, there is NO "*two" here !
        enddo
      enddo
      do ic=1,nxhp
        EnergySpecX(ic,jc,n3)= EnergySpecX(ic,jc,n3) +EgyX(ic)*normEgy
      enddo
    enddo

    deallocate(arrYplane)
  end subroutine clcLTS_x

  !******************************************************************
  ! clcLTS_z
  !******************************************************************
  subroutine clcLTS_z(arrz1,arrz2,ArrTemp,n1,n2,n3)
    implicit none
    integer,intent(in)::n1,n2,n3
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(in) ::arrz1
    real(RK),dimension(z1size(1),z1size(2),z1size(3)),intent(out)::arrz2
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::ArrTemp

    ! locals
    integer::ic,jc,kc,jt,kcCnter
    real(RK)::normEgy,EgyZ(nzhp)
    real(RK),allocatable,dimension(:,:,:)::arrYplane
    
    if(FlowType==FT_HC) then
      allocate(arrYplane(y1size(1),y1size(3),1))
    else
      allocate(arrYplane(y1size(1),y1size(3),2))        
    endif
    normEgy= one/(real(nxc,RK)*real(nzc,RK)*real(nzc,RK))
        
    ! LCS -Real -Ref1
    call transpose_z1_to_y1(arrz1,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(1),kc)
        enddo
      enddo
    ELSE
      jt=nyc+1-jForLCS(1)
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=ArrTemp(ic,jt,kc)
        enddo
      enddo        
    ENDIF
    
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_z1(ArrTemp,arrz2)
    do jc=1,z1size(2)
      EgyZ=zero
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz2(ic,jc,1)   
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz2(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz2(ic,jc,kc)+arrz2(ic,jc,kcCnter)) ! Note, there is NO "*two" here !
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,n1)= EnergySpecZ(kc,jc,n1) +EgyZ(kc)*normEgy
      enddo
    enddo
    
    ! LCS -Imag -Ref1
    arrz2(:,:,1)=arrz1(:,:,1)
    do kc=2,z1size(3)
      kcCnter=nzc+2-kc
      do jc=1,z1size(2)
        do ic=1,z1size(1)
          arrz2(ic,jc,kc)=arrz1(ic,jc,kcCnter)
        enddo
      enddo
    enddo
    call transpose_z1_to_y1(arrz2,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_z1(ArrTemp,arrz2)
    do jc=1,z1size(2)
      EgyZ=zero
      ! EgyZ(1)   = EgyZ(nzhp)= zero
      do ic=1,z1size(1)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz2(ic,jc,kc)-arrz2(ic,jc,kcCnter)) ! Note here !
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,n2)= EnergySpecZ(kc,jc,n2) +EgyZ(kc)*normEgy
      enddo
    enddo

    ! LCS -Real -Ref2
    call transpose_z1_to_y1(arrz1,ArrTemp)
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(2),kc)
        enddo
      enddo
    ELSE
      jt=nyc+1-jForLCS(2)
      do kc=1,y1size(3)
        do ic=1,y1size(1)
          arrYplane(ic,kc,1)=ArrTemp(ic,jForLCS(2),kc)
          arrYplane(ic,kc,2)=ArrTemp(ic,jt,kc)
        enddo
      enddo        
    ENDIF
    
    IF(FlowType==FT_HC) THEN
      do kc=1,y1size(3)
        do jc=1,y1size(2)
           do ic=1,y1size(1)
             ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
           enddo
        enddo
      enddo            
    ELSE
      do kc=1,y1size(3)
        do jc=1,nyc/2
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,y1size(1)
            ArrTemp(ic,jc,kc)=ArrTemp(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo
      enddo            
    ENDIF
    call transpose_y1_to_z1(ArrTemp,arrz2)
    do jc=1,z1size(2)
      EgyZ=zero
      do ic=1,z1size(1)
        EgyZ(1)   = EgyZ(1)   + arrz2(ic,jc,1)   
        EgyZ(nzhp)= EgyZ(nzhp)+ arrz2(ic,jc,nzhp)
        do kc=2,nzh
          kcCnter=nzc+2-kc
          EgyZ(kc)= EgyZ(kc) + (arrz2(ic,jc,kc)+arrz2(ic,jc,kcCnter)) ! Note, there is NO "*two" here !
        enddo
      enddo
      do kc=1,nzhp
        EnergySpecZ(kc,jc,n3)= EnergySpecZ(kc,jc,n3) +EgyZ(kc)*normEgy
      enddo
    enddo
    deallocate(arrYplane)
  end subroutine clcLTS_z

  !******************************************************************
  ! clcSpectraAndLCSR2D
  !******************************************************************
  subroutine clcSpectraAndLCSR2D(ArrIN,m1,m2)
    implicit none
    integer,intent(in)::m1,m2
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in)::ArrIn
    
    ! locals
    real(RK)::normEgy
    integer::ic,jc,kc,ict,jct,kct,jcs
    real(RK),allocatable,dimension(:,:,:)::arrYplane
    real(RK),dimension(:,:,:),allocatable::arrx,arrx_xhzf
    real(RK),dimension(:,:,:),allocatable::arry_xhzh,arry_xhzf
    real(RK),dimension(:,:,:),allocatable::arrz_xhzf,arrz_xhzh,arrz_xhzh_LCS
#ifndef OverWriteFFT
    real(RK),dimension(:,:,:),allocatable::ArrFFT
#endif
    
    normEgy= one/(real(nxc,RK)*real(nzc,RK))
    normEgy= normEgy*normEgy
    IF(FlowType==FT_CH)normEgy=normEgy*half
    allocate(arrx(nxc,x1size(2),x1size(3)))
#ifdef OverWriteFFT
    call transpose_y1_to_x1(ArrIN,arrx)
    call dfftw_execute_r2r(fft_plan_x, arrx, arrx)
#else
    allocate(ArrFFT(nxc,x1size(2),x1size(3)),arrx(nxc,x1size(2),x1size(3)))
    call transpose_y1_to_x1(ArrIN,ArrFFT)
    call dfftw_execute_r2r(fft_plan_x, ArrFFT, arrx)
    deallocate(ArrFFT)
#endif
    allocate(arrx_xhzf(nxhp,x1size(2),x1size(3)),   &
             arrz_xhzf(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    
    ! Real part
    arrx_xhzf=arrx(1:nxhp,:,:)
#ifdef OverWriteFFT
    call transpose_x1_to_z2(arrx_xhzf,arrz_xhzf,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrz_xhzf,arrz_xhzf)
#else
    allocate(arrFFT(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    call transpose_x1_to_z2(arrx_xhzf,arrFFT,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrFFT,arrz_xhzf)
    deallocate(arrFFT)
#endif
    arrx_xhzf(1,:,:)=zero; arrx_xhzf(nxhp,:,:)=zero
    do ic=2,nxh
      ict=nxc+2-ic
      arrx_xhzf(ic,:,:)=arrx(ict,:,:)
    enddo
    deallocate(arrx)
    allocate(arrz_xhzh(decomp_xhzh%z2sz(1),decomp_xhzh%z2sz(2),decomp_xhzh%z2sz(3)))
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc )*arrz_xhzf(ic,jc,kc )+ &
                               arrz_xhzf(ic,jc,kct)*arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc) 
      enddo
    enddo
    
    ! LCS
    allocate(arry_xhzf(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(2),decomp_xhzf%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzf,arry_xhzf,decomp_xhzf)
    IF(FlowType==FT_HC) THEN
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),1))
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,decomp_xhzf%y2sz(2)
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
      enddo
    else
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y1sz(3),2))
      jct=nyc+1-jForLCS(1)
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=arry_xhzf(ic,jct,kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,nyc/2
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo        
      enddo
    ENDIF
    call transpose_y2_to_z2(arry_xhzf,arrz_xhzf,decomp_xhzf)
    deallocate(arry_xhzf,arrYplane)
    allocate(arrz_xhzh_LCS(decomp_xhzh%z2sz(1),decomp_xhzh%z2sz(2),decomp_xhzh%z2sz(3)))
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc )+ arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc)
      enddo
    enddo
    
    ! Imag part
#ifdef OverWriteFFT
    call transpose_x1_to_z2(arrx_xhzf,arrz_xhzf,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrz_xhzf,arrz_xhzf)
#else
    allocate(arrFFT(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    call transpose_x1_to_z2(arrx_xhzf,arrFFT,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrFFT,arrz_xhzf)
    deallocate(arrFFT)
#endif
    deallocate(arrx_xhzf)
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc )*arrz_xhzf(ic,jc,kc )+ &
                               arrz_xhzf(ic,jc,kct)*arrz_xhzf(ic,jc,kct)+ arrz_xhzh(ic,jc,kc)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
      enddo
    enddo
    
    ! LCS
    allocate(arry_xhzf(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(2),decomp_xhzf%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzf,arry_xhzf,decomp_xhzf)
    IF(FlowType==FT_HC) THEN
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),1))
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,decomp_xhzf%y2sz(2)
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
      enddo
    else
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y1sz(3),2))
      jct=nyc+1-jForLCS(1)
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=arry_xhzf(ic,jct,kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,nyc/2
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo        
      enddo
    ENDIF
    call transpose_y2_to_z2(arry_xhzf,arrz_xhzf,decomp_xhzf)
    deallocate(arry_xhzf,arrYplane)
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzh_LCS(ic,jc,kc)+arrz_xhzf(ic,jc,kc )+ arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh_LCS(ic,jc,kc)=arrz_xhzh_LCS(ic,jc,kc)+ arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh_LCS(ic,jc,kc)=arrz_xhzh_LCS(ic,jc,kc)+ arrz_xhzf(ic,jc,kc)
      enddo
    enddo
    
    deallocate(arrz_xhzf)
    allocate(arry_xhzh(decomp_xhzh%y2sz(1),decomp_xhzh%y2sz(2),decomp_xhzh%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzh,arry_xhzh,decomp_xhzh)
    deallocate(arrz_xhzh)

    IF(FlowType==FT_HC) THEN
      do jc=1,nyc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m1)=EnergySpec2D(ic,jc,kc,m1)+normEgy*arry_xhzh(ic,jc,kc)
          enddo
        enddo
      enddo
    ELSE 
      do jc=1,nyc/2
        jcs=nyc+1-jc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m1)=EnergySpec2D(ic,jc,kc,m1)+ &
                                      normEgy*(arry_xhzh(ic,jc,kc)+arry_xhzh(ic,jcs,kc))
          enddo
        enddo
      enddo
    ENDIF
    
    ! LCS
    call transpose_z2_to_y2(arrz_xhzh_LCS,arry_xhzh,decomp_xhzh)
    deallocate(arrz_xhzh_LCS)
    IF(FlowType==FT_HC) THEN
      do jc=1,nyc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m2)=EnergySpec2D(ic,jc,kc,m2)+normEgy*arry_xhzh(ic,jc,kc)
          enddo
        enddo
      enddo
    ELSE 
      do jc=1,nyc/2
        jcs=nyc+1-jc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m2)=EnergySpec2D(ic,jc,kc,m2)+ &
                                      normEgy*(arry_xhzh(ic,jc,kc)+arry_xhzh(ic,jcs,kc))
          enddo
        enddo
      enddo
    ENDIF
    deallocate(arry_xhzh)
  end subroutine clcSpectraAndLCSR2D

  !******************************************************************
  ! clcSpecAndLCSR2DFrom1D
  !******************************************************************
  subroutine clcSpecAndLCSR2DFrom1D(arrx,m1,m2)
    implicit none
    integer,intent(in)::m1,m2
    real(RK),dimension(x1size(1),x1size(2),x1size(3)),intent(in)::arrx
    
    ! locals
    real(RK)::normEgy
    integer::ic,jc,kc,ict,jct,kct,jcs
    real(RK),allocatable,dimension(:,:,:)::arrYplane
    real(RK),dimension(:,:,:),allocatable::arrx_xhzf,arry_xhzh,arry_xhzf
    real(RK),dimension(:,:,:),allocatable::arrz_xhzf,arrz_xhzh,arrz_xhzh_LCS
#ifndef OverWriteFFT
    real(RK),dimension(:,:,:),allocatable::ArrFFT
#endif
    
    normEgy= one/(real(nxc,RK)*real(nzc,RK))
    normEgy= normEgy*normEgy
    IF(FlowType==FT_CH)normEgy=normEgy*half
    allocate(arrx_xhzf(nxhp,x1size(2),x1size(3)),   &
             arrz_xhzf(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    
    ! Real part
    arrx_xhzf=arrx(1:nxhp,:,:)
#ifdef OverWriteFFT
    call transpose_x1_to_z2(arrx_xhzf,arrz_xhzf,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrz_xhzf,arrz_xhzf)
#else
    allocate(arrFFT(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    call transpose_x1_to_z2(arrx_xhzf,arrFFT,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrFFT,arrz_xhzf)
    deallocate(arrFFT)
#endif
    arrx_xhzf(1,:,:)=zero; arrx_xhzf(nxhp,:,:)=zero
    do ic=2,nxh
      ict=nxc+2-ic
      arrx_xhzf(ic,:,:)=arrx(ict,:,:)
    enddo
    allocate(arrz_xhzh(decomp_xhzh%z2sz(1),decomp_xhzh%z2sz(2),decomp_xhzh%z2sz(3)))
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc )*arrz_xhzf(ic,jc,kc )+ &
                               arrz_xhzf(ic,jc,kct)*arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc) 
      enddo
    enddo
    
    ! LCS
    allocate(arry_xhzf(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(2),decomp_xhzf%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzf,arry_xhzf,decomp_xhzf)
    IF(FlowType==FT_HC) THEN
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),1))
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,decomp_xhzf%y2sz(2)
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
      enddo
    else
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y1sz(3),2))
      jct=nyc+1-jForLCS(1)
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=arry_xhzf(ic,jct,kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,nyc/2
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo        
      enddo
    ENDIF
    call transpose_y2_to_z2(arry_xhzf,arrz_xhzf,decomp_xhzf)
    deallocate(arry_xhzf,arrYplane)
    allocate(arrz_xhzh_LCS(decomp_xhzh%z2sz(1),decomp_xhzh%z2sz(2),decomp_xhzh%z2sz(3)))
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc )+ arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzf(ic,jc,kc)
      enddo
    enddo
    
    ! Imag part
#ifdef OverWriteFFT
    call transpose_x1_to_z2(arrx_xhzf,arrz_xhzf,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrz_xhzf,arrz_xhzf)
#else
    allocate(arrFFT(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    call transpose_x1_to_z2(arrx_xhzf,arrFFT,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrFFT,arrz_xhzf)
    deallocate(arrFFT)
#endif
    deallocate(arrx_xhzf)
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc )*arrz_xhzf(ic,jc,kc )+ &
                               arrz_xhzf(ic,jc,kct)*arrz_xhzf(ic,jc,kct)+ arrz_xhzh(ic,jc,kc)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf(ic,jc,kc)*arrz_xhzf(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
      enddo
    enddo
    
    ! LCS
    allocate(arry_xhzf(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(2),decomp_xhzf%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzf,arry_xhzf,decomp_xhzf)
    IF(FlowType==FT_HC) THEN
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y2sz(3),1))
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,decomp_xhzf%y2sz(2)
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
      enddo
    else
      allocate(arrYplane(decomp_xhzf%y2sz(1),decomp_xhzf%y1sz(3),2))
      jct=nyc+1-jForLCS(1)
      do kc=1,decomp_xhzf%y2sz(3)
        do ic=1,decomp_xhzf%y2sz(1)
          arrYplane(ic,kc,1)=arry_xhzf(ic,jForLCS(1),kc)
          arrYplane(ic,kc,2)=arry_xhzf(ic,jct,kc)
        enddo
      enddo
      do kc=1,decomp_xhzf%y2sz(3)
        do jc=1,nyc/2
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,1)
          enddo
        enddo
        do jc=nyc/2+1,nyc
          do ic=1,decomp_xhzf%y2sz(1)      
            arry_xhzf(ic,jc,kc)=arry_xhzf(ic,jc,kc)*arrYplane(ic,kc,2)
          enddo
        enddo        
      enddo
    ENDIF
    call transpose_y2_to_z2(arry_xhzf,arrz_xhzf,decomp_xhzf)
    deallocate(arry_xhzf,arrYplane)
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh_LCS(ic,jc,kc)= arrz_xhzh_LCS(ic,jc,kc)+arrz_xhzf(ic,jc,kc )+ arrz_xhzf(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh_LCS(ic,jc,kc)=arrz_xhzh_LCS(ic,jc,kc)+ arrz_xhzf(ic,jc,kc)
        kc=nzhp; arrz_xhzh_LCS(ic,jc,kc)=arrz_xhzh_LCS(ic,jc,kc)+ arrz_xhzf(ic,jc,kc)
      enddo
    enddo
    
    deallocate(arrz_xhzf)
    allocate(arry_xhzh(decomp_xhzh%y2sz(1),decomp_xhzh%y2sz(2),decomp_xhzh%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzh,arry_xhzh,decomp_xhzh)
    deallocate(arrz_xhzh)

    IF(FlowType==FT_HC) THEN
      do jc=1,nyc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m1)=EnergySpec2D(ic,jc,kc,m1)+normEgy*arry_xhzh(ic,jc,kc)
          enddo
        enddo
      enddo
    ELSE 
      do jc=1,nyc/2
        jcs=nyc+1-jc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m1)=EnergySpec2D(ic,jc,kc,m1)+ &
                                      normEgy*(arry_xhzh(ic,jc,kc)+arry_xhzh(ic,jcs,kc))
          enddo
        enddo
      enddo
    ENDIF
    
    ! LCS
    call transpose_z2_to_y2(arrz_xhzh_LCS,arry_xhzh,decomp_xhzh)
    deallocate(arrz_xhzh_LCS)
    IF(FlowType==FT_HC) THEN
      do jc=1,nyc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m2)=EnergySpec2D(ic,jc,kc,m2)+normEgy*arry_xhzh(ic,jc,kc)
          enddo
        enddo
      enddo
    ELSE 
      do jc=1,nyc/2
        jcs=nyc+1-jc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m2)=EnergySpec2D(ic,jc,kc,m2)+ &
                                      normEgy*(arry_xhzh(ic,jc,kc)+arry_xhzh(ic,jcs,kc))
          enddo
        enddo
      enddo
    ENDIF 
    deallocate(arry_xhzh)
  end subroutine clcSpecAndLCSR2DFrom1D
  
  !******************************************************************
  ! clcCosSpectra2D
  !******************************************************************
  subroutine clcCosSpectra2D(ArrIN1,ArrIN2,m,IsSymmetry)
    implicit none
    integer,intent(in)::m
    logical,intent(in)::IsSymmetry
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in)::ArrIn1,ArrIN2
    
    ! locals
    real(RK)::normEgy,rcoe
    integer::ic,jc,kc,ict,kct,jcs
    real(RK),dimension(:,:,:),allocatable::arry_xhzh
    real(RK),dimension(:,:,:),allocatable::arrx,arrx_xhzf1,arrx_xhzf2
    real(RK),dimension(:,:,:),allocatable::arrz_xhzf1,arrz_xhzf2,arrz_xhzh
#ifndef OverWriteFFT
    real(RK),dimension(:,:,:),allocatable::ArrFFT
#endif
    
    normEgy= one/(real(nxc,RK)*real(nzc,RK))
    normEgy= normEgy*normEgy
    IF(FlowType==FT_CH)normEgy=normEgy*half
    allocate(arrx(nxc,x1size(2),x1size(3)))
#ifdef OverWriteFFT
    call transpose_y1_to_x1(ArrIN1,arrx)
    call dfftw_execute_r2r(fft_plan_x, arrx, arrx)
#else
    allocate(ArrFFT(nxc,x1size(2),x1size(3)),arrx(nxc,x1size(2),x1size(3)))
    call transpose_y1_to_x1(ArrIN1,ArrFFT)
    call dfftw_execute_r2r(fft_plan_x, ArrFFT, arrx)
    deallocate(ArrFFT)
#endif
    allocate(arrx_xhzf1(nxhp,x1size(2),x1size(3)),   &
             arrz_xhzf1(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
                          
    ! Real part
    arrx_xhzf1=arrx(1:nxhp,:,:)
#ifdef OverWriteFFT
    call transpose_x1_to_z2(arrx_xhzf1,arrz_xhzf1,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrz_xhzf1,arrz_xhzf1)
#else
    allocate(arrFFT(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    call transpose_x1_to_z2(arrx_xhzf1,arrFFT,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrFFT,arrz_xhzf1)
    deallocate(arrFFT)
#endif   
    arrx_xhzf1(1,:,:)=zero; arrx_xhzf1(nxhp,:,:)=zero
    do ic=2,nxh
      ict=nxc+2-ic
      arrx_xhzf1(ic,:,:)=arrx(ict,:,:)
    enddo
#ifdef OverWriteFFT
    call transpose_y1_to_x1(ArrIN2,arrx)
    call dfftw_execute_r2r(fft_plan_x, arrx, arrx)
#else
    allocate(ArrFFT(nxc,x1size(2),x1size(3)),arrx(nxc,x1size(2),x1size(3)))
    call transpose_y1_to_x1(ArrIN2,ArrFFT)
    call dfftw_execute_r2r(fft_plan_x, ArrFFT, arrx)
    deallocate(ArrFFT)
#endif
    allocate(arrx_xhzf2(nxhp,x1size(2),x1size(3)),   &
             arrz_xhzf2(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    arrx_xhzf2=arrx(1:nxhp,:,:)
#ifdef OverWriteFFT
    call transpose_x1_to_z2(arrx_xhzf2,arrz_xhzf2,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrz_xhzf2,arrz_xhzf2)
#else
    allocate(arrFFT(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    call transpose_x1_to_z2(arrx_xhzf2,arrFFT,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrFFT,arrz_xhzf2)
    deallocate(arrFFT)
#endif        
    arrx_xhzf2(1,:,:)=zero; arrx_xhzf2(nxhp,:,:)=zero
    do ic=2,nxh
      ict=nxc+2-ic
      arrx_xhzf2(ic,:,:)=arrx(ict,:,:)
    enddo
    deallocate(arrx)
    allocate(arrz_xhzh(decomp_xhzh%z2sz(1),decomp_xhzh%z2sz(2),decomp_xhzh%z2sz(3)))
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc )*arrz_xhzf2(ic,jc,kc )+ &
                               arrz_xhzf1(ic,jc,kct)*arrz_xhzf2(ic,jc,kct)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc)*arrz_xhzf2(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc)*arrz_xhzf2(ic,jc,kc) 
      enddo
    enddo    
 
    ! Imag part   
#ifdef OverWriteFFT
    call transpose_x1_to_z2(arrx_xhzf1,arrz_xhzf1,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrz_xhzf1,arrz_xhzf1)  
    deallocate(arrx_xhzf1)
    call transpose_x1_to_z2(arrx_xhzf2,arrz_xhzf2,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrz_xhzf2,arrz_xhzf2)  
    deallocate(arrx_xhzf2)
#else
    allocate(arrFFT(decomp_xhzf%z2sz(1),decomp_xhzf%z2sz(2),decomp_xhzf%z2sz(3)))
    call transpose_x1_to_z2(arrx_xhzf1,arrFFT,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrFFT,arrz_xhzf1)  
    deallocate(arrx_xhzf1)
    call transpose_x1_to_z2(arrx_xhzf2,arrFFT,decomp_xhzf)
    call dfftw_execute_r2r(fft_plan_z2,arrFFT,arrz_xhzf2)  
    deallocate(arrx_xhzf2,arrFFT)
#endif    
    
    do kc=2,nzh
      kct=nzc+2-kc
      do jc=1,decomp_xhzh%z2sz(2)
        do ic=1,decomp_xhzh%z2sz(1)
          arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc )*arrz_xhzf2(ic,jc,kc )+ &
                               arrz_xhzf1(ic,jc,kct)*arrz_xhzf2(ic,jc,kct)+ arrz_xhzh(ic,jc,kc)
        enddo
      enddo
    enddo
    do jc=1,decomp_xhzh%z2sz(2)
      do ic=1,decomp_xhzh%z2sz(1)
        kc=1;    arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc)*arrz_xhzf2(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
        kc=nzhp; arrz_xhzh(ic,jc,kc)= arrz_xhzf1(ic,jc,kc)*arrz_xhzf2(ic,jc,kc)+ arrz_xhzh(ic,jc,kc)
      enddo
    enddo
    deallocate(arrz_xhzf1,arrz_xhzf2)    
    allocate(arry_xhzh(decomp_xhzh%y2sz(1),decomp_xhzh%y2sz(2),decomp_xhzh%y2sz(3)))
    call transpose_z2_to_y2(arrz_xhzh,arry_xhzh,decomp_xhzh)
    deallocate(arrz_xhzh)
            
    IF(FlowType==FT_HC) THEN
      do jc=1,nyc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m)=EnergySpec2D(ic,jc,kc,m)+normEgy*arry_xhzh(ic,jc,kc)
          enddo
        enddo
      enddo
    ELSE
      if(IsSymmetry) then
        rCoe= one      
      else
        rCoe=-one
      endif 
      do jc=1,nyc/2
        jcs=nyc+1-jc
        do kc=1,decomp_xhzh%y2sz(3)
          do ic=1,decomp_xhzh%y2sz(1)
            EnergySpec2D(ic,jc,kc,m)=EnergySpec2D(ic,jc,kc,m)+ &
                                     normEgy*(arry_xhzh(ic,jc,kc)+rCoe*arry_xhzh(ic,jcs,kc))
          enddo
        enddo
      enddo
    ENDIF
    deallocate(arry_xhzh)     
  end subroutine clcCosSpectra2D
  
  !******************************************************************
  ! transform_uy
  !******************************************************************
  subroutine transform_uy(uyIn,uyOut)
    implicit none
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in) ::uyIn
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::uyOut
    
    ! locals
    integer::ic,jc,kc,jp
    
    if(FlowType==FT_CH) then
      DO kc=1,y1size(3)
        do jc=1,nyc
          jp=jc+1
          if(jc>nyc/2)jp=jc
          do ic=1,y1size(1)
            uyOut(ic,jc,kc)=uyIn(ic,jp,kc)
          enddo
        enddo
      ENDDO    
    else
      DO kc=1,y1size(3)
        do jc=1,nyc-1
          jp=jc+1
          do ic=1,y1size(1)
            uyOut(ic,jc,kc)=uyIn(ic,jp,kc)
          enddo
        enddo
      ENDDO
      uyOut(:,nyc,:)=zero    
    endif
  end subroutine transform_uy

#undef iSpec1DUU
#undef iSpec1DVV
#undef iSpec1DWW
#undef iSpec1DPP
#undef iSpec1DUV
#undef iSpec1DUV2
#undef iLCSR1DUU
#undef iLCSI1DUU
#undef iLCSR1DUU2
#undef iLCSR1DVV
#undef iLCSI1DVV
#undef iLCSR1DVV2
#undef iLCSR1DWW
#undef iLCSI1DWW
#undef iLCSR1DWW2
#undef iLCSR1DPP
#undef iLCSI1DPP
#undef iLCSR1DPP2
#undef iSpec1DCC
#undef iSpec1DUC
#undef iSpec1DVC
#undef iLCSR1DCC
#undef iLCSI1DCC
#undef iLCSR1DCC2

#undef iSpec2DUU
#undef iSpec2DVV
#undef iSpec2DWW
#undef iSpec2DPP
#undef iSpec2DUV
#undef iLCSR2DUU
#undef iLCSR2DVV
#undef iLCSR2DWW
#undef iLCSR2DPP
#undef iSpec2DCC
#undef iLCSR2DCC

#undef NCHASTAT
#undef NEnergySpec1D
#undef NEnergySpec2D
#ifdef SAVE_SINGLE_Spec2D
#undef SAVE_SINGLE_Spec2D
#endif

#ifdef HighOrderGradStat
#undef HighOrderGradStat
#undef nGradStat
#endif
end module m_FlowCase
