module m_FlowType_Channel
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use iso_c_binding
  use m_MeshAndMetries
  use m_Variables,only:mb1
  use m_Tools,only:CalcUxAver
#ifdef CFDLPT_TwoWay
  use m_Variables,only:FpForce_x,FpForce_y,FpForce_z
#endif
#ifdef CFDACM
  use m_Variables,only:FluidIndicator
#endif
  implicit none
  private
  include "fftw3.f03"
#define HighOrderGradStat
#define SAVE_SINGLE_Spec2D

  ! statistics variabls
  integer:: nfstime
  real(RK):: PrGradsum
  real(RK),allocatable,dimension(:,:):: SumStat
#ifdef HighOrderGradStat 
  real(RK),allocatable,dimension(:,:)::SumGrad
#endif

  ! SpectraOptions
  integer:: ivSpec,jForLCS(2)
  logical:: clcSpectra1D,clcSpectra2D
  
  ! Spectra variables
  type(C_PTR)::fft_plan_x,fft_plan_z
  integer:: nxh,nxhp,nzh,nzhp,nSpectime
  type(decomp_info),allocatable::decomp_xhzf,decomp_xhzh
  real(RK),allocatable,dimension(:,:,:,:)::EnergySpec2D
  real(RK),allocatable,dimension(:,:,:)::EnergySpecX,EnergySpecZ
    
  public:: InitVelocity_CH,InitStatVar_CH,clcStat_CH

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
#define iImag1DUC  21
#define iSpec1DVC  22
#define iImag1DVC  23
#define iLCSR1DCC  24
#define iLCSI1DCC  25
#define iLCSR1DCC2 26

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

#ifdef CFDLPT_TwoWay
#define NCHASTAT 45
#elif CFDACM
#define NCHASTAT 51
#else
#define NCHASTAT 35
#endif

#define NEnergySpec1D 18
#define NEnergySpec2D 9

#if defined(HighOrderGradStat)
#define nGradStat     98
#endif

contains
#include "my_FFTW_inc.f90"
#define EnergySpectra_staggered_2nd
#include "EnergySpectra_calcu_fun_inc.f90"
#undef  EnergySpectra_staggered_2nd

  !******************************************************************
  ! InitVelocity_CH
  !******************************************************************
  subroutine InitVelocity_CH(ux,uy,uz,Deviation,chFile)
    implicit none
    character(*),intent(in)::chFile
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
  
    ! locals
    integer::ierror,ic,jc,kc,m1,m2,nUnit,iTV(8)
    real(RK)::uBulkTemp,ybulk1,ybulk2,twopi
    real(RK)::xplus,yplus,zplus,yct,ybar,xpt,zpt,ratiot
    real(RK)::retau_guass,utau_guass,height,rem,wx,wz,xlxPlus,zlzPlus
    NAMELIST/uBulk_Param/ybulk1,ybulk2

    open(newunit=nUnit, file=chFile, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitVelocity_CH","Cannot open file: "//trim(chFile))
    endif
    read(nUnit,nml=uBulk_Param)
    close(nUnit,IOSTAT=ierror)
    height=ybulk2-ybulk1
    uBulkTemp=uBulk*yly/height
    if(FlowType==FT_CH) then
      height=0.5_RK*height
    endif
    rem=uBulkTemp*height/xnu
    ux=0.0_RK; uy=0.0_RK; uz=0.0_RK
    if(abs(uBulk)<1.0E-12)return
    
    retau_guass = 0.1538_RK*rem**0.887741_RK
    utau_guass  = retau_guass*xnu/height
    if(nrank==0)print*,'************** retau_gauss=',retau_guass
    if(nrank==0)print*,'************** utau_gauss= ',utau_guass

    call date_and_time(values=iTV); !iTV=0
    call random_seed(size= ic)
    call random_seed(put = iTV(7)*iTV(8)+[(jc,jc=1,ic)])
    call random_number(Deviation)
    Deviation= 0.2_RK* Deviation + 0.9_RK ! [0.8, 1.2]

    ! modulation of the random noise + initial velocity profile
    twopi=2.0_RK*PI
    wx=twopi/500.0_RK; wz=twopi/200.0_RK
    xlxPlus=xlx*utau_guass/xnu;   zlzPlus=zlz*utau_guass/xnu;
    m1=floor(xlxPlus*wx/twopi)+1; wx=real(m1,RK)*twopi/xlxPlus
    m2=floor(zlzPlus*wz/twopi)+1; wz=real(m2,RK)*twopi/zlzPlus
    do jc=y1start(2),y1end(2)
      yct = height-abs(height-(yc(jc)-ybulk1))
      if(yc(jc)<ybulk1 .or. yc(jc)>=ybulk2)yct=0.0_RK
      ybar= yct/height; yplus=utau_guass*yct/xnu
      do kc=y1start(3),y1end(3)
        zpt  =real(kc-1,kind=RK)*dz+dz*0.5_RK
        zplus=utau_guass*zpt/xnu
        do ic=y1start(1),y1end(1)
          xpt  =real(ic-1,kind=RK)*dx+dx*0.5_RK
          xplus=utau_guass*xpt/xnu
          !ux(ic,jc,kc) = 0.0052_RK*uBulkTemp*yplus*exp(-yplus*yplus/1800.0_RK)*cos(wz*zplus)*Deviation(ic,jc,kc) ! original expression
          !uz(ic,jc,kc) = 0.0050_RK*uBulkTemp*yplus*exp(-yplus*yplus/1800.0_RK)*sin(wx*xplus)*Deviation(ic,jc,kc) ! original expression
          ux(ic,jc,kc) = uBulkTemp*ybar*exp(-4.5_RK*ybar*ybar)*cos(wz*zplus)*Deviation(ic,jc,kc)
          uz(ic,jc,kc) = uBulkTemp*ybar*exp(-4.5_RK*ybar*ybar)*sin(wx*xplus)*Deviation(ic,jc,kc)
          ux(ic,jc,kc) = ux(ic,jc,kc)+ 3.0_RK*uBulkTemp*(ybar-0.5_RK*ybar*ybar)
        enddo
      enddo
    enddo
    ratiot=ubulk/CalcUxAver(ux)
    call MPI_Bcast(ratiot,1,real_type,0,MPI_COMM_WORLD,ierror)
    ux= ux*ratiot
    Deviation=0.0_RK
  end subroutine InitVelocity_CH

  !******************************************************************
  ! InitStatVar_CH
  !******************************************************************
  subroutine InitStatVar_CH(chFile)
    implicit none
    character(*),intent(in)::chFile

    ! locals
    real(RK)::ybulk1,ybulk2
    character(len=128)::filename
    NAMELIST/uBulk_Param/ybulk1,ybulk2
    integer::ierror,nUnit,plan_type,nySpec2D
    real(RK),dimension(:),allocatable::Vec1,Vec2
    NAMELIST/SpectraOptions/clcSpectra1D,clcSpectra2D,ivSpec,jForLCS

    open(newunit=nUnit, file=chFile, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","Cannot open file: "//trim(chFile))
    endif
    read(nUnit,nml=uBulk_Param)
    rewind(nUnit)
    read(nUnit, nml=SpectraOptions)
    close(nUnit,IOSTAT=ierror)
        
    if(nrank==0) then
      write(MainLog%nUnit, nml=uBulk_Param)
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
    if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","Allocation failed")  
    nfstime=0; SumStat=0.0_RK; PrGradsum=0.0_RK
#ifdef HighOrderGradStat
    allocate(SumGrad(nGradStat,nyp),Stat=ierror)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar: ","Allocation failed 2")  
    SumGrad=0.0_RK
#endif

    if(FFTW_plan_type == 1) then
      plan_type=FFTW_MEASURE
    else
      plan_type=FFTW_ESTIMATE
    endif    
    nxh=nxc/2; nxhp=nxh+1
    nzh=nzc/2; nzhp=nzh+1
    
    ! FFT_plan_x
    allocate(Vec1(x1size(1)),Vec2(x1size(1)))
    fft_plan_x = fftw_plan_r2r_1d(x1size(1),Vec1,Vec2,FFTW_R2HC,plan_type)         
    deallocate(Vec1,Vec2)

    ! FFT_plan_z
    allocate(Vec1(z1size(3)),Vec2(z1size(3)))
    fft_plan_z = fftw_plan_r2r_1d(z1size(3),Vec1,Vec2,FFTW_R2HC,plan_type)         
    deallocate(Vec1,Vec2)
    
    IF(clcSpectra1D) THEN
      allocate(EnergySpecX(nxhp,x1size(2),NEnergySpec1D));EnergySpecX=0.0_RK
      allocate(EnergySpecZ(nzhp,z1size(2),NEnergySpec1D));EnergySpecZ=0.0_RK
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
 
      if(FlowType==FT_CH) then
        nySpec2D=nyc/2
      else
        nySpec2D=nyc
      endif
      allocate(EnergySpec2D(decomp_xhzh%y2sz(1),nySpec2D,decomp_xhzh%y2sz(3),NEnergySpec2D),Stat=ierror)
      if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar","Allocation failed For Spectra2D")
      EnergySpec2D=0.0_RK
    ENDIF
  end subroutine InitStatVar_CH

  !******************************************************************
  ! clcStat_CH
  !******************************************************************
  subroutine clcStat_CH(ux,uy,uz,pr,ArrTemp1,ArrTemp2)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pr
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::ArrTemp1,ArrTemp2
   
    ! locals
    character(len=128)::filename
    integer(kind=8)::disp,disp_inc
    real(RK),allocatable,dimension(:,:,:)::arrx1,arrx2,arrz1,arrz2
    real(RK),allocatable,dimension(:,:,:)::EnergySpecXR,EnergySpecZR
    real(RK)::dudyU,dvdyM,dudzC,dvdzC,dudxC,dvdxU,SumVec(NCHASTAT),rdxh
    integer::ic,jc,kc,im,jm,km,ip,jp,kp,it,jt,kt,ierror,ids,nrankX,nrankZ,nUnit
    real(RK)::dudxx,dudyy,dudzz,dvdxx,dvdyy,dvdzz,dwdxx,dwdyy,dwdzz,SumStatR(NCHASTAT,nyp)
    real(RK)::uxloc,uyloc,uzloc,prCell,uxCell,uyCell,uzCell,prloc2,inxz,infstime,cac,caj,cacU
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vor_x,vor_y,vor_z,InterpY1,InterpY2
#ifdef HighOrderGradStat
    real(RK)::dpdx,dpdy,dpdz,dpdxx,dpdyy,dpdzz,SumGradR(nGradStat,nyp),SumVec2(nGradStat)
#endif
#ifdef CFDLPT_TwoWay
    real(RK)::uyCellm
#endif

    rdxh=0.5_RK*rdx
    inxz=1.0_RK/(real(nxc,RK)*real(nzc,RK))
    DO jc=1,nyc
      jp=jc+1; jm=jc-1;
      InterpY1= 0.5_RK*YinterpCoe(jc); InterpY2=0.5_RK-InterpY1
      cac=rdyc(jc);cacU=rdyc(jp); caj=rdyp(jc); SumVec=0.0_RK
#ifdef HighOrderGradStat
      SumVec2=0.0_RK
#endif
      do kc=y1start(3),y1end(3)
        km=kc-1;kp=kc+1
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1
          
          uxloc = ux(ic,jc,kc)
          uyloc = uy(ic,jc,kc)
          uzloc = uz(ic,jc,kc)
          prCell= pr(ic,jc,kc)
          uxCell= 0.5_RK*(ux(ic,jc,kc)+ux(ip,jc,kc))
          uyCell= 0.5_RK*(uy(ic,jc,kc)+uy(ic,jp,kc))
          uzCell= 0.5_RK*(uz(ic,jc,kc)+uz(ic,jc,kp))
          prloc2= InterpY1*(pr(im,jm,kc)+pr(ic,jm,kc))+ InterpY2*(pr(im,jc,kc)+pr(ic,jc,kc))
 
          dudx= (ux(ip,jc,kc)-ux(ic,jc,kc))*rdx
          dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
          dudz= (ux(ic,jc,kc)-ux(ic,jc,km))*rdz
          dvdx= (uy(ic,jc,kc)-uy(im,jc,kc))*rdx
          dvdy= (uy(ic,jp,kc)-uy(ic,jc,kc))*caj
          dvdz= (uy(ic,jc,kc)-uy(ic,jc,km))*rdz
          dwdx= (uz(ic,jc,kc)-uz(im,jc,kc))*rdx
          dwdy= (uz(ic,jc,kc)-uz(ic,jm,kc))*cac
          dwdz= (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz

          dudxx= (ux(ip,jc,kc)-2.0_RK*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          dudyy= ap2c(jc)*ux(ic,jp,kc)+ac2c(jc)*ux(ic,jc,kc)+am2c(jc)*ux(ic,jm,kc)
          dudzz= (ux(ic,jc,kp)-2.0_RK*ux(ic,jc,kc)+ux(ic,jc,km))*rdz2
          dvdxx= (uy(ip,jc,kc)-2.0_RK*uy(ic,jc,kc)+uy(im,jc,kc))*rdx2
          dvdyy= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          dvdzz= (uy(ic,jc,kp)-2.0_RK*uy(ic,jc,kc)+uy(ic,jc,km))*rdz2
          dwdxx= (uz(ip,jc,kc)-2.0_RK*uz(ic,jc,kc)+uz(im,jc,kc))*rdx2
          dwdyy= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)
          dwdzz= (uz(ic,jc,kp)-2.0_RK*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          vor_x= dwdy-dvdz
          vor_y= dudz-dwdx
          vor_z= dvdx-dudy

          dudxC= (ux(ip,jc,kc)-ux(im,jc,kc))*rdxh
          dvdxU= (uy(ic,jp,kc)-uy(im,jp,kc))*rdx
          dudyU= (ux(ic,jp,kc)-ux(ic,jc,kc))*cacU
          dvdyM= (uy(im,jp,kc)-uy(im,jc,kc))*caj
          dudzC= (InterpY1*(ux(ic,jm,kp)-ux(ic,jm,km)) +InterpY2*(ux(ic,jc,kp)-ux(ic,jc,km)))*rdz
          dvdzC= (uy(im,jc,kp)-uy(im,jc,km)+ uy(ic,jc,kp)-uy(ic,jc,km))*rdz*0.25_RK

          ids=1
          SumVec(ids)=SumVec(ids)+ uxloc;                     ids=ids+1 ! 01
          SumVec(ids)=SumVec(ids)+ uyloc;                     ids=ids+1 ! 02, *yp*
          SumVec(ids)=SumVec(ids)+ uzloc;                     ids=ids+1 ! 03
          SumVec(ids)=SumVec(ids)+ prCell;                    ids=ids+1 ! 04
          SumVec(ids)=SumVec(ids)+ uxloc*uxloc;               ids=ids+1 ! 05
          SumVec(ids)=SumVec(ids)+ uyloc*uyloc;               ids=ids+1 ! 06, *yp*  
          SumVec(ids)=SumVec(ids)+ uzloc*uzloc;               ids=ids+1 ! 07
          SumVec(ids)=SumVec(ids)+ prCell*prCell;             ids=ids+1 ! 08
          SumVec(ids)=SumVec(ids)+ uxCell*uyCell;             ids=ids+1 ! 09 
          SumVec(ids)=SumVec(ids)+ uyCell*uzCell;             ids=ids+1 ! 10
          SumVec(ids)=SumVec(ids)+ uxCell*uzCell;             ids=ids+1 ! 11
          SumVec(ids)=SumVec(ids)+ uxCell*prCell;             ids=ids+1 ! 12
          SumVec(ids)=SumVec(ids)+ uyCell*prCell;             ids=ids+1 ! 13
          SumVec(ids)=SumVec(ids)+ uzCell*prCell;             ids=ids+1 ! 14
          SumVec(ids)=SumVec(ids)+ uxCell*uxCell*uyCell;      ids=ids+1 ! 15
          SumVec(ids)=SumVec(ids)+ uyloc *uyloc *uyloc;       ids=ids+1 ! 16, *yp*
          SumVec(ids)=SumVec(ids)+ uzCell*uzCell*uyCell;      ids=ids+1 ! 17
          SumVec(ids)=SumVec(ids)+ uxCell*uyCell*uyCell;      ids=ids+1 ! 18
          SumVec(ids)=SumVec(ids)+ uxloc*uxloc*uxloc;         ids=ids+1 ! 19 
          SumVec(ids)=SumVec(ids)+ uzloc*uzloc*uzloc;         ids=ids+1 ! 20
          SumVec(ids)=SumVec(ids)+ uxloc*uxloc*uxloc*uxloc;   ids=ids+1 ! 21
          SumVec(ids)=SumVec(ids)+ uyloc*uyloc*uyloc*uyloc;   ids=ids+1 ! 22, *yp*
          SumVec(ids)=SumVec(ids)+ uzloc*uzloc*uzloc*uzloc;   ids=ids+1 ! 23
          SumVec(ids)=SumVec(ids)+ prCell*dudx;               ids=ids+1 ! 24 
          SumVec(ids)=SumVec(ids)+ prCell*dvdy;               ids=ids+1 ! 25 
          SumVec(ids)=SumVec(ids)+ prCell*dwdz;               ids=ids+1 ! 26 
          SumVec(ids)=SumVec(ids)+ prloc2*(dudy+dvdx);        ids=ids+1 ! 27, *yp*
          SumVec(ids)=SumVec(ids)+ uxloc*(dudxx+dudyy+dudzz); ids=ids+1 ! 28 
          SumVec(ids)=SumVec(ids)+ uyloc*(dvdxx+dvdyy+dvdzz); ids=ids+1 ! 29, *yp*
          SumVec(ids)=SumVec(ids)+ uzloc*(dwdxx+dwdyy+dwdzz); ids=ids+1 ! 30
          SumVec(ids)=SumVec(ids)+ dudxC*(dvdxU+dvdx)*0.5_RK  &
                     +(dudyU+dudy)*(dvdyM+dvdy)*0.25_RK;      ids=ids+1 ! 31
          SumVec(ids)=SumVec(ids)+ dudzC*dvdzC;               ids=ids+1 ! 32, *yp*
          SumVec(ids)=SumVec(ids)+ vor_x*vor_x;               ids=ids+1 ! 33, *yp*
          SumVec(ids)=SumVec(ids)+ vor_y*vor_y;               ids=ids+1 ! 34
          SumVec(ids)=SumVec(ids)+ vor_z*vor_z;               ids=ids+1 ! 35, *yp*
#ifdef CFDLPT_TwoWay
          SumVec(ids)=SumVec(ids)+ FpForce_x(ic,jc,kc);       ids=ids+1 ! 36
          SumVec(ids)=SumVec(ids)+ FpForce_y(ic,jc,kc);       ids=ids+1 ! 37, *yp*
          SumVec(ids)=SumVec(ids)+ FpForce_z(ic,jc,kc);       ids=ids+1 ! 38
          SumVec(ids)=SumVec(ids)+ FpForce_x(ic,jc,kc)*FpForce_x(ic,jc,kc); ids=ids+1 ! 39
          SumVec(ids)=SumVec(ids)+ FpForce_y(ic,jc,kc)*FpForce_y(ic,jc,kc); ids=ids+1 ! 40, *yp*
          SumVec(ids)=SumVec(ids)+ FpForce_z(ic,jc,kc)*FpForce_z(ic,jc,kc); ids=ids+1 ! 41          
          SumVec(ids)=SumVec(ids)+ FpForce_x(ic,jc,kc)*uxloc; ids=ids+1 ! 42
          SumVec(ids)=SumVec(ids)+ FpForce_y(ic,jc,kc)*uyloc; ids=ids+1 ! 43, *yp*
          SumVec(ids)=SumVec(ids)+ FpForce_z(ic,jc,kc)*uzloc; ids=ids+1 ! 44
          
          SumVec(ids)=SumVec(ids)+ 0.5_RK*(FpForce_y(ic,jc,kc)+FpForce_y(ic,jp,kc))*uxCell
          uyCellm=(uy(im,jc,kc)+uy(im,jp,kc))*0.5_RK
          SumVec(ids)=SumVec(ids)+FpForce_x(ic,jc,kc)*(uyCellm+uyCell)*0.5_RK; ids=ids+1 ! 45
#endif

#ifdef CFDACM
          SumVec(ids)=SumVec(ids)+ uxCell;                    ids=ids+1 ! 36
          SumVec(ids)=SumVec(ids)+ uyCell;                    ids=ids+1 ! 37
          SumVec(ids)=SumVec(ids)+ uzCell;                    ids=ids+1 ! 38
          SumVec(ids)=SumVec(ids)+ uxCell*uxCell;             ids=ids+1 ! 39
          SumVec(ids)=SumVec(ids)+ uyCell*uyCell;             ids=ids+1 ! 40
          SumVec(ids)=SumVec(ids)+ uzCell*uzCell;             ids=ids+1 ! 41
          if(FluidIndicator(ic,jc,kc)=='P')cycle
          SumVec(ids)=SumVec(ids)+ 1.0_RK;                    ids=ids+1 ! 42          
          SumVec(ids)=SumVec(ids)+ uxCell;                    ids=ids+1 ! 43
          SumVec(ids)=SumVec(ids)+ uyCell;                    ids=ids+1 ! 44
          SumVec(ids)=SumVec(ids)+ uzCell;                    ids=ids+1 ! 45
          SumVec(ids)=SumVec(ids)+ uxCell*uxCell;             ids=ids+1 ! 46
          SumVec(ids)=SumVec(ids)+ uyCell*uyCell;             ids=ids+1 ! 47
          SumVec(ids)=SumVec(ids)+ uzCell*uzCell;             ids=ids+1 ! 48
          SumVec(ids)=SumVec(ids)+ uxCell*uyCell;             ids=ids+1 ! 49
          SumVec(ids)=SumVec(ids)+ prCell;                    ids=ids+1 ! 50
          SumVec(ids)=SumVec(ids)+ prCell*prCell;             ids=ids+1 ! 51
#endif

#ifdef HighOrderGradStat
          dpdx=(pr(ic,jc,kc)-pr(im,jc,kc))*rdx
          dpdy=(pr(ic,jc,kc)-pr(ic,jm,kc))*cac          ! yp ! 
          dpdz=(pr(ic,jc,kc)-pr(ic,jc,km))*rdz
          dpdxx= (pr(ip,jc,kc)-2.0_RK*pr(ic,jc,kc)+pr(im,jc,kc))*rdx2
          dpdyy= ap2c(jc)*pr(ic,jp,kc)+ac2c(jc)*pr(ic,jc,kc)+am2c(jc)*pr(ic,jm,kc)
          dpdzz= (pr(ic,jc,kp)-2.0_RK*pr(ic,jc,kc)+pr(ic,jc,km))*rdz2
          
          ids=1
          SumVec2(ids)=SumVec2(ids)+ prCell*prCell*prCell;        ids=ids+1 ! 01
          SumVec2(ids)=SumVec2(ids)+ prCell*prCell*prCell*prCell; ids=ids+1 ! 02
          
          SumVec2(ids)=SumVec2(ids)+ dudx;                        ids=ids+1 ! 03
          SumVec2(ids)=SumVec2(ids)+ dudx*dudx;                   ids=ids+1 ! 04
          SumVec2(ids)=SumVec2(ids)+ dudx*dudx*dudx;              ids=ids+1 ! 05
          SumVec2(ids)=SumVec2(ids)+ dudx*dudx*dudx*dudx;         ids=ids+1 ! 06
          SumVec2(ids)=SumVec2(ids)+ dudy;                        ids=ids+1 ! 07, *yp*
          SumVec2(ids)=SumVec2(ids)+ dudy*dudy;                   ids=ids+1 ! 08, *yp*
          SumVec2(ids)=SumVec2(ids)+ dudy*dudy*dudy;              ids=ids+1 ! 09, *yp*
          SumVec2(ids)=SumVec2(ids)+ dudy*dudy*dudy*dudy;         ids=ids+1 ! 10, *yp*
          SumVec2(ids)=SumVec2(ids)+ dudz;                        ids=ids+1 ! 11
          SumVec2(ids)=SumVec2(ids)+ dudz*dudz;                   ids=ids+1 ! 12
          SumVec2(ids)=SumVec2(ids)+ dudz*dudz*dudz;              ids=ids+1 ! 13
          SumVec2(ids)=SumVec2(ids)+ dudz*dudz*dudz*dudz;         ids=ids+1 ! 14

          SumVec2(ids)=SumVec2(ids)+ dvdx;                        ids=ids+1 ! 15, *yp*       
          SumVec2(ids)=SumVec2(ids)+ dvdx*dvdx;                   ids=ids+1 ! 16, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdx*dvdx*dvdx;              ids=ids+1 ! 17, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdx*dvdx*dvdx*dvdx;         ids=ids+1 ! 18, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdy;                        ids=ids+1 ! 19                  
          SumVec2(ids)=SumVec2(ids)+ dvdy*dvdy;                   ids=ids+1 ! 20          
          SumVec2(ids)=SumVec2(ids)+ dvdy*dvdy*dvdy;              ids=ids+1 ! 21      
          SumVec2(ids)=SumVec2(ids)+ dvdy*dvdy*dvdy*dvdy;         ids=ids+1 ! 22   
          SumVec2(ids)=SumVec2(ids)+ dvdz;                        ids=ids+1 ! 23, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdz*dvdz;                   ids=ids+1 ! 24, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdz*dvdz*dvdz;              ids=ids+1 ! 25, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dvdz*dvdz*dvdz*dvdz;         ids=ids+1 ! 26, *yp* 
          
          SumVec2(ids)=SumVec2(ids)+ dwdx;                        ids=ids+1 ! 27
          SumVec2(ids)=SumVec2(ids)+ dwdx*dwdx;                   ids=ids+1 ! 28
          SumVec2(ids)=SumVec2(ids)+ dwdx*dwdx*dwdx;              ids=ids+1 ! 29
          SumVec2(ids)=SumVec2(ids)+ dwdx*dwdx*dwdx*dwdx;         ids=ids+1 ! 30
          SumVec2(ids)=SumVec2(ids)+ dwdy;                        ids=ids+1 ! 31, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dwdy*dwdy;                   ids=ids+1 ! 32, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dwdy*dwdy*dwdy;              ids=ids+1 ! 33, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dwdy*dwdy*dwdy*dwdy;         ids=ids+1 ! 34, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dwdz;                        ids=ids+1 ! 35
          SumVec2(ids)=SumVec2(ids)+ dwdz*dwdz;                   ids=ids+1 ! 36
          SumVec2(ids)=SumVec2(ids)+ dwdz*dwdz*dwdz;              ids=ids+1 ! 37
          SumVec2(ids)=SumVec2(ids)+ dwdz*dwdz*dwdz*dwdz;         ids=ids+1 ! 38
          
          SumVec2(ids)=SumVec2(ids)+ dpdx;                        ids=ids+1 ! 39
          SumVec2(ids)=SumVec2(ids)+ dpdx*dpdx;                   ids=ids+1 ! 40
          SumVec2(ids)=SumVec2(ids)+ dpdx*dpdx*dpdx;              ids=ids+1 ! 41
          SumVec2(ids)=SumVec2(ids)+ dpdx*dpdx*dpdx*dpdx;         ids=ids+1 ! 42
          SumVec2(ids)=SumVec2(ids)+ dpdy;                        ids=ids+1 ! 43, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dpdy*dpdy;                   ids=ids+1 ! 44, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dpdy*dpdy*dpdy;              ids=ids+1 ! 45, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dpdy*dpdy*dpdy*dpdy;         ids=ids+1 ! 46, *yp* 
          SumVec2(ids)=SumVec2(ids)+ dpdz;                        ids=ids+1 ! 47
          SumVec2(ids)=SumVec2(ids)+ dpdz*dpdz;                   ids=ids+1 ! 48
          SumVec2(ids)=SumVec2(ids)+ dpdz*dpdz*dpdz;              ids=ids+1 ! 49
          SumVec2(ids)=SumVec2(ids)+ dpdz*dpdz*dpdz*dpdz;         ids=ids+1 ! 50
          
          SumVec2(ids)=SumVec2(ids)+ dudxx;                       ids=ids+1 ! 51
          SumVec2(ids)=SumVec2(ids)+ dudxx*dudxx;                 ids=ids+1 ! 52
          SumVec2(ids)=SumVec2(ids)+ dudxx*dudxx*dudxx;           ids=ids+1 ! 53
          SumVec2(ids)=SumVec2(ids)+ dudxx*dudxx*dudxx*dudxx;     ids=ids+1 ! 54
          SumVec2(ids)=SumVec2(ids)+ dudyy;                       ids=ids+1 ! 55
          SumVec2(ids)=SumVec2(ids)+ dudyy*dudyy;                 ids=ids+1 ! 56
          SumVec2(ids)=SumVec2(ids)+ dudyy*dudyy*dudyy;           ids=ids+1 ! 57
          SumVec2(ids)=SumVec2(ids)+ dudyy*dudyy*dudyy*dudyy;     ids=ids+1 ! 58
          SumVec2(ids)=SumVec2(ids)+ dudzz;                       ids=ids+1 ! 59
          SumVec2(ids)=SumVec2(ids)+ dudzz*dudzz;                 ids=ids+1 ! 60
          SumVec2(ids)=SumVec2(ids)+ dudzz*dudzz*dudzz;           ids=ids+1 ! 61
          SumVec2(ids)=SumVec2(ids)+ dudzz*dudzz*dudzz*dudzz;     ids=ids+1 ! 62
          
          SumVec2(ids)=SumVec2(ids)+ dvdxx;                       ids=ids+1 ! 63, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdxx*dvdxx;                 ids=ids+1 ! 64, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdxx*dvdxx*dvdxx;           ids=ids+1 ! 65, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdxx*dvdxx*dvdxx*dvdxx;     ids=ids+1 ! 66, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdyy;                       ids=ids+1 ! 67, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdyy*dvdyy;                 ids=ids+1 ! 68, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdyy*dvdyy*dvdyy;           ids=ids+1 ! 69, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdyy*dvdyy*dvdyy*dvdyy;     ids=ids+1 ! 70, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdzz;                       ids=ids+1 ! 71, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdzz*dvdzz;                 ids=ids+1 ! 72, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdzz*dvdzz*dvdzz;           ids=ids+1 ! 73, *yp*
          SumVec2(ids)=SumVec2(ids)+ dvdzz*dvdzz*dvdzz*dvdzz;     ids=ids+1 ! 74, *yp*                                 

          SumVec2(ids)=SumVec2(ids)+ dwdxx;                       ids=ids+1 ! 75
          SumVec2(ids)=SumVec2(ids)+ dwdxx*dwdxx;                 ids=ids+1 ! 76
          SumVec2(ids)=SumVec2(ids)+ dwdxx*dwdxx*dwdxx;           ids=ids+1 ! 77
          SumVec2(ids)=SumVec2(ids)+ dwdxx*dwdxx*dwdxx*dwdxx;     ids=ids+1 ! 78
          SumVec2(ids)=SumVec2(ids)+ dwdyy;                       ids=ids+1 ! 79
          SumVec2(ids)=SumVec2(ids)+ dwdyy*dwdyy;                 ids=ids+1 ! 80
          SumVec2(ids)=SumVec2(ids)+ dwdyy*dwdyy*dwdyy;           ids=ids+1 ! 81
          SumVec2(ids)=SumVec2(ids)+ dwdyy*dwdyy*dwdyy*dwdyy;     ids=ids+1 ! 82
          SumVec2(ids)=SumVec2(ids)+ dwdzz;                       ids=ids+1 ! 83
          SumVec2(ids)=SumVec2(ids)+ dwdzz*dwdzz;                 ids=ids+1 ! 84
          SumVec2(ids)=SumVec2(ids)+ dwdzz*dwdzz*dwdzz;           ids=ids+1 ! 85
          SumVec2(ids)=SumVec2(ids)+ dwdzz*dwdzz*dwdzz*dwdzz;     ids=ids+1 ! 86
          
          SumVec2(ids)=SumVec2(ids)+ dpdxx;                       ids=ids+1 ! 87
          SumVec2(ids)=SumVec2(ids)+ dpdxx*dpdxx;                 ids=ids+1 ! 88
          SumVec2(ids)=SumVec2(ids)+ dpdxx*dpdxx*dpdxx;           ids=ids+1 ! 89
          SumVec2(ids)=SumVec2(ids)+ dpdxx*dpdxx*dpdxx*dpdxx;     ids=ids+1 ! 90
          SumVec2(ids)=SumVec2(ids)+ dpdyy;                       ids=ids+1 ! 91
          SumVec2(ids)=SumVec2(ids)+ dpdyy*dpdyy;                 ids=ids+1 ! 92
          SumVec2(ids)=SumVec2(ids)+ dpdyy*dpdyy*dpdyy;           ids=ids+1 ! 93
          SumVec2(ids)=SumVec2(ids)+ dpdyy*dpdyy*dpdyy*dpdyy;     ids=ids+1 ! 94
          SumVec2(ids)=SumVec2(ids)+ dpdzz;                       ids=ids+1 ! 95
          SumVec2(ids)=SumVec2(ids)+ dpdzz*dpdzz;                 ids=ids+1 ! 96
          SumVec2(ids)=SumVec2(ids)+ dpdzz*dpdzz*dpdzz;           ids=ids+1 ! 97
          SumVec2(ids)=SumVec2(ids)+ dpdzz*dpdzz*dpdzz*dpdzz;     ids=ids+1 ! 98
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
    jc=nyp; jm=jc-1; cac=rdyc(jc); SumVec=0.0_RK
#ifdef HighOrderGradStat
    SumVec2=0.0_RK
#endif
    InterpY1= 0.5_RK*YinterpCoe(jc); InterpY2=0.5_RK-InterpY1
    do kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do ic=y1start(1),y1end(1)
        im=ic-1;ip=ic+1
        prloc2= InterpY1*(pr(im,jm,kc)+pr(ic,jm,kc))+ InterpY2*(pr(im,jc,kc)+pr(ic,jc,kc))
        dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
        dwdy= (uz(ic,jc,kc)-uz(ic,jm,kc))*cac
        vor_x=  dwdy
        vor_z= -dudy
        SumVec(27)=SumVec(27)+ prloc2*(dudy+0.0_RK)   ! yp !
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

    ! shear stress and pressure gradient
    PrGradsum   = PrGradsum+ PrGradData(1)
    if(nrank==0 .and. IsUxConst) then
      write(filename,'(A,I10.10)')trim(ResultsDir)//"PrGrad",ilast
      open(newunit=nUnit,file=filename,status='old',position='append',form='formatted',IOSTAT=ierror)
      if(ierror/=0) then
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      else
        write(nUnit,'(I7,2ES24.15)')itime,SimTime,PrGradData(1)
      endif
      close(nUnit,IOSTAT=ierror)
    endif
    nfstime= nfstime + 1
#define EnergySpectra_staggered_2nd
#include "EnergySpectra_staggered_inc.f90"
#undef  EnergySpectra_staggered_2nd
    if(mod(itime,SaveStat)/=0) return

    ! Write statistics
    call MPI_REDUCE(SumStat,SumStatR,NCHASTAT*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      infstime = 1.0_RK/real(nfstime,RK)
      write(filename,"(A,I10.10)") trim(ResultsDir)//'stats',itime
      open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
      IF(ierror/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      ELSE
        write(nUnit,'(a,I7,a,I7,a,I7)')'  The time step range for this fluid statistics is ', &
                                    itime-(nfstime-1)*ivstats, ':', ivstats, ':', itime
        write(nUnit,'(A)')'  '
        if(IsUxConst) then
          write(nUnit,'(A)')'  Constant velocity in x-dir by adding a pressure gradient.'
          write(nUnit,'(A, ES24.15)')'    time averaged pressure gradient is: ',PrGradsum*infstime
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
        End block
      ENDIF
      close(nUnit,IOSTAT=ierror)
    endif
    
#ifdef HighOrderGradStat
    call MPI_REDUCE(SumGrad,SumGradR,nGradStat*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      infstime = 1.0_RK/real(nfstime,RK)
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
    endif
    SumGrad=0.0_RK
#endif

#include "EnergySpectra_write_fun_inc.f90"

    nfstime=0; SumStat=0.0_RK; PrGradsum=0.0_RK; nSpectime=0; 
    if(clcSpectra1D) then
      EnergySpecX=0.0_RK; EnergySpecZ=0.0_RK
    endif
    if(clcSpectra2D) EnergySpec2D=0.0_RK
  end subroutine clcStat_CH

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
#undef iImag1DUC
#undef iSpec1DVC
#undef iImag1DVC
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
end module m_FlowType_Channel
