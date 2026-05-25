#include "definitions_inc.f90"
module f2_FlowType_HIT
  use MPI
  use iso_c_binding
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d
  use f2_Parameters
  use f2_MeshAndMetries
  use f2_Variables,only: mb1
  implicit none
  private
  include "fftw3.f03"

  real(RK):: F0_intensity
  integer:: nfstime, nLenSpec3D 

  ! Spectra variables
  integer,allocatable,dimension(:)::IndXYZ
  type(C_PTR):: fft_plan_x, fft_plan_y, fft_plan_z
  real(RK),allocatable,dimension(:,:)::EnergySpec3D
  
  public:: InitVelocity_HIT, InitStatVar_HIT, clcStat_HIT, add_FlowType_Forcing_HIT

#define NEnergySpec3D 4
contains
#define my_FFTW_inc_add_y
#include "my_FFTW_inc.f90"
#undef  my_FFTW_inc_add_y

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitStatVar_HIT
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitStatVar_HIT(chFile); implicit none
    character(len=*),intent(in)::chFile

    ! locals
    NAMELIST/HIT_Options/F0_intensity
    character(len=:),allocatable::filename
    real(RK),dimension(:),allocatable::Vec1,Vec2
    integer:: ierr,iUnit,plan_type,ic,jc,kc,ix,iy,iz

    open(newunit=iUnit, file=chFile, status='old',form='formatted',IOSTAT=ierr )
    if(ierr/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitVelocity_CH","Cannot open file: "//strip(chFile))
    endif
    read(iUnit,nml=HIT_Options)
    close(iUnit,IOSTAT=ierr)
    
    if(FFTW_plan_type == 1) then
      plan_type=FFTW_MEASURE
    else
      plan_type=FFTW_ESTIMATE
    endif
    
    ! FFT_plan_x
    allocate(Vec1(nxc), Vec2(nxc))
    fft_plan_x = fftw_plan_r2r_1d(nxc, Vec1,Vec2,FFTW_R2HC,plan_type)         
    deallocate(Vec1,Vec2)

    ! FFT_plan_y
    allocate(Vec1(nyc), Vec2(nyc))
    fft_plan_y = fftw_plan_r2r_1d(nyc, Vec1,Vec2,FFTW_R2HC,plan_type)         
    deallocate(Vec1,Vec2)
    
    ! FFT_plan_z
    allocate(Vec1(nzc), Vec2(nzc))
    fft_plan_z = fftw_plan_r2r_1d(nzc, Vec1,Vec2,FFTW_R2HC,plan_type)         
    deallocate(Vec1,Vec2)

    nLenSpec3D = ceiling(real(nxc/2,RK) * sqrt(3.0_RK)) +2
    allocate(IndXYZ(nxc))
    IndXYZ(1)=0
    do ic=2,nxc/2+1
      IndXYZ(ic) = ic-1
      IndXYZ(nxc+2-ic) = ic-1
    enddo
    do kc=y2start(3),y2end(3)
      iz = IndXYZ(kc)*IndXYZ(kc)
      do jc=y2start(2),y2end(2)
        iy = IndXYZ(jc)*IndXYZ(jc) +iz
        do ic=y2start(1),y2end(1)
          ix = IndXYZ(ic)*IndXYZ(ic) +iy
          ix = nint(sqrt(1.0_RK*ix)) +1
          if(ix<1 .or. ix>nLenSpec3D) then
            call MainLog%CheckForError(ErrT_Abort,"InitStatVar_HIT","Error for nLenSpec3D")
          endif
        enddo
      enddo
    enddo
    allocate(EnergySpec3D(NEnergySpec3D, nLenSpec3D),Stat=ierr)
    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_HIT","Allocation failed For Spectra3D")
    nfstime=0; EnergySpec3D=0.0_RK
        
    if(nrank/=0) return
    filename = strip(ResultsDir_) // "HITStat" // int2str(ilast,10)
    open(newunit=iUnit, file=filename,status='replace',form='formatted',IOSTAT=ierr)
    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_HIT","Cannot open file: "//filename)
    close(iUnit,IOSTAT=ierr)
  end subroutine InitStatVar_HIT
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitVelocity_HIT
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitVelocity_HIT(ux,uy,uz,Deviation,chFile); implicit none
    character(len=*),intent(in)::chFile
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
  
    ! locals
    NAMELIST/HIT_Options/F0_intensity
    type(decomp_info)::decomp_HIT
    logical,dimension(6)::initializeIn
    complex(RK)::cexp1,cexp2,alpha,beta
    real(RK),allocatable,dimension(:)::Ekk
    real(RK),allocatable,dimension(:,:,:)::ArrX
    real(RK)::WaveNumX(nxc),WaveNumY(nyc),WaveNumZ(nzc),rndnum(3)
    complex(RK),dimension(:,:,:),allocatable::cux,cuy,cuz,cMatX,cMatZ
    integer::ic,jc,kc,nxh,nxhp,nyh,nyhp,nzh,nzhp,kmax,kWave,iUnit,ierr
    real(RK)::uu0,rkp,rkxyz,rkxy,rkx,rky,rkz,rA,rTmp,phi,Theta1,Theta2,VelRatio

    open(newunit=iUnit, file=chFile, status='old',form='formatted',IOSTAT=ierr )
    if(ierr/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitVelocity_CH","Cannot open file: "//strip(chFile))
    endif
    read(iUnit,nml=HIT_Options)
    close(iUnit,IOSTAT=ierr)
    VelRatio = 0.5_RK*xlx/Pi
    VelRatio = sqrt(F0_intensity *VelRatio)
        
    rkp=20.0_RK
    uu0=1.0_RK/3.0_RK
    nxh=nxc/2; nxhp=nxh+1
    nyh=nyc/2; nyhp=nyh+1
    nzh=nzc/2; nzhp=nzh+1
    do ic=1,nxh
      WaveNumX(ic)= real(ic-1,RK)
    enddo
    do ic=nxhp,nxc
      WaveNumX(ic)=-real(nxc-ic+1,RK)
    enddo
    do jc=1,nyh
      WaveNumY(jc)= real(jc-1,RK)
    enddo
    do jc=nyhp,nyc
      WaveNumY(jc)=-real(nyc-jc+1,RK)
    enddo    
    do kc=1,nzh
      WaveNumZ(kc)= real(kc-1,RK)
    enddo
    do kc=nzhp,nzc
      WaveNumZ(kc)=-real(nzc-kc+1,RK)
    enddo
    kmax=nint(sqrt(1.0_RK*nxh*nxh+nyh*nyh+nzh*nzh))+1
    allocate(Ekk(kmax)); Ekk=0.0_RK
    rA=16.0_RK*sqrt(2.0_RK/Pi)*uu0
    do kc=1,kmax
      rTmp=real(kc,RK)/rkp
      Ekk(kc)=(rA/rkp)*(rTmp**4)*exp(-2.0_RK*rTmp*rTmp)
    enddo
    initializeIn=.false.; initializeIn(5:6)=.true.
    call decomp_info_init(nxhp,nyc,nzc,decomp_HIT,initialize=initializeIn)
    allocate(cux(decomp_HIT%y2st(1):decomp_HIT%y2en(1),decomp_HIT%y2st(2):decomp_HIT%y2en(2),decomp_HIT%y2st(3):decomp_HIT%y2en(3)))
    allocate(cuy(decomp_HIT%y2st(1):decomp_HIT%y2en(1),decomp_HIT%y2st(2):decomp_HIT%y2en(2),decomp_HIT%y2st(3):decomp_HIT%y2en(3)))
    allocate(cuz(decomp_HIT%y2st(1):decomp_HIT%y2en(1),decomp_HIT%y2st(2):decomp_HIT%y2en(2),decomp_HIT%y2st(3):decomp_HIT%y2en(3)))
    cux=cmplx(0.0_RK,0.0_RK,kind=RK); cuy=cmplx(0.0_RK,0.0_RK,kind=RK); cuz=cmplx(0.0_RK,0.0_RK,kind=RK); 
    
    call system_clock(count=ic); ic=0
    call random_seed(size  =jc)
    call random_seed(put   =ic+63946*(/(kc-1,kc=1,jc)/))
    do kc=decomp_HIT%y2st(3),decomp_HIT%y2en(3)
      rkz=WaveNumZ(kc)
      do jc=decomp_HIT%y2st(2),decomp_HIT%y2en(2)
        rky=WaveNumY(jc)
        do ic=decomp_HIT%y2st(1),decomp_HIT%y2en(1)
          rkx=WaveNumX(ic)
          if(ic+jc+kc==3) cycle
          rkxy =sqrt(rkx*rkx+rky*rky)
          rkxyz=sqrt(rkx*rkx+rky*rky+rkz*rkz)
          kWave=nint(rkxyz)
          call random_number(rndnum)
          Phi   =2.0_RK*Pi*rndnum(1)
          Theta1=2.0_RK*Pi*rndnum(2)
          Theta2=2.0_RK*Pi*rndnum(3)
          cexp1=cmplx(cos(Theta1),sin(Theta1),kind=RK)
          cexp2=cmplx(cos(Theta2),sin(Theta2),kind=RK)
          rTmp=sqrt(Ekk(kWave)/(2.0_RK*Pi))/rkxyz
          alpha=rTmp*cexp1*cos(Phi)
          beta =rTmp*cexp2*sin(Phi)
          if(ic+jc==2) then
            cux(ic,jc,kc)= alpha
            cuy(ic,jc,kc)= beta       
          else
            cux(ic,jc,kc)=( alpha*rkxyz*rky +beta*rkx*rkz)/(rkxyz*rkxy)
            cuy(ic,jc,kc)=(-alpha*rkxyz*rkx +beta*rky*rkz)/(rkxyz*rkxy)
          endif
          cuz(ic,jc,kc)=-beta*rkxy/rkxyz
        enddo
      enddo
    enddo
    deallocate(Ekk)
    
    ! ux
    call ifft3d_y(cux,decomp_HIT%y2sz(1),decomp_HIT%y2sz(2),decomp_HIT%y2sz(3))   ! inverse complex FFT in y-dir
    allocate(cMatZ(decomp_HIT%z2sz(1),decomp_HIT%z2sz(2),decomp_HIT%z2sz(3)))
    call cmplx_transpose_y2_to_z2(cux,cMatZ,decomp_HIT)
    deallocate(cux)
    
    call ifft3d_z(cMatZ,decomp_HIT%z2sz(1),decomp_HIT%z2sz(2),decomp_HIT%z2sz(3)) ! inverse complex FFT in z-dir
    allocate(cMatX(decomp_HIT%x1sz(1),decomp_HIT%x1sz(2),decomp_HIT%x1sz(3)))
    call cmplx_transpose_z2_to_x1(cMatZ,cMatX,decomp_HIT)
    deallocate(cMatZ)
    
    allocate(ArrX(x1size(1),x1size(2),x1size(3)))
    call ifft3d_x(cMatX,ArrX,x1size(1),x1size(2),x1size(3))                       ! inverse real FFT in x-dir
    call transpose_x1_to_y1(ArrX,Deviation)
    deallocate(cMatX,ArrX)
    ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))=Deviation*VelRatio
       
    ! uy
    call ifft3d_y(cuy,decomp_HIT%y2sz(1),decomp_HIT%y2sz(2),decomp_HIT%y2sz(3))   ! inverse complex FFT in y-dir
    allocate(cMatZ(decomp_HIT%z2sz(1),decomp_HIT%z2sz(2),decomp_HIT%z2sz(3)))
    call cmplx_transpose_y2_to_z2(cuy,cMatZ,decomp_HIT)
    deallocate(cuy)
    
    call ifft3d_z(cMatZ,decomp_HIT%z2sz(1),decomp_HIT%z2sz(2),decomp_HIT%z2sz(3)) ! inverse complex FFT in z-dir
    allocate(cMatX(decomp_HIT%x1sz(1),decomp_HIT%x1sz(2),decomp_HIT%x1sz(3)))
    call cmplx_transpose_z2_to_x1(cMatZ,cMatX,decomp_HIT)
    deallocate(cMatZ)
    
    allocate(ArrX(x1size(1),x1size(2),x1size(3)))
    call ifft3d_x(cMatX,ArrX,x1size(1),x1size(2),x1size(3))                       ! inverse real FFT in x-dir
    call transpose_x1_to_y1(ArrX,Deviation)
    deallocate(cMatX,ArrX)
    uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))=Deviation*VelRatio
    
    ! uz
    call ifft3d_y(cuz,decomp_HIT%y2sz(1),decomp_HIT%y2sz(2),decomp_HIT%y2sz(3))   ! inverse complex FFT in y-dir
    allocate(cMatZ(decomp_HIT%z2sz(1),decomp_HIT%z2sz(2),decomp_HIT%z2sz(3)))
    call cmplx_transpose_y2_to_z2(cuz,cMatZ,decomp_HIT)
    deallocate(cuz)
    
    call ifft3d_z(cMatZ,decomp_HIT%z2sz(1),decomp_HIT%z2sz(2),decomp_HIT%z2sz(3)) ! inverse complex FFT in z-dir
    allocate(cMatX(decomp_HIT%x1sz(1),decomp_HIT%x1sz(2),decomp_HIT%x1sz(3)))
    call cmplx_transpose_z2_to_x1(cMatZ,cMatX,decomp_HIT)
    deallocate(cMatZ)
    
    allocate(ArrX(x1size(1),x1size(2),x1size(3)))
    call ifft3d_x(cMatX,ArrX,x1size(1),x1size(2),x1size(3))                       ! inverse real FFT in x-dir
    call transpose_x1_to_y1(ArrX,Deviation)
    deallocate(cMatX,ArrX)
    uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))=Deviation*VelRatio
    
    Deviation=0.0_RK
    call decomp_info_finalize(decomp_HIT)
  end subroutine InitVelocity_HIT 
        
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ifft3d_y
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ifft3d_y(cMatY,im1,im2,im3); implicit none
    include "fftw3.f03"
    integer,intent(in)::im1,im2,im3
    complex(RK),dimension(im1,im2,im3),intent(inout)::cMatY
    
    ! locals
    integer::ic,jc,kc
    type(C_PTR)::FFTW_plan
    type(fftw_iodim)::iodim(1),iodim_howmany(1) 
    complex(RK),dimension(:,:),allocatable::MatIn,MatOut
    
    allocate(MatIn (im2,im1))
    allocate(MatOut(im2,im1))
    iodim(1)%n  = im2
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n  = im1
    iodim_howmany(1)%is = im2
    iodim_howmany(1)%os = im2
    FFTW_plan = fftw_plan_guru_dft(1,iodim,1,iodim_howmany,MatIn,MatOut,+1,FFTW_ESTIMATE)
    do kc=1,im3
      do jc=1,im2
        do ic=1,im1  
          MatIn(jc,ic)=cMatY(ic,jc,kc)
        enddo
      enddo
      call dfftw_execute_dft(FFTW_plan,MatIn,MatOut)
      do jc=1,im2
        do ic=1,im1
          cMatY(ic,jc,kc)=MatOut(jc,ic)
        enddo
      enddo
    enddo
    call fftw_destroy_plan(FFTW_plan)
    deallocate(MatIn,MatOut)    
  end subroutine ifft3d_y
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! cmplx_transpose_y2_to_z2
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90  
  subroutine cmplx_transpose_y2_to_z2(cMatY,cMatZ,decomp_HIT); implicit none
    type(decomp_info),intent(in)::decomp_HIT
    complex(RK),dimension(decomp_HIT%y2sz(1),decomp_HIT%y2sz(2),decomp_HIT%y2sz(3)),intent(in) ::cMatY
    complex(RK),dimension(decomp_HIT%z2sz(1),decomp_HIT%z2sz(2),decomp_HIT%z2sz(3)),intent(out)::cMatZ
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(decomp_HIT%y2sz(1),decomp_HIT%y2sz(2),decomp_HIT%y2sz(3))::rMatY
    real(RK),dimension(decomp_HIT%z2sz(1),decomp_HIT%z2sz(2),decomp_HIT%z2sz(3))::rMatZ
        
    rMatY=real(cMatY,RK)
    call transpose_y2_to_z2(rMatY,rMatZ,decomp_HIT)
    do kc=1,decomp_HIT%z2sz(3)
      do jc=1,decomp_HIT%z2sz(2)
        do ic=1,decomp_HIT%z2sz(1)
          cMatZ(ic,jc,kc)=cmplx(rMatZ(ic,jc,kc),0.0_RK,RK)
        enddo
      enddo
    enddo
    
    rMatY=aimag(cMatY)
    call transpose_y2_to_z2(rMatY,rMatZ,decomp_HIT)
    do kc=1,decomp_HIT%z2sz(3)
      do jc=1,decomp_HIT%z2sz(2)
        do ic=1,decomp_HIT%z2sz(1)
          cMatZ(ic,jc,kc)=cmplx(0.0_RK,rMatZ(ic,jc,kc),RK)+cMatZ(ic,jc,kc)
        enddo
      enddo
    enddo    
  end subroutine cmplx_transpose_y2_to_z2

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ifft3d_z
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ifft3d_z(cMatZ,im1,im2,im3); implicit none
    include "fftw3.f03"
    integer,intent(in)::im1,im2,im3
    complex(RK),dimension(im1,im2,im3),intent(inout)::cMatZ
    
    ! locals
    integer::ic,jc,kc
    type(C_PTR)::FFTW_plan
    type(fftw_iodim)::iodim(1),iodim_howmany(1) 
    complex(RK),dimension(:,:),allocatable::MatIn,MatOut
    
    allocate(MatIn (im3,im1))
    allocate(MatOut(im3,im1))
    iodim(1)%n  = im3
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n  = im1
    iodim_howmany(1)%is = im3
    iodim_howmany(1)%os = im3
    FFTW_plan = fftw_plan_guru_dft(1,iodim,1,iodim_howmany,MatIn,MatOut,+1,FFTW_ESTIMATE)
    do jc=1,im2
      do kc=1,im3
        do ic=1,im1
          MatIn(kc,ic)=cMatZ(ic,jc,kc)
        enddo
      enddo
      call dfftw_execute_dft(FFTW_plan,MatIn,MatOut)
      do kc=1,im3
        do ic=1,im1
          cMatZ(ic,jc,kc)=MatOut(kc,ic)
        enddo
      enddo
    enddo
    call fftw_destroy_plan(FFTW_plan)
    deallocate(MatIn,MatOut)    
  end subroutine ifft3d_z

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! cmplx_transpose_z2_to_x1
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90  
  subroutine cmplx_transpose_z2_to_x1(cMatZ,cMatX,decomp_HIT); implicit none
    type(decomp_info),intent(in)::decomp_HIT
    complex(RK),dimension(decomp_HIT%z2sz(1),decomp_HIT%z2sz(2),decomp_HIT%z2sz(3)),intent(in) ::cMatZ
    complex(RK),dimension(decomp_HIT%x1sz(1),decomp_HIT%x1sz(2),decomp_HIT%x1sz(3)),intent(out)::cMatX
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(decomp_HIT%z2sz(1),decomp_HIT%z2sz(2),decomp_HIT%z2sz(3))::rMatZ
    real(RK),dimension(decomp_HIT%x1sz(1),decomp_HIT%x1sz(2),decomp_HIT%x1sz(3))::rMatX
        
    rMatZ=real(cMatZ,RK)
    call transpose_z2_to_x1(rMatZ,rMatX,decomp_HIT)
    do kc=1,decomp_HIT%x1sz(3)
      do jc=1,decomp_HIT%x1sz(2)
        do ic=1,decomp_HIT%x1sz(1)
          cMatX(ic,jc,kc)=cmplx(rMatX(ic,jc,kc),0.0_RK,RK)
        enddo
      enddo
    enddo
    
    rMatZ=aimag(cMatZ)
    call transpose_z2_to_x1(rMatZ,rMatX,decomp_HIT)
    do kc=1,decomp_HIT%x1sz(3)
      do jc=1,decomp_HIT%x1sz(2)
        do ic=1,decomp_HIT%x1sz(1)
          cMatX(ic,jc,kc)=cmplx(0.0_RK,rMatX(ic,jc,kc),RK)+cMatX(ic,jc,kc)
        enddo
      enddo
    enddo    
  end subroutine cmplx_transpose_z2_to_x1

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ifft3d_x
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ifft3d_x(cMatX,ArrX,im1,im2,im3); implicit none
    include "fftw3.f03"
    integer,intent(in)::im1,im2,im3
    real(RK),dimension(im1,im2,im3),intent(out)::ArrX
    complex(RK),dimension(im1/2+1,im2,im3),intent(in)::cMatX
    
    ! locals
    type(C_PTR)::FFTW_plan
    integer::ic,jc,kc,imc,it
    type(fftw_iodim)::iodim(1),iodim_howmany(1) 
    real(RK),dimension(:,:),allocatable::MatIn,MatOut
    
    imc=im1/2+1
    allocate(MatIn (im1,im2))
    allocate(MatOut(im1,im2))
    iodim(1)%n  = im1
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n  = im2
    iodim_howmany(1)%is = im1
    iodim_howmany(1)%os = im1
    FFTW_plan= fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,MatIn,MatOut,[FFTW_HC2R],FFTW_ESTIMATE)
    do kc=1,im3
      do jc=1,im2
        do ic=1,imc
          MatIn(ic,jc)=real(cMatX(ic,jc,kc),RK)
        enddo
        do ic=imc+1,im1
          it=2*imc-ic
          MatIn(ic,jc)=aimag(cMatX(it,jc,kc))
        enddo
      enddo
      call dfftw_execute_r2r(FFTW_plan,MatIn,MatOut)
      do jc=1,im2
        do ic=1,im1
          ArrX(ic,jc,kc)=MatOut(ic,jc)
        enddo
      enddo
    enddo
    call fftw_destroy_plan(FFTW_plan)
    deallocate(MatIn,MatOut)  
  end subroutine ifft3d_x

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clcStat_HIT
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine clcStat_HIT(ux,uy,uz,pr,ArrTemp); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pr
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2), y1start(3):y1end(3)),intent(out)::ArrTemp
    
    ! locals
    integer::iUnit,ierr,ic,jc,kc
    character(len=:),allocatable::filename
    real(Rk)::sumStat(2),sumr(2),tau,epsi,ug,Lambda,ReLam
    
    filename = ' '
    !================
    call CalcDissipationAndTKE(ux,uy,uz,sumr)
    call MPI_REDUCE(sumr,sumStat,2,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(nrank==0) then
      epsi=xnu*sumStat(1)/real(nxc*nyc*nzc,RK)
      ug=sqrt(sumStat(2)/real(nxc*nyc*nzc,RK)/3.0_RK)
      Lambda=sqrt(15.0_RK*xnu/epsi)*ug
      ReLam=Lambda*ug/xnu
      tau=ug*ug/epsi
      filename = strip(ResultsDir_) // "HITStat" // int2str(ilast,10)
      open(newunit=iUnit, file=filename, status='old',position='append',form='formatted',IOSTAT=ierr )
      IF(ierr/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_HIT","Cannot open file: "// filename)
      ELSE
        write(iUnit,'(6ES24.15)')SimTime,tau,ug*ug,epsi,ReLam
      ENDIF
      close(iUnit,IOSTAT=ierr)
    endif
    
    !======================
    do kc=y1start(3),y1end(3); do jc=y1start(2),y1end(2); do ic=y1start(1),y1end(1)
    ArrTemp(ic,jc,kc) =ux(ic, jc, kc); enddo; enddo; enddo
    call clcSpectra3D(ArrTemp,1)

    do kc=y1start(3),y1end(3); do jc=y1start(2),y1end(2); do ic=y1start(1),y1end(1)
    ArrTemp(ic,jc,kc) =uy(ic, jc, kc); enddo; enddo; enddo
    call clcSpectra3D(ArrTemp,2)
    
    do kc=y1start(3),y1end(3); do jc=y1start(2),y1end(2); do ic=y1start(1),y1end(1)
    ArrTemp(ic,jc,kc) =uz(ic, jc, kc); enddo; enddo; enddo
    call clcSpectra3D(ArrTemp,3)
    
    do kc=y1start(3),y1end(3); do jc=y1start(2),y1end(2); do ic=y1start(1),y1end(1)
    ArrTemp(ic,jc,kc) =pr(ic, jc, kc); enddo; enddo; enddo
    call clcSpectra3D(ArrTemp,4)
        
    nfstime= nfstime + 1
    if(mod(itime,SaveStat)/=0) return
    block
    real(RK)::infstime
    character(len=128)::FormatStr
    real(RK),allocatable,dimension(:,:)::EnergySpec3DR
    allocate(EnergySpec3DR(NEnergySpec3D, nLenSpec3D),Stat=ierr)
    call MPI_REDUCE(EnergySpec3D,EnergySpec3DR, nLenSpec3D*NEnergySpec3D,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(nrank==0) then
      filename = strip(ResultsDir_) // 'EnergySpectra3D' // int2str(itime,10)
      open(newunit=iUnit,file=filename,status='replace',form='formatted',IOSTAT=ierr)
      IF(ierr/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_HIT","Cannot open file: "//filename)
      ELSE
        write(iUnit,'(a,I7,a,I7,a,I7)')'  The time step range for this fluid statistics is ', &
                                    itime-(nfstime-1)*ivstats, ':', ivstats, ':', itime
        write(iUnit,'(A)')'  '
        write(iUnit,'(A)')' k, E_u, E_v, E_w, E_p'        
        infstime = 1.0_RK/real(nfstime,RK)
        write(FormatStr,'(A,I3,A)')'(I8,',NEnergySpec3D,'ES24.15)'
        do ic=1,nLenSpec3D
          write(iUnit,FormatStr) ic-1, EnergySpec3DR(1:NEnergySpec3D, ic)*infstime
        enddo
      ENDIF
      close(iUnit,IOSTAT=ierr)
    endif
    end block
    
    nfstime=0; EnergySpec3D=0.0_RK
  end subroutine clcStat_HIT

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CalcDissipationAndTKE
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CalcDissipationAndTKE(ux,uy,uz,sumStat); implicit none
    real(RK),intent(out)::sumStat(2)
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz

    ! locals
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::st1,st2,st3,st4,st5,st6,caj,cac
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    
    ! epsilon_ij= (2*xnu*<S_ij*S_ij>), where S_ij = 0.5*(dui/dxj + duj/dxi)
    sumStat=0.0_RK
    DO kc=y1start(3),y1end(3)
      kp=kc+1; km=kc-1
      do jc=y1start(2),y1end(2)
        jp=jc+1; jm=jc-1
        caj= rdyp(jc)
        cac= rdyc(jc)
        do ic=y1start(1),y1end(1)
          ip=ic+1; im=ic-1
          dudx=(ux(ip,jc,kc)-ux(ic,jc,kc))*rdx
          dudy=(ux(ic,jc,kc)-ux(ic,jm,kc))*cac
          dudz=(ux(ic,jc,kc)-ux(ic,jc,km))*rdz
          dvdx=(uy(ic,jc,kc)-uy(im,jc,kc))*rdx
          dvdy=(uy(ic,jp,kc)-uy(ic,jc,kc))*caj
          dvdz=(uy(ic,jc,kc)-uy(ic,jc,km))*rdz
          dwdx=(uz(ic,jc,kc)-uz(im,jc,kc))*rdx
          dwdy=(uz(ic,jc,kc)-uz(ic,jm,kc))*cac
          dwdz=(uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
          st1=dudx
          st2=dvdy
          st3=dwdz
          st4=0.5_RK*(dudy+dvdx)
          st5=0.5_RK*(dudz+dwdx)
          st6=0.5_RK*(dvdz+dwdy)
          sumStat(1)=sumStat(1)+2.0_RK*(st1*st1+ st2*st2+ st3*st3)+ 4.0_RK*(st4*st4 + st5*st5 + st6*st6)
          sumStat(2)=sumStat(2)+ux(ic,jc,kc)*ux(ic,jc,kc)+uy(ic,jc,kc)*uy(ic,jc,kc)+uz(ic,jc,kc)*uz(ic,jc,kc)
        enddo
      enddo
    ENDDO
  end subroutine CalcDissipationAndTKE

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! add_FlowType_Forcing_HIT
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90 
  subroutine add_FlowType_Forcing_HIT(RhsX,RhsY,RhsZ); implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX,RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ
    
    ! locals
    integer::ic,jc,kc
    real(RK)::xcTmp, ycTmp, zcTmp, A_, B_, C_, twoPiDL, F0
    
    F0 = F0_intensity*pmAlpha
    A_ = 1.0_RK *F0
    B_ = 1.0_RK *F0
    C_ = 1.0_RK *F0
    twoPiDL = 2.0_RK*Pi/xlx
    DO kc=y1start(3),y1end(3)
      zcTmp = (real(kc, RK)-0.5_RK)*dz
      do jc=y1start(2),y1end(2)
        ycTmp = yc(jc)
        do ic=y1start(1),y1end(1)
          xcTmp = (real(ic, RK)-0.5_RK)*dx
          RhsX(ic,jc,kc) = RhsX(ic,jc,kc) +A_ *sin(twoPiDL *zcTmp) +C_ *cos(twoPiDL *ycTmp)
          RhsY(ic,jc,kc) = RhsY(ic,jc,kc) +B_ *sin(twoPiDL *xcTmp) +A_ *cos(twoPiDL *zcTmp)
          RhsZ(ic,jc,kc) = RhsZ(ic,jc,kc) +C_ *sin(twoPiDL *ycTmp) +B_ *cos(twoPiDL *xcTmp)
        enddo
      enddo
    ENDDO 
  end subroutine add_FlowType_Forcing_HIT

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clcSpectra3D
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine clcSpectra3D(ArrIN,m1); implicit none
    integer,intent(in)::m1
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in)::ArrIn
    
    ! locals
    integer::ic,jc,kc,ix,iy,iz
    real(RK)::Ratio1, Ratio2, Ratio3
    real(RK),allocatable,dimension(:)::RatioXYZ
    real(RK),allocatable,dimension(:,:,:)::arr1, arr2
    
    allocate(arr1(nxc,x1size(2),x1size(3)))
    call transpose_y1_to_x1(ArrIN, arr1)
    call my_execute_FFTW_r2r_x(fft_plan_x,x1size(1),x1size(2),x1size(3),arr1)
    
    allocate(arr2(z2size(1),z2size(2), nzc))
    call transpose_x1_to_z2(arr1, arr2)
    call my_execute_FFTW_r2r_z(fft_plan_z,z2size(1),z2size(2),z2size(3),arr2)
    deallocate(arr1)
    
    allocate(arr1(y2start(1):y2end(1), nyc, y2start(3):y2end(3)))
    call transpose_z2_to_y2(arr2, arr1)
    call my_execute_FFTW_r2r_y(fft_plan_y,y2size(1),y2size(2),y2size(3),arr1)
    deallocate(arr2)
    
    allocate(RatioXYZ(nxc))
    RatioXYZ =  2.0_RK/real(nxc, RK)/real(nxc, RK)
    ic=1;       RatioXYZ(ic) =0.5_RK*RatioXYZ(ic)
    ic=nxc/2+1; RatioXYZ(ic) =0.5_RK*RatioXYZ(ic)
    
    do kc=y2start(3),y2end(3)
      Ratio3 = RatioXYZ(kc)
      iz = IndXYZ(kc)*IndXYZ(kc)
      do jc=y2start(2),y2end(2)
        Ratio2 = RatioXYZ(jc)*Ratio3
        iy = IndXYZ(jc)*IndXYZ(jc) +iz
        do ic=y2start(1),y2end(1)
          Ratio1 = RatioXYZ(ic)*Ratio2
          ix = IndXYZ(ic)*IndXYZ(ic) +iy
          ix = nint(sqrt(1.0_RK*ix)) +1
          EnergySpec3D(m1, ix) =EnergySpec3D(m1, ix) +Ratio1 *arr1(ic,jc,kc)*arr1(ic,jc,kc)
        enddo
      enddo
    enddo
    deallocate(arr1)
  end subroutine clcSpectra3D
      
#undef NEnergySpec2D    
end module f2_FlowType_HIT
