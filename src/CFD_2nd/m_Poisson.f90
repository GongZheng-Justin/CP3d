module m_Poisson
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_decomp2d
  use m_Parameters
  use iso_c_binding
  use m_MeshAndMetries
  use m_Variables,only: mb1
  use m_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal
  implicit none
  private
  include "fftw3.f03"
  
  real(RK)::normfft
  procedure(),pointer::clcPPE=>null()
  type(decomp_info),allocatable::decomp_PPE
  real(RK),allocatable,dimension(:)::WaveNumX,WaveNumZ
  type(C_PTR)::fwd_plan_x,bwd_plan_x,fwd_plan_z,bwd_plan_z
  real(RK),allocatable,dimension(:,:,:)::a_reduce,c_reduce 
  
  public:: InitPoissonSolver,clcPPE,Destory_Poisson_FFT_Plan
contains
    
  !******************************************************************
  ! InitPoissonSolver
  !******************************************************************     
  subroutine InitPoissonSolver(prsrc,prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc,prphiHalo
        
    ! locals
    real(RK)::WaveCoe,wa1,wa3,best_time,t1,t2,norm
    integer::i,k,iErr01,iErr02,iChoice,IsReduce,IsReduceR

    normfft = one
    
    ! Modified wave number in x-dir
    allocate(WaveNumX(nxc), Stat =iErr01)
    IF(BcOption(xm_dir)==BC_PERIOD ) THEN    
      WaveCoe=twopi
      norm = one
    ELSE
      WaveCoe=pi
      norm = two
    ENDIF
    do i=1,nxc
      wa1= WaveCoe*real(i-1,RK)/real(nxc,RK)
      WaveNumX(i)=two*rdx2*(cos(wa1)-one)
    enddo
    normfft = normfft*norm*real(nxc,RK)

    ! Modified wave number in z-dir    
    allocate(WaveNumZ(nzc), Stat =iErr02)
    IF(BcOption(zm_dir)==BC_PERIOD ) THEN    
      WaveCoe=twopi
      norm = one
    ELSE
      WaveCoe=pi
      norm = two
    ENDIF   
    do k=1,nzc
      wa3= WaveCoe*real(k-1,RK)/real(nzc,RK)
      WaveNumZ(k)=two*rdz2*(cos(wa3)-one)
    enddo
    normfft = normfft*norm*real(nzc,RK); normfft=one/normfft
        
    best_time=huge(t1)
    if(nrank==0) call MainLog%OutInfo("Auto-tuning mode for Poisson Solver......",1)
    
    IF(BcOption(ym_dir)==BC_PERIOD ) THEN
      ! Choice-1
      call Create_Poisson_FFT_Plan(x1size,z2size)
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      t1=MPI_WTIME()
      do k=1,15
        call clcPPE_x1_periodic(prsrc,prphiHalo)
      enddo
      t2=MPI_WTIME()-t1
      call MPI_ALLREDUCE(t2,t1,1,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0)call MainLog%OutInfo("Choice-1, time="//num2str(t1),2)
      if(best_time>t1) then
        best_time=t1
        iChoice=1
      endif
      call Destory_Poisson_FFT_Plan()
   
      ! Choice-2    
      call Create_Poisson_FFT_Plan(x2size,z1size)
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      t1=MPI_WTIME()
      do k=1,15
        call clcPPE_z1_periodic(prsrc,prphiHalo)
      enddo
      t2=MPI_WTIME()-t1
      call MPI_ALLREDUCE(t2,t1,1,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0)call MainLog%OutInfo("Choice-2, time="//num2str(t1),2)
      if(best_time>t1) then
        best_time=t1
        iChoice=2
      endif
      call Destory_Poisson_FFT_Plan()    
    ELSE
      ! Choice-1
      call Create_Poisson_FFT_Plan(x1size,z2size)
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      t1=MPI_WTIME()
      do k=1,15
        call clcPPE_x1(prsrc,prphiHalo)
      enddo
      t2=MPI_WTIME()-t1
      call MPI_ALLREDUCE(t2,t1,1,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0)call MainLog%OutInfo("Choice-1, time="//num2str(t1),2)
      if(best_time>t1) then
        best_time=t1
        iChoice=1
      endif
      call Destory_Poisson_FFT_Plan()
   
      ! Choice-2    
      call Create_Poisson_FFT_Plan(x2size,z1size)
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      t1=MPI_WTIME()
      do k=1,15
        call clcPPE_z1(prsrc,prphiHalo)
      enddo
      t2=MPI_WTIME()-t1
      call MPI_ALLREDUCE(t2,t1,1,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0)call MainLog%OutInfo("Choice-2, time="//num2str(t1),2)
      if(best_time>t1) then
        best_time=t1
        iChoice=2
      endif
      call Destory_Poisson_FFT_Plan()

      ! Choice-3
      call Create_Poisson_FFT_Plan(x1size,z2size)
      IsReduce=1
      if(z2size(2)<4)IsReduce=0
      call MPI_ALLREDUCE(IsReduce,IsReduceR,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD,iErr01)
      if(IsReduceR==1) then
        allocate(decomp_PPE)
        allocate(a_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
        allocate(c_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
        call Initialize_ReduceMatrix(z2start,z2end,'z')
        call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
        t1=MPI_WTIME()
        do k=1,15
          call clcPPE_x1_reduce(prsrc,prphiHalo)
        enddo
        t2=MPI_WTIME()-t1
        call MPI_ALLREDUCE(t2,t1,1,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
        if(nrank==0)call MainLog%OutInfo("Choice-3, time="//num2str(t1),2)
        if(best_time>t1) then
          best_time=t1
          iChoice=3
        endif
        call decomp_info_finalize(decomp_PPE)
        deallocate(decomp_PPE,a_reduce,c_reduce)
      else
        if(nrank==0)call MainLog%OutInfo("Choice-3, z2size<4, ignore",2)
      endif
      call Destory_Poisson_FFT_Plan()
        
      ! Choice-4 
      call Create_Poisson_FFT_Plan(x2size,z1size)
      IsReduce=1
      if(x2size(2)<4)IsReduce=0
      call MPI_ALLREDUCE(IsReduce,IsReduceR,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD,iErr01)
      if(IsReduceR==1) then   
      allocate(decomp_PPE)
        allocate(a_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
        allocate(c_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
        call Initialize_ReduceMatrix(x2start,x2end,'x')
        call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
        t1=MPI_WTIME()
        do k=1,15
          call clcPPE_z1_reduce(prsrc,prphiHalo)
        enddo
        t2=MPI_WTIME()-t1
        call MPI_ALLREDUCE(t2,t1,1,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
        if(nrank==0)call MainLog%OutInfo("Choice-4, time="//num2str(t1),2)
        if(best_time>t1) then
          best_time=t1
          iChoice=4
        endif
        call decomp_info_finalize(decomp_PPE)
        deallocate(decomp_PPE,a_reduce,c_reduce)
      else
        if(nrank==0)call MainLog%OutInfo("Choice-4, x2size<4, ignore",2)
      endif
      call Destory_Poisson_FFT_Plan()
    ENDIF

#ifdef FFTW_1D
    iChoice=2 ! Caution, For testing here !
#endif
    if(nrank==0) then
      call MainLog%OutInfo("The best Poisson Solver choice is probably Choice-"//num2str(iChoice),2)
      call MainLog%OutInfo("Corresponding Global Data Transpose is:",2)
    endif
    if(iChoice==1) then
      IF(BcOption(ym_dir)==BC_PERIOD ) THEN
        clcPPE => clcPPE_x1_periodic      
      ELSE
        clcPPE => clcPPE_x1
      ENDIF
      call Create_Poisson_FFT_Plan(x1size,z2size)
      if(nrank==0) call MainLog%OutInfo("y1 -> x1 -> z2 -> y2 -> z2 -> x1 -> y1",3)
    elseif(iChoice==2) then
      IF(BcOption(ym_dir)==BC_PERIOD ) THEN
        clcPPE => clcPPE_z1_periodic      
      ELSE
        clcPPE => clcPPE_z1
      ENDIF
      call Create_Poisson_FFT_Plan(x2size,z1size)
      if(nrank==0) call MainLog%OutInfo("y1 -> z1 -> x2 -> y2 -> x2 -> z1 -> y1",3)
    elseif(iChoice==3) then
      clcPPE => clcPPE_x1_reduce 
      call Create_Poisson_FFT_Plan(x1size,z2size)
      if(nrank==0) call MainLog%OutInfo("y1 -> x1 -> z2 -> x1 -> y1",3)
      allocate(decomp_PPE)
      allocate(a_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
      allocate(c_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
      call Initialize_ReduceMatrix(z2start,z2end,'z')
    elseif(iChoice==4) then
      clcPPE => clcPPE_z1_reduce 
      call Create_Poisson_FFT_Plan(x2size,z1size)
      if(nrank==0) call MainLog%OutInfo("y1 -> z1 -> x2 -> z1 -> y1",3)
      allocate(decomp_PPE)
      allocate(a_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
      allocate(c_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
      call Initialize_ReduceMatrix(x2start,x2end,'x')
    endif
    prsrc=zero; prphiHalo=zero
  end subroutine InitPoissonSolver

  !******************************************************************
  ! Create_Poisson_FFT_Plan
  !******************************************************************
  subroutine Create_Poisson_FFT_Plan(xsizeIn,zsizeIn)
    implicit none
    integer,dimension(3),intent(in)::xsizeIn,zsizeIn
   
    ! locals
    integer::plan_type
#ifdef FFTW_1D
    integer(C_FFTW_R2R_KIND)::kind_fwd,kind_bwd
    real(RK),dimension(:),allocatable::Vec1,Vec2
#else
    type(fftw_iodim)::iodim(1),iodim_howmany(1)
#ifdef OverWriteFFT
    real(RK),dimension(:,:,:),allocatable::Arr
#else
    real(RK),dimension(:,:,:),allocatable::Arr1,Arr2
#endif
    integer(C_FFTW_R2R_KIND)::kind_fwd(1),kind_bwd(1)
#endif
    
    if(FFTW_plan_type == 1) then
      plan_type=FFTW_MEASURE
    else
      plan_type=FFTW_ESTIMATE
    endif
        
    ! FFT in x
    IF(BcOption(xm_dir)==BC_PERIOD ) THEN
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
    ELSE
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01    
    ENDIF
#ifdef FFTW_1D
    allocate(Vec1(xsizeIn(1)),Vec2(xsizeIn(1)))
    fwd_plan_x= fftw_plan_r2r_1d(xsizeIn(1),Vec1,Vec2,kind_fwd,plan_type)
    bwd_plan_x= fftw_plan_r2r_1d(xsizeIn(1),Vec1,Vec2,kind_bwd,plan_type)
    deallocate(Vec1,Vec2)
#else
    iodim(1)%n  = xsizeIn(1)
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n  = xsizeIn(2)*xsizeIn(3)
    iodim_howmany(1)%is = xsizeIn(1)
    iodim_howmany(1)%os = xsizeIn(1)
#ifdef OverWriteFFT
    allocate(Arr(xsizeIn(1),xsizeIn(2),xsizeIn(3)))
    fwd_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr, Arr, kind_fwd,plan_type)
    bwd_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr, Arr, kind_bwd,plan_type)
    deallocate(Arr)
#else
    allocate(Arr1(xsizeIn(1),xsizeIn(2),xsizeIn(3)))
    allocate(Arr2(xsizeIn(1),xsizeIn(2),xsizeIn(3)))
    fwd_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,kind_fwd,plan_type)
    bwd_plan_x = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,kind_bwd,plan_type)
    deallocate(Arr1,Arr2)
#endif
#endif
   
    ! FFT in z
    IF(BcOption(zm_dir)==BC_PERIOD ) THEN
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
    ELSE
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01    
    ENDIF
#ifdef FFTW_1D
    allocate(Vec1(zsizeIn(3)),Vec2(zsizeIn(3)))
    fwd_plan_z= fftw_plan_r2r_1d(zsizeIn(3),Vec1,Vec2,kind_fwd,plan_type)
    bwd_plan_z= fftw_plan_r2r_1d(zsizeIn(3),Vec1,Vec2,kind_bwd,plan_type)
    deallocate(Vec1,Vec2)
#else
    iodim(1)%n  = zsizeIn(3)
    iodim(1)%is = zsizeIn(1)*zsizeIn(2)
    iodim(1)%os = zsizeIn(1)*zsizeIn(2)
    iodim_howmany(1)%n  = zsizeIn(1)*zsizeIn(2)
    iodim_howmany(1)%is = 1
    iodim_howmany(1)%os = 1
#ifdef OverWriteFFT
    allocate(Arr(zsizeIn(1),zsizeIn(2),zsizeIn(3)))
    fwd_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr, Arr, kind_fwd,plan_type)
    bwd_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr, Arr, kind_bwd,plan_type)
    deallocate(Arr)
#else
    allocate(Arr1(zsizeIn(1),zsizeIn(2),zsizeIn(3)))
    allocate(Arr2(zsizeIn(1),zsizeIn(2),zsizeIn(3)))    
    fwd_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,kind_fwd,plan_type)
    bwd_plan_z = fftw_plan_guru_r2r(1,iodim,1,iodim_howmany,Arr1,Arr2,kind_bwd,plan_type)
    deallocate(Arr1,Arr2)
#endif
#endif
  end subroutine Create_Poisson_FFT_Plan

  !******************************************************************
  ! Initialize_ReduceMatrix
  !******************************************************************
  subroutine Initialize_ReduceMatrix(startIn,endIn,x_or_z)
    implicit none
    character(len=1),intent(in)::x_or_z
    integer,dimension(3),intent(in)::startIn,endIn

    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc
    logical,dimension(6)::initializeIn
    real(RK),dimension(:,:,:),allocatable::ArrM,ArrC,ArrP
    
    initializeIn=.false.
    if(x_or_z=='x') then
      initializeIn(4)=.true.
      call decomp_info_init(nxc,2*p_col,nzc,decomp_PPE,initialize=initializeIn)
    elseif(x_or_z=='z') then
      initializeIn(5)=.true.
      call decomp_info_init(nxc,2*p_row,nzc,decomp_PPE,initialize=initializeIn)
    endif
    allocate(ArrM(startIn(1):endIn(1),startIn(2):endIn(2),startIn(3):endIn(3)))
    allocate(ArrC(startIn(1):endIn(1),startIn(2):endIn(2),startIn(3):endIn(3)))
    allocate(ArrP(startIn(1):endIn(1),startIn(2):endIn(2),startIn(3):endIn(3)))
    do kc=startIn(3),endIn(3)
      do jc=startIn(2),endIn(2)
        do ic=startIn(1),endIn(1)
          ArrM(ic,jc,kc)=am2Pr(jc)
          ArrC(ic,jc,kc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          ArrP(ic,jc,kc)=ap2Pr(jc)
        enddo
      enddo
    enddo
    IF(nrank==0) THEN
      ArrM(1,1,1)=zero
      ArrC(1,1,1)=one
      ArrP(1,1,1)=zero
    ENDIF
    do kc=startIn(3),endIn(3)
      do ic=startIn(1),endIn(1)
        jc=startIn(2)
        ArrM(ic,jc,kc)=ArrM(ic,jc,kc)/ArrC(ic,jc,kc)
        ArrP(ic,jc,kc)=ArrP(ic,jc,kc)/ArrC(ic,jc,kc)
        jc=startIn(2)+1
        ArrM(ic,jc,kc)=ArrM(ic,jc,kc)/ArrC(ic,jc,kc)
        ArrP(ic,jc,kc)=ArrP(ic,jc,kc)/ArrC(ic,jc,kc)
        do jc=startIn(2)+2,endIn(2)
          rTemp=one/(ArrC(ic,jc,kc)-ArrM(ic,jc,kc)*ArrP(ic,jc-1,kc))
          ArrP(ic,jc,kc)= rTemp*ArrP(ic,jc,kc)
          ArrM(ic,jc,kc)=-rTemp*ArrM(ic,jc,kc)*ArrM(ic,jc-1,kc)
        enddo
        do jc=endIn(2)-2,startIn(2)+1,-1
          ArrM(ic,jc,kc)= ArrM(ic,jc,kc)-ArrP(ic,jc,kc)*ArrM(ic,jc+1,kc)
          ArrP(ic,jc,kc)=-ArrP(ic,jc,kc)*ArrP(ic,jc+1,kc)
        enddo
        jc=startIn(2)
        rTemp=one/(one-ArrM(ic,jc+1,kc)*ArrP(ic,jc,kc))
        ArrM(ic,jc,kc)= rTemp*ArrM(ic,jc,kc)
        ArrP(ic,jc,kc)=-rTemp*ArrP(ic,jc,kc)*ArrP(ic,jc+1,kc)        
      enddo
    enddo    
    deallocate(ArrC)
    allocate(ArrC(startIn(1):endIn(1),2,startIn(3):endIn(3)))
    
    do kc=startIn(3),endIn(3)
      do ic=startIn(1),endIn(1)
        ArrC(ic,1,kc)=ArrM(ic,startIn(2),kc)
        ArrC(ic,2,kc)=ArrM(ic,endIn(2),  kc)
      enddo
    enddo
    if(x_or_z=='x') then
      call transpose_x2_to_y2(ArrC,a_reduce,decomp_PPE)
    elseif(x_or_z=='z') then
      call transpose_z2_to_y2(ArrC,a_reduce,decomp_PPE)
    endif
    do kc=startIn(3),endIn(3)
      do ic=startIn(1),endIn(1)
        ArrC(ic,1,kc)=ArrP(ic,startIn(2),kc)
        ArrC(ic,2,kc)=ArrP(ic,endIn(2),  kc)
      enddo
    enddo
    if(x_or_z=='x') then
      call transpose_x2_to_y2(ArrC,c_reduce,decomp_PPE)
    elseif(x_or_z=='z') then
      call transpose_z2_to_y2(ArrC,c_reduce,decomp_PPE)
    endif  
    deallocate(ArrM,ArrC,ArrP)         
  end subroutine Initialize_ReduceMatrix

  !******************************************************************
  ! Destory_Poisson_FFT_Plan
  !******************************************************************
  subroutine Destory_Poisson_FFT_Plan()
    implicit none
    call dfftw_destroy_plan(fwd_plan_x,bwd_plan_x) 
    call dfftw_destroy_plan(fwd_plan_z,bwd_plan_z)
    if(allocated(decomp_PPE)) then
      call decomp_info_finalize(decomp_PPE)
      deallocate(decomp_PPE)
    endif
  end subroutine Destory_Poisson_FFT_Plan

  !******************************************************************
  ! clcPPE_x1_periodic
  !******************************************************************
  subroutine clcPPE_x1_periodic(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    integer::ic,jc,kc
#ifdef FFTW_1D
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx
    real(RK),dimension(z2size(1),z2size(2),z2size(3))::arrz
#else
#ifdef OverWriteFFT
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx
    real(RK),dimension(z2size(1),z2size(2),z2size(3))::arrz
#else
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1,arrx2
    real(RK),dimension(z2size(1),z2size(2),z2size(3))::arrz1,arrz2
#endif
#endif
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2),y2start(3):y2end(3))::prTemp
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2))::tridmj,tridcj,tridpj,tridfj

#ifdef FFTW_1D
    call transpose_y1_to_x1(prsrc,arrx)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_z2(arrx,arrz)
    call my_execute_FFTW_r2r_z(fwd_plan_z,z2size(1),z2size(2),z2size(3),arrz)
    call transpose_z2_to_y2(arrz,prTemp)    
#else
#ifdef OverWriteFFT
    call transpose_y1_to_x1(prsrc,arrx)
    call dfftw_execute_r2r(fwd_plan_x,arrx,arrx)
    call transpose_x1_to_z2(arrx,arrz)
    call dfftw_execute_r2r(fwd_plan_z,arrz,arrz)
    call transpose_z2_to_y2(arrz,prTemp)
#else
    call transpose_y1_to_x1(prsrc,arrx1)
    call dfftw_execute_r2r(fwd_plan_x,arrx1,arrx2)
    call transpose_x1_to_z2(arrx2,arrz1)
    call dfftw_execute_r2r(fwd_plan_z,arrz1,arrz2)
    call transpose_z2_to_y2(arrz2,prTemp)
#endif
#endif

    do jc=y2start(2),y2end(2)
      do ic=y2start(1),y2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridpj(ic,jc)=ap2Pr(jc)
      enddo
    enddo
    IF(nrank==0) THEN
      do kc=2,y2end(3)
        do jc=1,nyc
          do ic=1,y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo

      ! for nrank=0, y1start(1)=y1start(3)=1
      ic=1; jc=1; kc=1;
      tridpj(ic,jc)=zero
      tridcj(ic,jc)=one
      tridmj(ic,jc)=zero
      tridfj(ic,jc)=zero
      do ic=2,y2end(1)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridfj(ic,jc)=prTemp(ic,jc,kc)
      enddo
      do jc=2,nyc
        do ic=1,y2end(1)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridfj(ic,jc)=prTemp(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,tridfj,y2size(1),nyc)
      do jc=1,nyc
        do ic=1,y2end(1)
          prTemp(ic,jc,kc)= tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=y2start(3),y2end(3)
        do jc=y2start(2),y2end(2)
          do ic=y2start(1),y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo       
    ENDIF
#ifdef FFTW_1D
    call transpose_y2_to_z2(prTemp,arrz)
    call my_execute_FFTW_r2r_z(bwd_plan_z,z2size(1),z2size(2),z2size(3),arrz)
    call transpose_z2_to_x1(arrz,arrx)
    call my_execute_FFTW_r2r_x(bwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_y1(arrx,prsrc)
#else
#ifdef OverWriteFFT
    call transpose_y2_to_z2(prTemp,arrz)
    call dfftw_execute_r2r(bwd_plan_z,arrz, arrz)
    call transpose_z2_to_x1(arrz,arrx)
    call dfftw_execute_r2r(bwd_plan_x,arrx, arrx)
    call transpose_x1_to_y1(arrx,prsrc)
#else
    call transpose_y2_to_z2(prTemp,arrz2)
    call dfftw_execute_r2r(bwd_plan_z,arrz2, arrz1)
    call transpose_z2_to_x1(arrz1,arrx2)
    call dfftw_execute_r2r(bwd_plan_x,arrx2, arrx1)
    call transpose_x1_to_y1(arrx1,prsrc)
#endif
#endif
    
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_x1_periodic
        
  !******************************************************************
  ! clcPPE_x1
  !******************************************************************
  subroutine clcPPE_x1(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    integer::ic,jc,kc
#ifdef FFTW_1D
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx
    real(RK),dimension(z2size(1),z2size(2),z2size(3))::arrz
#else
#ifdef OverWriteFFT
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx
    real(RK),dimension(z2size(1),z2size(2),z2size(3))::arrz
#else
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1,arrx2
    real(RK),dimension(z2size(1),z2size(2),z2size(3))::arrz1,arrz2
#endif
#endif
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2),y2start(3):y2end(3))::prTemp
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2))::tridmj,tridcj,tridpj,tridfj
    
#ifdef FFTW_1D
    call transpose_y1_to_x1(prsrc,arrx)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_z2(arrx,arrz)
    call my_execute_FFTW_r2r_z(fwd_plan_z,z2size(1),z2size(2),z2size(3),arrz)
    call transpose_z2_to_y2(arrz,prTemp)    
#else
#ifdef OverWriteFFT
    call transpose_y1_to_x1(prsrc,arrx)
    call dfftw_execute_r2r(fwd_plan_x,arrx,arrx)
    call transpose_x1_to_z2(arrx,arrz)
    call dfftw_execute_r2r(fwd_plan_z,arrz,arrz)
    call transpose_z2_to_y2(arrz,prTemp)
#else
    call transpose_y1_to_x1(prsrc,arrx1)
    call dfftw_execute_r2r(fwd_plan_x,arrx1,arrx2)
    call transpose_x1_to_z2(arrx2,arrz1)
    call dfftw_execute_r2r(fwd_plan_z,arrz1,arrz2)
    call transpose_z2_to_y2(arrz2,prTemp)
#endif
#endif

    do jc=y2start(2),y2end(2)
      do ic=y2start(1),y2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridpj(ic,jc)=ap2Pr(jc)
      enddo
    enddo
    IF(nrank==0) THEN
      do kc=2,y2end(3)
        do jc=1,nyc
          do ic=1,y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InverseTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo

      ! for nrank=0, y1start(1)=y1start(3)=1
      ic=1; jc=1; kc=1;
      tridpj(ic,jc)=zero
      tridcj(ic,jc)=one
      tridmj(ic,jc)=zero
      tridfj(ic,jc)=zero
      do ic=2,y2end(1)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridfj(ic,jc)=prTemp(ic,jc,kc)
      enddo
      do jc=2,nyc
        do ic=1,y2end(1)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridfj(ic,jc)=prTemp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y2size(1),nyc)
      do jc=1,nyc
        do ic=1,y2end(1)
          prTemp(ic,jc,kc)= tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=y2start(3),y2end(3)
        do jc=y2start(2),y2end(2)
          do ic=y2start(1),y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InverseTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo       
    ENDIF
#ifdef FFTW_1D
    call transpose_y2_to_z2(prTemp,arrz)
    call my_execute_FFTW_r2r_z(bwd_plan_z,z2size(1),z2size(2),z2size(3),arrz)
    call transpose_z2_to_x1(arrz,arrx)
    call my_execute_FFTW_r2r_x(bwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_y1(arrx,prsrc)
#else
#ifdef OverWriteFFT
    call transpose_y2_to_z2(prTemp,arrz)
    call dfftw_execute_r2r(bwd_plan_z,arrz, arrz)
    call transpose_z2_to_x1(arrz,arrx)
    call dfftw_execute_r2r(bwd_plan_x,arrx, arrx)
    call transpose_x1_to_y1(arrx,prsrc)
#else
    call transpose_y2_to_z2(prTemp,arrz2)
    call dfftw_execute_r2r(bwd_plan_z,arrz2, arrz1)
    call transpose_z2_to_x1(arrz1,arrx2)
    call dfftw_execute_r2r(bwd_plan_x,arrx2, arrx1)
    call transpose_x1_to_y1(arrx1,prsrc)
#endif
#endif
    
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_x1

  !******************************************************************
  ! clcPPE_x1_reduce
  !******************************************************************
  subroutine clcPPE_x1_reduce(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc,jp,jm
#ifdef FFTW_1D
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx
#else
#ifdef OverWriteFFT
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx
#else
    real(RK),dimension(x1size(1),x1size(2),x1size(3))::arrx1,arrx2
#endif
#endif
    real(RK),dimension(z2start(1):z2end(1),z2start(2):z2end(2),z2start(3):z2end(3))::arrz1,arrz2,arrz3
    real(RK),dimension(z2start(1):z2end(1),2,z2start(3):z2end(3))::arrz2_reduce
    real(RK),dimension(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3))::arry2_reduce    
    real(RK),dimension(z2start(1):z2end(1),z2start(2):z2end(2))::tridmj,tridcj,tridpj,tridfj    

#ifdef FFTW_1D
    call transpose_y1_to_x1(prsrc,arrx)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_z2(arrx,arrz2)
    call my_execute_FFTW_r2r_z(fwd_plan_z,z2size(1),z2size(2),z2size(3),arrz2)
#else
#ifdef OverWriteFFT
    call transpose_y1_to_x1(prsrc,arrx)  
    call dfftw_execute_r2r(fwd_plan_x,arrx,arrx)
    call transpose_x1_to_z2(arrx,arrz2)
    call dfftw_execute_r2r(fwd_plan_z,arrz2,arrz2)
#else
    call transpose_y1_to_x1(prsrc,arrx1)  
    call dfftw_execute_r2r(fwd_plan_x,arrx1,arrx2)
    call transpose_x1_to_z2(arrx2,arrz1)
    call dfftw_execute_r2r(fwd_plan_z,arrz1,arrz2)
#endif
#endif

    IF(nrank==0) THEN
      do kc=2,z2end(3)
        do jc=z2start(2),z2end(2)
          do ic=z2start(1),z2end(1)
            tridmj(ic,jc)=am2Pr(jc)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
            tridpj(ic,jc)=ap2Pr(jc)
            tridfj(ic,jc)=arrz2(ic,jc,kc)
          enddo
        enddo      
        do ic=z2start(1),z2end(1)
          jc=z2start(2)
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
          jc=z2start(2)+1
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        enddo
        do jc=z2start(2)+2,z2end(2)
          jm=jc-1
          do ic=z2start(1),z2end(1)
            rTemp=one/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
            tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
            tridpj(ic,jc)= rTemp*tridpj(ic,jc)
            tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
          enddo
        enddo
        do jc=z2end(2)-2,z2start(2)+1,-1
          jp=jc+1
          do ic=z2start(1),z2end(1)
            tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
            tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
            tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
          enddo
        enddo
        do ic=z2start(1),z2end(1)
          jc=z2start(2); jp=jc+1
          rTemp=one/(one-tridmj(ic,jp)*tridpj(ic,jc))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
          tridmj(ic,jc)= rTemp*tridmj(ic,jc)
          tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
          arrz2_reduce(ic,1,kc)=tridfj(ic,z2start(2))
          arrz2_reduce(ic,2,kc)=tridfj(ic,z2end(2))      
        enddo
        do jc=z2start(2),z2end(2)
          do ic=z2start(1),z2end(1)
            arrz1(ic,jc,kc)=tridmj(ic,jc)
            arrz3(ic,jc,kc)=tridpj(ic,jc)
            arrz2(ic,jc,kc)=tridfj(ic,jc)
          enddo
        enddo
      enddo
      
      ! for nrank=0, z2start(1)=z2start(2)=1
      ic=1; jc=1; kc=1 
      tridmj(ic,jc)=zero
      tridcj(ic,jc)=one
      tridpj(ic,jc)=zero
      tridfj(ic,jc)=zero
      do ic=2,z2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridpj(ic,jc)=ap2Pr(jc)
        tridfj(ic,jc)=arrz2(ic,jc,kc)
      enddo
      do jc=2,z2end(2)
        do ic=z2start(1),z2end(1)
          tridmj(ic,jc)=am2Pr(jc)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridpj(ic,jc)=ap2Pr(jc)
          tridfj(ic,jc)=arrz2(ic,jc,kc)
        enddo
      enddo            
      do ic=z2start(1),z2end(1)
        jc=z2start(2)
        tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
        tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
        tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        jc=z2start(2)+1
        tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
        tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
        tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
      enddo
      do jc=z2start(2)+2,z2end(2)
        jm=jc-1
        do ic=z2start(1),z2end(1)
          rTemp=one/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
          tridpj(ic,jc)= rTemp*tridpj(ic,jc)
          tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
        enddo
      enddo
      do jc=z2end(2)-2,z2start(2)+1,-1
        jp=jc+1
        do ic=z2start(1),z2end(1)
          tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
          tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
          tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
        enddo
      enddo
      do ic=z2start(1),z2end(1)
        jc=z2start(2); jp=jc+1
        rTemp=one/(one-tridmj(ic,jp)*tridpj(ic,jc))
        tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
        tridmj(ic,jc)= rTemp*tridmj(ic,jc)
        tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
        arrz2_reduce(ic,1,kc)=tridfj(ic,z2start(2))
        arrz2_reduce(ic,2,kc)=tridfj(ic,z2end(2))      
      enddo
      do jc=z2start(2),z2end(2)
        do ic=z2start(1),z2end(1)
          arrz1(ic,jc,kc)=tridmj(ic,jc)
          arrz3(ic,jc,kc)=tridpj(ic,jc)
          arrz2(ic,jc,kc)=tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=z2start(3),z2end(3)
        do jc=z2start(2),z2end(2)
          do ic=z2start(1),z2end(1)
            tridmj(ic,jc)=am2Pr(jc)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
            tridpj(ic,jc)=ap2Pr(jc)
            tridfj(ic,jc)=arrz2(ic,jc,kc)
          enddo
        enddo      
        do ic=z2start(1),z2end(1)
          jc=z2start(2)
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
          jc=z2start(2)+1
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        enddo
        do jc=z2start(2)+2,z2end(2)
          jm=jc-1
          do ic=z2start(1),z2end(1)
            rTemp=one/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
            tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
            tridpj(ic,jc)= rTemp*tridpj(ic,jc)
            tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
          enddo
        enddo
        do jc=z2end(2)-2,z2start(2)+1,-1
          jp=jc+1
          do ic=z2start(1),z2end(1)
            tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
            tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
            tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
          enddo
        enddo
        do ic=z2start(1),z2end(1)
          jc=z2start(2); jp=jc+1
          rTemp=one/(one-tridmj(ic,jp)*tridpj(ic,jc))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
          tridmj(ic,jc)= rTemp*tridmj(ic,jc)
          tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
          arrz2_reduce(ic,1,kc)=tridfj(ic,z2start(2))
          arrz2_reduce(ic,2,kc)=tridfj(ic,z2end(2))      
        enddo
        do jc=z2start(2),z2end(2)
          do ic=z2start(1),z2end(1)
            arrz1(ic,jc,kc)=tridmj(ic,jc)
            arrz3(ic,jc,kc)=tridpj(ic,jc)
            arrz2(ic,jc,kc)=tridfj(ic,jc)
          enddo
        enddo
      enddo
    ENDIF
    call transpose_z2_to_y2(arrz2_reduce,arry2_reduce,decomp_PPE)
    do kc=y2start(3),y2end(3)
      call InverseTridiagonal_reduce(a_reduce(:,:,kc),c_reduce(:,:,kc),arry2_reduce(:,:,kc),y2size(1),2*p_row)
    enddo
    call transpose_y2_to_z2(arry2_reduce,arrz2_reduce,decomp_PPE)
    do kc=z2start(3),z2end(3)
      do ic=z2start(1),z2end(1)
        arrz2(ic,z2start(2),kc)=arrz2_reduce(ic,1,kc)
        arrz2(ic,z2end(2),  kc)=arrz2_reduce(ic,2,kc)
      enddo
      do jc=z2start(2)+1,z2end(2)-1
        do ic=z2start(1),z2end(1)    
          arrz2(ic,jc,kc)=arrz2(ic,jc,kc)-arrz1(ic,jc,kc)*arrz2_reduce(ic,1,kc)-arrz3(ic,jc,kc)*arrz2_reduce(ic,2,kc)
        enddo
      enddo
    enddo

#ifdef FFTW_1D
    call my_execute_FFTW_r2r_z(bwd_plan_z,z2size(1),z2size(2),z2size(3),arrz2)
    call transpose_z2_to_x1(arrz2,arrx)
    call my_execute_FFTW_r2r_x(bwd_plan_x,x1size(1),x1size(2),x1size(3),arrx)
    call transpose_x1_to_y1(arrx,prsrc)
#else
#ifdef OverWriteFFT
    call dfftw_execute_r2r(bwd_plan_z,arrz2, arrz2)
    call transpose_z2_to_x1(arrz2,arrx)
    call dfftw_execute_r2r(bwd_plan_x,arrx, arrx)
    call transpose_x1_to_y1(arrx,prsrc)
#else
    call dfftw_execute_r2r(bwd_plan_z,arrz2, arrz1)
    call transpose_z2_to_x1(arrz1,arrx2)
    call dfftw_execute_r2r(bwd_plan_x,arrx2, arrx1)
    call transpose_x1_to_y1(arrx1,prsrc)
#endif
#endif
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_x1_reduce

  !******************************************************************
  ! clcPPE_z1_periodic
  !******************************************************************
  subroutine clcPPE_z1_periodic(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    integer::ic,jc,kc
#ifdef FFTW_1D
    real(RK),dimension(x2size(1),x2size(2),x2size(3))::arrx
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz
#else
#ifdef OverWriteFFT
    real(RK),dimension(x2size(1),x2size(2),x2size(3))::arrx
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz
#else
    real(RK),dimension(x2size(1),x2size(2),x2size(3))::arrx1,arrx2
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1,arrz2
#endif
#endif
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2),y2start(3):y2end(3))::prTemp
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2))::tridmj,tridcj,tridpj,tridfj

#ifdef FFTW_1D
    call transpose_y1_to_z1(prsrc,arrz)
    call my_execute_FFTW_r2r_z(fwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_x2(arrz,arrx)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x2size(1),x2size(2),x2size(3),arrx)
    call transpose_x2_to_y2(arrx,prTemp)
#else
#ifdef OverWriteFFT
    call transpose_y1_to_z1(prsrc,arrz)
    call dfftw_execute_r2r(fwd_plan_z,arrz,arrz)
    call transpose_z1_to_x2(arrz,arrx)
    call dfftw_execute_r2r(fwd_plan_x,arrx,arrx)
    call transpose_x2_to_y2(arrx,prTemp)
#else
    call transpose_y1_to_z1(prsrc,arrz1)
    call dfftw_execute_r2r(fwd_plan_z,arrz1,arrz2)
    call transpose_z1_to_x2(arrz2,arrx1)
    call dfftw_execute_r2r(fwd_plan_x,arrx1,arrx2)
    call transpose_x2_to_y2(arrx2,prTemp)
#endif
#endif
    do jc=y2start(2),y2end(2)
      do ic=y2start(1),y2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridpj(ic,jc)=ap2Pr(jc)
      enddo
    enddo
    IF(nrank==0) THEN
      do kc=2,y2end(3)
        do jc=1,nyc
          do ic=1,y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo

      ! for nrank=0, y1start(1)=y1start(3)=1
      ic=1; jc=1; kc=1;
      tridpj(ic,jc)=zero
      tridcj(ic,jc)=one
      tridmj(ic,jc)=zero
      tridfj(ic,jc)=zero
      do ic=2,y2end(1)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridfj(ic,jc)=prTemp(ic,jc,kc)
      enddo
      do jc=2,nyc
        do ic=1,y2end(1)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridfj(ic,jc)=prTemp(ic,jc,kc)
        enddo
      enddo
      call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,tridfj,y2size(1),nyc)
      do jc=1,nyc
        do ic=1,y2end(1)
          prTemp(ic,jc,kc)= tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=y2start(3),y2end(3)
        do jc=y2start(2),y2end(2)
          do ic=y2start(1),y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InversePeriodicTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo       
    ENDIF
#ifdef FFTW_1D
    call transpose_y2_to_x2(prTemp,arrx)
    call my_execute_FFTW_r2r_x(bwd_plan_x,x2size(1),x2size(2),x2size(3),arrx)
    call transpose_x2_to_z1(arrx,arrz)
    call my_execute_FFTW_r2r_z(bwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_y1(arrz,prsrc)
#else
#ifdef OverWriteFFT
    call transpose_y2_to_x2(prTemp,arrx)
    call dfftw_execute_r2r(bwd_plan_x,arrx, arrx)
    call transpose_x2_to_z1(arrx,arrz)
    call dfftw_execute_r2r(bwd_plan_z,arrz, arrz)
    call transpose_z1_to_y1(arrz,prsrc)
#else
    call transpose_y2_to_x2(prTemp,arrx2)
    call dfftw_execute_r2r(bwd_plan_x,arrx2, arrx1)
    call transpose_x2_to_z1(arrx1,arrz2)
    call dfftw_execute_r2r(bwd_plan_z,arrz2, arrz1)
    call transpose_z1_to_y1(arrz1,prsrc)
#endif
#endif
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_z1_periodic
    
  !******************************************************************
  ! clcPPE_z1
  !******************************************************************
  subroutine clcPPE_z1(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    integer::ic,jc,kc
#ifdef FFTW_1D
    real(RK),dimension(x2size(1),x2size(2),x2size(3))::arrx
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz
#else
#ifdef OverWriteFFT
    real(RK),dimension(x2size(1),x2size(2),x2size(3))::arrx
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz
#else
    real(RK),dimension(x2size(1),x2size(2),x2size(3))::arrx1,arrx2
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1,arrz2
#endif
#endif
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2),y2start(3):y2end(3))::prTemp
    real(RK),dimension(y2start(1):y2end(1),y2start(2):y2end(2))::tridmj,tridcj,tridpj,tridfj

#ifdef FFTW_1D
    call transpose_y1_to_z1(prsrc,arrz)
    call my_execute_FFTW_r2r_z(fwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_x2(arrz,arrx)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x2size(1),x2size(2),x2size(3),arrx)
    call transpose_x2_to_y2(arrx,prTemp)
#else
#ifdef OverWriteFFT
    call transpose_y1_to_z1(prsrc,arrz)
    call dfftw_execute_r2r(fwd_plan_z,arrz,arrz)
    call transpose_z1_to_x2(arrz,arrx)
    call dfftw_execute_r2r(fwd_plan_x,arrx,arrx)
    call transpose_x2_to_y2(arrx,prTemp)
#else
    call transpose_y1_to_z1(prsrc,arrz1)
    call dfftw_execute_r2r(fwd_plan_z,arrz1,arrz2)
    call transpose_z1_to_x2(arrz2,arrx1)
    call dfftw_execute_r2r(fwd_plan_x,arrx1,arrx2)
    call transpose_x2_to_y2(arrx2,prTemp)
#endif
#endif
    do jc=y2start(2),y2end(2)
      do ic=y2start(1),y2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridpj(ic,jc)=ap2Pr(jc)
      enddo
    enddo
    IF(nrank==0) THEN
      do kc=2,y2end(3)
        do jc=1,nyc
          do ic=1,y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InverseTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo

      ! for nrank=0, y1start(1)=y1start(3)=1
      ic=1; jc=1; kc=1;
      tridpj(ic,jc)=zero
      tridcj(ic,jc)=one
      tridmj(ic,jc)=zero
      tridfj(ic,jc)=zero
      do ic=2,y2end(1)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridfj(ic,jc)=prTemp(ic,jc,kc)
      enddo
      do jc=2,nyc
        do ic=1,y2end(1)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridfj(ic,jc)=prTemp(ic,jc,kc)
        enddo
      enddo
      call InverseTridiagonal(tridmj,tridcj,tridpj,tridfj,y2size(1),nyc)
      do jc=1,nyc
        do ic=1,y2end(1)
          prTemp(ic,jc,kc)= tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=y2start(3),y2end(3)
        do jc=y2start(2),y2end(2)
          do ic=y2start(1),y2end(1)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          enddo
        enddo
        call InverseTridiagonal(tridmj,tridcj,tridpj,prTemp(:,:,kc),y2size(1),nyc)
      enddo       
    ENDIF
#ifdef FFTW_1D
    call transpose_y2_to_x2(prTemp,arrx)
    call my_execute_FFTW_r2r_x(bwd_plan_x,x2size(1),x2size(2),x2size(3),arrx)
    call transpose_x2_to_z1(arrx,arrz)
    call my_execute_FFTW_r2r_z(bwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_y1(arrz,prsrc)
#else    
#ifdef OverWriteFFT
    call transpose_y2_to_x2(prTemp,arrx)
    call dfftw_execute_r2r(bwd_plan_x,arrx, arrx)
    call transpose_x2_to_z1(arrx,arrz)
    call dfftw_execute_r2r(bwd_plan_z,arrz, arrz)
    call transpose_z1_to_y1(arrz,prsrc)
#else
    call transpose_y2_to_x2(prTemp,arrx2)
    call dfftw_execute_r2r(bwd_plan_x,arrx2, arrx1)
    call transpose_x2_to_z1(arrx1,arrz2)
    call dfftw_execute_r2r(bwd_plan_z,arrz2, arrz1)
    call transpose_z1_to_y1(arrz1,prsrc)
#endif
#endif
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_z1

  !******************************************************************
  ! clcPPE_z1_reduce
  !******************************************************************
  subroutine clcPPE_z1_reduce(prsrc, prphiHalo)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::prphiHalo
      
    ! locals
    real(RK)::rTemp
    integer::ic,jc,kc,jp,jm
#ifdef FFTW_1D
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz
#else
#ifdef OverWriteFFT
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz
#else
    real(RK),dimension(z1size(1),z1size(2),z1size(3))::arrz1,arrz2
#endif
#endif
    real(RK),dimension(x2start(1):x2end(1),2,x2start(3):x2end(3))::arrx2_reduce
    real(RK),dimension(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3))::arry2_reduce    
    real(RK),dimension(x2start(1):x2end(1),x2start(2):x2end(2))::tridmj,tridcj,tridpj,tridfj
    real(RK),dimension(x2start(1):x2end(1),x2start(2):x2end(2),x2start(3):x2end(3))::arrx1,arrx2,arrx3

#ifdef FFTW_1D
    call transpose_y1_to_z1(prsrc,arrz)  
    call my_execute_FFTW_r2r_z(fwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_x2(arrz,arrx2)
    call my_execute_FFTW_r2r_x(fwd_plan_x,x2size(1),x2size(2),x2size(3),arrx2)
#else
#ifdef OverWriteFFT
    call transpose_y1_to_z1(prsrc,arrz)  
    call dfftw_execute_r2r(fwd_plan_z,arrz,arrz)
    call transpose_z1_to_x2(arrz,arrx2)
    call dfftw_execute_r2r(fwd_plan_x,arrx2,arrx2)
#else
    call transpose_y1_to_z1(prsrc,arrz1)  
    call dfftw_execute_r2r(fwd_plan_z,arrz1,arrz2)
    call transpose_z1_to_x2(arrz2,arrx1)
    call dfftw_execute_r2r(fwd_plan_x,arrx1,arrx2)
#endif    
#endif
    IF(nrank==0) THEN      
      do kc=2,x2end(3)
        do jc=x2start(2),x2end(2)
          do ic=x2start(1),x2end(1)
            tridmj(ic,jc)=am2Pr(jc)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
            tridpj(ic,jc)=ap2Pr(jc)
            tridfj(ic,jc)=arrx2(ic,jc,kc)
          enddo
        enddo
        do ic=x2start(1),x2end(1)
          jc=x2start(2)
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
          jc=x2start(2)+1
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        enddo
        do jc=x2start(2)+2,x2end(2)
          jm=jc-1
          do ic=x2start(1),x2end(1)
            rTemp=one/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
            tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
            tridpj(ic,jc)= rTemp*tridpj(ic,jc)
            tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
          enddo
        enddo
        do jc=x2end(2)-2,x2start(2)+1,-1
          jp=jc+1
          do ic=x2start(1),x2end(1)
            tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
            tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
            tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
          enddo
        enddo
        do ic=x2start(1),x2end(1)
          jc=x2start(2); jp=jc+1
          rTemp=one/(one-tridmj(ic,jp)*tridpj(ic,jc))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
          tridmj(ic,jc)= rTemp*tridmj(ic,jc)
          tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
          arrx2_reduce(ic,1,kc)=tridfj(ic,x2start(2))
          arrx2_reduce(ic,2,kc)=tridfj(ic,x2end(2))       
        enddo
        do jc=x2start(2),x2end(2)
          do ic=x2start(1),x2end(1)
            arrx1(ic,jc,kc)=tridmj(ic,jc)
            arrx3(ic,jc,kc)=tridpj(ic,jc)
            arrx2(ic,jc,kc)=tridfj(ic,jc)
          enddo
        enddo
      enddo  

      ! for nrank=0, x2start(2)=x2start(3)=1
      ic=1; jc=1; kc=1
      tridmj(ic,jc)=zero
      tridcj(ic,jc)=one
      tridpj(ic,jc)=zero
      tridfj(ic,jc)=zero
      do ic=2,x2end(1)
        tridmj(ic,jc)=am2Pr(jc)
        tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
        tridpj(ic,jc)=ap2Pr(jc)
        tridfj(ic,jc)=arrx2(ic,jc,kc)
      enddo
      do jc=2,x2end(2)
        do ic=1,x2end(1)
          tridmj(ic,jc)=am2Pr(jc)
          tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
          tridpj(ic,jc)=ap2Pr(jc)
          tridfj(ic,jc)=arrx2(ic,jc,kc)
        enddo
      enddo
      do ic=x2start(1),x2end(1)
        jc=x2start(2)
        tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
        tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
        tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        jc=x2start(2)+1
        tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
        tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
        tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
      enddo
      do jc=x2start(2)+2,x2end(2)
        jm=jc-1
        do ic=x2start(1),x2end(1)
          rTemp=one/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
          tridpj(ic,jc)= rTemp*tridpj(ic,jc)
          tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
        enddo
      enddo
      do jc=x2end(2)-2,x2start(2)+1,-1
        jp=jc+1
        do ic=x2start(1),x2end(1)
          tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
          tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
          tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
        enddo
      enddo
      do ic=x2start(1),x2end(1)
        jc=x2start(2); jp=jc+1
        rTemp=one/(one-tridmj(ic,jp)*tridpj(ic,jc))
        tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
        tridmj(ic,jc)= rTemp*tridmj(ic,jc)
        tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
        arrx2_reduce(ic,1,kc)=tridfj(ic,x2start(2))
        arrx2_reduce(ic,2,kc)=tridfj(ic,x2end(2))       
      enddo
      do jc=x2start(2),x2end(2)
        do ic=x2start(1),x2end(1)
          arrx1(ic,jc,kc)=tridmj(ic,jc)
          arrx3(ic,jc,kc)=tridpj(ic,jc)
          arrx2(ic,jc,kc)=tridfj(ic,jc)
        enddo
      enddo
    ELSE
      do kc=x2start(3),x2end(3)
        do jc=x2start(2),x2end(2)
          do ic=x2start(1),x2end(1)
            tridmj(ic,jc)=am2Pr(jc)
            tridcj(ic,jc)=ac2Pr(jc)+WaveNumX(ic)+WaveNumZ(kc)
            tridpj(ic,jc)=ap2Pr(jc)
            tridfj(ic,jc)=arrx2(ic,jc,kc)
          enddo
        enddo
        do ic=x2start(1),x2end(1)
          jc=x2start(2)
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
          jc=x2start(2)+1
          tridmj(ic,jc)=tridmj(ic,jc)/tridcj(ic,jc)
          tridpj(ic,jc)=tridpj(ic,jc)/tridcj(ic,jc)
          tridfj(ic,jc)=tridfj(ic,jc)/tridcj(ic,jc)
        enddo
        do jc=x2start(2)+2,x2end(2)
          jm=jc-1
          do ic=x2start(1),x2end(1)
            rTemp=one/(tridcj(ic,jc)-tridmj(ic,jc)*tridpj(ic,jm))
            tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridmj(ic,jc)*tridfj(ic,jm))
            tridpj(ic,jc)= rTemp*tridpj(ic,jc)
            tridmj(ic,jc)=-rTemp*tridmj(ic,jc)*tridmj(ic,jm)
          enddo
        enddo
        do jc=x2end(2)-2,x2start(2)+1,-1
          jp=jc+1
          do ic=x2start(1),x2end(1)
            tridfj(ic,jc)= tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp)
            tridmj(ic,jc)= tridmj(ic,jc)-tridpj(ic,jc)*tridmj(ic,jp)
            tridpj(ic,jc)=-tridpj(ic,jc)*tridpj(ic,jp)
          enddo
        enddo
        do ic=x2start(1),x2end(1)
          jc=x2start(2); jp=jc+1
          rTemp=one/(one-tridmj(ic,jp)*tridpj(ic,jc))
          tridfj(ic,jc)= rTemp*(tridfj(ic,jc)-tridpj(ic,jc)*tridfj(ic,jp))
          tridmj(ic,jc)= rTemp*tridmj(ic,jc)
          tridpj(ic,jc)=-rTemp*tridpj(ic,jc)*tridpj(ic,jp) 
          arrx2_reduce(ic,1,kc)=tridfj(ic,x2start(2))
          arrx2_reduce(ic,2,kc)=tridfj(ic,x2end(2))       
        enddo
        do jc=x2start(2),x2end(2)
          do ic=x2start(1),x2end(1)
            arrx1(ic,jc,kc)=tridmj(ic,jc)
            arrx3(ic,jc,kc)=tridpj(ic,jc)
            arrx2(ic,jc,kc)=tridfj(ic,jc)
          enddo
        enddo
      enddo
    ENDIF
    call transpose_x2_to_y2(arrx2_reduce,arry2_reduce,decomp_PPE)
    do kc=y2start(3),y2end(3)
      call InverseTridiagonal_reduce(a_reduce(:,:,kc),c_reduce(:,:,kc),arry2_reduce(:,:,kc),y2size(1),2*p_col)
    enddo
    call transpose_y2_to_x2(arry2_reduce,arrx2_reduce,decomp_PPE)
    do kc=x2start(3),x2end(3)
      do ic=x2start(1),x2end(1)
        arrx2(ic,x2start(2),kc)=arrx2_reduce(ic,1,kc)
        arrx2(ic,x2end(2),  kc)=arrx2_reduce(ic,2,kc)
      enddo
      do jc=x2start(2)+1,x2end(2)-1
        do ic=x2start(1),x2end(1)    
          arrx2(ic,jc,kc)=arrx2(ic,jc,kc)-arrx1(ic,jc,kc)*arrx2_reduce(ic,1,kc)-arrx3(ic,jc,kc)*arrx2_reduce(ic,2,kc)
        enddo
      enddo
    enddo

#ifdef FFTW_1D
    call my_execute_FFTW_r2r_x(bwd_plan_x,x2size(1),x2size(2),x2size(3),arrx2)
    call transpose_x2_to_z1(arrx2,arrz)
    call my_execute_FFTW_r2r_z(bwd_plan_z,z1size(1),z1size(2),z1size(3),arrz)
    call transpose_z1_to_y1(arrz,prsrc)
#else
#ifdef OverWriteFFT
    call dfftw_execute_r2r(bwd_plan_x,arrx2,arrx2)
    call transpose_x2_to_z1(arrx2,arrz)
    call dfftw_execute_r2r(bwd_plan_z,arrz,arrz)
    call transpose_z1_to_y1(arrz,prsrc)
#else
    call dfftw_execute_r2r(bwd_plan_x,arrx2,arrx1)
    call transpose_x2_to_z1(arrx1,arrz2)
    call dfftw_execute_r2r(bwd_plan_z,arrz2,arrz1)
    call transpose_z1_to_y1(arrz1,prsrc)
#endif
#endif
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          prphiHalo(ic,jc,kc) = prsrc(ic,jc,kc)*normfft
        enddo
      enddo
    enddo
  end subroutine clcPPE_z1_reduce

  !******************************************************************
  ! InverseTridiagonal_reduce
  !******************************************************************  
  subroutine InverseTridiagonal_reduce(aj,cj,fj,m,n)
    implicit none
    integer,intent(in)::m,n
    real(RK),dimension(m,n),intent(in):: aj,cj
    real(RK),dimension(m,n),intent(inout)::fj
    
    ! locals
    integer:: i,j
    real(RK),dimension(m):: vecm
    real(RK),dimension(m,n)::arrmn
      
    do i=1,m
      vecm(i)=one
    enddo
    do j=2,n
      do i=1,m
        arrmn(i,j)=cj(i,j-1)/vecm(i)
        vecm(i)=one-aj(i,j)*arrmn(i,j)
        fj(i,j)=(fj(i,j)-aj(i,j)*fj(i,j-1))/vecm(i)
      enddo
    enddo
    do j=n-1,1,-1
      do i=1,m
        fj(i,j)=fj(i,j)-arrmn(i,j+1)*fj(i,j+1)
       enddo
    enddo    
  end subroutine InverseTridiagonal_reduce
  
#ifdef FFTW_1D
  !******************************************************************
  ! my_execute_FFTW_r2r_x
  !************************************************************
  subroutine my_execute_FFTW_r2r_x(plan_flag,nx_FFT,ny_FFT,nz_FFT,arrFFT)
    implicit none
    type(C_PTR),intent(in)::plan_flag
    integer,intent(in)::nx_FFT,ny_FFT,nz_FFT
    real(RK),dimension(nx_FFT,ny_FFT,nz_FFT),intent(inout)::arrFFT
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(nx_FFT)::VecIn,VecOut
    
    do kc=1,nz_FFT
      do jc=1,ny_FFT
        do ic=1,nx_FFT
          VecIn(ic)=arrFFT(ic,jc,kc)
        enddo
        call dfftw_execute_r2r(plan_flag,VecIn,VecOut)
        do ic=1,nx_FFT
          arrFFT(ic,jc,kc)=VecOut(ic)
        enddo
      enddo
    enddo
  end subroutine my_execute_FFTW_r2r_x

  !******************************************************************
  ! my_execute_FFTW_r2r_z
  !************************************************************
  subroutine my_execute_FFTW_r2r_z(plan_flag,nx_FFT,ny_FFT,nz_FFT,arrFFT)
    implicit none
    type(C_PTR),intent(in)::plan_flag
    integer,intent(in)::nx_FFT,ny_FFT,nz_FFT
    real(RK),dimension(nx_FFT,ny_FFT,nz_FFT),intent(inout)::arrFFT
    
    ! locals
    integer::ic,jc,kc
    real(RK),dimension(nz_FFT)::VecIn,VecOut
    
    do jc=1,ny_FFT
      do ic=1,nx_FFT
        do kc=1,nz_FFT
          VecIn(kc)=arrFFT(ic,jc,kc)
        enddo
        call dfftw_execute_r2r(plan_flag,VecIn,VecOut)
        do kc=1,nz_FFT
          arrFFT(ic,jc,kc)=VecOut(kc)
        enddo        
      enddo
    enddo
  end subroutine my_execute_FFTW_r2r_z
#endif
end module m_Poisson
