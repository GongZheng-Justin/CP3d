#include "definitions_inc.f90"
module f2_Poisson
  use MPI
  use iso_c_binding
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d
  use f2_Parameters
  use f2_MeshAndMetries
  use f2_Variables,only: mb1
  use f2_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal
  implicit none
  private
  include "myfftw3.f03"

  real(RK)::normfft  
  type(decomp_info),allocatable::decomp_PPE
  real(RK),allocatable,dimension(:,:,:)::a_reduce,c_reduce
  real(RK),allocatable,dimension(:)::WaveNumX,WaveNumY,WaveNumZ
  procedure(),pointer::clcPPE=>null(),execute_FFTW_r2r_z=>null()
  type(C_PTR)::fwd_plan_x,bwd_plan_x,fwd_plan_y,bwd_plan_y,fwd_plan_z,bwd_plan_z
  
  public:: InitPoissonSolver,clcPPE, Destory_Poisson_FFT_Plan
contains
#define nTime_FFT_Test 5

#define my_FFTW_inc_add_y
#define my_FFTW_inc_add_z2
#include "my_FFTW_inc.f90"
#undef  my_FFTW_inc_add_z2
#undef  my_FFTW_inc_add_y

#define my_Poisson_inc_add_Periodic_2d
#define my_Poisson_inc_add_UniformYmesh
#include "my_Poisson_inc.f90"
#undef  my_Poisson_inc_add_UniformYmesh
#undef  my_Poisson_inc_add_Periodic_2d
    
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitPoissonSolver
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90     
  subroutine InitPoissonSolver(); implicit none
    
    ! locals
    logical::IsUniformYmesh
    integer(C_FFTW_R2R_KIND)::kind_fwd,kind_bwd
    real(RK),dimension(:),allocatable::Vec1,Vec2
    integer::i,j,k,iErr01,iErr02,iErr03,iChoice,iFFTz,IsReduce,IsReduceR,plan_type
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm)::prphiHalo
    real(RK)::normfftX,normfftY,normfftZ,WaveCoe,wa1,wa3,best_time,t1(2),t2(2),normTmp
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::prsrc

    if(FFTW_plan_type == 1) then
      plan_type = FFTW_PATIENT
    else
      plan_type = FFTW_ESTIMATE
    endif
            
    ! Modified wave number in x-dir
    allocate(WaveNumX(nxc), Stat =iErr01)
    IF(BcOption(xm_dir)==BC_PERIOD) THEN    
      WaveCoe = PI*2.0_RK
      normTmp = 1.0_RK
    ELSE
      WaveCoe = PI
      normTmp = 2.0_RK
    ENDIF
    do i=1,nxc
      wa1= WaveCoe*real(i-1,RK)/real(nxc,RK)
      WaveNumX(i)=2.0_RK*rdx2*(cos(wa1)-1.0_RK)
    enddo
    normfftX= 1.0_RK/(normTmp*real(nxc,RK))
    ! FFT in x
    IF(BcOption(xm_dir)==BC_PERIOD) THEN
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
    ELSE
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01
    ENDIF
    allocate(Vec1(nxc), Vec2(nxc))
    fwd_plan_x= fftw_plan_r2r_1d(nxc,Vec1,Vec2,kind_fwd,plan_type)
    bwd_plan_x= fftw_plan_r2r_1d(nxc,Vec1,Vec2,kind_bwd,plan_type)
    deallocate(Vec1,Vec2)
    
    ! Modified wave number in y-dir, maynot be used.
    allocate(WaveNumY(nyc), Stat =iErr02)
    IF(BcOption(ym_dir)==BC_PERIOD) THEN    
      WaveCoe = PI*2.0_RK
      normTmp = 1.0_RK
    ELSE
      WaveCoe = PI
      normTmp = 2.0_RK
    ENDIF
    do j=1,nyc
      wa1= WaveCoe*real(j-1,RK)/real(nyc,RK)
      WaveNumY(j)=2.0_RK*rdy2*(cos(wa1)-1.0_RK)
    enddo
    normfftY= 1.0_RK/(normTmp*real(nyc,RK))
    ! FFT in y
    IF(BcOption(ym_dir)==BC_PERIOD) THEN
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
    ELSE
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01
    ENDIF
    allocate(Vec1(nyc), Vec2(nyc))
    fwd_plan_y= fftw_plan_r2r_1d(nyc,Vec1,Vec2,kind_fwd,plan_type)
    bwd_plan_y= fftw_plan_r2r_1d(nyc,Vec1,Vec2,kind_bwd,plan_type)
    deallocate(Vec1,Vec2)
        
    ! Modified wave number in z-dir
    allocate(WaveNumZ(nzc), Stat =iErr03)
    IF(BcOption(zm_dir)==BC_PERIOD) THEN    
      WaveCoe = PI*2.0_RK
      normTmp = 1.0_RK
    ELSE
      WaveCoe = PI
      normTmp = 2.0_RK
    ENDIF   
    do k=1,nzc
      wa3= WaveCoe*real(k-1,RK)/real(nzc,RK)
      WaveNumZ(k)=2.0_RK*rdz2*(cos(wa3)-1.0_RK)
    enddo
    normFFTZ= 1.0_RK/(normTmp*real(nzc,RK))
    ! FFT in z
    IF(BcOption(zm_dir)==BC_PERIOD ) THEN
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
    ELSE
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01    
    ENDIF
    allocate(Vec1(nzc), Vec2(nzc))
    fwd_plan_z= fftw_plan_r2r_1d(nzc,Vec1,Vec2,kind_fwd,plan_type)
    bwd_plan_z= fftw_plan_r2r_1d(nzc,Vec1,Vec2,kind_bwd,plan_type)
    deallocate(Vec1,Vec2)
    
    best_time=max(huge(wa1),1.0E+20)
    if(nrank==0) call MainLog%OutInfo("Auto-tuning mode for Poisson Solver......",1)
    IF(BcOption(ym_dir)==BC_PERIOD) THEN
      ! Choice-1
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      !      
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
      t1(1)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_x1_periodic(prsrc,prphiHalo)
      enddo
      t2(1)=MPI_WTIME()-t1(1)
      nullify(execute_FFTW_r2r_z)
      ! 
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
      t1(2)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_x1_periodic(prsrc,prphiHalo)
      enddo
      t2(2)=MPI_WTIME()-t1(2)
      nullify(execute_FFTW_r2r_z)
      !       
      call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0)call MainLog%OutInfo("Choice-1, time= "//strip(num2str(t1(1)))//", "//strip(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=1; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=1; iFFTz=2
      endif
   
      ! Choice-2
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      !      
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
      t1(1)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_z1_periodic(prsrc,prphiHalo)
      enddo
      t2(1)=MPI_WTIME()-t1(1)
      nullify(execute_FFTW_r2r_z)
      !      
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
      t1(2)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_z1_periodic(prsrc,prphiHalo)
      enddo
      t2(2)=MPI_WTIME()-t1(2)
      nullify(execute_FFTW_r2r_z)
      !       
      call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0)call MainLog%OutInfo("Choice-2, time= "//strip(num2str(t1(1)))//", "//strip(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=2; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=2; iFFTz=2
      endif
    ELSE
      ! Choice-1
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      !
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
      t1(1)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_x1(prsrc,prphiHalo)
      enddo
      t2(1)=MPI_WTIME()-t1(1)
      nullify(execute_FFTW_r2r_z)
      ! 
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
      t1(2)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_x1(prsrc,prphiHalo)
      enddo
      t2(2)=MPI_WTIME()-t1(2)
      nullify(execute_FFTW_r2r_z)
      !      
      call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0) call MainLog%OutInfo("Choice-1, time= "//strip(num2str(t1(1)))//", "//strip(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=1; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=1; iFFTz=2
      endif
   
      ! Choice-2
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      !
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
      t1(1)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_z1(prsrc,prphiHalo)
      enddo
      t2(1)=MPI_WTIME()-t1(1)
      nullify(execute_FFTW_r2r_z)
      ! 
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
      t1(2)=MPI_WTIME()
      do k=1,nTime_FFT_Test
        call clcPPE_z1(prsrc,prphiHalo)
      enddo
      t2(2)=MPI_WTIME()-t1(2)
      nullify(execute_FFTW_r2r_z)
      !
      call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0) call MainLog%OutInfo("Choice-2, time= "//strip(num2str(t1(1)))//", "//strip(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=2; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=2; iFFTz=2
      endif

      ! Choice-3
      IsReduce=1
      if(z2size(2)<4) IsReduce=0
      call MPI_ALLREDUCE(IsReduce,IsReduceR,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD,iErr01)
      if(IsReduceR==1) then
        allocate(decomp_PPE)
        allocate(a_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
        allocate(c_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
        call Initialize_ReduceMatrix(z2start,z2end,'z')
        call MPI_BARRIER(MPI_COMM_WORLD,iErr01)        
        !
        execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
        t1(1)=MPI_WTIME()
        do k=1,nTime_FFT_Test
          call clcPPE_x1_reduce(prsrc,prphiHalo)
        enddo
        t2(1)=MPI_WTIME()-t1(1)
        nullify(execute_FFTW_r2r_z)
        ! 
        execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
        t1(2)=MPI_WTIME()
        do k=1,nTime_FFT_Test
          call clcPPE_x1_reduce(prsrc,prphiHalo)
        enddo
        t2(2)=MPI_WTIME()-t1(2)
        nullify(execute_FFTW_r2r_z)
        !        
        call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
        if(nrank==0)call MainLog%OutInfo("Choice-3, time= "//strip(num2str(t1(1)))//", "//strip(num2str(t1(2))),2)
        if(best_time>t1(1)) then
          best_time=t1(1); iChoice=3; iFFTz=1
        endif
        if(best_time>t1(2)) then
          best_time=t1(2); iChoice=3; iFFTz=2
        endif
        call decomp_info_finalize(decomp_PPE)
        deallocate(decomp_PPE,a_reduce,c_reduce)
      else
        if(nrank==0) call MainLog%OutInfo("Choice-3, z2size<4, ignore",2)
      endif
        
      ! Choice-4
      IsReduce=1
      if(x2size(2)<4) IsReduce=0
      call MPI_ALLREDUCE(IsReduce,IsReduceR,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD,iErr01)
      if(IsReduceR==1) then   
        allocate(decomp_PPE)
        allocate(a_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
        allocate(c_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
        call Initialize_ReduceMatrix(x2start,x2end,'x')
        call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
        ! 
        execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
        t1(1)=MPI_WTIME()
        do k=1,nTime_FFT_Test
          call clcPPE_z1_reduce(prsrc,prphiHalo)
        enddo
        t2(1)=MPI_WTIME()-t1(1)
        nullify(execute_FFTW_r2r_z)
        ! 
        execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
        t1(2)=MPI_WTIME()
        do k=1,nTime_FFT_Test
          call clcPPE_z1_reduce(prsrc,prphiHalo)
        enddo
        t2(2)=MPI_WTIME()-t1(2)
        nullify(execute_FFTW_r2r_z)
        !
        call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
        if(nrank==0)call MainLog%OutInfo("Choice-4, time= "//strip(num2str(t1(1)))//", "//strip(num2str(t1(2))),2)
        if(best_time>t1(1)) then
          best_time=t1(1); iChoice=4; iFFTz=1
        endif
        if(best_time>t1(2)) then
          best_time=t1(2); iChoice=4; iFFTz=2
        endif
        call decomp_info_finalize(decomp_PPE)
        deallocate(decomp_PPE,a_reduce,c_reduce)
      else
        if(nrank==0) call MainLog%OutInfo("Choice-4, x2size<4, ignore",2)
      endif
    ENDIF
    IsUniformYmesh = .true.
    do j=1,nyc
      if(abs(dyp(j)/dyp(1) -1.0_RK) > 1.0E-11_RK) then
        IsUniformYmesh = .false.; exit
      endif
    enddo
    IF(IsUniformYmesh) THEN
      ! Choice-5
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      !      
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
      t1(1)=MPI_WTIME()
      if(BcOption(xm_dir)==BC_PERIOD) then      
        do k=1,nTime_FFT_Test
          call clcPPE_z1_UniformYmesh_periodic(prsrc,prphiHalo)
        enddo
      else        
        do k=1,nTime_FFT_Test
          call clcPPE_z1_UniformYmesh(prsrc,prphiHalo)
        enddo             
      endif
      t2(1)=MPI_WTIME()-t1(1) 
      nullify(execute_FFTW_r2r_z)
      !      
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
      t1(2)=MPI_WTIME()
      if(BcOption(xm_dir)==BC_PERIOD) then      
        do k=1,nTime_FFT_Test
          call clcPPE_z1_UniformYmesh_periodic(prsrc,prphiHalo)
        enddo
      else        
        do k=1,nTime_FFT_Test
          call clcPPE_z1_UniformYmesh(prsrc,prphiHalo)
        enddo             
      endif
      t2(2)=MPI_WTIME()-t1(2)
      nullify(execute_FFTW_r2r_z)
      !       
      call MPI_ALLREDUCE(t2,t1,2,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0) call MainLog%OutInfo("Choice-5, time= "//strip(num2str(t1(1)))//", "//strip(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=5; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=5; iFFTz=2
      endif
        
      ! Choice-6
      call MPI_BARRIER(MPI_COMM_WORLD,iErr01)
      !
      t1(1)=MPI_WTIME()
      if(BcOption(zm_dir)==BC_PERIOD) then      
        do k=1,nTime_FFT_Test
          call clcPPE_x1_UniformYmesh_periodic(prsrc,prphiHalo)
        enddo
      else        
        do k=1,nTime_FFT_Test
          call clcPPE_x1_UniformYmesh(prsrc,prphiHalo)
        enddo      
      endif
      t2(1)=MPI_WTIME()-t1(1)
      !
      call MPI_ALLREDUCE(t2,t1,1,real_type,MPI_SUM,MPI_COMM_WORLD,iErr01)
      if(nrank==0) call MainLog%OutInfo("Choice-6, time= "//strip(num2str(t1(1))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=6; iFFTz=-99
      endif           
    ENDIF

    !
    if(nrank==0) then
      call MainLog%OutInfo("The best Poisson Solver choice is probably Choice-"//num2str(iChoice),2)
      call MainLog%OutInfo("Corresponding Global Data Transpose is:",2)
    endif
    if(iFFTz==1) then
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z
    else
      execute_FFTW_r2r_z => my_execute_FFTW_r2r_z_2
    endif
    IF(BcOption(ym_dir) == BC_PERIOD) THEN
      if(iChoice==1) then
        clcPPE => clcPPE_x1_periodic
        normfft = normfftX*normfftZ
        if(nrank==0) call MainLog%OutInfo("y1 -> x1 -> z2 -> y2 -> z2 -> x1 -> y1",3)
      elseif(iChoice==2) then
        clcPPE => clcPPE_z1_periodic
        normfft = normfftX*normfftZ
        if(nrank==0) call MainLog%OutInfo("y1 -> z1 -> x2 -> y2 -> x2 -> z1 -> y1",3)
      elseif(iChoice==5) then
        if(BcOption(xm_dir)==BC_PERIOD) then
          clcPPE => clcPPE_z1_UniformYmesh_periodic
        else
          clcPPE => clcPPE_z1_UniformYmesh
        endif
        normfft = normfftY*normfftZ
        if(nrank==0) call MainLog%OutInfo("y1 -> z1 -> x2 -> z1 -> y1",3)
      elseif(iChoice==6) then
        if(BcOption(zm_dir)==BC_PERIOD) then
          clcPPE => clcPPE_x1_UniformYmesh_periodic
        else
          clcPPE => clcPPE_x1_UniformYmesh
        endif
        normfft = normfftX*normfftY
        if(nrank==0) call MainLog%OutInfo("y1 -> x1 -> z2 -> x1 -> y1",3)
      else
        if(nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitPoissonSolver", "iChoice Wrong-1 !!!")
      endif
    ELSE
      if(iChoice==1) then
        clcPPE => clcPPE_x1
        normfft = normfftX*normfftZ
        if(nrank==0) call MainLog%OutInfo("y1 -> x1 -> z2 -> y2 -> z2 -> x1 -> y1",3)
      elseif(iChoice==2) then
        clcPPE => clcPPE_z1
        normfft = normfftX*normfftZ
        if(nrank==0) call MainLog%OutInfo("y1 -> z1 -> x2 -> y2 -> x2 -> z1 -> y1",3)
      elseif(iChoice==3) then
        clcPPE => clcPPE_x1_reduce
        normfft = normfftX*normfftZ
        if(nrank==0) call MainLog%OutInfo("y1 -> x1 -> z2 -> x1 -> y1",3)
        allocate(decomp_PPE)
        allocate(a_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
        allocate(c_reduce(y2start(1):y2end(1),2*p_row,y2start(3):y2end(3)))
        call Initialize_ReduceMatrix(z2start,z2end,'z')      
      elseif(iChoice==4) then
        clcPPE => clcPPE_z1_reduce
        normfft = normfftX*normfftZ
        if(nrank==0) call MainLog%OutInfo("y1 -> z1 -> x2 -> z1 -> y1",3)
        allocate(decomp_PPE)
        allocate(a_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
        allocate(c_reduce(y2start(1):y2end(1),2*p_col,y2start(3):y2end(3)))
        call Initialize_ReduceMatrix(x2start,x2end,'x')
      elseif(iChoice==5) then
        if(BcOption(xm_dir)==BC_PERIOD) then
          clcPPE => clcPPE_z1_UniformYmesh_periodic
        else
          clcPPE => clcPPE_z1_UniformYmesh
        endif
        normfft = normfftY*normfftZ
        if(nrank==0) call MainLog%OutInfo("y1 -> z1 -> x2 -> z1 -> y1",3)
      elseif(iChoice==6) then
        if(BcOption(zm_dir)==BC_PERIOD) then
          clcPPE => clcPPE_x1_UniformYmesh_periodic
        else
          clcPPE => clcPPE_x1_UniformYmesh
        endif
        normfft = normfftX*normfftY
        if(nrank==0) call MainLog%OutInfo("y1 -> x1 -> z2 -> x1 -> y1",3)
      else
        if(nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitPoissonSolver", "iChoice Wrong-2 !!!")
      endif  
    ENDIF
    if(nrank==0) print*," "
    prsrc=0.0_RK; prphiHalo=0.0_RK
  end subroutine InitPoissonSolver

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Destory_Poisson_FFT_Plan
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Destory_Poisson_FFT_Plan(); implicit none
    nullify(execute_FFTW_r2r_z)
    call dfftw_destroy_plan(fwd_plan_x,bwd_plan_x)
    call dfftw_destroy_plan(fwd_plan_y,bwd_plan_y)
    call dfftw_destroy_plan(fwd_plan_z,bwd_plan_z)
    if(allocated(decomp_PPE)) then
      call decomp_info_finalize(decomp_PPE)
      deallocate(decomp_PPE)
    endif
  end subroutine Destory_Poisson_FFT_Plan
  
#undef nTime_FFT_Test
end module f2_Poisson
