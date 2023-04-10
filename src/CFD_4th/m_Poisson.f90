module m_Poisson
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_decomp2d
  use m_Parameters
  use iso_c_binding
  use m_MeshAndMetries
  use m_Variables,only:mb1
  use m_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal
  implicit none
  private
  include "fftw3.f03"
  
  real(RK)::normfft
  type(decomp_info),allocatable::decomp_PPE
  real(RK),allocatable,dimension(:)::WaveNumX,WaveNumZ
  type(C_PTR)::fwd_plan_x,bwd_plan_x,fwd_plan_z,bwd_plan_z
  real(RK),allocatable,dimension(:,:,:)::a_reduce,c_reduce
  procedure(),pointer::clcPPE=>null(),execute_FFTW_r2r_z=>null()
  
  public:: InitPoissonSolver,clcPPE,Destory_Poisson_FFT_Plan
contains
#define nTime_FFT_Test 5

#define my_FFTW_inc_add_z2
#include "../Common/my_FFTW_inc.f90"
#undef  my_FFTW_inc_add_z2

!#define my_Poisson_inc_add_Periodic_2d
#include "../Common/my_Poisson_inc.f90"
!#undef  my_Poisson_inc_add_Periodic_2d
    
  !******************************************************************
  ! InitPoissonSolver
  !******************************************************************     
  subroutine InitPoissonSolver()
    implicit none
        
    ! locals
    real(RK)::wa1,wa3,best_time,t1(2),t2(2)
    integer::i,k,iErr01,iErr02,iChoice,iFFTz,IsReduce,IsReduceR
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm)::prphiHalo
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::prsrc
    
    prsrc=0.0_RK
    
    ! Modified wave number in x1 and x3 direct.
    allocate(WaveNumX(nxc), Stat =iErr01)
    allocate(WaveNumZ(nzc), Stat =iErr02)
    if(abs(iErr01)+abs(iErr02) /= 0) call MainLog%CheckForError(ErrT_Abort,"InitPoissonSolver","Allocation failed")    
    do i=1,nxc
      wa1= 2.0_RK*PI*real(i-1,RK)/real(nxc,RK)
      WaveNumX(i)=rdx2/288.0_RK*(cos(3.0_RK*wa1)-54.0_RK*cos(2.0_RK*wa1)+783.0_RK*cos(wa1)-730.0_RK)
    enddo
    do k=1,nzc
      wa3= 2.0_RK*PI*real(k-1,RK)/real(nzc,RK)
      WaveNumZ(k)=rdz2/288.0_RK*(cos(3.0_RK*wa3)-54.0_RK*cos(2.0_RK*wa3)+783.0_RK*cos(wa3)-730.0_RK)
    enddo
    normfft = 1.0_RK/(real(nxc,RK)*real(nzc,RK))
        
    best_time=max(huge(wa1),1.0E+20)
    if(nrank==0) call MainLog%OutInfo("Auto-tuning mode for Poisson Solver......",1)

      ! Choice-1
      call Create_Poisson_FFT_Plan(x1size,z2size)
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
      if(nrank==0)call MainLog%OutInfo("Choice-1, time= "//trim(num2str(t1(1)))//", "//trim(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=1; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=1; iFFTz=2
      endif
      call Destory_Poisson_FFT_Plan()
   
      ! Choice-2    
      call Create_Poisson_FFT_Plan(x2size,z1size)
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
      if(nrank==0)call MainLog%OutInfo("Choice-2, time= "//trim(num2str(t1(1)))//", "//trim(num2str(t1(2))),2)
      if(best_time>t1(1)) then
        best_time=t1(1); iChoice=2; iFFTz=1
      endif
      if(best_time>t1(2)) then
        best_time=t1(2); iChoice=2; iFFTz=2
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
        if(nrank==0)call MainLog%OutInfo("Choice-3, time= "//trim(num2str(t1(1)))//", "//trim(num2str(t1(2))),2)
        if(best_time>t1(1)) then
          best_time=t1(1); iChoice=3; iFFTz=1
        endif
        if(best_time>t1(2)) then
          best_time=t1(2); iChoice=3; iFFTz=2
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
        if(nrank==0)call MainLog%OutInfo("Choice-4, time= "//trim(num2str(t1(1)))//", "//trim(num2str(t1(2))),2)
        if(best_time>t1(1)) then
          best_time=t1(1); iChoice=4; iFFTz=1
        endif
        if(best_time>t1(2)) then
          best_time=t1(2); iChoice=4; iFFTz=2
        endif
        call decomp_info_finalize(decomp_PPE)
        deallocate(decomp_PPE,a_reduce,c_reduce)
      else
        if(nrank==0)call MainLog%OutInfo("Choice-4, x2size<4, ignore",2)
      endif
      call Destory_Poisson_FFT_Plan()

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
    if(iChoice==1) then
      clcPPE => clcPPE_x1
      call Create_Poisson_FFT_Plan(x1size,z2size)
      if(nrank==0) call MainLog%OutInfo("y1 -> x1 -> z2 -> y2 -> z2 -> x1 -> y1",3)
    elseif(iChoice==2) then
      clcPPE => clcPPE_z1
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
    if(nrank==0) print*," "
    prsrc=0.0_RK; prphiHalo=0.0_RK
  end subroutine InitPoissonSolver

  !******************************************************************
  ! Create_Poisson_FFT_Plan
  !******************************************************************
  subroutine Create_Poisson_FFT_Plan(xsizeIn,zsizeIn)
    implicit none
    integer,dimension(3),intent(in)::xsizeIn,zsizeIn
   
    ! locals
    integer::plan_type
    integer(C_FFTW_R2R_KIND)::kind_fwd,kind_bwd
    real(RK),dimension(:),allocatable::Vec1,Vec2

    if(FFTW_plan_type == 1) then
      plan_type=FFTW_MEASURE
    else
      plan_type=FFTW_ESTIMATE
    endif
        
    ! FFT in x
    kind_fwd = FFTW_R2HC
    kind_bwd = FFTW_HC2R
    allocate(Vec1(xsizeIn(1)),Vec2(xsizeIn(1)))
    fwd_plan_x= fftw_plan_r2r_1d(xsizeIn(1),Vec1,Vec2,kind_fwd,plan_type)
    bwd_plan_x= fftw_plan_r2r_1d(xsizeIn(1),Vec1,Vec2,kind_bwd,plan_type)
    deallocate(Vec1,Vec2)

    ! FFT in z
    kind_fwd = FFTW_R2HC
    kind_bwd = FFTW_HC2R
    allocate(Vec1(zsizeIn(3)),Vec2(zsizeIn(3)))
    fwd_plan_z= fftw_plan_r2r_1d(zsizeIn(3),Vec1,Vec2,kind_fwd,plan_type)
    bwd_plan_z= fftw_plan_r2r_1d(zsizeIn(3),Vec1,Vec2,kind_bwd,plan_type)
    deallocate(Vec1,Vec2)
  end subroutine Create_Poisson_FFT_Plan
  
#undef nTime_FFT_Test
end module m_Poisson
