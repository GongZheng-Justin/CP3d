#include "definitions_inc.f90"
module f4_CFDSystem
  use MPI
  use mc_Timer
  use mc_TypeDef
  use mc_LogInfo
  use mc_decomp2d
  use mc_FileOperator,only:mkdir
  use f4_Tools
  use f4_TScheme
  use f4_Poisson
  use f4_FlowCase
  use f4_Variables
  use f4_IOAndVisu
  use f4_Parameters
  use f4_BC_and_Halo
  use f4_MeshAndMetries
  implicit none
  private
    
  !// timers
  type(timer):: total_timer
  public:: CFDInitialize, CFDIterate
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CFDInitialize
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CFDInitialize(CFDPrm); implicit none
    character(len=*),intent(in)::CFDPrm
    
    ! locals
    integer::ierr
   
    ! Initializing main log info
    if(nrank==0) then
      call mkdir(ResultsDir_, ierr); if(ierr<0) then; print*,'Cannot create folder ResultsDir'; stop; endif
      call mkdir(RestartDir_, ierr); if(ierr<0) then; print*,'Cannot create folder RestartDir'; stop; endif
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MainLog%InitLog(ResultsDir_, RunName_, LF_file_lvl,LF_cmdw_lvl)
    if(nrank==0) call MainLog%CreateFile(RunName_)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MainLog%OpenFile()    
    if(nrank==0) call DumpReadedParam()

    call InitMeshAndMetries(CFDPrm)
    call InitVisu(CFDPrm)

    call AllocateVariables()
    call InitPoissonSolver()
    call Init_BC_and_Halo()
    call InitTimeScheme()
    call InitStatVar(CFDPrm)
    
    if(.not. RestartFlag) then
      assoDevia: associate(Deviation=>RealArr1)
#ifdef ScalarFlow
      call InitFlowField(ux,uy,uz,Deviation,scalar)
#else
      call InitFlowField(ux,uy,uz,Deviation)
#endif
      end associate assoDevia
      if(nrank==0) call MainLog%OutInfo("Initializing all the needed variables ...", 1 )
    else
#ifdef ScalarFlow
      call read_restart(ux,uy,uz,pressure,scalar,HistXOld,HistYOld,HistZOld,HistCOld)
#else
      call read_restart(ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
#endif
      if(nrank==0) call MainLog%OutInfo("Reading all the needed variables ...", 1 )
    endif
#ifdef ScalarFlow
      if(nrank==0) call MainLog%OutInfo("Scalar ON  ...", 1 )
#else
      if(nrank==0) call MainLog%OutInfo("Scalar OFF ...", 1 )
#endif
    call Update_uy_ym(uy_ym, duy_ym)
    SimTime=0.0_RK; dt = dtMax
    call SetBC_and_UpdateHalo(ux,uy,uz,uy_ym)
    call SetBC_and_UpdateHalo_pr( pressure )
#ifdef ScalarFlow
    call SetBC_and_UpdateHalo_scalar(scalar)
    call dump_visu(ifirst-1,ux,uy,uz,pressure,scalar,RealArr1)
#else
    call dump_visu(ifirst-1,ux,uy,uz,pressure,RealArr1)
#endif
    RealArr1=0.0_RK; RealArr2=0.0_RK
        
    ! Timers
    call total_timer%reset()
   end subroutine CFDInitialize

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CFDIterate
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CFDIterate(); implicit none

    !locals
    integer:: ns
    real(RK)::uddxmax,cflmp,divmax1,divmax2,vmaxabs(3)
#ifdef ScalarFlow
    real(RK)::MeanValue(3)
#else
    real(RK)::UxMean
#endif

    call total_timer%start()
    call CalcMaxCFL(ux,uy,uz,uddxmax)
    if( icfl==1 ) then
      dt = CFLc/uddxmax
      dt = min(dt, dtMax)
    else
      dt = dtMax  
    endif
    cflmp=uddxmax*dt

    do ns=1, iadvance
      ! step0: Update the Projection Method coefficients.
      call PMcoeUpdate(ns)
      call Update_uy_ym(uy_ym, duy_ym)

      ! step1: Calculate the right hand side of the three velocity equations
      asso_RHS123: associate( RhsX=>RealArr1, RhsY=>RealArr2,  RhsZ=>RealHalo)
#ifdef ScalarFlow
      call clcRhsXYZ(ux,uy,uz,scalar,RhsX,RhsY,RhsZ,RealArrc,HistXOld,HistYOld,HistZOld,HistCOld,pressure)
      call clcScalar(scalar,RealArrc)
#else
      call clcRhsXYZ(ux,uy,uz,RhsX,RhsY,RhsZ,HistXOld,HistYOld,HistZOld,pressure)
#endif

      ! step2: Calculate the Uhat
      call clcU1Hat(ux,RhsX)
      call clcU2Hat(uy,RhsY,duy_ym)
      call clcU3Hat(uz,RhsZ)
      end associate asso_RHS123

      ! step3: Calculate the source term of the PPE 
      call SetBC_and_UpdateHaloForPrSrc( ux,uy,uz, uy_ym )
      asso_Pr: associate(prsrc =>RealArr1, prphi =>RealArr2, prphiHalo =>RealHalo  )
      call clcPrSrc(ux,uy,uz,prsrc,divmax1)
      call clcPPE(prsrc,prphiHalo)
      call SetBC_and_UpdateHalo_pr(prphiHalo)
            
      ! step4: Update the velocity field to get the final real velocity.
      call FluidVelUpdate(prphiHalo,ux,uy,uz)
            
      ! step5: Update the real pressure field to get the final pressure.
      call PressureUpdate(pressure, prphiHalo)
      end associate asso_Pr
      call SetBC_and_UpdateHalo( ux,uy,uz,uy_ym )
      call SetBC_and_UpdateHalo_pr( pressure )
#ifdef ScalarFlow
      call SetBC_and_UpdateHalo_scalar(scalar)
#endif
    enddo
#ifdef ScalarFlow
    if(mod(itime,ivstats)==0)   call clcStat(ux,uy,uz,pressure,scalar,RealArr1,RealArr2)
    if(mod(itime,SaveVisu)== 0) call dump_visu(itime,ux,uy,uz,pressure,scalar,RealArr1)
    if(mod(itime,BackupFreq)== 0 .or. itime==ilast) then
      call Write_Restart(itime,ux,uy,uz,pressure,scalar,HistXOld,HistYOld,HistZOld,HistCOld)
      call Delete_Prev_Restart(itime)
    endif
#else
    if(mod(itime,ivstats)==0)   call clcStat(ux,uy,uz,pressure,RealArr1,RealArr2)
    if(mod(itime,SaveVisu)== 0) call dump_visu(itime,ux,uy,uz,pressure,RealArr1)
    if(mod(itime,BackupFreq)== 0 .or. itime==ilast) then
      call Write_Restart(itime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
      call Delete_Prev_Restart(itime)
    endif
#endif
    call total_timer%finish()

    ! command window and log file output
    IF((itime==ifirst .or. mod(itime, Cmd_LFile_Freq)==0) ) THEN
      call CheckDivergence(ux,uy,uz, divmax2)
      if(nrank==0 .and. divmax2>div_limit) call MainLog%CheckForError(ErrT_Abort,"CFDIterate","too big div: "//strip(num2str(divmax2)))
      vmaxabs = CalcVmax(ux,uy,uz)
      if(nrank==0 .and. minval(vmaxabs)>vel_limit) call MainLog%CheckForError(ErrT_Abort,"CFDIterate","too big velocity: "//strip(num2str(vmaxabs(1)))//", "//strip(num2str(vmaxabs(2)))//", "//strip(num2str(vmaxabs(3))) )
#ifdef ScalarFlow
      MeanValue = ClcMeanValue(ux,scalar)
      !MeanValue(3)=MeanValue(3)/(Scalar_InitValue*Scalar_InitValue)-1.0_RK
#else
      UxMean = CalcUxAver(ux)
#endif
      if(nrank==0) then
        call MainLog%OutInfo("CFD_4th performed "//strip(num2str(itime))//" iterations up to here!",1)
        call MainLog%OutInfo("Execution time [tot, last, ave] [sec]: "//strip(num2str(total_timer%tot_time))//", "// &
        strip(num2str(total_timer%last_time ))//", "//strip(num2str(total_timer%average())),2)
        call MainLog%OutInfo("SimTime | dt | CFL : "//strip(num2str(SimTime))//' | '//strip(num2str(dt))//' | '//strip(num2str(cflmp)),3)
        call MainLog%OutInfo("Max Abs Div: "//strip(num2str(divmax1))//" | "//strip(num2str(divmax2)) ,3)
        call MainLog%OutInfo("Max Abs Vel: "//strip(num2str(vmaxabs(1)))//" | "//strip(num2str(vmaxabs(2)))//" | "//strip(num2str(vmaxabs(3))), 3)
#ifdef ScalarFlow
        call MainLog%OutInfo("Mean streamwise velocity | concentration | c^2: "//strip(num2str(MeanValue(1)))//' | '//strip(num2str(MeanValue(2)))//' | '//strip(num2str(MeanValue(3))),3)
#else
        call MainLog%OutInfo("Mean streamwise velocity: "//strip(num2str(UxMean)),3)
#endif
      endif
    ENDIF

  end subroutine CFDIterate

end module f4_CFDSystem
