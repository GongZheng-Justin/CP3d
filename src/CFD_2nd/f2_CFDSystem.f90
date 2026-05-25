#include "definitions_inc.f90"
module f2_CFDSystem
  use MPI
  use mc_Timer
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d
  use mc_FileOperator
  use f2_Tools
  use f2_Poisson
  use f2_TScheme
  use f2_FlowCase
  use f2_IOAndVisu
  use f2_Variables
  use f2_DumpPlane
  use f2_Stat_User
  use f2_Parameters
  use f2_BC_and_Halo
  use f2_MeshAndMetries
  implicit none
  private
    
  !// timers
  type(timer):: total_timer

  public:: CFDInitialize, CFDIterate
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CFDInitialize
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CFDInitialize(chFile); implicit none
    character(len=*),intent(in)::chFile
    
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
    
    call InitMeshAndMetries(chFile)
    call InitVisu(chFile)
    call Initialize_DumpPlane(chFile)

    call AllocateVariables()
    call InitPoissonSolver()
    call Init_BC_and_Halo()
    call InitTimeScheme()

    if(.not. RestartFlag) then
      assoDevia: associate(Deviation=>RealArr1)
      call InitVelocity(ux,uy,uz,Deviation,chFile)
      end associate assoDevia
      if(nrank==0) call MainLog%OutInfo("Initializing all the needed variables into CFDSystem ...", 1 )
    else
      call read_restart(ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
      if(nrank==0) call MainLog%OutInfo("Reading all the needed variables into CFDSystem ...", 1 )
    endif
    call InitStatVar(chFile)
    call InitStatVar_user(chFile)
    
    SimTime=0.0_RK; dt = dtMax
    call SetBC_and_UpdateHalo(ux,uy,uz)
    call SetBC_and_UpdateHalo_pr( pressure )
    RealArr1=0.0_RK; RealArr2=0.0_RK
       
    ! Timers
    call total_timer%reset()
    call dump_visu(ifirst-1,ux,uy,uz,pressure,RealArr1)
   end subroutine CFDInitialize

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CFDIterate
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CFDIterate(); implicit none

    !locals
    integer:: ns
    real(RK)::uddxmax,cflmp,divmax1,divmax2,uxm,vmaxabs(3)

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
      
      ! step1: Calculate the right hand side of the three velocity equations.
      asso_RHS123: associate(RhsX=>RealArr1, RhsY=>RealArr2,  RhsZ=>RealHalo)
      call clcRhsX(ux,uy,uz,RhsX,HistXOld,pressure)
      call clcRhsY(ux,uy,uz,RhsY,HistYOld,pressure)
      call clcRhsZ(ux,uy,uz,RhsZ,HistZOld,pressure)
      call add_FlowType_Forcing(RhsX,RhsY,RhsZ)     

      ! step2: update the out flow
      call clcOutFlowVelocity(ux,uy,uz)
      
      ! step3: Calculate the Uhat
      call clcU1Hat(ux,RhsX)
      call clcU2Hat(uy,RhsY)
      call clcU3Hat(uz,RhsZ)
      end associate asso_RHS123

      ! step4: Calculate the source term of the PPE 
      call SetBC_and_UpdateHaloForPrSrc( ux,uy,uz)
      call correctOutFlowFaceVelocity(ux,uy,uz)
      asso_Pr: associate(prsrc =>RealArr1, prphi =>RealArr2, prphiHalo =>RealHalo  )
      call clcPrSrc(ux,uy,uz,prsrc,pressure,divmax1)
      call clcPPE(prsrc,prphiHalo)
      call SetBC_and_UpdateHalo_pr(prphiHalo)
            
      ! step5: Update the velocity field to get the final real velocity.
      call FluidVelUpdate(prphiHalo,ux,uy,uz)
            
      ! step6: Update the real pressure field to get the final pressure.
      call PressureUpdate(pressure, prphiHalo)
      end associate asso_Pr
      call SetBC_and_UpdateHalo( ux,uy,uz)
      call SetBC_and_UpdateHalo_pr( pressure )
      !if(nrank==0) print*, itime,'PrGradData=', PrGradData
    enddo
    if(mod(itime,BackupFreq)== 0 .or. itime==ilast) then
      call Write_Restart(itime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
      call Delete_Prev_Restart(itime)
    endif

    if(mod(itime,ivstats)==0) then
      call clcStat(ux,uy,uz,pressure,RealArr1,RealArr2)
      call clcStat_User(ux,uy,uz,pressure)
    endif
    call dump_plane(itime,ux,uy,uz,pressure)
    if(mod(itime,SaveVisu)== 0)  call dump_visu(itime,ux,uy,uz,pressure,RealArr1)
    call total_timer%finish()

    ! command window and log file output
    IF((itime==ifirst .or. mod(itime, Cmd_LFile_Freq)==0) ) THEN
      call CheckDivergence(ux,uy,uz, divmax2)
      if(nrank==0 .and. divmax2>div_limit) call MainLog%CheckForError(ErrT_Abort,"CFDIterate","too big div: "//strip(num2str(divmax2)))
      vmaxabs = CalcVmax(ux,uy,uz)
      if(nrank==0 .and. minval(vmaxabs)>vel_limit) call MainLog%CheckForError(ErrT_Abort,"CFDIterate","too big velocity: "//strip(num2str(vmaxabs(1)))//", "//strip(num2str(vmaxabs(2)))//", "//strip(num2str(vmaxabs(3))) )

      uxm = CalcUxAver(ux)
      if(nrank==0) then
        call MainLog%OutInfo("CFD_2nd performed "//strip(num2str(itime))//" iterations up to here!",1)
        call MainLog%OutInfo("Execution time [tot, last, ave] [sec]: "//strip(num2str(total_timer%tot_time))//", "// &
        strip(num2str(total_timer%last_time ))//", "//strip(num2str(total_timer%average())),2)
        call MainLog%OutInfo("SimTime | dt | CFL : "//strip(num2str(SimTime))//' | '//strip(num2str(dt))//' | '//strip(num2str(cflmp)),3)
        call MainLog%OutInfo("Max Abs Div: "//strip(num2str(divmax1))//" | "//strip(num2str(divmax2)) ,3)
        call MainLog%OutInfo("Max Abs Vel: "//strip(num2str(vmaxabs(1)))//" | "//strip(num2str(vmaxabs(2)))//" | "//strip(num2str(vmaxabs(3))), 3)
        call MainLog%OutInfo("Mean Velocity in streamwise: "//strip(num2str(uxm)),3)
      endif
    ENDIF

  end subroutine CFDIterate

end module f2_CFDSystem
