module m_ChannelSystem
  use m_Tools
  use m_Timer
  use m_TScheme
  use m_Typedef
  use m_LogInfo
  use m_Poisson
  use m_FlowCase
  use m_decomp2d
  use m_Variables
  use m_IOAndVisu
  use m_Parameters
  use m_BC_and_Halo
  use m_MeshAndMetries
  implicit none
  private
    
  !// timers
  type(timer):: total_timer
  public:: ChannelInitialize, ChannelIterate
contains

  !******************************************************************
  ! ChannelInitialize
  !******************************************************************
  subroutine ChannelInitialize(ChannelPrm)
    implicit none
    character(*),intent(in)::ChannelPrm
    character(256):: chStr
   
    !// Initializing main log info
    write(chStr,"(A)") 'mkdir -p '//ResultsDir//' '//RestartDir//' 2> /dev/null'
    if(nrank==0) call system(trim(adjustl(chStr)))
    call MainLog%InitLog(ResultsDir,RunName,LF_file_lvl,LF_cmdw_lvl)
    if(nrank==0) call DumpReadedParam()

    call InitMeshAndMetries(ChannelPrm)
    call InitVisu(ChannelPrm)

    call AllocateVariables() 
    call Init_Halo()
    call InitPoissonSolver(RealArr1,RealHalo)
    call InitTimeScheme()
    call InitStatVar(ChannelPrm)
    
    if(.not. RestartFlag) then
      assoDevia: associate(Deviation=>RealArr1)
      call InitVelocity(ux,uy,uz,Deviation)
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
    SimTime=zero; dt = dtMax
    call SetBC_and_UpdateHalo(ux,uy,uz,uy_ym)
    call SetBC_and_UpdateHalo_pr( pressure )
#ifdef ScalarFlow
    call SetBC_and_UpdateHalo_scalar(scalar)
    call dump_visu(ifirst-1,ux,uy,uz,pressure,scalar,RealArr1)
#else
    call dump_visu(ifirst-1,ux,uy,uz,pressure,RealArr1)
#endif
    RealArr1=zero; RealArr2=zero
        
    ! Timers
    call total_timer%reset()
   end subroutine  ChannelInitialize

  !******************************************************************
  ! ChannelIterate
  !******************************************************************
  subroutine ChannelIterate()
    implicit none

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
    cflmp=uddxmax*dt
    if( icfl==1 ) then
      dt = CFLc/uddxmax
      dt = min(dt, dtMax)
    else
      dt = dtMax  
    endif

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
    
    call CheckDivergence(ux,uy,uz, divmax2)
    if(nrank==0 .and. divmax2>div_limit) call MainLog%CheckForError(ErrT_Abort,"ChannelIterate","too big div: "//trim(num2str(divmax2)))
    vmaxabs = CalcVmax(ux,uy,uz)
    if(nrank==0 .and. minval(vmaxabs)>vel_limit) call MainLog%CheckForError(ErrT_Abort,"ChannelIterate","too big velocity: "//trim(num2str(vmaxabs(1)))//", "//trim(num2str(vmaxabs(2)))//", "//trim(num2str(vmaxabs(3))) )

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
#ifdef ScalarFlow
      MeanValue = ClcMeanValue(ux,scalar)
      MeanValue(3)=MeanValue(3)/(Scalar_InitValue*Scalar_InitValue)-one
#else
      UxMean = CalcUxAver(ux)
#endif
      if(nrank==0) then
        call MainLog%OutInfo("Channel3d_4th performed "//trim(num2str(itime))//" iterations up to here!",1)
        call MainLog%OutInfo("Execution time [tot, last, ave] [sec]: "//trim(num2str(total_timer%tot_time))//", "// &
        trim(num2str(total_timer%last_time ))//", "//trim(num2str(total_timer%average())),2)
        call MainLog%OutInfo("SimTime | dt | CFL : "//trim(num2str(SimTime))//' | '//trim(num2str(dt))//' | '//trim(num2str(cflmp)),3)
        call MainLog%OutInfo("Max Abs Div: "//trim(num2str(divmax1))//" | "//trim(num2str(divmax2)) ,3)
        call MainLog%OutInfo("Max Abs Vel: "//trim(num2str(vmaxabs(1)))//" | "//trim(num2str(vmaxabs(2)))//" | "//trim(num2str(vmaxabs(3))), 3)
#ifdef ScalarFlow
        call MainLog%OutInfo("Mean streamwise velocity | concentration | c^2: "//trim(num2str(MeanValue(1)))//' | '//trim(num2str(MeanValue(2)))//' | '//trim(num2str(MeanValue(3))),3)
#else
        call MainLog%OutInfo("Mean streamwise velocity: "//trim(num2str(UxMean)),3)
#endif
      endif
    ENDIF

  end subroutine ChannelIterate

end module m_ChannelSystem
