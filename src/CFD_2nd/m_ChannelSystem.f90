module m_ChannelSystem
  use m_Timer
  use m_Tools
  use m_Typedef
  use m_LogInfo
  use m_Poisson
  use m_TScheme
  use m_FlowCase
  use m_Decomp2d
  use m_LESModels
  use m_IOAndVisu
  use m_Variables
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
    call MainLog%InitLog(ResultsDir, RunName,LF_file_lvl,LF_cmdw_lvl)
    if(nrank==0) call DumpReadedParam()

    call InitMeshAndMetries(ChannelPrm)
    call InitVisu(ChannelPrm)

    call AllocateVariables()
    call Init_Halo()
    call InitPoissonSolver(RealArr1,RealHalo)
    call InitTimeScheme(ChannelPrm)
    if(LES_type>0) call InitFilterCoe()

    if(.not. RestartFlag) then
      assoDevia: associate(Deviation=>RealArr1)
      call InitVelocity(ux,uy,uz,Deviation,ChannelPrm)
      end associate assoDevia
      if(nrank==0) call MainLog%OutInfo("Initializing all the needed variables into ChannelSystem ...", 1 )
    else
      call read_restart(ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
      if(nrank==0) call MainLog%OutInfo("Reading all the needed variables into ChannelSystem ...", 1 )
    endif
    call Update_uy_ym(uy_ym, duy_ym, zero)
    call InitStatVar(ChannelPrm) 
    
    SimTime=zero; dt = dtMax
    call SetBC_and_UpdateHalo(ux,uy,uz,uy_ym)
    call SetBC_and_UpdateHalo_pr( pressure )
    RealArr1=zero; RealArr2=zero
       
    ! Timers
    call total_timer%reset()
    call dump_visu(ifirst-1,ux,uy,uz,pressure,RealArr1)
   end subroutine  ChannelInitialize

  !******************************************************************
  ! ChannelIterate
  !******************************************************************
  subroutine ChannelIterate()
    implicit none

    !locals
    integer:: ns
    real(RK)::uddxmax,cflmp,divmax1,divmax2,uxm,vmaxabs(3)

    call total_timer%start()
    call CalcMaxCFL(ux,uy,uz,uddxmax)
    cflmp=uddxmax*dt
    if( icfl==1 ) then
      dt = CFLc/uddxmax
      dt = min(dt, dtMax)
    else
      dt = dtMax  
    endif
    if(LES_type>0) call Clc_SGS_Vis(ux,uy,uz,nut)

    do ns=1, iadvance
      ! step0: Update the Projection Method coefficients.
      call PMcoeUpdate(ns)
      call Update_uy_ym(uy_ym, duy_ym, SimTime)
      
      ! step1: Calculate the right hand side of the three velocity equations.
      asso_RHS123: associate( RhsX=>RealArr1, RhsY=>RealArr2,  RhsZ=>RealHalo)
      call clcRhsX(ux,uy,uz,RhsX,HistXOld,pressure)
      call clcRhsY(ux,uy,uz,RhsY,HistYOld,pressure)
      call clcRhsZ(ux,uy,uz,RhsZ,HistZOld,pressure)

      ! step2: Calculate the Uhat
      call clcOutFlowVelocity(ux,uy,uz)
      call clcU1Hat(ux,RhsX)
      call clcU2Hat(uy,RhsY,duy_ym)
      call clcU3Hat(uz,RhsZ)
      end associate asso_RHS123

      ! step3: Calculate the source term of the PPE 
      call SetBC_and_UpdateHaloForPrSrc( ux,uy,uz, uy_ym )
      call correctOutFlowFaceVelocity(ux,uy,uz)
      asso_Pr: associate(prsrc =>RealArr1, prphi =>RealArr2, prphiHalo =>RealHalo  )
      call clcPrSrc(ux,uy,uz,prsrc,pressure,divmax1)
      call clcPPE(prsrc,prphiHalo)
      call SetBC_and_UpdateHalo_pr(prphiHalo)
            
      ! step4: Update the velocity field to get the final real velocity.
      call FluidVelUpdate(prphiHalo,ux,uy,uz)
            
      ! step5: Update the real pressure field to get the final pressure.
      call PressureUpdate(pressure, prphiHalo)
      end associate asso_Pr
      call SetBC_and_UpdateHalo( ux,uy,uz,uy_ym )
      call SetBC_and_UpdateHalo_pr( pressure )
    enddo
    
    call CheckDivergence(ux,uy,uz, divmax2)  
    if(nrank==0 .and. divmax2>div_limit) call MainLog%CheckForError(ErrT_Abort,"ChannelIterate","too big div: "//trim(num2str(divmax2)))
    vmaxabs = CalcVmax(ux,uy,uz)
    if(nrank==0 .and. minval(vmaxabs)>vel_limit) call MainLog%CheckForError(ErrT_Abort,"ChannelIterate","too big velocity: "//trim(num2str(vmaxabs(1)))//", "//trim(num2str(vmaxabs(2)))//", "//trim(num2str(vmaxabs(3))) )

    if(mod(itime,ivstats)==0)    call clcStat(ux,uy,uz,pressure)
    if(mod(itime,SaveVisu)== 0)  call dump_visu(itime,ux,uy,uz,pressure,RealArr1)
    if(mod(itime,BackupFreq)== 0 .or. itime==ilast) then
      call Write_Restart(itime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
      call Delete_Prev_Restart(itime)
    endif
    call total_timer%finish()

    ! command window and log file output
    IF((itime==ifirst .or. mod(itime, Cmd_LFile_Freq)==0) ) THEN
      uxm = CalcUxAver(ux)
      if(nrank==0) then
        call MainLog%OutInfo("Channel3d performed "//trim(num2str(itime))//" iterations up to here!",1)
        call MainLog%OutInfo("Execution time [tot, last, ave] [sec]: "//trim(num2str(total_timer%tot_time))//", "// &
        trim(num2str(total_timer%last_time ))//", "//trim(num2str(total_timer%average())),2)
        call MainLog%OutInfo("SimTime | dt | CFL : "//trim(num2str(SimTime))//' | '//trim(num2str(dt))//' | '//trim(num2str(cflmp)),3)
        call MainLog%OutInfo("Max Abs Div: "//trim(num2str(divmax1))//" | "//trim(num2str(divmax2)) ,3)
        call MainLog%OutInfo("Max Abs Vel: "//trim(num2str(vmaxabs(1)))//" | "//trim(num2str(vmaxabs(2)))//" | "//trim(num2str(vmaxabs(3))), 3)
        call MainLog%OutInfo("Mean Velocity in streamwise: "//trim(num2str(uxm)),3)
      endif
    ENDIF

  end subroutine ChannelIterate

end module m_ChannelSystem
