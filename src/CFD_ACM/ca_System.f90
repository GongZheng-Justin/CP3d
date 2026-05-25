#include "definitions_inc.f90"
module ca_system
  use MPI
  use mc_Timer
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d
  use f2_Tools
  use f2_TScheme
  use f2_Poisson
  use f2_FlowCase
  use f2_Variables
  use f2_IOAndVisu
  use f2_DumpPlane
  use f2_Stat_User
  use f2_Parameters
  use f2_BC_and_Halo  
  use sp_Comm
  use sp_System
  use sp_Property
  use sp_variables
  use sp_Parameters
  use ca_Statistics
  use ca_IBM
  use ca_IBM_implicit
  use ca_BC_and_Halo,only: SetBC_and_UpdateHalo_VelIBM
  implicit none
  private
 
  !// timers
  type(timer):: total_timer
  public::CFDACM_Iterate, CFDACM_Iterate_imp
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CFDACM_Iterate
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CFDACM_Iterate(); implicit none

    ! locals
    real(RK)::uddxmax,cflmp,divmax1,divmax2,uxm,vmaxabs(3)
    integer::ns,iForcingExtra,idem_start,idem_end,idem,nForcingExtra_ns(3)

#ifdef SeveralSphereInfo
    character(len=:),allocatable::chFile
    if(DEM_Opt%numPrtcl==1 .and. GPrtcl_list%nlocal==1) then
      chFile = strip(ResultsDir_) // "SphereInfo_" // strip(RunName_) // "_" // int2str(GPrtcl_id(1),2) // ".txt"
      open(newunit=idem,file=chFile,status='replace',form='formatted')
      write(idem,'(A)')' t, x, y, z, r, vx, vy, vz, wx, wy, wz fx fy fz'
      close(idem)    
    else
      do ns=1,GPrtcl_list%nlocal
        chFile = strip(ResultsDir_) // "SphereInfo_" // strip(RunName_) // "_" // int2str(GPrtcl_id(ns),2) // ".txt"
        open(newunit=idem,file=chFile,status='replace',form='formatted')
        write(idem,'(A)')' t, x, y, z, r, vx, vy, vz, wx, wy, wz'
        close(idem)
      enddo    
    endif
#endif
#ifdef ChanBraunJFM2011
    integer::nfstime,iUnit,ierr
    character(len=:),allocatable::filename    
    real(RK),dimension(8,3)::ForceAndTorque
    type(real3),dimension(:),allocatable::Ave_FpForce
    type(real3),dimension(:),allocatable::Ave_FpTorque
    
    nfstime=0
    ForceAndTorque=0.0_RK
    if(GPrtcl_list%nlocal>0) then
      allocate(Ave_FpForce(1:max(GPrtcl_list%nlocal,1)))
      allocate(Ave_FpTorque(1:max(GPrtcl_list%nlocal,1)))
    endif

    if(nrank==0) then
      filename = strip(DEM_opt%ResultsDir) // "FpForce" // int2str(ilast,10) // '.txt'
      open(newunit=iUnit,file=filename,status='replace',form='formatted',IOSTAT=ierr)
      if(ierr /= 0) call MainLog%CheckForError(ErrT_Abort,"CFDACM_Iterate","Cannot open file: "// filename)
      close(iUnit,IOSTAT=ierr)
    endif
#endif

    nForcingExtra_ns = nForcingExtra
    !if(ischeme==FI_RK3 .and. nForcingExtra>1) nForcingExtra_ns(2) = nForcingExtra-1
    if(.not. DEM_Opt%RestartFlag) then
      dt= dtMax
      call PMcoeUpdate(iadvance)
      SimTime=0.0_RK
      associate_VolForce: associate(VolForce_x=>IBMArr1,VolForce_y=>IBMArr2,VolForce_z=>IBMArr3)
      call PrepareIBM_interp(VolForce_x,VolForce_y,VolForce_z)
      end associate associate_VolForce
      call IntegrateFluidPrtclForce(ux,uy,uz) 
    endif
    DO itime=ifirst, ilast
      call total_timer%start()
      call CalcMaxCFL(ux,uy,uz,uddxmax)
      if( icfl==1 ) then
        dt= CFLc/uddxmax
        dt= min(dt,dtMax)
        DEM_Opt%dt= dt/real(icouple,RK)
      else
        dt= dtMax  
      endif
      cflmp=uddxmax*dt

      do ns=1,iadvance
        ! step0: Update the Projection Method coefficients.
        call PMcoeUpdate(ns)
        idem_start=(itime-1)*icouple+ idem_advance_start(ns)
        idem_end  =(itime-1)*icouple+ idem_advance_end(ns)

        ! DEM Part for explicit coupling =========================
        do idem=idem_start,idem_end
          call DEM%iterate(idem)
        enddo
        
        ! step1: Calculate the right hand side of the three velocity equations.
        asso_RHS123: associate( RhsX=>RealArr1, RhsY=>RealArr2, RhsZ=>RealHalo)
        call clcRhsX(ux,uy,uz,RhsX,HistXOld,pressure)
        call clcRhsY(ux,uy,uz,RhsY,HistYOld,pressure)
        call clcRhsZ(ux,uy,uz,RhsZ,HistZOld,pressure)

        ! IBM Part1: compute the predicted velocity without the effect of boundary
!#define Double_FIMP_IBM
#ifdef Double_FIMP_IBM
        uxStar=ux; uyStar=uy; uzStar=uz
        call clcU1Hat(uxStar,RhsX)
        call clcU2Hat(uyStar,RhsY)
        call clcU3Hat(uzStar,RhsZ)
#else
        call clc_uStar(uxStar,uyStar,uzStar,ux,uy,uz,RhsX,RhsY,RhsZ)
#endif
    
        ! step2: update the out flow
        call clcOutFlowVelocity(ux,uy,uz)

        ! IBM Part2: interpolate and spread the IBM force
        associate_VolForce: associate(VolForce_x=>IBMArr1,VolForce_y=>IBMArr2,VolForce_z=>IBMArr3)
        call PrepareIBM_interp(VolForce_x,VolForce_y,VolForce_z)
        do iForcingExtra=0, nForcingExtra_ns(ns)
          call SetBC_and_UpdateHalo_VelIBM(uxStar,uyStar,uzStar)
          call AdditionalForceIBM(uxStar,uyStar,uzStar,VolForce_x,VolForce_y,VolForce_z)
        enddo
        end associate associate_VolForce
#ifdef NoSlipTest
        call SetBC_and_UpdateHalo_VelIBM(uxStar,uyStar,uzStar)
        call Clc_NoSlipErr(uxStar,uyStar,uzStar)
#endif

        ! step3: Calculate the Uhat
        if(IsImplicit==0) then
          call clcU123Hat(uxStar,uyStar,uzStar,ux,uy,uz)
        else
          call updateRhsIBM(uxStar,uyStar,uzStar,ux,uy,uz,RhsX,RhsY,RhsZ)
          call clcU1Hat(ux,RhsX)
          call clcU2Hat(uy,RhsY)
          call clcU3Hat(uz,RhsZ)        
        endif
        end associate asso_RHS123

        ! step4: Calculate the source term of the PPE 
        call SetBC_and_UpdateHaloForPrSrc( ux,uy,uz)
        call correctOutFlowFaceVelocity(ux,uy,uz)
        asso_Pr: associate(prsrc =>RealArr1, prphi =>RealArr2, prphiHalo =>RealHalo)
        call clcPrSrc(ux,uy,uz,prsrc,pressure,divmax1)
        call clcPPE(prsrc,prphiHalo)
        call SetBC_and_UpdateHalo_pr(prphiHalo)
            
        ! step5: Update the velocity field to get the final real velocity.
        call FluidVelUpdate(prphiHalo,ux,uy,uz)
            
        ! step6: Update the real pressure field  to get the final pressure.
        call PressureUpdate(pressure, prphiHalo)
        end associate asso_Pr
        call SetBC_and_UpdateHalo(ux,uy,uz)
        call SetBC_and_UpdateHalo_pr( pressure )
        
        ! IBM Part3: compute the IBM force for every single particle
        call IntegrateFluidPrtclForce(ux,uy,uz)
        !if(nrank==0) print*, itime,'PrGradData=', PrGradData
        
#ifdef ChanBraunJFM2011
#include "ca_Force_inc.f90"
#endif
#ifdef SeveralSphereInfo
        if(DEM_Opt%numPrtcl==1 .and. GPrtcl_list%nlocal==1) then
          chFile = strip(ResultsDir_) //"SphereInfo_" // strip(RunName_) // "_" // int2str(GPrtcl_id(1),2) // ".txt"
          open(newunit=idem,file=chFile,status='old',position='append',form='formatted')
          write(idem,'(15ES24.15)')SimTime,GPrtcl_PosR(1),GPrtcl_linVel(1,1),GPrtcl_rotVel(1,1),GPrtcl_FpForce(1)
          close(idem) 
        else
          do iForcingExtra=1,GPrtcl_list%nlocal
            chFile = strip(ResultsDir_) //"SphereInfo_" // strip(RunName_) // "_" // int2str(GPrtcl_id(iForcingExtra),2) // ".txt"
            open(newunit=idem,file=chFile,status='old',position='append',form='formatted')
            write(idem,'(15ES24.15)')SimTime,GPrtcl_PosR(iForcingExtra),GPrtcl_linVel(1,iForcingExtra),GPrtcl_rotVel(1,iForcingExtra)
            close(idem)
          enddo        
        endif
#endif
      enddo
      if(mod(itime,ivstats)==0) then
        call Update_FluidIndicator(FluidIndicator) 
        call clcStat(ux,uy,uz,pressure,RealArr1,RealArr2)
        call clcStat_User(ux,uy,uz,pressure)
        call ClcCAStatistics()
        call dump_plane(itime,ux,uy,uz,pressure)   
      endif
      if(mod(itime,SaveVisu)== 0) call dump_visu(itime,ux,uy,uz,pressure,RealArr1)
      if(mod(itime,BackupFreq)== 0 .or. itime==ilast) then
        call Write_Restart(itime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
        call Delete_Prev_Restart(itime)
        call DEM%Write_Restart(itime*icouple)
      endif
      call total_timer%finish()

      ! command window and log file output
      IF((itime==ifirst .or. mod(itime, Cmd_LFile_Freq)==0) ) THEN
        call CheckDivergence(ux,uy,uz, divmax2)  
        if(nrank==0 .and. divmax2>div_limit) call MainLog%CheckForError(ErrT_Abort,"CFDACM_Iterate","too big div: "//strip(num2str(divmax2)))
        vmaxabs = CalcVmax(ux,uy,uz)
        if(nrank==0 .and. minval(vmaxabs)>vel_limit) call MainLog%CheckForError(ErrT_Abort,"CFDACM_Iterate","too big velocity: "//strip(num2str(vmaxabs(1)))//", "//strip(num2str(vmaxabs(2)))//", "//strip(num2str(vmaxabs(3))) )
        uxm = CalcUxAver(ux)
        if(nrank==0) then
          call MainLog%OutInfo("CFD_ACM performed "//strip(num2str(itime))//" iterations up to here!",1)
          call MainLog%OutInfo("Execution time [tot, last, ave] [sec]: "//strip(num2str(total_timer%tot_time))//", "// &
          strip(num2str(total_timer%last_time ))//", "//strip(num2str(total_timer%average())),2)
          call MainLog%OutInfo("SimTime | dt | CFL : "//strip(num2str(SimTime))//' | '//strip(num2str(dt))//' | '//strip(num2str(cflmp)),3)
          call MainLog%OutInfo("Max Abs Div: "//strip(num2str(divmax1))//" | "//strip(num2str(divmax2)) ,3)
          call MainLog%OutInfo("Max Abs Vel: "//strip(num2str(vmaxabs(1)))//" | "//strip(num2str(vmaxabs(2)))//" | "//strip(num2str(vmaxabs(3))), 3)
          call MainLog%OutInfo("Mean Velocity in streamwise: "//strip(num2str(uxm)),3)
        endif
      ENDIF
    ENDDO
  end subroutine CFDACM_Iterate

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CFDACM_Iterate_imp
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CFDACM_Iterate_imp(); implicit none

    ! locals
    integer::ns,idem_start,idem_end,idem
    real(RK)::uddxmax,cflmp,divmax1,divmax2,uxm,vmaxabs(3)

#ifdef SeveralSphereInfo
    character(len=:),allocatable::chFile
    if(DEM_Opt%numPrtcl==1 .and. GPrtcl_list%nlocal==1) then
      chFile = strip(ResultsDir_) // "SphereInfo_" // strip(RunName_) // "_" // int2str(GPrtcl_id(1),2) // ".txt"
      open(newunit=idem,file=chFile,status='replace',form='formatted')
      write(idem,'(A)')' t, x, y, z, r, vx, vy, vz, wx, wy, wz fx fy fz'
      close(idem)    
    else
      do ns=1,GPrtcl_list%nlocal
        chFile = strip(ResultsDir_) // "SphereInfo_" // strip(RunName_) // "_" // int2str(GPrtcl_id(ns),2) // ".txt"
        open(newunit=idem,file=chFile,status='replace',form='formatted')
        write(idem,'(A)')' t, x, y, z, r, vx, vy, vz, wx, wy, wz'
        close(idem)
      enddo    
    endif
#endif
#ifdef ChanBraunJFM2011
    integer::nfstime,iUnit,ierr
    character(len=:),allocatable::filename    
    real(RK),dimension(8,3)::ForceAndTorque
    type(real3),dimension(:),allocatable::Ave_FpForce
    type(real3),dimension(:),allocatable::Ave_FpTorque
    
    nfstime=0
    ForceAndTorque=0.0_RK
    if(GPrtcl_list%nlocal>0) then
      allocate(Ave_FpForce(1:max(GPrtcl_list%nlocal,1)))
      allocate(Ave_FpTorque(1:max(GPrtcl_list%nlocal,1)))
    endif

    if(nrank==0) then
      filename = strip(DEM_opt%ResultsDir) // "FpForce" // int2str(ilast,10) // '.txt'
      open(newunit=iUnit,file=filename,status='replace',form='formatted',IOSTAT=ierr)
      if(ierr /= 0) call MainLog%CheckForError(ErrT_Abort,"CFDACM_Iterate","Cannot open file: "// filename)
      close(iUnit,IOSTAT=ierr)
    endif
#endif

    if(.not. DEM_Opt%RestartFlag) then
      dt= dtMax
      call PMcoeUpdate(iadvance)
      SimTime=0.0_RK
      call PrepareIBM_interp_imp1()
      call IntegrateFluidPrtclForce_imp(ux,uy,uz)
    endif
    DO itime=ifirst, ilast
      call total_timer%start()
      call CalcMaxCFL(ux,uy,uz,uddxmax)
      if( icfl==1 ) then
        dt= CFLc/uddxmax
        dt= min(dt,dtMax)
        DEM_Opt%dt= dt/real(icouple,RK)
      else
        dt= dtMax  
      endif
      cflmp=uddxmax*dt

      do ns=1,iadvance
        ! step0: Update the Projection Method coefficients.
        call PMcoeUpdate(ns)
        idem_start=(itime-1)*icouple+ idem_advance_start(ns)
        idem_end  =(itime-1)*icouple+ idem_advance_end(ns)
        
        ! step1: Calculate the right hand side of the three velocity equations.
        asso_RHS123: associate( RhsX=>RealArr1, RhsY=>RealArr2, RhsZ=>RealHalo)
        call clcRhsX(ux,uy,uz,RhsX,HistXOld,pressure)
        call clcRhsY(ux,uy,uz,RhsY,HistYOld,pressure)
        call clcRhsZ(ux,uy,uz,RhsZ,HistZOld,pressure)

        ! IBM Part1 for semi-implicit coupling =================
        call clc_uStar(uxStar,uyStar,uzStar,ux,uy,uz,RhsX,RhsY,RhsZ)
        call SetBC_and_UpdateHalo_VelIBM( uxStar,uyStar,uzStar)
        call PrepareIBM_interp_imp1()
        call InterpolationIBM_imp(uxStar,uyStar,uzStar)
        call IntegrateFluidPrtclForce_imp(uxStar,uyStar,uzStar)
        associate_VolForce_imp: associate(VolForce_x=>IBMArr1,VolForce_y=>IBMArr2,VolForce_z=>IBMArr3)
        call SpreadIbpForce_imp(VolForce_x,VolForce_y,VolForce_z,1)
        call UpdateRhsIBM_imp1(RhsX,RhsY,RhsZ,VolForce_x,VolForce_y,VolForce_z)

        ! DEM Part for semi-implicit coupling ==================
        do idem=idem_start,idem_end
          call DEM%iterate(idem)
        enddo

        ! IBM Part2 for semi-implicit coupling =================
        call PrepareIBM_interp_imp2()
        call SpreadIbpForce_imp(VolForce_x,VolForce_y,VolForce_z,2)
        call UpdateRhsIBM_imp2(RhsX,RhsY,RhsZ,VolForce_x,VolForce_y,VolForce_z)
        end associate associate_VolForce_imp
          
        ! step2: Calculate the Uhat
        call clcOutFlowVelocity(ux,uy,uz)
        call clcU1Hat(ux,RhsX)
        call clcU2Hat(uy,RhsY)
        call clcU3Hat(uz,RhsZ)
        end associate asso_RHS123

        ! step3: Calculate the source term of the PPE 
        call SetBC_and_UpdateHaloForPrSrc( ux,uy,uz)
        call correctOutFlowFaceVelocity(ux,uy,uz)
        asso_Pr: associate(prsrc =>RealArr1, prphi =>RealArr2, prphiHalo =>RealHalo)
        call clcPrSrc(ux,uy,uz,prsrc,pressure,divmax1)
        call clcPPE(prsrc,prphiHalo)
        call SetBC_and_UpdateHalo_pr(prphiHalo)
            
        ! step4: Update the velocity field to get the final real velocity.
        call FluidVelUpdate(prphiHalo,ux,uy,uz)
            
        ! step5: Update the real pressure field  to get the final pressure.
        call PressureUpdate(pressure, prphiHalo)
        end associate asso_Pr
        call SetBC_and_UpdateHalo( ux,uy,uz)
        call SetBC_and_UpdateHalo_pr( pressure )
        
#ifdef ChanBraunJFM2011
#include "ca_Force_inc.f90"
#endif
#ifdef SeveralSphereInfo
        if(DEM_Opt%numPrtcl==1 .and. GPrtcl_list%nlocal==1) then
          chFile = strip(ResultsDir_) //"SphereInfo_" // strip(RunName_) // "_" // int2str(GPrtcl_id(1),2) // ".txt"
          open(newunit=idem,file=chFile,status='old',position='append',form='formatted')
          write(idem,'(15ES24.15)')SimTime,GPrtcl_PosR(1),GPrtcl_linVel(1,1),GPrtcl_rotVel(1,1),GPrtcl_FpForce(1)
          close(idem) 
        else
          do idem_end=1,GPrtcl_list%nlocal
            chFile = strip(ResultsDir_) //"SphereInfo_" // strip(RunName_) // "_" // int2str(GPrtcl_id(idem_end),2) // ".txt"
            open(newunit=idem,file=chFile,status='old',position='append',form='formatted')
            write(idem,'(15ES24.15)')SimTime,GPrtcl_PosR(idem_end),GPrtcl_linVel(1,idem_end),GPrtcl_rotVel(1,idem_end)
            close(idem)
          enddo        
        endif
#endif
      enddo
      if(mod(itime,ivstats)==0) then
        call Update_FluidIndicator_imp(FluidIndicator) 
        call clcStat(ux,uy,uz,pressure,RealArr1,RealArr2)
        call clcStat_User(ux,uy,uz,pressure)
        call ClcCAStatistics()
        call dump_plane(itime,ux,uy,uz,pressure)   
      endif
      if(mod(itime,SaveVisu)== 0) call dump_visu(itime,ux,uy,uz,pressure,RealArr1)
      if(mod(itime,BackupFreq)== 0 .or. itime==ilast) then
        call Write_Restart(itime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
        call Delete_Prev_Restart(itime)
        call DEM%Write_Restart(itime*icouple)
      endif
      call total_timer%finish()

      ! command window and log file output
      IF((itime==ifirst .or. mod(itime, Cmd_LFile_Freq)==0) ) THEN
        call CheckDivergence(ux,uy,uz, divmax2)  
        if(nrank==0 .and. divmax2>div_limit) call MainLog%CheckForError(ErrT_Abort,"CFDACM_Iterate","too big div: "//strip(num2str(divmax2)))
        vmaxabs = CalcVmax(ux,uy,uz)
        if(nrank==0 .and. minval(vmaxabs)>vel_limit) call MainLog%CheckForError(ErrT_Abort,"CFDACM_Iterate","too big velocity: "//strip(num2str(vmaxabs(1)))//", "//strip(num2str(vmaxabs(2)))//", "//strip(num2str(vmaxabs(3))) )
        uxm = CalcUxAver(ux)
        if(nrank==0) then
          call MainLog%OutInfo("CFDACM performed "//strip(num2str(itime))//" iterations up to here!",1)
          call MainLog%OutInfo("Execution time [tot, last, ave] [sec]: "//strip(num2str(total_timer%tot_time))//", "// &
          strip(num2str(total_timer%last_time ))//", "//strip(num2str(total_timer%average())),2)
          call MainLog%OutInfo("SimTime | dt | CFL : "//strip(num2str(SimTime))//' | '//strip(num2str(dt))//' | '//strip(num2str(cflmp)),3)
          call MainLog%OutInfo("Max Abs Div: "//strip(num2str(divmax1))//" | "//strip(num2str(divmax2)) ,3)
          call MainLog%OutInfo("Max Abs Vel: "//strip(num2str(vmaxabs(1)))//" | "//strip(num2str(vmaxabs(2)))//" | "//strip(num2str(vmaxabs(3))), 3)
          call MainLog%OutInfo("Mean Velocity in streamwise: "//strip(num2str(uxm)),3)
        endif
      ENDIF
    ENDDO
  end subroutine CFDACM_Iterate_imp
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clc_uStar
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90 
  subroutine clc_uStar(uxStar,uyStar,uzStar,ux,uy,uz,RhsX,RhsY,RhsZ); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::uxStar,uyStar,uzStar
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in) ::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(in)::RhsX,RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::RhsZ

    ! locals
    integer:: ic,jc,kc
    
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          uxStar(ic,jc,kc)= ux(ic,jc,kc)+ RhsX(ic,jc,kc)
          uyStar(ic,jc,kc)= uy(ic,jc,kc)+ RhsY(ic,jc,kc)
          uzStar(ic,jc,kc)= uz(ic,jc,kc)+ RhsZ(ic,jc,kc)         
        enddo
      enddo
    ENDDO
  end subroutine clc_uStar

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clcU123Hat
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90 
  subroutine clcU123Hat(uxStar,uyStar,uzStar,ux,uy,uz); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::uxStar,uyStar,uzStar
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz

    ! locals
    integer:: ic,jc,kc
    
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)= uxStar(ic,jc,kc)
          uy(ic,jc,kc)= uyStar(ic,jc,kc)
          uz(ic,jc,kc)= uzStar(ic,jc,kc)     
        enddo
      enddo
    ENDDO
  end subroutine clcU123Hat
end module ca_system
