module ca_system
  use MPI
  use ca_IBM
  use m_Timer
  use m_Tools
  use m_TScheme
  use m_Typedef
  use m_LogInfo
  use m_Poisson
  use m_Decomp2d
  use m_FlowCase
  use Prtcl_Comm
  use m_Variables
  use m_LESModels
  use m_IOAndVisu
  use m_Parameters
  use Prtcl_System 
  use m_BC_and_Halo
  use ca_Statistics
  use Prtcl_Property
  use Prtcl_variables
  use ca_IBM_implicit
  use Prtcl_Parameters
  use ca_BC_and_Halo,only: SetBC_and_UpdateHalo_VelIBM
  implicit none
  private
 
  !// timers
  type(timer):: total_timer
  public::ChannelACM_Iterate
contains

  !******************************************************************
  ! ChannelACM_Iterate
  !******************************************************************
  subroutine ChannelACM_Iterate()
    implicit none

    ! locals
    integer::ns,iForcingExtra,idem_start,idem_end,idem
    real(RK)::uddxmax,cflmp,divmax1,divmax2,uxm,vmaxabs(3)

#ifdef SeveralSphereInfo
    character(128)::chFile
    type(real3)::Ave_FpForce
    if(DEM_Opt%numPrtcl==1 .and. GPrtcl_list%nlocal==1) then
      write(chFile,"(A,I2.2,A)") trim(ResultsDir)//"SphereInfo_"//trim(RunName)//"_",GPrtcl_id(1),".txt"
      open(newunit=idem,file=chFile,status='replace',form='formatted')
      write(idem,'(A)')' t, x, y, z, r, vx, vy, vz, wx, wy, wz fx fy fz'
      close(idem)    
    else
      do ns=1,GPrtcl_list%nlocal
        write(chFile,"(A,I2.2,A)") trim(ResultsDir)//"SphereInfo_"//trim(RunName)//"_",GPrtcl_id(ns),".txt"
        open(newunit=idem,file=chFile,status='replace',form='formatted')
        write(idem,'(A)')' t, x, y, z, r, vx, vy, vz, wx, wy, wz'
        close(idem)
      enddo    
    endif
#endif
#ifdef ChanBraunJFM2011
    character(len=128)::filename
    integer::nfstime,nUnit,ierror
    real(RK),dimension(8,3)::ForceAndTorque
    type(real3),dimension(:),allocatable::Ave_FpForce
    type(real3),dimension(:),allocatable::Ave_FpTorque
    
    nfstime=0
    ForceAndTorque=0.0_RK
    if(GPrtcl_list%nlocal>0) then
      allocate(Ave_FpForce(GPrtcl_list%nlocal))
      allocate(Ave_FpTorque(GPrtcl_list%nlocal))
    endif

    if(nrank==0) then
      write(filename,'(A,I10.10,A)')trim(DEM_opt%ResultsDir)//"FpForce",ilast,'.txt'
      open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
      if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"ChannelACM_Iterate","Cannot open file: "//trim(filename))
      close(nUnit,IOSTAT=ierror)
    endif
#endif

    if(.not. DEM_Opt%RestartFlag) then
      dt= dtMax
      call PMcoeUpdate(iadvance)
      SimTime=0.0_RK
      IF(IBM_Scheme<2) THEN
        call PrepareIBM_interp()
        call IntegrateFluidPrtclForce(ux,uy,uz)
      ELSE
        call PrepareIBM_interp_imp1()
        call IntegrateFluidPrtclForce_imp(ux,uy,uz)
      ENDIF 
    endif
    DO itime=ifirst, ilast
      call total_timer%start()
      call CalcMaxCFL(ux,uy,uz,uddxmax)
      cflmp=uddxmax*dt
      if( icfl==1 ) then
        dt= CFLc/uddxmax
        dt= min(dt,dtMax)
        DEM_Opt%dt= dt/real(icouple,RK)
      else
        dt= dtMax  
      endif
      if(LES_type>0) call Clc_SGS_Vis(ux,uy,uz,nut)

      do ns=1,iadvance
        ! step0: Update the Projection Method coefficients.
        call PMcoeUpdate(ns)
        call Update_uy_ym(uy_ym, duy_ym, SimTime)
        idem_start=(itime-1)*icouple+ idem_advance_start(ns)
        idem_end  =(itime-1)*icouple+ idem_advance_end(ns)

        ! step1: Calculate the right hand side of the three velocity equations.
        asso_RHS123: associate( RhsX=>RealArr1, RhsY=>RealArr2, RhsZ=>RealHalo)
        call clcRhsX(ux,uy,uz,RhsX,HistXOld,pressure)
        call clcRhsY(ux,uy,uz,RhsY,HistYOld,pressure)
        call clcRhsZ(ux,uy,uz,RhsZ,HistZOld,pressure)

        IF(IBM_Scheme<2) THEN
          ! IBM Part1 for explicit coupling ======================
          associate_uStar: associate(uxStar=>IBMArr1, uyStar=>IBMArr2, uzStar=>IBMArr3)
          call clc_uStar(uxStar,uyStar,uzStar,ux,uy,uz,RhsX,RhsY,RhsZ)
          call SetBC_and_UpdateHalo_VelIBM( uxStar,uyStar,uzStar, uy_ym )
          call PrepareIBM_interp()
          call InterpolationAndForcingIBM(uxStar,uyStar,uzStar)
          end associate associate_uStar

          associate_VolForce: associate(VolForce_x=>IBMArr1,VolForce_y=>IBMArr2,VolForce_z=>IBMArr3)
          call SpreadIbpForce(VolForce_x,VolForce_y,VolForce_z)
          call updateRhsIBM(RhsX,RhsY,RhsZ,VolForce_x,VolForce_y,VolForce_z)
          end associate associate_VolForce

        ELSE
          ! IBM Part1 for semi-implicit coupling =================
          associate_uStar_imp: associate(uxStar=>IBMArr1, uyStar=>IBMArr2, uzStar=>IBMArr3)
          call clc_uStar(uxStar,uyStar,uzStar,ux,uy,uz,RhsX,RhsY,RhsZ)
          call SetBC_and_UpdateHalo_VelIBM( uxStar,uyStar,uzStar, uy_ym )
          call PrepareIBM_interp_imp1()
          call InterpolationIBM_imp(uxStar,uyStar,uzStar)
          call IntegrateFluidPrtclForce_imp(uxStar,uyStar,uzStar)
          end associate associate_uStar_imp
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
        ENDIF
         
        ! step2: Calculate the Uhat
        call clcOutFlowVelocity(ux,uy,uz)
        call clcU1Hat(ux,RhsX)
        call clcU2Hat(uy,RhsY,duy_ym)
        call clcU3Hat(uz,RhsZ) 
        end associate asso_RHS123

        ! IBM Part2 for explicit coupling ========================
        IF(IBM_Scheme<2) THEN
          associate_VolForce2: associate(VolForce_x=>IBMArr1,VolForce_y=>IBMArr2,VolForce_z=>IBMArr3)
          do iForcingExtra=1, nForcingExtra
            call SetBC_and_UpdateHalo_VelIBM( ux,uy,uz, uy_ym )
            call AdditionalForceIBM(ux,uy,uz,VolForce_x,VolForce_y,VolForce_z)
          enddo
          end associate associate_VolForce2
        ENDIF

        ! step3: Calculate the source term of the PPE 
        call SetBC_and_UpdateHaloForPrSrc( ux,uy,uz, uy_ym )
        call correctOutFlowFaceVelocity(ux,uy,uz)
        asso_Pr: associate(prsrc =>RealArr1, prphi =>RealArr2, prphiHalo =>RealHalo  )
        call clcPrSrc(ux,uy,uz,prsrc,pressure,divmax1)
        call clcPPE(prsrc,prphiHalo)
        call SetBC_and_UpdateHalo_pr(prphiHalo)
            
        ! step4: Update the velocity field to get the final real velocity.
        call FluidVelUpdate(prphiHalo,ux,uy,uz)
            
        ! step5: Update the real pressure field  to get the final pressure.
        call PressureUpdate(pressure, prphiHalo)
        end associate asso_Pr
        call SetBC_and_UpdateHalo( ux,uy,uz,uy_ym )
        call SetBC_and_UpdateHalo_pr( pressure )

        ! DEM Part for explicit coupling ========================= 
        IF(IBM_Scheme<2) THEN
          call IntegrateFluidPrtclForce(ux,uy,uz)
          do idem=idem_start,idem_end
            call DEM%iterate(idem)
          enddo
        ENDIF
        
#ifdef ChanBraunJFM2011
#include "ca_Force_inc.f90"
#endif
#ifdef SeveralSphereInfo
        if(ns==1 .and. DEM_Opt%numPrtcl==1) Ave_FpForce=zero_r3
        if(DEM_Opt%numPrtcl==1 .and. GPrtcl_list%nlocal==1) then
          Ave_FpForce =Ave_FpForce +pmAlphaC*GPrtcl_FpForce(1)
          write(chFile,"(A,I2.2,A)") trim(ResultsDir)//"SphereInfo_"//trim(RunName)//"_",GPrtcl_id(1),".txt"
          open(newunit=idem,file=chFile,status='old',position='append',form='formatted')
          write(idem,'(15ES24.15)')SimTime,GPrtcl_PosR(1),GPrtcl_linVel(1,1),GPrtcl_rotVel(1,1),Ave_FpForce
          close(idem) 
        else
          do iForcingExtra=1,GPrtcl_list%nlocal
            write(chFile,"(A,I2.2,A)") trim(ResultsDir)//"SphereInfo_"//trim(RunName)//"_",GPrtcl_id(iForcingExtra),".txt"
            open(newunit=idem,file=chFile,status='old',position='append',form='formatted')
            write(idem,'(15ES24.15)')SimTime,GPrtcl_PosR(iForcingExtra),GPrtcl_linVel(1,iForcingExtra),GPrtcl_rotVel(1,iForcingExtra)
            close(idem)
          enddo        
        endif
#endif
      enddo
      if(mod(itime,ivstats)==0) then
        IF(IBM_Scheme<2) THEN
          call Update_FluidIndicator(FluidIndicator)
        ELSE
          call Update_FluidIndicator_imp(FluidIndicator)
        ENDIF  
        call clcStat(ux,uy,uz,pressure)
        call ClcCAStatistics()      
      endif
      if(mod(itime,SaveVisu)== 0)   call dump_visu(itime,ux,uy,uz,pressure,RealArr1)
      if(mod(itime,BackupFreq)== 0 .or. itime==ilast) then
        call Write_Restart(itime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
        call Delete_Prev_Restart(itime)
      endif
      call total_timer%finish()

      ! command window and log file output
      IF((itime==ifirst .or. mod(itime, Cmd_LFile_Freq)==0) ) THEN
        call CheckDivergence(ux,uy,uz, divmax2)  
        if(nrank==0 .and. divmax2>div_limit) call MainLog%CheckForError(ErrT_Abort,"ChannelACM_Iterate","too big div: "//trim(num2str(divmax2)))
        vmaxabs = CalcVmax(ux,uy,uz)
        if(nrank==0 .and. minval(vmaxabs)>vel_limit) call MainLog%CheckForError(ErrT_Abort,"ChannelACM_Iterate","too big velocity: "//trim(num2str(vmaxabs(1)))//", "//trim(num2str(vmaxabs(2)))//", "//trim(num2str(vmaxabs(3))) )
        uxm = CalcUxAver(ux)
        if(nrank==0) then
          call MainLog%OutInfo("ChannelACM performed "//trim(num2str(itime))//" iterations up to here!",1)
          call MainLog%OutInfo("Execution time [tot, last, ave] [sec]: "//trim(num2str(total_timer%tot_time))//", "// &
          trim(num2str(total_timer%last_time ))//", "//trim(num2str(total_timer%average())),2)
          call MainLog%OutInfo("SimTime | dt | CFL : "//trim(num2str(SimTime))//' | '//trim(num2str(dt))//' | '//trim(num2str(cflmp)),3)
          call MainLog%OutInfo("Max Abs Div: "//trim(num2str(divmax1))//" | "//trim(num2str(divmax2)) ,3)
          call MainLog%OutInfo("Max Abs Vel: "//trim(num2str(vmaxabs(1)))//" | "//trim(num2str(vmaxabs(2)))//" | "//trim(num2str(vmaxabs(3))), 3)
          call MainLog%OutInfo("Mean Velocity in streamwise: "//trim(num2str(uxm)),3)
        endif
      ENDIF
    ENDDO
  end subroutine ChannelACM_Iterate

  !******************************************************************
  ! clc_uStar
  !****************************************************************** 
  subroutine clc_uStar(uxStar,uyStar,uzStar,ux,uy,uz,RhsX,RhsY,RhsZ)
    implicit none
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

end module ca_system
