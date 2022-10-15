module cd_System
  use MPI
  use m_Timer
  use m_Typedef
  use m_LogInfo
  use m_Decomp2d
  use cd_FpForce
  use m_Variables
  use m_Parameters
  use cd_Statistics
  use m_ChannelSystem
  use Prtcl_Variables
  use Prtcl_DEMSystem
  use Prtcl_DumpPrtcl
  use Prtcl_Parameters
  implicit none
  private
  type(timer):: CoupleTimer
  
  public::ChannelDEM_Iterate
contains

  !******************************************************************
  ! ChannelDEM_Iterate
  !******************************************************************
  subroutine ChannelDEM_Iterate()
    implicit none

    !locals
    integer:: idem,idem_start,idem_end

#ifdef SeveralSphereInfo
    integer::pid,nUnit
    character(128)::chFile
    do pid=1,GPrtcl_list%nlocal
      write(chFile,"(A,I2.2,A)") trim(ResultsDir)//"SphereInfo_"//trim(RunName)//"_",GPrtcl_id(pid),".txt"
      open(newunit=nUnit,file=chFile,status='replace',form='formatted')
      write(nUnit,'(A)')' t, x, y, z, r, vx, vy, vz, wx, wy, wz'
      close(nUnit)
    enddo 
#endif
    call CoupleTimer%reset()
    do itime=ifirst, ilast
      ! CFD-DEM coupling part
      call CoupleTimer%start()
      asso_Fpforce: associate(RatioYp_interp =>GPrtcl_cntctForce, RatioYc_interp =>GPrtcl_torque)
      call PrepareDistribute(RatioYp_interp,RatioYc_interp)
      call clc_FpForce(ux,uy,uz,pressure,RatioYp_interp,RatioYc_interp)
      end associate asso_Fpforce
      call distribute_FpForce()
      call FinalDistribute()
      call CoupleTimer%finish()

      ! DEM Iterate
      idem_start=(itime-1)*icouple+1
      idem_end=itime*icouple
      do idem=idem_start,idem_end
        call DEM%iterate(idem)
      enddo
      if(mod(itime,ivstats)==0)call ClcCDStatistics()

      ! CFD Iterate
      call ChannelIterate()

#ifdef SeveralSphereInfo
      do pid=1,GPrtcl_list%nlocal
        write(chFile,"(A,I2.2,A)") trim(ResultsDir)//"SphereInfo_"//trim(RunName)//"_",GPrtcl_id(pid),".txt"
        open(newunit=nUnit,file=chFile,status='old',position='append',form='formatted')
        write(nUnit,'(15ES24.15)')SimTime,GPrtcl_PosR(pid),GPrtcl_linVel(1,pid),GPrtcl_rotVel(1,pid)
        close(nUnit)
      enddo
#endif
      if(nrank==0 .and. mod(itime, Cmd_LFile_Freq)==0) then
        call MainLog%OutInfo("Coupling time  [tot, last, ave] [sec]: "//trim(num2str(CoupleTimer%tot_time))//", "// &
            trim(num2str(CoupleTimer%last_time ))//", "//trim(num2str(CoupleTimer%average())),2)
      endif
    enddo
  end subroutine ChannelDEM_Iterate

end module cd_System
