#include "definitions_inc.f90"
module cd_System
  use MPI
  use mc_Timer
  use mc_TypeDef
  use mc_Decomp2d
  use f2_Variables
  use f2_Parameters
  use f2_CFDSystem
  use sp_System
  use sp_Variables
  use sp_DumpPrtcl
  use sp_Parameters  
  use cd_FpForce
  use cd_Statistics
  implicit none
  private
  type(timer):: CoupleTimer
  
  public::CFDDEM_Iterate
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CFDDEM_Iterate
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CFDDEM_Iterate(); implicit none

    !locals
    integer:: idem,idem_start,idem_end

#ifdef SeveralSphereInfo
    integer::pid,iUnit
    character(len=:),allocatable::chFile
    do pid=1,GPrtcl_list%nlocal
      chFile = strip(ResultsDir_) // "SphereInfo_" // strip(RunName_) // "_" // int2str(GPrtcl_id(pid),2) // ".txt"
      open(newunit=iUnit,file=chFile,status='replace',form='formatted')
      write(iUnit,'(A)')' t, x, y, z, r, vx, vy, vz, wx, wy, wz'
      close(iUnit)
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
      if(mod(itime,ivstats)==0) call ClcCDStatistics()

      ! CFD Iterate
      call CFDIterate()

#ifdef SeveralSphereInfo
      do pid=1,GPrtcl_list%nlocal
        chFile = strip(ResultsDir_)//"SphereInfo_"//strip(RunName_)//"_"//int2str(GPrtcl_id(pid),2)//".txt"
        open(newunit=iUnit,file=chFile,status='old',position='append',form='formatted')
        write(iUnit,'(15ES24.15)') SimTime,GPrtcl_PosR(pid),GPrtcl_linVel(1,pid),GPrtcl_rotVel(1,pid)
        close(iUnit)
      enddo
#endif
      if(nrank==0 .and. mod(itime, Cmd_LFile_Freq)==0) then
        call MainLog%OutInfo("Coupling time  [tot, last, ave] [sec]: "//strip(num2str(CoupleTimer%tot_time))//", "// &
            strip(num2str(CoupleTimer%last_time ))//", "//strip(num2str(CoupleTimer%average())),2)
      endif
    enddo
  end subroutine CFDDEM_Iterate

end module cd_System
