module PPD_Parameters
  use PPD_TypeDef
  implicit none
  private

  integer,save,public::nrank,nproc


  character(64),public,save:: Stat_Dir      ! Statistic directory
  character(64),save,public:: StatSuffix

  logical,save,public::IsReadPrtclDumpOut
  logical,save,public::IsRewritePrtclDump
  character(64),save,public::PrtclDumpDir
  character(64),save,public::PrtclDumpOutDir
  character(64),save,public::PrtclDumpFilter

  logical,save,public:: IsSortPrtclDumpOut
  logical,save,public:: IsAddPeriodicLen
  real(RK),save,public::xlx,yly,zlz

  logical,save,public:: IsShrinkPrtclDump
  logical,save,public::IsRewriteShrink
  integer,save,public::ShrinkDumpFreq
  integer,save,public::ShrinkDumpNeglect
  character(64),save,public::ShrinkDumpDir

  logical,save,public:: IsClcHopStat
  character(64),save,public:: ReadHopDir
  integer,save,public:: iDetectHopSpan
  integer,save,public:: iDetectHopSet
  real(RK),save,public::dtInterval
  real(RK),save,public::moveDistLimit
  real(RK),save,public::yLimit
  real(RK),save,public::VelBinSetX
  real(RK),save,public::VelBinDeltaX
  integer,save,public:: VelBinNumX
  real(RK),save,public::VelBinSetZ
  real(RK),save,public::VelBinDeltaZ
  integer,save,public:: VelBinNumZ
  real(RK),save,public::AccBinSetX
  real(RK),save,public::AccBinDeltaX
  integer,save,public:: AccBinNumX
  real(RK),save,public::AccBinSetZ
  real(RK),save,public::AccBinDeltaZ
  integer,save,public:: AccBinNumZ
  real(RK),save,public::TimeBinSet
  real(RK),save,public::TimeBinDelta
  integer,save,public:: TimeBinNum
  real(RK),save,public::LenBinSetX
  real(RK),save,public::LenBinDeltaX
  integer,save,public:: LenBinNumX
  real(RK),save,public::LenBinSetZ
  real(RK),save,public::LenBinDeltaZ
  integer,save,public:: LenBinNumZ

  logical,save,public:: IsClcCrossDiffusion
  character(64),save,public:: ReadDiffuseDir
  real(RK),save,public::CrossDiffBinSet
  real(RK),save,public::CrossDiffBinDelta
  integer,save,public::CrossDiffusionBinNum
  integer,save,public::iDetectCrossSet

  public:: ReadParameters
contains

  !******************************************************************
  ! ReadParameters
  !******************************************************************
  subroutine ReadParameters(chFile)
    implicit none
    character(*),intent(in)::chFile

    ! locals
    integer:: myistat
    NAMELIST/BasicParam/Stat_Dir,StatSuffix,IsReadPrtclDumpOut,PrtclDumpDir,PrtclDumpOutDir,PrtclDumpFilter,  &
                        IsRewritePrtclDump,IsSortPrtclDumpOut,IsAddPeriodicLen,xlx,yly,zlz,                   &
                        IsShrinkPrtclDump,IsRewriteShrink,ShrinkDumpDir,ShrinkDumpFreq,ShrinkDumpNeglect,     &
                        IsClcHopStat,ReadHopDir,iDetectHopSpan,iDetectHopSet,dtInterval,moveDistLimit,yLimit, &
                          VelBinSetX,VelBinDeltaX,VelBinNumX,VelBinSetZ,VelBinDeltaZ,VelBinNumZ,AccBinSetX,   &
                          AccBinDeltaX,AccBinNumX,AccBinSetZ,AccBinDeltaZ,AccBinNumZ,TimeBinSet,TimeBinDelta, &
                          TimeBinNum,LenBinSetX,LenBinDeltaX,LenBinNumX,LenBinSetZ,LenBinDeltaZ,LenBinNumZ,   &
                        IsClcCrossDiffusion,ReadDiffuseDir,CrossDiffBinSet,CrossDiffBinDelta,CrossDiffusionBinNum,iDetectCrossSet
                                 
    open(25,file=chFile,status='old',form='formatted',IOSTAT=myistat)
    if(myistat /= 0) then
      print*,"Cannot open file: "//trim(adjustl(chFile)); STOP
    endif
    read(25, nml=BasicParam)
    close(25)

  end subroutine ReadParameters

end module PPD_Parameters
