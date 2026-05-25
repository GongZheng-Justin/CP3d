module PPD_Statistics
  use MPI
  use PPD_Tools
  use PPD_TypeDef
  use PPD_Parameters
  implicit none
  private

  public:: ClcHopStat, ClcCrossDiffusion
contains

  !******************************************************************
  ! ClcHopStat
  !******************************************************************
  subroutine ClcHopStat()
    implicit none

    ! locals
    character::Dummy
    character(14)::PrtclStr
    type(DumpPrtclVarOut)::DumpVarT
    character(128)::FileName1,FileName2,WriteHopDir
    integer,dimension(:),allocatable::SetVec,EndVec
    type(DumpPrtclVarOut),allocatable,dimension(:)::PrtclDataIn
    integer::pid,fid,j,k,FileByte,SinglePrtclByte,nPrtclSize,ierror,iSet,iEnd,SumHop,ContinueFlag,mn
    real(RK)::dispDist,HopDistX,HopDistZ,HopTime,yPosMin
    logical::IsFree

    !===== allocate variables for PDF file =====!
    integer::BinId,SumBinCount(6),SumBinCountR(6)
    real(RK)::VelX,VelZ,AccX,AccZ,SumBinValue(8),SumBinValueR(8)
    integer,dimension(:),allocatable::CountVelBinX,CountVelBinXR,CountVelBinZ,CountVelBinZR
    integer,dimension(:),allocatable::CountAccBinX,CountAccBinXR,CountAccBinZ,CountAccBinZR
    integer:: SumHopCount(4),CountTemp
    real(RK)::SumHopStat(5)
    integer,dimension(:),allocatable::CountHopTime,CountHopLenX,COuntHopLenZ
    
    SumBinCount=0; SumBinCountR=0; SumBinValue=0.0_RK; SumBinValueR=0.0_RK;
    allocate(CountVelBinX(VelBinNumX),CountVelBinXR(VelBinNumX));CountVelBinX=0;CountVelBinXR=0
    allocate(CountVelBinZ(VelBinNumZ),CountVelBinZR(VelBinNumZ));CountVelBinZ=0;CountVelBinZR=0
    allocate(CountAccBinX(AccBinNumX),CountAccBinXR(AccBinNumX));CountAccBinX=0;CountAccBinXR=0
    allocate(CountAccBinZ(AccBinNumZ),CountAccBinZR(AccBinNumZ));CountAccBinZ=0;CountAccBinZR=0
    SumHopCount=0;
    SumHopStat=0.0_RK
    allocate(CountHopTime(TimeBinNum)); CountHopTime=0
    allocate(CountHopLenX(LenBinNumX)); CountHopLenX=0
    allocate(CountHopLenZ(LenBinNumZ)); CountHopLenZ=0
    !===== allocate variables for PDF file =====!

    WriteHopDir="./HoStatTemp/"
    SinglePrtclByte=sizeof(DumpVarT)
    if(nrank==0) then
      call CreateDir(.true.,WriteHopDir)
      call FileSystem%WriteFileNameList(ReadHopDir,"pid*")
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call FileSystem%GetFileName(ReadHopDir)
    do fid=1,FileSystem%nFile
      if(mod(fid,nproc)/=nrank) cycle
      FileName1=trim(adjustl(ReadHopDir))//trim(adjustl(FileSystem%FileName(fid)))
      FileName2=trim(adjustl(WriteHopDir))//trim(adjustl(FileSystem%FileName(fid)))//'.txt'

      FileByte=GetFileByte(FileName1)
      if(mod(FileByte,SinglePrtclByte)/=0) then
        print*,"ClcHopStat: FileByte wrong";stop
      endif
      nPrtclSize=FileByte/SinglePrtclByte
      if(nPrtclSize<1) cycle
      if(iDetectHopSet+1>=nPrtclSize) cycle

      allocate(PrtclDataIn(nPrtclSize),SetVec(nPrtclSize),EndVec(nPrtclSize))
      open(35,file=FileName1,status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
      read(35)PrtclDataIn
      close(35)

      !========== Start to Calculate Hop Statistics ==========!
      iSet=iDetectHopSet+1;iEnd=iSet;SumHop=0;ContinueFlag=1;
      if(PrtclDataIn(nPrtclSize)%Pos(1)-PrtclDataIn(iSet)%Pos(1)<1.1E-2) then
        deallocate(PrtclDataIn,SetVec,EndVec);cycle
      endif
      dispDist=ClcDisplacement(PrtclDataIn,iSet)
      if(dispDist>=moveDistLimit) then
        do j=iSet,nPrtclSize-1
         dispDist=ClcDisplacement(PrtclDataIn,j)
         if(dispDist<moveDistLimit) then
           iSet=j;ContinueFlag=1;exit
         endif
         ContinueFlag=0
        enddo
      endif
      do
        if(ContinueFlag==0)exit
        do j=iSet,nPrtclSize-1
          dispDist=ClcDisplacement(PrtclDataIn,j)
          if(dispDist>=moveDistLimit) then
            iSet=j; ContinueFlag=1;exit
          endif
          ContinueFlag=0
        enddo

        if(ContinueFlag==0)exit
        do j=iSet,nPrtclSize-1
          dispDist=ClcDisplacement(PrtclDataIn,j)
          if(dispDist<moveDistLimit) then
            iEnd=j
            yPosMin=ClcYmin(PrtclDataIn,iSet,iEnd)
            IsFree=.false.
            do mn=iSet,iEnd
              if(PrtclDataIn(mn)%cntctFlag==0) then
                IsFree=.true.;exit
              endif
            enddo
            if(iEnd-iSet>iDetectHopSpan .and. yPosMin>=yLimit ) then
              if(PrtclDataIn(iEnd)%Pos(1)-PrtclDataIn(iSet)%Pos(1) >5.7E-3 .and. PrtclDataIn(iEnd)%Pos(1)-PrtclDataIn(iSet)%Pos(1)<=73E-3) then
                SumHop=SumHop+1
                SetVec(SumHop)=iSet
                EndVec(SumHop)=iEnd
              endif
            endif
            iSet=j;ContinueFlag=1;exit
          endif
          ContinueFlag=0
        enddo
      enddo
      if(SumHop>0) then
        open(36,file=FileName2,status='replace',form='formatted',action='write',IOSTAT=ierror)
        do j=1,SumHop
          iSet=SetVec(j); iEnd=EndVec(j)
          HopTime=real(iEnd-iSet,kind=RK)*dtInterval
          HopDistX=PrtclDataIn(iEnd)%Pos(1)-PrtclDataIn(iSet)%Pos(1)
          HopDistZ=PrtclDataIn(iEnd)%Pos(3)-PrtclDataIn(iSet)%Pos(3)
          write(36,'(2I10,9ES18.9)')iSet,iEnd,HopTime,HopDistX,HopDistZ

          do pid=iSet,iEnd-1
            VelX=(PrtclDataIn(pid+1)%Pos(1)-PrtclDataIn(pid)%Pos(1))/dtInterval
            BinId=floor((VelX-VelBinSetX)/VelBinDeltaX)+1
            if(BinId>0 .and. BinId<=VelBinNumX) then
              CountVelBinX(BinId)=CountVelBinX(BinId)+1
              SumBinCount(1)=SumBinCount(1)+1
            endif

            VelZ=(PrtclDataIn(pid+1)%Pos(3)-PrtclDataIn(pid)%Pos(3))/dtInterval
            BinId=floor((VelZ-VelBinSetZ)/VelBinDeltaZ)+1
            if(BinId>0 .and. BinId<=VelBinNumZ) then
              CountVelBinZ(BinId)=CountVelBinZ(BinId)+1
              SumBinCount(2)=SumBinCount(2)+1
            endif

            SumBinCount(3)=SumBinCount(3)+1
            SumBinValue(1)=SumBinValue(1)+VelX
            SumBinValue(2)=SumBinValue(2)+abs(VelX)
            SumBinValue(3)=SumBinValue(3)+VelZ
            SumBinValue(4)=SumBinValue(4)+abs(VelZ)
          enddo
          do pid=iSet,iEnd-2
            AccX=(PrtclDataIn(pid+2)%Pos(1)-2.0*PrtclDataIn(pid+1)%Pos(1)+PrtclDataIn(pid)%Pos(1))/(dtInterval*dtInterval)
            BinId=floor((AccX-AccBinSetX)/AccBinDeltaX)+1
            if(BinId>0 .and. BinId<=AccBinNumX) then
              CountAccBinX(BinId)=CountAccBinX(BinId)+1
              SumBinCount(4)=SumBinCount(4)+1
            endif

            AccZ=(PrtclDataIn(pid+2)%Pos(3)-2.0*PrtclDataIn(pid+1)%Pos(3)+PrtclDataIn(pid)%Pos(3))/(dtInterval*dtInterval)
            BinId=floor((AccZ-AccBinSetZ)/AccBinDeltaZ)+1
            if(BinId>0 .and. BinId<=AccBinNumZ) then
              CountAccBinZ(BinId)=CountAccBinZ(BinId)+1
              SumBinCount(5)=SumBinCount(5)+1
            endif

            SumBinCount(6)=SumBinCount(6)+1
            SumBinValue(5)=SumBinValue(5)+AccX
            SumBinValue(6)=SumBinValue(6)+abs(AccX)
            SumBinValue(7)=SumBinValue(7)+AccZ
            SumBinValue(8)=SumBinValue(8)+abs(AccZ)
          enddo
        enddo
        close(36,IOSTAT=ierror)
      endif

      deallocate(PrtclDataIn,SetVec,EndVec)
      write(*,'(A,I4,A)')'rank',nrank,' ClcHopStat: '//trim(adjustl(FileName1))//' successfully'
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(SumBinCount, SumBinCountR, 6,         MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(SumBinValue, SumBinValueR, 8,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(CountVelBinX,CountVelBinXR,VelBinNumX,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(CountVelBinZ,CountVelBinZR,VelBinNumZ,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(CountAccBinX,CountAccBinXR,AccBinNumX,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(CountAccBinZ,CountAccBinZR,AccBinNumZ,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank/=0) return
    call FileSystem%DeleteFileNameList()

    ! Gather and dump the Hop statistics
    write(FileName1,"(A,A,A)")trim(adjustl(Stat_Dir))//'LinVelXPDF_',trim(StatSuffix),'.txt'
    open(37,file=FileName1,status='replace',form='formatted',action='write',IOSTAT=ierror)
    write(37,'(A,I18)')' Total Sample :    ', SumBinCountR(3)
    write(37,'(A,I18)')' Total PDF Sample :', SumBinCountR(1)
    write(37,'(A,ES18.9)') ' Mean Velocity (cm/s): ', SumBinValueR(1)/real(SumBinCountR(3),RK)*100_RK
    write(37,'(A,ES18.9)') ' Mean Abs Velo (cm/s): ', SumBinValueR(2)/real(SumBinCountR(3),RK)*100_RK
    write(37,'(A)')' VelX(cm/s) nP PD'
    do pid=1,VelBinNumX
      write(37,'(ES18.9,I18,ES18.9)')(VelBinSetX+(real(pid,RK)-0.5_RK)*VelBinDeltaX)*100_RK,CountVelBinXR(pid),real(CountVelBinXR(pid),RK)/real(SumBinCountR(1),RK)
    enddo
    close(37,IOSTAT=ierror)
    write(FileName1,"(A,A,A)")trim(adjustl(Stat_Dir))//'LinVelZPDF_',trim(StatSuffix),'.txt'
    open(37,file=FileName1,status='replace',form='formatted',action='write',IOSTAT=ierror)
    write(37,'(A,I18)')' Total Sample :    ', SumBinCountR(3)
    write(37,'(A,I18)')' Total PDF Sample :', SumBinCountR(2)
    write(37,'(A,ES18.9)') ' Mean Velocity (cm/s): ', SumBinValueR(3)/real(SumBinCountR(3),RK)*100_RK
    write(37,'(A,ES18.9)') ' Mean Abs Velo (cm/s): ', SumBinValueR(4)/real(SumBinCountR(3),RK)*100_RK
    write(37,'(A)')' VelZ(cm/s) nP PD'
    do pid=1,VelBinNumZ
      write(37,'(ES18.9,I18,ES18.9)')(VelBinSetZ+(real(pid,RK)-0.5_RK)*VelBinDeltaZ)*100_RK,CountVelBinZR(pid),real(CountVelBinZR(pid),RK)/real(SumBinCountR(2),RK)
    enddo
    close(37,IOSTAT=ierror)

    write(FileName1,"(A,A,A)")trim(adjustl(Stat_Dir))//'LinAccXPDF_',trim(StatSuffix),'.txt'
    open(38,file=FileName1,status='replace',form='formatted',action='write',IOSTAT=ierror)
    write(38,'(A,I18)')' Total Sample :    ', SumBinCountR(6)
    write(38,'(A,I18)')' Total PDF Sample :', SumBinCountR(4)
    write(38,'(A,ES18.9)') ' Mean Acceleration (cm2/s): ', SumBinValueR(5)/real(SumBinCountR(6),RK)*100_RK
    write(38,'(A,ES18.9)') ' Mean Abs Accelera (cm2/s): ', SumBinValueR(6)/real(SumBinCountR(6),RK)*100_RK
    write(38,'(A)')' AccX(cm/s2) nP PD'
    do pid=1,AccBinNumX
      write(38,'(ES18.9,I18,ES18.9)')(AccBinSetX+(real(pid,RK)-0.5_RK)*AccBinDeltaX)*100_RK,CountAccBinXR(pid),real(CountAccBinXR(pid),RK)/real(SumBinCountR(4),RK)
    enddo
    close(38,IOSTAT=ierror)
    write(FileName1,"(A,A,A)")trim(adjustl(Stat_Dir))//'LinAccZPDF_',trim(StatSuffix),'.txt'
    open(38,file=FileName1,status='replace',form='formatted',action='write',IOSTAT=ierror)
    write(38,'(A,I18)')' Total Sample :    ', SumBinCountR(6)
    write(38,'(A,I18)')' Total PDF Sample :', SumBinCountR(5)
    write(38,'(A,ES18.9)') ' Mean Acceleration (cm2/s): ', SumBinValueR(7)/real(SumBinCountR(6),RK)*100_RK
    write(38,'(A,ES18.9)') ' Mean Abs Accelera (cm2/s): ', SumBinValueR(8)/real(SumBinCountR(6),RK)*100_RK
    write(38,'(A)')' AccZ(cm/s2) nP PD'
    do pid=1,AccBinNumZ
      write(38,'(ES18.9,I18,ES18.9)')(AccBinSetZ+(real(pid,RK)-0.5_RK)*AccBinDeltaZ)*100_RK,CountAccBinZR(pid),real(CountAccBinZR(pid),RK)/real(SumBinCountR(5),RK)
    enddo
    close(38,IOSTAT=ierror)

    call FileSystem%WriteFileNameList(WriteHopDir,"pid*")
    call FileSystem%GetFileName(WriteHopDir)
    write(FileName1,"(A,A,A)")trim(adjustl(Stat_Dir))//'HopGather_',trim(StatSuffix),'.txt'
    open(39,file=FileName1,status='replace',form='formatted',action='write',IOSTAT=ierror)
    write(39,*)' pid iSet iEnd HopTime(s) HopDistX(cm) HopDistZ(cm)'
    do fid=1,FileSystem%nFile
      PrtclStr=trim(adjustl(FileSystem%FileName(fid)))
      read(PrtclStr,'(A3,I7,A4)')Dummy,pid,Dummy
      FileName2=trim(adjustl(WriteHopDir))//PrtclStr
      open(40,file=FileName2,status='old',form='formatted',action='read',IOSTAT=ierror)
      do
        read(40,'(2I10,9ES18.9)',IOSTAT=ierror)iSet,iEnd,HopTime,HopDistX,HopDistZ
        if(ierror/=0) exit
        write(39,'(I7,A,2I10,9ES18.9)')pid,' ',iSet,iEnd,HopTime,HopDistX*100_RK,HopDistZ*100_RK

        SumHopCount(1)=SumHopCount(1)+1
        SumHopStat(1) =SumHopStat(1) +HopTime
        SumHopStat(2) =SumHopStat(2) +HopDistX
        SumHopStat(3) =SumHopStat(3) +abs(HopDistX)
        SumHopStat(4) =SumHopStat(4) +HopDistZ
        SumHopStat(5) =SumHopStat(5) +abs(HopDistZ)

        BinId=floor((HopTime-TimeBinSet)/TimeBinDelta)+1
        if(BinId>0 .and. BinId<=TimeBinNum) then
          CountHopTime(BinId)=CountHopTime(BinId)+1
          SumHopCount(2)=SumHopCount(2)+1
        endif

        BinId=floor((HopDistX-LenBinSetX)/LenBinDeltaX)+1
        if(BinId>0 .and. BinId<=LenBinNumX) then
          CountHopLenX(BinId)=CountHopLenX(BinId)+1
          SumHopCount(3)=SumHopCount(3)+1
        endif

        BinId=floor((HopDistZ-LenBinSetZ)/LenBinDeltaZ)+1
        if(BinId>0 .and. BinId<=LenBinNumZ) then
          CountHopLenZ(BinId)=CountHopLenZ(BinId)+1
          SumHopCount(4)=SumHopCount(4)+1
        endif
      enddo
      close(40,IOSTAT=ierror)
    enddo
    close(39,IOSTAT=ierror)
    call system("rm -rf "//trim(adjustl(WriteHopDir))//" 2>/dev/null")
    call FileSystem%DeleteFileNameList()

    write(FileName1,"(A,A,A)")trim(adjustl(Stat_Dir))//'HopTimePDF_',trim(StatSuffix),'.txt'
    open(41,file=FileName1,status='replace',form='formatted',action='write',IOSTAT=ierror)
    write(41,'(A,I18)')' Total Sample :    ', SumHopCount(1)
    write(41,'(A,I18)')' Total PDF Sample :', SumHopCount(2)
    write(41,'(A,ES18.9)') ' Mean Hop Time (s): ', SumHopStat(1)/real(SumHopCount(1),RK)
    write(41,'(A)')' Time(s) nP PD'
    CountTemp=sum(CountHopTime)
    do pid=1,TimeBinNum
      write(41,'(ES18.9,I18,ES18.9)')TimeBinSet+(real(pid,RK)-0.5_RK)*TimeBinDelta,CountHopTime(pid),real(CountHopTime(pid),RK)/real(CountTemp,RK)
    enddo
    close(41,IOSTAT=ierror)

    write(FileName1,"(A,A,A)")trim(adjustl(Stat_Dir))//'HopLenXPDF_',trim(StatSuffix),'.txt'
    open(41,file=FileName1,status='replace',form='formatted',action='write',IOSTAT=ierror)
    write(41,'(A,I18)')' Total Sample :    ', SumHopCount(1)
    write(41,'(A,I18)')' Total PDF Sample :', SumHopCount(3)
    write(41,'(A,ES18.9)') ' Mean Hop Length (cm): ', SumHopStat(2)/real(SumHopCount(1),RK)*100.0_RK
    write(41,'(A,ES18.9)') ' Mean Abs Length (cm): ', SumHopStat(3)/real(SumHopCount(1),RK)*100.0_RK
    write(41,'(A)')' LenX(cm) nP PD'
    CountTemp=sum(CountHopLenX)
    do pid=1,LenBinNumX
      write(41,'(ES18.9,I18,ES18.9)')(LenBinSetX+(real(pid,RK)-0.5_RK)*LenBinDeltaX)*100.0_RK,CountHopLenX(pid),real(CountHopLenX(pid),RK)/real(CountTemp,RK)
    enddo
    close(41,IOSTAT=ierror)

    write(FileName1,"(A,A,A)")trim(adjustl(Stat_Dir))//'HopLenZPDF_',trim(StatSuffix),'.txt'
    open(41,file=FileName1,status='replace',form='formatted',action='write',IOSTAT=ierror)
    write(41,'(A,I18)')' Total Sample :    ', SumHopCount(1)
    write(41,'(A,I18)')' Total PDF Sample :', SumHopCount(4)
    write(41,'(A,ES18.9)') ' Mean Hop Length (cm): ', SumHopStat(4)/real(SumHopCount(1),RK)*100.0_RK
    write(41,'(A,ES18.9)') ' Mean Abs Length (cm): ', SumHopStat(5)/real(SumHopCount(1),RK)*100.0_RK
    write(41,'(A)')' LenZ(cm) nP PD'
    CountTemp=sum(CountHopLenZ)
    do pid=1,LenBinNumZ
      write(41,'(ES18.9,I18,ES18.9)')(LenBinSetZ+(real(pid,RK)-0.5_RK)*LenBinDeltaZ)*100.0_RK,CountHopLenZ(pid),real(CountHopLenZ(pid),RK)/real(CountTemp,RK)
    enddo
    close(41,IOSTAT=ierror)

  end subroutine ClcHopStat

  !******************************************************************
  ! ClcDisplacement
  !******************************************************************
  function ClcDisplacement(PrtclDataIn,pid) result(disp)
    implicit none
    integer,intent(in)::pid
    type(DumpPrtclVarOut),dimension(:),intent(in)::PrtclDataIn

    ! locals
    real(RK)::disp,dispX,dispZ

    dispX=PrtclDataIn(pid+1)%Pos(1)-PrtclDataIn(pid)%Pos(1)
    dispZ=PrtclDataIn(pid+1)%Pos(3)-PrtclDataIn(pid)%Pos(3)
    !disp=sqrt(dispX*dispX +dispZ*dispZ)
    disp=dispX
  end function ClcDisplacement

  !******************************************************************
  ! ClcYmin
  !******************************************************************
  function ClcYmin(PrtclDataIn,iSet,iEnd) result(yMin)
    implicit none
    integer,intent(in)::iSet,iEnd
    type(DumpPrtclVarOut),dimension(:),intent(in)::PrtclDataIn    

    ! locals
    integer::k
    real(RK)::yMin

    yMin=1.0E20_RK
    do k=iSet,iEnd
      if(PrtclDataIn(k)%Pos(2)<yMin)yMin=PrtclDataIn(k)%Pos(2)
    enddo
  end function ClcYmin

  !******************************************************************
  ! ClcCrossDiffusion
  !******************************************************************
  subroutine ClcCrossDiffusion()
    implicit none
 
    ! locals
    character(128)::FileName1
    type(DumpPrtclVarOut)::DumpVarT
    integer,allocatable,dimension(:)::CrossCount,CrossCountR
    real(RK),allocatable,dimension(:,:)::MeanVariance,MeanVarianceR
    type(DumpPrtclVarOut),allocatable,dimension(:)::PrtclDataIn
    real(RK)::PosX0,PosZ0,MoveX1,MoveX2,MoveZ1,MoveZ2,MoveX,MoveZ,MoveXj
    integer::pid,fid,j,FileByte,SinglePrtclByte,nPrtclSize,ierror,nSample

    allocate(CrossCount(CrossDiffusionBinNum));      CrossCount=0
    allocate(CrossCountR(CrossDiffusionBinNum));     CrossCountR=0
    allocate(MeanVariance(2,CrossDiffusionBinNum));  MeanVariance=0.0_RK
    allocate(MeanVarianceR(2,CrossDiffusionBinNum)); MeanVarianceR=0.0_RK

    SinglePrtclByte=sizeof(DumpVarT)
    if(nrank==0) then
      call FileSystem%WriteFileNameList(ReadDiffuseDir,"pid*")
    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call FileSystem%GetFileName(ReadDiffuseDir)
    do fid=1,FileSystem%nFile
      if(mod(fid,nproc)/=nrank) cycle
      FileName1=trim(adjustl(ReadDiffuseDir))//trim(adjustl(FileSystem%FileName(fid)))

      FileByte=GetFileByte(FileName1)
      if(mod(FileByte,SinglePrtclByte)/=0) then
        print*,"ClcCrossDiffusion: FileByte wrong";stop
      endif
      nPrtclSize=FileByte/SinglePrtclByte
      if(nPrtclSize<1) cycle
      if(iDetectCrossSet+1>=nPrtclSize) cycle

      allocate(PrtclDataIn(nPrtclSize))
      open(42,file=FileName1,status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
      read(42)PrtclDataIn
      close(42)

      !========== Start to Calculate Cross-Diffusion Statistics ==========!
      PosX0=PrtclDataIn(iDetectCrossSet+1)%Pos(1)
      PosZ0=PrtclDataIn(iDetectCrossSet+1)%Pos(3)
      MoveX=PrtclDataIn(nPrtclSize)%Pos(1)- PosX0
      nSample=floor(MoveX/CrossDiffBinDelta)

      pid=iDetectCrossSet+1
      do j=1,nSample
        MoveXj=real(j,RK)*CrossDiffBinDelta
        do
          MoveX1=PrtclDataIn(pid  )%Pos(1)- PosX0
          MoveX2=PrtclDataIn(pid+1)%Pos(1)- PosX0
          if(MoveXj>=MoveX1 .and. MoveXj<MoveX2) then
            MoveZ1=PrtclDataIn(pid  )%Pos(3)- PosZ0
            MoveZ2=PrtclDataIn(pid+1)%Pos(3)- PosZ0
            MoveZ = (MoveZ1*(MoveX2-MoveXj) + MoveZ2*(MoveXj-MoveX1))/(MoveX2-MoveX1)
            CrossCount(j)=CrossCount(j)+1
            MeanVariance(1,j)=MeanVariance(1,j)+MoveZ
            MeanVariance(2,j)=MeanVariance(2,j)+MoveZ*MoveZ            
            exit
          endif
          pid=pid+1
          if(pid>nPrtclSize) exit
        enddo
      enddo

      deallocate(PrtclDataIn)
      write(*,'(A,I4,A)')'rank',nrank,' ClcCrossDiffusion: '//trim(adjustl(FileName1))//' successfully'
    enddo
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(CrossCount,  CrossCountR,   CrossDiffusionBinNum, MPI_INTEGER,         MPI_SUM,0,MPI_COMM_WORLD,ierror)
    call MPI_REDUCE(MeanVariance,MeanVarianceR,2*CrossDiffusionBinNum,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank/=0) return
    call FileSystem%DeleteFileNameList()

    ! Dump the Cross-Diffusion statistics
    write(FileName1,"(A,A,A)")trim(adjustl(Stat_Dir))//'CrossDuffusion_',trim(StatSuffix),'.txt'
    open(43,file=FileName1,status='replace',form='formatted',action='write',IOSTAT=ierror)
    write(43,'(A)')' StreamPos(cm) Mean(cm) Variance(cm^2) nP'
    write(43,'(3ES18.9,I10)')0.0,0.0,0.0,0
    do pid=1,CrossDiffusionBinNum
      nSample=CrossCountR(pid)
      if(nSample>0) then
        write(43,'(3ES18.9,I10)')real(pid,RK)*CrossDiffBinDelta*100_RK,MeanVarianceR(1,pid)*100_RK/real(nSample,RK),MeanVarianceR(2,pid)*10000_RK/real(nSample,RK),nSample
      else
        write(43,'(3ES18.9,I10)')real(pid,RK)*CrossDiffBinDelta*100_RK,0.0,0.0,nSample
      endif
    enddo
    close(37,IOSTAT=ierror)

  end subroutine ClcCrossDiffusion

end module PPD_Statistics
