module PPD_ReadAndSortDump
  use MPI
  use PPD_Tools
  use PPD_TypeDef
  use PPD_Parameters
  implicit none
  private

  public:: ReadPrtclDump, SortPrtclDump, ShrinkPrtclDump
contains

  !******************************************************************
  ! ReadPrtclDump
  !******************************************************************
  subroutine ReadPrtclDump()
    implicit none

    ! locals
    type(DumpPrtclVar)::DumpVarT
    type(DumpPrtclVarOut)::DumpOutTemp
    character(128)::FileName,FileNameOut
    type(DumpPrtclVar),allocatable,dimension(:)::DumpPrtclIn
    integer::i,pid,FileByte,nDumpPrtclSize,SinglePrtclByte,ierror

    SinglePrtclByte=sizeof(DumpVarT)
    call CreateDir(IsRewritePrtclDump,PrtclDumpOutDir)

    call FileSystem%WriteFileNameList(PrtclDumpDir,PrtclDumpFilter)
    call FileSystem%GetFileName(PrtclDumpDir)
    do i=1,FileSystem%nFile
      FileName=trim(adjustl(PrtclDumpDir))//trim(adjustl(FileSystem%FileName(i)))

      FileByte=GetFileByte(FileName)
      if(mod(FileByte,SinglePrtclByte)/=0) then
        print*,"ReadPrtclDump: FileByte wrong",FileByte;stop
      endif
      nDumpPrtclSize=FileByte/SinglePrtclByte
      if(nDumpPrtclSize<1) cycle
      allocate(DumpPrtclIn(nDumpPrtclSize))
      open(26,file=FileName,status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
      read(26)DumpPrtclIn
      close(26)

      do pid=1,nDumpPrtclSize
        DumpOutTemp=DumpPrtclIn(pid)
        write(FileNameOut,"(A,A,I7.7)")trim(adjustl(PrtclDumpOutDir)),"pid",DumpPrtclIn(pid)%id
        open(27,file=FileNameOut,form='unformatted',access='stream',action='write',position='append',IOSTAT=ierror)
        write(27)DumpOutTemp
        close(27)
      enddo

      deallocate(DumpPrtclIn)
      print*,'read: '//trim(adjustl(FileName))//' successfully'
    enddo
    call FileSystem%DeleteFileNameList()
    
  end subroutine ReadPrtclDump

  !******************************************************************
  ! SortPrtclDump
  !******************************************************************
  subroutine SortPrtclDump()
    implicit none

    ! locals
    character(128)::FileName
    type(DumpPrtclVarOut)::DumpVarT
    real(RK),dimension(:),allocatable::xlenInc,zlenInc
    integer,allocatable,dimension(:)::iTimeVec,iTimeSeq
    integer::i,j,FileByte,nDumpPrtclSize,SinglePrtclByte,ierror
    type(DumpPrtclVarOut),allocatable,dimension(:)::DumpPrtclOut1,DumpPrtclOut2

    SinglePrtclByte=sizeof(DumpVarT)
    if(nrank==0) call FileSystem%WriteFileNameList(PrtclDumpOutDir,"pid*")
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    call FileSystem%GetFileName(PrtclDumpOutDir)
    do i=1,FileSystem%nFile
      if(mod(i,nproc)/=nrank) cycle
      FileName=trim(adjustl(PrtclDumpOutDir))//trim(adjustl(FileSystem%FileName(i)))

      FileByte=GetFileByte(FileName)
      if(mod(FileByte,SinglePrtclByte)/=0) then
        print*,"SortPrtclDump: FileByte wrong";stop
      endif
      nDumpPrtclSize=FileByte/SinglePrtclByte
      if(nDumpPrtclSize<1) cycle

      allocate(iTimeVec(nDumpPrtclSize),iTimeSeq(nDumpPrtclSize))
      allocate(xlenInc(nDumpPrtclSize), zlenInc(nDumpPrtclSize))
      allocate(DumpPrtclOut1(nDumpPrtclSize),DumpPrtclOut2(nDumpPrtclSize))
      open(28,file=FileName,status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
      read(28)DumpPrtclOut1
      close(28)
     
      do j=1,nDumpPrtclSize
        iTimeSeq(j)=j
        iTimeVec(j)=DumpPrtclOut1(j)%iTime
      enddo
      call QuickSort(iTimeVec,iTimeSeq,1,nDumpPrtclSize)
      do j=1,nDumpPrtclSize
        DumpPrtclOut2(j)=DumpPrtclOut1(iTimeSeq(j))
      enddo


      if(IsAddPeriodicLen) then 
        xlenInc=0.0_RK; zlenInc=0.0_RK
        do j=1,nDumpPrtclSize-1
           xlenInc(j+1)=xlenInc(j); zlenInc(j+1)=zlenInc(j)
           if(DumpPrtclOut2(j)%Pos(1)>0.8_RK*xlx .and. DumpPrtclOut2(j+1)%Pos(1)<0.2_RK*xlx)xlenInc(j+1)=xlenInc(j+1)+xlx
           if(DumpPrtclOut2(j)%Pos(1)<0.2_RK*xlx .and. DumpPrtclOut2(j+1)%Pos(1)>0.8_RK*xlx)xlenInc(j+1)=xlenInc(j+1)-xlx
           if(DumpPrtclOut2(j)%Pos(3)>0.8_RK*zlz .and. DumpPrtclOut2(j+1)%Pos(3)<0.2_RK*zlz)zlenInc(j+1)=zlenInc(j+1)+zlz
           if(DumpPrtclOut2(j)%Pos(3)<0.2_RK*zlz .and. DumpPrtclOut2(j+1)%Pos(3)>0.8_RK*zlz)zlenInc(j+1)=zlenInc(j+1)-zlz
        enddo
        do j=1,nDumpPrtclSize
          DumpPrtclOut2(j)%Pos(1)=DumpPrtclOut2(j)%Pos(1)+real(xlenInc(j),RKP)
          DumpPrtclOut2(j)%Pos(3)=DumpPrtclOut2(j)%Pos(3)+real(zlenInc(j),RKP)
        enddo
      endif

      open(29,file=FileName,status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
      write(29)DumpPrtclOut2
      close(29)
      deallocate(iTimeVec,iTimeSeq,xlenInc,zlenInc,DumpPrtclOut1,DumpPrtclOut2)

      write(*,'(A,I4,A)')'rank',nrank,' sort: '//trim(adjustl(FileName))//' successfully'
    enddo

    if(nrank==0) call FileSystem%DeleteFileNameList()
  end subroutine SortPrtclDump

  !******************************************************************
  ! ShrinkPrtclDump
  !******************************************************************
  subroutine ShrinkPrtclDump()
    implicit none

    ! locals
    type(DumpPrtclVarOut)::DumpVarT
    character(128)::FileName1,FileName2
    type(DumpPrtclVarOut),allocatable,dimension(:)::DumpPrtclOut1,DumpPrtclOut2
    integer::i,j,FileByte,nDumpPrtclSize,nDumpPrtclSizeSmall,SinglePrtclByte,ierror
  
    SinglePrtclByte=sizeof(DumpVarT)
    if(nrank==0) then
      call CreateDir(IsRewriteShrink,ShrinkDumpDir)
      call FileSystem%WriteFileNameList(PrtclDumpOutDir,"pid*")
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    call FileSystem%GetFileName(PrtclDumpOutDir)
    do i=1,FileSystem%nFile
      if(mod(i,nproc)/=nrank) cycle
      FileName1=trim(adjustl(PrtclDumpOutDir))//trim(adjustl(FileSystem%FileName(i)))
      FileName2=trim(adjustl(ShrinkDumpDir))//trim(adjustl(FileSystem%FileName(i)))

      FileByte=GetFileByte(FileName1)
      if(mod(FileByte,SinglePrtclByte)/=0) then
        print*,"ShrinkPrtclDump: FileByte wrong";stop
      endif
      nDumpPrtclSize=FileByte/SinglePrtclByte
      nDumpPrtclSizeSmall= (nDumpPrtclSize-ShrinkDumpNeglect)/ShrinkDumpFreq
      if(nDumpPrtclSizeSmall<1) cycle

      allocate(DumpPrtclOut1(nDumpPrtclSize),DumpPrtclOut2(nDumpPrtclSizeSmall))
      open(30,file=FileName1,status='old',form='unformatted',access='stream',action='read',IOSTAT=ierror)
      read(30)DumpPrtclOut1
      close(30)

      do j=1,nDumpPrtclSizeSmall
        DumpPrtclOut2(j)=DumpPrtclOut1(j*ShrinkDumpFreq+ShrinkDumpNeglect)
      enddo
      open(31,file=FileName2,status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
      write(31)DumpPrtclOut2
      close(31)
      deallocate(DumpPrtclOut1,DumpPrtclOut2)

      write(*,'(A,I4,A)')'rank',nrank,' shrink: '//trim(adjustl(FileName1))//' successfully'
    enddo

    if(nrank==0) call FileSystem%DeleteFileNameList()
  end subroutine ShrinkPrtclDump

end module PPD_ReadAndSortDump
