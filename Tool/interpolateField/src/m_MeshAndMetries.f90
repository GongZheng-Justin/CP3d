module m_MeshAndMetries
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  implicit none
  private
    
  real(RK),public,allocatable,dimension(:):: xpOld,xpNew,xcOld,xcNew
  real(RK),public,allocatable,dimension(:):: ypOld,ypNew,ycOld,ycNew
  real(RK),public,allocatable,dimension(:):: zpOld,zpNew,zcOld,zcNew
  public:: InitMeshAndMetries
contains

  !******************************************************************
  ! InitMeshAndMetries
  !****************************************************************** 
  subroutine InitMeshAndMetries(interpPrm)
    implicit none 
    character(*),intent(in)::interpPrm
    
    ! locals
    integer::nSection
    character(128)::chFile
    NAMELIST/MeshSectionOld/nSection
    NAMELIST/MeshSectionNew/nSection
    integer:: j,ierrTemp,ierror=0,nUnit,myistat
    real(RK),allocatable,dimension(:)::SectionLength,SectioncStret,ySectionCoord
    integer,allocatable,dimension(:)::nycSection,StretType,StretOption,nidYSection
    NAMELIST/MeshOptionsOld/SectionLength,SectioncStret,nycSection,StretType,StretOption
    NAMELIST/MeshOptionsNew/SectionLength,SectioncStret,nycSection,StretType,StretOption

    ! x mesh =============
    allocate(xpOld(1:nxpOld), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    allocate(xpNew(1:nxpNew), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    allocate(xcOld(0:nxpOld), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    allocate(xcNew(0:nxpNew), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries: ","Allocation failed: x-dir")
    do j=1,nxpOld
      xpOld(j)=real(j-1,RK)/real(nxcOld,RK)
    enddo
    do j=1,nxpNew
      xpNew(j)=real(j-1,RK)/real(nxcNew,RK)
    enddo
    xpOld(1)=zero; xpOld(nxpOld)=one
    xpNew(1)=zero; xpNew(nxpNew)=one

    ! xc, center coordinate interval in x-dir
    do j=1,nxcOld
      xcOld(j) = half*(xpOld(j)+xpOld(j+1))
    enddo
    do j=1,nxcNew
      xcNew(j) = half*(xpNew(j)+xpNew(j+1))
    enddo
    xcOld(0)     =two*xpOld(1)     -xcOld(1)
    xcNew(0)     =two*xpNew(1)     -xcNew(1)
    xcOld(nxpOld)=two*xpOld(nxpOld)-xcOld(nxcOld)
    xcNew(nxpNew)=two*xpNew(nxpNew)-xcNew(nxcNew)

    ! z mesh =============
    allocate(zpOld(1:nzpOld), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    allocate(zpNew(1:nzpNew), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    allocate(zcOld(0:nzpOld), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    allocate(zcNew(0:nzpNew), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries: ","Allocation failed: z-dir")
    do j=1,nzpOld
      zpOld(j)=real(j-1,RK)/real(nzcOld,RK)
    enddo
    do j=1,nzpNew
      zpNew(j)=real(j-1,RK)/real(nzcNew,RK)
    enddo
    zpOld(1)=zero; zpOld(nzpOld)=one
    zpNew(1)=zero; zpNew(nzpNew)=one

    ! zc, center coordinate interval in z-dir
    do j=1,nzcOld
      zcOld(j) = half*(zpOld(j)+zpOld(j+1))
    enddo
    do j=1,nzcNew
      zcNew(j) = half*(zpNew(j)+zpNew(j+1))
    enddo
    zcOld(0)     =two*zpOld(1)     -zcOld(1)
    zcNew(0)     =two*zpNew(1)     -zcNew(1)
    zcOld(nzpOld)=two*zpOld(nzpOld)-zcOld(nzcOld)
    zcNew(nzpNew)=two*zpNew(nzpNew)-zcNew(nzcNew)
 
    ! y mesh =============
    allocate(ypOld(1:nypOld), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    allocate(ypNew(1:nypNew), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    allocate(ycOld(0:nypOld), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    allocate(ycNew(0:nypNew), Stat=ierrTemp); ierror=ierror+abs(ierrTemp)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries: ","Allocation failed: y-dir")
        
    ! Read mesh options and claculate ypOld for every section 
    open(newunit=nUnit, file=interpPrm, status='old',form='formatted',IOSTAT=myistat )
    if(myistat/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "Cannot open file: "//trim(interpPrm))
    read(nUnit, nml=MeshSectionOld)
    if(nSection<0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nSection Wrong !!!")
    allocate(nycSection(nSection),StretType(nSection),StretOption(nSection))
    allocate(SectionLength(nSection),SectioncStret(nSection))
    allocate(nidYSection(nSection+1),ySectionCoord(nSection+1))
    read(nUnit, nml=MeshOptionsOld)
    if(nrank==0) then
      do j=1,nSection
        if(SectionLength(j)<zero)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "SectionLength Wrong !!!")
        if(SectioncStret(j)<zero)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "SectioncStret Wrong !!!")
        if(nycSection(j)<1)      call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nycSection Wrong 1 !!!")
        if(StretType(j)<0  .or. StretType(j)>3 ) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "StretType Wrong !!!")
        if(StretOption(j)<0.or. StretOption(j)>1)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "StretOption Wrong !!!")
      enddo
      if(sum(nycSection) /= nycOld) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nycSection Wrong 2 !!!")
      write(MainLog%nUnit, nml=MeshSectionOld)
      write(MainLog%nUnit, nml=MeshOptionsOld)
    endif
    close(nUnit,IOSTAT=myistat)
    nidYSection=1
    ySectionCoord=zero
    Block 
      real(RK)::SumLength
      SumLength=sum(SectionLength)
      do j=1,nSection
        nidYSection(j+1)  =nidYSection(j)   + nycSection(j)
        ySectionCoord(j+1)=ySectionCoord(j) + SectionLength(j)/SumLength*yly
      enddo
    End Block
    ySectionCoord(nSection+1)=yly
    deallocate(nycSection,SectionLength)
    call MPI_BARRIER(MPI_COMM_WORLD,myistat)
    do j=1,nSection
      call clcYcoord(ypOld(1:nypOld),nidYSection(j:j+1),ySectionCoord(j:j+1),StretType(j),SectioncStret(j),StretOption(j))
    enddo
    deallocate(StretType,StretOption,SectioncStret,nidYSection,ySectionCoord)
    write(chFile,"(A)") trim(ResultsDir)//"yMeshOldFor"//trim(RunName)//".txt"
    if(nrank==0) then
      open(newunit=nUnit, file=chfile,status='replace',form='formatted',IOSTAT=myistat)
      if(myistat/=0.and.nrank==0)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Cannot open file: "//trim(chFile))
      do j=1,nypOld
        write(nUnit,*)j,ypOld(j)
      enddo
      close(nUnit,IOSTAT=myistat)
    endif

    ! Read mesh options and claculate ypNew for every section
    open(newunit=nUnit, file=interpPrm, status='old',form='formatted',IOSTAT=myistat )
    if(myistat/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "Cannot open file: "//trim(interpPrm))
    read(nUnit, nml=MeshSectionNew)
    if(nSection<0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nSection Wrong !!!")
    allocate(nycSection(nSection),StretType(nSection),StretOption(nSection))
    allocate(SectionLength(nSection),SectioncStret(nSection))
    allocate(nidYSection(nSection+1),ySectionCoord(nSection+1))
    read(nUnit, nml=MeshOptionsNew)
    if(nrank==0) then
      do j=1,nSection
        if(SectionLength(j)<zero)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "SectionLength Wrong !!!")
        if(SectioncStret(j)<zero)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "SectioncStret Wrong !!!")
        if(nycSection(j)<1)      call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nycSection Wrong 1 !!!")
        if(StretType(j)<0  .or. StretType(j)>3 ) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "StretType Wrong !!!")
        if(StretOption(j)<0.or. StretOption(j)>1)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "StretOption Wrong !!!")
      enddo
      if(sum(nycSection) /= nycNew) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nycSection Wrong 2 !!!")
      write(MainLog%nUnit, nml=MeshSectionNew)
      write(MainLog%nUnit, nml=MeshOptionsNew)
    endif
    close(nUnit,IOSTAT=myistat)
    nidYSection=1
    ySectionCoord=zero
    Block 
      real(RK)::SumLength
      SumLength=sum(SectionLength)
      do j=1,nSection
        nidYSection(j+1)  =nidYSection(j)   + nycSection(j)
        ySectionCoord(j+1)=ySectionCoord(j) + SectionLength(j)/SumLength*yly
      enddo
    End Block
    ySectionCoord(nSection+1)=yly
    deallocate(nycSection,SectionLength)
    call MPI_BARRIER(MPI_COMM_WORLD,myistat)
    do j=1,nSection
      call clcYcoord(ypNew(1:nypNew),nidYSection(j:j+1),ySectionCoord(j:j+1),StretType(j),SectioncStret(j),StretOption(j))
    enddo
    deallocate(StretType,StretOption,SectioncStret,nidYSection,ySectionCoord)
    write(chFile,"(A)") trim(ResultsDir)//"yMeshNewFor"//trim(RunName)//".txt"
    if(nrank==0) then
      open(newunit=nUnit, file=chfile,status='replace',form='formatted',IOSTAT=myistat)
      if(myistat/=0.and.nrank==0)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Cannot open file: "//trim(chFile))
      do j=1,nypNew
        write(nUnit,*)j,ypNew(j)
      enddo
      close(nUnit,IOSTAT=myistat)
    endif
    
    ! yc, center coordinate interval in y-dir
    do j=1,nycOld
      ycOld(j) = half*(ypOld(j)+ypOld(j+1))
    enddo
    do j=1,nycNew
      ycNew(j) = half*(ypNew(j)+ypNew(j+1))
    enddo
    ycOld(0)     =two*ypOld(1)     -ycOld(1)
    ycNew(0)     =two*ypNew(1)     -ycNew(1)
    ycOld(nypOld)=two*ypOld(nypOld)-ycOld(nycOld)
    ycNew(nypNew)=two*ypNew(nypNew)-ycNew(nycNew)
  end subroutine InitMeshAndMetries

  !******************************************************************
  ! clcYcoord
  !******************************************************************
  subroutine clcYcoord(ycoord,nidY,ySectionCoord,StretType,cStret,StretOption)
    implicit none
    real(RK),dimension(:),intent(inout)::ycoord
    real(RK),intent(in)::ySectionCoord(2),cStret
    integer,intent(in)::nidY(2),StretType,StretOption
    
    ! locals
    integer::j,jt,m
    real(RK)::yLen,tstr,xi,y0
    
    m=nidY(2)-nidY(1)
    yLen=ySectionCoord(2)-ySectionCoord(1)
    ycoord(nidY(1))=ySectionCoord(1)
    ycoord(nidY(2))=ySectionCoord(2)
    SELECT CASE(StretType)
    CASE(0) ! Uniform
      do j=0,m
        jt=nidY(1)+j
        ycoord(jt)=ySectionCoord(1)+real(j,RK)/real(m,RK)*yLen
      enddo
    CASE(1) ! Tangent hyperbolic function
      tstr= tanh(cStret)
      do j=0,m
        xi= real(j,kind=RK)/real(m)           ! For j: [0,m],  xi: [0, 1]
        y0= tanh(cStret*(xi-one))/tstr + one  ! For j: [0,m],  y0: [0, 1]
        if(StretOption==0) then  ! bottom
          jt=nidY(1)+j
          ycoord(jt)=ySectionCoord(1) + y0*yLen
        else                     ! top
          jt=nidY(2)-j
          ycoord(jt)=ySectionCoord(2) - y0*yLen
        endif        
      enddo
    CASE(2) ! Sine/cosine function
      tstr=sin(half*cStret*PI)
      do j=0,m
        xi= real(j,RK)/real(m,RK)-one         ! For j: [0,m],  xi: [-1,0]
        y0= sin(cStret*xi*PI*half)/tstr + one ! For j: [0,m],  y0: [0, 1]
        if(StretOption==0) then  ! bottom
          jt=nidY(1)+j
          ycoord(jt)=ySectionCoord(1) + y0*yLen
        else                     ! top
          jt=nidY(2)-j
          ycoord(jt)=ySectionCoord(2) - y0*yLen
        endif
      enddo
    CASE(3) ! Proportional sequence
      if(cStret==one) then
        do j=0,m
          jt=j+nidY(1)
          ycoord(jt)=ySectionCoord(1)+real(j,RK)/real(m,RK)*yLen
        enddo
      else
        tstr=(cStret**m -one)/(cStret-one)
        do j=0,m-1
          y0=(cStret**j)/tstr
          if(StretOption==0) then  ! bottom
            jt=nidY(1)+j
            ycoord(jt+1)=ycoord(jt) + y0*yLen
          else                     ! top
            jt=nidY(2)-j
            ycoord(jt-1)=ycoord(jt) - y0*yLen
          endif
        enddo
      endif
    END SELECT
    ycoord(nidY(1))=ySectionCoord(1)
    ycoord(nidY(2))=ySectionCoord(2)  
  end subroutine clcYcoord

end module m_MeshAndMetries
