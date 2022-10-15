module m_MeshAndMetries
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  implicit none
  private 
  
  real(RK),public::dx,  dy,  dz   ! average mesh intervals in three directions
  real(RK),public::dx2, dy2, dz2  ! square of average mesh intervals in three directions  
  real(RK),public::rdx, rdy, rdz  ! inverse average mesh intervals in three directions
  real(RK),public::rdx2,rdy2,rdz2 ! square of the inverse average mesh intervals in three directions
#if defined CFDACM
  real(RK),public:: dyUniform,rdyUniform
#endif

  real(RK),public,allocatable,dimension(:):: yp   ! point coordinate in y-dir. Suffix 'v' means 'vector'
  real(RK),public,allocatable,dimension(:):: xc   ! center coordinate in x-dir    
  real(RK),public,allocatable,dimension(:):: yc   ! center coordinate in y-dir 
  real(RK),public,allocatable,dimension(:):: zc   ! center coordinate in z-dir
  real(RK),public,allocatable,dimension(:):: VolCell   ! cell volume
  real(RK),public,allocatable,dimension(:):: DeltaCell ! (dx*dy*dz)^(1/3)

  real(RK),public,allocatable,dimension(:):: dyp  ! point coordinate interval in y-dir
  real(RK),public,allocatable,dimension(:):: dyc  ! center coordinate interval in y-dir
  real(RK),public,allocatable,dimension(:):: rdyp ! inverse point coordinate interval in y-dir
  real(RK),public,allocatable,dimension(:):: rdyc ! inverse center coordinate interval in y-dir

  ! Pressure Laplacian metries in y-dir
  real(RK),public,allocatable,dimension(:):: ap2Pr
  real(RK),public,allocatable,dimension(:):: ac2Pr  
  real(RK),public,allocatable,dimension(:):: am2Pr 

  ! uy/uz Laplacian metries in x-dir (STAGGERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am1c
  real(RK),public,allocatable,dimension(:):: ac1c
  real(RK),public,allocatable,dimension(:):: ap1c  
  
  ! ux/uz Laplacian metries in y-dir (STAGGERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am2c
  real(RK),public,allocatable,dimension(:):: ac2c
  real(RK),public,allocatable,dimension(:):: ap2c
  
  ! uy Laplacian metries in y-dir (CENTERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am2p
  real(RK),public,allocatable,dimension(:):: ac2p  
  real(RK),public,allocatable,dimension(:):: ap2p
  
  ! ux/uy Laplacian metries in z-dir (STAGGERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am3c
  real(RK),public,allocatable,dimension(:):: ac3c
  real(RK),public,allocatable,dimension(:):: ap3c 

  ! for linear interpolation in y-dir
  real(RK),public,allocatable,dimension(:):: YinterpCoe
  
  public:: InitMeshAndMetries
contains

  !******************************************************************
  ! InitMeshAndMetries
  !****************************************************************** 
  subroutine InitMeshAndMetries(ChannelPrm)
    implicit none 
    character(*),intent(in)::ChannelPrm
    
    ! locals
    integer::nSection
    character(128)::chFile
    NAMELIST/MeshSection/nSection
    integer::j,nUnitFile,ierrTmp,ierror=0
    real(RK),allocatable,dimension(:)::SectionLength,SectioncStret,ySectionCoord
    integer,allocatable,dimension(:)::nycSection,StretType,StretOption,nidYSection
    NAMELIST/MeshOptions/SectionLength,SectioncStret,nycSection,StretType,StretOption
    
    dx = xlx/real(nxc,kind=RK)
    dy = yly/real(nyc,kind=RK)
    dz = zlz/real(nzc,kind=RK)
    rdx= real(nxc,kind=RK)/xlx
    rdy= real(nyc,kind=RK)/yly
    rdz= real(nzc,kind=RK)/zlz 
    dx2= dx*dx;  rdx2= rdx*rdx
    dy2= dy*dy;  rdy2= rdy*rdy
    dz2= dz*dz;  rdz2= rdz*rdz
    
    allocate(yp(0:nyp),    Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(yc(0:nyp),    Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(dyp(0:nyp),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(rdyp(0:nyp),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(dyc(0:nyp),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(rdyc(0:nyp),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)    
    
    allocate(ap2Pr(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac2Pr(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)   
    allocate(am2Pr(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    allocate(ap1c(1:nxc),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac1c(1:nxc),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(am1c(1:nxc),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    allocate(ap2c(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac2c(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(am2c(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    allocate(ap2p(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac2p(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)   
    allocate(am2p(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    allocate(ap3c(1:nzc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac3c(1:nzc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(am3c(1:nzc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)

    allocate(xc(0:nxp),        Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(zc(0:nzp),        Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(VolCell(0:nyp),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(DeltaCell(0:nyp), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)

    allocate(YinterpCoe(1:nyp), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Allocation failed")

    ! Read mesh options and claculate yp for every section
    open(newunit=nUnitFile, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "Cannot open file: "//trim(ChannelPrm))
    read(nUnitFile, nml=MeshSection)
    if(nSection<0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nSection Wrong !!!")
    allocate(nycSection(nSection),StretType(nSection),StretOption(nSection))
    allocate(SectionLength(nSection),SectioncStret(nSection))
    allocate(nidYSection(nSection+1),ySectionCoord(nSection+1))
    read(nUnitFile, nml=MeshOptions)
    if(nrank==0) then
      do j=1,nSection
        if(SectionLength(j)<zero)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "SectionLength Wrong !!!")
        if(SectioncStret(j)<zero)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "SectioncStret Wrong !!!")
        if(nycSection(j)<1)      call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nycSection Wrong 1 !!!")
        if(StretType(j)<0  .or. StretType(j)>3 ) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "StretType Wrong !!!")
        if(StretOption(j)<0.or. StretOption(j)>1)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "StretOption Wrong !!!")
      enddo
      if(sum(nycSection) /= nyc) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries", "nycSection Wrong 2 !!!")
      write(MainLog%nUnit, nml=MeshSection)
      write(MainLog%nUnit, nml=MeshOptions)
    endif
    close(nUnitFile,IOSTAT=ierror)
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
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    do j=1,nSection
      call clcCoord(yp(1:nyp),nidYSection(j:j+1),ySectionCoord(j:j+1),StretType(j),SectioncStret(j),StretOption(j))
    enddo
    deallocate(StretType,StretOption,SectioncStret,nidYSection,ySectionCoord)
    write(chFile,"(A)") trim(ResultsDir)//"yMeshFor"//trim(RunName)//".txt"
    if(nrank==0) then
      open(newunit=nUnitFile, file=chfile,status='replace',form='formatted',IOSTAT=ierror)
      if(ierror/=0.and.nrank==0)call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries","Cannot open file: "//trim(chFile))
      do j=1,nyp
        write(nUnitFile,*)j,yp(j)
      enddo
      close(nUnitFile,IOSTAT=ierror)
    endif

    ! yc, center coordinate interval in y-dir
    do j=1,nyc
      yc(j) = half*(yp(j)+yp(j+1))
    enddo
    yc(0)=two*yp(1)-yc(1)
    yc(nyp)=two*yp(nyp)-yc(nyc)

    ! xc,zc, center coordinate interval in x-dir and z-dir
    do j=0,nxp
      xc(j)=dx*(real(j,RK)-half)
    enddo
    do j=0,nzp
      zc(j)=dz*(real(j,RK)-half)
    enddo

    ! dyp, point coordinate interval in y-dir 
    do j=1,nyc
      dyp(j) = yp(j+1) - yp(j)    
    enddo
    dyp(0)  =dyp(1)
    dyp(nyp)=dyp(nyc)
#ifdef CFDACM
    dyUniform=dyp(1)
    rdyUniform=one/dyUniform
#endif

    ! VolCell,volume of the cell
    do j=0,nyp
      VolCell(j)=dyp(j)*dx*dz
    enddo
  
    ! DeltaCell, (VolCell)^(1/3)
    do j=0,nyp
      DeltaCell(j)=(dyp(j)*dx*dz)**(one/three)
    enddo

    ! dyc, center coordinate interval in y-dir
    do j=2,nyc
      dyc(j) = yc(j)-yc(j-1)
    enddo
    dyc(1)=dyp(1)
    dyc(nyp)=dyp(nyc)
    
    ! rdyp and rdyc, the reverse of the dyp and dyc, respectively.
    do j=0, nyp
      rdyp(j) = one/dyp(j)    
    enddo    
    do j=1,nyp
      rdyc(j) = one/dyc(j)    
    enddo

    ! for linear interpolation in y-dir
    do j=1,nyp
      YinterpCoe(j)= dyp(j)/(dyp(j)+dyp(j-1))
    enddo
    
    ! Pressure Laplacian metries in y-dir
    do j=1,nyc
      am2Pr(j)= rdyp(j)*rdyc(j)
      ap2Pr(j)= rdyp(j)*rdyc(j+1)
    enddo
    if(BcOption(ym_dir)/=BC_Period) then
      am2Pr(1)= zero
      ap2Pr(1)= rdyp(1)*rdyc(2) 
    endif
    if(BcOption(yp_dir)/=BC_Period) then
      am2Pr(nyc)= rdyp(nyc)*rdyc(nyc)
      ap2Pr(nyc)= zero 
    endif    
    ac2Pr = -(am2Pr+ap2Pr)
    
    ! uy/uz Laplacian metries in x-dir (STAGGERED VARIABLE)=====================
    am1c = rdx2;  ap1c=rdx2
    if(BcOption(xm_dir)==BC_NoSlip ) then
      am1c(1)= four/three*rdx2
      ap1c(1)= four/three*rdx2
    endif    
    if(BcOption(xp_dir)==BC_NoSlip ) then
      am1c(nxc)= four/three*rdx2
      ap1c(nxc)= four/three*rdx2 
    endif 
    ac1c= -(am1c+ap1c)
    
    ! ux/uz Laplacian metries in y-dir (STAGGERED VARIABLE)=====================
    do j=1,nyc
      am2c(j)= rdyp(j)*rdyc(j)
      ap2c(j)= rdyp(j)*rdyc(j+1)
    enddo
    if(BcOption(ym_dir)==BC_NoSlip ) then
      am2c(1)= four*rdyc(1)/( dyc(1)+two*dyc(2) )
      ap2c(1)= four*rdyc(2)/( dyc(1)+two*dyc(2) )
    endif
    if(BcOption(yp_dir)==BC_NoSlip ) then
      am2c(nyc)= four*rdyc(nyc)/( dyc(nyp)+two*dyc(nyc) )
      ap2c(nyc)= four*rdyc(nyp)/( dyc(nyp)+two*dyc(nyc) )
    endif
    ac2c= -(am2c+ap2c)
    
    ! ux/uy Laplacian metries in z-dir (STAGGERED VARIABLE)=====================
    am3c = rdz2;  ap3c=rdz2;
    if(BcOption(zm_dir)==BC_NoSlip ) then
      am3c(1)= four/three*rdz2
      ap3c(1)= four/three*rdz2
    endif
    if(BcOption(zp_dir)==BC_NoSlip ) then
      am3c(nzc)= four/three*rdz2
      ap3c(nzc)= four/three*rdz2
    endif
    ac3c= -(am3c+ap3c)

    ! uy Laplacian metries in y-dir (CENTERED VARIABLE)
    do j=1,nyc
      am2p(j)= rdyc(j)*rdyp(j-1)        
      ap2p(j)= rdyc(j)*rdyp(j)        
    enddo
    ac2p= -(am2p+ap2p)    
  end subroutine InitMeshAndMetries

  !******************************************************************
  ! clcCoord
  !******************************************************************
  subroutine clcCoord(coordinate,ncoorId,SectionCoord,StretType,cStret,StretOption)
    implicit none
    real(RK),dimension(:),intent(inout)::coordinate
    real(RK),intent(in)::SectionCoord(2),cStret
    integer,intent(in)::ncoorId(2),StretType,StretOption
    
    ! locals
    integer::j,jt,m
    real(RK)::secLen,tstr,xi,setCrd
    
    m=ncoorId(2)-ncoorId(1)
    secLen=SectionCoord(2)-SectionCoord(1)
    coordinate(ncoorId(1))=SectionCoord(1)
    coordinate(ncoorId(2))=SectionCoord(2)
    SELECT CASE(StretType)
    CASE(0) ! Uniform
      do j=0,m
        jt=ncoorId(1)+j
        coordinate(jt)=SectionCoord(1)+real(j,RK)/real(m,RK)*secLen
      enddo
    CASE(1) ! Tangent hyperbolic function
      tstr= tanh(cStret)
      do j=0,m
        xi= real(j,kind=RK)/real(m)               ! For j: [0,m],  xi: [0, 1]
        setCrd= tanh(cStret*(xi-one))/tstr + one  ! For j: [0,m],  setCrd: [0, 1]
        if(StretOption==0) then  ! bottom
          jt=ncoorId(1)+j
          coordinate(jt)=SectionCoord(1) + setCrd*secLen
        else                     ! top
          jt=ncoorId(2)-j
          coordinate(jt)=SectionCoord(2) - setCrd*secLen
        endif        
      enddo
    CASE(2) ! Sine/cosine function
      tstr=sin(half*cStret*PI)
      do j=0,m
        xi= real(j,RK)/real(m,RK)-one             ! For j: [0,m],  xi: [-1,0]
        setCrd= sin(cStret*xi*PI*half)/tstr + one ! For j: [0,m],  setCrd: [0, 1]
        if(StretOption==0) then  ! bottom
          jt=ncoorId(1)+j
          coordinate(jt)=SectionCoord(1) + setCrd*secLen
        else                     ! top
          jt=ncoorId(2)-j
          coordinate(jt)=SectionCoord(2) - setCrd*secLen
        endif
      enddo
    CASE(3) ! Proportional sequence
      if(cStret==one) then
        do j=0,m
          jt=j+ncoorId(1)
          coordinate(jt)=SectionCoord(1)+real(j,RK)/real(m,RK)*secLen
        enddo
      else
        tstr=(cStret**m -one)/(cStret-one)
        do j=0,m-1
          setCrd=(cStret**j)/tstr
          if(StretOption==0) then  ! bottom
            jt=ncoorId(1)+j
            coordinate(jt+1)=coordinate(jt) + setCrd*secLen
          else                     ! top
            jt=ncoorId(2)-j
            coordinate(jt-1)=coordinate(jt) - setCrd*secLen
          endif
        enddo
      endif
    END SELECT
    coordinate(ncoorId(1))=SectionCoord(1)
    coordinate(ncoorId(2))=SectionCoord(2)  
  end subroutine clcCoord

end module m_MeshAndMetries
