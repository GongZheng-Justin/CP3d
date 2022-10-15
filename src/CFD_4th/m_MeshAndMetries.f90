module m_MeshAndMetries
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_decomp2d
  use m_Parameters
  implicit none
  private 
  
  real(RK),public::dx,  dy,  dz   ! average mesh intervals in three directions
  real(RK),public::dx2, dy2, dz2  ! square of average mesh intervals in three directions  
  real(RK),public::rdx, rdy, rdz  ! inverse average mesh intervals in three directions
  real(RK),public::rdx2,rdy2,rdz2 ! square of the inverse average mesh intervals in three directions

  real(RK),public,allocatable,dimension(:):: yp   ! point coordinate in y-dir. Suffix 'v' means 'vector'   
  real(RK),public,allocatable,dimension(:):: yc   ! center coordinate in y-dir 
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
  
  ! ux/uz Laplacian metries in y-dir (STAGGERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am2c
  real(RK),public,allocatable,dimension(:):: ac2c
  real(RK),public,allocatable,dimension(:):: ap2c  

  ! uy Laplacian metries in y-dir (CENTERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am2p
  real(RK),public,allocatable,dimension(:):: ac2p  
  real(RK),public,allocatable,dimension(:):: ap2p

  ! for linear interpolation in y-dir
  real(RK),public,allocatable,dimension(:):: YinterpCoe

  ! for 4-order FD
  real(RK),public::dxCoe1,dxCoe2,dxCoe3,dxCoe4
  real(RK),public::dzCoe1,dzCoe2,dzCoe3,dzCoe4
  real(RK),public::dxxCoe1,dxxCoe2,dxxCoe3,dxxCoe4,dxxCoe5
  real(RK),public::dzzCoe1,dzzCoe2,dzzCoe3,dzzCoe4,dzzCoe5
  real(RK),public::InterpCoe1,InterpCoe2,InterpCoe3,InterpCoe4
  
#ifdef ScalarFlow
  ! scalar Laplacian metries in y-dir (STAGGERED VARIABLE)
  real(RK),public,allocatable,dimension(:):: am2cSc
  real(RK),public,allocatable,dimension(:):: ac2cSc
  real(RK),public,allocatable,dimension(:):: ap2cSc
#endif

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
    integer:: j,nUnitFile,ierrTmp,ierror=0
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

    ! for 4-order FD 
    dxCoe1=  rdx/24.00_RK; dxCoe2= -rdx*1.125_RK; dxCoe3=  rdx*1.125_RK; dxCoe4= -rdx/24.00_RK;
    dzCoe1=  rdz/24.00_RK; dzCoe2= -rdz*1.125_RK; dzCoe3=  rdz*1.125_RK; dzCoe4= -rdz/24.00_RK
    dxxCoe1=-rdx2/12.0_RK; dxxCoe2=4.0_RK/3.0_RK*rdx2; dxxCoe3=-2.5_RK*rdx2; dxxCoe4=4.0_RK/3.0_RK*rdx2; dxxCoe5=-rdx2/12.0_RK;
    dzzCoe1=-rdz2/12.0_RK; dzzCoe2=4.0_RK/3.0_RK*rdz2; dzzCoe3=-2.5_RK*rdz2; dzzCoe4=4.0_RK/3.0_RK*rdz2; dzzCoe5=-rdz2/12.0_RK;
    InterpCoe1= -1.0_RK/16.00_RK; InterpCoe2=  9.0_RK/16.00_RK; InterpCoe3=  9.0_RK/16.00_RK; InterpCoe4= -1.0_RK/16.00_RK
    
    allocate(yp(0:nyp),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(yc(0:nyp),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(dyp(0:nyp),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(rdyp(0:nyp), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(dyc(0:nyp),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(rdyc(0:nyp), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)    
    
    allocate(ap2Pr(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac2Pr(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)  
    allocate(am2Pr(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    allocate(ap2c(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac2c(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(am2c(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    allocate(ap2p(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac2p(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)    
    allocate(am2p(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)

    allocate(VolCell(0:nyp),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(DeltaCell(0:nyp), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(YinterpCoe(1:nyp),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
#ifdef ScalarFlow
    allocate(ap2cSc(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(ac2cSc(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(am2cSc(1:nyc), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
#endif
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitMeshAndMetries: ","Allocation failed")
        
    ! Read mesh options and claculate yp for every section
    open(newunit=nUnitFile, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierror )
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
      call clcYcoord(yp(1:nyp),nidYSection(j:j+1),ySectionCoord(j:j+1),StretType(j),SectioncStret(j),StretOption(j))
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

    ! dyp, point coordinate interval in y-dir 
    do j=1,nyc
      dyp(j) = yp(j+1) - yp(j)    
    enddo
    dyp(0)  =dyp(1)
    dyp(nyp)=dyp(nyc)

    !VolCell,volume of the cell
    do j=0,nyp
      VolCell(j)=dyp(j)*dx*dz
    enddo
  
    !DeltaCell, (VolCell)^(1/3)
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
    am2Pr(1)= zero
    ap2Pr(1)= rdyp(1)*rdyc(2) 
    am2Pr(nyc)= rdyp(nyc)*rdyc(nyc)
    ap2Pr(nyc)= zero 
    ac2Pr = -(am2Pr+ap2Pr)
    
    ! ux/uz Laplacian metries in y-dir (STAGGERED VARIABLE)=====================
    do j=1,nyc
      am2c(j)= rdyp(j)*rdyc(j)
      ap2c(j)= rdyp(j)*rdyc(j+1)
    enddo
    am2c(1)= four*rdyc(1)/( dyc(1)+two*dyc(2) )
    ap2c(1)= four*rdyc(2)/( dyc(1)+two*dyc(2) )
    if(FlowType==FT_CH ) then
      am2c(nyc)= four*rdyc(nyc)/( dyc(nyp)+two*dyc(nyc) )
      ap2c(nyc)= four*rdyc(nyp)/( dyc(nyp)+two*dyc(nyc) )
    elseif(FlowType==FT_HC) then
      am2c(nyc)= rdyp(nyc)*rdyc(nyc)
      ap2c(nyc)= zero
    endif
    ac2c= -(am2c+ap2c)

    ! uy Laplacian metries in y-dir (CENTERED VARIABLE)
    do j=1,nyc
      am2p(j)= rdyc(j)*rdyp(j-1)        
      ap2p(j)= rdyc(j)*rdyp(j)        
    enddo
    ac2p= -(am2p+ap2p)

#ifdef ScalarFlow
    ! scalar Laplacian metries in y-dir (STAGGERED VARIABLE)=====================
    do j=1,nyc
      am2cSc(j)= rdyp(j)*rdyc(j)
      ap2cSc(j)= rdyp(j)*rdyc(j+1)
    enddo
    if(ScalarBcOption(1)==-1) then
      am2cSc(1)= four*rdyc(1)/( dyc(1)+two*dyc(2) )
      ap2cSc(1)= four*rdyc(2)/( dyc(1)+two*dyc(2) )
    endif
    if(ScalarBcOption(2)==-1) then
      am2cSc(nyc)= four*rdyc(nyc)/( dyc(nyp)+two*dyc(nyc) )
      ap2cSc(nyc)= four*rdyc(nyp)/( dyc(nyp)+two*dyc(nyc) )
    endif
    ac2cSc= -(am2cSc+ap2cSc)
#endif
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
