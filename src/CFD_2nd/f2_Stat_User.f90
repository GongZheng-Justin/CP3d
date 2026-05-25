#include "definitions_inc.f90"
module f2_Stat_User
!#define BedForm_Along_X
!#define BedForm_Along_Z
#define nStat_User 11
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d
  use f2_Parameters
  use f2_MeshAndMetries
  use f2_Variables,only:mb1
  implicit none

  integer:: nfstime

#ifdef BedForm_Along_X
  real(RK),allocatable,dimension(:,:,:):: SumStat_plane_x
#endif
#ifdef BedForm_Along_Z
  real(RK),allocatable,dimension(:,:,:):: SumStat_plane_z
#endif
  
  public:: InitStatVar_User, clcStat_User
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitStatVar_User
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitStatVar_User(FileStr); implicit none
    character(len=*),intent(in)::FileStr
    
    ! locals
    integer::iUnit, ierr
     
    open(newunit=iUnit, file=FileStr, status='old',form='formatted',IOSTAT=ierr)
    if(ierr/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitStatVar_User","Cannot open file: "//strip(FileStr))
    endif
    close(iUnit, IOSTAT=ierr)
    nfstime=0; 

#ifdef BedForm_Along_X
    allocate(SumStat_plane_x(y1start(2):y1end(2), y1start(3):y1end(3), nStat_User), Stat=ierr)
    SumStat_plane_x = 0.0_RK
    
    BLOCK
    integer::ifld
    integer(kind=MPI_OFFSET_KIND)::disp
    character(len=:),allocatable::XdmfFile
    XdmfFile = strip(ResultsDir_) // 'VisuFor_BedForm_Along_X_' // strip(RunName_) // '.xmf'
    open(newunit=iUnit, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierr)
    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,'InitStatVar_User','Cannot open file: '//strip(XdmfFile))
    call Write_XDMF_Head(iUnit,nyc,nzc,1)
    call Write_XDMF_yGrid(iUnit)       ! y-grid
    call Write_XDMF_zGrid(iUnit)       ! z-grid
    call Write_XDMF_0Grid(iUnit)       ! 0-grid
    call Write_XDMF_Time_Series(iUnit) ! Time series

    ! attribute
    do ifld=saveStat,ilast,saveStat
      disp = int(0, MPI_OFFSET_KIND)
      write(iUnit,'(A,I10.10,A)')'    <Grid Name="T',ifld,'" GridType="Uniform">'
      write(iUnit,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(iUnit,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'      
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'unm_x','BedForm_Along_X_For_') ! 01
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'vnm_x','BedForm_Along_X_For_') ! 02
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'wnm_x','BedForm_Along_X_For_') ! 03
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'pnm_x','BedForm_Along_X_For_') ! 04
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'ucm_x','BedForm_Along_X_For_') ! 05
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'vcm_x','BedForm_Along_X_For_') ! 06
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'wcm_x','BedForm_Along_X_For_') ! 07
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'uuc_x','BedForm_Along_X_For_') ! 08
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'vvc_x','BedForm_Along_X_For_') ! 09
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'wwc_x','BedForm_Along_X_For_') ! 10
      call Write_XDMF_One(iUnit,ifld,nyc,nzc,1,disp,'uvc_x','BedForm_Along_X_For_') ! 11
      write(iUnit,'(A)')'    </Grid>'
    enddo
    call Write_XDMF_Tail(iUnit)                       
    close(unit=iUnit,IOSTAT=ierr)
    ENDBLOCK 
#endif

#ifdef BedForm_Along_Z
    allocate(SumStat_plane_z(y1start(1):y1end(1), y1start(2):y1end(2), nStat_User), Stat=ierr)
    SumStat_plane_z = 0.0_RK

    BLOCK
    integer::ifld
    integer(kind=MPI_OFFSET_KIND)::disp
    character(len=:),allocatable::XdmfFile
    XdmfFile = strip(ResultsDir_) // 'VisuFor_BedForm_Along_Z_' // strip(RunName_) // '.xmf'
    open(newunit=iUnit, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierr)
    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,'InitStatVar_User','Cannot open file: '//strip(XdmfFile))
    call Write_XDMF_Head(iUnit,nxc,nyc,1)
    call Write_XDMF_xGrid(iUnit)       ! x-grid
    call Write_XDMF_yGrid(iUnit)       ! y-grid
    call Write_XDMF_0Grid(iUnit)       ! 0-grid
    call Write_XDMF_Time_Series(iUnit) ! Time series

    ! attribute
    do ifld=saveStat,ilast,saveStat
      disp = int(0, MPI_OFFSET_KIND)
      write(iUnit,'(A,I10.10,A)')'    <Grid Name="T',ifld,'" GridType="Uniform">'
      write(iUnit,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(iUnit,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'      
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'unm_z','BedForm_Along_Z_For_') ! 01
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'vnm_z','BedForm_Along_Z_For_') ! 02
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'wnm_z','BedForm_Along_Z_For_') ! 03
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'pnm_z','BedForm_Along_Z_For_') ! 04
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'ucm_z','BedForm_Along_Z_For_') ! 05
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'vcm_z','BedForm_Along_Z_For_') ! 06
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'wcm_z','BedForm_Along_Z_For_') ! 07
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'uuc_z','BedForm_Along_Z_For_') ! 08
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'vvc_z','BedForm_Along_Z_For_') ! 09
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'wwc_z','BedForm_Along_Z_For_') ! 10
      call Write_XDMF_One(iUnit,ifld,nxc,nyc,1,disp,'uvc_z','BedForm_Along_Z_For_') ! 11
      write(iUnit,'(A)')'    </Grid>'
    enddo
    call Write_XDMF_Tail(iUnit)                       
    close(unit=iUnit,IOSTAT=ierr)    
    ENDBLOCK 
#endif

  CONTAINS
    subroutine Write_XDMF_Head(iUnit,i1,i2,i3); implicit none
      integer,intent(in)::iUnit,i1,i2,i3
      write(iUnit,'(A)') '<?xml version="1.0" ?>'
      write(iUnit,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
      write(iUnit,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
      write(iUnit,'(A)') '<Domain>'    
      ! grid
      write(iUnit,'(A,3I7,A)')'  <Topology name="TOPO" TopologyType="3DRectMesh" Dimensions="',i3,i2,i1,'"/>'
      write(iUnit,'(A)')'  <Geometry name="GEO" GeometryType="VXVYVZ">'
    end subroutine Write_XDMF_Head

    subroutine Write_XDMF_xGrid(iUnit); implicit none
      integer,intent(in)::iUnit
      integer::i
      write(iUnit,'(A,I1,A,I5,A)') '    <DataItem Format="XML" DataType="Float" Precision="',8,'" Endian="Native" Dimensions="',nxc,'">'
      write(iUnit,'(A)',advance='no') '    '
      do i=1,nxc
        if(mod(i,10)==0) then
          write(iUnit,'(ES15.7)') (i-1)*dx+dx*0.5_RK
          if(i<nxc) write(iUnit,'(A)',advance='no') '    '
        else
          write(iUnit,'(ES15.7)',advance='no') (i-1)*dx+dx*0.5_RK
        endif
      enddo
      if(mod(nxc,10) /=0) write(iUnit,*)' '
      write(iUnit,'(A)')'    </DataItem>'
    end subroutine Write_XDMF_xGrid

    subroutine Write_XDMF_yGrid(iUnit); implicit none
      integer,intent(in)::iUnit
      integer::j
      write(iUnit,'(A,I1,A,I5,A)') '    <DataItem Format="XML" DataType="Float" Precision="',8,'" Endian="Native" Dimensions="',nyc,'">'
      write(iUnit,'(A)',advance='no') '    '
      do j=1,nyc
        if(mod(j,10)==0) then
          write(iUnit,'(ES15.7)') yc(j)
          if(j<nyc) write(iUnit,'(A)',advance='no') '    '
        else
          write(iUnit,'(ES15.7)',advance='no') yc(j)
        endif
      enddo
      if(mod(nyc,10) /=0) write(iUnit,*)' '
      write(iUnit,'(A)')'    </DataItem>'
    end subroutine Write_XDMF_yGrid
    
    subroutine Write_XDMF_zGrid(iUnit); implicit none
      integer,intent(in)::iUnit
      integer::k
      write(iUnit,'(A,I1,A,I5,A)') '    <DataItem Format="XML" DataType="Float" Precision="',8,'" Endian="Native" Dimensions="',nzc,'">'
      write(iUnit,'(A)',advance='no') '    '
      do k=1,nzc
        if(mod(k,10)==0) then
          write(iUnit,'(ES15.7)') (k-1)*dz+dz*0.5_RK
          if(k<nzc) write(iUnit,'(A)',advance='no') '    '
        else
          write(iUnit,'(ES15.7)',advance='no') (k-1)*dz+dz*0.5_RK
        endif
      enddo
      if(mod(nzc,10) /=0) write(iUnit,*)' '
      write(iUnit,'(A)')'    </DataItem>'
    end subroutine Write_XDMF_zGrid
    
    subroutine Write_XDMF_0Grid(iUnit); implicit none
      integer,intent(in)::iUnit
      write(iUnit,'(A,I1,A,I5,A)') '    <DataItem Format="XML" DataType="Float" Precision="',8,'" Endian="Native" Dimensions="',1,'">'
      write(iUnit,'(A,E15.7)')'    ',0.0_RK
      write(iUnit,'(A)')'    </DataItem>'
    end subroutine Write_XDMF_0Grid
          
    subroutine Write_XDMF_Time_Series(iUnit); implicit none
      integer,intent(in)::iUnit
      integer::ifld,nfld
      nfld = (ilast - ifirst +1)/saveStat
      write(iUnit,'(A)')'  </Geometry>'
      write(iUnit,'(A)')'  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
      write(iUnit,'(A)')'    <Time TimeType="List">'
      write(iUnit,'(A,I6,A)')'      <DataItem Format="XML" NumberType="Int" Dimensions="',nfld,'">' 
      write(iUnit,'(A)',advance='no')'        '
      do ifld = 1, nfld
        if(mod(ifld,10)==0) then
          write(iUnit,'(I10)') ifld*saveStat +(ifirst-1)
          if(ifld < nfld) write(iUnit,'(A)',advance='no') '        '
        else
          write(iUnit,'(I10)',advance='no') ifld*saveStat +(ifirst-1)
        endif     
      enddo
      if(mod(nfld,10) /=0) write(iUnit,*)' '
      write(iUnit,'(A)') '      </DataItem>'
      write(iUnit,'(A)') '    </Time>'    
    end subroutine Write_XDMF_Time_Series

    subroutine Write_XDMF_Tail(iUnit); implicit none
      integer,intent(in)::iUnit
      write(iUnit,'(A)')'  </Grid>'
      write(iUnit,'(A)')'</Domain>'
      write(iUnit,'(A)')'</Xdmf>'     
    end subroutine Write_XDMF_Tail
  end subroutine InitStatVar_User
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Write_XDMF_One
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Write_XDMF_One(iUnit,ifld,i1,i2,i3,disp,chAttribute,chPrefix); implicit none
    integer,intent(in)::iUnit,ifld,i1,i2,i3
    character(len=*),intent(in)::chAttribute,chPrefix
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

    ! locals
    integer::iprec
    character(len=:),allocatable::chFile
    
    iprec = mytype
    chFile= strip(chPrefix) // strip(RunName_) // int2str(ifld,10) // '.bin'
    write(iUnit,'(A)')'      <Attribute Name="'//strip(chAttribute)//'" Center="Node">'
    write(iUnit,'(A,I1,A,3I7,A,I15,A)')'        <DataItem Format="Binary" DataType="Float" Precision="',iprec, &
       '" Endian="Native" Dimensions="',i3,i2,i1,'" Seek="',disp,'">'
    write(iUnit,'(A)')'          '//strip(chFile)
    write(iUnit,'(A)')'        </DataItem>'
    write(iUnit,'(A)')'      </Attribute>'
    disp = disp +int(i1,8)*int(i2,8)*int(i3,8)*int(iprec,8)
  end subroutine Write_XDMF_One
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clcStat_User
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine clcStat_User(ux,uy,uz,pr); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pr

#ifdef BedForm_Along_X
    ! locals
    integer::ic,jc,kc,ip,jp,kp
    character(len=:),allocatable::filename
    real(RK)::infstime,uxCell,uyCell,uzCell
    
    do kc=y1start(3),y1end(3)
      kp=kc+1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        do ic=y1start(1),y1end(1)
          ip=ic+1
          uxCell= 0.5_RK*(ux(ic,jc,kc)+ux(ip,jc,kc))
          uyCell= 0.5_RK*(uy(ic,jc,kc)+uy(ic,jp,kc))
          uzCell= 0.5_RK*(uz(ic,jc,kc)+uz(ic,jc,kp))
          SumStat_plane_x(jc,kc, 1) = SumStat_plane_x(jc,kc, 1) +ux(ic,jc,kc)
          SumStat_plane_x(jc,kc, 2) = SumStat_plane_x(jc,kc, 2) +uy(ic,jc,kc)
          SumStat_plane_x(jc,kc, 3) = SumStat_plane_x(jc,kc, 3) +uz(ic,jc,kc)
          SumStat_plane_x(jc,kc, 4) = SumStat_plane_x(jc,kc, 4) +pr(ic,jc,kc)
          SumStat_plane_x(jc,kc, 5) = SumStat_plane_x(jc,kc, 5) +uxCell
          SumStat_plane_x(jc,kc, 6) = SumStat_plane_x(jc,kc, 6) +uyCell
          SumStat_plane_x(jc,kc, 7) = SumStat_plane_x(jc,kc, 7) +uzCell
          SumStat_plane_x(jc,kc, 8) = SumStat_plane_x(jc,kc, 8) +uxCell*uxCell
          SumStat_plane_x(jc,kc, 9) = SumStat_plane_x(jc,kc, 9) +uyCell*uyCell
          SumStat_plane_x(jc,kc,10) = SumStat_plane_x(jc,kc,10) +uzCell*uzCell
          SumStat_plane_x(jc,kc,11) = SumStat_plane_x(jc,kc,11) +uxCell*uyCell
        enddo
      enddo
    enddo
    nfstime= nfstime + 1
    if(mod(itime,SaveStat)/=0) return
    
    infstime = 1.0_RK/real(nfstime,RK)/real(nxc,RK)
    SumStat_plane_x  = infstime*SumStat_plane_x
    filename = strip(ResultsDir_) // 'BedForm_Along_X_For_' // strip(RunName_) // int2str(itime,10) // '.bin'
    call decomp_2d_write_plane_reduce(SumStat_plane_x, 2, 1, nStat_User, MPI_SUM, filename)
    SumStat_plane_x=0.0_RK
#endif

#ifdef BedForm_Along_Z  
    ! locals
    integer::ic,jc,kc,ip,jp,kp
    character(len=:),allocatable::filename
    real(RK)::infstime,uxCell,uyCell,uzCell
    
    do kc=y1start(3),y1end(3)
      kp=kc+1
      do jc=y1start(2),y1end(2)
        jp=jc+1
        do ic=y1start(1),y1end(1)
          ip=ic+1
          uxCell= 0.5_RK*(ux(ic,jc,kc)+ux(ip,jc,kc))
          uyCell= 0.5_RK*(uy(ic,jc,kc)+uy(ic,jp,kc))
          uzCell= 0.5_RK*(uz(ic,jc,kc)+uz(ic,jc,kp))
          SumStat_plane_z(ic,jc, 1) = SumStat_plane_z(ic,jc, 1) +ux(ic,jc,kc)
          SumStat_plane_z(ic,jc, 2) = SumStat_plane_z(ic,jc, 2) +uy(ic,jc,kc)
          SumStat_plane_z(ic,jc, 3) = SumStat_plane_z(ic,jc, 3) +uz(ic,jc,kc)
          SumStat_plane_z(ic,jc, 4) = SumStat_plane_z(ic,jc, 4) +pr(ic,jc,kc)
          SumStat_plane_z(ic,jc, 5) = SumStat_plane_z(ic,jc, 5) +uxCell
          SumStat_plane_z(ic,jc, 6) = SumStat_plane_z(ic,jc, 6) +uyCell
          SumStat_plane_z(ic,jc, 7) = SumStat_plane_z(ic,jc, 7) +uzCell
          SumStat_plane_z(ic,jc, 8) = SumStat_plane_z(ic,jc, 8) +uxCell*uxCell
          SumStat_plane_z(ic,jc, 9) = SumStat_plane_z(ic,jc, 9) +uyCell*uyCell
          SumStat_plane_z(ic,jc,10) = SumStat_plane_z(ic,jc,10) +uzCell*uzCell
          SumStat_plane_z(ic,jc,11) = SumStat_plane_z(ic,jc,11) +uxCell*uyCell
        enddo
      enddo
    enddo
    nfstime= nfstime + 1
    if(mod(itime,SaveStat)/=0) return
    
    infstime = 1.0_RK/real(nfstime,RK)/real(nzc,RK)
    SumStat_plane_z  = infstime*SumStat_plane_z
    filename = strip(ResultsDir_) // 'BedForm_Along_Z_For_' //strip(RunName_) // int2str(itime,10) // '.bin'
    call decomp_2d_write_plane_reduce(SumStat_plane_z, 2, 3, nStat_User, MPI_SUM, filename)
    SumStat_plane_z=0.0_RK
#endif
    
    nfstime=0
  end subroutine clcStat_User

end module f2_Stat_User
#ifdef nStat_User
#undef nStat_User
#endif
