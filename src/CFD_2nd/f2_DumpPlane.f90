#include "definitions_inc.f90"
module f2_DumpPlane
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d
  use f2_Parameters
  use f2_Variables,only: mb1
  use mc_FileOperator, only: mkdir
  use f2_MeshAndMetries
  implicit none
  private

#define RKP_Dump 4

  integer :: nxPlane = 0
  integer :: nyPlane = 0
  integer :: nzPlane = 0
  integer :: nTimePlane, iTimePlane

  logical :: DumpPlaneFlag = .false.          ! Output plane data or not
  character(len=:),allocatable:: PlaneDir_    ! Output plane directory

  integer :: WritePlaneFileFreq = 100000000   ! Output plane frequency
  integer :: RecordPlaneFreq = 100000000      ! Record plane data frequency
  
  logical :: save_ux = .false.
  logical :: save_uy = .false.
  logical :: save_uz = .false.
  logical :: save_pr = .false.

  integer, allocatable :: ixPlane(:)
  integer, allocatable :: iyPlane(:)
  integer, allocatable :: izPlane(:)
  real(RKP_Dump), dimension(:, :, :, :), allocatable :: uxPlane_x, uyPlane_x, uzPlane_x, prPlane_x
  real(RKP_Dump), dimension(:, :, :, :), allocatable :: uxPlane_y, uyPlane_y, uzPlane_y, prPlane_y
  real(RKP_Dump), dimension(:, :, :, :), allocatable :: uxPlane_z, uyPlane_z, uzPlane_z, prPlane_z


  public :: Initialize_DumpPlane, dump_plane


contains


  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Initialize_DumpPlane
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Initialize_DumpPlane(chFile); implicit none
     character(len = *), intent(in) :: chFile

    ! locals
    character(512):: PlaneDir
    integer(8) :: disp, disp_max
    integer, allocatable :: iPlane_tmp(:)
    character(len = :), allocatable :: XdmfFile, attri_str
    integer :: ierrTmp, ierr, iunit, ijk, isum, idtmp, i_dir, len_str, ifld, itime_tail, nd1, nd2, nd3

    namelist /DumpPlane_x_Options/ ixPlane
    namelist /DumpPlane_y_Options/ iyPlane
    namelist /DumpPlane_z_Options/ izPlane
    namelist /DumpPlane_Options/ DumpPlaneFlag, PlaneDir, RecordPlaneFreq, WritePlaneFileFreq, &
            nxPlane, nyPlane, nzPlane, save_ux, save_uy, save_uz, save_pr
    
    open(newunit=iUnit, file=chFile, status='old',form='formatted', iostat = ierr)
    if(ierr/=0 .and. nrank == 0) call MainLog%CheckForError(ErrT_Abort, "Initialize_DumpPlane", "Cannot openfile: " // chFile)
    read(iUnit, nml = DumpPlane_Options, iostat = ierr)
    close(iUnit, iostat = ierrTmp)
    if (ierr /= 0) return
    if (nxPlane < 1 .and. nyPlane < 1 .and. nzPlane < 1) DumpPlaneFlag = .false.
    if ((.not. save_ux) .and. (.not. save_uy) .and. (.not. save_uz) .and. (.not. save_pr)) DumpPlaneFlag = .false.
    if(.not. DumpPlaneFlag) return

    PlaneDir_ = strip(PlaneDir);  len_str = len(PlaneDir_)
    if(PlaneDir_(len_str : len_str) /= '/') PlaneDir_ = PlaneDir_ // '/'
    if (nrank == 0) then
       write(MainLog%iUnit, nml = DumpPlane_Options)
       call mkdir(PlaneDir_, ierr); if(ierr<0) then; print*,'Cannot create folder PlaneDir'; stop; endif
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    if(.not. DumpPlaneFlag) return
    if(mod(WritePlaneFileFreq, RecordPlaneFreq) /= 0 .and. nrank == 0) then
      call MainLog%CheckForError(ErrT_Abort,"Initialize_DumpPlane","nTimePlane WRONG")
    endif
    nTimePlane = WritePlaneFileFreq / RecordPlaneFreq;  iTimePlane = 0

    ! nxPlane
    i_dir = 1;  nd1 = 1;  nd2 = nyc;  nd3 = nzc
    if (nxPlane > 0) then
       allocate( ixPlane(nxPlane), stat = ierr, source = 0)
       if (ierr /= 0) call MainLog%CheckForError(ErrT_Abort, "Initialize_DumpPlane", 'Allocation failed-1')
       open(newunit=iUnit, file=chFile, status='old',form='formatted', iostat = ierr)
       read(iUnit, nml = DumpPlane_x_Options)
       close(iUnit, iostat = ierr)
       call move_alloc(ixPlane, iPlane_tmp)

       ! nxPlane xdmf =====================
       if (nrank == 0) then
          XdmfFile = strip(PlaneDir_) // 'xPlane_' // strip(RunName_) // '.xmf'
          open(newunit=iUnit, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierr)
          if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,'InitStatVar_User','Cannot open file: '//strip(XdmfFile))
          call Write_XDMF_Head(iUnit, nd1, nd2, nd3)
          call Write_XDMF_0Grid(iUnit)       ! 0-grid
          call Write_XDMF_yGrid(iUnit)       ! y-grid
          call Write_XDMF_zGrid(iUnit)       ! z-grid
          call Write_XDMF_Time_Series(iUnit) ! Time series

          ! mesh attribute
          disp = int(0, MPI_OFFSET_KIND)
          disp_max = int(WritePlaneFileFreq / RecordPlaneFreq, 8) * int(nd1,8)*int(nd2,8)*int(nd3,8)*int(RKP_Dump, 8)
          do ifld = RecordPlaneFreq, ilast, RecordPlaneFreq
             write(iUnit,'(A,I10.10,A)')'    <Grid Name="T',ifld,'" GridType="Uniform">'
             write(iUnit,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
             write(iUnit,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
             if (mod(ifld, WritePlaneFileFreq) /= 0) then 
                itime_tail = min((ifld / WritePlaneFileFreq + 1) * WritePlaneFileFreq, ilast)
             else
                itime_tail = min((ifld / WritePlaneFileFreq ) * WritePlaneFileFreq, ilast)
             endif
             do ijk = 1, nxPlane
                idtmp = iPlane_tmp(ijk)

                if (save_ux) then
                   attri_str = 'ux_xPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 01
                endif

                if (save_uy) then
                   attri_str = 'uy_xPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 02
                endif

                if (save_uz) then
                   attri_str = 'uz_xPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 03
                endif

                if (save_pr) then
                   attri_str = 'pr_xPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 04
                endif
             enddo
             disp = disp +int(nd1,8)*int(nd2,8)*int(nd3,8)*int(RKP_Dump, 8)
             if (disp == disp_max) disp = int(0, MPI_OFFSET_KIND)
             write(iUnit,'(A)')'    </Grid>'
          enddo
          call Write_XDMF_Tail(iUnit)                       
          close(unit=iUnit,IOSTAT=ierr)
       endif
       ! nxPlane xdmf =====================

       do ijk = 1, nxPlane
          idtmp = iPlane_tmp(ijk)
          if (idtmp > nxc .or. idtmp < 1 .and. nrank == 0) then
             call MainLog%CheckForError(ErrT_Abort, "Initialize_DumpPlane", 'ixPlane wrong')
          endif
       enddo

       isum = 0
       do ijk = 1, nxPlane
          idtmp = iPlane_tmp(ijk)
          if (idtmp >= y1start(i_dir) .and. idtmp <= y1end(i_dir)) isum = isum + 1
       enddo
       if (isum > 0) then
          allocate( ixPlane(isum), stat = ierr, source = 0)
          if (ierr /= 0) call MainLog%CheckForError(ErrT_Abort, "Initialize_DumpPlane", 'Allocation failed-2')
          isum = 0
          do ijk = 1, nxPlane
             idtmp = iPlane_tmp(ijk)    
             if (idtmp >= y1start(i_dir) .and. idtmp <= y1end(i_dir)) then
                isum = isum + 1;  ixPlane(isum) = idtmp
             endif
         enddo
       endif
       nxPlane = isum
       deallocate(iPlane_tmp, stat = ierr)
    endif

    ! nyPlane
    i_dir = 2;  nd1 = nxc;  nd2 = 1;  nd3 = nzc
    if (nyPlane > 0) then
       allocate( iyPlane(nyPlane), stat = ierr, source = 0)
       if (ierr /= 0) call MainLog%CheckForError(ErrT_Abort, "Initialize_DumpPlane", 'Allocation failed-3')
       open(newunit=iUnit, file=chFile, status='old',form='formatted', iostat = ierr)
       read(iUnit, nml = DumpPlane_y_Options)
       close(iUnit, iostat = ierr)
       call move_alloc(iyPlane, iPlane_tmp)

       ! nyPlane xdmf =====================
       if (nrank == 0) then
          XdmfFile = strip(PlaneDir_) // 'yPlane_' // strip(RunName_) // '.xmf'
          open(newunit=iUnit, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierr)
          if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,'InitStatVar_User','Cannot open file: '//strip(XdmfFile))
          call Write_XDMF_Head(iUnit, nd1, nd2, nd3)
          call Write_XDMF_xGrid(iUnit)       ! x-grid
          call Write_XDMF_0Grid(iUnit)       ! 0-grid
          call Write_XDMF_zGrid(iUnit)       ! z-grid
          call Write_XDMF_Time_Series(iUnit) ! Time series

          ! mesh attribute
          disp = int(0, MPI_OFFSET_KIND)
          disp_max = int(WritePlaneFileFreq / RecordPlaneFreq, 8) * int(nd1,8)*int(nd2,8)*int(nd3,8)*int(RKP_Dump, 8)
          do ifld = RecordPlaneFreq, ilast, RecordPlaneFreq
             write(iUnit,'(A,I10.10,A)')'    <Grid Name="T',ifld,'" GridType="Uniform">'
             write(iUnit,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
             write(iUnit,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
             if (mod(ifld, WritePlaneFileFreq) /= 0) then 
                itime_tail = min((ifld / WritePlaneFileFreq + 1) * WritePlaneFileFreq, ilast)
             else
                itime_tail = min((ifld / WritePlaneFileFreq ) * WritePlaneFileFreq, ilast)
             endif
             do ijk = 1, nyPlane
                idtmp = iPlane_tmp(ijk)

                if (save_ux) then
                   attri_str = 'ux_yPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 01
                endif

                if (save_uy) then
                   attri_str = 'uy_yPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 02
                endif

                if (save_uz) then
                   attri_str = 'uz_yPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 03
                endif

                if (save_pr) then
                   attri_str = 'pr_yPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 04
                endif
             enddo
             disp = disp +int(nd1,8)*int(nd2,8)*int(nd3,8)*int(RKP_Dump, 8)
             if (disp == disp_max) disp = int(0, MPI_OFFSET_KIND)
             write(iUnit,'(A)')'    </Grid>'
          enddo
          call Write_XDMF_Tail(iUnit)                       
          close(unit=iUnit,IOSTAT=ierr)
       endif
       ! nyPlane xdmf =====================

       do ijk = 1, nyPlane
          idtmp = iPlane_tmp(ijk)
          if (idtmp > nyc .or. idtmp < 1 .and. nrank == 0) then
             call MainLog%CheckForError(ErrT_Abort, "Initialize_DumpPlane", 'iyPlane wrong')
          endif
       enddo

       isum = 0
       do ijk = 1, nyPlane
          idtmp = iPlane_tmp(ijk)
          if (idtmp >= y1start(i_dir) .and. idtmp <= y1end(i_dir)) isum = isum + 1
       enddo
       if (isum > 0) then
          allocate( iyPlane(isum), stat = ierr, source = 0)
          if (ierr /= 0) call MainLog%CheckForError(ErrT_Abort, "Initialize_DumpPlane", 'Allocation failed-4')
          isum = 0
          do ijk = 1, nyPlane
             idtmp = iPlane_tmp(ijk)    
             if (idtmp >= y1start(i_dir) .and. idtmp <= y1end(i_dir)) then
                isum = isum + 1;  iyPlane(isum) = idtmp
             endif
         enddo
       endif
       nyPlane = isum
       deallocate(iPlane_tmp, stat = ierr)
    endif

    ! nzPlane
    i_dir = 3;  nd1 = nxc;  nd2 = nyc;  nd3 = 1
    if (nzPlane > 0) then
       allocate( izPlane(nzPlane), stat = ierr, source = 0)
       if (ierr /= 0) call MainLog%CheckForError(ErrT_Abort, "Initialize_DumpPlane", 'Allocation failed-5')
       open(newunit=iUnit, file=chFile, status='old',form='formatted', iostat = ierr)
       read(iUnit, nml = DumpPlane_z_Options)
       close(iUnit, iostat = ierr)
       call move_alloc(izPlane, iPlane_tmp)

       ! nzPlane xdmf =====================
       if (nrank == 0) then
          XdmfFile = strip(PlaneDir_) // 'zPlane_' // strip(RunName_) // '.xmf'
          open(newunit=iUnit, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierr)
          if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,'InitStatVar_User','Cannot open file: '//strip(XdmfFile))
          call Write_XDMF_Head(iUnit, nd1, nd2, nd3)
          call Write_XDMF_xGrid(iUnit)       ! x-grid
          call Write_XDMF_yGrid(iUnit)       ! y-grid
          call Write_XDMF_0Grid(iUnit)       ! 0-grid
          call Write_XDMF_Time_Series(iUnit) ! Time series

          ! mesh attribute
          disp = int(0, MPI_OFFSET_KIND)
          disp_max = int(WritePlaneFileFreq / RecordPlaneFreq, 8) * int(nd1,8)*int(nd2,8)*int(nd3,8)*int(RKP_Dump, 8)
          do ifld = RecordPlaneFreq, ilast, RecordPlaneFreq
             write(iUnit,'(A,I10.10,A)')'    <Grid Name="T',ifld,'" GridType="Uniform">'
             write(iUnit,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
             write(iUnit,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
             if (mod(ifld, WritePlaneFileFreq) /= 0) then 
                itime_tail = min((ifld / WritePlaneFileFreq + 1) * WritePlaneFileFreq, ilast)
             else
                itime_tail = min((ifld / WritePlaneFileFreq ) * WritePlaneFileFreq, ilast)
             endif
             do ijk = 1, nzPlane
                idtmp = iPlane_tmp(ijk)

                if (save_ux) then
                   attri_str = 'ux_zPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 01
                endif

                if (save_uy) then
                   attri_str = 'uy_zPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 02
                endif

                if (save_uz) then
                   attri_str = 'uz_zPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 03
                endif

                if (save_pr) then
                   attri_str = 'pr_zPlane' // int2str(idtmp, 5);  XdmfFile = attri_str // '_' // int2str(itime_tail, 10)                
                   call Write_XDMF_One(iUnit, RKP_Dump, nd1, nd2, nd3, disp, attri_str, XdmfFile) ! 04
                endif
             enddo
             disp = disp +int(nd1,8)*int(nd2,8)*int(nd3,8)*int(RKP_Dump, 8)
             if (disp == disp_max) disp = int(0, MPI_OFFSET_KIND)
             write(iUnit,'(A)')'    </Grid>'
          enddo
          call Write_XDMF_Tail(iUnit)                       
          close(unit=iUnit,IOSTAT=ierr)
       endif
       ! nzPlane xdmf =====================

       do ijk = 1, nzPlane
          idtmp = iPlane_tmp(ijk)
          if (idtmp > nzc .or. idtmp < 1 .and. nrank == 0) then
             call MainLog%CheckForError(ErrT_Abort, "Initialize_DumpPlane", 'izPlane wrong')
          endif
       enddo

       isum = 0
       do ijk = 1, nzPlane
          idtmp = iPlane_tmp(ijk)
          if (idtmp >= y1start(i_dir) .and. idtmp <= y1end(i_dir)) isum = isum + 1
       enddo
       if (isum > 0) then
          allocate( izPlane(isum), stat = ierr, source = 0)
          if (ierr /= 0) call MainLog%CheckForError(ErrT_Abort, "Initialize_DumpPlane", 'Allocation failed-6')
          isum = 0
          do ijk = 1, nzPlane
             idtmp = iPlane_tmp(ijk)    
             if (idtmp >= y1start(i_dir) .and. idtmp <= y1end(i_dir)) then
                isum = isum + 1;  izPlane(isum) = idtmp
             endif
         enddo
       endif
       nzPlane = isum
       deallocate(iPlane_tmp, stat = ierr)
    endif
    if (nxPlane < 1 .and. nyPlane < 1 .and. nzPlane < 1) DumpPlaneFlag = .false.
    if(.not. DumpPlaneFlag) return


    if (nxPlane > 0) then
       ierr = 0
       if (save_ux) allocate(uxPlane_x(y1start(2):y1end(2), y1start(3):y1end(3), nTimePlane, nxPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if (save_uy) allocate(uyPlane_x(y1start(2):y1end(2), y1start(3):y1end(3), nTimePlane, nxPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if (save_uz) allocate(uzPlane_x(y1start(2):y1end(2), y1start(3):y1end(3), nTimePlane, nxPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if (save_pr) allocate(prPlane_x(y1start(2):y1end(2), y1start(3):y1end(3), nTimePlane, nxPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if(ierr /= 0) call MainLog%CheckForError(ErrT_Abort,"Initialize_DumpPlane","Allocation failed-7")
    endif
    if (nyPlane > 0) then
       ierr = 0
       if (save_ux) allocate(uxPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane, nyPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if (save_uy) allocate(uyPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane, nyPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if (save_uz) allocate(uzPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane, nyPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if (save_pr) allocate(prPlane_y(y1start(1):y1end(1), y1start(3):y1end(3), nTimePlane, nyPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp) 
       if(ierr /= 0) call MainLog%CheckForError(ErrT_Abort,"Initialize_DumpPlane","Allocation failed-8")
    endif
    if (nzPlane > 0) then
       ierr = 0
       if (save_ux) allocate(uxPlane_z(y1start(1):y1end(1), y1start(2):y1end(2), nTimePlane, nzPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if (save_uy) allocate(uyPlane_z(y1start(1):y1end(1), y1start(2):y1end(2), nTimePlane, nzPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if (save_uz) allocate(uzPlane_z(y1start(1):y1end(1), y1start(2):y1end(2), nTimePlane, nzPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if (save_pr) allocate(prPlane_z(y1start(1):y1end(1), y1start(2):y1end(2), nTimePlane, nzPlane),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
       if(ierr /= 0) call MainLog%CheckForError(ErrT_Abort,"Initialize_DumpPlane","Allocation failed-9")
    endif

  contains

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
      nfld = (ilast - ifirst +1)/RecordPlaneFreq
      write(iUnit,'(A)')'  </Geometry>'
      write(iUnit,'(A)')'  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
      write(iUnit,'(A)')'    <Time TimeType="List">'
      write(iUnit,'(A,I6,A)')'      <DataItem Format="XML" NumberType="Int" Dimensions="',nfld,'">' 
      write(iUnit,'(A)',advance='no')'        '
      do ifld = 1, nfld
        if(mod(ifld,10)==0) then
          write(iUnit,'(I10)') ifld*RecordPlaneFreq +(ifirst-1)
          if(ifld < nfld) write(iUnit,'(A)',advance='no') '        '
        else
          write(iUnit,'(I10)',advance='no') ifld*RecordPlaneFreq +(ifirst-1)
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

  endsubroutine Initialize_DumpPlane


  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Write_XDMF_One
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Write_XDMF_One(iUnit,iprec, i1,i2,i3,disp,chAttribute,chFile); implicit none
    integer,intent(in)::iUnit, iprec, i1, i2, i3
    character(len=*),intent(in)::chAttribute,chFile
    integer(kind=8),intent(in)::disp

    write(iUnit,'(A)')'      <Attribute Name="'//strip(chAttribute)//'" Center="Node">'
    write(iUnit,'(A,I1,A,3I7,A,I15,A)')'        <DataItem Format="Binary" DataType="Float" Precision="',iprec, &
       '" Endian="Native" Dimensions="',i3,i2,i1,'" Seek="',disp,'">'
    write(iUnit,'(A)')'          '//strip(chFile)
    write(iUnit,'(A)')'        </DataItem>'
    write(iUnit,'(A)')'      </Attribute>'
  end subroutine Write_XDMF_One

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! dump_plane
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine dump_plane(ntime,ux,uy,uz,pr); implicit none
    integer,intent(in)::ntime
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pr
    
    ! locals
    integer,dimension(3)::sizes,subsizes,starts
    integer::ic,jc,kc, iplane, jplane, kplane, data_type,ierr,newtype,iUnit
        
    if((.not. DumpPlaneFlag) .or. mod(itime, RecordPlaneFreq) /= 0) return
    
    iTimePlane = iTimePlane + 1
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    if (save_ux) then
       if (nxPlane > 0) then
          do kc = y1start(3), y1end(3)
            do jc = y1start(2), y1end(2)
              do ic = 1, nxPlane
                iplane = ixPlane(ic)
                uxPlane_x(jc, kc, iTimePlane, ic) = real(ux(iplane, jc, kc), RKP_Dump)
              enddo
            enddo
          enddo
       endif
       if (nyPlane > 0) then
          do kc = y1start(3), y1end(3)
            do jc = 1, nyPlane
              jplane = iyPlane(jc)
              do ic = y1start(1), y1end(1)
                uxPlane_y(ic, kc, iTimePlane, jc) = real(ux(ic, jplane, kc), RKP_Dump)
              enddo
            enddo
          enddo
       endif
       if (nzPlane > 0) then
          do kc = 1, nzPlane
            kplane = izPlane(kc)
            do jc = y1start(2), y1end(2)
              do ic = y1start(1), y1end(1)
                uxPlane_z(ic, jc, iTimePlane, kc) = real(ux(ic, jc, kplane), RKP_Dump)
              enddo
            enddo
          enddo
       endif
    
    endif


    if (save_uy) then
       if (nxPlane > 0) then
          do kc = y1start(3), y1end(3)
            do jc = y1start(2), y1end(2)
              do ic = 1, nxPlane
                iplane = ixPlane(ic)
                uyPlane_x(jc, kc, iTimePlane, ic) = real(uy(iplane, jc, kc), RKP_Dump)
              enddo
            enddo
          enddo
       endif
       if (nyPlane > 0) then
          do kc = y1start(3), y1end(3)
            do jc = 1, nyPlane
              jplane = iyPlane(jc)
              do ic = y1start(1), y1end(1)
                uyPlane_y(ic, kc, iTimePlane, jc) = real(uy(ic, jplane, kc), RKP_Dump)
              enddo
            enddo
          enddo
       endif
       if (nzPlane > 0) then
          do kc = 1, nzPlane
            kplane = izPlane(kc)
            do jc = y1start(2), y1end(2)
              do ic = y1start(1), y1end(1)
                uyPlane_z(ic, jc, iTimePlane, kc) = real(uy(ic, jc, kplane), RKP_Dump)
              enddo
            enddo
          enddo
       endif
    endif


    if (save_uz) then
       if (nxPlane > 0) then
          do kc = y1start(3), y1end(3)
            do jc = y1start(2), y1end(2)
              do ic = 1, nxPlane
                iplane = ixPlane(ic)
                uzPlane_x(jc, kc, iTimePlane, ic) = real(uz(iplane, jc, kc), RKP_Dump)
              enddo
            enddo
          enddo
       endif
       if (nyPlane > 0) then
          do kc = y1start(3), y1end(3)
            do jc = 1, nyPlane
              jplane = iyPlane(jc)
              do ic = y1start(1), y1end(1)
                uzPlane_y(ic, kc, iTimePlane, jc) = real(uz(ic, jplane, kc), RKP_Dump)
              enddo
            enddo
          enddo
       endif
       if (nzPlane > 0) then
          do kc = 1, nzPlane
            kplane = izPlane(kc)
            do jc = y1start(2), y1end(2)
              do ic = y1start(1), y1end(1)
                uzPlane_z(ic, jc, iTimePlane, kc) = real(uz(ic, jc, kplane), RKP_Dump)
              enddo
            enddo
          enddo
       endif
    endif


    if (save_pr) then
       if (nxPlane > 0) then
          do kc = y1start(3), y1end(3)
            do jc = y1start(2), y1end(2)
              do ic = 1, nxPlane
                iplane = ixPlane(ic)
                prPlane_x(jc, kc, iTimePlane, ic) = real(pr(iplane, jc, kc), RKP_Dump)
              enddo
            enddo
          enddo
       endif
       if (nyPlane > 0) then
          do kc = y1start(3), y1end(3)
            do jc = 1, nyPlane
              jplane = iyPlane(jc)
              do ic = y1start(1), y1end(1)
                prPlane_y(ic, kc, iTimePlane, jc) = real(pr(ic, jplane, kc), RKP_Dump)
              enddo
            enddo
          enddo
       endif
       if (nzPlane > 0) then
          do kc = 1, nzPlane
            kplane = izPlane(kc)
            do jc = y1start(2), y1end(2)
              do ic = y1start(1), y1end(1)
                prPlane_z(ic, jc, iTimePlane, kc) = real(pr(ic, jc, kplane), RKP_Dump)
              enddo
            enddo
          enddo
       endif
    endif


    if(iTimePlane /= nTimePlane .and. itime /= ilast) return

    ! Write plane info ========================
    if(RKP_Dump == 4) then
      data_type = MPI_REAL
    else    
      data_type = MPI_DOUBLE_PRECISION    
    endif

    if (nxPlane > 0) then
       sizes(1)= nyc
       sizes(2)= nzc
       sizes(3)= iTimePlane
       subsizes(1)= y1size(2)
       subsizes(2)= y1size(3)
       subsizes(3)= iTimePlane
       starts(1)= y1start(2) - 1
       starts(2)= y1start(3) - 1
       starts(3)= 0
       call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierr)
       call MPI_TYPE_COMMIT(newtype,ierr)
       if (save_ux) call dump_one_plane('ux_xPlane', uxPlane_x, nxPlane, ixPlane, DECOMP_2D_COMM_ROW)
       if (save_uy) call dump_one_plane('uy_xPlane', uyPlane_x, nxPlane, ixPlane, DECOMP_2D_COMM_ROW)
       if (save_uz)  call dump_one_plane('uz_xPlane', uzPlane_x, nxPlane, ixPlane, DECOMP_2D_COMM_ROW)
       if (save_pr)  call dump_one_plane('pr_xPlane', prPlane_x, nxPlane, ixPlane, DECOMP_2D_COMM_ROW)
       call MPI_TYPE_FREE(newtype,ierr)
    endif

    if (nyPlane > 0) then
       sizes(1)= nxc
       sizes(2)= nzc
       sizes(3)= iTimePlane
       subsizes(1)= y1size(1)
       subsizes(2)= y1size(3)
       subsizes(3)= iTimePlane
       starts(1)= y1start(1) - 1
       starts(2)= y1start(3) - 1
       starts(3)= 0
       call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierr)
       call MPI_TYPE_COMMIT(newtype,ierr)
       if (save_ux) call dump_one_plane('ux_yPlane', uxPlane_y, nyPlane, iyPlane, MPI_COMM_WORLD)
       if (save_uy) call dump_one_plane('uy_yPlane', uyPlane_y, nyPlane, iyPlane, MPI_COMM_WORLD)
       if (save_uz) call dump_one_plane('uz_yPlane', uzPlane_y, nyPlane, iyPlane, MPI_COMM_WORLD)    
       if (save_pr) call dump_one_plane('pr_yPlane', prPlane_y, nyPlane, iyPlane, MPI_COMM_WORLD)
       call MPI_TYPE_FREE(newtype,ierr)
    endif

    if (nzPlane > 0) then
       sizes(1)= nxc
       sizes(2)= nyc
       sizes(3)= iTimePlane
       subsizes(1)= y1size(1)
       subsizes(2)= y1size(2)
       subsizes(3)= iTimePlane
       starts(1)= y1start(1) - 1
       starts(2)= y1start(2) - 1
       starts(3)= 0    
       call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierr)
       call MPI_TYPE_COMMIT(newtype,ierr)
       if (save_ux) call dump_one_plane('ux_zPlane', uxPlane_z, nzPlane, izPlane, DECOMP_2D_COMM_COL)
       if (save_uy) call dump_one_plane('uy_zPlane', uyPlane_z, nzPlane, izPlane, DECOMP_2D_COMM_COL)
       if (save_uz) call dump_one_plane('uz_zPlane', uzPlane_z, nzPlane, izPlane, DECOMP_2D_COMM_COL)
       if (save_pr) call dump_one_plane('pr_zPlane', prPlane_z, nzPlane, izPlane, DECOMP_2D_COMM_COL)
       call MPI_TYPE_FREE(newtype,ierr)
    endif

    iTimePlane=0

  contains

    subroutine dump_one_plane(file_str, plane_data, nPlane, iPlane, comm)
       implicit none
       character(len=*), intent(in) :: file_str
       real(RKP_Dump), intent(in) :: plane_data(:, :, :, :)
       integer, intent(in) :: nPlane
       integer, intent(in) :: iPlane(nPlane)
       integer, intent(in) :: comm

       integer :: ids, id_p
       character(len=:),allocatable::filename
       do ids = 1, nPlane
          id_p = iPlane(ids)
          filename = strip(PlaneDir_)// file_str // int2str(id_p, 5) // '_' // int2str(ntime,10)
          call MPI_FILE_OPEN(comm, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, iUnit, ierr)
          call MPI_BARRIER(comm, ierr)
          call my_mpi_file_set_size(iUnit, 0_MPI_OFFSET_KIND, ierr)  ! guarantee overwriting
          call MPI_BARRIER(comm, ierr)
          call MPI_FILE_SET_VIEW(iUnit, 0_MPI_OFFSET_KIND, data_type, newtype, 'native', MPI_INFO_NULL, ierr)
          call MPI_FILE_WRITE_ALL(iUnit, plane_data(:,:,1:iTimePlane,ids), subsizes(1)*subsizes(2)*subsizes(3), data_type,MPI_STATUS_IGNORE,ierr)
          call MPI_FILE_CLOSE(iUnit,ierr)
       enddo
    endsubroutine dump_one_plane

  endsubroutine dump_plane


end module f2_DumpPlane

#undef RKP_Dump
