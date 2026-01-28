module m_IOAndVisu
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only:mb1 
  use m_Tools,only: Clc_Q_vor,Clc_lamda2
  implicit none
  private

  ! VisuOption
  integer:: iskip,jskip,kskip
  integer:: Prev_BackUp_itime  = 53456791
  logical:: save_ux,save_uy,save_uz,save_wx,save_wy,save_wz,save_wMag
  logical:: save_pr,save_Q_vor,save_lamda2,WriteHistOld,ReadHistOld
#ifdef ScalarFlow
  logical:: save_scalar
#endif
  public:: InitVisu, dump_visu, read_restart, write_restart, Delete_Prev_Restart

contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitVisu
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitVisu(ChannelPrm); implicit none
    character(len=*),intent(in)::ChannelPrm

    ! locals
    character(len=:),allocatable::XdmfFile
    integer::iUnit,ierr,nfld,ifld,iprec,i,j,k, int_t, nxc_o, nyc_o, nzc_o, line_num
    NAMELIST /IO_Options/ save_ux,save_uy,save_uz,save_pr,save_wx,save_wy,save_wz,save_wMag,save_Q_vor, &
                          save_lamda2,WriteHistOld,ReadHistOld,iskip,jskip,kskip
#ifdef ScalarFlow
    NAMELIST /SaveScalarOption/save_scalar
#endif
 
    open(newunit=iUnit, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierr )
    if(ierr/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitVisu", "Cannot open file: "//strip(ChannelPrm))
    read(iUnit, nml=IO_Options)
#ifdef ScalarFlow
    rewind(iUnit)
    read(iUnit, nml=SaveScalarOption)
#endif
    close(iUnit,IOSTAT=ierr)
    if(nrank==0) then
      write(MainLog%iUnit, nml=IO_Options)
#ifdef ScalarFlow
      write(MainLog%iUnit, nml=SaveScalarOption)
#endif
    endif

    ! Write XDMF file
    if(nrank/=0) return
    xdmfFile = strip(ResultsDir_)//"VisuFor"//strip(RunName_)//".xmf"
    open(newunit=iUnit, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierr)
    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"InitVisu","Cannot open file: "//trim(XdmfFile))
    write(iUnit,'(A)') '<?xml version="1.0" ?>'
    write(iUnit,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(iUnit,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(iUnit,'(A)') '<Domain>'

    ! grid
    nxc_o = nxc / iskip
    int_t = mod(nxc, iskip)
    if (int_t == 1) nxc_o = nxc_o + 1

    nyc_o = nyc / jskip
    int_t = mod(nyc, jskip)
    if (int_t == 1) nyc_o = nyc_o + 1

    nzc_o = nzc / kskip
    int_t = mod(nzc, kskip)
    if (int_t == 1) nzc_o = nzc_o + 1

    line_num = 8
    iprec=mytype_save
    
    write(iUnit,'(A,3I7,A)')'  <Topology name="TOPO" TopologyType="3DRectMesh" Dimensions="',nzc_o, nyc_o, nxc_o, '"/>'
    write(iUnit,'(A)')'  <Geometry name="GEO" GeometryType="VXVYVZ">'
    ! x-grid
    write(iUnit,'(A,I1,A,I5,A)') '    <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nxc_o,'">'
    write(iUnit,'(A)',advance='no') '    '
    int_t = 0
    do i = 1, nxc, iskip
      int_t = int_t + 1
      if(mod(int_t, line_num) == 0) then
        write(iUnit,'(ES15.7)') (i - 1)*dx+dx*0.5_RK
        if(int_t < nxc_o) write(iUnit,'(A)',advance='no') '    '
      else
        write(iUnit,'(ES15.7)',advance='no') (i - 1)*dx+dx*0.5_RK
      endif
    enddo
    if(mod(nxc_o, line_num) /= 0) write(iUnit, *) ' '    
    write(iUnit,'(A)')'    </DataItem>'

    ! y-grid
    write(iUnit,'(A,I1,A,I5,A)') '    <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nyc_o,'">'
    write(iUnit,'(A)',advance='no') '    '
    int_t = 0
    do j = 1, nyc, jskip
      int_t = int_t + 1
      if(mod(int_t, line_num) == 0) then
        write(iUnit,'(ES15.7)') yc(j)
        if(int_t < nyc_o) write(iUnit,'(A)',advance='no') '    '
      else
        write(iUnit,'(ES15.7)',advance='no') yc(j)
      endif
    enddo
    if(mod(nyc_o, line_num) /=0) write(iUnit,*)' '
    write(iUnit,'(A)')'    </DataItem>'

    ! z-grid
    write(iUnit,'(A,I1,A,I5,A)') '    <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nzc_o,'">'
    write(iUnit,'(A)',advance='no') '    '
    int_t = 0
    do k = 1, nzc, kskip
      int_t = int_t + 1
      if(mod(int_t, line_num) == 0) then
        write(iUnit,'(ES15.7)') (k-1)*dz+dz*0.5_RK
        if(int_t < nzc_o) write(iUnit,'(A)',advance='no') '    '
      else
        write(iUnit,'(ES15.7)',advance='no') (k-1)*dz+dz*0.5_RK
      endif
    enddo
    if(mod(nzc_o, line_num) /= 0) write(iUnit,*) ' '
    write(iUnit,'(A)')'    </DataItem>'
    write(iUnit,'(A)')'  </Geometry>'

    ! Time series
    nfld = (ilast - ifirst +1)/SaveVisu  + 1
    write(iUnit,'(A)')'  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    write(iUnit,'(A)')'    <Time TimeType="List">'
    write(iUnit,'(A,I6,A)')'      <DataItem Format="XML" NumberType="Int" Dimensions="',nfld,'">' 
    write(iUnit,'(A)',advance='no')'        '
    do ifld = 1, nfld
      if(mod(ifld, line_num)==0) then
        write(iUnit,'(I10)') (ifld-1)*SaveVisu +(ifirst-1)
        if(ifld < nfld) write(iUnit,'(A)',advance='no') '        '
      else
        write(iUnit,'(I10)',advance='no') (ifld-1)*SaveVisu +(ifirst-1)
      endif     
    enddo
    if(mod(nfld, line_num) /=0) write(iUnit,*)' '
    write(iUnit,'(A)') '      </DataItem>'
    write(iUnit,'(A)') '    </Time>'

    ! attribute
    do ifld=ifirst-1,ilast,SaveVisu
      write(iUnit,'(A,I10.10,A)')'    <Grid Name="T',ifld,'" GridType="Uniform">'
      write(iUnit,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(iUnit,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
      if(save_ux)     call Write_XDMF_One(iUnit,ifld,'ux', nxc_o, nyc_o, nzc_o)
      if(save_uy)     call Write_XDMF_One(iUnit,ifld,'uy', nxc_o, nyc_o, nzc_o)
      if(save_uz)     call Write_XDMF_One(iUnit,ifld,'uz', nxc_o, nyc_o, nzc_o)
      if(save_pr)     call Write_XDMF_One(iUnit,ifld,'pr', nxc_o, nyc_o, nzc_o)
      if(save_wx)     call Write_XDMF_One(iUnit,ifld,'wx', nxc_o, nyc_o, nzc_o)
      if(save_wy)     call Write_XDMF_One(iUnit,ifld,'wy', nxc_o, nyc_o, nzc_o)
      if(save_wz)     call Write_XDMF_One(iUnit,ifld,'wz', nxc_o, nyc_o, nzc_o)
      if(save_wMag)   call Write_XDMF_One(iUnit,ifld,'wMag', nxc_o, nyc_o, nzc_o)
      if(save_Q_vor)  call Write_XDMF_One(iUnit,ifld,'Q' , nxc_o, nyc_o, nzc_o)
      if(save_lamda2) call Write_XDMF_One(iUnit,ifld,'lambda2', nxc_o, nyc_o, nzc_o)
#ifdef ScalarFlow
     if(save_scalar) call Write_XDMF_One(iUnit,ifld,'scalar', nxc_o, nyc_o, nzc_o)
#endif
      write(iUnit,'(A)')'    </Grid>'
    enddo

    write(iUnit,'(A)')'  </Grid>'
    write(iUnit,'(A)')'</Domain>'
    write(iUnit,'(A)')'</Xdmf>'
    close(iUnit,IOSTAT=ierr)
#ifdef SaveNode
    call MainLog%OutInfo("Choose to save the visualizing file at grid node",2)
#else
    call MainLog%OutInfo("Choose to save the visualizing file at cell center",2)
#endif
  end subroutine InitVisu

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Write_XDMF_One
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Write_XDMF_One(iUnit, ifld,chAttribute, nxc_o, nyc_o, nzc_o); implicit none
    integer,intent(in)::iUnit, ifld, nxc_o, nyc_o, nzc_o
    character(len=*),intent(in)::chAttribute

    ! locals
    integer::iprec
    character(len=:),allocatable::chFile
    
    iprec = mytype_save
    chFile = "VisuFor" // strip(RunName_) // "_" // strip(chAttribute) // "_" // int2str(ifld,10)
    write(iUnit,'(A)')'      <Attribute Name="'//trim(chAttribute)//'" Center="Node">'
    write(iUnit,'(A,I1,A,3I7,A)')'        <DataItem Format="Binary" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="', &
                                 nzc_o, nyc_o, nxc_o, '">'
    write(iUnit,'(A)')'          '//trim(chFile)
    write(iUnit,'(A)')'        </DataItem>'
    write(iUnit,'(A)')'      </Attribute>'
  end subroutine Write_XDMF_One

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! dump_visu
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
#ifdef ScalarFlow
  subroutine dump_visu(ntime,ux,uy,uz,pressure,scalar,ArrTemp); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure,scalar
#else
  subroutine dump_visu(ntime,ux,uy,uz,pressure,ArrTemp); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
#endif
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::ArrTemp
    integer,intent(in)::ntime

    ! locals
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    character(len=:),allocatable::chFile
    real(RK)::dudy,dudz,dvdx,dvdz,dwdx,dwdy
    real(RK)::caj,cac1,cac2,cac12,vor_x,vor_y,vor_z
 
    chFile = ' '
    ! ux
    if(save_ux) then
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
#ifdef SaveNode
             ArrTemp(ic,jc,kc)=ux(ic,jc,kc)+uCRF
#else
             ArrTemp(ic,jc,kc)=0.5_RK*(ux(ic+1,jc,kc)+ux(ic,jc,kc))+uCRF
#endif
          enddo
        enddo
      enddo
      chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_ux_" //int2str(ntime,10)
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! uy
    if(save_uy) then
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
#ifdef SaveNode
             ArrTemp(ic,jc,kc)=uy(ic,jc,kc)
#else
             ArrTemp(ic,jc,kc)=0.5_RK*(uy(ic,jc+1,kc)+uy(ic,jc,kc))
#endif
          enddo
        enddo
      enddo
      chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_uy_" //int2str(ntime,10)
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! uz
    if(save_uz) then
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
#ifdef SaveNode
             ArrTemp(ic,jc,kc)=uz(ic,jc,kc)
#else
             ArrTemp(ic,jc,kc)=0.5_RK*(uz(ic,jc,kc+1)+uz(ic,jc,kc))
#endif
          enddo
        enddo
      enddo
      chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_uz_" //int2str(ntime,10)
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! pressure
    if(save_pr) then
      chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_pr_" //int2str(ntime,10)
      call decomp_2d_write_every(y_pencil,pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),iskip,jskip,kskip,chFile,from1=.true.)
    endif

#ifdef ScalarFlow
    if(save_scalar) then
      chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_scalar_" //int2str(ntime,10)
      call decomp_2d_write_every(y_pencil,scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),iskip,jskip,kskip,chFile,from1=.true.)
    endif
#endif

    ! wx
    if(save_wx) then
      do kc=y1start(3),y1end(3)
        kp=kc+1
        km=kc-1
        do jc=y1start(2),y1end(2)
          jp=jc+1
          jm=jc-1
          cac1 = rdyc(jc)
          cac2 = rdyc(jp)    
          cac12= cac1 - cac2
          do ic=y1start(1),y1end(1)
            dvdz=  (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz *0.25_RK
            dwdy= ((uz(ic,jp,kp) +uz(ic,jp,kc))*cac2    &
                  +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12   &
                  -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1    )    *0.25_RK
            ArrTemp(ic,jc,kc)= dwdy -dvdz
          enddo
        enddo
      enddo
      chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_wx_" //int2str(ntime,10)
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wy
    if(save_wy) then
      do kc=y1start(3),y1end(3)
        kp=kc+1
        km=kc-1
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
            ip=ic+1
            im=ic-1
            dudz=  (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz *0.25_RK
            dwdx=  (uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx *0.25_RK
            ArrTemp(ic,jc,kc)= dudz -dwdx
          enddo
        enddo
      enddo
      chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_wy_" //int2str(ntime,10)
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wz
    if(save_wz) then
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          jp=jc+1
          jm=jc-1
          caj  = rdyp(jc)
          cac1 = rdyc(jc)
          cac2 = rdyc(jp)    
          cac12= cac1 - cac2
          do ic=y1start(1),y1end(1)
            ip=ic+1
            im=ic-1
            dudy= ((ux(ip,jp,kc) +ux(ic,jp,kc))*cac2  &
                  +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12 &
                  -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1  )  *0.25_RK
            dvdx=  (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx *0.25_RK
            ArrTemp(ic,jc,kc)= dvdx -dudy
          enddo
        enddo
      enddo
      chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_wz_" //int2str(ntime,10)
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! wMag
    if(save_wMag) then
      do kc=y1start(3),y1end(3)
        kp=kc+1
        km=kc-1
        do jc=y1start(2),y1end(2)
          jp=jc+1
          jm=jc-1
          caj  = rdyp(jc)
          cac1 = rdyc(jc)
          cac2 = rdyc(jp)    
          cac12= cac1 - cac2
          do ic=y1start(1),y1end(1)
            ip=ic+1
            im=ic-1

            dudy= ((ux(ip,jp,kc) +ux(ic,jp,kc))*cac2    &
                  +(ux(ip,jc,kc) +ux(ic,jc,kc))*cac12   &
                  -(ux(ip,jm,kc) +ux(ic,jm,kc))*cac1    )    *0.25_RK
            dudz=  (ux(ip,jc,kp) +ux(ic,jc,kp) -ux(ip,jc,km) -ux(ic,jc,km))*rdz *0.25_RK
      
            dvdx=  (uy(ip,jp,kc) -uy(im,jp,kc) +uy(ip,jc,kc) -uy(im,jc,kc))*rdx *0.25_RK
            dvdz=  (uy(ic,jp,kp) +uy(ic,jc,kp) -uy(ic,jp,km) -uy(ic,jc,km))*rdz *0.25_RK

            dwdx=  (uz(ip,jc,kp) -uz(im,jc,kp) +uz(ip,jc,kc) -uz(im,jc,kc))*rdx *0.25_RK
            dwdy= ((uz(ic,jp,kp) +uz(ic,jp,kc))*cac2    &
                  +(uz(ic,jc,kp) +uz(ic,jc,kc))*cac12   &
                  -(uz(ic,jm,kp) +uz(ic,jm,kc))*cac1    )    *0.25_RK

            vor_x= dwdy-dvdz
            vor_y= dudz-dwdx
            vor_z= dvdx-dudy
            ArrTemp(ic,jc,kc)= sqrt(vor_x*vor_x +vor_y*vor_y +vor_z*vor_z)
          enddo
        enddo
      enddo
      chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_wmag_" //int2str(ntime,10)
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! Q
    if(save_Q_vor) then
      call Clc_Q_vor(ux,uy,uz,ArrTemp)
      chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_Q_" //int2str(ntime,10) 
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif  
  
   ! lamda2
   if(save_lamda2) then
     call Clc_lamda2(ux,uy,uz,ArrTemp)
     chFile = strip(ResultsDir_) // "VisuFor" // strip(RunName_) // "_lambda2_" //int2str(ntime,10)
     call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
   endif

  end subroutine dump_visu

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Delete_Prev_Restart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Delete_Prev_Restart(ntime); implicit none
    integer,intent(in):: ntime

    ! locals
    integer::iUnit,ierr
    character(len=:),allocatable::chFile

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if(nrank/=0) return
    !
    chFile = strip(RestartDir_) // "RestartFor" // strip(RunName_) // int2str(Prev_BackUp_itime,10)
    open(newunit=iUnit,file=chFile,IOSTAT=ierr)
    close(unit=iUnit,status='delete',IOSTAT=ierr)
    !
    chFile = strip(RestartDir_) // "PrDataFor" // strip(RunName_) // int2str(Prev_BackUp_itime,10)
    open(newunit=iUnit,file=chFile,IOSTAT=ierr)
    close(unit=iUnit,status='delete',IOSTAT=ierr)
    ! 
    Prev_BackUp_itime = ntime
  end subroutine Delete_Prev_Restart

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! write_restart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
#ifdef ScalarFlow
  subroutine write_restart(ntime,ux,uy,uz,pressure,scalar,HistXOld,HistYOld,HistZOld,HistCOld); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure,scalar
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in):: HistXOld,HistYOld,HistZOld,HistCOld
#else
  subroutine write_restart(ntime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in):: HistXOld,HistYOld,HistZOld
#endif
    integer,intent(in)::ntime
    
    ! locals
    integer::fh,ierr
    character(len=:),allocatable::chFile
    integer(kind=MPI_OFFSET_KIND)::disp

    ! begin to write restart file
    chFile = strip(RestartDir_) // "RestartFor" // strip(RunName_) // int2str(ntime,10)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_File_set_size(fh,0_MPI_OFFSET_KIND,ierr)  ! guarantee overwriting
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    disp = 0_MPI_OFFSET_KIND

    call decomp_2d_write_var(fh,disp,y_pencil,      ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))+uCRF)
    call decomp_2d_write_var(fh,disp,y_pencil,      uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_write_var(fh,disp,y_pencil,      uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_write_var(fh,disp,y_pencil,pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
#ifdef ScalarFlow
    call decomp_2d_write_var(fh,disp,y_pencil,scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
#endif
    if(ischeme==FI_AB2 .and. WriteHistOld) then
      call decomp_2d_write_var(fh,disp,y_pencil,HistXOld)
      call decomp_2d_write_var(fh,disp,y_pencil,HistYOld)
      call decomp_2d_write_var(fh,disp,y_pencil,HistZOld)
#ifdef ScalarFlow
      call decomp_2d_write_var(fh,disp,y_pencil,HistCOld)
#endif
    endif
    call MPI_FILE_CLOSE(fh,ierr)

    ! Write PrGradData
    if(IsUxConst .and. nrank==0) then
      chFile = strip(RestartDir_) // "PrDataFor" // strip(RunName_) // int2str(ntime,10)
      open(newunit=fh, file=chFile, status='replace', action='write', IOSTAT=ierr)
      if(ierr/=0) then
        call MainLog%CheckForError(ErrT_Abort,"write_restart","Cannot open file: "//chFile)
      else
        write(fh,'(4ES26.17)') PrGradData(1:4)
      endif
      close(fh,IOSTAT=ierr)
    endif
  end subroutine write_restart

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! read_restart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
#ifdef ScalarFlow
  subroutine read_restart(ux,uy,uz,pressure,scalar,HistXOld,HistYOld,HistZOld,HistCOld); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz,pressure,scalar
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out):: HistXOld,HistYOld,HistZOld,HistCOld
#else
  subroutine read_restart(ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz,pressure
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out):: HistXOld,HistYOld,HistZOld
#endif

    ! locals
    character(len=:),allocatable::chFile
    integer::fh,ierr,ntime
    integer(kind=MPI_OFFSET_KIND)::disp,byte_total1,byte_total2,filebyte

    ! begin to write restart file
    ntime= ifirst - 1
    chFile = strip(RestartDir_) // "RestartFor" // strip(RunName_) // int2str(ntime,10)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    if(ierr/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"Read_Restart","Cannot open file: "//chFile)

    call MPI_FILE_GET_SIZE(fh,filebyte,ierr)
#ifdef ScalarFlow
    byte_total1=int(mytype_bytes,8)*int(nxc,8)*int(nyc,8)*int(nzc,8)*9_MPI_OFFSET_KIND
    byte_total2=int(mytype_bytes,8)*int(nxc,8)*int(nyc,8)*int(nzc,8)*5_MPI_OFFSET_KIND
#else
    byte_total1=int(mytype_bytes,8)*int(nxc,8)*int(nyc,8)*int(nzc,8)*7_MPI_OFFSET_KIND
    byte_total2=int(mytype_bytes,8)*int(nxc,8)*int(nyc,8)*int(nzc,8)*4_MPI_OFFSET_KIND
#endif
    if(ischeme==FI_AB2 .and. ReadHistOld) then
      if(filebyte /= byte_total1 .and. nrank==0) then
        call MainLog%CheckForError(ErrT_Abort,"Read_Restart","file byte wrong1")
      endif      
    else
      if((filebyte /= byte_total1 .and. filebyte /= byte_total2) .and. nrank==0) then
        call MainLog%CheckForError(ErrT_Abort,"Read_Restart","file byte wrong2")
      endif      
    endif

    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_read_var(fh,disp,y_pencil,      ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)));ux=ux-uCRF
    call decomp_2d_read_var(fh,disp,y_pencil,      uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_read_var(fh,disp,y_pencil,      uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_read_var(fh,disp,y_pencil,pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
#ifdef ScalarFlow
    call decomp_2d_read_var(fh,disp,y_pencil,scalar(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
#endif
    if(ischeme==FI_AB2 .and. ReadHistOld) then
      call decomp_2d_read_var(fh,disp,y_pencil,HistXOld)
      call decomp_2d_read_var(fh,disp,y_pencil,HistYOld)
      call decomp_2d_read_var(fh,disp,y_pencil,HistZOld)
#ifdef ScalarFlow
      call decomp_2d_read_var(fh,disp,y_pencil,HistCOld)
#endif
    endif
    call MPI_FILE_CLOSE(fh,ierr)

    ! Read PrGradData
    if(IsUxConst) then
      chFile = strip(RestartDir_) // "PrDataFor" // strip(RunName_) // int2str(ntime,10)
      open(newunit=fh, file=chFile, status='old', action='read', IOSTAT=ierr)
      if(ierr/=0) then
        if(nrank==0) then
          call MainLog%OutInfo("read_restart: Cannot open file "// chFile,1)
          call MainLog%OutInfo(" PrGradData=0.0 will be used ! ",2)
        endif
        PrGradData=0.0_RK
      else
        read(fh,*) PrGradData(1:4)
      endif
      close(fh,IOSTAT=ierr)
    endif
  end subroutine read_restart

end module m_IOAndVisu
