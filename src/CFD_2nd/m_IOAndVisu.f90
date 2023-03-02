module m_IOAndVisu
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Tools,only: Clc_Q_vor,Clc_lamda2
  use m_Variables,only:mb1,OutFlowInfoX,OutFlowInfoY
  implicit none
  private

  ! VisuOption
  integer:: iskip,jskip,kskip
  integer:: Prev_BackUp_itime  = 53456791
  logical:: save_ux,save_uy,save_uz,save_wx,save_wy,save_wz,save_wMag
  logical:: save_pr,save_Q_vor,save_lamda2,WriteHistOld,ReadHistOld

  public:: InitVisu, dump_visu, read_restart, write_restart, Delete_Prev_Restart

contains

  !******************************************************************
  ! InitVisu
  !******************************************************************
  subroutine InitVisu(ChannelPrm)
    implicit none
    character(*),intent(in)::ChannelPrm

    ! locals
    character(128)::XdmfFile
    integer::nUnitFile,ierror,nflds,ifld,iprec,i,j,k
    NAMELIST /IO_Options/ save_ux,save_uy,save_uz,save_pr,save_wx,save_wy,save_wz,save_wMag,save_Q_vor, &
                          save_lamda2,WriteHistOld,ReadHistOld,iskip,jskip,kskip
 
    open(newunit=nUnitFile, file=ChannelPrm, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitVisu", "Cannot openfile: "//trim(ChannelPrm))
    read(nUnitFile, nml=IO_Options)
    close(nUnitFile,IOSTAT=ierror)
    if(nrank==0)write(MainLog%nUnit,nml=IO_Options)

    ! write XDMF file
    if(nrank/=0) return
    write(xdmfFile,"(A)") trim(ResultsDir)//"VisuFor"//trim(RunName)//".xmf"
    open(newunit=nUnitFile, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitVisu","Cannot open file: "//trim(XdmfFile))
    ! XDMF/XMF Title
    write(nUnitFile,'(A)') '<?xml version="1.0" ?>'
    write(nUnitFile,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(nUnitFile,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(nUnitFile,'(A)') '<Domain>'

    ! grid
    iprec=mytype_save
    write(nUnitFile,'(A,3I7,A)')'    <Topology name="TOPO" TopologyType="3DRectMesh" Dimensions="',nzc,nyc,nxc,'"/>'
    write(nUnitFile,'(A)')'    <Geometry name="GEO" GeometryType="VXVYVZ">'
    ! x-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nxc,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do i=1,nxc
      write(nUnitFile,'(E14.7)',advance='no') (i-1)*dx+dx*0.5_RK
    enddo
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'        </DataItem>'
    ! y-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nyc,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do j=1,nyc
      write(nUnitFile,'(E14.7)',advance='no') yc(j)
    enddo
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'        </DataItem>'
    ! z-grid
    write(nUnitFile,'(A,I1,A,I5,A)') '        <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nzc,'">'
    write(nUnitFile,'(A)',advance='no') '        '
    do k=1,nzc
      write(nUnitFile,'(E14.7)',advance='no') (k-1)*dz+dz*0.5_RK
    enddo
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'        </DataItem>'
    write(nUnitFile,'(A)')'    </Geometry>'

    ! Time series
    nflds = (ilast - ifirst +1)/SaveVisu  + 1
    write(nUnitFile,'(A)')'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    write(nUnitFile,'(A)')'        <Time TimeType="List">'
    write(nUnitFile,'(A,I6,A)')'        <DataItem Format="XML" NumberType="Int" Dimensions="',nflds,'">' 
    write(nUnitFile,'(A)',advance='no')'        '
    do ifld = ifirst-1,ilast,SaveVisu
      write(nUnitFile,'(I10)',advance='no') ifld
    enddo
    write(nUnitFile,'(A)')'        </DataItem>'
    write(nUnitFile,'(A)')' '; write(nUnitFile,'(A)')'       </Time>'

    ! attribute
    do  ifld=ifirst-1,ilast,SaveVisu
      write(nUnitFile,'(A,I10.10,A)')'        <Grid Name="T',ifld,'" GridType="Uniform">'
      write(nUnitFile,'(A)')'            <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(nUnitFile,'(A)')'            <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
      if(save_ux)    call Write_XDMF_One(nUnitFile,ifld,'ux')
      if(save_uy)    call Write_XDMF_One(nUnitFile,ifld,'uy')
      if(save_uz)    call Write_XDMF_One(nUnitFile,ifld,'uz')
      if(save_pr)    call Write_XDMF_One(nUnitFile,ifld,'pr')
      if(save_wx)    call Write_XDMF_One(nUnitFile,ifld,'wx')
      if(save_wy)    call Write_XDMF_One(nUnitFile,ifld,'wy')
      if(save_wz)    call Write_XDMF_One(nUnitFile,ifld,'wz')
      if(save_wMag)  call Write_XDMF_One(nUnitFile,ifld,'wMag')
      if(save_Q_vor) call Write_XDMF_One(nUnitFile,ifld,'Q' )
      if(save_lamda2)call Write_XDMF_One(nUnitFile,ifld,'lamda2')
      write(nUnitFile,'(A)')'        </Grid>'
    enddo

    write(nUnitFile,'(A)')'    </Grid>'
    write(nUnitFile,'(A)')'</Domain>'
    write(nUnitFile,'(A)')'</Xdmf>'
    close(nUnitFile,IOSTAT=ierror)
  end subroutine InitVisu

  !******************************************************************
  ! Write_XDMF_One
  !******************************************************************
  subroutine Write_XDMF_One(nUnitFile, ifld,chAttribute)
    implicit none
    integer,intent(in)::nUnitFile,ifld
    character(*),intent(in)::chAttribute

    ! locals
    character(128)::chFile
    integer::iprec=mytype_save
    
    write(chFile,'(A,A,I10.10)')"VisuFor"//trim(RunName),"_"//trim(adjustl(chAttribute))//"_",ifld
    write(nUnitFile,'(A)')'            <Attribute Name="'//trim(chAttribute)//'" Center="Node">'
    write(nUnitFile,'(A,I1,A,3I7,A)')'                <DataItem Format="Binary" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nzc,nyc,nxc,'">'
    write(nUnitFile,'(A)')'                    '//trim(chFile)
    write(nUnitFile,'(A)')'                </DataItem>'
    write(nUnitFile,'(A)')'            </Attribute>'
  end subroutine Write_XDMF_One

  !******************************************************************
  ! dump_visu
  !******************************************************************
  subroutine dump_visu(ntime,ux,uy,uz,pressure,ArrTemp)
    implicit none
    integer,intent(in)::ntime
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::ArrTemp

    ! locals
    character(128)::chFile
    integer::ic,jc,kc,ip,jp,kp,im,jm,km
    real(RK)::dudy,dudz,dvdx,dvdz,dwdx,dwdy
    real(RK)::caj,cac1,cac2,cac12,vor_x,vor_y,vor_z
 
    ! ux
    if(save_ux) then
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          do ic=y1start(1),y1end(1)
#ifdef SaveNode
            ArrTemp(ic,jc,kc)=ux(ic,jc,kc)
#else
            ArrTemp(ic,jc,kc)=0.5_RK*(ux(ic+1,jc,kc)+ux(ic,jc,kc))
#endif
          enddo
        enddo
      enddo
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_ux_",ntime
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
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_uy_",ntime
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
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_uz_",ntime
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! pressure
    if(save_pr) then
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_pr_",ntime
      call decomp_2d_write_every(y_pencil,pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),iskip,jskip,kskip,chFile,from1=.true.)
    endif

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
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_wx_",ntime
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
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_wy_",ntime
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
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_wz_",ntime
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
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_wMag_",ntime 
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif

    ! Q
    if(save_Q_vor) then
      call Clc_Q_vor(ux,uy,uz,ArrTemp)
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_Q_",ntime 
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif  
  
    ! lamda2
    if(save_lamda2) then
      call Clc_lamda2(ux,uy,uz,ArrTemp)
      write(chFile,"(A,I10.10)") trim(ResultsDir)//"VisuFor"//trim(RunName)//"_lamda2_",ntime 
      call decomp_2d_write_every(y_pencil,ArrTemp,iskip,jskip,kskip,chFile,from1=.true.)
    endif
  end subroutine dump_visu

  !**********************************************************************
  ! Delete_Prev_Restart
  !**********************************************************************
  subroutine Delete_Prev_Restart(ntime)
    implicit none
    integer,intent(in)::ntime

    ! locals
    integer::nUnit,ierror
    character(128)::chFile

    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    if(nrank/=0) return
    ! 
    write(chFile,"(A,I10.10)") trim(RestartDir)//"RestartFor"//trim(RunName),Prev_BackUp_itime
    open(newunit=nUnit,file=trim(chFile),IOSTAT=ierror)
    close(unit=nUnit,status='delete',IOSTAT=ierror)
    ! 
    write(chFile,"(A,I10.10)") trim(RestartDir)//"OutFlowFor"//trim(RunName),Prev_BackUp_itime 
    open(newunit=nUnit,file=trim(chFile),IOSTAT=ierror)
    close(unit=nUnit,status='delete',IOSTAT=ierror)
    !
    write(chFile,"(A,I10.10)") trim(RestartDir)//"PrDataFor"//trim(RunName),Prev_BackUp_itime 
    open(newunit=nUnit,file=trim(chFile),IOSTAT=ierror)
    close(unit=nUnit,status='delete',IOSTAT=ierror)
    !
    Prev_BackUp_itime = ntime
  end subroutine Delete_Prev_Restart

  !******************************************************************
  ! write_restart
  !******************************************************************
  subroutine write_restart(ntime,ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
    implicit none
    integer,intent(in)::ntime
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(in):: HistXOld,HistYOld,HistZOld

    ! locals
    character(128)::chFile
    integer(MPI_OFFSET_KIND)::disp
    integer::fh,ierror,ic,jc,kc,newtype,sizes(3),subsizes(3),starts(3)

    ! begin to write restart file
    write(chFile,"(A,I10.10)") trim(RestartDir)//"RestartFor"//trim(RunName),ntime
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    disp = 0_MPI_OFFSET_KIND
    call decomp_2d_write_var(fh,disp,y_pencil,      ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_write_var(fh,disp,y_pencil,      uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_write_var(fh,disp,y_pencil,      uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_write_var(fh,disp,y_pencil,pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    if(ischeme==FI_AB2 .and. WriteHistOld) then
      call decomp_2d_write_var(fh,disp,y_pencil,HistXOld)
      call decomp_2d_write_var(fh,disp,y_pencil,HistYOld)
      call decomp_2d_write_var(fh,disp,y_pencil,HistZOld)
    endif
    call MPI_FILE_CLOSE(fh,ierror)

    ! Write PrGradData
    if(IsUxConst .and. nrank==0) then
      write(chFile,"(A,I10.10)") trim(RestartDir)//"PrDataFor"//trim(RunName),ntime
      open(newunit=fh, file=chFile, status='replace', action='write', IOSTAT=ierror)
      if(ierror/=0) then
        call MainLog%CheckForError(ErrT_Abort,"write_restart","Cannot open file: "//trim(chFile))
      else
        write(fh,'(4ES26.17)') PrGradData(1:4)
      endif
      close(fh,IOSTAT=ierror)
    endif
    
    ! Write OutFlow
    if(BcOption(xp_dir)/=BC_OutFlow .and. BcOption(yp_dir)/=BC_OutFlow) return
    write(chFile,"(A,I10.10)") trim(RestartDir)//"OutFlowFor"//trim(RunName),ntime
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
    call MPI_FILE_SET_SIZE(fh,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    disp = 0_MPI_OFFSET_KIND 
    if(BcOption(yp_dir)==BC_OutFlow) then
      jc=nyp
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          OutFlowInfoY(4,ic,kc)=ux(ic,jc,kc)
          OutFlowInfoY(5,ic,kc)=uy(ic,jc,kc)
          OutFlowInfoY(6,ic,kc)=uz(ic,jc,kc)         
        enddo
      enddo
      sizes=[3, nxc, nzc]
      subsizes=[3,y1size(1),y1size(3)]
      starts=[0,y1start(1)-1,y1start(3)-1]
      call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,real_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(fh,OutFlowInfoY(4:6,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
      disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      if(ischeme==FI_AB2 .and. WriteHistOld) then
        call MPI_FILE_WRITE_ALL(fh,OutFlowInfoY(1:3,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      endif
      call MPI_TYPE_FREE(newtype,ierror)
      call MPI_FILE_CLOSE(fh,ierror)
    endif
    if(myProcNghBC(y_pencil,3)==BC_OutFlow) then
      ic=nxp
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          OutFlowInfoX(4,jc,kc)=ux(ic,jc,kc)
          OutFlowInfoX(5,jc,kc)=uy(ic,jc,kc)
          OutFlowInfoX(6,jc,kc)=uz(ic,jc,kc)          
        enddo
      enddo
      sizes=[3, nyc, nzc]
      subsizes=[3,y1size(2),y1size(3)]
      starts=[0,y1start(2)-1,y1start(3)-1]
      call MPI_FILE_OPEN(DECOMP_2D_COMM_ROW, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierror)
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,real_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_WRITE_ALL(fh,OutFlowInfoX(4:6,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
      disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      if(ischeme==FI_AB2 .and. WriteHistOld) then
        call MPI_FILE_WRITE_ALL(fh,OutFlowInfoX(1:3,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)    
      endif
      call MPI_TYPE_FREE(newtype,ierror)
      call MPI_FILE_CLOSE(fh,ierror)
    endif    
  end subroutine write_restart

  !******************************************************************
  ! read_restart
  !******************************************************************
  subroutine read_restart(ux,uy,uz,pressure,HistXOld,HistYOld,HistZOld)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz,pressure
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out):: HistXOld,HistYOld,HistZOld

    ! locals
    character(128)::chFile
    integer(MPI_OFFSET_KIND)::disp,byte_total1,byte_total2,filebyte
    integer::fh,ierror,ntime,ic,jc,kc,newtype,sizes(3),subsizes(3),starts(3)
    
    ! begin to write restart file
    ntime= ifirst - 1
    write(chFile,"(A,I10.10)") trim(RestartDir)//"RestartFor"//trim(RunName),ntime
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror)
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"Read_Restart","Cannot open file: "//trim(chFile))
    call MPI_FILE_GET_SIZE(fh,filebyte,ierror)
    byte_total1=int(mytype_bytes,8)*int(nxc,8)*int(nyc,8)*int(nzc,8)*7_MPI_OFFSET_KIND
    byte_total2=int(mytype_bytes,8)*int(nxc,8)*int(nyc,8)*int(nzc,8)*4_MPI_OFFSET_KIND
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
    call decomp_2d_read_var(fh,disp,y_pencil,      ux(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_read_var(fh,disp,y_pencil,      uy(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_read_var(fh,disp,y_pencil,      uz(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    call decomp_2d_read_var(fh,disp,y_pencil,pressure(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)))
    if(ischeme==FI_AB2 .and. ReadHistOld) then
      call decomp_2d_read_var(fh,disp,y_pencil,HistXOld)
      call decomp_2d_read_var(fh,disp,y_pencil,HistYOld)
      call decomp_2d_read_var(fh,disp,y_pencil,HistZOld)
    endif
    call MPI_FILE_CLOSE(fh,ierror)

    ! Read PrGradData
    if(IsUxConst) then
      write(chFile,"(A,I10.10)") trim(RestartDir)//"PrDataFor"//trim(RunName),ntime
      open(newunit=fh, file=chFile, status='old', action='read', IOSTAT=ierror)
      if(ierror/=0) then
        if(nrank==0) then
          call MainLog%OutInfo("read_restart: Cannot open file "//trim(chFile),1)
          call MainLog%OutInfo(" PrGradData=0.0 will be used ! ",2)
        endif
        PrGradData=0.0_RK
      else
        read(fh,*) PrGradData(1:4)
      endif
      close(fh,IOSTAT=ierror)
    endif
    
    ! Read OutFlow
    if(BcOption(xp_dir)/=BC_OutFlow .and. BcOption(yp_dir)/=BC_OutFlow) return
    write(chFile,"(A,I10.10)") trim(RestartDir)//"OutFlowFor"//trim(RunName),ntime
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror)
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"Read_Restart","Cannot open file: "//trim(chFile))
    call MPI_FILE_GET_SIZE(fh,filebyte,ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    byte_total1=0; byte_total2=0;
    if(BcOption(xp_dir)==BC_OutFlow) then
      byte_total1=byte_total1+int(mytype_bytes,8)*int(nyc,8)*int(nzc,8)*6_8
      byte_total2=byte_total2+int(mytype_bytes,8)*int(nyc,8)*int(nzc,8)*3_8
    endif
    if(BcOption(yp_dir)==BC_OutFlow) then
      byte_total1=byte_total1+int(mytype_bytes,8)*int(nxc,8)*int(nzc,8)*6_8
      byte_total2=byte_total2+int(mytype_bytes,8)*int(nxc,8)*int(nzc,8)*3_8
    endif
    if(ischeme==FI_AB2 .and. ReadHistOld) then
      if(filebyte /= byte_total1 .and. nrank==0) then
        call MainLog%CheckForError(ErrT_Abort,"Read_Restart","file byte wrong3")
      endif      
    else
      if((filebyte /= byte_total1 .and. filebyte /= byte_total2) .and. nrank==0) then
        call MainLog%CheckForError(ErrT_Abort,"Read_Restart","file byte wrong4")
      endif      
    endif
    disp = 0_MPI_OFFSET_KIND 
    if(BcOption(yp_dir)==BC_OutFlow) then
      sizes=[3, nxc, nzc]
      subsizes=[3,y1size(1),y1size(3)]
      starts=[0,y1start(1)-1,y1start(3)-1]
      call MPI_FILE_OPEN(MPI_COMM_WORLD, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror)
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,real_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_READ_ALL(fh,OutFlowInfoY(4:6,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
      disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      if(ischeme==FI_AB2 .and. WriteHistOld) then
        call MPI_FILE_READ_ALL(fh,OutFlowInfoY(1:3,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      endif
      call MPI_TYPE_FREE(newtype,ierror)
      call MPI_FILE_CLOSE(fh,ierror)
      jc=nyp
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,jc,kc)=OutFlowInfoY(4,ic,kc)
          uy(ic,jc,kc)=OutFlowInfoY(5,ic,kc)
          uz(ic,jc,kc)=OutFlowInfoY(6,ic,kc)       
        enddo
      enddo   
    endif    
    if(myProcNghBC(y_pencil,3)==BC_OutFlow) then
      sizes=[3, nyc, nzc]
      subsizes=[3,y1size(2),y1size(3)]
      starts=[0,y1start(2)-1,y1start(3)-1]
      call MPI_FILE_OPEN(DECOMP_2D_COMM_ROW, chFile, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierror)
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,real_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
      call MPI_FILE_READ_ALL(fh,OutFlowInfoX(4:6,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
      disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      if(ischeme==FI_AB2 .and. WriteHistOld) then
        call MPI_FILE_READ_ALL(fh,OutFlowInfoX(1:3,:,:),subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
      endif
      call MPI_TYPE_FREE(newtype,ierror)
      call MPI_FILE_CLOSE(fh,ierror)
      ic=nxp
      do kc=y1start(3),y1end(3)
        do jc=y1start(2),y1end(2)
          ux(ic,jc,kc)=OutFlowInfoX(4,jc,kc)
          uy(ic,jc,kc)=OutFlowInfoX(5,jc,kc)
          uz(ic,jc,kc)=OutFlowInfoX(6,jc,kc)         
        enddo
      enddo                      
    endif
  end subroutine read_restart
end module m_IOAndVisu
