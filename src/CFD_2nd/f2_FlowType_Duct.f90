#include "definitions_inc.f90"
module f2_FlowType_Duct
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d
  use f2_Parameters
  use f2_MeshAndMetries
  use f2_Variables,only:mb1
  use f2_Tools,only:CalcUxAver
  implicit none
  private

  ! statistics variabls
  integer:: nfstime
  real(RK),allocatable,dimension(:,:,:):: SumStat_plane_x
  
  public:: InitVelocity_Duct, InitStatVar_Duct, clcStat_Duct

#define _nDuct_Stat_ 11

contains


  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitVelocity_Duct
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitVelocity_Duct(ux,uy,uz,Deviation, input_file); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
    character(len = *), intent(in) :: input_file

    ! locals
    real(RK) :: ratiot, retau_guass,utau_guass,height,rem,twopi
    integer :: iUnit, ierr, vertical_direction,iTV(8), ic,jc,kc,m1,m2
    namelist /Duct_Param/ vertical_direction
    real(RK)::xplus,yplus,zplus, ybar,zbar, xct,yct,zct, wx,wy,wz, xlxPlus,ylyPlus,zlzPlus

    
    vertical_direction = 0
    open(newunit=iUnit, file=input_file, status='old',form='formatted',IOSTAT=ierr )
    if(ierr/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitVelocity_Duct","Cannot open file: "//strip(input_file))
    endif
    read(iUnit,nml=Duct_Param)
    close(iUnit,IOSTAT=ierr)
 
    call date_and_time(values=iTV); !iTV=0
    call random_seed(size= ic)
    call random_seed(put = iTV(7)*iTV(8)+[(jc,jc=1,ic)])
    call random_number(Deviation)
    Deviation= 0.2_RK* Deviation + 0.9_RK ! [0.8, 1.2] 

    height = 0.0_RK
    IF (vertical_direction == 2) THEN ! y is the vertical direction  
       height = yly
       if(nrank==0) print*, 'Here y is the vertical direction'
    ELSEIF (vertical_direction == 3) THEN ! z is the vertical direction
       height = zlz
       if(nrank==0) print*, 'Here z is the vertical direction'
    ELSE
       if (nrank == 0) call MainLog%CheckForError(ErrT_Abort,"InitVelocity_Duct","vertical direction of duct is wrong") 
    ENDIF

    if(FlowType==FT_CH) height=0.5_RK*height
    rem=uBulk*height/xnu
    ux=0.0_RK; uy=0.0_RK; uz=0.0_RK
    if(abs(uBulk)<1.0E-12)return
    
    retau_guass = 0.1538_RK*rem**0.887741_RK
    utau_guass  = retau_guass*xnu/height
    if(nrank==0) print*,' retau_gauss=',retau_guass
    if(nrank==0) print*,' utau_gauss= ',utau_guass
    twopi=2.0_RK*PI

    IF (vertical_direction == 2) THEN ! y is the vertical direction

    ! modulation of the random noise + initial velocity profile
    wx=twopi/500.0_RK; wz=twopi/200.0_RK
    xlxPlus=xlx*utau_guass/xnu;   zlzPlus=zlz*utau_guass/xnu;
    m1=floor(xlxPlus*wx/twopi)+1; wx=real(m1,RK)*twopi/xlxPlus
    m2=floor(zlzPlus*wz/twopi)+1; wz=real(m2,RK)*twopi/zlzPlus
    do jc=y1start(2),y1end(2)
      yct = height-abs(height-yc(jc))
      ybar= yct/height; yplus=utau_guass*yct/xnu
      do kc=y1start(3),y1end(3)
        zct  =real(kc-1,kind=RK)*dz+dz*0.5_RK
        zplus=utau_guass*zct/xnu
        do ic=y1start(1),y1end(1)
          xct  =real(ic-1,kind=RK)*dx+dx*0.5_RK
          xplus=utau_guass*xct/xnu
          !ux(ic,jc,kc) = 0.0052_RK*uBulk*yplus*exp(-yplus*yplus/1800.0_RK)*cos(wz*zplus)*Deviation(ic,jc,kc) ! original expression
          !uz(ic,jc,kc) = 0.0050_RK*uBulk*yplus*exp(-yplus*yplus/1800.0_RK)*sin(wx*xplus)*Deviation(ic,jc,kc) ! original expression
          ux(ic,jc,kc) = uBulk*ybar*exp(-4.5_RK*ybar*ybar)*cos(wz*zplus)*Deviation(ic,jc,kc)
          uz(ic,jc,kc) = uBulk*ybar*exp(-4.5_RK*ybar*ybar)*sin(wx*xplus)*Deviation(ic,jc,kc)
          ux(ic,jc,kc) = ux(ic,jc,kc)+ 3.0_RK*uBulk*(ybar-0.5_RK*ybar*ybar)
        enddo
      enddo
    enddo

    ELSEIF (vertical_direction == 3) THEN ! z is the vertical direction
    
    ! modulation of the random noise + initial velocity profile
    wx=twopi/500.0_RK; wy=twopi/200.0_RK
    xlxPlus=xlx*utau_guass/xnu;   ylyPlus=yly*utau_guass/xnu;
    m1=floor(xlxPlus*wx/twopi)+1; wx=real(m1,RK)*twopi/xlxPlus
    m2=floor(ylyPlus*wy/twopi)+1; wy=real(m2,RK)*twopi/ylyPlus

    do kc=y1start(3),y1end(3)
      zct = real(kc-1, RK)*dz+dz*0.5_RK
      zct = height-abs(height-zct)
      zplus = utau_guass*zct/xnu
      zbar = zct / height
      do jc=y1start(2),y1end(2)
        yct = yc(jc); yplus=utau_guass * yct / xnu
        do ic=y1start(1),y1end(1)
          xct = real(ic-1, RK)*dx + dx * 0.5_RK
          xplus = utau_guass * xct / xnu
          ux(ic,jc,kc) = uBulk*zbar*exp(-4.5_RK*zbar*zbar)*cos(wy*yplus)*Deviation(ic,jc,kc)
          uy(ic,jc,kc) = uBulk*zbar*exp(-4.5_RK*zbar*zbar)*sin(wx*xplus)*Deviation(ic,jc,kc)
          ux(ic,jc,kc) = ux(ic,jc,kc)+ 3.0_RK*uBulk*(zbar-0.5_RK*zbar*zbar)
        enddo
      enddo
    enddo

    ENDIF

    if(abs(ubulk)<1.0E-20_RK) then
      ratiot=0.0_RK
    else
      ratiot=ubulk/CalcUxAver(ux)
      call MPI_Bcast(ratiot,1,real_type,0,MPI_COMM_WORLD,ierr)    
    endif
    ux= ux*ratiot
    Deviation=0.0_RK

  end subroutine InitVelocity_Duct

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitStatVar_Duct
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitStatVar_Duct(); implicit none
    
    ! locals
    integer(8) :: offset
    character(len=:),allocatable::FileStr
    integer::ierr, iUnit, line_num, j, k, int_t, iprec, nfld, ifld


    allocate(SumStat_plane_x(y1start(2):y1end(2), y1start(3):y1end(3), _nDuct_Stat_), Stat=ierr)
    nfstime=0; SumStat_plane_x = 0.0_RK

    if(nrank/=0) return

    FileStr = strip(ResultsDir_) // "statFor" // strip(RunName_) // ".xmf"
    open(newunit=iUnit, file=FileStr,status='replace',form='formatted',IOSTAT=ierr)
    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"InitVisu","Cannot open file: "//trim(FileStr))
    write(iUnit,'(A)') '<?xml version="1.0" ?>'
    write(iUnit,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(iUnit,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(iUnit,'(A)') '<Domain>'

    iprec = RK
    line_num = 8

    write(iUnit,'(A,3I7,A)')'  <Topology name="TOPO" TopologyType="3DRectMesh" Dimensions="',nzc, nyc, 1, '"/>'
    write(iUnit,'(A)')'  <Geometry name="GEO" GeometryType="VXVYVZ">'    
    ! x-grid
    write(iUnit,'(A,I1,A,I5,A)') '    <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',1,'">'
    write(iUnit,'(A, ES15.7)') '    ', 0.0_RK
    write(iUnit,'(A)')'    </DataItem>'

    ! y-grid
    write(iUnit,'(A,I1,A,I5,A)') '    <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nyc,'">'
    write(iUnit,'(A)',advance='no') '    '
    int_t = 0
    do j = 1, nyc
      int_t = int_t + 1
      if(mod(int_t, line_num) == 0) then
        write(iUnit,'(ES15.7)') yc(j)
        if(int_t < nyc) write(iUnit,'(A)',advance='no') '    '
      else
        write(iUnit,'(ES15.7)',advance='no') yc(j)
      endif
    enddo
    if(mod(nyc, line_num) /=0) write(iUnit,*)' '
    write(iUnit,'(A)')'    </DataItem>'

    ! z-grid
    write(iUnit,'(A,I1,A,I5,A)') '    <DataItem Format="XML" DataType="Float" Precision="',iprec,'" Endian="Native" Dimensions="',nzc,'">'
    write(iUnit,'(A)',advance='no') '    '
    int_t = 0
    do k = 1, nzc
      int_t = int_t + 1
      if(mod(int_t, line_num) == 0) then
        write(iUnit,'(ES15.7)') (k-1)*dz+dz*0.5_RK
        if(int_t < nzc) write(iUnit,'(A)',advance='no') '    '
      else
        write(iUnit,'(ES15.7)',advance='no') (k-1)*dz+dz*0.5_RK
      endif
    enddo
    if(mod(nzc, line_num) /= 0) write(iUnit,*) ' '
    write(iUnit,'(A)')'    </DataItem>'
    write(iUnit,'(A)')'  </Geometry>'


    ! Time series
    nfld = (ilast - ifirst + 1)/SaveStat
    write(iUnit,'(A)')'  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    write(iUnit,'(A)')'    <Time TimeType="List">'
    write(iUnit,'(A,I6,A)')'      <DataItem Format="XML" NumberType="Int" Dimensions="',nfld,'">' 
    write(iUnit,'(A)',advance='no')'        '
    do ifld = 1, nfld
      if(mod(ifld, line_num)==0) then
        write(iUnit,'(I10)') ifld*SaveStat +(ifirst-1)
        if(ifld < nfld) write(iUnit,'(A)',advance='no') '        '
      else
        write(iUnit,'(I10)',advance='no') ifld*SaveStat +(ifirst-1)
      endif     
    enddo
    if(mod(nfld, line_num) /=0) write(iUnit,*)' '
    write(iUnit,'(A)') '      </DataItem>'
    write(iUnit,'(A)') '    </Time>'


    ! attribute
    do ifld=ifirst-1 + SaveStat, ilast, SaveStat
      write(iUnit,'(A,I10.10,A)')'    <Grid Name="T',ifld,'" GridType="Uniform">'
      write(iUnit,'(A)')'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      write(iUnit,'(A)')'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

      offset = 0_8
      FileStr = 'stats_Duct_' // int2str(ifld, 10) // '.bin'
      call Write_XDMF_One('u',  offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 01
      call Write_XDMF_One('v',  offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 02
      call Write_XDMF_One('w',  offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 03
      call Write_XDMF_One('p',  offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 04
      call Write_XDMF_One('uu', offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 05
      call Write_XDMF_One('vv', offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 06
      call Write_XDMF_One('ww', offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 07
      call Write_XDMF_One('pp', offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 08
      call Write_XDMF_One('uw', offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 09
      call Write_XDMF_One('uv', offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 10
      call Write_XDMF_One('vw', offset);  offset = offset + int(nyc * nzc * RK, 8)   ! 11
      write(iUnit,'(A)')'    </Grid>'
    enddo

    write(iUnit,'(A)') '  </Grid>'
    write(iUnit,'(A)') '</Domain>'
    write(iUnit,'(A)') '</Xdmf>'
    close(iUnit,IOSTAT=ierr)

  contains

     subroutine Write_XDMF_One(chAttribute, offset)
       implicit none
       character(len=*), intent(in) :: chAttribute
       integer(8), intent(in) :: offset

       write(iUnit,'(A)')'      <Attribute Name="'//trim(chAttribute)//'" Center="Node">'
       write(iUnit,'(A,I1,A,3I7,A,I12,A)') '        <DataItem Format="Binary" DataType="Float" Precision="',iprec, &
                                           '" Endian="Native" Dimensions="', nzc, nyc, 1, '" Seek="', offset, '">'
       write(iUnit,'(A)')'          ' // trim(FileStr)
       write(iUnit,'(A)')'        </DataItem>'
       write(iUnit,'(A)')'      </Attribute>'
     endsubroutine Write_XDMF_One
  endsubroutine InitStatVar_Duct

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clcStat_Duct
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine clcStat_Duct(ux,uy,uz, pr); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz, pr
  
    ! locals
    integer::ic, jc, kc, ids, ip, jp, kp
    character(len=:),allocatable::filename
    real(RK)::infstime, uxp, uyp, uzp, prc, uxc, uyc, uzc
   
    do kc = y1start(3), y1end(3)
      kp = kc + 1
      do jc = y1start(2), y1end(2)
        jp = jc + 1
        do ic = y1start(1), y1end(1)
          ip = ic + 1
          uxp = ux(ic, jc, kc)
          uyp = uy(ic, jc, kc)
          uzp = uz(ic, jc, kc)
          prc = pr(ic, jc, kc)
          uxc = 0.5_RK * (uxp + ux(ip,jc,kc))
          uyc = 0.5_RK * (uyp + uy(ic,jp,kc))
          uzc = 0.5_RK * (uzp + uz(ic,jc,kp))

          ids = 1
          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + uxp;  ids = ids + 1  ! 01
          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + uyp;  ids = ids + 1  ! 02
          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + uzp;  ids = ids + 1  ! 03
          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + prc;  ids = ids + 1  ! 04
          
          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + uxp * uxp;  ids = ids + 1  ! 05
          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + uyp * uyp;  ids = ids + 1  ! 06
          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + uzp * uzp;  ids = ids + 1  ! 07
          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + prc * prc;  ids = ids + 1  ! 08

          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + uxc * uzc;  ids = ids + 1  ! 09
          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + uxc * uyc;  ids = ids + 1  ! 10
          SumStat_plane_x(jc, kc, ids) = SumStat_plane_x(jc, kc, ids) + uyc * uzc;  ids = ids + 1  ! 11
        enddo
      enddo
    enddo
    nfstime= nfstime + 1
    if(mod(itime,SaveStat)/=0) return
    
    infstime = 1.0_RK/real(nfstime,RK)/real(nxc,RK)
    SumStat_plane_x  = infstime*SumStat_plane_x
    filename = strip(ResultsDir_) // 'stats_Duct_' // int2str(itime,10) // '.bin'
    call decomp_2d_write_plane_reduce(SumStat_plane_x, 2, 1, _nDuct_Stat_, MPI_SUM, filename)

    nfstime=0; SumStat_plane_x=0.0_RK
  end subroutine clcStat_Duct

#undef _nDuct_Stat_
end module f2_FlowType_Duct
