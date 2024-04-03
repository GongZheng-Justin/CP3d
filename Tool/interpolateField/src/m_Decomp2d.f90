module m_Decomp2d
  use MPI
  implicit none
  private

  integer:: mytype_bytes  
  integer,parameter:: mytype = KIND(0.0D0)
  integer,parameter,public:: real_type = MPI_DOUBLE_PRECISION
  
  ! some key global variables
  integer::p_row,p_col
  integer,public:: nrank  ! local MPI rank 
  integer,public:: nproc  ! total number of processors

  ! Define neighboring blocks (including the Bc INFO)
  !   second dimension 4 neighbour numProcessors
  integer,public,dimension(4):: myProcNghBC

  ! staring/ending index and size of data held by current numProcessor
  type,public::DECOMP_INFO
    integer, dimension(3)::yst,yen,ysz
  end type DECOMP_INFO
  TYPE(DECOMP_INFO),public::decomp_main
  integer,dimension(3),public:: ystart,yend,ysize

  type,public:: MatBound
    integer:: xmr, xmm, xme  ! real index, memory index, expand size in xm_dir respectively, xmm=xmr-xme
    integer:: xpr, xpm, xpe  ! real index, memory index, expand size in xp_dir respectively, xpm=xpr+xpe 
    integer:: ymr, ymm, yme  ! real index, memory index, expand size in ym_dir respectively, ymm=ymr-yme
    integer:: ypr, ypm, ype  ! real index, memory index, expand size in yp_dir respectively, ypm=ypr+ype
    integer:: zmr, zmm, zme  ! real index, memory index, expand size in zm_dir respectively, zmm=zmr-zme
    integer:: zpr, zpm, zpe  ! real index, memory index, expand size in zp_dir respectively, zpm=zpr+zpe
  end type MatBound

  type,public:: HaloINFO
    integer :: xmh, xph   ! xm_dir halo, xp_dir halo
    integer :: ymh, yph   ! ym_dir halo, yp_dir halo
    integer :: zmh, zph   ! zm_dir halo, zp_dir halo
  end type HaloINFO

  public:: decomp2d_alloc,decomp2d_updateHalo
  public:: decomp2d_init,decomp2d_readVar,decomp2d_writeVar
contains

  !**********************************************************************
  ! decomp2d_init
  !**********************************************************************
  subroutine decomp2d_init(nxc,nyc,nzc,row,col,BcOpt)
    implicit none
    integer,intent(in)::nxc,nyc,nzc,row,col,BcOpt(6)

    ! locals
    integer::coord1,coord2,idTop,idBottom,idLeft,idRight,ierror

    ! If we have 6 numProcessors = 3 row * 2 col
    !
    !  nrank, coord1, coord2, and neighbor_index  are as follows:
    !   x               x               x               x             
    !   |        4 5    |        2 2    |        0 1    |        6 3 5
    !   |        2 3    |        1 1    |        0 1    |        2 0 1
    !   |_ _ _z  0 1    |_ _ _z  0 0    |_ _ _z  0 1    |_ _ _z  7 4 8

    p_row=row; p_col=col
    if(nproc/=p_row*p_col) call decomp_2d_abort(1,'Invalid 2D processor grid - nproc /= p_row*p_col')
    coord1 = int(nrank/col)
    coord2 = mod(nrank,col)

    ! Firstly, all the boundaries are assumed to be periodic
    idTop     = mod(coord1+1,    row)  ! top
    idBottom  = mod(coord1+row-1,row)  ! bottom 
    idLeft    = mod(coord2+col-1,col)  ! left
    idRight   = mod(coord2+1,    col)  ! right
    myProcNghBC(1) = coord1   * col + idRight; ! myProcNghBC(5) = idTop    * col + idRight
    myProcNghBC(2) = coord1   * col + idLeft ; ! myProcNghBC(6) = idTop    * col + idLeft
    myProcNghBC(3) = idTop    * col + coord2 ; ! myProcNghBC(7) = idBottom * col + idLeft
    myProcNghBC(4) = idBottom * col + coord2 ; ! myProcNghBC(8) = idBottom * col + idRight

    ! Secondly, modify the edge neighbour ids
    IF(BcOpt(1)<0 .and. coord1==0    ) myProcNghBC(4)=BcOpt(1)
    IF(BcOpt(2)<0 .and. coord1==row-1) myProcNghBC(3)=BcOpt(2)
    IF(BcOpt(5)<0 .and. coord2==0    ) myProcNghBC(2)=BcOpt(5)
    IF(BcOpt(6)<0 .and. coord2==col-1) myProcNghBC(1)=BcOpt(6)

    ! Thirdly, determine decomp_main
    decomp_main%yst(2)=1
    decomp_main%yen(2)=nyc
    decomp_main%ysz(2)=nyc
    call distribute_proc(nxc,row,coord1,decomp_main%yst(1),decomp_main%yen(1),decomp_main%ysz(1))
    call distribute_proc(nzc,col,coord2,decomp_main%yst(3),decomp_main%yen(3),decomp_main%ysz(3))
    ystart=decomp_main%yst; yend=decomp_main%yen; ysize=decomp_main%ysz
    
    call MPI_TYPE_SIZE(real_type,mytype_bytes,ierror)
  end subroutine decomp2d_init

  !**********************************************************************
  ! decomp2d_alloc
  !**********************************************************************
  subroutine decomp2d_alloc(var,mb,opt_decomp,opt_global)
    implicit none
    type(MatBound),intent(inout)::mb
    logical,intent(in),optional::opt_global
    TYPE(DECOMP_INFO),intent(in),optional :: opt_decomp
    real(mytype), allocatable, dimension(:,:,:),intent(out) :: var
 
    ! locals
    logical:: global
    integer:: ierror
    TYPE(DECOMP_INFO)::decomp

    if (present(opt_decomp)) then
      decomp=opt_decomp
    else
      decomp=decomp_main
    endif
    if (present(opt_global)) then
      global = opt_global
    else
      global = .false.
    endif

    ! first update the MatBound
    if (global) then
      mb%xmr=decomp%yst(1); mb%xpr=decomp%yen(1)
      mb%ymr=decomp%yst(2); mb%ypr=decomp%yen(2)
      mb%zmr=decomp%yst(3); mb%zpr=decomp%yen(3)
    else
      mb%xmr=1; mb%xpr=decomp%ysz(1)
      mb%ymr=1; mb%ypr=decomp%ysz(2)
      mb%zmr=1; mb%zpr=decomp%ysz(3)
    end if
    mb%xmm= mb%xmr- mb%xme;  mb%xpm= mb%xpr+ mb%xpe
    mb%ymm= mb%ymr- mb%yme;  mb%ypm= mb%ypr+ mb%ype
    mb%zmm= mb%zmr- mb%zme;  mb%zpm= mb%zpr+ mb%zpe 

    allocate(var(mb%xmm:mb%xpm,mb%ymm:mb%ypm,mb%zmm:mb%zpm),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(15,'myallocate-allocate wrong')
  end subroutine decomp2d_alloc

  !**********************************************************************
  ! decomp2d_updateHalo
  !**********************************************************************
  subroutine decomp2d_updateHalo(mat,mb,hi)
    implicit none
    type(MatBound),intent(in)::mb
    type(HaloINFO),intent(in)::hi
    real(mytype), dimension(mb%xmm:mb%xpm,mb%ymm:mb%ypm,mb%zmm:mb%zpm),intent(inout):: mat

    ! locals
    integer:: j,halo_type,ierror,numProcNgh(4) 
    integer:: s1,s2,s3,icount,ilength,ijump,requests(2)
    integer:: xs1,xs2,ys1,ys2,zs1,zs2   ! index for halo sending
    integer:: xr1,xr2,yr1,yr2,zr1,zr2   ! index for halo receiving
    integer,dimension(MPI_STATUS_SIZE,2) :: SRstatus
  
    if(mb%xme<hi%xmh .or. mb%yme<hi%ymh .or. mb%zme<hi%zmh .or.  &
       mb%xpe<hi%xph .or. mb%ype<hi%yph .or. mb%zpe<hi%zph .or.  &
       hi%xmh<0      .or. hi%ymh<0      .or. hi%zmh<0      .or.  &
       hi%xph<0      .or. hi%yph<0      .or. hi%zph<0       ) then
       call decomp_2d_abort(16,'decomp2d_updateHalo-HaloSize error')
    endif

    s1= mb%xpm-mb%xmm+1
    s2= mb%ypm-mb%ymm+1
    s3= mb%zpm-mb%zmm+1
    do j=1,4
      if(myProcNghBC(j)<0) then
        numProcNgh(j)=MPI_PROC_NULL
      else
        numProcNgh(j)=myProcNghBC(j)
      endif
    enddo

    ! step1: receive from ym_dir, and send to yp_dir
    IF(hi%ymh >0 ) THEN
      xr1=mb%xmr;             xr2=mb%xpr 
      yr1=mb%ymr-hi%ymh;      yr2=mb%ymr-1  
      zr1=mb%zmr;             zr2=mb%zpr

      xs1=xr1;                xs2=xr2
      ys1=mb%ypr-hi%ymh+1;    ys2=mb%ypr    
      zs1=zr1;                zs2=zr2
      mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
    ENDIF

    ! step2: receive from yp_dir, and send to ym_dir
    IF(hi%yph >0 ) THEN
      xr1=mb%xmr;             xr2=mb%xpr 
      yr1=mb%ypr+1;           yr2=mb%ypr+hi%yph  
      zr1=mb%zmr;             zr2=mb%zpr

      xs1=xr1;                xs2=xr2
      ys1=mb%ymr;             ys2=mb%ymr+hi%yph-1   
      zs1=zr1;                zs2=zr2
      mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
    ENDIF

    ! step3: receive from xm_dir, and send to xp_dir
    IF(hi%xmh > 0) THEN
      xr1=mb%xmr-hi%xmh;      xr2=mb%xmr-1
      yr1=mb%ymm;             yr2=mb%ypm
      zr1=mb%zmm;             zr2=mb%zpm
        
      xs1=mb%xpr-hi%xmh+1;    xs2=mb%xpr
      ys1=yr1;                ys2=yr2    
      zs1=zr1;                zs2=zr2
      if(numProcNgh(3)==nrank) then  ! neighbour is nrank itself, and numProcNgh(4)==nrank also !!!
        mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
      else
        icount = s2*s3
        ilength= hi%xmh
        ijump  = s1
        call MPI_TYPE_VECTOR(icount,ilength,ijump, real_type, halo_type, ierror)
        call MPI_TYPE_COMMIT(halo_type, ierror)
        call MPI_IRECV( mat(xr1,yr1,zr1),1,halo_type,numProcNgh(4),1,MPI_COMM_WORLD,requests(1),ierror)
        call MPI_ISSEND(mat(xs1,ys1,zs1),1,halo_type,numProcNgh(3),1,MPI_COMM_WORLD,requests(2),ierror)
        call MPI_WAITALL(2,requests,SRstatus,ierror)
        call MPI_TYPE_FREE(halo_type, ierror)
      endif
    ENDIF

    ! step4: receive from xp_dir, and send to xm_dir
    IF(hi%xph > 0) THEN
      xr1=mb%xpr+1;           xr2=mb%xpr+hi%xph
      yr1=mb%ymm;             yr2=mb%ypm
      zr1=mb%zmm;             zr2=mb%zpm
        
      xs1=mb%xmr;             xs2=mb%xmr+hi%xph-1
      ys1=yr1;                ys2=yr2    
      zs1=zr1;                zs2=zr2
      if(numProcNgh(4)==nrank) then  ! neighbour is nrank itself, and numProcNgh(3)==nrank also !!!
        mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
      else
        icount = s2*s3
        ilength= hi%xph
        ijump  = s1
        call MPI_TYPE_VECTOR(icount, ilength, ijump, real_type, halo_type, ierror)
        call MPI_TYPE_COMMIT(halo_type, ierror)
        call MPI_IRECV( mat(xr1,yr1,zr1),1,halo_type,numProcNgh(3),2,MPI_COMM_WORLD,requests(1),ierror)
        call MPI_ISSEND(mat(xs1,ys1,zs1),1,halo_type,numProcNgh(4),2,MPI_COMM_WORLD,requests(2),ierror)
        call MPI_WAITALL(2,requests,SRstatus,ierror)
        call MPI_TYPE_FREE(halo_type, ierror)
      endif
    ENDIF

    ! step5: receive from zm_dir, and send to zp_dir
    IF(hi%zmh > 0) THEN
      xr1=mb%xmm;             xr2=mb%xpm
      yr1=mb%ymm;             yr2=mb%ypm
      zr1=mb%zmr-hi%zmh;      zr2=mb%zmr-1
        
      xs1=xr1;                xs2=xr2
      ys1=yr1;                ys2=yr2    
      zs1=mb%zpr-hi%zmh+1;    zs2=mb%zpr
      if(numProcNgh(1)==nrank) then  ! neighbour is nrank itself, and numProcNgh(2)==nrank also !!!
        mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
      else
        icount=s1*s2*hi%zmh
        call MPI_IRECV( mat(xr1,yr1,zr1),icount,real_type,numProcNgh(2),3,MPI_COMM_WORLD,requests(1),ierror)
        call MPI_ISSEND(mat(xs1,ys1,zs1),icount,real_type,numProcNgh(1),3,MPI_COMM_WORLD,requests(2),ierror)
        call MPI_WAITALL(2,requests,SRstatus,ierror)
      endif
    ENDIF

    ! step6: receive from zp_dir, and send to zm_dir
    IF(hi%zph > 0) THEN
      xr1=mb%xmm;             xr2=mb%xpm
      yr1=mb%ymm;             yr2=mb%ypm
      zr1=mb%zpr+1;           zr2=mb%zpr+hi%zph
        
      xs1=xr1;                xs2=xr2
      ys1=yr1;                ys2=yr2    
      zs1=mb%zmr;             zs2=mb%zmr+hi%zph-1
      if(numProcNgh(2)==nrank) then  ! neighbour is nrank itself, and numProcNgh(1)==nrank also !!!
        mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
      else
        icount=s1*s2*hi%zph
        call MPI_IRECV( mat(xr1,yr1,zr1),icount,real_type,numProcNgh(1),4,MPI_COMM_WORLD,requests(1),ierror)
        call MPI_ISSEND(mat(xs1,ys1,zs1),icount,real_type,numProcNgh(2),4,MPI_COMM_WORLD,requests(2),ierror)
        call MPI_WAITALL(2,requests,SRstatus,ierror)
      endif
     ENDIF
  end subroutine decomp2d_updateHalo

  !**********************************************************************
  !   - distibutes grid points in one dimension
  !   - handles uneven distribution properly 
  !**********************************************************************
  subroutine distribute_proc(nDataLen,numProc,pcoord,st,en,sz)
    implicit none
    integer,intent(in) ::nDataLen,numProc,pcoord
    integer,intent(out)::st,en,sz

    ! locals
    integer:: i,size1,nl,nu
    integer,dimension(:),allocatable::stArr,enArr,szArr 

    allocate(stArr(0:numProc-1),enArr(0:numProc-1),szArr(0:numProc-1))
    size1=nDataLen/numProc
    nu = nDataLen-size1 * numProc
    nl = numProc -nu
    stArr(0)= 1
    szArr(0)= size1
    enArr(0)= size1
    do i=1,nl-1
     stArr(i)= stArr(i-1) + size1
     szArr(i)= size1
     enArr(i)= enArr(i-1) + size1
    end do
    size1 = size1 + 1
    do i=nl,numProc-1
     stArr(i) = enArr(i-1) + 1
     szArr(i) = size1
     enArr(i) = enArr(i-1) + size1
    end do
    enArr(numProc-1)= nDataLen 
    szArr(numProc-1)= nDataLen-stArr(numProc-1)+1

    st=stArr(pcoord);en=enArr(pcoord);sz=szArr(pcoord)
    deallocate(stArr,enArr,szArr)
  end subroutine distribute_proc

  !**********************************************************************
  ! decomp2d_readVar
  !**********************************************************************
  subroutine decomp2d_readVar(fh,disp,var,opt_decomp)
    implicit none
    integer,intent(in)::fh
    integer(MPI_OFFSET_KIND),intent(inout)::disp
    real(mytype), dimension(:,:,:), intent(out)::var
    TYPE(DECOMP_INFO),intent(in),optional::opt_decomp

    ! locals
    TYPE(DECOMP_INFO)::decomp
    integer::ierror,newtype,sizes(3),subsizes(3),starts(3)

    if(present(opt_decomp)) then
      decomp=opt_decomp
    else
      decomp=decomp_main
    endif

    ! Create file type and set file view
    sizes= decomp%ysz; subsizes= decomp%ysz; starts= decomp%yst-1
    call MPI_ALLREDUCE(decomp%ysz(1),sizes(1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror); sizes(1)=sizes(1)/p_col
    call MPI_ALLREDUCE(decomp%ysz(3),sizes(3),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror); sizes(3)=sizes(3)/p_row
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,real_type,newtype,ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh,var,subsizes(1)*subsizes(2)*subsizes(3),real_type,MPI_STATUS_IGNORE,ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for the next read operation
    disp = disp + int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)
  end subroutine decomp2d_readVar

  !**********************************************************************
  ! decomp2d_writeVar
  !**********************************************************************
  subroutine decomp2d_writeVar(fh,disp,var,opt_decomp)
    implicit none
    integer,intent(in)::fh
    integer(MPI_OFFSET_KIND),intent(inout)::disp
    real(mytype), dimension(:,:,:), intent(in)::var
    TYPE(DECOMP_INFO),intent(in),optional::opt_decomp

    ! locals
    TYPE(DECOMP_INFO)::decomp
    integer::ierror,newtype,sizes(3),subsizes(3),starts(3)

    if(present(opt_decomp)) then
      decomp=opt_decomp
    else
      decomp=decomp_main
    endif

    ! Create file type and set file view
    sizes= decomp%ysz; subsizes= decomp%ysz; starts= decomp%yst-1
    call MPI_ALLREDUCE(decomp%ysz(1),sizes(1),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror);sizes(1)=sizes(1)/p_col
    call MPI_ALLREDUCE(decomp%ysz(3),sizes(3),1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror);sizes(3)=sizes(3)/p_row
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,real_type,newtype,ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,real_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,var,subsizes(1)*subsizes(2)*subsizes(3),real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for the next write operation
    disp = disp + int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8) ! modified by gongzheng, 2021-02-18
  end subroutine decomp2d_writeVar

  !******************************************************************
  ! Error handling
  !******************************************************************
  subroutine decomp_2d_abort(errorcode,msg)
    implicit none
    integer,intent(in)::errorcode
    character(len=*),intent(in)::msg

    ! locals
    integer::ierror
    if(nrank==0) then
      write(*,*) '2DECOMP ERROR - errorcode: ', errorcode
      write(*,*) 'ERROR MESSAGE: ' // msg
    endif
    call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
  end subroutine decomp_2d_abort
end module m_Decomp2d
