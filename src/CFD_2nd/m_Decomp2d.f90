module m_decomp2d
  use MPI
  implicit none
  private
#define DECOMP_TRANSPOSE_X1Y1 1
#define DECOMP_TRANSPOSE_Y1Z1 2
#define DECOMP_TRANSPOSE_Z1X2 3
#define DECOMP_TRANSPOSE_X2Y2 4
#define DECOMP_TRANSPOSE_Y2Z2 5
#define DECOMP_TRANSPOSE_Z2X1 6
       
  integer,parameter,public:: mytype = KIND(0.0D0)
  integer,parameter,public:: real_type = MPI_DOUBLE_PRECISION
#ifdef SAVE_SINGLE
  integer,parameter,public:: mytype_save = KIND(0.0)
  integer,parameter,public:: real_type_save = MPI_REAL
#else
  integer,parameter,public:: mytype_save = KIND(0.0D0)
  integer,parameter,public:: real_type_save = MPI_DOUBLE_PRECISION
#endif
  integer,public:: mytype_bytes

  ! some key global variables
  integer,public:: nrank  ! local MPI rank 
  integer,public:: nproc  ! total number of processors

  ! parameters for 2D Cartesian topology
  logical::IsInitTopoEnd=.false.
  integer::np_row,np_col,pcoord(2)
  integer,allocatable,dimension(:)::tran_count,tran_disp
  integer,public::DECOMP_2D_COMM_ROW,DECOMP_2D_COMM_COL

  ! staring/ending index and size of data held by current processor
  ! duplicate 'decomp_main', needed by apps to define data structure 
  integer,dimension(3),public::x1start,x1end,x1size  ! x-pencil
  integer,dimension(3),public::y1start,y1end,y1size  ! y-pencil
  integer,dimension(3),public::z1start,z1end,z1size  ! z-pencil
  integer,dimension(3),public::x2start,x2end,x2size  ! x-pencil
  integer,dimension(3),public::y2start,y2end,y2size  ! y-pencil
  integer,dimension(3),public::z2start,z2end,z2size  ! z-pencil
  
  ! derived type to store decomposition info for a given global data size
  type,public::decomp_info
    integer::nrow,ncol
    logical,dimension(6)::Initialize
    integer,dimension(3)::x1st,x1en,x1sz  ! x-pencil
    integer,dimension(3)::y1st,y1en,y1sz  ! y-pencil
    integer,dimension(3)::z1st,z1en,z1sz  ! z-pencil
    integer,dimension(3)::x2st,x2en,x2sz  ! x-pencil
    integer,dimension(3)::y2st,y2en,y2sz  ! y-pencil
    integer,dimension(3)::z2st,z2en,z2sz  ! z-pencil
  
    integer,allocatable,dimension(:)::xtype_x1y1,ytype_x1y1
    integer,allocatable,dimension(:)::ytype_y1z1,ztype_y1z1
    integer,allocatable,dimension(:)::ztype_z1x2,xtype_z1x2
    integer,allocatable,dimension(:)::xtype_x2y2,ytype_x2y2
    integer,allocatable,dimension(:)::ytype_y2z2,ztype_y2z2
    integer,allocatable,dimension(:)::ztype_z2x1,xtype_z2x1    
  END TYPE decomp_info
  type(decomp_info),target::decomp_main ! main (default) decomposition information

  ! Define neighboring blocks (including the Bc info)
  !   first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
  !   second dimension 4 neighbour processors
  integer,public,dimension(3,4):: myProcNghBC
  type,public:: MatBound
    integer:: pencil
    integer:: xmr, xmm, xme  ! real index, memory index, expand size in xm_dir respectively, xmm=xmr-xme
    integer:: xpr, xpm, xpe  ! real index, memory index, expand size in xp_dir respectively, xpm=xpr+xpe 
    integer:: ymr, ymm, yme  ! real index, memory index, expand size in ym_dir respectively, ymm=ymr-yme
    integer:: ypr, ypm, ype  ! real index, memory index, expand size in yp_dir respectively, ypm=ypr+ype
    integer:: zmr, zmm, zme  ! real index, memory index, expand size in zm_dir respectively, zmm=zmr-zme
    integer:: zpr, zpm, zpe  ! real index, memory index, expand size in zp_dir respectively, zpm=zpr+zpe
  end type MatBound

  type,public:: HaloInfo
    integer :: pencil
    integer :: xmh, xph   ! xm_dir halo, xp_dir halo
    integer :: ymh, yph   ! ym_dir halo, yp_dir halo
    integer :: zmh, zph   ! zm_dir halo, zp_dir halo
  end type HaloInfo

  ! public user routines
  public::decomp_2d_init,decomp_2d_finalize,decomp_info_init,myallocate, &
          decomp_info_finalize, decomp_2d_write_var, decomp_2d_read_var, &
          decomp_2d_write_plane,decomp_2d_write_every, myupdate_halo,    &
          transpose_x1_to_y1, transpose_y1_to_z1, transpose_z1_to_x2,    &
          transpose_x2_to_y2, transpose_y2_to_z2, transpose_z2_to_x1,    &
          transpose_x1_to_z2, transpose_z2_to_y2, transpose_y2_to_x2,    &
          transpose_x2_to_z1, transpose_z1_to_y1, transpose_y1_to_x1
contains

  !******************************************************************
  ! Routine to be called by applications to initialise this library
  !******************************************************************
  subroutine decomp_2d_init(nx,ny,nz,iproc,p_row,p_col,mainPencil,BcOptions)
    implicit none
    integer,intent(inout)::p_row,p_col
    integer,intent(in)::nx,ny,nz,mainPencil,iproc,BcOptions(6)

    ! locals
    integer::errorcode,ierror,nrow,ncol

    nrow=p_row
    ncol=p_col
    nproc=iproc
    IsInitTopoEnd=.false.
    if(p_row==0 .and. p_col==0) then
      call best_2d_grid(nx,ny,nz,mainPencil,iproc,nrow,ncol,BcOptions)
    else
      if(nproc/=p_row*p_col) then
        errorcode = 1
        call decomp_2d_abort(errorcode,'Invalid 2D processor grid - nproc /= p_row*p_col')
      endif
    endif
    call decomp_info_init(nx,ny,nz,decomp_main,[nrow,ncol])
    call Init_ProcNeighbour(nrow,ncol,BcOptions)
    IsInitTopoEnd=.true.
        
    ! make a copy of the decomposition information 
    np_row = nrow; p_row = nrow
    np_col = ncol; p_col = ncol
    x1start= decomp_main%x1st
    y1start= decomp_main%y1st
    z1start= decomp_main%z1st
    x1end  = decomp_main%x1en
    y1end  = decomp_main%y1en
    z1end  = decomp_main%z1en
    x1size = decomp_main%x1sz
    y1size = decomp_main%y1sz
    z1size = decomp_main%z1sz
    x2start= decomp_main%x2st
    y2start= decomp_main%y2st
    z2start= decomp_main%z2st
    x2end  = decomp_main%x2en
    y2end  = decomp_main%y2en
    z2end  = decomp_main%z2en
    x2size = decomp_main%x2sz
    y2size = decomp_main%y2sz
    z2size = decomp_main%z2sz
    call MPI_TYPE_SIZE(real_type,mytype_bytes,ierror)
  end subroutine decomp_2d_init

  !******************************************************************
  ! Function to define globle domain of any size, distribute it, and
  ! then transpose data among pencils.
  ! - generate 2D decomposition details as defined in decomp_info
  !   using its own decomp_info object
  !******************************************************************
  subroutine decomp_info_init(nx,ny,nz,decomp,np_row_col,initialize)
    implicit none
    integer,intent(in)::nx,ny,nz
    type(decomp_info),intent(inout)::decomp
    integer,dimension(2),intent(in),optional::np_row_col
    logical,dimension(6),intent(in),optional::initialize    

    ! locals
    integer::sizes(3),subsizes(3),starts(3),istart
    integer::k,ksize,ierror,DECOMP_2D_COMM_CART,p_row,p_col
    integer,allocatable,dimension(:)::xrstart,xrend,xrsize
    integer,allocatable,dimension(:)::xcstart,xcend,xcsize
    integer,allocatable,dimension(:)::yrstart,yrend,yrsize
    integer,allocatable,dimension(:)::ycstart,ycend,ycsize
    integer,allocatable,dimension(:)::zrstart,zrend,zrsize
    integer,allocatable,dimension(:)::zcstart,zcend,zcsize
        
    if(present(np_row_col)) then
      if(IsInitTopoEnd) then
        call decomp_2d_abort(1,'initialize topology end. New np_row_col is forbidden !!!')      
      endif
      p_row=np_row_col(1)
      p_col=np_row_col(2)
    else
      p_row=np_row
      p_col=np_col
    endif
    decomp%nrow=p_row
    decomp%ncol=p_col
    if(.not.IsInitTopoEnd) then
      call MPI_CART_CREATE(MPI_COMM_WORLD,2,[p_row,p_col],[.false.,.false.],.false.,DECOMP_2D_COMM_CART, ierror)
      call MPI_CART_COORDS(DECOMP_2D_COMM_CART,nrank,2,pcoord,ierror)
      call MPI_CART_SUB(DECOMP_2D_COMM_CART,[.true.,.false.],DECOMP_2D_COMM_COL,ierror)
      call MPI_CART_SUB(DECOMP_2D_COMM_CART,[.false.,.true.],DECOMP_2D_COMM_ROW,ierror)
      call MPI_COMM_FREE(DECOMP_2D_COMM_CART,ierror)
    endif
    if(present(initialize)) then
      decomp%initialize=initialize
    else
      decomp%initialize=.true.
    endif

    ! Verify the global size can actually be distributed as pencils 
    ierror=0
    if(decomp%initialize(DECOMP_TRANSPOSE_X1Y1)) then
      if(min(nx,ny)<p_row .or. nz<p_col) ierror=1
    endif
    if(decomp%initialize(DECOMP_TRANSPOSE_Y1Z1)) then
      if(min(ny,nz)<p_col .or. nx<p_row) ierror=2  
    endif
    if(decomp%initialize(DECOMP_TRANSPOSE_Z1X2)) then
      if(min(nz,nx)<p_row .or. ny<p_col) ierror=3     
    endif    
    if(decomp%initialize(DECOMP_TRANSPOSE_X2Y2)) then
      if(min(nx,ny)<p_col .or. nz<p_row) ierror=4 
    endif
    if(decomp%initialize(DECOMP_TRANSPOSE_Y2Z2)) then
      if(min(ny,nz)<p_row .or. nx<p_col) ierror=5   
    endif
    if(decomp%initialize(DECOMP_TRANSPOSE_Z2X1)) then
      if(min(nz,nx)<p_col .or. ny<p_row) ierror=6
    endif            
    if(ierror/=0) then
      call decomp_2d_abort(ierror,'Invalid 2D processor grid !!!')
    endif

    ! distribute mesh points
    allocate(xrstart(p_row), xrend(p_row), xrsize(p_row))
    allocate(xcstart(p_col), xcend(p_col), xcsize(p_col))
    allocate(yrstart(p_row), yrend(p_row), yrsize(p_row))
    allocate(ycstart(p_col), ycend(p_col), ycsize(p_col))
    allocate(zrstart(p_row), zrend(p_row), zrsize(p_row))
    allocate(zcstart(p_col), zcend(p_col), zcsize(p_col)) 

    call distribute(nx,p_row,xrstart,xrend); xrsize=xrend-xrstart+1
    call distribute(nx,p_col,xcstart,xcend); xcsize=xcend-xcstart+1
    call distribute(ny,p_row,yrstart,yrend); yrsize=yrend-yrstart+1
    call distribute(ny,p_col,ycstart,ycend); ycsize=ycend-ycstart+1
    call distribute(nz,p_row,zrstart,zrend); zrsize=zrend-zrstart+1
    call distribute(nz,p_col,zcstart,zcend); zcsize=zcend-zcstart+1
    decomp%x1st= [1,  yrstart(pcoord(1)+1), zcstart(pcoord(2)+1)]
    decomp%x1en= [nx, yrend(pcoord(1)+1),   zcend(pcoord(2)+1)  ]
    decomp%y1st= [xrstart(pcoord(1)+1), 1,  zcstart(pcoord(2)+1)]
    decomp%y1en= [xrend(pcoord(1)+1),   ny, zcend(pcoord(2)+1)  ]
    decomp%z1st= [xrstart(pcoord(1)+1), ycstart(pcoord(2)+1), 1 ]
    decomp%z1en= [xrend(pcoord(1)+1),   ycend(pcoord(2)+1),   nz]
    decomp%x2st= [1,  ycstart(pcoord(2)+1), zrstart(pcoord(1)+1)]
    decomp%x2en= [nx, ycend(pcoord(2)+1),   zrend(pcoord(1)+1)  ]
    decomp%y2st= [xcstart(pcoord(2)+1), 1,  zrstart(pcoord(1)+1)]
    decomp%y2en= [xcend(pcoord(2)+1),   ny, zrend(pcoord(1)+1)  ]
    decomp%z2st= [xcstart(pcoord(2)+1), yrstart(pcoord(1)+1), 1 ]
    decomp%z2en= [xcend(pcoord(2)+1),   yrend(pcoord(1)+1),   nz]
    decomp%x1sz= decomp%x1en-decomp%x1st+1
    decomp%y1sz= decomp%y1en-decomp%y1st+1
    decomp%z1sz= decomp%z1en-decomp%z1st+1
    decomp%x2sz= decomp%x2en-decomp%x2st+1
    decomp%y2sz= decomp%y2en-decomp%y2st+1
    decomp%z2sz= decomp%z2en-decomp%z2st+1
    deallocate(xrend,  xcend,  yrend,  ycend,  zrend,  zcend)
    deallocate(xrstart,xcstart,yrstart,ycstart,zrstart,zcstart)
    
    ! x1 <-> y1
    if(decomp%initialize(DECOMP_TRANSPOSE_X1Y1)) then
      allocate(decomp%xtype_x1y1(p_row),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 2 !!!')
      sizes=decomp%x1sz
      subsizes=decomp%x1sz
      starts=0; istart=0
      do k=1,p_row
        ksize=xrsize(k)
        starts(1)=istart
        subsizes(1)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%xtype_x1y1(k),ierror)
        call MPI_TYPE_COMMIT(decomp%xtype_x1y1(k),ierror)      
      enddo
      allocate(decomp%ytype_x1y1(p_row),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 3 !!!')
      sizes=decomp%y1sz
      subsizes=decomp%y1sz
      starts=0; istart=0
      do k=1,p_row
        ksize=yrsize(k)
        starts(2)=istart
        subsizes(2)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%ytype_x1y1(k),ierror)
        call MPI_TYPE_COMMIT(decomp%ytype_x1y1(k),ierror) 
      enddo
    endif
    
    ! y1 <-> z1
    if(decomp%initialize(DECOMP_TRANSPOSE_Y1Z1)) then
      allocate(decomp%ytype_y1z1(p_col),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 4 !!!')
      sizes=decomp%y1sz
      subsizes=decomp%y1sz
      starts=0; istart=0
      do k=1,p_col
        ksize=ycsize(k)
        starts(2)=istart
        subsizes(2)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%ytype_y1z1(k),ierror)
        call MPI_TYPE_COMMIT(decomp%ytype_y1z1(k),ierror) 
      enddo
      allocate(decomp%ztype_y1z1(p_col),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 5 !!!')
      sizes=decomp%z1sz
      subsizes=decomp%z1sz
      starts=0; istart=0  
      do k=1,p_col
        ksize=zcsize(k)
        starts(3)=istart
        subsizes(3)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%ztype_y1z1(k),ierror)
        call MPI_TYPE_COMMIT(decomp%ztype_y1z1(k),ierror)       
      enddo
    endif 
        
    ! z1 <-> x2
    if(decomp%initialize(DECOMP_TRANSPOSE_Z1X2)) then
      allocate(decomp%ztype_z1x2(p_row),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 6 !!!')
      sizes=decomp%z1sz
      subsizes=decomp%z1sz
      starts=0; istart=0
      do k=1,p_row
        ksize=zrsize(k)
        starts(3)=istart
        subsizes(3)=ksize
        istart=istart+ksize  
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%ztype_z1x2(k),ierror)
        call MPI_TYPE_COMMIT(decomp%ztype_z1x2(k),ierror) 
      enddo
      allocate(decomp%xtype_z1x2(p_row),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 7 !!!')
      sizes=decomp%x2sz
      subsizes=decomp%x2sz
      starts=0; istart=0
      do k=1,p_row
        ksize=xrsize(k)
        starts(1)=istart
        subsizes(1)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%xtype_z1x2(k),ierror)
        call MPI_TYPE_COMMIT(decomp%xtype_z1x2(k),ierror)      
      enddo
    endif
    
    ! x2 <-> y2
    if(decomp%initialize(DECOMP_TRANSPOSE_X2Y2)) then
      allocate(decomp%xtype_x2y2(p_col),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 8 !!!')
      sizes=decomp%x2sz
      subsizes=decomp%x2sz
      starts=0; istart=0
      do k=1,p_col
        ksize=xcsize(k)
        starts(1)=istart
        subsizes(1)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%xtype_x2y2(k),ierror)
        call MPI_TYPE_COMMIT(decomp%xtype_x2y2(k),ierror)       
      enddo
      allocate(decomp%ytype_x2y2(p_col),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 9 !!!')
      sizes=decomp%y2sz
      subsizes=decomp%y2sz
      starts=0; istart=0
      do k=1,p_col
        ksize=ycsize(k)
        starts(2)=istart
        subsizes(2)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%ytype_x2y2(k),ierror)
        call MPI_TYPE_COMMIT(decomp%ytype_x2y2(k),ierror)       
      enddo 
    endif
       
    ! y2 <-> z2
    if(decomp%initialize(DECOMP_TRANSPOSE_Y2Z2)) then
      allocate(decomp%ytype_y2z2(p_row),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 10 !!!')
      sizes=decomp%y2sz
      subsizes=decomp%y2sz
      starts=0; istart=0  
      do k=1,p_row
        ksize=yrsize(k)
        starts(2)=istart
        subsizes(2)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%ytype_y2z2(k),ierror)
        call MPI_TYPE_COMMIT(decomp%ytype_y2z2(k),ierror)       
      enddo
      allocate(decomp%ztype_y2z2(p_row),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 11 !!!')
      sizes=decomp%z2sz
      subsizes=decomp%z2sz
      starts=0; istart=0  
      do k=1,p_row
        ksize=zrsize(k)
        starts(3)=istart
        subsizes(3)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%ztype_y2z2(k),ierror)
        call MPI_TYPE_COMMIT(decomp%ztype_y2z2(k),ierror)       
      enddo 
    endif   
            
    ! z2 <-> x1
    if(decomp%initialize(DECOMP_TRANSPOSE_Z2X1)) then
      allocate(decomp%ztype_z2x1(p_col),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 12 !!!')
      sizes=decomp%z2sz
      subsizes=decomp%z2sz
      starts=0; istart=0  
      do k=1,p_col
        ksize=zcsize(k)
        starts(3)=istart
        subsizes(3)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%ztype_z2x1(k),ierror)
        call MPI_TYPE_COMMIT(decomp%ztype_z2x1(k),ierror)       
      enddo    
      allocate(decomp%xtype_z2x1(p_col),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(ierror,'decomp_info_init- Allocate Wrong 13 !!!')
      sizes=decomp%x1sz
      subsizes=decomp%x1sz
      starts=0; istart=0   
      do k=1,p_col
        ksize=xcsize(k)
        starts(1)=istart
        subsizes(1)=ksize
        istart=istart+ksize
        call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, &
          real_type,decomp%xtype_z2x1(k),ierror)
        call MPI_TYPE_COMMIT(decomp%xtype_z2x1(k),ierror)       
      enddo
    endif
  end subroutine decomp_info_init

  !******************************************************************
  ! Routine to be called by applications to clean things up
  !******************************************************************
  subroutine decomp_2d_finalize()
    implicit none
    IsInitTopoEnd=.false.
    call decomp_info_finalize(decomp_main)
  end subroutine decomp_2d_finalize
  
  !******************************************************************
  ! Release memory associated with a decomp_info object
  !******************************************************************
  subroutine decomp_info_finalize(decomp)
    implicit none
    type(decomp_info),intent(inout)::decomp

    ! locals
    integer::k,p_row,p_col,ierror

    p_row=decomp%nrow
    p_col=decomp%ncol
    
    ! x1 <-> y1
    if(decomp%initialize(DECOMP_TRANSPOSE_X1Y1)) then
      do k=1,p_row
        call MPI_TYPE_FREE(decomp%xtype_x1y1(k),ierror) 
        call MPI_TYPE_FREE(decomp%ytype_x1y1(k),ierror) 
      enddo
      deallocate(decomp%xtype_x1y1,decomp%ytype_x1y1)
    endif
    
    ! y1 <-> z1
    if(decomp%initialize(DECOMP_TRANSPOSE_Y1Z1)) then
      do k=1,p_col
        call MPI_TYPE_FREE(decomp%ytype_y1z1(k),ierror) 
        call MPI_TYPE_FREE(decomp%ztype_y1z1(k),ierror)
      enddo
      deallocate(decomp%ytype_y1z1,decomp%ztype_y1z1) 
    endif 
        
    ! z1 <-> x2
    if(decomp%initialize(DECOMP_TRANSPOSE_Z1X2)) then
      do k=1,p_row
        call MPI_TYPE_FREE(decomp%ztype_z1x2(k),ierror)
        call MPI_TYPE_FREE(decomp%xtype_z1x2(k),ierror)     
      enddo
      deallocate(decomp%ztype_z1x2,decomp%xtype_z1x2)
    endif
    
    ! x2 <-> y2
    if(decomp%initialize(DECOMP_TRANSPOSE_X2Y2)) then
      do k=1,p_col
        call MPI_TYPE_FREE(decomp%xtype_x2y2(k),ierror)       
        call MPI_TYPE_FREE(decomp%ytype_x2y2(k),ierror)       
      enddo 
      deallocate(decomp%xtype_x2y2,decomp%ytype_x2y2)
    endif
       
    ! y2 <-> z2
    if(decomp%initialize(DECOMP_TRANSPOSE_Y2Z2)) then
      do k=1,p_row
        call MPI_TYPE_FREE(decomp%ytype_y2z2(k),ierror)
        call MPI_TYPE_FREE(decomp%ztype_y2z2(k),ierror)       
      enddo 
      deallocate(decomp%ytype_y2z2,decomp%ztype_y2z2)
    endif   
            
    ! z2 <-> x1
    if(decomp%initialize(DECOMP_TRANSPOSE_Z2X1)) then
      do k=1,p_col
        call MPI_TYPE_FREE(decomp%ztype_z2x1(k),ierror)
        call MPI_TYPE_FREE(decomp%xtype_z2x1(k),ierror)       
      enddo
      deallocate(decomp%ztype_z2x1,decomp%xtype_z2x1)
    endif

    if(.not.IsInitTopoEnd) then
      call MPI_COMM_FREE(DECOMP_2D_COMM_ROW,ierror)
      call MPI_COMM_FREE(DECOMP_2D_COMM_COL,ierror)   
    endif
  end subroutine decomp_info_finalize
 
  !******************************************************************
  ! Distibutes grid points in one dimension
  !******************************************************************
  subroutine distribute(ndata,mproc,istart,iend)
    implicit none
    integer::ndata,mproc,istart(mproc),iend(mproc)

    ! locals
    integer::i,isize,nl

    isize=ndata/mproc
    nl=mproc-mod(ndata,mproc)
    do i=1,mproc
      istart(i)=(i-1)*isize+1+max(i-1-nl,0)
      iend(i)=i*isize+max(i-nl,0)
    enddo
  end subroutine distribute

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

  !******************************************************************
  ! findfactor
  !******************************************************************
  subroutine findfactor(num, factors, nfact)
    implicit none
    integer,intent(in)::num
    integer,intent(out),dimension(*)::factors
    integer,intent(out)::nfact

    ! locals
    integer::i,m

    nfact = 1
    m = int(sqrt(real(num)))
    do i=1,m
      if(num/i*i /= num)cycle
      factors(nfact)= i
      nfact = nfact + 1
    enddo
    nfact = nfact - 1
    if(factors(nfact)**2/=num) then
      do i=nfact+1, 2*nfact
        factors(i) = num/factors(2*nfact-i+1)
      enddo
      nfact = nfact * 2
    else
      do i=nfact+1, 2*nfact-1
        factors(i)= num/factors(2*nfact-i)
      enddo
      nfact = nfact * 2 - 1
    endif
  end subroutine findfactor

  !******************************************************************
  ! transpose_x1_to_y1
  !******************************************************************  
  subroutine transpose_x1_to_y1(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer :: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_X1Y1)) then
      call decomp_2d_abort(1, 'transpose_x1_to_y1 Have Not been initialized')
    endif   
    
    np=decomp%nrow
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(1, 'transpose_x1_to_y1-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(1)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%xtype_x1y1, &
                       dst,tran_count,tran_disp,decomp%ytype_x1y1, DECOMP_2D_COMM_COL, ierror)
    dst(:,decomp%x1st(2):decomp%x1en(2),:)=src(decomp%y1st(1):decomp%y1en(1),:,:)
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_x1_to_y1

  !******************************************************************
  ! transpose_y1_to_x1
  !******************************************************************
  subroutine transpose_y1_to_x1(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_X1Y1)) then
      call decomp_2d_abort(2, 'transpose_y1_to_x1 Have Not been initialized')
    endif

    np=decomp%nrow
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(2, 'transpose_y1_to_x1-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(1)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%ytype_x1y1, &
                       dst,tran_count,tran_disp,decomp%xtype_x1y1, DECOMP_2D_COMM_COL, ierror)
    dst(decomp%y1st(1):decomp%y1en(1),:,:)=src(:,decomp%x1st(2):decomp%x1en(2),:)
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_y1_to_x1

  !******************************************************************
  ! transpose_y1_to_z1
  !******************************************************************
  subroutine transpose_y1_to_z1(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_Y1Z1)) then
      call decomp_2d_abort(3, 'transpose_y1_to_z1 Have Not been initialized')
    endif 
    
    np=decomp%ncol 
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(3, 'transpose_y1_to_z1-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(2)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%ytype_y1z1, &
                       dst,tran_count,tran_disp,decomp%ztype_y1z1, DECOMP_2D_COMM_ROW, ierror)
    dst(:,:,decomp%y1st(3):decomp%y1en(3))=src(:,decomp%z1st(2):decomp%z1en(2),:)
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_y1_to_z1

  !******************************************************************
  ! transpose_z1_to_y1
  !******************************************************************
  subroutine transpose_z1_to_y1(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_Y1Z1)) then
      call decomp_2d_abort(4, 'transpose_z1_to_y1 Have Not been initialized')
    endif
    
    np=decomp%ncol
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(4, 'transpose_z1_to_y1-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(2)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%ztype_y1z1, &
                       dst,tran_count,tran_disp,decomp%ytype_y1z1, DECOMP_2D_COMM_ROW, ierror)
    dst(:,decomp%z1st(2):decomp%z1en(2),:)=src(:,:,decomp%y1st(3):decomp%y1en(3))
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_z1_to_y1

  !******************************************************************
  ! transpose_z1_to_x2
  !******************************************************************
  subroutine transpose_z1_to_x2(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_Z1X2)) then
      call decomp_2d_abort(5, 'transpose_z1_to_x2 Have Not been initialized')
    endif 
    
    np=decomp%nrow
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(5, 'transpose_z1_to_x2-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(1)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%ztype_z1x2, &
                       dst,tran_count,tran_disp,decomp%xtype_z1x2, DECOMP_2D_COMM_COL, ierror)
    dst(decomp%z1st(1):decomp%z1en(1),:,:)=src(:,:,decomp%x2st(3):decomp%x2en(3))
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_z1_to_x2    
 
  !******************************************************************
  ! transpose_x2_to_z1
  !******************************************************************
  subroutine transpose_x2_to_z1(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_Z1X2)) then
      call decomp_2d_abort(6, 'transpose_x2_to_z1 Have Not been initialized')
    endif 
  
    np=decomp%nrow
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(6, 'transpose_x2_to_z1-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(1)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%xtype_z1x2, &
                       dst,tran_count,tran_disp,decomp%ztype_z1x2, DECOMP_2D_COMM_COL, ierror)
    dst(:,:,decomp%x2st(3):decomp%x2en(3))=src(decomp%z1st(1):decomp%z1en(1),:,:)
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_x2_to_z1
    
  !******************************************************************
  ! transpose_x2_to_y2
  !******************************************************************
  subroutine transpose_x2_to_y2(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_X2Y2)) then
      call decomp_2d_abort(7, 'transpose_x2_to_y2 Have Not been initialized')
    endif 
 
    np=decomp%ncol
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(7, 'transpose_x2_to_y2-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(2)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%xtype_x2y2, &
                       dst,tran_count,tran_disp,decomp%ytype_x2y2, DECOMP_2D_COMM_ROW, ierror)
    dst(:,decomp%x2st(2):decomp%x2en(2),:)=src(decomp%y2st(1):decomp%y2en(1),:,:)
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_x2_to_y2

  !******************************************************************
  ! transpose_y2_to_x2
  !******************************************************************
  subroutine transpose_y2_to_x2(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_X2Y2)) then
      call decomp_2d_abort(8, 'transpose_y2_to_x2 Have Not been initialized')
    endif 
 
    np=decomp%ncol
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(8, 'transpose_y2_to_x2-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(2)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%ytype_x2y2, &
                       dst,tran_count,tran_disp,decomp%xtype_x2y2, DECOMP_2D_COMM_ROW, ierror)
    dst(decomp%y2st(1):decomp%y2en(1),:,:)=src(:,decomp%x2st(2):decomp%x2en(2),:)
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_y2_to_x2 
  
  !******************************************************************
  ! transpose_y2_to_z2
  !******************************************************************
  subroutine transpose_y2_to_z2(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_Y2Z2)) then
      call decomp_2d_abort(9, 'transpose_y2_to_z2 Have Not been initialized')
    endif
    
    np=decomp%nrow
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(9, 'transpose_y2_to_z2-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(1)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%ytype_y2z2, &
                       dst,tran_count,tran_disp,decomp%ztype_y2z2, DECOMP_2D_COMM_COL, ierror)
    dst(:,:,decomp%y2st(3):decomp%y2en(3))=src(:,decomp%z2st(2):decomp%z2en(2),:)
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_y2_to_z2

  !******************************************************************
  ! transpose_z2_to_y2
  !******************************************************************
  subroutine transpose_z2_to_y2(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_Y2Z2)) then
      call decomp_2d_abort(10,'transpose_z2_to_y2 Have Not been initialized')
    endif
    
    np=decomp%nrow
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(10, 'transpose_z2_to_y2-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(1)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%ztype_y2z2, &
                       dst,tran_count,tran_disp,decomp%ytype_y2z2, DECOMP_2D_COMM_COL, ierror)
    dst(:,decomp%z2st(2):decomp%z2en(2),:)=src(:,:,decomp%y2st(3):decomp%y2en(3))
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_z2_to_y2 
    
  !******************************************************************
  ! transpose_z2_to_x1
  !******************************************************************
  subroutine transpose_z2_to_x1(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_Z2X1)) then
      call decomp_2d_abort(11,'transpose_z2_to_x1 Have Not been initialized')
    endif
    
    np=decomp%ncol
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(11, 'transpose_z2_to_x1-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(2)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%ztype_z2x1, &
                       dst,tran_count,tran_disp,decomp%xtype_z2x1, DECOMP_2D_COMM_ROW, ierror)
    dst(decomp%z2st(1):decomp%z2en(1),:,:)=src(:,:,decomp%x1st(3):decomp%x1en(3))
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_z2_to_x1

  !******************************************************************
  ! transpose_x1_to_z2
  !******************************************************************
  subroutine transpose_x1_to_z2(src, dst, opt_decomp)
    implicit none
    real(mytype),dimension(:,:,:),intent(in) ::src
    real(mytype),dimension(:,:,:),intent(out)::dst
    type(decomp_info),intent(in),optional,target::opt_decomp

    ! locals
    integer::np,ierror
    type(decomp_info),pointer:: decomp
    
    if(present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if(.not.decomp%Initialize(DECOMP_TRANSPOSE_Z2X1)) then
      call decomp_2d_abort(12,'transpose_x1_to_z2 Have Not been initialized')
    endif
  
    np=decomp%ncol
    allocate(tran_count(np),tran_disp(np),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(12, 'transpose_x1_to_z2-allocate wrong')
    tran_count=1; tran_disp=0
    tran_count(pcoord(2)+1)=0
    
    call MPI_ALLTOALLW(src,tran_count,tran_disp,decomp%xtype_z2x1, &
                       dst,tran_count,tran_disp,decomp%ztype_z2x1, DECOMP_2D_COMM_ROW, ierror)
    dst(:,:,decomp%x1st(3):decomp%x1en(3))=src(decomp%z2st(1):decomp%z2en(1),:,:)
    
    deallocate(tran_count,tran_disp)
    nullify(decomp)
  end subroutine transpose_x1_to_z2
    
  !******************************************************************
  ! Write a 3D array as part of a big MPI-IO file.
  !  'disp' will be updated after the writing operation.
  !******************************************************************
  subroutine decomp_2d_write_var(fh,disp,ipencil,var,opt_decomp)
    implicit none
    integer,intent(in)::fh
    integer,intent(in)::ipencil
    integer(MPI_OFFSET_KIND),intent(inout)::disp
    real(mytype),dimension(:,:,:),intent(in)::var
    type(decomp_info),intent(in),target,optional::opt_decomp

    ! locals
    type(decomp_info),pointer::decomp
    integer::ierror,data_type,newtype
    integer,dimension(3)::sizes,subsizes,starts

    data_type = real_type
    if(present(opt_decomp))then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    sizes(1)= decomp%x1sz(1)
    sizes(2)= decomp%y1sz(2)
    sizes(3)= decomp%z1sz(3)
    if(ipencil==1)then
      subsizes= decomp%x1sz
      starts= decomp%x1st-1
    elseif(ipencil==2)then
      subsizes= decomp%y1sz
      starts= decomp%y1st-1
    elseif(ipencil==3)then
      subsizes= decomp%z1sz
      starts= decomp%z1st-1
    endif
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,data_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,var,subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8) ! Update displacement
  end subroutine decomp_2d_write_var

  !******************************************************************
  ! Read a 3D array as part of a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated after the reading
  !  operation to prepare the reading of next chunk of data.
  !******************************************************************
  subroutine decomp_2d_read_var(fh,disp,ipencil,var,opt_decomp)
    implicit none
    integer,intent(in)::fh
    integer,intent(in)::ipencil
    integer(MPI_OFFSET_KIND),intent(inout)::disp
    real(mytype),dimension(:,:,:),intent(inout)::var
    type(decomp_info),intent(in),target,optional::opt_decomp

    ! locals
    type(decomp_info),pointer::decomp
    integer::ierror,data_type,newtype
    integer,dimension(3)::sizes,subsizes,starts

    data_type = real_type
    if(present(opt_decomp))then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    sizes(1)= decomp%x1sz(1)
    sizes(2)= decomp%y1sz(2)
    sizes(3)= decomp%z1sz(3)
    if(ipencil==1)then
      subsizes= decomp%x1sz
      starts= decomp%x1st-1
    elseif(ipencil==2)then
      subsizes= decomp%y1sz
      starts= decomp%y1st-1
    elseif(ipencil==3)then
      subsizes= decomp%z1sz
      starts= decomp%z1st-1
    endif
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,data_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh,var,subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    disp=disp + int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(mytype_bytes,8)  ! Update displacement
  end subroutine decomp_2d_read_var

  !******************************************************************
  ! Write a 2D slice of the 3D data to a file
  ! It is much easier to implement if all mpi ranks participate I/O.
  ! Transpose the 3D data if necessary.
  !******************************************************************
  subroutine decomp_2d_write_plane(ipencil,var,iplane,n,filename,opt_decomp)
    implicit none
    integer,intent(in)::n        !which plane to write (global coordinate)
    integer,intent(in)::iplane   !x-plane= 1; y-plane= 2; z-plane= 3
    integer,intent(in)::ipencil  !x-pencil=1; y-pencil=2; z-pencil=3
    character(len=*),intent(in)::filename
    real(mytype),dimension(:,:,:),intent(in) :: var
    type(decomp_info),intent(in),target,optional::opt_decomp

    ! lcoals
    type(decomp_info),pointer::decomp
    integer::sizes(3),subsizes(3),starts(3)
    integer::i,j,k,ierror,newtype,fh,data_type
    real(mytype),allocatable,dimension(:,:,:)::wk3d
    real(mytype_save),allocatable,dimension(:,:,:)::wk2d
    
    data_type= real_type_save
    if(present(opt_decomp))then
      decomp=> opt_decomp
    else
      decomp=> decomp_main
    endif
    if(iplane==1)then
      allocate(wk2d(1,decomp%x1sz(2),decomp%x1sz(3)),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(13, 'decomp_2d_write_plane-allocate wrong-1')
      IF(ipencil==1)THEN
        do k=1,decomp%x1sz(3)
          do j=1,decomp%x1sz(2)
            wk2d(1,j,k)=real(var(n,j,k),mytype_save)
          enddo
        enddo
        sizes   =[1,decomp%y1sz(2),decomp%z1sz(3)]
        subsizes=[1,decomp%x1sz(2),decomp%x1sz(3)]
        starts  =[1,decomp%x1st(2),decomp%x1st(3)]-1
      ELSEIF(ipencil==2)then
        allocate(wk3d(decomp%x1sz(1),decomp%x1sz(2),decomp%x1sz(3)),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(13, 'decomp_2d_write_plane-allocate wrong-2')
        call transpose_y1_to_x1(var,wk3d,decomp)
        do k=1,decomp%x1sz(3)
          do j=1,decomp%x1sz(2)
            wk2d(1,j,k)=real(wk3d(n,j,k),mytype_save)
          enddo
        enddo
        deallocate(wk3d)
        sizes   =[1,decomp%y1sz(2),decomp%z1sz(3)]
        subsizes=[1,decomp%x1sz(2),decomp%x1sz(3)]
        starts  =[1,decomp%x1st(2),decomp%x1st(3)]-1
      ELSE
        allocate(wk3d(decomp%x2sz(1),decomp%x2sz(2),decomp%x2sz(3)),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(13, 'decomp_2d_write_plane-allocate wrong-3')
        call transpose_z1_to_x2(var,wk3d,decomp)
        do k=1,decomp%x2sz(3)
          do j=1,decomp%x2sz(2)
            wk2d(1,j,k)=real(wk3d(n,j,k),mytype_save)
          enddo
        enddo
        deallocate(wk3d)
        sizes   =[1,decomp%y1sz(2),decomp%z1sz(3)]
        subsizes=[1,decomp%x2sz(2),decomp%x2sz(3)]
        starts  =[1,decomp%x2st(2),decomp%x2st(3)]-1
      ENDIF
    elseif(iplane==2) then
      allocate(wk2d(decomp%y1sz(1),1,decomp%y1sz(3)),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(13, 'decomp_2d_write_plane-allocate wrong-4')
      IF(ipencil==2)THEN
        do k=1,decomp%y1sz(3)
          do i=1,decomp%y1sz(1)
            wk2d(i,1,k)=real(var(i,n,k),mytype_save)
          enddo
        enddo
      ELSE
        allocate(wk3d(decomp%y1sz(1),decomp%y1sz(2),decomp%y1sz(3)),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(13, 'decomp_2d_write_plane-allocate wrong-5')
        if(ipencil==1)then
          call transpose_x1_to_y1(var,wk3d,decomp)
        elseif(ipencil==3)then
          call transpose_z1_to_y1(var,wk3d,decomp)
        endif
        do k=1,decomp%y1sz(3)
          do i=1,decomp%y1sz(1)
            wk2d(i,1,k)=real(wk3d(i,n,k),mytype_save)
          enddo
        enddo
        deallocate(wk3d)
      ENDIF
      sizes   =(/decomp%x1sz(1),1,decomp%z1sz(3)/)
      subsizes=(/decomp%y1sz(1),1,decomp%y1sz(3)/)
      starts  =(/decomp%y1st(1),1,decomp%y1st(3)/)-1
    elseif(iplane==3)then
      IF(ipencil==1) THEN
        allocate(wk2d(decomp%z2sz(1),decomp%z2sz(2),1),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(13, 'decomp_2d_write_plane-allocate wrong-6')
        allocate(wk3d(decomp%z2sz(1),decomp%z2sz(2),decomp%z2sz(3)),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(13, 'decomp_2d_write_plane-allocate wrong-7')
        call transpose_x1_to_z2(var,wk3d,decomp)
        do j=1,decomp%z2sz(2)
          do i=1,decomp%z2sz(1) 
            wk2d(i,j,1)=real(wk3d(i,j,n),mytype_save)
          enddo
        enddo 
        deallocate(wk3d)
        sizes   =[decomp%x1sz(1),decomp%y1sz(2),1]
        subsizes=[decomp%z2sz(1),decomp%z2sz(2),1]
        starts  =[decomp%z2st(1),decomp%z2st(2),1]-1
      ELSEIF(ipencil==2) THEN
        allocate(wk2d(decomp%z1sz(1),decomp%z1sz(2),1),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(13, 'decomp_2d_write_plane-allocate wrong-8')
        allocate(wk3d(decomp%z1sz(1),decomp%z1sz(2),decomp%z1sz(3)),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(13, 'decomp_2d_write_plane-allocate wrong-9')
        call transpose_y1_to_z1(var,wk3d,decomp)
        do j=1,decomp%z1sz(2)
          do i=1,decomp%z1sz(1) 
            wk2d(i,j,1)=real(wk3d(i,j,n),mytype_save)
          enddo
        enddo
        deallocate(wk3d)
        sizes   =[decomp%x1sz(1),decomp%y1sz(2),1]
        subsizes=[decomp%z1sz(1),decomp%z1sz(2),1]
        starts  =[decomp%z1st(1),decomp%z1st(2),1]-1 
      ELSE
        allocate(wk2d(decomp%z1sz(1),decomp%z1sz(2),1),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(13, 'decomp_2d_write_plane-allocate wrong-10')
        do j=1,decomp%z1sz(2)
          do i=1,decomp%z1sz(1) 
            wk2d(i,j,1)=real(var(i,j,n),mytype_save)
          enddo
        enddo
        sizes   =[decomp%x1sz(1),decomp%y1sz(2),1]
        subsizes=[decomp%z1sz(1),decomp%z1sz(2),1]
        starts  =[decomp%z1st(1),decomp%z1st(2),1]-1
      ENDIF
    endif
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh,ierror)
    call MPI_FILE_SET_SIZE(fh,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,wk2d,subsizes(1)*subsizes(2)*subsizes(3),data_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    deallocate(wk2d)
  end subroutine decomp_2d_write_plane

  !******************************************************************
  ! Write 3D array data for every specified mesh point
  !******************************************************************
  subroutine decomp_2d_write_every(ipencil,var,iskip,jskip,kskip,filename,from1)
    implicit none
    integer,intent(in)::ipencil
    integer,intent(in)::iskip,jskip,kskip 
    character(len=*),intent(in)::filename
    real(mytype),dimension(:,:,:),intent(in)::var
    logical,intent(in)::from1  ! T:save 1,n+1,2n+1...; F:save n,2n,3n...

    ! locals
    real(mytype_save),allocatable,dimension(:,:,:)::wk
    integer::i,j,k,id,jd,kd,ierror,newtype,fh,key,color,newcomm,data_type
    integer,dimension(3)::sizes,subsizes,starts,xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip

    data_type = real_type_save
    skip=(/iskip,jskip,kskip/)
    do i=1,3
      if(from1)then
        xst(i)= (x1start(i)+skip(i)-1)/skip(i)
        yst(i)= (y1start(i)+skip(i)-1)/skip(i)
        zst(i)= (z1start(i)+skip(i)-1)/skip(i)
        if(mod(x1start(i)+skip(i)-1,skip(i))/=0) xst(i)=xst(i)+1
        if(mod(y1start(i)+skip(i)-1,skip(i))/=0) yst(i)=yst(i)+1
        if(mod(z1start(i)+skip(i)-1,skip(i))/=0) zst(i)=zst(i)+1
        xen(i)= (x1end(i)+skip(i)-1)/skip(i)
        yen(i)= (y1end(i)+skip(i)-1)/skip(i)
        zen(i)= (z1end(i)+skip(i)-1)/skip(i)
      else
        xst(i)= x1start(i)/skip(i)
        yst(i)= y1start(i)/skip(i)
        zst(i)= z1start(i)/skip(i)
        if(mod(x1start(i),skip(i))/=0) xst(i)=xst(i)+1
        if(mod(y1start(i),skip(i))/=0) yst(i)=yst(i)+1
        if(mod(z1start(i),skip(i))/=0) zst(i)=zst(i)+1
        xen(i)= x1end(i)/skip(i)
        yen(i)= y1end(i)/skip(i)
        zen(i)= z1end(i)/skip(i)
      endif
      xsz(i)= xen(i)-xst(i)+1
      ysz(i)= yen(i)-yst(i)+1
      zsz(i)= zen(i)-zst(i)+1
    enddo

    ! if 'skip' value is large it is possible that some ranks do not 
    ! contain any points to be written. Subarray constructor requires 
    ! nonzero size so it is not possible to use MPI_COMM_WORLD for IO.
    ! Create a sub communicator for this...
    color=1; key=0  ! rank order doesn't matter
    if(ipencil==1)then
      if(xsz(1)==0.or.xsz(2)==0.or.xsz(3)==0)color=2
    elseif(ipencil==2)then
      if(ysz(1)==0.or.ysz(2)==0.or.ysz(3)==0)color=2
    elseif(ipencil==3)then
      if(zsz(1)==0.or.zsz(2)==0.or.zsz(3)==0)color=2
    endif
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,newcomm,ierror)

    if(color/=1) return ! only ranks in this group do IO collectively
    sizes(1)= xsz(1)
    sizes(2)= ysz(2)
    sizes(3)= zsz(3)
    if(ipencil==1) then
      subsizes= xsz
      starts= xst-1
    elseif(ipencil==2) then
      subsizes= ysz
      starts= yst-1
    elseif(ipencil==3) then
      subsizes= zsz
      starts= zst-1
    endif

    ! Copy data from original array
    IF(ipencil==1)THEN
      allocate(wk(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(14, 'decomp_2d_write_every-allocate wrong-1')
      if(from1)then
        do k=xst(3),xen(3)
          kd=(k-1)*kskip-x1start(3)+2
          do j=xst(2),xen(2)
            jd=(j-1)*jskip-x1start(2)+2
            do i=xst(1),xen(1)
              id=(i-1)*iskip-x1start(1)+2
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
        enddo
      else
        do k=xst(3),xen(3)
          kd=k*kskip-x1start(3)+1
          do j=xst(2),xen(2)
            jd=j*jskip-x1start(2)+1
            do i=xst(1),xen(1)
              id=i*iskip-x1start(1)+1
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
        enddo
      endif   
    ELSEIF(ipencil==2)THEN
      allocate(wk(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(14, 'decomp_2d_write_every-allocate wrong-2')
      if(from1)then
        do k=yst(3),yen(3)
          kd=(k-1)*kskip-y1start(3)+2
          do j=yst(2),yen(2)
            jd=(j-1)*jskip-y1start(2)+2
            do i=yst(1),yen(1)
              id=(i-1)*iskip-y1start(1)+2
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
        enddo
      else
        do k=yst(3),yen(3)
          kd=k*kskip-y1start(3)+1
          do j=yst(2),yen(2)
            jd=j*jskip-y1start(2)+1
            do i=yst(1),yen(1)
              id=i*iskip-y1start(1)+1
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
        enddo
      endif
    ELSEIF(ipencil==3)THEN
      allocate(wk(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)),stat=ierror)
      if(ierror/=0) call decomp_2d_abort(14, 'decomp_2d_write_every-allocate wrong-3')
      if(from1)then
        do k=zst(3),zen(3)
          kd=(k-1)*kskip-z1start(3)+2
          do j=zst(2),zen(2)
            jd=(j-1)*jskip-z1start(2)+2
            do i=zst(1),zen(1)
              id=(i-1)*iskip-z1start(1)+2
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
         enddo
      else
        do k=zst(3),zen(3)
          kd=k*kskip-z1start(3)+1
          do j=zst(2),zen(2)
            jd=j*jskip-z1start(2)+1
            do i=zst(1),zen(1)
              id=i*iskip-z1start(1)+1
              wk(i,j,k)=real(var(id,jd,kd),mytype_save)
            enddo
          enddo
        enddo
      endif
    ENDIF
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN, data_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(newcomm,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,fh,ierror)
    call MPI_FILE_SET_SIZE(fh,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
    call MPI_BARRIER(newcomm,ierror)
    call MPI_FILE_SET_VIEW(fh,0_MPI_OFFSET_KIND,data_type,newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh,wk,subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    deallocate(wk)
  end subroutine decomp_2d_write_every

  !**********************************************************************
  ! Init_ProcNeighbour
  !**********************************************************************
  subroutine Init_ProcNeighbour(row,col,BcOpt)
    implicit none
    integer,intent(in)::row,col,BcOpt(6)

    !locals
    integer :: i,coord1,coord2,idTop,idBottom,idLeft,idRight

    ! ------------------------- neigbor information begins------------------------
    ! 
    ! 2D domain decomposition method in Decomp2D
    ! 
    !   (from left to right:  x1, y1, z1, x2, y2, z2 pencil respectively)
    ! 
    !   If we have 6 processors = 3 row * 2 col
    ! 
    !     The arrangement of the smubomains(nrank) is as follow:
    !       y             x             x             y              x              x
    !       |       4 5   |       4 5   |       4 5   |              |              |
    !       |       2 3   |       2 3   |       2 3   |       1 3 5  |       1 3 5  |       1 3 5
    !       |_ _ _z 0 1   |_ _ _z 0 1   |_ _ _y 0 1   |_ _ _z 0 2 4  |_ _ _z 0 2 4  |_ _ _y 0 2 4
    !
    !     the arrangement of the coord1 is as follow:
    !       y             x             x
    !       |       2 2   |       2 2   |       2 2
    !       |       1 1   |       1 1   |       1 1
    !       |_ _ _z 0 0   |_ _ _z 0 0   |_ _ _y 0 0
    ! 
    !     the arrangement of the coord2 is as follow:
    !       y             x             x
    !       |       0 1   |       0 1   |       0 1
    !       |       0 1   |       0 1   |       0 1
    !       |_ _ _z 0 1   |_ _ _z 0 1   |_ _ _y 0 1
    ! 
    !     neighbor index:
    ! 
    !       y             x             x
    !       |       6 3 5 |       6 3 5 |       6 3 5
    !       |       2 * 1 |       2 * 1 |       2 * 1
    !       |_ _ _z 7 4 8 |_ _ _z 7 4 8 |_ _ _y 7 4 8
    ! 
    !        Here * means the center subdomain, and 1-8 stands for the reduative location of the eight neighbors

    coord1 = int ( nrank / col)
    coord2 = mod ( nrank,  col)
    !print*,'====',nrank,coord1,coord2

    ! Firstly, all the boundaries are assumed to be periodic
    idTop     = mod(coord1+1,    row)  ! top
    idBottom  = mod(coord1+row-1,row)  ! bottom 
    idLeft    = mod(coord2+col-1,col)  ! left
    idRight   = mod(coord2+1,    col)  ! right
    do i=1,3
      myProcNghBC(i,1) = coord1   * col + idRight; ! myProcNghBC(i,5) = idTop    * col + idRight
      myProcNghBC(i,2) = coord1   * col + idLeft ; ! myProcNghBC(i,6) = idTop    * col + idLeft
      myProcNghBC(i,3) = idTop    * col + coord2 ; ! myProcNghBC(i,7) = idBottom * col + idLeft
      myProcNghBC(i,4) = idBottom * col + coord2 ; ! myProcNghBC(i,8) = idBottom * col + idRight
    enddo

    ! Secondly, modify the edge neighbour ids
    IF(coord1==0) THEN
      if(BcOpt(3)<0) myProcNghBC(1,4)=BcOpt(3)
      if(BcOpt(1)<0) myProcNghBC(2,4)=BcOpt(1)
      if(BcOpt(1)<0) myProcNghBC(3,4)=BcOpt(1)
    ENDIF
    IF(coord1==row-1) THEN
      if(BcOpt(4)<0) myProcNghBC(1,3)=BcOpt(4)
      if(BcOpt(2)<0) myProcNghBC(2,3)=BcOpt(2)
      if(BcOpt(2)<0) myProcNghBC(3,3)=BcOpt(2)
    ENDIF
    IF(coord2==0) THEN
      if(BcOpt(5)<0) myProcNghBC(1,2)=BcOpt(5)
      if(BcOpt(5)<0) myProcNghBC(2,2)=BcOpt(5)
      if(BcOpt(3)<0) myProcNghBC(3,2)=BcOpt(3)
    ENDIF
    IF(coord2==col-1) THEN
      if(BcOpt(6)<0) myProcNghBC(1,1)=BcOpt(6)
      if(BcOpt(6)<0) myProcNghBC(2,1)=BcOpt(6)
      if(BcOpt(4)<0) myProcNghBC(3,1)=BcOpt(4)
    ENDIF
  end subroutine Init_ProcNeighbour

  !**********************************************************************
  ! myallocate
  !**********************************************************************
  subroutine myallocate(var, mb, opt_decomp, opt_global)
    implicit none
    type(MatBound),intent(inout)::mb
    logical,intent(in),optional :: opt_global
    type(decomp_info),intent(in),target,optional :: opt_decomp
    real(mytype), allocatable, dimension(:,:,:),intent(out) :: var
 
    ! locals
    logical::global
    integer::ierror
    TYPE(decomp_info),pointer::decomp

    if(mb%pencil<1 .or. mb%pencil>3) then
      call decomp_2d_abort(15,'wrong input pencil when creating new arrays')
    endif

    if (present(opt_decomp)) then
      decomp => opt_decomp
    else
      decomp => decomp_main
    endif
    if (present(opt_global)) then
      global = opt_global
    else
      global = .false.
    endif

    ! first update the MatBound
    select case(mb%pencil)
    case(1) ! x-pencil
      if (global) then
        mb%xmr = decomp%x1st(1); mb%xpr = decomp%x1en(1)
        mb%ymr = decomp%x1st(2); mb%ypr = decomp%x1en(2)
        mb%zmr = decomp%x1st(3); mb%zpr = decomp%x1en(3)
      else
        mb%xmr = 1; mb%xpr = decomp%x1sz(1)
        mb%ymr = 1; mb%ypr = decomp%x1sz(2)
        mb%zmr = 1; mb%zpr = decomp%x1sz(3)
      end if
    case(2)  ! y-pencil
      if (global) then
        mb%xmr = decomp%y1st(1); mb%xpr = decomp%y1en(1)
        mb%ymr = decomp%y1st(2); mb%ypr = decomp%y1en(2)
        mb%zmr = decomp%y1st(3); mb%zpr = decomp%y1en(3)
      else
        mb%xmr = 1; mb%xpr = decomp%y1sz(1)
        mb%ymr = 1; mb%ypr = decomp%y1sz(2)
        mb%zmr = 1; mb%zpr = decomp%y1sz(3)
      end if
    case(3)  ! z-pencil
      if (global) then
        mb%xmr = decomp%z1st(1); mb%xpr = decomp%z1en(1)
        mb%ymr = decomp%z1st(2); mb%ypr = decomp%z1en(2)
        mb%zmr = decomp%z1st(3); mb%zpr = decomp%z1en(3)
      else
        mb%xmr = 1; mb%xpr = decomp%z1sz(1)
        mb%ymr = 1; mb%ypr = decomp%z1sz(2)
        mb%zmr = 1; mb%zpr = decomp%z1sz(3)
      end if
    end select
    mb%xmm= mb%xmr- mb%xme;  mb%xpm= mb%xpr+ mb%xpe
    mb%ymm= mb%ymr- mb%yme;  mb%ypm= mb%ypr+ mb%ype
    mb%zmm= mb%zmr- mb%zme;  mb%zpm= mb%zpr+ mb%zpe 

    allocate(var(mb%xmm:mb%xpm,mb%ymm:mb%ypm,mb%zmm:mb%zpm),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(15,'myallocate-allocate wrong')
  end subroutine myallocate

  !**********************************************************************
  ! myupdate_halo
  !**********************************************************************
  subroutine myupdate_halo(mat,mb,hi)
    implicit none
    type(MatBound),intent(in)::mb
    type(HaloInfo),intent(in)::hi
    real(mytype), dimension(mb%xmm:mb%xpm,mb%ymm:mb%ypm,mb%zmm:mb%zpm),intent(INout):: mat

    ! locals
    integer:: j,data_type,halo_type,ierror,ProcNgh(4) 
    integer:: s1,s2,s3,icount,ilength,ijump,requests(2)
    integer:: xs1,xs2,ys1,ys2,zs1,zs2   ! index for halo sending
    integer:: xr1,xr2,yr1,yr2,zr1,zr2   ! index for halo receiving
    integer,dimension(MPI_STATUS_SIZE,2) :: SRstatus
  
    data_type=real_type
    if(mb%pencil /= hi%pencil ) then
      ierror = 14
      call decomp_2d_abort(ierror,'myupdate_halo wrong: pencil error, Gong Zheng')
    endif
    if(mb%xme<hi%xmh .or. mb%yme<hi%ymh .or. mb%zme<hi%zmh .or.  &
       mb%xpe<hi%xph .or. mb%ype<hi%yph .or. mb%zpe<hi%zph .or.  &
       hi%xmh<0      .or. hi%ymh<0      .or. hi%zmh<0      .or.  &
       hi%xph<0      .or. hi%yph<0      .or. hi%zph<0       ) then
      ierror = 14
      call decomp_2d_abort(ierror,'myupdate_halo wrong: HaloSize error, Gong Zheng')
    endif

    s1= mb%xpm-mb%xmm+1
    s2= mb%ypm-mb%ymm+1
    s3= mb%zpm-mb%zmm+1
    do j=1,4
      if(myProcNghBC(mb%pencil,j)<0) then
        ProcNgh(j)=MPI_PROC_NULL
      else
        ProcNgh(j)=myProcNghBC(mb%pencil,j)
      endif
    enddo

    SELECT CASE(mb%pencil)
    CASE(1)  ! x-pencil
      ! To DO
    CASE(2)  ! y-pencil

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
        if(ProcNgh(3)==nrank) then  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
          mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
        else
          icount = s2*s3
          ilength= hi%xmh
          ijump  = s1
          call MPI_TYPE_VECTOR(icount,ilength,ijump, data_type, halo_type, ierror)
          call MPI_TYPE_COMMIT(halo_type, ierror)
          call MPI_IRECV( mat(xr1,yr1,zr1),1,halo_type,ProcNgh(4),1,MPI_COMM_WORLD,requests(1),ierror)
          call MPI_ISSEND(mat(xs1,ys1,zs1),1,halo_type,ProcNgh(3),1,MPI_COMM_WORLD,requests(2),ierror)
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
        if(ProcNgh(4)==nrank) then  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
          mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
        else
          icount = s2*s3
          ilength= hi%xph
          ijump  = s1
          call MPI_TYPE_VECTOR(icount, ilength, ijump, data_type, halo_type, ierror)
          call MPI_TYPE_COMMIT(halo_type, ierror)
          call MPI_IRECV( mat(xr1,yr1,zr1),1,halo_type,ProcNgh(3),2,MPI_COMM_WORLD,requests(1),ierror)
          call MPI_ISSEND(mat(xs1,ys1,zs1),1,halo_type,ProcNgh(4),2,MPI_COMM_WORLD,requests(2),ierror)
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
        if(ProcNgh(1)==nrank) then  ! neighbour is nrank itself, and ProcNgh(2)==nrank also !!!
          mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
        else
          icount=s1*s2*hi%zmh
          call MPI_IRECV( mat(xr1,yr1,zr1),icount,data_type,ProcNgh(2),3,MPI_COMM_WORLD,requests(1),ierror)
          call MPI_ISSEND(mat(xs1,ys1,zs1),icount,data_type,ProcNgh(1),3,MPI_COMM_WORLD,requests(2),ierror)
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
        if(ProcNgh(2)==nrank) then  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
          mat(xr1:xr2, yr1:yr2, zr1:zr2) = mat(xs1:xs2, ys1:ys2, zs1:zs2)
        else
          icount=s1*s2*hi%zph
          call MPI_IRECV( mat(xr1,yr1,zr1),icount,data_type,ProcNgh(1),4,MPI_COMM_WORLD,requests(1),ierror)
          call MPI_ISSEND(mat(xs1,ys1,zs1),icount,data_type,ProcNgh(2),4,MPI_COMM_WORLD,requests(2),ierror)
          call MPI_WAITALL(2,requests,SRstatus,ierror)
        endif
      ENDIF
    CASE(3)  ! z-pencil
      ! To DO
    END SELECT

  end subroutine myupdate_halo

  !**********************************************************************
  ! best_2d_grid
  !**********************************************************************
  subroutine best_2d_grid(nx,ny,nz,pencil,iproc,best_p_row,best_p_col,BcOptions)
    implicit none
    integer,intent(in)::  nx,ny,nz,pencil,iproc,BcOptions(6)
    integer,intent(out):: best_p_row, best_p_col

    ! locals
    TYPE(MatBound)::mb_test
    TYPE(HaloInfo)::hi_test
    type(decomp_info)::decomp
    real(mytype)::t1,t2,best_time
    integer,allocatable,dimension(:)::factors
    integer::nfact,i,k,row,col,ierror,errorcode
    real(mytype),allocatable, dimension(:,:,:):: Arr1,Arr2

    if (nrank==0) write(*,*) 'Auto-tuning mode......'
    best_time = huge(t1)
    best_p_row= -1
    best_p_col= -1

    i = int(sqrt(real(iproc))) + 10  ! enough space to save all factors 
    allocate(factors(i),stat=ierror)
    if(ierror/=0) call decomp_2d_abort(16,'best_2d_grid-allocate wrong-1')
    call findfactor(iproc, factors, nfact)
    if(nrank==0) write(*,*) 'factors: ', (factors(i), i=1,nfact)

    mb_test%pencil = pencil  
    mb_test%xme=2;  mb_test%xpe=2
    mb_test%yme=2;  mb_test%ype=2
    mb_test%zme=2;  mb_test%zpe=2

    hi_test%pencil = pencil
    hi_test%xmh=2;  hi_test%xph=2
    hi_test%ymh=0;  hi_test%yph=0
    hi_test%zmh=2;  hi_test%zph=2

    do i=1,nfact
      row= factors(i)
      col= iproc / row
      if(min(nx,ny,nz)>=row .and. min(ny,nz,nz)>=col) then
        call decomp_info_init(nx,ny,nz,decomp,[row,col])
        call Init_ProcNeighbour(row,col,BcOptions)

        allocate(Arr1(decomp%x1sz(1),decomp%x1sz(2),decomp%x1sz(3)),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(16,'best_2d_grid-allocate wrong-2')
        allocate(Arr2(decomp%y1sz(1),decomp%y1sz(2),decomp%y1sz(3)),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(16,'best_2d_grid-allocate wrong-3')
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        t1 = MPI_WTIME()
        do k=1,5
          call transpose_x1_to_y1(Arr1,Arr2,decomp)
          call transpose_y1_to_x1(Arr2,Arr1,decomp)
        enddo
        t2 = MPI_WTIME() -t1
        
        deallocate(Arr1)
        allocate(Arr1(decomp%z1sz(1),decomp%z1sz(2),decomp%z1sz(3)),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(16,'best_2d_grid-allocate wrong-4')        
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        t1 = MPI_WTIME()
        do k=1,5
          call transpose_y1_to_z1(Arr2,Arr1,decomp)
          call transpose_z1_to_y1(Arr1,Arr2,decomp)
        enddo
        t2 = MPI_WTIME() -t1 +t2
        
        deallocate(Arr2)     
        allocate(Arr2(decomp%x2sz(1),decomp%x2sz(2),decomp%x2sz(3)),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(16,'best_2d_grid-allocate wrong-5')
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        t1 = MPI_WTIME()
        do k=1,5
          call transpose_z1_to_x2(Arr1,Arr2,decomp)
          call transpose_x2_to_z1(Arr2,Arr1,decomp)
        enddo
        t2 = MPI_WTIME() -t1 +t2           

        deallocate(Arr1)     
        allocate(Arr1(decomp%y2sz(1),decomp%y2sz(2),decomp%y2sz(3)),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(16,'best_2d_grid-allocate wrong-6')
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        t1 = MPI_WTIME()
        do k=1,5
          call transpose_x2_to_y2(Arr2,Arr1,decomp)
          call transpose_y2_to_x2(Arr1,Arr2,decomp)
        enddo
        t2 = MPI_WTIME() -t1 +t2 

        deallocate(Arr2)     
        allocate(Arr2(decomp%z2sz(1),decomp%z2sz(2),decomp%z2sz(3)),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(16,'best_2d_grid-allocate wrong-7')
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        t1 = MPI_WTIME()
        do k=1,5
          call transpose_y2_to_z2(Arr1,Arr2,decomp)
          call transpose_z2_to_y2(Arr2,Arr1,decomp)
        enddo
        t2 = MPI_WTIME() -t1 +t2
        
        deallocate(Arr1)     
        allocate(Arr1(decomp%x1sz(1),decomp%x1sz(2),decomp%x1sz(3)),stat=ierror)
        if(ierror/=0) call decomp_2d_abort(16,'best_2d_grid-allocate wrong-8')
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        t1 = MPI_WTIME()
        do k=1,5
          call transpose_z2_to_x1(Arr2,Arr1,decomp)
          call transpose_x1_to_z2(Arr1,Arr2,decomp)
        enddo
        t2 = MPI_WTIME() -t1 +t2
        
        deallocate(Arr1,Arr2)
        call myallocate(Arr1,mb_test)
        call MPI_BARRIER(MPI_COMM_WORLD,ierror)
        t1 = MPI_WTIME()
        do k=1,40
          call myupdate_halo(Arr1,mb_test,hi_test)
        enddo
        t2 = MPI_WTIME() - t1 + t2
        
        deallocate(Arr1)
        call decomp_info_finalize(decomp)
        call MPI_ALLREDUCE(t2,t1,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)

        if(nrank==0)write(*,*) 'processor grid', row, ' by ', col, ' time=', t1
        if(best_time > t1) then
          best_time = t1
          best_p_row= row
          best_p_col= col
        endif
      endif
    enddo
    deallocate(factors)

    if(best_p_row/=-1) then
      if(nrank==0)write(*,*)'the best processor grid is probably ',best_p_row,' by ',best_p_col
    else
      errorcode = 9
      call decomp_2d_abort(errorcode,'The processor-grid auto-tuning code failed. ' // &
        'The number of processes requested is probably too large.')
    endif
  end subroutine best_2d_grid

#undef DECOMP_TRANSPOSE_X1Y1
#undef DECOMP_TRANSPOSE_Y1Z1
#undef DECOMP_TRANSPOSE_Z1X2
#undef DECOMP_TRANSPOSE_X2Y2
#undef DECOMP_TRANSPOSE_Y2Z2
#undef DECOMP_TRANSPOSE_Z2X1
end module m_decomp2d
