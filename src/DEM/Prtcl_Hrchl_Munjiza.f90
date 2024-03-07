module Prtcl_Hrchl_Munjiza
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Property
  use Prtcl_CL_and_CF
  use Prtcl_Variables
  use Prtcl_decomp_2d
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
#endif
  implicit none
  private

  integer::nlocal,nlocalp,nghost
  real(RK)::xst_cs,yst_cs,zst_cs
    
  type(integer3),dimension(:),allocatable:: box_index       ! integer coordinate of box
  integer,dimension(:), allocatable :: NextX
  integer,dimension(:), allocatable :: NextY
  integer,dimension(:), allocatable :: NextZ  
    
  type::lvl_Munjiza
    integer:: lvl              ! level number
    integer:: lvl_multiple     ! for lvl=1, lvl_multiple=1
    integer:: numPrtcl_lvl = 0 ! number of  particles in this level        
    real(RK):: minD_lvl        ! the minimum diameter of bonnding boxes in this level
    real(RK):: maxD_lvl        ! the maximum diameter of bounding boxes in this level

    real(RK):: cell_len_lvl      ! length of cell
    integer::  nx_lvl      ! number of divisions in x direction
    integer::  ny_lvl      ! number of divisions in y direction
    integer::  nz_lvl      ! number of divisions in z direction
        
    integer,dimension(:), allocatable :: HeadY
    integer,dimension(:), allocatable :: HeadX0
    integer,dimension(:), allocatable :: HeadX
    integer,dimension(:), allocatable :: HeadX2        
    integer:: curr_xList_ind
    integer:: curr_xList2_ind        
        
    integer,dimension(:,:), allocatable :: HeadZ0 ! head list for (iy-1), lower row
    integer,dimension(:,:), allocatable :: HeadZ  ! head list for iy, current row 
    integer,dimension(:,:), allocatable :: HeadZ2
    integer,dimension(0:2):: curr_zList0_ind        
    integer,dimension(0:2):: curr_zList_ind
    integer,dimension(0:2):: curr_zList2_ind
  end type lvl_Munjiza
  type(lvl_Munjiza),allocatable,dimension(:):: lvls_Munjiza ! level

  type::NBS_Munjiza_Hrchl
    integer:: mbox
    integer:: num_lvls = 1            ! number of levels
    integer:: num_Cnsv_cntct = 0      ! number of conservative contacts in the broad search phase
    integer:: lvl_num_cnsv_cntct = 0  ! number of conservative contact in this level        
  contains
    procedure:: Init_Munjiza_Hrchl
        
    ! performing a full contact search 
    procedure:: ContactSearch => NBSMH_ContactSearch
    ! calculating integer coordinates of all boxes
    procedure:: clcBoxIndexAndBuildYList   => NBSMH_clcBoxIndex_and_BuildYList 
            
    procedure:: BuildXList    => NBSMH_BuildXList
    procedure:: BuildXList2   => NBSMH_BuildXList2
    procedure:: BuildZList0   => NBSMH_BuildZList0
    procedure:: BuildZList    => NBSMH_BuildZList
    procedure:: BuildZList2   => NBSMH_BuildZList2
        
    procedure:: LoopNBSMask   => NBSMH_LoopNBSMask
    procedure:: Loop_CrossMask=> NBSMH_Loop_CrossMask
    procedure:: FineSearch    => NBSMH_FineSearch
        
    procedure,private:: OneLevelBroadSearch
    procedure:: Grow_Box_And_Next => NBSMH_Grow_Box_And_Next        
  end type NBS_Munjiza_Hrchl
  type(NBS_Munjiza_Hrchl),public,allocatable :: m_NBS_Munjiza_Hrchl
    
contains
 
  !*********************************************************************
  ! NBS_Munjiza_Hrchl
  !*********************************************************************
  subroutine Init_Munjiza_Hrchl(this )
    implicit none
    class(NBS_Munjiza_Hrchl):: this
        
    integer:: numLevels, idh, nx,ny,nz, dmax_min,numPrtcl_lvl,i, id_level,nLevel_final
    integer:: iErr1, iErr2, iErr3, iErr4, iErr5, iErr6, iErr7,iErrSum
    real(RK):: minD, maxD, minD_lvl, maxD_lvl, cell_len,Diam
    real(RK):: xed_cs,yed_cs,zed_cs
    type(integer3)::numCell
        
    maxD   = 2.0_RK*maxval( DEMProperty%Prtcl_PureProp%Radius )
    minD   = 2.0_RK*minval( DEMProperty%Prtcl_PureProp%Radius )
    dmax_min = int(maxD/minD)
    if( DEM_opt%CS_numlvls <=0) then
      select case( dmax_min )
      case(0:1)
        numLevels = 1
      case(2:3)
        numLevels = 2
      case(4:7)
        numLevels = 3
      case(8:15)
        numLevels = 4
      case(16:31)
        numLevels = 5
      case default
        numLevels = 6
      end select
    else
      numLevels = DEM_opt%CS_numlvls            
    endif

    maxD_lvl = maxD
    nLevel_final = 0
    DO idh = 1, numLevels
      numPrtcl_lvl = 0
      minD_lvl = maxD_lvl/2.0_RK; if(idh==numLevels) minD_lvl=0.0_RK

      do i = 1, DEM_opt%numPrtcl_Type
        Diam = 2.0_RK * DEMProperty%Prtcl_PureProp(i)%Radius
        if(Diam >minD_lvl .and. Diam <= maxD_lvl) then
          numPrtcl_lvl = numPrtcl_lvl + DEMProperty%nPrtcl_in_Bin(i)
        endif
      enddo
      if(numPrtcl_lvl>0) nLevel_final=nLevel_final+1
      maxD_lvl = minD_lvl
    ENDDO
    this%num_lvls = nLevel_final
    this%mbox = GPrtcl_list%mlocal + GPrtcl_list%mGhost_CS
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Contact search method is NBS Munjiza Hierarchy", 3 )
      call DEMLogInfo%OutInfo(" Number of levels is :" // trim( num2str(nLevel_final)),3)
    endif

    allocate( lvls_Munjiza(nLevel_final),  Stat = iErr1 )
    allocate( box_index(this%mbox), Stat = iErr2 ) 
    allocate( NextX( this%mbox ),   STAT = iErr3 )
    allocate( NextY( this%mbox ),   STAT = iErr4 )
    allocate( NextZ( this%mbox ),   STAT = iErr5 ) 
    iErrSum =  abs(iErr1) + abs(iErr2) + abs(iErr3) + abs(iErr4) + abs(iErr5) 
    if(iErrSum /= 0) then
      call DEMLogInfo%CheckForError( ErrT_Abort, "Init_Munjiza_Hrchl", "Allocations failed 1" )
    endif
    NextX = -1
    NextY = -1
    NextZ = -1

    cell_len =  DEM_Opt%Prtcl_cs_ratio * maxD
    xst_cs = DEM_decomp%xSt - cell_len*1.05_RK
    yst_cs = DEM_decomp%ySt - cell_len*1.05_RK
    zst_cs = DEM_decomp%zSt - cell_len*1.05_RK
    xed_cs = DEM_decomp%xEd + cell_len*1.05_RK
    yed_cs = DEM_decomp%yEd + cell_len*1.05_RK 
    zed_cs = DEM_decomp%zEd + cell_len*1.05_RK
    nx = int((xed_cs-xst_cs)/cell_len)+1
    ny = int((yed_cs-yst_cs)/cell_len)+1
    nz = int((zed_cs-zst_cs)/cell_len)+1

    id_level = 0    
    maxD_lvl = maxD ! setting maximum diameter of the first level equal to maximum diameter of bounding boxes
    DO idh = 1, numLevels
      ! the minimum diameter of the level is half of the maximum diameter
      ! modifying the minimum diameter of the last level and sets it to a very small value
      minD_lvl = maxD_lvl/2.0_RK 
      if( idh == numLevels) minD_lvl = 0.0_RK

      numPrtcl_lvl = 0
      do i = 1, DEM_opt%numPrtcl_Type
        Diam = 2.0_RK * DEMProperty%Prtcl_PureProp(i)%Radius
        if(Diam >minD_lvl .and. ((idh==1) .or. (idh>1 .and. Diam <= maxD_lvl))) then
          numPrtcl_lvl = numPrtcl_lvl + DEMProperty%nPrtcl_in_Bin(i)
        endif
      enddo
            
      if(numPrtcl_lvl>0) then
        id_level = id_level + 1
        do i = 1, DEM_opt%numPrtcl_Type
          Diam = 2.0_RK * DEMProperty%Prtcl_PureProp(i)%Radius
          if(Diam >minD_lvl .and. ((idh==1) .or. (idh>1 .and. Diam <= maxD_lvl))) then
            DEMProperty%CS_Hrchl_level(i) = id_level
          endif
        enddo
        lvls_Munjiza(id_level)%lvl = id_level
        lvls_Munjiza(id_level)%lvl_multiple = idh
        lvls_Munjiza(id_level)%nx_lvl = nx
        lvls_Munjiza(id_level)%ny_lvl = ny
        lvls_Munjiza(id_level)%nz_lvl = nz
        lvls_Munjiza(id_level)%cell_len_lvl = cell_len
        lvls_Munjiza(id_level)%numPrtcl_lvl = numPrtcl_lvl

        allocate(lvls_Munjiza(id_level)%HeadY( 0:ny+1), STAT= iErr1)
        allocate(lvls_Munjiza(id_level)%HeadX0(0:nx+1), STAT= iErr2)
        allocate(lvls_Munjiza(id_level)%HeadX( 0:nx+1), STAT= iErr3)
        allocate(lvls_Munjiza(id_level)%HeadX2(0:nx+1), STAT= iErr4)
        allocate(lvls_Munjiza(id_level)%HeadZ0(0:2,0:nz+1), STAT= iErr5)
        allocate(lvls_Munjiza(id_level)%HeadZ( 0:2,0:nz+1), STAT= iErr6)
        allocate(lvls_Munjiza(id_level)%HeadZ2(0:2,0:nz+1), STAT= iErr7)
        iErrSum=abs(iErr1)+abs(iErr2)+abs(iErr3)+abs(iErr4)+abs(iErr5)+ abs(iErr6)+abs(iErr7)
        if(iErrSum/= 0) then
          call DEMLogInfo%CheckForError(ErrT_Abort,"Init_Munjiza_Hrchl","Allocation failed 2")
        endif
                
        lvls_Munjiza(id_level)%HeadY  = -1
        lvls_Munjiza(id_level)%HeadX0 = -1
        lvls_Munjiza(id_level)%HeadX  = -1
        lvls_Munjiza(id_level)%HeadX2 = -1
                
        lvls_Munjiza(id_level)%curr_xList_ind  = -1
        lvls_Munjiza(id_level)%curr_XList2_ind = -1
                
        lvls_Munjiza(id_level)%HeadZ0 = -1
        lvls_Munjiza(id_level)%HeadZ  = -1
        lvls_Munjiza(id_level)%HeadZ2 = -1
        lvls_Munjiza(id_level)%curr_zList0_ind = -1
        lvls_Munjiza(id_level)%curr_zList_ind  = -1
        lvls_Munjiza(id_level)%curr_ZList2_ind = -1

        !>>> log file
        numCell = integer3(nx,ny,nz)
        call DEMLogInfo%OutInfo("Level "//trim( num2str(id_level)), 4)
        call DEMLogInfo%OutInfo("   Cell size is [m]: "// trim(num2str(cell_len)), 4)
        call DEMLogInfo%OutInfo("   number of cells (x,y,z) :"//trim(num2str(numCell)),4,.true.)
        call DEMLogInfo%OutInfo("   number of particles : "//trim(num2str(numPrtcl_lvl)),4,.true.)
      endif
              
      ! halving the maximum diameter of the level to be the maximum diameter  for the next level
      nx= nx*2
      ny= ny*2
      nz= nz*2            
      maxD_lvl = minD_lvl
      cell_len = cell_len/2.0_RK
    ENDDO

  end subroutine Init_Munjiza_Hrchl
    
  !******************************************************************
  ! NBSMH_Grow_Box_And_Next
  !******************************************************************
  subroutine NBSMH_Grow_Box_And_Next(this,nbox)
    implicit none
    class(NBS_Munjiza_Hrchl):: this 
    integer,intent(in)::nbox
    
    ! lcoals
    integer:: sizen
    
    sizen= int(1.2_RK*real(this%mbox,kind=RK))
    sizen= max(sizen, nbox+1)

    deallocate(box_index); allocate(box_index(sizen))
    deallocate(NextX); allocate(NextX(sizen))
    deallocate(NextY); allocate(NextY(sizen))
    deallocate(NextZ); allocate(NextZ(sizen))

    call DEMLogInfo%CheckForError(ErrT_Pass," NBSMH_Grow_Box_And_Next"," Need to reallocate Box_And_Next")
    call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),    3)
    call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(this%mbox)),3)
    call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),    3)

    this%mbox=sizen
  end subroutine NBSMH_Grow_Box_And_Next

  !*********************************************************************
  ! a full contact search for all levels
  !*********************************************************************
  subroutine NBSMH_ContactSearch(this)
    implicit none
    class(NBS_Munjiza_Hrchl) :: this
    integer :: idh
    
    nlocal=GPrtcl_list%nlocal
    if(nlocal == 0) return
    nlocalp = nlocal + 1
    nghost=GPrtcl_list%nGhost_CS
    this%num_Cnsv_cntct = 0
    this%lvl_num_cnsv_cntct = 0

    ! grid index and Ylists of all levels
    call this%clcBoxIndexAndBuildYList()
    
    ! contact search of all levels
    do idh = this%num_lvls, 1, -1
      call this%OneLevelBroadSearch( lvls_Munjiza(idh) )            
    end do
  end subroutine NBSMH_ContactSearch

  !******************************************************************
  ! calculating integer coordinates of all boxes
  !******************************************************************    
  subroutine NBSMH_clcBoxIndex_and_BuildYList(this)
    implicit none
    class(NBS_Munjiza_Hrchl)::this
    integer::i,idh,m,n12, iy
    real(RK):: rpdx,rpdy, rpdz,cell_len

    n12=nlocal + nghost
    if(n12>this%mbox)call this%Grow_Box_And_Next(n12)

    ! nullifying list Y
    do idh=1,this%num_lvls
      lvls_Munjiza(idh)%HeadY = -1
    enddo 
        
    do i= 1,nlocal
      rpdx = GPrtcl_PosR(i)%x - xst_cs
      rpdy = GPrtcl_PosR(i)%y - yst_cs
      rpdz = GPrtcl_PosR(i)%z - zst_cs
      idh= DEMProperty%CS_Hrchl_level(GPrtcl_pType(i))
      cell_len = lvls_Munjiza(idh)%cell_len_lvl
      box_index(i)%x = floor(rpdx/cell_len)+1
      iy = floor(rpdy/cell_len)+1; box_index(i)%y = iy
      box_index(i)%z = floor(rpdz/cell_len)+1

      NextY(i) = lvls_Munjiza(idh)%HeadY(iy)
      lvls_Munjiza(idh)%HeadY(iy) = i
    enddo   
 
    m=1
    do i=nlocal+1,n12
      rpdx = GhostP_PosR(m)%x - xst_cs
      rpdy = GhostP_PosR(m)%y - yst_cs
      rpdz = GhostP_PosR(m)%z - zst_cs
      idh= DEMProperty%CS_Hrchl_level(GhostP_pType(m))
      cell_len = lvls_Munjiza(idh)%cell_len_lvl
      box_index(i)%x = floor(rpdx/cell_len)+1
      iy = floor(rpdy/cell_len)+1; box_index(i)%y = iy
      box_index(i)%z = floor(rpdz/cell_len)+1

      NextY(i) = lvls_Munjiza(idh)%HeadY(iy)
      lvls_Munjiza(idh)%HeadY(iy) = i
      m=m+1    
    enddo

  end subroutine NBSMH_clcBoxIndex_and_BuildYList
 
  !******************************************************************
  ! constructing Xlist of row iy 
  !******************************************************************
  subroutine NBSMH_BuildXList(this,base_lvl,iy)
  implicit none
    class(NBS_Munjiza_Hrchl):: this
    type(lvl_Munjiza):: base_lvl
    integer,intent(in)::iy ! row index
  integer:: n,ix
        
    ! first checking if the xlist of row iy has been constructed previously 
    if(iy==base_lvl%curr_xList_ind) return

    ! nullifying the xlist of current row but keeps the previous raw
    base_lvl%HeadX= -1
    base_lvl%curr_xList_ind = iy
        
    n = base_lvl%HeadY(iy)
    do while (n.ne.-1)
    ix = box_index(n)%x
    NextX( n ) = base_lvl%HeadX(ix)
    base_lvl%HeadX(ix) = n  
    n = NextY(n)
  enddo
  end subroutine NBSMH_BuildXList
    
  !******************************************************************
  ! constructing Xlist of row iy+1
  !******************************************************************
  subroutine NBSMH_BuildXList2(this,base_lvl,iy)
  implicit none
    class(NBS_Munjiza_Hrchl):: this
    type(lvl_Munjiza):: base_lvl
  integer,intent(in)::iy ! row index
  integer:: n,ix
        
  ! first checking if the xlist of row iy has been constructed previously 
    if(iy == base_lvl%curr_xList2_ind) return

    ! nullifying the xlist of current row but keeps the previous raw
  base_lvl%HeadX2= -1
    base_lvl%curr_xList_ind = iy
        
  n = base_lvl%HeadY(iy)
  do while(n.ne.-1)
    ix = box_index(n)%x
    NextX(n) = base_lvl%HeadX2(ix)
    base_lvl%HeadX2(ix) = n  
    n=NextY(n)
  enddo
  end subroutine NBSMH_BuildXList2    

  !*********************************************************************
  !   Constructing the ZList0 of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine NBSMH_BuildZList0(this,base_lvl, ix,m,lcheck)
    implicit none
    class(NBS_Munjiza_Hrchl ):: this
    type(lvl_Munjiza):: base_lvl
  integer,intent(in):: ix ! column index
    integer,intent(in):: m  ! the column location with respect to ix  
    logical,optional,intent(in):: lcheck
    integer:: n,iz
  
  if(present(lcheck)) then
      if(lcheck.and.ix == base_lvl%curr_zList0_ind(m)) return
    endif
        
    ! 0: previous (left) column, 1: current column, 2: next (right) column
  base_lvl%HeadZ0(m,:) = -1
  base_lvl%curr_ZList0_ind(m) = ix
  n = base_lvl%HeadX0(ix+m-1) ! reading from the row below iy (or iy-1)
  do while (n.ne.-1)
      iz = box_index(n)%z
    NextZ(n) = base_lvl%HeadZ0(m,iz)
      base_lvl%HeadZ0(m,iz) = n
      n = NextX(n)
  enddo
  end subroutine NBSMH_BuildZList0    

  !*********************************************************************
  !   Constructing the ZList of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine NBSMH_BuildZList(this, base_lvl,ix,m,lcheck)
    implicit none
    class(NBS_Munjiza_Hrchl ):: this
    type(lvl_Munjiza):: base_lvl
  integer,intent(in):: ix ! col index
    integer,intent(in):: m  ! the column location with resect to ix 
    logical,optional,intent(in):: lcheck
    integer::n,iz
      
    if(present(lcheck)) then
      if(lcheck.and.ix == base_lvl%curr_zList_ind(m)) return
    endif
        
  ! nullifying the zlist of current col ix, but keeps the previous and next cols
    ! 0: previous (left) column, 1: current column, 2: next (right) column
  base_lvl%HeadZ(m,:) = -1
  base_lvl%curr_ZList_ind(m) = ix
  n = base_lvl%HeadX(ix+m-1) ! reading from current row iy
        
    do while (n.ne.-1)
      iz = box_index(n)%z
    NextZ(n) = base_lvl%HeadZ(m,iz)
      base_lvl%HeadZ(m,iz) = n
      n = NextX(n)
  enddo
  end subroutine NBSMH_BuildZList

  !*********************************************************************
  !   Constructing the ZList2 of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine NBSMH_BuildZList2(this,base_lvl,ix, m,lcheck )
    implicit none
    class(NBS_Munjiza_Hrchl ):: this
    type(lvl_Munjiza):: base_lvl
  integer,intent(in) :: ix ! col index
    integer,intent(in) :: m  ! the column location with resect to ix 
    logical,optional,intent(in):: lcheck
    integer::n,iz
      
    if(present(lcheck)) then
      if(lcheck.and.ix == base_lvl%curr_zList2_ind(m)) return
    endif
        
  ! nullifying the zlist of current col ix, but keeps the previous and next cols
    ! 0: previous (left) column, 1: current column, 2: next (right) column
  base_lvl%HeadZ2(m,:) = -1
  base_lvl%curr_ZList2_ind(m) = ix
  n = base_lvl%HeadX2(ix+m-1) ! reading from current row iy
        
    do while (n.ne.-1)
      iz = box_index(n)%z
    NextZ(n) = base_lvl%HeadZ2(m,iz)
      base_lvl%HeadZ2(m,iz) = n
      n = NextX(n)
  enddo
  end subroutine NBSMH_BuildZList2
    
  !******************************************************************************
  !One level contact search (intra-level and cross level with next levels)
  !******************************************************************************
  subroutine OneLevelBroadSearch(this,base_lvl)
    implicit none
    class(NBS_Munjiza_Hrchl):: this
    type(lvl_Munjiza) base_lvl
    integer:: idh,lvlCoe
    integer:: ix,iy,iz,crs_indx,crs_indy,crs_indz
    
    ! same level
    base_lvl%HeadX0(:) = -1
    ! starting the loop over all rows in base level
    DO iy=1,base_lvl%ny_lvl
      call this%BuildXList(base_lvl, iy)
      ! constructing the xLists of lvl at above rows
      DO idh=base_lvl%lvl-1, 1, -1
        lvlCoe= 2**(base_lvl%lvl_multiple-lvls_Munjiza(idh)%lvl_multiple)
        crs_indy=(iy -1)/lvlCoe + 1
        if(crs_indy == 1) then
          ! for the first row, the xLists of below and top rows should be constructed
          lvls_Munjiza(idh)%HeadX0(:) = -1
          call this%BuildXlist(lvls_Munjiza(idh), 1)
        endif
        call this%BuildXList2(lvls_Munjiza(idh),  crs_indy+1 )
      ENDDO
                        
      ! if row is non-empty
      IF(base_lvl%HeadY(iy) .ne. -1 ) THEN
        base_lvl%HeadZ(0,:) = -1
        base_lvl%HeadZ0(0,:)= -1
               
        ! Creating the zlist of the lower row and current ix (column)
        call this%BuildZList0(base_lvl, 1, 1)
        do ix = 1,base_lvl%nx_lvl
          call this%BuildZList(base_lvl, ix, 1)
          call this%BuildZList0(base_lvl,ix, 2)
                
          ! constructing the zLists of the lvl at the right column
          do idh =base_lvl%lvl-1, 1, -1
                
            lvlCoe= 2**(base_lvl%lvl_multiple-lvls_Munjiza(idh)%lvl_multiple)
            crs_indx = (ix-1)/lvlCoe + 1
            if(crs_indx == 1 ) then
              ! for the first column, the zLists of current and left columns should be constructed
              lvls_Munjiza(idh)%HeadZ0(0,:) = -1
              lvls_Munjiza(idh)%HeadZ(0,:)  = -1
              lvls_Munjiza(idh)%HeadZ2(0,:) = -1
              call this%BuildZList0(lvls_Munjiza(idh),crs_indx,1, .true. )    
              call this%BuildZList(lvls_Munjiza(idh), crs_indx,1, .true. )    
              call this%BuildZList2(lvls_Munjiza(idh),crs_indx,1, .true. )
            endif
                        
            call this%BuildZList0(lvls_Munjiza(idh),crs_indx,2, .true. )    
            call this%BuildZList(lvls_Munjiza(idh), crs_indx,2, .true. )    
            call this%BuildZList2(lvls_Munjiza(idh),crs_indx,2, .true. )  
          enddo                        
                              
          if(base_lvl%Headx(ix).ne.-1) then
                    
            do iz= 1, base_lvl%nz_lvl
              ! same level NBS mask check
              call this%LoopNBSMask(base_lvl, iz)
              do idh = base_lvl%lvl-1, 1, -1
                lvlCoe= 2**(base_lvl%lvl_multiple-lvls_Munjiza(idh)%lvl_multiple)
                crs_indz = (iz-1)/lvlCoe + 1
                call this%Loop_CrossMask(base_lvl,iz,crs_indz,lvls_Munjiza(idh))
              enddo
            enddo
          endif
                
          ! same row, subs
          base_lvl%HeadZ(0,:) = base_lvl%HeadZ(1,:)
            
          ! lower row, subs
          base_lvl%HeadZ0(0,:) = base_lvl%HeadZ0(1,:)
          base_lvl%HeadZ0(1,:) = base_lvl%HeadZ0(2,:)
                
          do idh = base_lvl%lvl-1, 1, -1
            lvlCoe= 2**(base_lvl%lvl_multiple-lvls_Munjiza(idh)%lvl_multiple)
            if(mod(ix, lvlCoe) == 0)then
              ! swap zlists
                        
              lvls_Munjiza(idh)%HeadZ0(0,:) = lvls_Munjiza(idh)%HeadZ0(1,:)
              lvls_Munjiza(idh)%HeadZ0(1,:) = lvls_Munjiza(idh)%HeadZ0(2,:)
              lvls_Munjiza(idh)%curr_ZList0_ind(0) = lvls_Munjiza(idh)%curr_ZList0_ind(1)
              lvls_Munjiza(idh)%curr_ZList0_ind(1) = lvls_Munjiza(idh)%curr_ZList0_ind(2)
    
              lvls_Munjiza(idh)%HeadZ(0,:) = lvls_Munjiza(idh)%HeadZ(1,:)
              lvls_Munjiza(idh)%HeadZ(1,:) = lvls_Munjiza(idh)%HeadZ(2,:)
              lvls_Munjiza(idh)%curr_ZList_ind(0) = lvls_Munjiza(idh)%curr_ZList_ind(1)
              lvls_Munjiza(idh)%curr_ZList_ind(1) = lvls_Munjiza(idh)%curr_ZList_ind(2)
    
              lvls_Munjiza(idh)%HeadZ2(0,:) = lvls_Munjiza(idh)%HeadZ2(1,:)
              lvls_Munjiza(idh)%HeadZ2(1,:) = lvls_Munjiza(idh)%HeadZ2(2,:)
              lvls_Munjiza(idh)%curr_ZList2_ind(0) = lvls_Munjiza(idh)%curr_ZList2_ind(1)
              lvls_Munjiza(idh)%curr_ZList2_ind(1) = lvls_Munjiza(idh)%curr_ZList2_ind(2)
                        
            endif 
          enddo
        enddo
      ENDIF
        
      base_lvl%Headx0(:) = base_lvl%Headx(:)
      DO idh= base_lvl%lvl-1, 1, -1
        lvlCoe= 2**(base_lvl%lvl_multiple-lvls_Munjiza(idh)%lvl_multiple)
        if(mod(iy, lvlCoe) == 0)then
          ! swap xlists
                
          lvls_Munjiza(idh)%HeadX0 = lvls_Munjiza(idh)%HeadX
          lvls_Munjiza(idh)%HeadX = lvls_Munjiza(idh)%HeadX2
          lvls_Munjiza(idh)%curr_xList_ind = lvls_Munjiza(idh)%curr_xList2_ind
        endif
      ENDDO
    ENDDO
    
  end subroutine OneLevelBroadSearch

  !**********************************************************************
  ! finding contacts between particles in the target cell and particles in
  ! cells determined by NBS mask.  
  !**********************************************************************
  subroutine NBSMH_LoopNBSMask( this, base_lvl, iz)
    implicit none
    class(NBS_Munjiza_Hrchl) :: this
    type(lvl_Munjiza):: base_lvl
    integer,intent(in)  :: iz
    integer m, n, i, lx
  
    m = base_lvl%HeadZ(1,iz)
    DO WHILE( m .ne. -1 )

      IF(m<nlocalp) THEN

        !over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          call this%FineSearch(n,m)
          this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = base_lvl%HeadZ(1,iz-1)
        DO WHILE (n .ne. -1)
          call this%FineSearch(n,m)
          this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = base_lvl%HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            call this%FineSearch(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = base_lvl%HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              call this%FineSearch(n,m)
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              n = NextZ(n)
            ENDDO
          enddo
        enddo
        m= NextZ(m)
      ELSE
        !over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          if(n<nlocalp) then
            call this%FineSearch(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = base_lvl%HeadZ(1,iz-1)
        DO WHILE (n .ne. -1)
          if(n<nlocalp) then
            call this%FineSearch(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = base_lvl%HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            if(n<nlocalp) then
              call this%FineSearch(n,m)
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            endif
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = base_lvl%HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              if(n<nlocalp) then
                call this%FineSearch(n,m)
                this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              endif
              n = NextZ(n)
            ENDDO
          enddo
        enddo
        m= NextZ(m)
      ENDIF
    ENDDO
  end subroutine  NBSMH_LoopNBSMask

  !********************************************************************** 
  ! particle  fine search
  !**********************************************************************    
  subroutine NBSMH_FineSearch(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza_Hrchl):: this
    integer,intent(in)  :: pid1, pid2
    integer::gid
    real(RK):: dx,dy,dz,dr,d2sum,dr2,ovrlp

    IF(pid1>nlocal) THEN
      gid=pid1-nlocal
      dr= GhostP_PosR(gid)%w + GPrtcl_PosR(pid2)%w
      dr2= dr*dr
      dx= GhostP_PosR(gid)%x - GPrtcl_PosR(pid2)%x
      d2sum=dx*dx;             if(d2sum>dr2) return
      dy= GhostP_PosR(gid)%y - GPrtcl_PosR(pid2)%y
      d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
      dz= GhostP_PosR(gid)%z - GPrtcl_PosR(pid2)%z
      d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
      ovrlp = dr-sqrt(d2sum)  
      call GPPW_CntctList%AddContactPPG(pid2,gid,ovrlp)
    ELSEIF(pid2>nlocal) THEN
      gid=pid2-nlocal
      dr= GPrtcl_PosR(pid1)%w + GhostP_PosR(gid)%w
      dr2= dr*dr
      dx= GPrtcl_PosR(pid1)%x - GhostP_PosR(gid)%x
      d2sum=dx*dx;             if(d2sum>dr2) return
      dy= GPrtcl_PosR(pid1)%y - GhostP_PosR(gid)%y
      d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
      dz= GPrtcl_PosR(pid1)%z - GhostP_PosR(gid)%z
      d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
      ovrlp = dr-sqrt(d2sum)  
      call GPPW_CntctList%AddContactPPG(pid1,gid,ovrlp)
    ELSE
      dr= GPrtcl_PosR(pid1)%w + GPrtcl_PosR(pid2)%w
      dr2= dr*dr
      dx= GPrtcl_PosR(pid1)%x - GPrtcl_PosR(pid2)%x
      d2sum=dx*dx;             if(d2sum>dr2) return
      dy= GPrtcl_PosR(pid1)%y - GPrtcl_PosR(pid2)%y
      d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
      dz= GPrtcl_PosR(pid1)%z - GPrtcl_PosR(pid2)%z
      d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
      ovrlp = dr-sqrt(d2sum)        

      ! this is a convention, the lower id should be the first item in the contact pair (particle & particle)
      if(GPrtcl_id(pid1) < GPrtcl_id(pid2) ) then
         call GPPW_CntctList%AddContactPP(pid1,pid2,ovrlp)
      else
         call GPPW_CntctList%AddContactPP(pid2,pid1,ovrlp)
      endif
    ENDIF
        
  end subroutine NBSMH_FineSearch 
    
  !******************************************************************************
  !   Performing a cross-level contact search between the current level (this)
  ! and CrossMunjiza level. A cross level contact search performed between 
  ! particles from target cell in current level and 27 neighbor cells from
  ! CrossMunjiza level
  !******************************************************************************
  subroutine NBSMH_Loop_CrossMask(this, base_lvl, iz, crs_iz, CrossMunjiza)
    implicit none
    class(NBS_Munjiza_Hrchl):: this
    type(lvl_Munjiza)::base_lvl
    integer, intent(in):: iz,crs_iz
    type(lvl_Munjiza),intent(in):: CrossMunjiza
    integer::m,n,lx,l
    
    m = base_lvl%HeadZ(1,iz)
    DO WHILE(m.ne.-1)
      IF(m<nlocalp) THEN
        ! first, looping all 9 cells of cross level located in the row below of the target cell
        do lx = 0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ0(lx,crs_iz+l)
            do while( n .ne. -1)
              ! performing the narrow phase search
              call this%FineSearch(n,m)
              this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              n = NextZ(n)
            enddo
          enddo
        enddo
         
        ! second, looping all 9 cells of cross level located in the same row as the target cell
        do lx=0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ(lx,crs_iz+l)
            do while(n.ne.-1)
              ! performing the narrow phase search
              call this%FineSearch(n,m)
              this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              n = NextZ(n)
            enddo
          enddo
        enddo
          
        ! third, looping all 9 cells of cross level located in the row above the target cell
        do lx = 0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ2(lx,crs_iz+l)
            do while(n.ne.-1)
              ! performing the narrow phase search
              call this%FineSearch(n,m)
              this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              n = NextZ(n)
            enddo
          enddo
        enddo
        m = NextZ(m)

      ELSE
        ! first, looping all 9 cells of cross level located in the row below of the target cell
        do lx = 0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ0(lx,crs_iz+l)
            do while( n .ne. -1)
              ! performing the narrow phase search
              if(n<nlocalp .or. m<nlocalp) then
                call this%FineSearch(n,m)
                this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              endif
              n = NextZ(n)
            enddo
          enddo
        enddo
         
        ! second, looping all 9 cells of cross level located in the same row as the target cell
        do lx=0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ(lx,crs_iz+l)
            do while(n.ne.-1)
              ! performing the narrow phase search
              if(n<nlocalp .or. m<nlocalp) then
                call this%FineSearch(n,m)
                this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              endif
              n = NextZ(n)
            enddo
          enddo
        enddo
          
        ! third, looping all 9 cells of cross level located in the row above the target cell
        do lx = 0,2
          do l=-1,1
            n = CrossMunjiza%HeadZ2(lx,crs_iz+l)
            do while(n.ne.-1)
              ! performing the narrow phase search
              if(n<nlocalp .or. m<nlocalp) then
                call this%FineSearch(n,m)
                this%lvl_num_cnsv_cntct = this%lvl_num_cnsv_cntct + 1
              endif
              n = NextZ(n)
            enddo
          enddo
        enddo
        m = NextZ(m)

      ENDIF
    ENDDO
  end subroutine NBSMH_Loop_CrossMask

end module Prtcl_Hrchl_Munjiza
