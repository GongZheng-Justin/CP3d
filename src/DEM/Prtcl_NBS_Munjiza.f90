module Prtcl_NBS_Munjiza
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_CL_and_CF
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_Parameters
  use Prtcl_decomp_2d
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
#endif
  implicit none
  private
  
  real(RK)::xst_cs,yst_cs,zst_cs
  integer::nFixed,nFixedP,nlocal,nFixedAndlocal,nFixedAndlocalP,nNeedCSTotoal
  
  type(integer3),dimension(:),allocatable:: box_index ! integer coordinate of box
  integer,dimension(:),allocatable :: NextX
  integer,dimension(:),allocatable :: NextY
  integer,dimension(:),allocatable :: NextZ
  integer,dimension(:),allocatable :: HeadY
  integer,dimension(:),allocatable :: HeadX
  integer,dimension(:),allocatable :: HeadX0
  integer,dimension(:,:),allocatable :: HeadZ  ! head list for iy, current row 
  integer,dimension(:,:),allocatable :: HeadZ0 ! head list for (iy-1), lower row
    
  type::NBS_Munjiza
    integer:: mbox
    real(RK)::maxDiam
    real(RK)::cell_len   ! length of cell
    integer:: nx         ! number of divisions in x direction
    integer:: ny         ! number of divisions in y direction
    integer:: nz         ! number of divisions in z direction       
    integer:: num_Cnsv_cntct = 0 !number of conservative contacts in the broad search phase
  contains
    procedure:: Init_NBSM
        
    ! performing contact search (includes all steps)
    procedure:: ContactSearch=> NBSM_ContactSearch
    procedure:: clcBoxIndex  => NBSM_clcBoxIndex ! calculating integer coordinates of all boxes
    procedure:: BuildYList   => NBSM_BuildYList  ! constructing YList
    procedure:: BuildXList   => NBSM_BuildXList  ! constructing XList
    procedure:: BuildZList   => NBSM_BuildZList  ! constructing ZList
    procedure:: BuildZList0  => NBSM_BuildZList0
    procedure:: FineSearch1  => NBSM_FineSearch1
    procedure:: FineSearch2  => NBSM_FineSearch2
    procedure:: FineSearch3  => NBSM_FineSearch3
    procedure:: LoopNBSMask  => NBSM_LoopNBSMask
  end type NBS_Munjiza
  type(NBS_Munjiza),public,allocatable:: m_NBS_Munjiza
    
contains
  !******************************************************************
  ! Initializing NBS_Munjiza object
  !******************************************************************
  subroutine Init_NBSM(this)
    implicit none
    class(NBS_Munjiza)::this

    ! locals
    type(integer3)::numCell
    real(RK)::xed_cs,yed_cs,zed_cs
    integer::i,iErr1,iErr2,iErr3,iErr4,iErr5,iErr6,iErr7,iErr8,iErr9,iErrSum

    this%maxDiam = two*maxval( DEMProperty%Prtcl_PureProp%Radius )
#ifndef CFDACM
    this%cell_len= DEM_Opt%Prtcl_cs_ratio*this%maxDiam
#else
    this%cell_len= DEM_Opt%Prtcl_cs_ratio*this%maxDiam +maxval(dlub_pp)
#endif

    xst_cs = DEM_decomp%xSt - this%cell_len*1.05_RK
    yst_cs = DEM_decomp%ySt - this%cell_len*1.05_RK
    zst_cs = DEM_decomp%zSt - this%cell_len*1.05_RK
    xed_cs = DEM_decomp%xEd + this%cell_len*1.05_RK
    yed_cs = DEM_decomp%yEd + this%cell_len*1.05_RK 
    zed_cs = DEM_decomp%zEd + this%cell_len*1.05_RK
    this%nx = int((xed_cs-xst_cs)/this%cell_len)+1
    this%ny = int((yed_cs-yst_cs)/this%cell_len)+1
    this%nz = int((zed_cs-zst_cs)/this%cell_len)+1
        
    numcell = integer3(this%nx, this%ny, this%nz)
    if(nrank==0) then
      call DEMLogInfo%OutInfo("Contact search method is NBS Munjiza", 3 )
      call DEMLogInfo%OutInfo("Cell size is [m]: "// trim(num2str(this%cell_len)), 4)
      call DEMLogInfo%OutInfo("Number of cells considered is (x,y,z) :"//trim(num2str(numCell)),4) 
    endif       
        
    nFixed= GPrtcl_list%mlocalFix+ GPrtcl_list%nGhostFix_CS; nFixedP=nFixed+1
    this%mbox = nFixed+ GPrtcl_list%mlocal + GPrtcl_list%mGhost_CS

    allocate(box_index(this%mbox),   STAT=iErr1) 
    allocate(HeadY(this%ny),         STAT=iErr2)
    allocate(HeadX(this%nx),         STAT=iErr3)
    allocate(HeadX0(this%nx),        STAT=iErr4)
    allocate(HeadZ(0:1,0:this%nz+1), STAT=iErr5)
    allocate(HeadZ0(0:2,0:this%nz+1),STAT=iErr6)
    allocate(NextY(this%mbox),       STAT=iErr7)
    allocate(NextX(this%mbox),       STAT=iErr8)
    allocate(NextZ(this%mbox),       STAT=iErr9)

    iErrSum=abs(iErr1)+abs(iErr2)+abs(iErr3)+abs(iErr4)+abs(iErr5)+abs(iErr6)+abs(iErr7)+abs(iErr8)+abs(iErr9)
    if(iErrSum/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"Init_NBSM: ","Allocation failed " )
        
    HeadY = -1; NextY = -1
    HeadX = -1; NextX = -1; HeadX0 = -1
    HeadZ = -1; NextZ = -1; HeadZ0 = -1

    ! Firstly, clculate the Box_index for fixed particles(including the relevant ghost fixed particles)
    do i=1,nFixed
      box_index(i)%x = floor(( GPFix_PosR(i)%x - xst_cs )/this%cell_len)+1
      box_index(i)%y = floor(( GPFix_PosR(i)%y - yst_cs )/this%cell_len)+1
      box_index(i)%z = floor(( GPFix_PosR(i)%z - zst_cs )/this%cell_len)+1
    enddo
  end subroutine Init_NBSM

  !******************************************************************
  ! performing a contact search on all particles (includes all steps)
  !******************************************************************
  subroutine NBSM_ContactSearch(this)    
    implicit none
    class(NBS_Munjiza):: this      

    ! locals
    integer::ix,iy,iz

    nlocal=GPrtcl_list%nlocal
    this%num_Cnsv_cntct = 0
    if(nlocal == 0) return
    nFixedAndlocal  = nFixed + nlocal
    nFixedAndlocalP = nFixedAndlocal + 1
    nNeedCSTotoal = nFixedAndlocal + GPrtcl_list%nGhost_CS

    call this%clcBoxIndex()
    call this%BuildYList()

    HeadX0 = -1
    DO iy = 1, this%ny
      call this%BuildXList(iy)
      if(HeadY(iy).ne. -1 ) then
        HeadZ(0,:)  = -1
        HeadZ0(0,:) = -1
        call this%BuildZlist0(1,1)
        do ix = 1, this%nx
          call this%BuildZlist(ix,1)
          call this%BuildZlist0(ix,2)

          if(HeadX(ix).ne.-1) then
            do iz=1,this%nz
              call this%LoopNBSMask(iz)
            enddo
          endif
          HeadZ(0,:) = HeadZ(1,:)  ! same row, subs
          HeadZ0(0,:)= HeadZ0(1,:) ! lower row, subs
          HeadZ0(1,:)= HeadZ0(2,:)
        enddo
      endif
      HeadX0 = HeadX
    ENDDO
  end subroutine NBSM_ContactSearch

  !******************************************************************
  ! calculating integer coordinates of all boxes
  !******************************************************************
  subroutine NBSM_clcBoxIndex(this)
    implicit none
    class(NBS_Munjiza):: this

    ! locals
    real(RK)::rpdx,rpdy,rpdz
    integer::i,m,sizen,ierrTmp,ierror=0
    type(integer3),allocatable,dimension(:)::Int3Vec
  
    if(nNeedCSTotoal>this%mbox) then
      sizen= int(1.2_RK*real(this%mbox,kind=RK))
      sizen= max(sizen, nNeedCSTotoal+1)

      call move_alloc(box_index,Int3Vec)
      allocate(box_index(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      if(nFixed>0)box_index(1:nFixed)= Int3Vec(1:nFixed)
      deallocate(Int3Vec)

      deallocate(NextX); allocate(NextX(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      deallocate(NextY); allocate(NextY(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      deallocate(NextZ); allocate(NextZ(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      if(ierror/=0) then
        call DEMLogInfo%CheckForError(ErrT_Abort," NBSM_clcBoxIndex"," Reallocate wrong!")
        call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
      endif      
      !call DEMLogInfo%CheckForError(ErrT_Pass," NBSM_clcBoxIndex"," Need to reallocate Box_And_Next")
      !call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),    3)
      !call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(this%mbox)),3)
      !call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),    3)
      this%mbox=sizen
    endif

    m=1
    do i=nFixed+1, nFixedAndlocal
      rpdx = GPrtcl_PosR(m)%x - xst_cs
      rpdy = GPrtcl_PosR(m)%y - yst_cs
      rpdz = GPrtcl_PosR(m)%z - zst_cs
      box_index(i)%x = floor(rpdx/this%cell_len)+1
      box_index(i)%y = floor(rpdy/this%cell_len)+1
      box_index(i)%z = floor(rpdz/this%cell_len)+1
      m=m+1
    enddo

    m=1
    do i=nFixedAndlocalP, nNeedCSTotoal
      rpdx = GhostP_PosR(m)%x - xst_cs
      rpdy = GhostP_PosR(m)%y - yst_cs
      rpdz = GhostP_PosR(m)%z - zst_cs
      box_index(i)%x = floor(rpdx/this%cell_len)+1
      box_index(i)%y = floor(rpdy/this%cell_len)+1
      box_index(i)%z = floor(rpdz/this%cell_len)+1
      m=m+1    
    enddo
  end subroutine NBSM_clcBoxIndex
    
  !******************************************************************
  ! constructing Ylist of particles
  !******************************************************************
  subroutine NBSM_BuildYList(this)
    implicit none
    class(NBS_Munjiza) this
    integer:: i,iy

    HeadY = -1  ! nullifying list Y
    do i=1,nNeedCSTotoal  
      iy=box_index(i)%y
      NextY(i) = HeadY(iy) 
      HeadY(iy)= i
    enddo
  end subroutine NBSM_BuildYList

  !******************************************************************
  ! constructing Xlist of row iy 
  !******************************************************************
  subroutine NBSM_BuildXList(this,iy)
    implicit none
    class(NBS_Munjiza)::this
    integer,intent(in)::iy ! row index
    integer::n,ix

    ! nullifying the xlist of current row but keeps the previous raw
    HeadX= -1
        
    n = HeadY(iy)
    do while (n .ne. -1)
      ix = box_index(n)%x
      NextX(n) = HeadX(ix)
      HeadX(ix)= n  
      n = NextY(n)
    enddo
  end subroutine NBSM_BuildXList

  !*********************************************************************
  !   Constructing the ZList of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine NBSM_BuildZList(this,ix,m)
    implicit none
    class(NBS_Munjiza ) this
    integer,intent(in) :: ix ! col index
    integer,intent(in) :: m  ! the column location with resect to ix
    integer:: n,iz
 
    ! nullifying the zlist of current col ix, but keeps the previous and next cols
    ! 0: previous (left) column, 1: current column, 2: next (right) column
    HeadZ(m,:) = -1
    n = HeadX(ix+m-1) ! reading from current row iy
    do while ( n .ne. -1 )
      iz = box_index(n)%z
      NextZ(n) = HeadZ(m,iz)
      HeadZ(m,iz) = n
      n = NextX(n)
    enddo
  end subroutine NBSM_BuildZList
    
  !*********************************************************************
  !   Constructing the ZList0 of column ix, ix-1, or ix+1 depending on the value of m
  !*********************************************************************
  subroutine NBSM_BuildZList0(this,ix,m)
    implicit none
    class(NBS_Munjiza ):: this
    integer,intent(in):: ix ! column index
    integer,intent(in):: m  ! the column location with respect to ix
    integer:: n,iz,ixm
        
    ! 0: previous (left) column, 1: current column, 2: next (right) column
    HeadZ0(m,:) = -1
    ixm=ix+m-1
    if(ixm>this%nx) return
   
    n = HeadX0(ixm) ! reading from the row below iy (or iy-1)
    do while ( n .ne. -1 )
      iz = box_index(n)%z
      NextZ(n) = HeadZ0(m,iz)
      HeadZ0(m,iz) = n
      n = NextX(n)
    enddo
  end subroutine NBSM_BuildZList0

  !**********************************************************************
  ! finding contacts between particles in the target cell and particles in
  ! cells determined by NBS mask.  
  !**********************************************************************
  subroutine NBSM_LoopNBSMask(this, iz)
    implicit none
    class(NBS_Munjiza) this
    integer,intent(in)  :: iz
    integer m, n, i, lx
  
    m = HeadZ(1,iz)
    DO WHILE(m.ne.-1)
      IF(m<nFixedP) THEN               !==================================

        ! over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          if(n<nFixedAndlocalP .and. n>nFixed ) then
            call this%FineSearch1(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = HeadZ(1,iz-1)
        DO WHILE (n.ne.-1)
          if(n<nFixedAndlocalP .and. n>nFixed ) then
            call this%FineSearch1(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            if(n<nFixedAndlocalP .and. n>nFixed ) then
              call this%FineSearch1(n,m)
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            endif
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              if(n<nFixedAndlocalP .and. n>nFixed ) then
                call this%FineSearch1(n,m)
                this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              endif
              n = NextZ(n)
            ENDDO
          enddo
        enddo

      ELSEIF(m>nFixedAndlocal) THEN    !==================================

        !over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          if(n<nFixedAndlocalP .and. n>nFixed ) then
            call this%FineSearch3(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = HeadZ(1,iz-1)
        DO WHILE (n.ne.-1)
          if(n<nFixedAndlocalP .and. n>nFixed ) then
            call this%FineSearch3(n,m)
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          endif
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            if(n<nFixedAndlocalP .and. n>nFixed ) then
              call this%FineSearch3(n,m)
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            endif
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              if(n<nFixedAndlocalP .and. n>nFixed ) then
                call this%FineSearch3(n,m)
                this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              endif
              n = NextZ(n)
            ENDDO
          enddo
        enddo

      ELSE                             !==================================

       !over particles in the same cell but not the same particle (to prevent self-contact)
        n = NextZ(m)
        DO WHILE(n.ne.-1)
          if(n<nFixedP) then
            call this%FineSearch1(m,n)
          elseif(n>nFixedAndlocal) then
            call this%FineSearch3(m,n)
          else
            call this%FineSearch2(m,n)
          endif
          this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          n = NextZ(n)
        ENDDO

        ! over particles in (ix, iy , iz-1)
        n = HeadZ(1,iz-1)
        DO WHILE (n.ne.-1)
          if(n<nFixedP) then
            call this%FineSearch1(m,n)
          elseif(n>nFixedAndlocal) then
            call this%FineSearch3(m,n)
          else
            call this%FineSearch2(m,n)
          endif
          this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
          n = NextZ(n)
        ENDDO

        ! over particles in all cells located at (ix-1) and (iy)
        do i = -1,1
          n = HeadZ(0,iz+i)
          DO WHILE(n.ne.-1)
            if(n<nFixedP) then
              call this%FineSearch1(m,n)
            elseif(n>nFixedAndlocal) then
              call this%FineSearch3(m,n)
            else
              call this%FineSearch2(m,n)
            endif
            this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
            n = NextZ(n)
          ENDDO
        enddo
                        
        ! over particles in all 9 cells located at row (iy-1)
        do lx = 0,2
          do i=-1,1
            n = HeadZ0(lx,iz+i)
            DO WHILE (n.ne.-1)
              if(n<nFixedP) then
                call this%FineSearch1(m,n)
              elseif(n>nFixedAndlocal) then
                call this%FineSearch3(m,n)
              else
                call this%FineSearch2(m,n)
              endif
              this%num_Cnsv_cntct = this%num_Cnsv_cntct + 1
              n = NextZ(n)
            ENDDO
          enddo
        enddo
      ENDIF                            !==================================

      m = NextZ(m)
    ENDDO
  end subroutine NBSM_LoopNBSMask

#ifdef CFDACM
  !********************************************************************** 
  ! particle fine search (Moving-particles with Fixed-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch1(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::pid1Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,ovrlp,drlub

    pid1Temp= pid1- nFixed
    dr= GPrtcl_PosR(pid1Temp)%w + GPFix_PosR(pid2)%w
    drlub= dr+ dlub_pp(GPrtcl_pType(pid1Temp),GPFix_pType(pid2))
    dr2= drlub*drlub

    dx= GPrtcl_PosR(pid1Temp)%x - GPFix_PosR(pid2)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GPFix_PosR(pid2)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GPFix_PosR(pid2)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return

    ovrlp = dr-sqrt(d2sum)  
    if(ovrlp>=zero) then
      call GPPW_CntctList%AddContactPPFix(pid1Temp,pid2,ovrlp)
    else
      call GPPW_CntctList%AddLubForcePPFix(pid1Temp,pid2,-ovrlp) 
    endif
  end subroutine NBSM_FineSearch1 

  !********************************************************************** 
  ! particle fine search(Moving-particles with Moving-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch2(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::pid1Temp,pid2Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,drlub,ovrlp

    pid1Temp= pid1- nFixed
    pid2Temp= pid2- nFixed
    dr= GPrtcl_PosR(pid1Temp)%w + GPrtcl_PosR(pid2Temp)%w
    drlub= dr+ dlub_pp(GPrtcl_pType(pid1Temp),GPrtcl_pType(pid2Temp))
    dr2= drlub*drlub

    dx= GPrtcl_PosR(pid1Temp)%x - GPrtcl_PosR(pid2Temp)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GPrtcl_PosR(pid2Temp)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GPrtcl_PosR(pid2Temp)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    ovrlp = dr-sqrt(d2sum)        

    if(ovrlp>=zero) then
      ! this is a convention, the lower id should be the first item in the contact pair (particle & particle)
      if(GPrtcl_id(pid1Temp) < GPrtcl_id(pid2Temp) ) then
        call GPPW_CntctList%AddContactPP(pid1Temp,pid2Temp,ovrlp)
      else
        call GPPW_CntctList%AddContactPP(pid2Temp,pid1Temp,ovrlp)
      endif
    else
      call GPPW_CntctList%AddLubForcePP(pid1Temp,pid2Temp,-ovrlp)
    endif
  end subroutine NBSM_FineSearch2

  !********************************************************************** 
  ! particle fine search (Moving-particles with Ghost-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch3(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::gid,pid1Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,drlub,ovrlp

    pid1Temp= pid1- nFixed
    gid     = pid2- nFixedAndlocal
    dr= GPrtcl_PosR(pid1Temp)%w + GhostP_PosR(gid)%w
    drlub= dr+ dlub_pp(GPrtcl_pType(pid1Temp),GhostP_pType(gid))
    dr2= drlub*drlub

    dx= GPrtcl_PosR(pid1Temp)%x - GhostP_PosR(gid)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GhostP_PosR(gid)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GhostP_PosR(gid)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    ovrlp = dr-sqrt(d2sum)
    
    if(ovrlp>=zero) then
      call GPPW_CntctList%AddContactPPG(pid1Temp,gid,ovrlp)
    else
      call GPPW_CntctList%AddLubForcePPG(pid1Temp,gid,-ovrlp) 
    endif
  end subroutine NBSM_FineSearch3
#else
  !********************************************************************** 
  ! particle fine search (Moving-particles with Fixed-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch1(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::pid1Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,ovrlp

    pid1Temp= pid1- nFixed
    dr= GPrtcl_PosR(pid1Temp)%w + GPFix_PosR(pid2)%w
    dr2= dr*dr
    
    dx= GPrtcl_PosR(pid1Temp)%x - GPFix_PosR(pid2)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GPFix_PosR(pid2)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GPFix_PosR(pid2)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    
    ovrlp = dr-sqrt(d2sum)  
    call GPPW_CntctList%AddContactPPFix(pid1Temp,pid2,ovrlp)
  end subroutine NBSM_FineSearch1 

  !********************************************************************** 
  ! particle fine search(Moving-particles with Moving-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch2(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::pid1Temp,pid2Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,ovrlp

    pid1Temp= pid1- nFixed
    pid2Temp= pid2- nFixed
    dr= GPrtcl_PosR(pid1Temp)%w + GPrtcl_PosR(pid2Temp)%w
    dr2= dr*dr
    dx= GPrtcl_PosR(pid1Temp)%x - GPrtcl_PosR(pid2Temp)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GPrtcl_PosR(pid2Temp)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GPrtcl_PosR(pid2Temp)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    ovrlp = dr-sqrt(d2sum)        

    ! this is a convention, the lower id should be the first item in the contact pair (particle & particle)
    if(GPrtcl_id(pid1Temp) < GPrtcl_id(pid2Temp) ) then
      call GPPW_CntctList%AddContactPP(pid1Temp,pid2Temp,ovrlp)
    else
      call GPPW_CntctList%AddContactPP(pid2Temp,pid1Temp,ovrlp)
    endif
  end subroutine NBSM_FineSearch2 

  !********************************************************************** 
  ! particle fine search (Moving-particles with Ghost-particles)
  !**********************************************************************    
  subroutine NBSM_FineSearch3(this,pid1,pid2)
    implicit none
    class(NBS_Munjiza):: this
    integer,intent(in):: pid1,pid2

    ! locals
    integer::gid,pid1Temp
    real(RK)::dx,dy,dz,dr,d2sum,dr2,ovrlp

    pid1Temp= pid1- nFixed
    gid     = pid2- nFixedAndlocal
    dr= GPrtcl_PosR(pid1Temp)%w + GhostP_PosR(gid)%w
    dr2= dr*dr
    dx= GPrtcl_PosR(pid1Temp)%x - GhostP_PosR(gid)%x
    d2sum=dx*dx;             if(d2sum>dr2) return
    dy= GPrtcl_PosR(pid1Temp)%y - GhostP_PosR(gid)%y
    d2sum=dy*dy+d2sum;       if(d2sum>dr2) return
    dz= GPrtcl_PosR(pid1Temp)%z - GhostP_PosR(gid)%z
    d2sum=dz*dz+d2sum;       if(d2sum>dr2) return
    ovrlp = dr-sqrt(d2sum)  
    call GPPW_CntctList%AddContactPPG(pid1Temp,gid,ovrlp)
  end subroutine NBSM_FineSearch3 
#endif
end module Prtcl_NBS_Munjiza
