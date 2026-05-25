#include "definitions_inc.f90"
#define DEM_ncv_Allowed 1000
module sp_CL_and_CF
  use MPI
  use mc_TypeDef
  use mc_LogInfo
#if defined(CFDDEM) || defined(CFDACM)
  use mc_Decomp2d,only: nrank,nproc
#else
  use sp_decomp_2d,only: nrank,nproc
#endif
  use sp_Property
  use sp_Geometry
  use sp_Variables
  use sp_Parameters
#ifdef CFDACM
  use f2_Parameters,only:dtMax
#endif
  implicit none
  private 
  real(RK),parameter::END_OF_PRTCL = -huge(Pi) + Pi   ! particle IO end flag
  
  integer,dimension(:),allocatable:: Bucket
  integer,dimension(:),allocatable:: id_j
  integer,dimension(:),allocatable:: Next
  integer,dimension(:),allocatable:: CntctStatus
  type(real3),dimension(:),allocatable:: TanDelta
  real(RK),dimension(:),allocatable::VelRel_Init

  integer,dimension(:),allocatable:: id_i
  integer,dimension(:),allocatable:: Head_Cp  ! cp: counterpart
  integer,dimension(:),allocatable:: Next_Cp
    
  type ContactList
    integer:: mBucket
    integer:: numCntcts(4)
    integer:: max_numCntcts
    integer:: NextInsert
  contains
    procedure:: InitContactList => CL_InitContactList
    procedure:: reallocateCL    => CL_reallocateCL
    procedure:: AddContactPP    => CL_AddContactPP
    procedure:: AddContactPPG   => CL_AddContactPPG
    procedure:: AddContactPPFix => CL_AddContactPPFix
    procedure:: AddContactPW    => CL_AddContactPW
    procedure,private:: Find_Insert
    procedure:: PreIteration    => CL_PreIteration
    procedure:: RemvReleased    => CL_RemvReleased
    procedure:: Add_Cntctlink   => CL_Add_Cntctlink
  
    procedure:: Prepare_Restart => CL_Prepare_Restart
    procedure:: GetNextTanDel_Un=> CL_GetNextTanDel_Un
    procedure:: Add_RestartCntctlink=> CL_Add_RestartCntctlink
#ifdef CFDACM
    procedure,nopass:: AddLubForcePP    => CL_AddLubForcePP
    procedure,nopass:: AddLubForcePPG   => CL_AddLubForcePPG
    procedure,nopass:: AddLubForcePPFix => CL_AddLubForcePPFix
    procedure,nopass:: AddLubForcePW    => CL_AddLubForcePW
#endif
  end type ContactList
  type(ContactList):: GPPW_CntctList 
    
  public:: GPPW_CntctList,END_OF_PRTCL
  public:: print_ContactList, resemble_ContactList, count_ContactList, copy_ContactList
  public:: gather_ContactList,get_numContacts, getPrtcl_nlink, IsPrtclContact
contains
    
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Initializing the contact list for ContactList class
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_InitContactList(this); implicit none
    class(ContactList) :: this
    integer:: max_numCntcts,i,ierr,ierrTmp
    
    this%numCntcts= 0
    max_numCntcts = DEM_opt%cntctList_Size * GPrtcl_list%mlocal

    this%Max_numCntcts = max_numCntcts
    this%mBucket = GPrtcl_list%mlocal
        
    ! initializing the linked list to stores id pairs
    ierr = 0
    allocate(Bucket(this%mBucket),      Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(id_j(max_numCntcts),       Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(CntctStatus(max_numCntcts),Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(Next(max_numCntcts),       Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(TanDelta(max_numCntcts),   Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(VelRel_Init(max_numCntcts),Stat=ierrTmp); ierr =ierr +abs(ierrTmp)

    allocate(id_i(max_numCntcts),       Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(Head_Cp(this%mBucket),     Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(Next_cp(max_numCntcts),    Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    if(ierr/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"CL_InitContactList","Allocation failed ")
        
    Bucket = 0
    do i = 1, max_numCntcts
      Next(i)=-i-1
      CntctStatus(i)=-1
    enddo
    Next(max_numCntcts) = 0
    this%NextInsert = 1
    Head_Cp = -1
    Next_cp = -1
  end subroutine CL_InitContactList

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! reallocate contact list
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_reallocateCL(this); implicit none
    class(ContactList):: this

    ! locals
    integer::sizep,sizen,ierrTmp,ierr
    integer,dimension(:),allocatable:: IntVec

    sizep= this%mBucket
    sizen= GPrtcl_list%mlocal    ! 2021-09-08, Zheng Gong
    this%mBucket = sizen

    ! ======= integer vector part =======
    ierr=0
    call move_alloc(Bucket, IntVec)
    allocate(Bucket(sizen),stat=ierrTmp);ierr=ierr+abs(ierrTmp)
    Bucket(1:sizep)=IntVec
    Bucket(sizeP+1:sizen)=0    ! Added at 11:24, 2020-09-11, Gong Zheng

    call move_alloc(Head_Cp,IntVec)
    allocate(Head_Cp(sizen),stat=ierrTmp);ierr=ierr+abs(ierrTmp)
    Head_Cp(1:sizep)=IntVec
    Head_Cp(sizep+1:sizen)=-1  ! Added at 11:24, 2020-09-11, Gong Zheng
    deallocate(IntVec)
    
    if(ierr/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort," CL_reallocateCL"," Reallocate wrong!")
      call DEMLogInfo%OutInfo("The present processor  is :"//strip(num2str(nrank)),3)
    endif   
   ! here nothing is done for id_j, Next, CntctStatus, TanDelta, VelRel_Init
   ! considering that Max_numCntcts is big enough.
   ! If not, a fatal error will occur in CL_AddContactPP/CL_AddContactPW/CL_AddContactPPG/CL_AddContactPPFix
  end subroutine CL_reallocateCL

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Adding a contact pair to the contact list (particle & particle)
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_AddContactPP(this,id1,id2,ovrlp); implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    real(RK),intent(in)::ovrlp
    integer::item1,item2,gid2

    gid2=GPrtcl_id(id2)
    call this%Find_Insert(id1,gid2,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(1) = this%numCntcts(1) + 1
      call ContactForce_PP(id1,id2,item2,ovrlp)

      id_i(item2)= GPrtcl_id(id1)
      Next_Cp(item2)=Head_Cp(id2)
      Head_Cp(id2)=item2

    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPP","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPP","The container is full and there is no space for new item" )      
    endif
  end subroutine CL_AddContactPP

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Adding a contact pair to the contact list (particle & ghost particle) 
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_AddContactPPG(this,id1,id2,ovrlp); implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    real(RK),intent(in)::ovrlp
    integer::item1,item2,gid2

    gid2=GhostP_id(id2)
    call this%Find_Insert(id1,gid2,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(2) = this%numCntcts(2) + 1
      call ContactForce_PPG(id1,id2,item2,ovrlp)

    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPG","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPG","The container is full and there is no space for new item" )      
    endif
  end subroutine CL_AddContactPPG

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Adding a contact pair to the contact list (particle & fixed particle) 
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_AddContactPPFix(this,id1,id2,ovrlp); implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    real(RK),intent(in)::ovrlp
    integer::item1,item2,gid2

    gid2=GPFix_id(id2)
    call this%Find_Insert(id1,gid2,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(3) = this%numCntcts(3) + 1
      call ContactForce_PPFix(id1,id2,item2,ovrlp)

    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPFix","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPPFix","The container is full and there is no space for new item" )      
    endif
  end subroutine CL_AddContactPPFix

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Adding a contact pair to the contact list 
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_AddContactPW(this,id1,id2,ovrlp,nv); implicit none
    class(ContactList):: this
    integer,intent(in):: id1,id2
    integer::wid,item1,item2
    real(RK)::ovrlp
    type(real3)::nv

    wid= DEMGeometry%pWall(id2)%wall_id
    call this%Find_Insert(id1,wid,item1,item2)
    
    ! item1 is the status of the item insertion. 1:old, 2:new 
    ! item2 is the container index
    if(item1>0) then
      CntctStatus(item2)=item1
      this%numCntcts(4) = this%numCntcts(4) + 1
      call ContactForce_PW(id1,id2,item2,ovrlp,nv)
    elseif(item1 == -1 ) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPW","The inserted item's id is greater than allowed value:"//num2str(id1))
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"CL_AddContactPW","The container is full and there is no space for new item")      
    endif
  end subroutine CL_AddContactPW
    
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  !* searching for an item, if it exists, it would return (/1, index of container that contain the item/),
  !     if it is new, pushing the item into the list and return (/2, index of container that
  !     contains the new item /), otherwise (-2,-2). 
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Find_insert(this,bkt_id,b_value,item1,item2); implicit none
    class(ContactList) this
    integer,intent(in) :: bkt_id,b_value
    integer,intent(out):: item1, item2

    ! locals
    integer::n,NextI
        
    ! The bucket id is not in the range and must return (-1,-1), this is an error for the linked list
    if(bkt_id>GPrtcl_list%mlocal) then
      item1 =-1; item2 = -1
      return
    end if
        
     n = Bucket(bkt_id) 
     ! the item already exists in the list; returning the proper code and the index of container in which the item exists
    do while(n>0)
      if(id_j(n)==b_value) then
        item1 = 1; item2 = n
        return
      end if
      n = Next(n)
    enddo
        
    !the container is full and there is no more space for this item
    if(this%NextInsert==0) then
      item1 =-2; item2 = -2
      return
    endif
        
    ! the item is new and it should be pushed into the list, inserting new item in the list 
    NextI = this%NextInsert     
    this%NextInsert = -Next(NextI)
    id_j(NextI)=b_value
    Next(NextI)=Bucket(bkt_id)
    Bucket(bkt_id) = NextI
    TanDelta(NextI)= zero_r3
    VelRel_Init(NextI)=-2.0_RK
    item1 =2; item2 = NextI
  end subroutine Find_insert

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! getPrtcl_nlink
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  function getPrtcl_nlink(pid) result(res); implicit none
    integer,intent(in)::pid

    ! locals
    integer::res,n

    res = 0
    n = Bucket(pid)
    do while(n>0)
      res = res + 1
      n = Next(n)
    enddo

    n = Head_cp(pid)
    do while(n .ne. -1)
      res = res + 1
      n = Next_Cp(n)
    enddo
  end function getPrtcl_nlink

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! IsPrtclContact
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  function IsPrtclContact(pid) result(res); implicit none
    integer,intent(in)::pid

    ! locals
    integer::res

    if(Bucket(pid)>0 .or. Head_cp(pid)>0) then
      res= 1
    else
      res= 0
    endif
  end function IsPrtclContact

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CL_Prepare_Restart
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_Prepare_Restart(this,nlink_ind); implicit none
    class(ContactList)::this
    integer,intent(in)::nlink_ind

    ! locals
    integer::i,CLid

    CLid = nlink_ind
    do i=1,this%Max_numCntcts
      if(CntctStatus(i)>0) then
        CLid = CLid + 1
        CntctStatus(i) = CLid 
      endif
    enddo
  end subroutine CL_Prepare_Restart

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! get_numContacts
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine get_numContacts(nCntct,nTanDel); implicit none
    integer,intent(inout):: nCntct,nTanDel

    ! locals
    integer:: pid,n

    nCntct=0; nTanDel=0
    DO pid=1,GPrtcl_list%nlocal
      n=Bucket(pid)
      do while(n>0)
        nTanDel=nTanDel+1; n=Next(n)
      enddo
     
      n=Head_cp(pid)
      do while(n>0)
        nCntct=nCntct+1; n=Next_cp(n)
      enddo
    ENDDO
    nCntct=nCntct+nTanDel
  end subroutine get_numContacts

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CL_GetNextTanDel_Un
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_GetNextTanDel_Un(this,TanDel_Un,prev,now); implicit none
    class(ContactList)::this 
    type(real4),intent(out)::TanDel_Un
    integer,intent(in)::prev
    integer,intent(out)::now

    ! locals
    integer::i
    
    do i=prev,this%Max_numCntcts
      if(CntctStatus(i)>0) then
        TanDel_Un   = TanDelta(i)
        TanDel_Un%w = VelRel_Init(i)
        now=i;exit
      endif
    enddo
  end subroutine CL_GetNextTanDel_Un

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! gather_ContactList
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine gather_ContactList(pid,buf_send,m); implicit none
    integer,intent(in)::pid
    real(RK),dimension(*),intent(out)::buf_send
    integer,intent(inout)::m

    ! locals
    integer::n

    n = Bucket(pid)
    do while(n>0)
      CntctStatus(n)=3  ! 2021-08-25,added
      buf_send(m)=real(id_j(n),RK);  m=m+1
      buf_send(m)=TanDelta(n)%x;     m=m+1
      buf_send(m)=TanDelta(n)%y;     m=m+1
      buf_send(m)=TanDelta(n)%z;     m=m+1
      buf_send(m)=VelRel_Init(n);    m=m+1
      n = Next(n)
    enddo        

    n = Head_cp(pid)
    do while(n .ne. -1)
      buf_send(m)=real(id_i(n),RK);  m=m+1
      buf_send(m)=TanDelta(n)%x;     m=m+1
      buf_send(m)=TanDelta(n)%y;     m=m+1
      buf_send(m)=TanDelta(n)%z;     m=m+1  
      buf_send(m)=VelRel_Init(n);    m=m+1   
      n = Next_Cp(n)
    enddo
  end subroutine gather_ContactList

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! copy_ContactList
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine copy_ContactList(id1,id2); implicit none
    integer,intent(in)::id1,id2

    bucket(id1) = bucket(id2)
    Head_Cp(id1)= Head_Cp(id2)
    bucket(id2) = 0
    Head_Cp(id2)=-1
  end subroutine copy_ContactList

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CL_Add_Cntctlink
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_Add_Cntctlink(this,pid,buf_recv,m); implicit none
    class(ContactList)::this
    integer,intent(in)::pid
    integer,intent(inout)::m
    real(RK),dimension(*),intent(in)::buf_recv

    ! locals
    integer::NextI
    real(RK)::realt

    do
      realt = buf_recv(m); m=m+1
      if(realt == END_OF_PRTCL) return

      NextI = this%NextInsert
      this%NextInsert = -Next(NextI)
      if(this%NextInsert==0) then
        call DEMLogInfo%CheckForError(ErrT_Abort,"CL_Add_Cntctlink","The container is full and there is no space for new item") 
      endif
      id_j(NextI) = nint(realt)
      Next(NextI) = Bucket(pid)
      Bucket(pid) = NextI
      CntctStatus(NextI)=2

      TanDelta(NextI)%x = buf_recv(m); m=m+1
      TanDelta(NextI)%y = buf_recv(m); m=m+1
      TanDelta(NextI)%z = buf_recv(m); m=m+1
      VelRel_Init(NextI)= buf_recv(m); m=m+1
    enddo
  end subroutine CL_Add_Cntctlink

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Removing all released contacts (those which are not in contact in this time step) 
  ! from contact list.
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_RemvReleased(this); implicit none
    class(ContactList)::this
    
    ! locals
    integer::pid,n,prev,NextI
    
    do pid=1,GPrtcl_list%nlocal ! Modify this part by Zheng Gong, at 2020-09-09
      prev=0
      n=Bucket(pid)
      do while(n>0)
        NextI=Next(n)
        if(CntctStatus(n)==-2) then
          if(prev==0) then
            Bucket(pid)=NextI
          else
            Next(prev)=NextI
          endif
          CntctStatus(n)=-1
        
          Next(n)=-this%nextInsert
          this%nextInsert=n
        else
          prev=n
        endif
        n=NextI
      enddo

      n= Head_cp(pid)
      do while(n>0)
        if(CntctStatus(n)==3) then
          Next(n)=Bucket(pid)
          Bucket(pid)=n
          id_j(n)=id_i(n)
          CntctStatus(n)=2
        endif
        n=Next_cp(n)
      enddo
    enddo
    
    do pid=1,this%Max_numCntcts
      if(CntctStatus(pid)==3) then
        CntctStatus(pid)=-1
        Next(pid)=-this%nextInsert
        this%nextInsert=pid
      endif
    enddo

#ifdef DEBUG_RemvReleased
    ! CntctStatus(pid)==-2 will only appear when some particle escapes the domain
    do pid=1,this%Max_numCntcts
      if(CntctStatus(pid)==-2) then
         print*,nrank,pid,CntctStatus(pid),'############';stop     
      endif
    enddo 
#endif
  end subroutine CL_RemvReleased

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CL_PreIteration
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_PreIteration(this); implicit none
    class(ContactList) this

    ! locals
    integer::i

    this%numCntcts = 0
    do i=1,this%Max_numCntcts
       if(CntctStatus(i)>0) CntctStatus(i) = -2  ! flag previous contact
    enddo
    Head_Cp = -1
    Next_Cp = -1
  end subroutine CL_PreIteration

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! CL_Add_RestartCntctlink
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_Add_RestartCntctlink(this,pid1,ncv,CntctVec,TanDel_Un); implicit none
    class(ContactList)::this
    integer,intent(in)::pid1,ncv
    integer,dimension(:),intent(in)::CntctVec
    type(real4),dimension(:),intent(in)::TanDel_Un

    ! locals
    real(RK)::dx,dy,dz,dist
    integer::i,j,nlocal,pid2,gid1,gid2,NextI
  
    dist  = 1.0_RK
    gid1  = GPrtcl_id(pid1)
    nlocal= GPrtcl_list%nlocal
    DO j=1,ncv
      pid2=0;  gid2=CntctVec(j)
      DO i=1,nlocal
        if(GPrtcl_id(i)==gid2) then
          pid2=i; exit
        endif
      ENDDO
      IF(pid2>0) THEN
        dx= GPrtcl_PosR(pid1)%x- GPrtcl_PosR(pid2)%x
        dy= GPrtcl_PosR(pid1)%y- GPrtcl_PosR(pid2)%y
        dz= GPrtcl_PosR(pid1)%z- GPrtcl_PosR(pid2)%z
        dist= sqrt(dx*dx+dy*dy+dz*dz)-1.5_RK*(GPrtcl_PosR(pid1)%w+GPrtcl_PosR(pid2)%w)
        if(gid1>gid2 .and. dist<=0.0_RK) cycle
      ENDIF

      NextI = this%NextInsert     
      this%NextInsert = -Next(NextI)
      if(this%NextInsert==0) then
        call DEMLogInfo%CheckForError(ErrT_Abort,"CL_Add_RestartCntctlink","The container is full and there is no space for new item") 
      endif        
      id_j(NextI) = gid2 ! here gid2 can be a particle_within_this_processor/ghost_particle/wall
      Next(NextI) = Bucket(pid1)
      Bucket(pid1)= NextI
      CntctStatus(NextI)=2
      TanDelta(NextI)= TanDel_Un(j)
      VelRel_Init(NextI)= TanDel_Un(j)%w
      if(pid2>0 .and. gid1<gid2 .and. dist<=0.0_RK) then  ! gid2 is also within this processor
        id_i(NextI)= gid1
        Next_Cp(NextI)=Head_Cp(pid2)
        Head_Cp(pid2)=NextI
      endif
    ENDDO
  end subroutine CL_Add_RestartCntctlink

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! count_ContactList
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine count_ContactList(pid,ncv); implicit none
    integer,intent(in)::pid
    integer,intent(out)::ncv

    ! locals
    integer::n
    
    ncv =0
    n = Bucket(pid)
    do while(n>0)
      ncv = ncv + 1
      ! for monosize particles, a particle can contact with NO MORE THAN 12 neighbor particles.
      if(ncv > DEM_ncv_Allowed) call DEMLogInfo%CheckForError(ErrT_Abort,"count_ContactList","so big ncv")
      n = Next(n)
    enddo        

    n = Head_cp(pid)
    do while(n .ne. -1)
      ncv = ncv + 1
      if(ncv > DEM_ncv_Allowed) call DEMLogInfo%CheckForError(ErrT_Abort,"count_ContactList","so big ncv") 
      n = Next_Cp(n)
    enddo
  end subroutine count_ContactList

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! resemble_ContactList
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine resemble_ContactList(pid,ncv,CntctVec); implicit none
    integer,intent(in)::pid
    integer,intent(out)::ncv
    integer,dimension(:),intent(out)::CntctVec

    ! locals
    integer::n
    
    ncv =0
    n = Bucket(pid)
    do while(n>0)
      ncv = ncv + 1
      CntctVec(2*ncv-1) =  id_j(n)
      CntctVec(2*ncv)   =  CntctStatus(n)
      n = Next(n)
    enddo        

    n = Head_cp(pid)
    do while(n .ne. -1)
      ncv = ncv + 1
      CntctVec(2*ncv-1) =  id_i(n)
      CntctVec(2*ncv)   =  CntctStatus(n)  
      n = Next_Cp(n)
    enddo
  end subroutine resemble_ContactList

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! print_ContactList
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine print_ContactList(itime); implicit none
    integer,intent(in)::itime
   
    ! locals
    integer::i,n,pid,ierr,iUnit
    character(len=:),allocatable::filename
    
    filename = "cntctlist" // int2str(itime,10) // ".txt"
    if(nrank==0) then  
      open(newunit=iUnit,file=filename,status='replace',form='formatted')
      close(iUnit)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    do i=0,nproc-1
      if(nrank==i) then
        open(newunit=iUnit,file=filename,status='old',position='append',form='formatted')
        DO pid=1,GPrtcl_list%nlocal
          n = Bucket(pid)
          do while(n>0)
            write(iUnit,'(3I10,10ES24.15)')nrank,GPrtcl_id(pid),id_j(n),TanDelta(n),VelRel_Init(n)
            n = Next(n)
          enddo
     
          n = Head_cp(pid)
          do while(n>0)
            write(iUnit,'(3I10,10ES24.15)')nrank,GPrtcl_id(pid),id_i(n),TanDelta(n),VelRel_Init(n)
            n = Next_cp(n)
          enddo
        ENDDO
        close(iUnit)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    enddo
  end subroutine print_ContactList
 
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ContactForce_PP
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ContactForce_PP(pid,pjd,ind,ovrlp); implicit none
    integer,intent(in)::pid,pjd,ind
    real(RK),intent(in)::ovrlp
    
    ! locals
    type(real4)::Posi,Posj
    type(BinaryProperty)::Prop_ij
    real(RK)::ri,rj,fn,ft,vrn,ft_fric,k_n,d_n,k_t,d_t,normTan1,normTan2,Vel_in!,w_hat_mag
    type(real3)::Norm_v,Veli,Velj,Rvei,Rvej,Vrij,Vel_w,Vij_n,Vij_t,Ovlp_t,Fnij,Ftij,Moment!,Mrij,W_hat
#ifdef CFDACM
    real(RK)::TCollision
#endif

    Prop_ij=DEMProperty%Prtcl_BnryProp(GPrtcl_pType(pid), GPrtcl_pType(pjd))
    Veli= GPrtcl_linVel(1,pid)
    Velj= GPrtcl_linVel(1,pjd)
    Rvei= GPrtcl_rotVel(1,pid)
    Rvej= GPrtcl_rotVel(1,pjd)
    Posi= GPrtcl_PosR(pid)
    Posj= GPrtcl_PosR(pjd)
    ri= Posi%w; rj= Posj%w
    Norm_v= Posj .nv. Posi  ! Normal vector, Posj-Posi

#define ContactForce_PP
#ifdef CFDACM
#include "ac_ContactForce_inc.f90"
#else
#include "sp_ContactForce_inc.f90"
#endif
#undef  ContactForce_PP
  end subroutine ContactForce_PP

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ContactForce_PPG
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ContactForce_PPG(pid,pjd,ind,ovrlp); implicit none
    integer,intent(in)::pid,pjd,ind
    real(RK),intent(in)::ovrlp
    
    ! locals
    type(real4)::Posi,Posj
    type(BinaryProperty)::Prop_ij
    real(RK)::ri,rj,fn,ft,vrn,ft_fric,k_n,d_n,k_t,d_t,normTan1,normTan2,Vel_in!,w_hat_mag
    type(real3)::Norm_v,Veli,Velj,Rvei,Rvej,Vrij,Vel_w,Vij_n,Vij_t,Ovlp_t,Fnij,Ftij,Moment!,Mrij,W_hat
#ifdef CFDACM
    real(RK)::TCollision
#endif

    Prop_ij=DEMProperty%Prtcl_BnryProp(GPrtcl_pType(pid), GhostP_pType(pjd))
    if(GPrtcl_id(pid)<GhostP_id(pjd)) then
      Veli= GPrtcl_linVel(1,pid)
      Velj= GhostP_linVel(pjd)
      Rvei= GPrtcl_rotVel(1,pid)
      Rvej= GhostP_rotVel(pjd)
      Posi= GPrtcl_PosR(pid)
      Posj= GhostP_PosR(pjd)
      ri= Posi%w; rj= Posj%w
      Norm_v= Posj .nv. Posi  ! Normal vector, Posj-Posi
#define ContactForce_PPG
#ifdef CFDACM
#include "ac_ContactForce_inc.f90"
#else
#include "sp_ContactForce_inc.f90"
#endif
#undef  ContactForce_PPG
    else
      Veli= GhostP_linVel(pjd)
      Velj= GPrtcl_linVel(1,pid)
      Rvei= GhostP_rotVel(pjd)
      Rvej= GPrtcl_rotVel(1,pid)
      Posi= GhostP_PosR(pjd)
      Posj= GPrtcl_PosR(pid)
      ri= Posi%w; rj= Posj%w
      Norm_v= Posj .nv. Posi  ! Normal vector, Posj-Posi
#define ContactForce_PGP
#ifdef CFDACM
#include "ac_ContactForce_inc.f90"
#else
#include "sp_ContactForce_inc.f90"
#endif
#undef  ContactForce_PGP
    endif
  end subroutine ContactForce_PPG

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ContactForce_PPFix
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ContactForce_PPFix(pid,pjd,ind,ovrlp); implicit none
    integer,intent(in)::pid,pjd,ind
    real(RK),intent(in)::ovrlp
    
    ! locals
    type(real4)::Posi,Posj
    type(BinaryProperty)::Prop_ij
    real(RK)::ri,rj,fn,ft,vrn,ft_fric,k_n,d_n,k_t,d_t,normTan1,normTan2,Vel_in!,w_hat_mag
    type(real3)::Norm_v,Veli,Velj,Rvei,Rvej,Vrij,Vel_w,Vij_n,Vij_t,Ovlp_t,Fnij,Ftij,Moment!,Mrij,W_hat
#ifdef CFDACM
    real(RK)::TCollision
#endif

    Prop_ij=DEMProperty%Prtcl_BnryProp(GPrtcl_pType(pid), GPFix_pType(pjd))
    Veli= GPrtcl_linVel(1,pid)
    Velj= zero_r3
    Rvei= GPrtcl_rotVel(1,pid)
    Rvej= zero_r3
    Posi= GPrtcl_PosR(pid)
    Posj= GPFix_PosR(pjd)
    ri= Posi%w; rj= Posj%w
    Norm_v= Posj .nv. Posi  ! Normal vector, Posj-Posi

#define ContactForce_PPFix_W
#ifdef CFDACM
#include "ac_ContactForce_inc.f90"
#else
#include "sp_ContactForce_inc.f90"
#endif
#undef  ContactForce_PPFix_W
  end subroutine ContactForce_PPFix

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! ContactForce_PW
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine ContactForce_PW(pid,mwi,ind,ovrlp,Norm_v); implicit none
    integer,intent(in)::pid,mwi,ind
    real(RK),intent(in):: ovrlp
    type(real3),intent(inout):: Norm_v
    
    ! locals
    type(BinaryProperty)::Prop_ij
    real(RK)::ri,rj,fn,ft,vrn,ft_fric,k_n,d_n,k_t,d_t,normTan1,normTan2,Vel_in!,w_hat_mag
    type(real3)::Veli,Velj,Rvei,Rvej,Vrij,Vel_w,Vij_n,Vij_t,Ovlp_t,Fnij,Ftij,Moment!,Mrij,W_hat
#ifdef CFDACM
    real(RK)::TCollision
#endif
        
    Prop_ij = DEMProperty%PrtclWall_BnryProp(GPrtcl_pType(pid),DEMGeometry%pWall(mwi)%wall_Type)
    Veli= GPrtcl_linVel(1,pid)
    Velj= DEMGeometry%pWall(mwi)%trans_vel
    Rvei= GPrtcl_rotVel(1,pid)
    Rvej= zero_r3
    ri= GPrtcl_PosR(pid)%w
    rj= 1.00E20_RK*ri
    ! since the normal vector points from particle i to j we must negate the normal vector
    Norm_v = (-1.0_RK)*Norm_v
    
#define ContactForce_PPFix_W
#ifdef CFDACM
#include "ac_ContactForce_inc.f90"
#else
#include "sp_ContactForce_inc.f90"
#endif
#undef  ContactForce_PPFix_W
  end subroutine ContactForce_PW

#ifdef CFDACM
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Adding lubrication force to the "contact list" (particle & particle) 
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_AddLubForcePP(pid,pjd,ovrlp); implicit none
    integer,intent(in):: pid,pjd
    real(RK),intent(in)::ovrlp
    
#ifdef DriftKissTumbleBreugem
    ! locals
    integer::pt_i,pt_j
    type(real3)::Norm_v,LubForce
    type(BinaryProperty)::Prop_ij
    real(RK)::Gravity_Norm,Mass,EpsValue,dlubDist,ratioD
    
    pt_i = GPrtcl_pType(pid)
    pt_j = GPrtcl_pType(pjd)
    dlubDist=dlub_pp(pt_i, pt_j)  
    Prop_ij=DEMProperty%Prtcl_BnryProp(pt_i, pt_j)
    
    EpsValue=1.0E-4
    Gravity_Norm=norm(DEM_Opt%Gravity)
    Mass=2.0_RK*Prop_ij%MassEff
    ratioD=ovrlp/dlubDist-1.0_RK
    ratioD=ratioD*ratioD
    
    ! normal vector, Posj-Posi
    Norm_v = (GPrtcl_PosR(pjd)) .nv. (GPrtcl_PosR(pid))  

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= (Mass*Gravity_Norm/EpsValue)*ratioD*Norm_v
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) - LubForce
    GPrtcl_cntctForce(pjd) = GPrtcl_cntctForce(pjd) + LubForce    
#else
    ! locals
    integer:: pt_i,pt_j
    type(real3)::Norm_v,Vij_n,LubForce

    ! normal vector, Posj-Posi
    pt_i = GPrtcl_pType(pid)
    pt_j = GPrtcl_pType(pjd)
    Norm_v = (GPrtcl_PosR(pjd)) .nv. (GPrtcl_PosR(pid))  

    ! normal relative velocity vectors
    Vij_n = ((GPrtcl_linVel(1,pid)-GPrtcl_linVel(1,pjd)) .dot. Norm_v)*Norm_v

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= LubCoe_pp(pt_i,pt_j)*Vij_n
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) - LubForce
    GPrtcl_cntctForce(pjd) = GPrtcl_cntctForce(pjd) + LubForce
#endif
  end subroutine CL_AddLubForcePP

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Adding lub force to the "contact list" (particle & ghost particle)
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_AddLubForcePPG(pid,gjd,ovrlp); implicit none
    integer,intent(in):: pid,gjd
    real(RK),intent(in)::ovrlp

#ifdef DriftKissTumbleBreugem
    ! locals
    integer::pt_i,pt_j
    type(real3)::Norm_v,LubForce
    type(BinaryProperty)::Prop_ij
    real(RK)::Gravity_Norm,Mass,EpsValue,dlubDist,ratioD
    
    pt_i = GPrtcl_pType(pid)
    pt_j = GhostP_pType(gjd)
    dlubDist=dlub_pp(pt_i, pt_j)  
    Prop_ij=DEMProperty%Prtcl_BnryProp(pt_i, pt_j)
    
    EpsValue=1.0E-4
    Gravity_Norm=norm(DEM_Opt%Gravity)
    Mass=2.0_RK*Prop_ij%MassEff
    ratioD=ovrlp/dlubDist-1.0_RK
    ratioD=ratioD*ratioD
    
    ! normal vector, Posj-Posi
    Norm_v = (GhostP_PosR(gjd)) .nv. (GPrtcl_PosR(pid))  

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= (Mass*Gravity_Norm/EpsValue)*ratioD*Norm_v
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) - LubForce
#else
    ! locals
    integer:: pt_i,pt_j
    type(real3)::Norm_v,Vij_n,LubForce

    ! normal vector, Posj-Posi
    pt_i = GPrtcl_pType(pid)
    pt_j = GhostP_pType(gjd)
    Norm_v = (GhostP_PosR(gjd)) .nv. (GPrtcl_PosR(pid))  

    ! normal relative velocity vectors
    Vij_n = ((GPrtcl_linVel(1,pid)-GhostP_linVel(gjd)) .dot. Norm_v)*Norm_v

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= LubCoe_pp(pt_i,pt_j)*Vij_n
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) - LubForce
#endif
  end subroutine CL_AddLubForcePPG

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Adding lub force to the "contact list" (particle & fixed particle)
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_AddLubForcePPFix(pid,fid,ovrlp); implicit none
    integer,intent(in):: pid,fid
    real(RK),intent(in)::ovrlp

    ! locals
    integer:: pt_i,pt_j
    type(real3)::Norm_v,Vij_n,LubForce

    ! normal vector, Posj-Posi
    pt_i = GPrtcl_pType(pid)
    pt_j = GPFix_pType(fid)
    Norm_v = (GPFix_PosR(fid) .nv. GPrtcl_PosR(pid))  

    ! normal relative velocity vectors
    Vij_n = (GPrtcl_linVel(1,pid) .dot. Norm_v)*Norm_v

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= LubCoe_pp(pt_i,pt_j)*Vij_n
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) -LubForce
  end subroutine CL_AddLubForcePPFix

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Adding lubrication force to the "contact list" (particle & wall)
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine CL_AddLubForcePW(pid,mwi,ovrlp); implicit none
    integer,intent(in)::pid,mwi
    real(RK),intent(in):: ovrlp
    
    ! locals
    type(real3):: Vij_n,LubForce,Norm_v
    
    ! normal and tangential velocities vectors 
    Norm_v= DEMGeometry%pWall(mwi)%n
    Vij_n = ((GPrtcl_linVel(1,pid)- DEMGeometry%pWall(mwi)%trans_vel) .dot. Norm_v)*Norm_v

    ! Here lubrication force is regarded as a special type of "contact force"
    LubForce= LubCoe_pw(GPrtcl_pType(pid)) *Vij_n
    GPrtcl_cntctForce(pid) = GPrtcl_cntctForce(pid) -LubForce
  end subroutine CL_AddLubForcePW
#endif
end module sp_CL_and_CF
#ifdef DEM_ncv_Allowed
#undef DEM_ncv_Allowed
#endif
