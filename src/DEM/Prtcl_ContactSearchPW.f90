module Prtcl_ContactSearchPW
  use MPI
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Geometry
  use Prtcl_Property
  use Prtcl_Variables
  use Prtcl_CL_and_CF
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
  use Prtcl_Decomp_2d,only:int_type,real_type
#else
  use Prtcl_Decomp_2d,only:nrank,int_type,real_type
#endif
  implicit none
  private
    
  integer,dimension(:),allocatable:: Bucket_PWCS ! Head for plane-wall contact search
  integer,dimension(:),allocatable:: id_Wall     ! particle id
  integer,dimension(:),allocatable:: Next_PWCS

  type::ContactSearchPW
    integer :: mHead
    integer :: NextInsert
    real(RK):: MaxWallVel  ! new added
    real(RK):: DeltaXmax1  ! new added
    real(RK):: DeltaXmax2  ! new added
    integer :: Next_iter_update = 0
    integer :: max_nearPrtcl_pWall
    integer :: num_nearPrtcl_pWall = 0
  contains
    procedure:: InitContactSearchPW
    procedure:: FindContactsPW
    procedure:: UpdateNearPrtclsPW
    procedure:: Init_PWCS_List
    procedure:: Reallocate_PWCS_List
    procedure:: reallocate_Bucket
    procedure:: InsertNearPW
    procedure:: copy => CSPW_copy
  end type ContactSearchPW
  type(ContactSearchPW),public::DEMContactSearchPW
    
contains

  !**********************************************************************
  ! Initializing the object
  !**********************************************************************
  subroutine InitContactSearchPW(this )
    implicit none
    class(ContactSearchPW)::this
    type(real3):: wallvel
    integer:: i,nw,iErr1,iErr2,iErr3,iErrSum,ierror
    real(RK)::MaxWallVel,maxvel1
        
    this%max_nearPrtcl_pWall = int(1.5_RK*GPrtcl_list%mlocal)+1
    this%mHead=GPrtcl_list%mlocal
    allocate(Bucket_PWCS(this%mHead), Stat=iErr1)
    allocate(id_Wall(this%max_nearPrtcl_pWall ), Stat=iErr2)
    allocate(Next_PWCS(this%max_nearPrtcl_pWall ), Stat=iErr3)
    iErrSum = abs(iErr1)+abs(iErr2)+abs(iErr3)
    if(iErrSum/=0 ) call DEMLogInfo%CheckForError(ErrT_Abort,"Initialize Contact List","Allocation failed ")
    call this%Init_PWCS_List()
        
    MaxWallVel = 0.0_RK
    nw=DEMGeometry%nPW_local
    do i=1,nw
      wallvel=DEMGeometry%pWall(i)%trans_vel
      MaxWallVel = max(MaxWallVel,norm(wallvel)) 
    enddo
    call MPI_ALLREDUCE(MaxWallVel, maxvel1,1, real_type,MPI_MAX,MPI_COMM_WORLD,ierror)

    this%MaxWallVel= maxvel1
    this%DeltaXmax2= maxval(DEMProperty%Prtcl_PureProp%Radius)*(DEM_opt%Wall_neighbor_ratio+1.0_RK)
#ifndef CFDACM
    this%DeltaXmax1= this%DeltaXmax2-maxval(DEMProperty%Prtcl_PureProp%Radius)
#else
    this%DeltaXmax1= this%DeltaXmax2-maxval(DEMProperty%Prtcl_PureProp%Radius)-maxval(dlub_pw)
#endif
  end subroutine InitContactSearchPW
    
  !**********************************************************************
  ! Init_PWCS_List
  !**********************************************************************
  subroutine Init_PWCS_List(this)
    implicit none
    class(ContactSearchPW) :: this 
    integer::i
       
    Bucket_PWCS = 0
    do i=1,this%max_nearPrtcl_pWall 
      Next_PWCS(i)=-i-1
    enddo
    Next_PWCS(this%max_nearPrtcl_pWall) = 0
    this%NextInsert = 1
  end subroutine Init_PWCS_List
  !**********************************************************************
  ! Reallocate_Bucket
  !********************************************************************** 
  subroutine Reallocate_Bucket(this,nB_new)
    implicit none
    class(ContactSearchPW)::this 
    integer,intent(in):: nB_new

    ! locals
    integer::sizep,sizen
    integer,dimension(:),allocatable:: IntVec

    sizep= this%mHead
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen, nB_new+1)
    sizen= min(sizen,DEM_Opt%numPrtcl)
    this%mHead = sizen

    call move_alloc(Bucket_PWCS, IntVec)
    allocate(Bucket_PWCS(sizen))
    Bucket_PWCS(1:sizep)=IntVec
    Bucket_PWCS(sizep+1:sizen)=0     ! Added at 22:37, 2020-11-06, Gong Zheng
    deallocate(IntVec)

  end subroutine Reallocate_Bucket 
    
  !**********************************************************************
  ! Reallocate_PWCS_List
  !**********************************************************************    
  subroutine Reallocate_PWCS_List(this)
    implicit none
    class(ContactSearchPW)::this 
      
    ! locals
    integer::i
    integer:: sizep,sizen
    integer,dimension(:),allocatable:: IntVec 
      
    sizep=this%max_nearPrtcl_pWall
    sizen=int(1.2_RK*real(sizep,kind=RK))
    this%max_nearPrtcl_pWall=sizen
      
    call move_alloc(id_Wall,IntVec)
    allocate(id_Wall(sizen))
    id_Wall(1:sizep)=IntVec
      
    call move_alloc(Next_PWCS,IntVec)
    allocate(Next_PWCS(sizen))
    Next_PWCS(1:sizep)=IntVec
    deallocate(IntVec)
    
    do i=sizep+1,sizen
      Next_PWCS(i)=-i-1
    enddo
    Next_PWCS(sizen)=0
    this%NextInsert=sizep+1
     
    call DEMLogInfo%CheckForError(ErrT_Pass,"Reallocate_PWCS_List","Need to reallocate particle variables")
    call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)    
  end subroutine Reallocate_PWCS_List
    
  !**********************************************************************
  ! Finding particles which are near all walls
  !**********************************************************************
  subroutine UpdateNearPrtclsPW(this,iterNumber)
    implicit none
    class(ContactSearchPW) :: this
    integer,intent(in):: iterNumber

    ! locals
    real(RK)::max_v,max_a,t,maxreal 
    integer:: nw,i,numIter,wid,nNear,nextI,ierror
    type(real3):: min_point, max_point, pmin_point, pmax_point

    ! checks if in this iteration the neighbour list should be updated
    if(iterNumber<this%Next_iter_update) return
    call this%Init_PWCS_List()
                
    nNear=0
    nw = DEMGeometry%nPW_local
    do wid = 1,nw
      min_point= DEMGeometry%pWall(wid)%min_point
      max_point= DEMGeometry%pWall(wid)%max_point
      pmin_point = min_point - this%DeltaXmax2*real3(1.0_RK,1.0_RK,1.0_RK)
      pmax_point = max_point + this%DeltaXmax2*real3(1.0_RK,1.0_RK,1.0_RK)

      do i=1,GPrtcl_list%nlocal
        IF(GPrtcl_PosR(i)%x >= pmin_point%x .and. GPrtcl_PosR(i)%x <= pmax_point%x .and. &
           GPrtcl_PosR(i)%y >= pmin_point%y .and. GPrtcl_PosR(i)%y <= pmax_point%y .and. &
           GPrtcl_PosR(i)%z >= pmin_point%z .and. GPrtcl_PosR(i)%z <= pmax_point%z) THEN
          nNear=nNear+1
          if(this%NextInsert==0) call this%Reallocate_PWCS_List()
          nextI=this%NextInsert
          this%NextInsert= -Next_PWCS(nextI)
          id_Wall(nextI)= wid
          Next_PWCS(nextI)=Bucket_PWCS(i)
          Bucket_PWCS(i)=nextI
        ENDIF
      enddo
    enddo
    this%num_nearPrtcl_pWall = nNear
    
    ! calculating number of iterations which should be performed
    ! untill the next update of the list of particles near the wall  
    
    ! calculating the maximum velocity in the system
    maxreal = 0.0_RK
    do i=1,GPrtcl_list%nlocal
      maxreal = max(maxreal, norm(GPrtcl_linVel(1,i)))
    enddo
    call MPI_ALLREDUCE(maxreal,max_v,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
    max_v = max_v + this%MaxWallVel
    max_v = max(max_v, 0.3_RK)
    max_v = max_v*1.5_RK
    
    !calculating the maximum acceleration in the system
    maxreal = 0.0_RK
    do i=1,GPrtcl_list%nlocal
      maxreal = max(maxreal, norm(GPrtcl_linAcc(1,i)))
    enddo
    call MPI_ALLREDUCE(maxreal,max_a,1,real_type,MPI_MAX,MPI_COMM_WORLD,ierror)
    max_a = max(max_a, norm(DEM_Opt%gravity))
    max_a = max_a * 1.5_RK
    
    ! solving the equation of motion for t
    ! or  0.5*a*t^2 + v*t - dx = 0 
    ! This equation has two roots, a positive and a negative
    ! we need the positive root
    if(nrank==0) then
      t = (sqrt(max_v**2+ 2.0_RK*max_a*this%DeltaXmax1)-max_v)/max_a
      numIter = max(int(t/DEM_opt%dt), 1)
      numIter = min(numIter, DEM_opt%Wall_max_update_iter)
    endif
    call MPI_BCAST(numIter,1,int_type,0,MPI_COMM_WORLD,ierror)
    this%Next_iter_update =  iterNumber + numIter
    
    if(nrank==0 .and. numIter<DEM_opt%Wall_max_update_iter ) then
      call DEMLogInfo%OutInfo("Neighbor particles of all walls are updated", 4)
      call DEMLogInfo%OutInfo("This process is repeated in Iteration :"//trim(num2str(DEMContactSearchPW%Next_iter_update)),4)
    endif    
  end subroutine UpdateNearPrtclsPW

  !**********************************************************************
  ! CSPW_copy
  !**********************************************************************
  subroutine CSPW_copy(this,id1,id2)
    implicit none
    class(ContactSearchPW):: this
    integer,intent(in)::id1,id2

    ! locals
    integer::n,NextI
   
    n=Bucket_PWCS(id1)
    do while(n>0)
      NextI = Next_PWCS(n)
      Next_PWCS(n) = -this%nextInsert
      this%nextInsert = n      
      n = NextI
    enddo

    Bucket_PWCS(id1) = Bucket_PWCS(id2)
    Bucket_PWCS(id2) = 0

  end subroutine CSPW_copy

  !**********************************************************************
  ! InsertNearPW
  !**********************************************************************
  subroutine InsertNearPW(this,pid)
    implicit none
    class(ContactSearchPW):: this
    integer,intent(in)::pid

    ! locals
    integer:: nw,wid,nextI
    type(real3):: min_point, max_point, pmin_point, pmax_point

    nw = DEMGeometry%nPW_local
    do wid= 1,nw
      min_point= DEMGeometry%pWall(wid)%min_point
      max_point= DEMGeometry%pWall(wid)%max_point
      pmin_point= min_point- this%DeltaXmax2*real3(1.0_RK,1.0_RK,1.0_RK)
      pmax_point= max_point+ this%DeltaXmax2*real3(1.0_RK,1.0_RK,1.0_RK) 

      IF(GPrtcl_PosR(pid)%x >=pmin_point%x .and. GPrtcl_PosR(pid)%x <=pmax_point%x .and. &
         GPrtcl_PosR(pid)%y >=pmin_point%y .and. GPrtcl_PosR(pid)%y <=pmax_point%y .and. &
         GPrtcl_PosR(pid)%z >=pmin_point%z .and. GPrtcl_PosR(pid)%z <=pmax_point%z) THEN

        if(this%NextInsert==0) call this%Reallocate_PWCS_List()
        nextI=this%NextInsert
        this%NextInsert= -Next_PWCS(nextI)
        id_Wall(nextI)= wid
        Next_PWCS(nextI)=Bucket_PWCS(pid)
        Bucket_PWCS(pid)=nextI
      ENDIF
    enddo
  end subroutine InsertNearPW

#ifdef CFDACM
  !**********************************************************************
  ! Performing contact search to determine particle-wall contacts 
  !**********************************************************************
  subroutine FindContactsPW(this)
    implicit none
    class(ContactSearchPW):: this

    ! locals
    logical:: Iscntct,clcLubFlag
    integer:: i,nid,wid
    real(RK):: ovrlp
    type(real3):: nv

    ! this is a convention, the particle id should be the first item in the contact pair (particle & wall)
    DO i=1,GPrtcl_list%nlocal
       nid=Bucket_PWCS(i)
       do while(nid>0)
         wid=id_Wall(nid)
         Iscntct= DEMGeometry%pWall(wid)%isInContact(GPrtcl_PosR(i),ovrlp,nv,clcLubFlag)
         if(Iscntct)then
           call GPPW_CntctList%AddContactPW(i,wid,ovrlp,nv)
         elseif(clcLubFlag) then
           ovrlp= -ovrlp
           if(ovrlp<=dlub_pw(GPrtcl_pType(i))) call GPPW_CntctList%AddLubForcePW(i,wid,ovrlp) 
         endif 
         nid=Next_PWCS(nid)
       enddo
    ENDDO
  end subroutine FindContactsPW
#else
  !**********************************************************************
  ! Performing contact search to determine particle-wall contacts 
  !**********************************************************************
  subroutine FindContactsPW(this)
    implicit none
    class(ContactSearchPW):: this

    ! locals
    integer:: i,nid,wid
    real(RK):: ovrlp
    type(real3):: nv

    ! this is a convention, the particle id should be the first item in the contact pair (particle & wall)
    DO i=1,GPrtcl_list%nlocal
       nid=Bucket_PWCS(i)
       do while(nid>0)
         wid=id_Wall(nid)
         if(DEMGeometry%pWall(wid)%isInContact(GPrtcl_PosR(i),ovrlp,nv))then
           call GPPW_CntctList%AddContactPW(i,wid,ovrlp,nv) 
         endif 
         nid=Next_PWCS(nid)
       enddo
    ENDDO
  end subroutine FindContactsPW
#endif
    
end module Prtcl_ContactSearchPW
