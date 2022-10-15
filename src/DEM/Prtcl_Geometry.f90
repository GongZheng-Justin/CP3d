module prtcl_Geometry
  use prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Property
  use Prtcl_decomp_2d
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
#endif
  implicit none
  private
    
  integer:: MaxWallSize=6
  type(real3):: pmin, pmax             ! global domain
  type(real3):: pmin_local, pmax_local ! local domain

  type PlaneWall
    type(real3):: P1       ! first point   
    type(real3):: P2       ! second point  
    type(real3):: P3       ! third point 
    type(real3):: P4       ! fourth point
    type(real3):: trans_vel = zero_r3 ! translational velocity
    integer:: user_id      ! user supplied wall id
    integer:: wall_Type    ! property type of wall material
        
    logical:: bothSide     ! checking if both sides are active 
    logical:: isInfinite
    
    integer:: wall_id      ! program generated wall id
    real (RK):: d          ! d in the implicit equation: ax+by+cz+d = 0   
    type(real3):: n        ! normal vector
    type(real3):: min_point
    type(real3):: max_point
  contains
    procedure:: isInContact => PW_isInContact
    procedure:: IsInDomain  => PW_IsInDomain
  end type PlaneWall 
    
  type Geometry
    integer:: num_pWall                              ! number of plane walls (total)
    integer:: nPW_local                              ! number of plane walls (local) 
    type(PlaneWall),allocatable,dimension(:):: pWall ! a vector that stores all plane walls 
  contains
    procedure:: InitAllocate  => G_InitAllocate
    procedure:: add_PlaneWall => G_add_PlaneWall
    procedure:: MakeGeometry  => G_MakeGeometry
  end type Geometry
  type(Geometry),public :: DEMGeometry
contains

#ifndef CFDACM
  !**********************************************************************
  ! determining if an sphere (box) has a contact with this wall
  !**********************************************************************
  logical function PW_isInContact(this, box,ovrlp,nv )
    implicit none
    class(PlaneWall)::this
    type(real4),intent(in):: box
    real(RK),intent(out):: ovrlp
    type(real3),intent(out)::nv
        
    !locals
    real(RK)::Radius, dist,t
    type(real3):: p,cp
    
    PW_isInContact = .false.
    Radius = box%w
    p = box

    dist = (this%n .dot. p) + this%d
    if(this%bothside) then
      ovrlp = Radius - abs(dist)
    else
      if( dist < zero ) return
      ovrlp = Radius - dist
    end if
    if(ovrlp<zero) return
        
    IF(this%isInfinite) THEN
      PW_isInContact = .true.
      nv = sign(one, dist) *this%n
      return
    ELSE
      t = -dist
      cp = t * this%n + p
      if(PW_IsInPlane(this%p1,this%p2,this%p3, this%p4,cp)) then
        PW_isInContact = .true.
        nv = sign(one, dist) *this%n
        return
      endif
               
      if(Line_point_check(this%p1,this%p2, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif
      if(Line_point_check(this%p2,this%p3, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif           
      if(Line_point_check(this%p3,this%p4, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif
      if(Line_point_check(this%p4,this%p1, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif           
    ENDIF
  end function PW_isInContact

#else
  !**********************************************************************
  ! determining if an sphere (box) has a contact with this wall
  !**********************************************************************
  logical function PW_isInContact(this, box,ovrlp,nv,clcLubFlag )
    implicit none
    class(PlaneWall)::this
    type(real4),intent(in):: box
    real(RK),intent(out):: ovrlp
    type(real3),intent(out)::nv
    logical,intent(out)::clcLubFlag
        
    !locals
    type(real3):: p,cp
    real(RK)::Radius, dist,t
    
    clcLubFlag= .false.
    PW_isInContact = .false.
    Radius = box%w
    p = box

    dist = (this%n .dot. p) + this%d
    if(this%bothside) then
      ovrlp = Radius - abs(dist)
    else
      if(dist < zero) return
      ovrlp = Radius - dist
    endif
    if(ovrlp<zero) then
      IF(this%isInfinite) THEN
        clcLubFlag= .true.
      ELSE
        t = -dist
        cp= t * this%n + p
        if(PW_IsInPlane(this%p1,this%p2,this%p3, this%p4,cp)) then
          clcLubFlag= .true.
        endif
      ENDIF
      return
    endif
        
    IF(this%isInfinite) THEN
      PW_isInContact = .true.
      nv = sign(one, dist) *this%n
      return
    ELSE
      t  = -dist
      cp = t * this%n + p
      if(PW_IsInPlane(this%p1,this%p2,this%p3, this%p4,cp)) then
        PW_isInContact = .true.
        nv = sign(one, dist) *this%n
        return
      endif
               
      if(Line_point_check(this%p1,this%p2, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif
      if(Line_point_check(this%p2,this%p3, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif           
      if(Line_point_check(this%p3,this%p4, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif
      if(Line_point_check(this%p4,this%p1, box, nv, ovrlp)) then
        PW_isInContact = .true.
        return
      endif           
    ENDIF
  end function PW_isInContact
#endif

  !**********************************************************************
  ! checking if the plane have some overlpa region with the domain
  !**********************************************************************
  function PW_IsInDomain(this,dpmin, dpmax) result(res)
    implicit none
    class(PlaneWall)::this
    type(real3),intent(in):: dpmin, dpmax
    logical:: res
      
    ! locals
    real(RK)::t
    type(real3),dimension(8):: DomPoint
    integer,dimension(8)::intarr
    integer::i,sumintarr
    type(real3)::fp1,fp2,fp3,fp4,sp1,sp2,sp3,sp4
      
    res= .false.
    DomPoint(1)=real3(dpmin%x, dpmin%y, dpmin%z)
    DomPoint(2)=real3(dpmax%x, dpmin%y, dpmin%z)
    DomPoint(3)=real3(dpmax%x, dpmax%y, dpmin%z)
    DomPoint(4)=real3(dpmin%x, dpmax%y, dpmin%z)
    DomPoint(5)=real3(dpmin%x, dpmin%y, dpmax%z)
    DomPoint(6)=real3(dpmax%x, dpmin%y, dpmax%z)
    DomPoint(7)=real3(dpmax%x, dpmax%y, dpmax%z)
    DomPoint(8)=real3(dpmin%x, dpmax%y, dpmax%z)
  
    ! Firstly, check whether there are some points within the domain or not
    if( PointIsInDomain(this%p1, dpmin, dpmax) .or. PointIsInDomain(this%p2, dpmin, dpmax)  .or. &
      PointIsInDomain(this%p3, dpmin, dpmax) .or. PointIsInDomain(this%p4, dpmin, dpmax)) then
      res= .true.
      return
    endif
         
    ! Secondly, check whether the reflections overlap or not
    ! reflection in x-y plane
    fp1= real3(this%p1%x,this%p1%y,zero); fp2= real3(this%p2%x,this%p2%y,zero)
    fp3= real3(this%p3%x,this%p3%y,zero); fp4= real3(this%p4%x,this%p4%y,zero)
    sp1= real3(dpmin%x,  dpmin%y,  zero); sp2= real3(dpmax%x,  dpmin%y,  zero)
    sp3= real3(dpmax%x,  dpmax%y,  zero); sp4= real3(dpmin%x,  dpmax%y,  zero) 
    if(.not.(IsReflectionOvlp(fp1,fp2,fp3,fp4, sp1,sp2,sp3,sp4))) return
      
    ! reflection in x-z plane
    fp1= real3(this%p1%x,zero,this%p1%z); fp2= real3(this%p2%x,zero,this%p2%z)
    fp3= real3(this%p3%x,zero,this%p3%z); fp4= real3(this%p4%x,zero,this%p4%z) 
    sp1= real3(dpmin%x,  zero,dpmin%z);   sp2= real3(dpmin%x,  zero,dpmax%z)
    sp3= real3(dpmax%x,  zero,dpmax%z);   sp4= real3(dpmax%x,  zero,dpmin%z)
    if(.not.(IsReflectionOvlp(fp1,fp2,fp3,fp4, sp1,sp2,sp3,sp4))) return
      
    ! reflection in y-z plane
    fp1= real3(zero,this%p1%y,this%p1%z); fp2= real3(zero,this%p2%y,this%p2%z)
    fp3= real3(zero,this%p3%y,this%p3%z); fp4= real3(zero,this%p4%y,this%p4%z) 
    sp1= real3(zero,dpmin%y,  dpmin%z);   sp2= real3(zero,dpmax%y,  dpmin%z)
    sp3= real3(zero,dpmax%y,  dpmax%z);   sp4= real3(zero,dpmin%y,  dpmax%z)
    if(.not.(IsReflectionOvlp(fp1,fp2,fp3,fp4, sp1,sp2,sp3,sp4))) return 
      
    ! Thirdly, check whether the domain is totally located in one side of the plane 
    intarr=-99
    do i=1,8
      t = ((this%n .dot. DomPoint(i))+ this%d )
      if(t>1.0E-8_RK) then
        intarr(i)=1
      elseif ( t< -1.0E-8_RK ) then
        intarr(i)=0
      endif
    enddo
    sumintarr=sum(intarr)
    if(sumintarr==0 .or. sumintarr==8) return

    res= .true.
  end function PW_IsInDomain
  
  !**********************************************************************
  ! checking if the point lays within the boundaries of the plane
  !**********************************************************************
  logical function PW_IsInPlane( p1,p2,p3,p4, cp )
    implicit none
  type(real3),intent(in)::p1,p2,p3,p4, cp
        
    !// locals
    real(RK):: p1p3, p2p4,p1p4,p2p3,p1p2,p3p4,p2p2
    type(real3) p1p, p2p, p3p, p4p
        
    !// body    
    ! Here (p1,p2,p3,p4) are the four point of the plane.
    PW_IsInPlane = .false.    
    p1p = P1-cp
  p2p = P2-cp
  p3p = P3-cp
    p4p = P4-cp
        
    ! first condition u.w<0
    ! u.w = [(p1-p)x(p2-p)].[(p3-p)x(p4-p)] = (p1p.p3p)(p2p.p4p) - (p1p.p4p)(p2p.p3p)
  p1p3 = p1p .dot. p3p
    p2p4 = p2p .dot. p4p
    p1p4 = p1p .dot. p4p
    p2p3 = p2p .dot. p3p
  if(p1p3*p2p4-p1p4*p2p3<zero) return
  
  ! second condition v.x < 0
    ! v.x = [(p2-p)x(p3-p)].[(p4-p)x(p1-p)] = (p1p.p3p)(p2p.p4p) - (p1p.p2p)(p3p.p4p)
  p1p2 = p1p .dot. p2p
  p3p4 = p3p .dot. p4p
  if(p1p3*p2p4-p1p2*p3p4<zero) return

  ! third condition u.v < 0
  ! u.v = [(p1-p)x(p2-p)].[(p2-p)x(p3-p)] = (p1p.p2p)(p2p.p3p) - (p1p.p3p)(p2p.p2p)
  p2p2 = p2p .dot. p2p
  if(p1p2*p2p3-p1p3*p2p2<zero) return 
  
    PW_IsInPlane = .true.
  end function PW_IsInPlane  
    
  !**********************************************************************
  ! checking if the point lays within the line
  !**********************************************************************
  logical function Line_point_check( lp1,lp2, dpos, nv,ovrlp )
    implicit none
    type(real3),intent(in) :: lp1,lp2
    type(real4),intent(in) :: dpos
    type(real3),intent(out):: nv
    real(RK),intent(out)::ovrlp
    
    real(RK):: t, r,length
    type(real3):: w,v,pos,cp
    
    pos = dpos
    w = pos-lp1 
    v = lp2-lp1
    length = (lp1.dist.lp2)
    r = dpos%w
    t = (w.dot.v )/(v.dot.v)
    
    Line_point_check = .false.
    if( t>= zero .and. t<= one ) then
      cp = (v * t) + lp1
    elseif( t >= (-r/length)  .and. t <zero )then
      cp = lp1
    elseif( t> one  .and. t>= (one+r/length) ) then
      cp = lp2
    else
      Line_point_check = .false.
      return
    endif
    
    ovrlp = r - (cp .dist. pos)
    if( ovrlp >= zero )then
      nv = (pos .nv. cp)   
      Line_point_check = .true.
      return 
    endif
  end function Line_point_check

  !**********************************************************************
  ! checking whether the point lays within the domain
  !**********************************************************************    
  function PointIsInDomain(point, dpmin, dpmax) result(res)
    implicit none
    type(real3),intent(in)::point, dpmin, dpmax
    logical::res
      
    !locals
    real(RK)::realeps=1.000E-10_RK
      
    ! here I slightly expand the domain, to make sure the point  on the domain surface can be regarded as "PointIsInDomain"
    res = .false.
    if(point%x<dpmin%x-realeps .or.  point%x>dpmax%x+realeps) return
    if(point%y<dpmin%y-realeps .or.  point%y>dpmax%y+realeps) return
    if(point%z<dpmin%z-realeps .or.  point%z>dpmax%z+realeps) return
    res = .true.
  end function PointIsInDomain
    
  !**********************************************************************
  ! checking whether the two reflection plane have some overlap region
  !**********************************************************************
  function IsReflectionOvlp(fp1,fp2,fp3,fp4,sp1,sp2,sp3,sp4) result(res)
    implicit none
    type(real3),intent(in)::fp1,fp2,fp3,fp4,sp1,sp2,sp3,sp4
    logical::res
       
    ! Here (fp1,fp2,fp3,fp4) represents first plane
    ! And (sp1,sp2,sp3,sp4) represent  second plane
    res = .true.
       
    ! Firstly, check whether there are some points from one plane located in the inner region of the other plane  
    if(PW_IsInPlane(fp1,fp2,fp3,fp4, sp1)) return
    if(PW_IsInPlane(fp1,fp2,fp3,fp4, sp2)) return
    if(PW_IsInPlane(fp1,fp2,fp3,fp4, sp3)) return
    if(PW_IsInPlane(fp1,fp2,fp3,fp4, sp4)) return
    if(PW_IsInPlane(sp1,sp2,sp3,sp4, fp1)) return
    if(PW_IsInPlane(sp1,sp2,sp3,sp4, fp2)) return
    if(PW_IsInPlane(sp1,sp2,sp3,sp4, fp3)) return
    if(PW_IsInPlane(sp1,sp2,sp3,sp4, fp4)) return
       
    ! secondly, check whether there are some inner insection points among the lines from different plane
    if(IsLineInsertInner(fp1,fp2,sp1,sp2)) return
    if(IsLineInsertInner(fp1,fp2,sp2,sp3)) return
    if(IsLineInsertInner(fp1,fp2,sp3,sp4)) return
    if(IsLineInsertInner(fp1,fp2,sp4,sp1)) return
       
    if(IsLineInsertInner(fp2,fp3,sp1,sp2)) return
    if(IsLineInsertInner(fp2,fp3,sp2,sp3)) return
    if(IsLineInsertInner(fp2,fp3,sp3,sp4)) return
    if(IsLineInsertInner(fp2,fp3,sp4,sp1)) return
     
    if(IsLineInsertInner(fp3,fp4,sp1,sp2)) return
    if(IsLineInsertInner(fp3,fp4,sp2,sp3)) return
    if(IsLineInsertInner(fp3,fp4,sp3,sp4)) return
    if(IsLineInsertInner(fp3,fp4,sp4,sp1)) return

    if(IsLineInsertInner(fp4,fp1,sp1,sp2)) return
    if(IsLineInsertInner(fp4,fp1,sp2,sp3)) return
    if(IsLineInsertInner(fp4,fp1,sp3,sp4)) return
    if(IsLineInsertInner(fp4,fp1,sp4,sp1)) return
    res = .false.
        
  end function IsReflectionOvlp
    
  !**********************************************************************
  ! checking whether the two line have inner insertion
  !**********************************************************************    
  function IsLineInsertInner(L1_p1, L1_p2, L2_p1, L2_p2) result(res)
    implicit none
    type(real3),intent(in)::L1_p1, L1_p2, L2_p1, L2_p2
    logical :: res
      
    ! locals
    real(RK)::a,b,c,d,e,dist,t1,t2
    type(real3):: v1,v2,L21_p1,v1_crs_v2, L21_crs_v1,L21_crs_v2,cp
      
    ! Here I use L1_p1 and L1_p2 to express the starting and ending points of the line L1
    !            L2_p1 and L2_p2 to express the starting and ending points of the line L2
    ! And
    !       direction vector of line L1 : v1 = L1_p2 - L1_p1
    !       direction vector of line L2 : v2 = L2_p2 - L2_p1
    ! So any point on line L1 and L2 can be determined by the following parameter equation:
    !       line L1:   P1 = L1_p1 + a * v1, where 'a' is the parameter
    !       line L2:   P2 = L2_p1 + b * v2, where 'b' is the parameter
    ! The inner insert point of line L1 and line L2 should satisfy the following expression:
    !       P_insertion = P1 = P2
    ! i.e.
    !       L1_p1 + a * v1  = L2_p1 + b * v2                                 (1)
    ! We use the cross production of v1 to both sides of Eq.(1):
    !       (L1_p1 + a * v1) x v1 = (L2_p1 + b * v2) x v1
    ! i.e.
    !       (L1_p1 - L2_p1) x v1 =  b* v2 x v1                               (2)    
    ! S
    !       b = sign_b*norm((L2_p1 - L1_p1) x v1 ) /norm( v1 x v2 )          (3)
    ! Where
    !       sign_b= 1,  if ((L2_p1 - L1_p1) x v1 ) .dot. ( v1 x v2 ) )>0     (4)
    !             =-1,  if ((L2_p1 - L1_p1) x v1 ) .dot. ( v1 x v2 ) )<0
    ! Similarly, we can get  that:
    !       a = sign_a*norm((L2_p1 - L1_p1) x v2 ) / norm( v1 x v2 )         (5)
    !       sign_a= 1,  if ((L2_p1 - L1_p1) x v2 ) .dot. ( v1 x v2 ) )>0
    !             =-1,  if ((L2_p1 - L1_p1) x v2 ) .dot. ( v1 x v2 ) )<0
    ! If    0 =<a <= 1, and 0 =<b <= 1, L1 and L2 have inner insection point.
    ! In this function, I use the following notes:
    !     L21_p1 = L2_p1 - L1_p1
    !     v1_crs_v2 = v1 x v2
    !     L21_crs_v1 = (L2_p1 - L1_p1) x v1
    !     L21_crs_v2 = (L2_p1 - L1_p1) x v2
    !     c = norm( (L2_p1 - L1_p1) x v1 ) =norm(L21_crs_v1)
    !     d = norm( (L2_p1 - L1_p1) x v2 ) =norm(L21_crs_v2)
    !     e = norm( v1 x v2 ) =norm(v1_crs_v2)
      
    v1 = L1_p2 - L1_p1
    v2 = L2_p2 - L2_p1
    L21_p1 = L2_p1 - L1_p1
      
    v1_crs_v2 = v1 .cross. v2
    L21_crs_v1= L21_p1 .cross. v1
    L21_crs_v2= L21_p1 .cross. v2
      
    c= norm(L21_crs_v1)
    d= norm(L21_crs_v2)
    e= norm(v1_crs_v2)
      
    ! L1 and L2 are parallel or they two are on the same line.
    IF(e<1.0E-10_RK) THEN
      t1 = (L21_p1 .dot. v1 )/ (v1 .dot. v1)
      cp = (v1 * t1) + L1_p1
      dist = cp .dist. L2_p1  ! the distance between two lines
      if(dist >1.0E-6_RK) then
        res = .false.
      else
        t2 = ((L2_p2-L1_p1).dot.v1)/(v1.dot.v1)
        if((t1<zero.and.t2<zero) .or. (t1>one.and.t2>one) ) then
          res = .false.
        else
          res = .true.
        endif
      endif
      return
    ENDIF
      
    if((L21_crs_v2 .dot. v1_crs_v2)>zero) then
      a = d/e
    else
      a= - d/e
    endif
    if((L21_crs_v1 .dot. v1_crs_v2)>zero) then
      b = c/e
    else
      b= - c/e
    endif      
    if(a>=zero-1.00E-10_RK .and. a<=one+1.00E-10_RK   .and. &
       b>=zero-1.00E-10_RK .and. b<=one+1.00E-10_RK)  then
      res = .true.
    else   
      res = .false.
    endif
  end function IsLineInsertInner
  
  !**********************************************************************
  ! Initializing the geometry object 
  !**********************************************************************
  subroutine G_InitAllocate(this)
    implicit none
    class(Geometry):: this
    integer:: nUnitFile,ierror
    character(128) :: chFile
        
    ! locals
    real(RK):: LenExp
         
    LenExp=  1.2_RK*maxval( DEMProperty%Prtcl_PureProp%Radius )
#ifdef CFDACM
    LenExp= LenExp+maxval(dlub_pw)
#endif
    pmin = DEm_opt%SimDomain_min - LenExp * one_r3
    pmax = DEm_opt%SimDomain_max + LenExp * one_r3
    pmin_local =real3(DEM_decomp%xSt,DEM_decomp%ySt,DEM_decomp%zSt)-LenExp*one_r3
    pmax_local =real3(DEM_decomp%xEd,DEM_decomp%yEd,DEM_decomp%zEd)+LenExp*one_r3
        
    this%num_pWall = 0
    this%nPW_local = 0
    allocate(this%pWall( MaxWallSize ))
        
    if(nrank/=0) return
    write(chFile,"(A)") trim(DEM_opt%ResultsDir)//"WallsFor"//trim(DEM_opt%RunName)//".backup"
    open(newunit=nUnitFile, file=chfile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0.and.nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"G_InitAllocate","Cannot open file: "//trim(chFile))
    close(nUnitFile,IOSTAT=ierror)
        
  end subroutine G_InitAllocate

  !**********************************************************************
  ! Add a plane to the geometry object 
  !**********************************************************************
  subroutine G_add_PlaneWall(this, p1, p2, p3, p4 , user_id ,prop_type,  both, infinite, t_vel_t)
    implicit none
    class(Geometry) this
    type(real3),intent(in) :: p1, p2, p3, p4   ! corner points 
    integer,intent(in) :: user_id, prop_type
    logical,optional,intent(in) :: both        !  both side active status
    logical,optional,intent(in) :: infinite
    type(real3),optional,intent(in) :: t_vel_t ! translational velocity
        
    ! locals
    logical ::lboth,linfinite
    type(PlaneWall)::wall
    type(PlaneWall),dimension(:),allocatable:: wall_temp
    type(real3)::t_vel, ln
    integer:: nUnitFile,ierror
    character(128) :: chFile        
 
    IF(this%nPW_local == MaxWallSize) THEN
      MaxWallSize  = int( real(MaxWallSize, RK) *1.2_RK) +1
      call move_alloc(this%pWall, wall_temp)
      allocate(this%pWall(MaxWallSize),Stat=ierror) 
      if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Reallocation failed, 2 ")
      this%pWall(1: this%nPW_local) = wall_temp
      deallocate(wall_temp)
    ENDIF        
        
    ! checking for both side spec. 
    lboth = .false.   !.false.
    linfinite=.true.
    t_vel=zero_r3
    if(present(both))    lboth = both
    if(present(infinite))linfinite=infinite
    if(present(t_vel_t)) t_vel = t_vel_t
    
    ! assignments 
    ln = (p2-p1).cross.(p3-p1)
    wall%p1 = p1
    wall%p2 = p2
    wall%p3 = p3
    wall%p4 = p4
    wall%min_point= min(p4,min(p3,min(p2,p1)))
    wall%max_point= max(p4,max(p3,max(p2,p1)))
    wall%user_id = user_id
    wall%wall_Type = prop_type
    wall%bothSide = lboth
    wall%IsInfinite=linfinite
    wall%n = ln/norm(ln)
    wall%d = - (wall%n .dot. p1)
    wall%trans_vel = t_vel
    if(abs((wall%n .dot. p4) +wall%d) >= 0.00001_RK .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Cannot create a plane wall, wall No. "//num2str(wall%wall_id))
    endif
        
    IF(wall%IsInDomain(pmin, pmax)) THEN
      this%num_pWall= this%num_pWall +1
      wall%wall_id  = this%num_pWall + DEM_opt%Base_wall_id
      if(wall%IsInDomain(pmin_local,pmax_local)) then
        this%nPW_local = this%nPW_local +1
        this%pWall(this%nPW_local) = wall 
      endif
        
      if(nrank/=0) return
      write(chFile,"(A)") trim(DEM_opt%ResultsDir)//"WallsFor"//trim(DEM_opt%RunName)//".backup"
      open(newunit=nUnitFile, file=chFile, status='old',position='append',form='formatted',IOSTAT=ierror )
      if(ierror/=0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Cannot open file: "//trim(chFile))
      write(nUnitFile,* ) p1
      write(nUnitFile,* ) p2
      write(nUnitFile,* ) p3
      write(nUnitFile,* ) p4
      close(nUnitFile,IOSTAT=ierror)
    ELSE
      if(nrank/=0) return
      call DEMLogInfo%CheckForError(ErrT_Pass,"G_add_PlaneWall","  The following plane ISNOT within the simulation domain: ")
      call DEMLogInfo%OutInfo("   It will be skipped :",3, .true.)
      call DEMLogInfo%OutInfo("   Point 1: "//trim(num2str(p1%x))//'  '//trim(num2str(p1%y))//'  '//trim(num2str(p1%z)), 3, .true.)
      call DEMLogInfo%OutInfo("   Point 2: "//trim(num2str(p2%x))//'  '//trim(num2str(p2%y))//'  '//trim(num2str(p2%z)), 3, .true.)
      call DEMLogInfo%OutInfo("   Point 3: "//trim(num2str(p3%x))//'  '//trim(num2str(p3%y))//'  '//trim(num2str(p3%z)), 3, .true.)
      call DEMLogInfo%OutInfo("   Point 4: "//trim(num2str(p4%x))//'  '//trim(num2str(p4%y))//'  '//trim(num2str(p4%z)), 3, .true.)
    ENDIF
  end subroutine G_add_PlaneWall

  !**********************************************************************
  ! MakeGeometry
  !**********************************************************************     
  subroutine G_MakeGeometry(this,chFile)
    implicit none
    class(Geometry)::this
    character(*),intent(in)::chFile
        
    !locals
    integer:: i,nplane,nUnitFile,ierror
    integer,allocatable,dimension(:):: user_id, wall_Type
    logical,allocatable,dimension(:):: BothSide, IsInfinite
    type(real3):: p01,p02,p03,p04,p05,p06,p07,p08,p09,p10,p11,p12
    type(real3):: p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24
    type(real3),allocatable,dimension(:)::Point1,Point2,Point3,Point4,TraVel
    NAMELIST/GeometryMakingNumPlane/nplane
    NAMELIST/GeometryMakingParam/Point1,Point2,Point3,Point4,TraVel,user_id, wall_Type,bothSide, IsInfinite
        
    call this%InitAllocate()
    if(DEM_Opt%GeometrySource ==0) then  ! add the geometry directly

      !sandbox
      p01= real3(-0.075_RK,  0.30_RK, -0.075_RK)
      p02= real3( 0.075_RK,  0.30_RK, -0.075_RK)
      p03= real3( 0.075_RK,  0.30_RK,  0.075_RK)
      p04= real3(-0.075_RK,  0.30_RK,  0.075_RK)

      p05= real3(-0.075_RK,  0.14_RK, -0.075_RK)
      p06= real3( 0.075_RK,  0.14_RK, -0.075_RK)
      p07= real3( 0.075_RK,  0.14_RK,  0.075_RK)
      p08= real3(-0.075_RK,  0.14_RK,  0.075_RK)
  
      p09= real3(-0.015_RK,  0.05_RK, -0.015_RK)
      p10= real3( 0.015_RK,  0.05_RK, -0.015_RK)
      p11= real3( 0.015_RK,  0.05_RK,  0.015_RK)
      p12= real3(-0.015_RK,  0.05_RK,  0.015_RK)

      p13= real3(-0.015_RK, -0.05_RK, -0.015_RK)
      p14= real3( 0.015_RK, -0.05_RK, -0.015_RK)
      p15= real3( 0.015_RK, -0.05_RK,  0.015_RK)
      p16= real3(-0.015_RK, -0.05_RK,  0.015_RK)

      p17= real3(-0.075_RK, -0.14_RK, -0.075_RK)
      p18= real3( 0.075_RK, -0.14_RK, -0.075_RK)
      p19= real3( 0.075_RK, -0.14_RK,  0.075_RK)
      p20= real3(-0.075_RK, -0.14_RK,  0.075_RK)

      p21= real3(-0.075_RK, -0.30_RK, -0.075_RK)
      p22= real3( 0.075_RK, -0.30_RK, -0.075_RK)
      p23= real3( 0.075_RK, -0.30_RK,  0.075_RK)
      p24= real3(-0.075_RK, -0.30_RK,  0.075_RK)            
      call this%add_PlaneWall( p01, p05, p06, p02, 1, 1, infinite=.false. ) !
      call this%add_PlaneWall( p05, p09, p10, p06, 1, 1, infinite=.false. ) !   
      call this%add_PlaneWall( p09, p13, p14, p10, 1, 1, infinite=.false. ) 
      call this%add_PlaneWall( p13, p17, p18, p14, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p17, p21, p22, p18, 1, 1, infinite=.false. )

      call this%add_PlaneWall( p01, p04, p08, p05, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p05, p08, p12, p09, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p09, p12, p16, p13, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p13, p16, p20, p17, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p17, p20, p24, p21, 1, 1, infinite=.false. )

      call this%add_PlaneWall( p03, p07, p08, p04, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p07, p11, p12, p08, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p11, p15, p16, p12, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p15, p19, p20, p16, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p19, p23, p24, p20, 1, 1, infinite=.false. )

      call this%add_PlaneWall( p02, p06, p07, p03, 1, 1, infinite=.false. )    
      call this%add_PlaneWall( p06, p10, p11, p07, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p10, p14, p15, p11, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p14, p18, p19, p15, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p18, p22, p23, p19, 1, 1, infinite=.false. )

      call this%add_PlaneWall( p01, p02, p03, p04, 1, 1, infinite=.false. )
      call this%add_PlaneWall( p21, p24, p23, p22, 1, 1, infinite=.false. )
          
    elseif(DEM_Opt%GeometrySource ==1) then  ! add the geometry from the NAMELIST "&GeometryMakingParam"
                
      open(newunit=nUnitFile, file=chFile,status='old', form='formatted', IOSTAT=ierror)
      if(ierror /= 0 .and. nrank==0) call DEMLogInfo%CheckForError(ErrT_Abort,"G_MakeGeometry","Cannot open file:"//trim(chFile))
      read(nUnitFile, nml=GeometryMakingNumPlane)
      if(nplane<1) then
        close(nUnitFile, IOSTAT=ierror);return
      endif
      allocate(user_id(nplane), wall_Type(nplane),bothSide(nplane),IsInfinite(nplane))
      allocate(Point1(nplane),Point2(nplane),Point3(nplane),Point4(nplane),TraVel(nplane))
      rewind(nUnitFile)
      read(nUnitFile, nml=GeometryMakingParam)
      close(nUnitFile, IOSTAT=ierror)
      do i=1,nplane
        call this%add_PlaneWall(Point1(i),Point2(i),Point3(i),Point4(i),user_id(i),wall_Type(i),bothSide(i),IsInfinite(i),TraVel(i))
      enddo
    elseif(DEM_Opt%GeometrySource ==2) then  ! add the geometry from the external STL file
            
    endif
  end subroutine G_MakeGeometry
 
end module prtcl_Geometry
