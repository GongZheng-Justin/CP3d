module LPT_Geometry
  use m_TypeDef
  use m_LogInfo
  use LPT_Property
  use LPT_decomp_2d
  use LPT_Parameters
  use m_Decomp2d,only:nrank
  use m_Parameters,only:xlx,yly,zlz
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
    real(RK):: d          ! d in the implicit equation: ax+by+cz+d = 0   
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
    procedure:: MakeGeometry =>G_MakeGeometry
  end type Geometry
    
  type(Geometry),public :: LPTGeometry
    
contains
    
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
      if( dist < 0.0_RK ) return
      ovrlp = Radius - dist
    end if
    if(ovrlp<0.0_RK) return
        
    IF(this%isInfinite) THEN
      PW_isInContact = .true.
      nv = sign(1.0_RK, dist) *this%n
      return
    ELSE
      t = -((this%n .dot. p)+ this%d )
      cp = t * this%n + p
      if(PW_IsInPlane(this%p1,this%p2,this%p3, this%p4,cp)) then
        PW_isInContact = .true.
        nv = sign(1.0_RK, dist) *this%n
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
    fp1= real3(this%p1%x,this%p1%y,0.0_RK); fp2= real3(this%p2%x,this%p2%y,0.0_RK)
    fp3= real3(this%p3%x,this%p3%y,0.0_RK); fp4= real3(this%p4%x,this%p4%y,0.0_RK)
    sp1= real3(dpmin%x,  dpmin%y,  0.0_RK); sp2= real3(dpmax%x,  dpmin%y,  0.0_RK)
    sp3= real3(dpmax%x,  dpmax%y,  0.0_RK); sp4= real3(dpmin%x,  dpmax%y,  0.0_RK) 
    if(.not.(IsReflectionOvlp(fp1,fp2,fp3,fp4, sp1,sp2,sp3,sp4))) return
      
    ! reflection in x-z plane
    fp1= real3(this%p1%x,0.0_RK,this%p1%z); fp2= real3(this%p2%x,0.0_RK,this%p2%z)
    fp3= real3(this%p3%x,0.0_RK,this%p3%z); fp4= real3(this%p4%x,0.0_RK,this%p4%z) 
    sp1= real3(dpmin%x,  0.0_RK,dpmin%z);   sp2= real3(dpmin%x,  0.0_RK,dpmax%z)
    sp3= real3(dpmax%x,  0.0_RK,dpmax%z);   sp4= real3(dpmax%x,  0.0_RK,dpmin%z)
    if(.not.(IsReflectionOvlp(fp1,fp2,fp3,fp4, sp1,sp2,sp3,sp4))) return
      
    ! reflection in y-z plane
    fp1= real3(0.0_RK,this%p1%y,this%p1%z); fp2= real3(0.0_RK,this%p2%y,this%p2%z)
    fp3= real3(0.0_RK,this%p3%y,this%p3%z); fp4= real3(0.0_RK,this%p4%y,this%p4%z) 
    sp1= real3(0.0_RK,dpmin%y,  dpmin%z);   sp2= real3(0.0_RK,dpmax%y,  dpmin%z)
    sp3= real3(0.0_RK,dpmax%y,  dpmax%z);   sp4= real3(0.0_RK,dpmin%y,  dpmax%z)
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
    if(p1p3*p2p4-p1p4*p2p3<0.0_RK) return
  
    ! second condition v.x < 0
    ! v.x = [(p2-p)x(p3-p)].[(p4-p)x(p1-p)] = (p1p.p3p)(p2p.p4p) - (p1p.p2p)(p3p.p4p)
    p1p2 = p1p .dot. p2p
    p3p4 = p3p .dot. p4p
    if(p1p3*p2p4-p1p2*p3p4<0.0_RK) return

    ! third condition u.v < 0
    ! u.v = [(p1-p)x(p2-p)].[(p2-p)x(p3-p)] = (p1p.p2p)(p2p.p3p) - (p1p.p3p)(p2p.p2p)
    p2p2 = p2p .dot. p2p
    if(p1p2*p2p3-p1p3*p2p2<0.0_RK) return 
  
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
    if( t>= 0.0_RK .and. t<= 1.0_RK ) then
      cp = (v * t) + lp1
    elseif( t >= (-r/length)  .and. t <0.0_RK )then
      cp = lp1
    elseif( t> 1.0_RK  .and. t>= (1.0_RK+r/length) ) then
      cp = lp2
    else
      Line_point_check = .false.
      return
    endif
    
    ovrlp = r - (cp .dist. pos)
    if( ovrlp >= 0.0_RK )then
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
        if((t1<0.0_RK.and.t2<0.0_RK) .or. (t1>1.0_RK .and. t2>1.0_RK) ) then
          res = .false.
        else
          res = .true.
        endif
      endif
      return
    ENDIF
      
    if((L21_crs_v2 .dot. v1_crs_v2)>0.0_RK) then
      a = d/e
    else
      a= - d/e
    endif
    if((L21_crs_v1 .dot. v1_crs_v2)>0.0_RK) then
      b = c/e
    else
      b= - c/e
    endif      
    if(a>=0.0_RK-1.00E-10_RK .and. a<=1.0_RK+1.00E-10_RK   .and. &
       b>=0.0_RK-1.00E-10_RK .and. b<=1.0_RK+1.00E-10_RK)  then
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
         
    LenExp=  1.2_RK*maxval( LPTProperty%Prtcl_PureProp%Radius )
    pmin = LPT_opt%SimDomain_min - LenExp *real3(1.0_RK,1.0_RK,1.0_RK)
    pmax = LPT_opt%SimDomain_max + LenExp *real3(1.0_RK,1.0_RK,1.0_RK)
    pmin_local =real3(LPT_decomp%xSt,LPT_decomp%ySt,LPT_decomp%zSt)-LenExp*real3(1.0_RK,1.0_RK,1.0_RK)
    pmax_local =real3(LPT_decomp%xEd,LPT_decomp%yEd,LPT_decomp%zEd)+LenExp*real3(1.0_RK,1.0_RK,1.0_RK)
        
    this%num_pWall = 0
    this%nPW_local = 0
    allocate(this%pWall( MaxWallSize ))
        
    if(nrank/=0) return
    write(chFile,"(A)") trim(LPT_opt%ResultsDir)//"WallsFor"//trim(LPT_opt%RunName)//".backup"
    open(newunit=nUnitFile, file=chfile,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0 .and. nrank==0) call LPTLogInfo%CheckForError(ErrT_Abort,"G_InitAllocate","Cannot open file: "//trim(chFile))
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
      if(ierror/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Reallocation failed, 2 ")
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
    wall%min_point%x= min(p4%x,min(p3%x,min(p2%x,p1%x)))
    wall%min_point%y= min(p4%y,min(p3%y,min(p2%y,p1%y)))
    wall%min_point%z= min(p4%z,min(p3%z,min(p2%z,p1%z)))
    wall%max_point%x= max(p4%x,max(p3%x,max(p2%x,p1%x)))
    wall%max_point%y= max(p4%y,max(p3%y,max(p2%y,p1%y)))
    wall%max_point%z= max(p4%z,max(p3%z,max(p2%z,p1%z)))
    wall%user_id = user_id
    wall%wall_Type = prop_type
    wall%bothSide = lboth
    wall%IsInfinite=linfinite
    wall%n = ln/norm(ln)
    wall%d = - (wall%n .dot. p1)
    wall%trans_vel = t_vel
    if(abs((wall%n .dot. p4) +wall%d) >= 0.00001_RK .and. nrank==0) then
      call LPTLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Cannot create a plane wall, wall No. "//num2str(wall%wall_id))
    endif
        
    IF(wall%IsInDomain(pmin, pmax)) THEN
      this%num_pWall= this%num_pWall +1
      wall%wall_id  = this%num_pWall
      if(wall%IsInDomain(pmin_local,pmax_local)) then
        this%nPW_local = this%nPW_local +1
        this%pWall(this%nPW_local) = wall 
      endif
        
      if(nrank/=0) return
      write(chFile,"(A)") trim(LPT_opt%ResultsDir)//"WallsFor"//trim(LPT_opt%RunName)//".backup"
      open(newunit=nUnitFile, file=chFile, status='old',position='append',form='formatted',IOSTAT=ierror )
      if(ierror/=0 .and. nrank==0) call LPTLogInfo%CheckForError(ErrT_Abort,"G_add_PlaneWall","Cannot open file: "//trim(chFile))
      write(nUnitFile,* ) p1
      write(nUnitFile,* ) p2
      write(nUnitFile,* ) p3
      write(nUnitFile,* ) p4
      close(nUnitFile,IOSTAT=ierror)
    ELSE
      if(nrank/=0) return
      call LPTLogInfo%CheckForError(ErrT_Pass,"G_add_PlaneWall","  The following plane ISNOT within the simulation domain: ")
      call LPTLogInfo%OutInfo("   It will be skipped :",3, .true.)
      call LPTLogInfo%OutInfo("   Point 1: "//trim(num2str(p1%x))//'  '//trim(num2str(p1%y))//'  '//trim(num2str(p1%z)), 3, .true.)
      call LPTLogInfo%OutInfo("   Point 2: "//trim(num2str(p2%x))//'  '//trim(num2str(p2%y))//'  '//trim(num2str(p2%z)), 3, .true.)
      call LPTLogInfo%OutInfo("   Point 3: "//trim(num2str(p3%x))//'  '//trim(num2str(p3%y))//'  '//trim(num2str(p3%z)), 3, .true.)
      call LPTLogInfo%OutInfo("   Point 4: "//trim(num2str(p4%x))//'  '//trim(num2str(p4%y))//'  '//trim(num2str(p4%z)), 3, .true.)
    ENDIF
  end subroutine G_add_PlaneWall

  !**********************************************************************
  ! MakeGeometry
  !**********************************************************************     
  subroutine G_MakeGeometry(this)
    implicit none
    class(Geometry)::this
        
    !locals
    type(real3):: p01,p02,p03,p04
        
    call this%InitAllocate()
    p01= real3( 0.0_RK,  0.0_RK,  0.0_RK)
    p02= real3( 0.0_RK,  0.0_RK,   zlz)
    p03= real3(  xlx,  0.0_RK,   zlz)
    p04= real3(  xlx,  0.0_RK,  0.0_RK)
    call this%add_PlaneWall( p01, p02, p03, p04, 1, 1, infinite=.true. )
    p01= real3( 0.0_RK,  yly,  0.0_RK)
    p02= real3(  xlx,  yly,  0.0_RK)
    p03= real3(  xlx,  yly,   zlz)
    p04= real3( 0.0_RK,  yly,   zlz)
    call this%add_PlaneWall( p01, p02, p03, p04, 1, 1, infinite=.true. )
  end subroutine G_MakeGeometry
 
end module LPT_Geometry
