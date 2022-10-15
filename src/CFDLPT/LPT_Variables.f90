module LPT_Variables
  use MPI
  use LPT_LogInfo
  use LPT_TypeDef
  use LPT_Property
  use LPT_decomp_2d
  use LPT_Parameters
  use m_Decomp2d,only: nrank
  implicit none
  private

  integer,dimension(:),allocatable,public:: GPrtcl_id
  integer,dimension(:),allocatable,public:: GPrtcl_pType
  integer,dimension(:),allocatable,public:: GPrtcl_usrMark
  type(real4),dimension(:),allocatable,public::   GPrtcl_PosR
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linVel
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linAcc
  type(real3),dimension(:),allocatable,public:: GPrtcl_PosOld
  type(real3),dimension(:),allocatable,public:: GPrtcl_VFluid
  
  type(real3),dimension(:),  allocatable,public:: GPrtcl_FpForce

  type VarList
    integer:: nlocal     ! number of particles in local processor
    integer:: mlocal     ! the possible maxium # of particles in local processor
    integer:: tsize      ! size for translational veloctity and acceleration
  contains
    procedure:: AllocateAllVar     => GL_AllocateAllVar
    procedure:: MakingAllPrtcl     => GL_MakingAllPrtcl
    procedure:: ReallocatePrtclVar => GL_ReallocatePrtclVar
    procedure:: copy               => GL_copy
  end type VarList
  type(VarList),public::GPrtcl_list
    
contains

  !**********************************************************************
  ! GL_AllocateAllVar
  !**********************************************************************
  subroutine GL_AllocateAllVar(this)
    implicit none
    class(VarList)::this

    ! locals
    type(real3)::SimLen
    integer::mlocal,numPrtcl,ierrTmp,ierror=0
    real(RK)::xst,xed,yst,yed,zst,zed,vol_tot,vol_local

    numPrtcl    = LPT_opt%numPrtcl

    ! step1: allocating memory for particles
    if (LPT_opt%PI_Method==PIM_AB2) then
      this%tsize=2
    elseif(LPT_opt%PI_Method==PIM_AB3) then
      this%tsize=3
    endif

    ! step0: determine initial msize
    xst=LPT_decomp%xSt; xed=LPT_decomp%xEd
    yst=LPT_decomp%ySt; yed=LPT_decomp%yEd
    zst=LPT_decomp%zSt; zed=LPT_decomp%zEd

    SimLen = LPT_Opt%SimDomain_max - LPT_Opt%SimDomain_min
    vol_tot= SimLen%x * SimLen%y * SimLen%z
    
    vol_local =(xed-xst)*(yed-yst)*(zed-zst)
    mlocal = int(vol_local/vol_tot*real(numPrtcl,kind=RK))
    mlocal = int(1.5_RK*mlocal)
    mlocal = min(mlocal, numPrtcl)
    mlocal = max(mlocal, 10)
        
    allocate(GPrtcl_id(mlocal),                Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_pType(mlocal),             Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_usrMark(mlocal),           Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_PosR(mlocal),              Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_linVel(this%tsize,mlocal), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_linAcc(this%tsize,mlocal), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_PosOld(mlocal),            Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_VFluid(mlocal),            Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar","Allocation failed 1")

    allocate(GPrtcl_FpForce(mlocal),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation failed 2")    
    
    GPrtcl_id = 0
    GPrtcl_pType = 1
    GPrtcl_usrMark = 1
    GPrtcl_PosR = zero_r4
    GPrtcl_linVel  = zero_r3
    GPrtcl_linAcc  = zero_r3
    GPrtcl_PosOld  = zero_r3
    GPrtcl_VFluid  = zero_r3
    GPrtcl_FpForce = zero_r3
    
    ! step2: Initialize this%mlocal,this%nlocal
    this%nlocal  = 0
    this%mlocal  = mlocal
  end subroutine GL_AllocateAllVar

  !**********************************************************************
  ! GL_MakingAllPrtcl
  !**********************************************************************
  subroutine GL_MakingAllPrtcl(this, chFile)
    implicit none
    class(VarList)::this
    character(*),intent(in)::chFile
        
    ! locals     
    logical::IsRandomDist
    character,dimension(3):: Fill_order
    integer,dimension(:),allocatable:: sum_bin    
    integer::i,j,k,code,bin_id,numPrtcl,nUnitFile,ierror,nlocal,nlocal_sum
    real(RK)::xst,xed,yst,yed,zst,zed,maxRad,randt,Distance_Ratio,MkPrtclMinpoint(3), MkPrtclMaxpoint(3),realt(3)
    type(real3):: MkPrtclMinpoint_real3, MkPrtclMaxpoint_real3,cntr,l1_vec,l2_vec,l3_vec,lmin_p,lmax_p,dx,pos1,pos2
    NAMELIST /ParticleMakingOption/  MkPrtclMinpoint,MkPrtclMaxpoint, Fill_order, Distance_Ratio,  IsRandomDist
        
    open(newunit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl", "Cannot open file: "//trim(chFile))
    read(nUnitFile, nml=ParticleMakingOption)
    if(nrank==0)write(LPTLogInfo%nUnit, nml=ParticleMakingOption)
    close(nUnitFile,IOSTAT=ierror)
    MkPrtclMinpoint_real3 = MkPrtclMinpoint
    MkPrtclMaxpoint_real3 = MkPrtclMaxpoint
    numPrtcl = LPT_opt%numPrtcl
        
    xst=LPT_decomp%xSt; xed=LPT_decomp%xEd
    yst=LPT_decomp%ySt; yed=LPT_decomp%yEd
    zst=LPT_decomp%zSt; zed=LPT_decomp%zEd

    ! step3: Assign Prtcl_Type, Prtcl_id and Prtcl_PosR
    call Fill_Vectors(Fill_order, l1_vec, l2_vec, l3_vec)
    call system_clock(count=code)
    if(IsRandomDist) code=0
    call random_seed(size = j)
    call random_seed(put = code+63946*(/ (i - 1, i = 1, j) /))
    if(IsRandomDist) Distance_Ratio=one
    maxRad = maxval( LPTProperty%Prtcl_PureProp%Radius )
    lmin_p = MkPrtclMinpoint_real3 + Distance_Ratio* maxRad*one_r3
    lmax_p = MkPrtclMaxpoint_real3 - Distance_Ratio* maxRad*one_r3
    dx = Distance_Ratio *two*maxRad* one_r3
        
    allocate(sum_bin(LPT_opt%numPrtcl_Type))
    sum_bin(1)=LPTProperty%nPrtcl_in_Bin(1)
    do j=2, LPT_opt%numPrtcl_Type
      sum_bin(j)=sum_bin(j-1)+LPTProperty%nPrtcl_in_Bin(j)
    enddo

    nlocal=0; bin_id=-1
    IF(IsRandomDist) THEN
      do i=1,numPrtcl
        do j=1,LPT_opt%numPrtcl_Type
          if(i<=sum_bin(j)) then
            bin_id = j; exit
          endif
        enddo
                
        call random_number(realt)
        cntr%x=  lmin_p%x +realt(1)*(lmax_p%x-lmin_p%x)
        cntr%y=  lmin_p%y +realt(2)*(lmax_p%y-lmin_p%y)
        cntr%z=  lmin_p%z +realt(3)*(lmax_p%z-lmin_p%z)         
        if(cntr%x>=xst.and.cntr%x<xed.and.cntr%y>=yst.and.cntr%y<yed.and.cntr%z>=zst.and.cntr%z<zed) then
          if(nlocal>=this%mlocal) call this%ReallocatePrtclVar(nlocal)
          nlocal=nlocal+1
          GPrtcl_id(nlocal) = i
          GPrtcl_pType(nlocal)= bin_id
          GPrtcl_PosR(nlocal) = cntr
          GPrtcl_PosR(nlocal)%w = LPTProperty%Prtcl_PureProp(bin_id)%Radius
        endif
      enddo    
    ELSE
      cntr = lmin_p       ! start point
      do i=1,numPrtcl
        do j=1,LPT_opt%numPrtcl_Type
          if(i<=sum_bin(j)) then
            bin_id = j; exit
          endif
        enddo

        if(cntr%x>=xst.and.cntr%x<xed.and.cntr%y>=yst.and.cntr%y<yed.and.cntr%z>=zst.and.cntr%z<zed) then
          if(nlocal>=this%mlocal) call this%ReallocatePrtclVar(nlocal)
          nlocal=nlocal+1
          GPrtcl_id(nlocal) = i
          GPrtcl_pType(nlocal)= bin_id
          GPrtcl_PosR(nlocal) = cntr
          GPrtcl_PosR(nlocal)%w = LPTProperty%Prtcl_PureProp(bin_id)%Radius
        endif

        cntr = cntr + l1_vec * dx
        if((l1_vec.dot.cntr)>=(l1_vec.dot.lmax_p))then
          cntr = (lmin_p*l1_vec) + ((cntr+dx)*l2_vec)+(cntr*l3_vec)
          if((l2_vec.dot.cntr)>=(l2_vec.dot.lmax_p))then
            cntr = (cntr*l1_vec)+(lmin_p*l2_vec)+((cntr+dx)*l3_vec)
            if((l3_vec.dot.cntr) >= (l3_vec.dot.lmax_p) .and. nrank==0) then
              call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Not enough space for positioning" ) 
            endif
          endif
        endif
      enddo
    ENDIF
    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    if(nlocal_sum/= numPrtcl .and. nrank==0) then
      call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: "," nlocal_sum/= numPrtcl " )
    endif
    this%nlocal = nlocal
    deallocate(sum_bin)

    ! randomize the GPrtcl_PosR
    do k=1,3
      do i=1,this%nlocal
        call random_number(randt)
        j=int(randt*this%nlocal)+1
        pos1=GPrtcl_PosR(i)
        pos2=GPrtcl_PosR(j)
        GPrtcl_PosR(i)=pos2
        GPrtcl_PosR(j)=pos1
      enddo
    enddo
  end subroutine GL_MakingAllPrtcl
  
  !**********************************************************************
  ! Fill order Vector 
  !**********************************************************************
  subroutine Fill_Vectors( fill_order , l1 , l2, l3 )
    implicit none
    character,dimension(3) :: fill_order
    type(real3),intent(out):: l1, l2, l3
    
    if(fill_order(1)=="x" .or. fill_order(1)=="X") then
      l1 = real3(one,zero,zero)
    elseif(fill_order(1)=="y" .or. fill_order(1)=="Y") then
      l1 = real3(zero,one,zero)
    elseif(fill_order(1)=="z" .or. fill_order(1)=="Z") then
      l1 = real3(zero,zero,one)
    else
      call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 1 ")
    endif
            
    if(fill_order(2)=="x".or.fill_order(2)=="X") then
      l2 = real3(one,zero,zero)
    elseif(fill_order(2)=="y".or.fill_order(2)=="Y") then
      l2 = real3(zero,one,zero)
    elseif(fill_order(2)=="z".or.fill_order(2)=="Z") then
      l2 = real3(zero,zero,one)
    else
      call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 2 ")
    endif            
            
    if(fill_order(3)=="x".or.fill_order(3)=="X") then
      l3 = real3(one,zero,zero)
    elseif(fill_order(3)=="y".or.fill_order(3)=="Y") then
      l3 = real3(zero,one,zero) 
    elseif(fill_order(3)=="z".or.fill_order(3)=="Z") then
      l3 = real3(zero,zero,one)
    else
      call LPTLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 3 ")
    endif
  end subroutine Fill_Vectors

  !**********************************************************************
  ! copy i2 to i1
  !**********************************************************************
  subroutine GL_copy(this,i1,i2)
    implicit none
    class(VarList)::this
    integer,intent(in)::i1,i2

    if(i1==i2) return
    GPrtcl_id(i1)       = GPrtcl_id(i2)
    GPrtcl_pType(i1)    = GPrtcl_pType(i2)
    GPrtcl_usrMark(i1)  = GPrtcl_usrMark(i2)
    GPrtcl_PosR(i1)     = GPrtcl_PosR(i2)
    GPrtcl_linVel(:,i1) = GPrtcl_linVel(:,i2)
    GPrtcl_linAcc(:,i1) = GPrtcl_linAcc(:,i2)
    GPrtcl_FpForce(i1)  = GPrtcl_FpForce(i2)
    
  end subroutine GL_copy

  !**********************************************************************
  ! Reallocate particle varaibles 
  !**********************************************************************
  subroutine GL_ReallocatePrtclVar(this,np_new)
    implicit none
    class(VarList)::this
    integer,intent(in)::np_new

    ! locals
    integer:: sizep,sizen
    integer,dimension(:),allocatable:: IntVec
    type(real3),dimension(:),allocatable::Real3Vec
    type(real4),dimension(:),allocatable::Real4Vec
    type(real3),dimension(:,:),allocatable::Real3Arr

    sizep= this%mlocal
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen, np_new+1)
    sizen= min(sizen,LPT_Opt%numPrtcl)
    this%mlocal=sizen

    ! ======= integer vector part =======
    call move_alloc(GPrtcl_id, IntVec)
    allocate(GPrtcl_id(sizen))
    GPrtcl_id(1:sizep)=IntVec
    GPrtcl_id(sizep+1:sizen)=0

    call move_alloc(GPrtcl_pType, IntVec)
    allocate(GPrtcl_pType(sizen))
    GPrtcl_pType(1:sizep)=IntVec
    GPrtcl_pType(sizep+1:sizen)=1

    call move_alloc(GPrtcl_usrMark, IntVec)
    allocate(GPrtcl_usrMark(sizen))
    GPrtcl_usrMark(1:sizep)=IntVec
    GPrtcl_usrMark(sizep+1:sizen)=1
    deallocate(IntVec)

    ! ======= real3 vector part =======
    deallocate(GPrtcl_PosOld)
    allocate(GPrtcl_PosOld(sizen))
    deallocate(GPrtcl_VFluid)
    allocate(GPrtcl_VFluid(sizen))
    call move_alloc(GPrtcl_FpForce,Real3Vec)
    allocate(GPrtcl_FpForce(sizen))
    GPrtcl_FpForce(1:sizep)=Real3Vec
    GPrtcl_FpForce(sizep+1:sizen)=zero_r3
    
    ! ======= real3 matrix part =======
    call move_alloc(GPrtcl_linVel,Real3Arr)
    allocate(GPrtcl_linVel(this%tsize,sizen))
    GPrtcl_linVel(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linVel(1:this%tsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_linAcc,Real3Arr)
    allocate(GPrtcl_linAcc(this%tsize,sizen))
    GPrtcl_linAcc(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linAcc(1:this%tsize,sizep+1:sizen)=zero_r3
    deallocate(Real3Arr) 

    ! ======= real4 vercor part =======
    call move_alloc(GPrtcl_PosR,Real4Vec)
    allocate(GPrtcl_PosR(sizen))
    GPrtcl_PosR(1:sizep)=Real4Vec
    GPrtcl_PosR(sizep+1:sizen)=zero_r4
    deallocate(Real4Vec)

    call LPTLogInfo%CheckForError(ErrT_Pass," ReallocatePrtclVar"," Need to reallocate particle variables")
    call LPTLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    call LPTLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    call LPTLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine GL_ReallocatePrtclVar

end module LPT_Variables
