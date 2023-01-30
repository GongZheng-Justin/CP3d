module ATP_Variables
  use MPI
  use m_TypeDef
  use m_LogInfo
  use ATP_Property
  use ATP_decomp_2d
  use ATP_Parameters
  use m_Decomp2d,only: nrank
  use m_Parameters,only:xlx,yly,zlz
  implicit none
  private
  integer,dimension(:),allocatable,public:: GPrtcl_id
  integer,dimension(:),allocatable,public:: GPrtcl_pType
  integer,dimension(:),allocatable,public:: GPrtcl_usrMark
  type(real3),dimension(:),allocatable,public::   GPrtcl_PosR
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linVel
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_SwimAcc
  type(real3),dimension(:),allocatable,public:: GPrtcl_PosOld
  type(real3),dimension(:),allocatable,public:: GPrtcl_SwimDir
  type(real3),dimension(:),allocatable,public:: GPrtcl_MoveDistance
  
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

    numPrtcl    = ATP_opt%numPrtcl

    ! step1: allocating memory for particles
    if (ATP_opt%PI_Method==PIM_AB2) then
      this%tsize=2
    elseif(ATP_opt%PI_Method==PIM_AB3) then
      this%tsize=3
    endif

    ! step0: determine initial msize
    xst=ATP_decomp%xSt; xed=ATP_decomp%xEd
    yst=ATP_decomp%ySt; yed=ATP_decomp%yEd
    zst=ATP_decomp%zSt; zed=ATP_decomp%zEd

    SimLen = ATP_Opt%SimDomain_max - ATP_Opt%SimDomain_min
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
    allocate(GPrtcl_SwimAcc(this%tsize,mlocal),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_PosOld(mlocal),            Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_SwimDir(mlocal),           Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_MoveDistance(mlocal),       Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call ATPLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar","Allocation failed 1")   
    
    GPrtcl_id = 0
    GPrtcl_pType = 1
    GPrtcl_usrMark = 1
    GPrtcl_PosR = zero_r3
    GPrtcl_linVel  = zero_r3
    GPrtcl_SwimAcc  = zero_r3
    GPrtcl_PosOld  = zero_r3
    GPrtcl_SwimDir  = zero_r3
    GPrtcl_MoveDistance = zero_r3
    
    ! step2: Initialize this%mlocal,this%nlocal
    this%nlocal  = 0
    this%mlocal  = mlocal
  end subroutine GL_AllocateAllVar

  !**********************************************************************
  ! GL_MakingAllPrtcl
  !**********************************************************************
  subroutine GL_MakingAllPrtcl(this)
    implicit none
    class(VarList)::this
        
    ! locals
    type(real3)::cntr
    integer,dimension(:),allocatable::sum_bin
    real(RK)::xst,xed,yst,yed,zst,zed,realt(3),MagTmp
    integer::i,j,k,iTV(8),bin_id,numPrtcl,nUnitFile,ierror,nlocal,nlocal_sum
    
    numPrtcl = ATP_opt%numPrtcl
    call date_and_time(values=iTV); !iTV=0
    call random_seed(size= i)
    call random_seed(put = iTV(7)*iTV(8)*[(j,j=1,i)])
    xst=ATP_decomp%xSt; xed=ATP_decomp%xEd
    yst=ATP_decomp%ySt; yed=ATP_decomp%yEd
    zst=ATP_decomp%zSt; zed=ATP_decomp%zEd
            
    allocate(sum_bin(ATP_opt%numPrtcl_Type))
    sum_bin(1)=ATPProperty%nPrtcl_in_Bin(1)
    do j=2, ATP_opt%numPrtcl_Type
      sum_bin(j)=sum_bin(j-1)+ATPProperty%nPrtcl_in_Bin(j)
    enddo

    nlocal=0; bin_id=-1
    do i=1,numPrtcl
      do j=1,ATP_opt%numPrtcl_Type
        if(i<=sum_bin(j)) then
          bin_id = j; exit
        endif
      enddo
      call random_number(realt)
      realt=2.0_RK*realt-1.0_RK
      MagTmp=sqrt(realt(1)*realt(1)+realt(2)*realt(2)+realt(3)*realt(3))
      realt=realt/MagTmp

      cntr%x=  0.5_RK*xlx
      cntr%y=  0.5_RK*yly
      cntr%z=  0.5_RK*zlz        
      if(cntr%x>=xst.and.cntr%x<xed.and.cntr%y>=yst.and.cntr%y<yed.and.cntr%z>=zst.and.cntr%z<zed) then
        if(nlocal>=this%mlocal) call this%ReallocatePrtclVar(nlocal)
        nlocal=nlocal+1
        GPrtcl_id(nlocal) = i
        GPrtcl_pType(nlocal)= bin_id
        GPrtcl_PosR(nlocal) = cntr
        GPrtcl_SwimDir(nlocal)=realt
      endif
    enddo
    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nlocal_sum/= numPrtcl .and. nrank==0) then
      call ATPLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: "," nlocal_sum/= numPrtcl " )
    endif
    this%nlocal = nlocal
    deallocate(sum_bin)
  end subroutine GL_MakingAllPrtcl

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
    GPrtcl_SwimDir(i1)  = GPrtcl_SwimDir(i2)
    GPrtcl_linVel(:,i1) = GPrtcl_linVel(:,i2)
    GPrtcl_SwimAcc(:,i1)= GPrtcl_SwimAcc(:,i2)
    GPrtcl_MoveDistance(i1)= GPrtcl_MoveDistance(i2)
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
    type(real3),dimension(:,:),allocatable::Real3Arr

    sizep= this%mlocal
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen, np_new+1)
    sizen= min(sizen,ATP_Opt%numPrtcl)
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
    call move_alloc(GPrtcl_SwimDir,Real3Vec)
    allocate(GPrtcl_SwimDir(sizen))
    GPrtcl_SwimDir(1:sizep)=Real3Vec
    GPrtcl_SwimDir(sizep+1:sizen)=zero_r3
    
    ! ======= real3 matrix part =======
    call move_alloc(GPrtcl_linVel,Real3Arr)
    allocate(GPrtcl_linVel(this%tsize,sizen))
    GPrtcl_linVel(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linVel(1:this%tsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_SwimAcc,Real3Arr)
    allocate(GPrtcl_SwimAcc(this%tsize,sizen))
    GPrtcl_SwimAcc(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_SwimAcc(1:this%tsize,sizep+1:sizen)=zero_r3
    deallocate(Real3Arr) 

    ! ======= real4 vercor part =======
    call move_alloc(GPrtcl_PosR,Real3Vec)
    allocate(GPrtcl_PosR(sizen))
    GPrtcl_PosR(1:sizep)=Real3Vec
    GPrtcl_PosR(sizep+1:sizen)=zero_r3
    call move_alloc(GPrtcl_MoveDistance,Real3Vec)
    allocate(GPrtcl_MoveDistance(sizen))
    GPrtcl_MoveDistance(1:sizep)=Real3Vec
    GPrtcl_MoveDistance(sizep+1:sizen)=zero_r3    
    deallocate(Real3Vec)

    call ATPLogInfo%CheckForError(ErrT_Pass," ReallocatePrtclVar"," Need to reallocate particle variables")
    call ATPLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    call ATPLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    call ATPLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine GL_ReallocatePrtclVar

end module ATP_Variables
