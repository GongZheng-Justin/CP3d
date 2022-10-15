!********************************************************************!
!*    file name  : Prtcl_Variables.f90                              *!
!*    module name: Prtcl_Variables                                  *!  
!*                                                                  *!
!*    purpose:                                                      *! 
!*      1) All datas required to represent spherical particles      *!
!*      2) Initialize all the particle variables                    *!
!*                                                                  *!
!*  Author: Zheng Gong           Date: 23:Feb:2020                  *!
!*                                                                  *!
!********************************************************************!

module Prtcl_Variables
  use MPI
  use Prtcl_TypeDef
  use Prtcl_LogInfo
  use Prtcl_Property
  use Prtcl_decomp_2d
  use Prtcl_Parameters
#if defined(CFDDEM) || defined(CFDACM)
  use m_Decomp2d,only: nrank
#endif
  implicit none
  private

  integer,dimension(:),allocatable,public:: GPrtcl_id
  integer,dimension(:),allocatable,public:: GPrtcl_pType
  integer,dimension(:),allocatable,public:: GPrtcl_usrMark
  type(real4),dimension(:),allocatable,public::   GPrtcl_PosR
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linVel
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_linAcc
  type(real3),dimension(:),allocatable,public::   GPrtcl_theta
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_rotVel
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_rotAcc
  type(real3),dimension(:),allocatable,public::   GPrtcl_cntctForce
  type(real3),dimension(:),allocatable,public::   GPrtcl_torque

  integer,dimension(:),allocatable,public::       GPFix_id
  integer,dimension(:),allocatable,public::       GPFix_pType
  type(real4),dimension(:),allocatable,public::   GPFix_PosR
#ifdef CFDDEM
  type(real3),dimension(:,:),allocatable,public:: GPFix_VFluid 
  type(real3),dimension(:),  allocatable,public:: GPrtcl_FpForce,GPrtcl_FpForce_old,GPrtcl_linVelOld
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_Vfluid
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_BassetData
  type(real3),dimension(:,:),allocatable,public:: GPFix_BassetData
#endif
#ifdef CFDACM
  type(real3),dimension(:),allocatable,public::   GPrtcl_FpForce
  type(real3),dimension(:),allocatable,public::   GPrtcl_FpTorque
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_FluidIntegrate
  type(real3),dimension(:,:),allocatable,public:: GPrtcl_FluidIntOld
  type(real4),dimension(:),allocatable,public::   GPrtcl_PosOld
  character(len=1),dimension(:),allocatable,public:: GPrtcl_HighSt
#endif

  type VarList
    integer:: nlocal     ! number of particles in local processor
    integer:: mlocal     ! the possible maxium # of particles in local processor
    integer:: mlocalFix  ! number of fixed particles in local processor

    integer:: nGhost_CS
    integer:: mGhost_CS
    integer:: nGhostFix_CS

    integer:: tsize  ! size for translational veloctity and acceleration
    integer:: rsize  ! size for rotational veloctity and acceleration
  contains
    procedure:: AllocateAllVar     => GL_AllocateAllVar
#if !defined(CFDDEM) && !defined(CFDACM)
    procedure:: MakingAllPrtcl     => GL_MakingAllPrtcl
#endif
    procedure:: ReallocatePrtclVar => GL_ReallocatePrtclVar
    procedure:: copy               => GL_copy
  end type VarList
  type(VarList),public::GPrtcl_list 

  ! Ghost particle variables
  integer,dimension(:),allocatable,public:: GhostP_id
  integer,dimension(:),allocatable,public:: GhostP_pType
  type(real4),dimension(:),allocatable,public:: GhostP_PosR
  type(real3),dimension(:),allocatable,public:: GhostP_linVel
  type(real3),dimension(:),allocatable,public:: GhostP_rotVel
    
contains

  !**********************************************************************
  ! GL_AllocateAllVar
  !**********************************************************************
  subroutine GL_AllocateAllVar(this)
    implicit none
    class(VarList)::this

    ! locals
    type(real3)::SimLen
    real(RK)::xst,xed,yst,yed,zst,zed,vol_tot,vol_local
    integer:: mlocal,numPrtcl,mlocalFix,numPrtclFix,ierrTmp,ierror=0

    numPrtcl    = DEM_opt%numPrtcl
    numPrtclFix = DEM_opt%numPrtclFix

    ! step0: determine initial misze
    xst=DEM_decomp%xSt; xed=DEM_decomp%xEd
    yst=DEM_decomp%ySt; yed=DEM_decomp%yEd
    zst=DEM_decomp%zSt; zed=DEM_decomp%zEd

    ! step1: allocating memory for particles
    if(DEM_opt%PI_Method==PIM_FE) then
      this%tsize=1
    elseif (DEM_opt%PI_Method==PIM_AB2) then
      this%tsize=2
    elseif(DEM_opt%PI_Method==PIM_AB3) then
      this%tsize=3
    endif
    if(DEM_opt%PRI_Method==PIM_FE) then
      this%rsize=1
    elseif (DEM_opt%PRI_Method==PIM_AB2) then
      this%rsize=2
    elseif(DEM_opt%PRI_Method==PIM_AB3) then
      this%rsize=3
    endif 

    SimLen = DEM_Opt%SimDomain_max - DEM_Opt%SimDomain_min
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
    allocate(GPrtcl_theta(mlocal),             Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_rotVel(this%rsize,mlocal), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_rotAcc(this%rsize,mlocal), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_cntctForce(mlocal),        Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_torque(mlocal),            Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation failed 1")

    mlocalFix = int(vol_local/vol_tot*real(numPrtclFix,kind=RK))
    mlocalFix = int(1.5_RK*mlocalFix)
    mlocalFix = max(mlocalFix, 10)
    mlocalFix = min(mlocalFix, numPrtclFix)
    allocate(GPFix_id(mlocalFix),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPFix_pType(mlocalFix),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPFix_PosR(mlocalFix), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation failed 2")
#ifdef CFDDEM
    allocate(GPrtcl_FpForce(mlocal),    Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_FpForce_old(mlocal),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_linVelOld(mlocal),  Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_Vfluid(2,mlocal),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation Fpforce failed")
    GPrtcl_FpForce     = zero_r3
    GPrtcl_FpForce_old = zero_r3
    GPrtcl_linVelOld   = zero_r3
    GPrtcl_Vfluid      = zero_r3
    if(Is_clc_Basset) then
      allocate(GPrtcl_BassetData(GPrtcl_BassetSeq%nDataLen, mlocal), Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"GL_AllocateAllVar: ","Allocation GPrtcl_BassetData failed")
      GPrtcl_BassetData= zero_r3
    endif
#endif
#ifdef CFDACM
    allocate(GPrtcl_FpForce(mlocal),         Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_FpTorque(mlocal),        Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_FluidIntegrate(2,mlocal),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_FluidIntOld(2,mlocal),   Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    allocate(GPrtcl_HighSt(mlocal),          Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(IBM_Scheme==2) allocate(GPrtcl_PosOld(mlocal),Stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"AllocateAllVar","Allocation Fpforce failed")
    GPrtcl_FpForce  = zero_r3
    GPrtcl_FpTorque = zero_r3
    GPrtcl_FluidIntegrate= zero_r3
    GPrtcl_FluidIntOld   = zero_r3
    GPrtcl_HighSt = "N"
    if(IBM_Scheme==2) GPrtcl_PosOld=zero_r4
#endif

    GPrtcl_id = 0
    GPrtcl_pType = 1
    GPrtcl_usrMark = 1
    GPrtcl_PosR = zero_r4
    GPrtcl_linVel  = zero_r3
    GPrtcl_linAcc  = zero_r3
    GPrtcl_theta   = zero_r3
    GPrtcl_rotVel  = zero_r3
    GPrtcl_rotAcc  = zero_r3
    GPrtcl_cntctForce = zero_r3
    GPrtcl_torque     = zero_r3
    GPFix_id = 0
    GPFix_pType = 1
    GPFix_PosR = zero_r4

    ! step2: Initialize this%mlocal,this%nlocal
    this%nlocal    = 0
    this%mlocal    = mlocal
    this%mlocalFix = mlocalFix
  end subroutine GL_AllocateAllVar

#if !defined(CFDDEM) && !defined(CFDACM)
  !**********************************************************************
  ! GL_MakingAllPrtcl
  !**********************************************************************
  subroutine GL_MakingAllPrtcl(this, chFile)
    implicit none
    class(VarList)::this
    character(*),intent(in)::chFile
        
    ! locals
    real(RK)::Distance_Ratio,VelMag
    character,dimension(3)::Fill_order
    real(RK),dimension(3)::MkPrtclMinpoint,MkPrtclMaxpoint
    type(real3):: MkPrtclMinpoint_real3,MkPrtclMaxpoint_real3
    NAMELIST /ParticleMakingOption/MkPrtclMinpoint,MkPrtclMaxpoint, Fill_order, Distance_Ratio, VelMag
    
    real(RK),dimension(3):: vel_dir
    integer,dimension(:),allocatable:: sum_bin
    real(RK)::xst,xed,yst,yed,zst,zed,maxRad,rvelt,randt
    type(real3)::cntr,l1_vec,l2_vec,l3_vec,lmin_p,lmax_p,dx,pos1,pos2
    integer::i,j,k,code,bin_id,numPrtcl,nUnitFile,ierror,nlocal,nlocal_sum
   
    open(newunit=nUnitFile, file=chFile,status='old',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl", "Cannot open file: "//trim(chFile))
    read(nUnitFile, nml=ParticleMakingOption)
    if(nrank==0)write(DEMLogInfo%nUnit, nml=ParticleMakingOption)
    close(nUnitFile,IOSTAT=ierror)
    MkPrtclMinpoint_real3%x=max(MkPrtclMinpoint(1),DEM_Opt%SimDomain_min%x)
    MkPrtclMinpoint_real3%y=max(MkPrtclMinpoint(2),DEM_Opt%SimDomain_min%y)
    MkPrtclMinpoint_real3%z=max(MkPrtclMinpoint(3),DEM_Opt%SimDomain_min%z)
    MkPrtclMaxpoint_real3%x=min(MkPrtclMaxpoint(1),DEM_Opt%SimDomain_max%x)
    MkPrtclMaxpoint_real3%y=min(MkPrtclMaxpoint(2),DEM_Opt%SimDomain_max%y)
    MkPrtclMaxpoint_real3%z=min(MkPrtclMaxpoint(3),DEM_Opt%SimDomain_max%z)
    numPrtcl = DEM_opt%numPrtcl
        
    xst=DEM_decomp%xSt; xed=DEM_decomp%xEd
    yst=DEM_decomp%ySt; yed=DEM_decomp%yEd
    zst=DEM_decomp%zSt; zed=DEM_decomp%zEd

    ! step3: Assign Prtcl_Type, Prtcl_id and Prtcl_PosR
    call Fill_Vectors(Fill_order, l1_vec, l2_vec, l3_vec)
    call system_clock(count=code); !code=0
    call random_seed(size = j)
    call random_seed(put = code+63946*(/ (i - 1, i = 1, j) /))
    maxRad = maxval( DEMProperty%Prtcl_PureProp%Radius )
    lmin_p = MkPrtclMinpoint_real3 + Distance_Ratio* maxRad*one_r3
    lmax_p = MkPrtclMaxpoint_real3 - Distance_Ratio* maxRad*one_r3
    dx = Distance_Ratio *two*maxRad* one_r3
        
    allocate(sum_bin(DEM_opt%numPrtcl_Type))
    sum_bin(1)=DEMProperty%nPrtcl_in_Bin(1)
    do j=2, DEM_opt%numPrtcl_Type
      sum_bin(j)=sum_bin(j-1)+DEMProperty%nPrtcl_in_Bin(j)
    enddo

    cntr = lmin_p       ! start point
    nlocal=0
    do i=1,numPrtcl
      do j=1,DEM_opt%numPrtcl_Type
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
        GPrtcl_PosR(nlocal)%w = DEMProperty%Prtcl_PureProp(bin_id)%Radius
      endif

      cntr = cntr + l1_vec * dx
      if((l1_vec.dot.cntr)>=(l1_vec.dot.lmax_p))then
        cntr = (lmin_p*l1_vec) + ((cntr+dx)*l2_vec)+(cntr*l3_vec)
        if((l2_vec.dot.cntr)>=(l2_vec.dot.lmax_p))then
          cntr = (cntr*l1_vec)+(lmin_p*l2_vec)+((cntr+dx)*l3_vec)
          if((l3_vec.dot.cntr) >= (l3_vec.dot.lmax_p) .and. nrank==0) then
            call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Not enough space for positioning" ) 
          endif
        endif
      endif
    enddo

    call MPI_REDUCE(nlocal,nlocal_sum,1,int_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    if(nlocal_sum/= numPrtcl .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: "," nlocal_sum/= numPrtcl " )
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
        
    ! step3: assign the remaining variables
    do i=1,this%nlocal
      call random_number(vel_dir)
      vel_dir=two*vel_dir-one
      rvelt=sqrt(vel_dir(1)**2+vel_dir(2)**2+vel_dir(3)**2)
      if(rvelt>1.0E-10_RK) then
        vel_dir(1)=vel_dir(1)/rvelt
        vel_dir(2)=vel_dir(2)/rvelt
        vel_dir(3)=vel_dir(3)/rvelt
      else
        vel_dir(1)= -one*l3_vec%x
        vel_dir(2)= -one*l3_vec%y
        vel_dir(3)= -one*l3_vec%z
      endif
      GPrtcl_linVel(1,i)%x=VelMag*vel_dir(1)
      GPrtcl_linVel(1,i)%y=VelMag*vel_dir(2)
      GPrtcl_linVel(1,i)%z=VelMag*vel_dir(3)
    enddo
#ifdef DEMObliqueCollideDry
    do i=1,this%nlocal
      GPrtcl_PosR(i)%x=(xSt+xEd)*half
      GPrtcl_PosR(i)%y=ySt+3.0_RK*DEMProperty%Prtcl_PureProp(GPrtcl_pType(i))%Radius
      GPrtcl_PosR(i)%z=(zSt+zEd)*half
      GPrtcl_linVel(1,i)= DEM_opt%Gravity
    enddo
    DEM_opt%Gravity=zero_r3
#endif
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
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 1 ")
    endif
            
    if(fill_order(2)=="x".or.fill_order(2)=="X") then
      l2 = real3(one,zero,zero)
    elseif(fill_order(2)=="y".or.fill_order(2)=="Y") then
      l2 = real3(zero,one,zero)
    elseif(fill_order(2)=="z".or.fill_order(2)=="Z") then
      l2 = real3(zero,zero,one)
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 2 ")
    endif            
            
    if(fill_order(3)=="x".or.fill_order(3)=="X") then
      l3 = real3(one,zero,zero)
    elseif(fill_order(3)=="y".or.fill_order(3)=="Y") then
      l3 = real3(zero,one,zero) 
    elseif(fill_order(3)=="z".or.fill_order(3)=="Z") then
      l3 = real3(zero,zero,one)
    else
      call DEMLogInfo%CheckForError(ErrT_Abort,"MakingAllPrtcl: ","Fill_order wrong: 3 ")
    endif
  end subroutine Fill_Vectors
#endif

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
    GPrtcl_theta(i1)    = GPrtcl_theta(i2)
    GPrtcl_rotVel(:,i1) = GPrtcl_rotVel(:,i2)
    GPrtcl_rotAcc(:,i1) = GPrtcl_rotAcc(:,i2)
#ifdef CFDDEM
    GPrtcl_FpForce(i1)     = GPrtcl_FpForce(i2)
    GPrtcl_FpForce_old(i1) = GPrtcl_FpForce_old(i2)
    GPrtcl_linVelOld(i1)   = GPrtcl_linVelOld(i2)
    GPrtcl_Vfluid(:,i1)    = GPrtcl_Vfluid(:,i2)
    if(is_clc_Basset)  GPrtcl_BassetData(:,i1)= GPrtcl_BassetData(:,i2)
#endif
#ifdef CFDACM
    GPrtcl_FpForce(i1) = GPrtcl_FpForce(i2)
    GPrtcl_FpTorque(i1)= GPrtcl_FpTorque(i2)
    GPrtcl_FluidIntOld(:,i1)= GPrtcl_FluidIntOld(:,i2)
    if(IBM_Scheme==2) GPrtcl_PosOld(i1) = GPrtcl_PosOld(i2)
#endif
  end subroutine GL_copy

  !**********************************************************************
  ! Reallocate particle varaibles 
  !**********************************************************************
  subroutine GL_ReallocatePrtclVar(this,np_new)
    implicit none
    class(VarList)::this
    integer,intent(in)::np_new

    ! locals
    integer:: sizep,sizen,ierrTmp,ierror=0
    integer,dimension(:),allocatable:: IntVec
    type(real3),dimension(:),allocatable::Real3Vec
    type(real4),dimension(:),allocatable::Real4Vec
    type(real3),dimension(:,:),allocatable::Real3Arr
#ifdef CFDACM
    character(len=1),dimension(:),allocatable::ChaVec
#endif

    sizep= this%mlocal
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= max(sizen, np_new+1)
    sizen= min(sizen,DEM_Opt%numPrtcl)
    this%mlocal=sizen

    ! ======= integer vector part =======
    call move_alloc(GPrtcl_id, IntVec)
    allocate(GPrtcl_id(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_id(1:sizep)=IntVec
    GPrtcl_id(sizep+1:sizen)=0

    call move_alloc(GPrtcl_pType, IntVec)
    allocate(GPrtcl_pType(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_pType(1:sizep)=IntVec
    GPrtcl_pType(sizep+1:sizen)=1

    call move_alloc(GPrtcl_usrMark, IntVec)
    allocate(GPrtcl_usrMark(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_usrMark(1:sizep)=IntVec
    GPrtcl_usrMark(sizep+1:sizen)=1
    deallocate(IntVec)

    ! ======= real3 vercor part =======
    call move_alloc(GPrtcl_theta,Real3Vec)
    allocate(GPrtcl_theta(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_theta(1:sizep)=Real3Vec
    GPrtcl_theta(sizep+1:sizen)=zero_r3
#ifdef CFDDEM
    call move_alloc(GPrtcl_FpForce,Real3Vec)
    allocate(GPrtcl_FpForce(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_FpForce(1:sizep)=Real3Vec
    GPrtcl_FpForce(sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_FpForce_old,Real3Vec)
    allocate(GPrtcl_FpForce_old(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_FpForce_old(1:sizep)=Real3Vec
    GPrtcl_FpForce_old(sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_linVelOld,Real3Vec)
    allocate(GPrtcl_linVelOld(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_linVelOld(1:sizep)=Real3Vec
    GPrtcl_linVelOld(sizep+1:sizen)=zero_r3
#endif
#ifdef CFDACM
    call move_alloc(GPrtcl_FpForce,Real3Vec)
    allocate(GPrtcl_FpForce(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_FpForce(1:sizep)=Real3Vec
    GPrtcl_FpForce(sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_FpTorque,Real3Vec)
    allocate(GPrtcl_FpTorque(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_FpTorque(1:sizep)=Real3Vec
    GPrtcl_FpTorque(sizep+1:sizen)=zero_r3
#endif
    deallocate(Real3Vec)
 
    deallocate(GPrtcl_cntctForce)
    allocate(GPrtcl_cntctForce(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)

    deallocate(GPrtcl_torque)
    allocate(GPrtcl_torque(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    
    ! ======= real3 matrix part =======
    call move_alloc(GPrtcl_linVel,Real3Arr)
    allocate(GPrtcl_linVel(this%tsize,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_linVel(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linVel(1:this%tsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_linAcc,Real3Arr)
    allocate(GPrtcl_linAcc(this%tsize,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_linAcc(1:this%tsize,1:sizep)=Real3Arr
    GPrtcl_linAcc(1:this%tsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_rotVel,Real3Arr)
    allocate(GPrtcl_rotVel(this%rsize,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_rotVel(1:this%rsize,1:sizep)=Real3Arr
    GPrtcl_rotVel(1:this%rsize,sizep+1:sizen)=zero_r3

    call move_alloc(GPrtcl_rotAcc,Real3Arr)
    allocate(GPrtcl_rotAcc(this%rsize,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_rotAcc(1:this%rsize,1:sizep)=Real3Arr
    GPrtcl_rotAcc(1:this%rsize,sizep+1:sizen)=zero_r3
#ifdef CFDDEM
    call move_alloc(GPrtcl_Vfluid,Real3Arr)
    allocate(GPrtcl_Vfluid(2,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_Vfluid(1:2,1:sizep)=Real3Arr
    GPrtcl_Vfluid(1:2,sizep+1:sizen)=zero_r3

    if(is_clc_Basset) then
      call move_alloc(GPrtcl_BassetData,Real3Arr)
      allocate(GPrtcl_BassetData(GPrtcl_BassetSeq%nDataLen ,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      GPrtcl_BassetData(1:GPrtcl_BassetSeq%nDataLen, 1:sizep)=Real3Arr
      GPrtcl_BassetData(1:GPrtcl_BassetSeq%nDataLen, sizep+1:sizen)=zero_r3
    endif
#endif
#ifdef CFDACM
    call move_alloc(GPrtcl_FluidIntOld,Real3Arr)
    allocate(GPrtcl_FluidIntOld(2,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_FluidIntOld(1:2,1:sizep)=Real3Arr    
    GPrtcl_FluidIntOld(1:2,sizep+1:sizen)=zero_r3    

    deallocate(GPrtcl_FluidIntegrate)
    allocate(GPrtcl_FluidIntegrate(2,sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
#endif   
    deallocate(Real3Arr) 

    ! ======= real4 vercor part =======
    call move_alloc(GPrtcl_PosR,Real4Vec)
    allocate(GPrtcl_PosR(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_PosR(1:sizep)=Real4Vec
    GPrtcl_PosR(sizep+1:sizen)=zero_r4
#ifdef CFDACM
    if(IBM_Scheme==2) then
      call move_alloc(GPrtcl_PosOld,Real4Vec)
      allocate(GPrtcl_PosOld(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
      GPrtcl_PosOld(1:sizep)=Real4Vec
      GPrtcl_PosOld(sizep+1:sizen)=zero_r4
    endif
#endif
    deallocate(Real4Vec)

#ifdef CFDACM
    ! ======= character vercor part =======
    call move_alloc(GPrtcl_HighSt,ChaVec)
    allocate(GPrtcl_HighSt(sizen),stat=ierrTmp);ierror=ierror+abs(ierrTmp)
    GPrtcl_HighSt(1:sizep)=ChaVec
    GPrtcl_HighSt(sizep+1:sizen)="N"
    deallocate(ChaVec)
#endif

    if(ierror/=0) then
      call DEMLogInfo%CheckForError(ErrT_Abort," GL_ReallocatePrtclVar"," Reallocate wrong!")
      call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    endif
    !call DEMLogInfo%CheckForError(ErrT_Pass," ReallocatePrtclVar"," Need to reallocate particle variables")
    !call DEMLogInfo%OutInfo("The present processor  is :"//trim(num2str(nrank)),3)
    !call DEMLogInfo%OutInfo("Previous matirx length is :"//trim(num2str(sizep)),3)
    !call DEMLogInfo%OutInfo("Updated  matirx length is :"//trim(num2str(sizen)),3)
  end subroutine GL_ReallocatePrtclVar

end module Prtcl_Variables
