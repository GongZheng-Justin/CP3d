module LPT_Property
  use LPT_TypeDef
  use LPT_LogInfo
  use LPT_Parameters
  use m_Decomp2d,only: nrank
  use m_Parameters,only: xnu,FluidDensity
  implicit none
  private     
    
  type PureProperty
    real(RK):: Radius
    real(RK):: Density
    real(RK):: Mass
    real(RK):: Volume
    real(RK):: RelaxionTime
  end type PureProperty
    
  type PhysicalProperty
    integer,allocatable,dimension(:) :: nPrtcl_in_Bin  
    type(pureProperty),allocatable,dimension(:) :: Prtcl_PureProp     
  contains
    procedure:: InitPrtclProperty
  end type PhysicalProperty
  type(PhysicalProperty),public::LPTProperty

contains

  !*******************************************************
  ! initializing the size distribution with property 
  !*******************************************************
  subroutine InitPrtclProperty( this,  chFile )
    implicit none
    class(PhysicalProperty)::this
    character(*),intent(in)::chFile
        
    ! locals
    real(RK):: sum_divided,rtemp 
    real(RK),dimension(:),allocatable:: Bin_Divided, Density, Diameter
    namelist/ParticlePhysicalProperty/Bin_Divided, Density, Diameter
    integer:: i,j,code,nPType,nUnitFile,ierror,sum_prtcl,bin_pnum,prdiff,bin_id
        
    nPType  = LPT_opt%numPrtcl_Type
    allocate( Bin_Divided(nPType))
    allocate( Density(nPType))
    allocate( Diameter(nPType))
    allocate( this%nPrtcl_in_Bin(nPType))
    allocate( this%Prtcl_PureProp(nPType))
          
    open(newunit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=ierror)
    if(ierror/=0.and.nrank==0) call LPTLogInfo%CheckForError(ErrT_Abort,"InitPrtclProperty","Cannot open file: "//trim(chFile))
    read(nUnitFile, nml=ParticlePhysicalProperty)
    if(nrank==0)write(LPTLogInfo%nUnit, nml=ParticlePhysicalProperty)
    close(nUnitFile, IOSTAT=ierror)
        
    ! calculate this%nPrtcl_in_Bin
    sum_divided=zero
    do i=1,nPType
      sum_divided = sum_divided + Bin_Divided(i)
    enddo
    sum_prtcl = 0
    do i=1,nPType
      bin_pnum = int(LPT_opt%numPrtcl*Bin_Divided(i)/sum_divided)
      this%nPrtcl_in_Bin(i)= bin_pnum
      sum_prtcl = sum_prtcl + bin_pnum
    enddo
    prdiff = LPT_opt%numPrtcl - sum_prtcl
    if( prdiff > 0 ) then
      call system_clock(count=code) !code=0
      call random_seed(size = j)
      call random_seed(put = code+63946*(/(i-1,i=1,j)/))
      do i=1, prdiff 
        call random_number(rtemp)
        bin_id = int(rtemp*nPType) + 1
        this%nPrtcl_in_Bin(bin_id) = this%nPrtcl_in_Bin(bin_id)  + 1
        enddo
     endif
        
     ! calculate particle properties
     do i = 1, nPType
       this%Prtcl_PureProp(i)%Density= Density(i)
       this%Prtcl_PureProp(i)%Radius = half * Diameter(i)
       this%Prtcl_PureProp(i)%Volume = four/three*Pi*(this%Prtcl_PureProp(i)%Radius)**3
       this%Prtcl_PureProp(i)%Mass   = Density(i)*this%Prtcl_PureProp(i)%Volume
       this%Prtcl_PureProp(i)%RelaxionTime= Density(i)*Diameter(i)*Diameter(i)/(18.00_RK*xnu*FluidDensity)
     end do       
  end subroutine InitPrtclProperty

end module LPT_Property
