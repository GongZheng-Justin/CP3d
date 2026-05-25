#include "definitions_inc.f90"
module lp_Property
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d,only: nrank
#ifdef CFDSecondOrder
  use f2_Parameters,only: xnu,FluidDensity
#else
  use f4_Parameters,only: xnu,FluidDensity
#endif
  use lp_Parameters
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

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! initializing the size distribution with property 
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitPrtclProperty(this, chFile); implicit none
    class(PhysicalProperty)::this
    character(len=*),intent(in)::chFile
        
    ! locals
    real(RK):: sum_divided,rtemp 
    real(RK),dimension(:),allocatable:: Bin_Divided, Density, Diameter
    namelist/ParticlePhysicalProperty/Bin_Divided, Density, Diameter
    integer:: i,j,iTV(8),nPType,iUnit,ierr,sum_prtcl,bin_pnum,prdiff,bin_id
        
    nPType  = LPT_opt%numPrtcl_Type
    allocate( Bin_Divided(nPType))
    allocate( Density(nPType))
    allocate( Diameter(nPType))
    allocate( this%nPrtcl_in_Bin(nPType))
    allocate( this%Prtcl_PureProp(nPType))
          
    open(newunit=iUnit, file=chFile, status='old', form='formatted', IOSTAT=ierr)
    if(ierr/=0.and.nrank==0) call LPTLogInfo%CheckForError(ErrT_Abort,"InitPrtclProperty","Cannot open file: "//strip(chFile))
    read(iUnit, nml=ParticlePhysicalProperty)
    if(nrank==0) write(LPTLogInfo%iUnit, nml=ParticlePhysicalProperty)
    close(iUnit, IOSTAT=ierr)
        
    ! calculate this%nPrtcl_in_Bin
    sum_divided=0.0_RK
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
    if(prdiff>0) then
      call date_and_time(values=iTV); !iTV=0
      call random_seed(size= i)
      call random_seed(put = iTV(7)*iTV(8)+[(j,j=1,i)])
      do i=1, prdiff 
        call random_number(rtemp)
        bin_id = int(rtemp*nPType) + 1
        this%nPrtcl_in_Bin(bin_id) = this%nPrtcl_in_Bin(bin_id)  + 1
      enddo
     endif
        
     ! calculate particle properties
     do i = 1, nPType
       this%Prtcl_PureProp(i)%Density= Density(i)
       this%Prtcl_PureProp(i)%Radius = 0.5_RK*Diameter(i)
       this%Prtcl_PureProp(i)%Volume = 1.333333333333333_RK*Pi*(this%Prtcl_PureProp(i)%Radius)**3
       this%Prtcl_PureProp(i)%Mass   = Density(i)*this%Prtcl_PureProp(i)%Volume
       this%Prtcl_PureProp(i)%RelaxionTime= Density(i)*Diameter(i)*Diameter(i)/(18.00_RK*xnu*FluidDensity)
     enddo       
  end subroutine InitPrtclProperty

end module lp_Property
