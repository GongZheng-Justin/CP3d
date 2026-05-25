#include "definitions_inc.f90"
module ap_Property
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d,only: nrank
#ifdef CFDSecondOrder
  use f2_Parameters,only: xnu,FluidDensity
#else
  use f4_Parameters,only: xnu,FluidDensity
#endif
  use ap_Parameters
  implicit none
  private     
    
  type PureProperty
    real(RK):: DiffuseR         ! Dr: Orientation diffusivity
    real(RK):: OneDtwoB         ! 1/2B
    real(RK):: SwimVelocityMag  ! Vs: Swim velocity
    real(RK):: Diameter
  end type PureProperty
    
  type PhysicalProperty
    integer,allocatable,dimension(:) :: nPrtcl_in_Bin  
    type(pureProperty),allocatable,dimension(:) :: Prtcl_PureProp     
  contains
    procedure:: InitPrtclProperty
  end type PhysicalProperty
  type(PhysicalProperty),public::ATPProperty

contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! initializing the size distribution with property 
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitPrtclProperty(this, chFile); implicit none
    class(PhysicalProperty)::this
    character(len=*),intent(in)::chFile
        
    ! locals
    real(RK)::sum_divided,rtemp 
    real(RK),dimension(:),allocatable::Bin_Divided,DiffuseR,OneDtwoB,SwimVelocityMag,Diameter
    namelist/SwimOptions/Bin_Divided,DiffuseR,OneDtwoB,SwimVelocityMag,Diameter
    integer:: i,j,iTV(8),nPType,iUnit,ierr,sum_prtcl,bin_pnum,prdiff,bin_id
        
    nPType  = ATP_opt%numPrtcl_Type
    allocate( Bin_Divided(nPType))
    allocate( DiffuseR(nPType))
    allocate( OneDtwoB(nPType))
    allocate( SwimVelocityMag(nPType))
    allocate( Diameter(nPType))
    allocate( this%nPrtcl_in_Bin(nPType))
    allocate( this%Prtcl_PureProp(nPType))
          
    open(newunit=iUnit, file=chFile, status='old', form='formatted', IOSTAT=ierr)
    if(ierr/=0.and.nrank==0) call ATPLogInfo%CheckForError(ErrT_Abort,"InitPrtclProperty","Cannot open file: "//strip(chFile))
    read(iUnit, nml=SwimOptions)
    if(nrank==0)write(ATPLogInfo%iUnit, nml=SwimOptions)
    close(iUnit, IOSTAT=ierr)
        
    ! calculate this%nPrtcl_in_Bin
    sum_divided=0.0_RK
    do i=1,nPType
      sum_divided = sum_divided + Bin_Divided(i)
    enddo
    sum_prtcl = 0
    do i=1,nPType
      bin_pnum = int(ATP_opt%numPrtcl*Bin_Divided(i)/sum_divided)
      this%nPrtcl_in_Bin(i)= bin_pnum
      sum_prtcl = sum_prtcl + bin_pnum
    enddo
    prdiff = ATP_opt%numPrtcl - sum_prtcl
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
       this%Prtcl_PureProp(i)%DiffuseR= DiffuseR(i)
       this%Prtcl_PureProp(i)%OneDtwoB= OneDtwoB(i)
       this%Prtcl_PureProp(i)%SwimVelocityMag= SwimVelocityMag(i)
       this%Prtcl_PureProp(i)%Diameter= Diameter(i)
     enddo       
  end subroutine InitPrtclProperty

end module ap_Property
