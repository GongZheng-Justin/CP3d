#include "definitions_inc.f90"
module f2_FlowCase
  use mc_TypeDef
  use mc_Decomp2d
  use mc_LogInfo
  use f2_Parameters
  use f2_FlowType_Channel
  use f2_FlowType_Duct
  use f2_FlowType_TGVortex
  use f2_FlowType_HIT
  use f2_FlowType_AddedNew
  use f2_Variables,only: mb1
  implicit none
  private
  
  public:: InitVelocity, InitStatVar, clcStat, add_FlowType_Forcing
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! SetBC_and_UpdateHalo
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitVelocity(ux,uy,uz,Deviation,InputPrmFile); implicit none
    character(len=*),intent(in)::InputPrmFile
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
      
    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call InitVelocity_CH(ux,uy,uz,Deviation,InputPrmFile)
    case(FT_Duct)
      call InitVelocity_Duct(ux,uy,uz,Deviation,InputPrmFile)
    case(FT_TG) ! Taylor-Green vortex
      call InitVelocity_TG(ux,uy,uz)         
    case(FT_HIT) ! Homogenerous isotropic turbulence
      call InitVelocity_HIT(ux,uy,uz,Deviation,InputPrmFile)    
    case(FT_AN) ! Added new
      call InitVelocity_AN(ux,uy,uz)  
    end select    
  end subroutine InitVelocity

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitStatVar
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitStatVar(chFile); implicit none
    character(len=*),intent(in)::chFile

    ! locals
    integer::iUnit, ierr
    character(len=:),allocatable::filename
    
    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call InitStatVar_CH(chFile)
    case(FT_Duct)
      call InitStatVar_Duct()
    case(FT_TG) ! Taylor-Green vortex
      call InitStatVar_TG()          
    case(FT_HIT) ! Homogenerous isotropic turbulence
      call InitStatVar_HIT(chFile)     
    case(FT_AN) ! Added new
          
    end select

    if(nrank==0) then
      if(IsUxConst) then
        filename = strip(ResultsDir_) // 'PrGrad' // int2str(ilast,10) // '_' //strip(RunName_)
        open(newunit=iUnit,file=filename,status='replace',form='formatted',IOSTAT=ierr)
        if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar","Cannot open file: "//filename)
        close(iUnit,IOSTAT=ierr)
      endif
    endif
  end subroutine InitStatVar

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clcStat
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90 
  subroutine clcStat(ux,uy,uz,pressure,ArrTemp1,ArrTemp2); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::ArrTemp1,ArrTemp2
    
    ! locals
    integer::iUnit, ierr
    character(len=:),allocatable::filename

    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call clcStat_CH(ux,uy,uz,pressure,ArrTemp1,ArrTemp2)
    case(FT_Duct)
      call clcStat_Duct(ux,uy,uz, pressure)
    case(FT_TG)  ! Taylor-Green vortex
      call clcStat_TG(ux,uy,uz)          
    case(FT_HIT) ! Homogenerous isotropic turbulence
      call clcStat_HIT(ux,uy,uz,pressure,ArrTemp1)
    case(FT_AN)  ! Added new
          
    end select
    
    if(nrank==0 .and. IsUxConst) then
      filename = strip(ResultsDir_) // 'PrGrad' // int2str(ilast,10) // '_' //strip(RunName_)
      open(newunit=iUnit,file=filename,status='old',position='append',form='formatted',IOSTAT=ierr)
      if(ierr/=0) then
        call MainLog%CheckForError(ErrT_Pass,"clcStat","Cannot open file: "//filename)
      else
        write(iUnit,'(I7,2ES24.15)')itime,SimTime,PrGradData(1)
      endif
      close(iUnit,IOSTAT=ierr)
    endif    
  end subroutine clcStat

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! add_FlowType_Forcing
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine add_FlowType_Forcing(RhsX,RhsY,RhsZ); implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(inout)::RhsX,RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::RhsZ

    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      
    case(FT_Duct)
      
    case(FT_TG) ! Taylor-Green vortex
              
    case(FT_HIT) ! Homogenerous isotropic turbulence
      call add_FlowType_Forcing_HIT(RhsX,RhsY,RhsZ)
    case(FT_AN) ! Added new
          
    end select
      
  end subroutine add_FlowType_Forcing
  
end module f2_FlowCase
