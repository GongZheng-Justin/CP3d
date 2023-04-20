module m_FlowCase
  use m_TypeDef
  use m_Decomp2d
  use m_Parameters
  use m_FlowType_Channel
  use m_FlowType_TGVortex
  use m_FlowType_AddedNew
  use m_Variables,only: mb1
  implicit none
  private
  
  public:: InitVelocity, InitStatVar, clcStat
contains

  !******************************************************************
  ! SetBC_and_UpdateHalo
  !******************************************************************
  subroutine InitVelocity(ux,uy,uz,Deviation,ChannelPrm)
    implicit none
    character(*),intent(in)::ChannelPrm
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
      
    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call InitVelocity_CH(ux,uy,uz,Deviation,ChannelPrm)
    case(FT_TG) ! Taylor-Green vortex
      call InitVelocity_TG(ux,uy,uz,Deviation)          
    case(FT_HI) ! Homogenerous isotropic turbulence
          
    case(FT_AN) ! Added new
      call InitVelocity_AN(ux,uy,uz,Deviation)          
    end select    
  end subroutine InitVelocity

  !******************************************************************
  ! InitStatVar
  !******************************************************************
  subroutine InitStatVar(chFile)
    implicit none
    character(*),intent(in)::chFile

    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call InitStatVar_CH(chFile)
    case(FT_TG) ! Taylor-Green vortex
      call InitStatVar_TG()          
    case(FT_HI) ! Homogenerous isotropic turbulence
          
    case(FT_AN) ! Added new
          
    end select

  end subroutine InitStatVar

  !******************************************************************
  ! clcStat
  !****************************************************************** 
  subroutine clcStat(ux,uy,uz,pressure,ArrTemp1,ArrTemp2)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
    real(RK),dimension(y1size(1),y1size(2),y1size(3)),intent(out)::ArrTemp1,ArrTemp2

    select case(FlowType)
    case(FT_CH,FT_HC) ! Channel
      call clcStat_CH(ux,uy,uz,pressure,ArrTemp1,ArrTemp2)
    case(FT_TG) ! Taylor-Green vortex
      call clcStat_TG(ux,uy,uz,pressure)          
    case(FT_HI) ! Homogenerous isotropic turbulence
          
    case(FT_AN) ! Added new
          
    end select
  end subroutine clcStat 
    
end module m_FlowCase
