#include "definitions_inc.f90"
module f2_TScheme
  use MPI
  use mc_LogInfo
  use mc_TypeDef
  use mc_Decomp2d
  use f2_Parameters
  use f2_MeshAndMetries
  use f2_Variables,only: mb1,OutFlowInfoX,OutFlowInfoY
#if defined CFDDEM || defined CFDLPT_TwoWay
  use f2_Variables,only: FpForce_x,FpForce_y,FpForce_z
#endif
  use f2_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal,InversePTriFixedCoe
  implicit none
  private
  
  ! uy/uz Laplacian metries in x-dir (for Crank-Nicolson scheme purpose)
  real(RK),allocatable,dimension(:)::am1cForCN,ap1cForCN
  
  ! ux/uz Laplacian metries in y-dir (for Crank-Nicolson scheme purpose)
  real(RK),allocatable,dimension(:)::am2cForCN,ap2cForCN  
  
  ! ux/uy Laplacian metries in z-dir (for Crank-Nicolson scheme purpose)
  real(RK),allocatable,dimension(:)::am3cForCN,ap3cForCN
  
  procedure(),pointer,public::clcRhsX, clcRhsY, clcRhsZ
  procedure(),pointer,public::clcU1Hat,clcU2Hat,clcU3Hat
  procedure(),pointer,public::clcPrSrc,PressureUpdate
  public:: InitTimeScheme,FluidVelUpdate,clcU1Hat_FEXP,clcU2Hat_FEXP,clcU3Hat_FEXP
contains    

#include "f2_TSchemeFEXP_inc.f90"
#include "f2_TSchemePIMP_inc.f90"
#include "f2_TSchemeFIMP_inc.f90"
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitTimeScheme
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitTimeScheme(); implicit none
    
    ! locals
    integer::ierr

    ! uy/uz Laplacian metries in x-dir (for Crank-Nicolson scheme purpose)
    allocate(ap1cForCN(1:nxc), am1cForCN(1:nxc),Stat=ierr)
    ap1cForCN=ap1c;   am1cForCN=am1c;
    if(BcOption(xm_dir)==BC_NoSlip) then
      am1cForCN(1)= 2.0_RK*am1c(1)
    elseif(BcOption(xm_dir)==BC_FreeSlip) then
      am1cForCN(1)= 0.0_RK      
    endif
    if(BcOption(xp_dir)==BC_NoSlip) then
      ap1cForCN(nxc)= 2.0_RK*ap1c(nxc)
    elseif(BcOption(xp_dir)==BC_FreeSlip) then
      ap1cForCN(nxc)= 0.0_RK    
    endif
        
    ! ux/uz Laplacian metries in y-dir (for Crank-Nicolson scheme purpose)
    allocate(ap2cForCN(1:nyc), am2cForCN(1:nyc),Stat=ierr)
    ap2cForCN = ap2c;   am2cForCN = am2c;
    if(BcOption(ym_dir)==BC_NoSlip) then
      am2cForCN(1)= 2.0_RK*am2c(1)    
    elseif(BcOption(ym_dir)==BC_FreeSlip) then
      am2cForCN(1)= 0.0_RK    
    endif
    if(BcOption(yp_dir)==BC_NoSlip) then
      ap2cForCN(nyc)= 2.0_RK*ap2c(nyc)
    elseif(BcOption(yp_dir)==BC_FreeSlip) then
      ap2cForCN(nyc)= 0.0_RK    
    endif
    
    ! ux/uy Laplacian metries in z-dir (for Crank-Nicolson scheme purpose)
    allocate(ap3cForCN(1:nzc), am3cForCN(1:nzc),Stat=ierr)
    ap3cForCN = ap3c;   am3cForCN = am3c;
    if(BcOption(zm_dir)==BC_NoSlip) then
      am3cForCN(1)= 2.0_RK*am3c(1)    
    elseif(BcOption(zm_dir)==BC_FreeSlip) then
      am3cForCN(1)= 0.0_RK
    endif
    if(BcOption(zp_dir)==BC_NoSlip) then
      ap3cForCN(nzc)= 2.0_RK*ap3c(nzc)
    elseif(BcOption(zp_dir)==BC_FreeSlip) then
      ap3cForCN(nzc)= 0.0_RK
    endif
    
    ! FEXP 0: full explicit 
    ! PIMP 1: partial implicit, only use C-N in y-dir 
    ! FIMP 2: full implicit, use C-N in all 3 dirs.
    if((BcOption(xm_dir)==BC_Period .and. BcOption(xp_dir)/=BC_Period) .or. &
       (BcOption(ym_dir)==BC_Period .and. BcOption(yp_dir)/=BC_Period) .or. &
       (BcOption(zm_dir)==BC_Period .and. BcOption(zp_dir)/=BC_Period)) then
      call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Periodic Bc Wrong")
    endif
    if( BcOption(xm_dir)==BC_OutFlow .or. BcOption(ym_dir)==BC_OutFlow .or. &
        BcOption(zm_dir)==BC_OutFlow .or. BcOption(zp_dir)==BC_OutFlow)     &
      call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","OutFlow Is ONLY Supported in xp-dir and yp-dir")

    SELECT CASE(IsImplicit)
    CASE(0)
      clcRhsX => clcRhsX_FEXP
      clcRhsY => clcRhsY_FEXP
      clcRhsZ => clcRhsZ_FEXP
      clcU1Hat=> clcU1Hat_FEXP
      clcU2Hat=> clcU2Hat_FEXP
      clcU3Hat=> clcU3Hat_FEXP
      clcPrSrc=> clcPrSrcOther
      PressureUpdate => PressureUpdate_FEXP
    CASE(1)
      clcRhsX => clcRhsX_PIMP
      clcRhsY => clcRhsY_PIMP
      clcRhsZ => clcRhsZ_PIMP
      if(BcOption(ym_dir)==BC_Period) then
        clcU1Hat=> clcU1Hat_PIMP_0
        clcU2Hat=> clcU2Hat_PIMP_0
        clcU3Hat=> clcU3Hat_PIMP_0
      else
        clcU1Hat=> clcU1Hat_PIMP
        clcU2Hat=> clcU2Hat_PIMP
        clcU3Hat=> clcU3Hat_PIMP
      endif
      clcPrSrc=> clcPrSrcOther
      PressureUpdate => PressureUpdate_PIMP
    CASE(2)
      clcRhsX => clcRhsX_FIMP
      clcRhsY => clcRhsY_FIMP
      clcRhsZ => clcRhsZ_FIMP
      clcU1Hat=> clcU1Hat_FIMP
      clcU2Hat=> clcU2Hat_FIMP
      clcU3Hat=> clcU3Hat_FIMP
      clcPrSrc=> clcPrSrc_FIMP
      PressureUpdate => PressureUpdate_FIMP
    END SELECT
  end subroutine InitTimeScheme

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clcPrSrcOther
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90  
  subroutine clcPrSrcOther(ux,uy,uz,prsrc,pressure,divmax); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
    real(RK),intent(out)::divmax
    
    ! locals
    integer::ic,jc,kc,ip,jp,kp,ierr
    real(RK)::sucaj,rdiv,divmax1

    divmax1=0.0_RK
    DO kc=y1start(3),y1end(3)
       kp=kc+1
       do jc=y1start(2),y1end(2)
         jp=jc+1
         sucaj=rdyp(jc)
         do ic=y1start(1),y1end(1)
           ip=ic+1
           rdiv= (ux(ip,jc,kc)-ux(ic,jc,kc))*rdx + (uy(ic,jp,kc)-uy(ic,jc,kc))*sucaj + &
                 (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
           divmax1=max(abs(rdiv),divmax1)
           prsrc(ic,jc,kc)= rdiv
         enddo
       enddo
     ENDDO
     call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  end subroutine clcPrSrcOther    

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clcPrSrc_FIMP
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90  
  subroutine clcPrSrc_FIMP(ux,uy,uz,prsrc,pressure,divmax); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
    real(RK),intent(out)::divmax
    
    ! locals
    integer::ic,jc,kc,ip,jp,kp,ierr
    real(RK)::sucaj,rdiv,divmax1,xnuhm

    divmax1=0.0_RK
    xnuhm= -0.5_RK*xnu
    DO kc=y1start(3),y1end(3)
       kp=kc+1
       do jc=y1start(2),y1end(2)
         jp=jc+1
         sucaj=rdyp(jc)
         do ic=y1start(1),y1end(1)
           ip=ic+1
           rdiv= (ux(ip,jc,kc)-ux(ic,jc,kc))*rdx + (uy(ic,jp,kc)-uy(ic,jc,kc))*sucaj + &
                 (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz
           divmax1=max(abs(rdiv),divmax1)
           prsrc(ic,jc,kc)= rdiv
           pressure(ic,jc,kc)= pressure(ic,jc,kc)+ xnuhm*rdiv         ! new added for full implicit scheme
         enddo
       enddo
     ENDDO
     call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  end subroutine clcPrSrc_FIMP
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! FluidVelUpdate
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine FluidVelUpdate(prphiHalo,ux,uy,uz); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out):: ux,uy,uz
    
    ! locals
    integer::ic,jc,kc,im,jm,km
    real(RK)::sucac,locphi

    DO kc=y1start(3),y1end(3)
      km=kc-1
      do jc=y1start(2),y1end(2)
        jm=jc-1
        sucac=rdyc(jc)
        do ic=y1start(1),y1end(1)
          im=ic-1
          locphi=prphiHalo(ic,jc,kc)
          ux(ic,jc,kc)=ux(ic,jc,kc) -(locphi-prphiHalo(im,jc,kc))*rdx
          uy(ic,jc,kc)=uy(ic,jc,kc) -(locphi-prphiHalo(ic,jm,kc))*sucac
          uz(ic,jc,kc)=uz(ic,jc,kc) -(locphi-prphiHalo(ic,jc,km))*rdz
        enddo
      enddo
    ENDDO
  end subroutine FluidVelUpdate

end module f2_TScheme
