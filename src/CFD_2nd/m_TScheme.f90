module m_TScheme
  use MPI
  use m_LogInfo
  use m_TypeDef
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1,nut,OutFlowInfoX,OutFlowInfoY
#if defined CFDDEM || defined CFDLPT_TwoWay
  use m_Variables,only: FpForce_x,FpForce_y,FpForce_z
#endif
  use m_Tools,only: InverseTridiagonal,InversePeriodicTridiagonal,InversePTriFixedCoe
  implicit none
  private 

  logical,public::IsUxConst
  real(RK)::ubulk,ybulk1,ybulk2
  
  ! uy/uz Laplacian metries in x-dir (for Crank-Nicolson scheme purpose)
  real(RK),allocatable,dimension(:)::am1cForCN,ap1cForCN
  
  ! ux/uz Laplacian metries in y-dir (for Crank-Nicolson scheme purpose)
  real(RK),allocatable,dimension(:)::am2cForCN,ap2cForCN  
  
  ! ux/uy Laplacian metries in z-dir (for Crank-Nicolson scheme purpose)
  real(RK),allocatable,dimension(:)::am3cForCN,ap3cForCN
  
  procedure(),pointer,public::clcRhsX, clcRhsY, clcRhsZ
  procedure(),pointer,public::clcU1Hat,clcU2Hat,clcU3Hat
  procedure(),pointer,public::clcPrSrc,PressureUpdate
  public:: InitTimeScheme,FluidVelUpdate
contains    

#include "m_TSchemeFEXP_inc.f90"
#include "m_TSchemePIMP_inc.f90"
#include "m_TSchemeFIMP_inc.f90"
  !******************************************************************
  ! InitTimeScheme
  !******************************************************************
  subroutine InitTimeScheme(chFile)
    implicit none
    character(*),intent(in)::chFile
    
    ! locals
    integer::nUnit,ierror
    logical::IsUbulkGlobal
    NAMELIST/ubulk_Param/IsUxConst,IsUbulkGlobal,ubulk,ybulk1,ybulk2

    open(newunit=nUnit, file=chFile, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme", "Cannot open file: "//trim(chFile))
    read(nUnit,nml=ubulk_Param)
    close(nUnit,IOSTAT=ierror)
    if(IsUbulkGlobal) then
      ybulk1=zero; ybulk2=yly
    endif
    if(IsUxConst .and. ybulk2<=ybulk1)  call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme", "ybulk1 OR ybulk2 wrong")
    if(nrank==0)write(MainLog%nUnit,nml=ubulk_Param)    

    ! uy/uz Laplacian metries in x-dir (for Crank-Nicolson scheme purpose)
    allocate(ap1cForCN(1:nxc),am1cForCN(1:nxc),Stat=ierror)
    ap1cForCN=ap1c;   am1cForCN=am1c;
    if(BcOption(xm_dir)==BC_NoSlip) then
      am1cForCN(1)= two*am1c(1)
    elseif(BcOption(xm_dir)==BC_FreeSlip) then
      am1cForCN(1)= zero      
    endif
    if(BcOption(xp_dir)==BC_NoSlip) then
      ap1cForCN(nxc)= two*ap1c(nxc)
    elseif(BcOption(xp_dir)==BC_FreeSlip) then
      ap1cForCN(nxc)= zero    
    endif
        
    ! ux/uz Laplacian metries in y-dir (for Crank-Nicolson scheme purpose)
    allocate(ap2cForCN(1:nyc),am2cForCN(1:nyc),Stat=ierror)
    ap2cForCN = ap2c;   am2cForCN = am2c;
    if(BcOption(ym_dir)==BC_NoSlip) then
      am2cForCN(1)= two*am2c(1)    
    elseif(BcOption(ym_dir)==BC_FreeSlip) then
      am2cForCN(1)= zero    
    endif
    if(BcOption(yp_dir)==BC_NoSlip) then
      ap2cForCN(nyc)= two*ap2c(nyc)
    elseif(BcOption(yp_dir)==BC_FreeSlip) then
      ap2cForCN(nyc)= zero    
    endif
    
    ! ux/uy Laplacian metries in z-dir (for Crank-Nicolson scheme purpose)
    allocate(ap3cForCN(1:nzc),am3cForCN(1:nzc),Stat=ierror)
    ap3cForCN = ap3c;   am3cForCN = am3c;
    if(BcOption(zm_dir)==BC_NoSlip) then
      am3cForCN(1)= two*am3c(1)    
    elseif(BcOption(zm_dir)==BC_FreeSlip) then
      am3cForCN(1)= zero
    endif
    if(BcOption(zp_dir)==BC_NoSlip) then
      ap3cForCN(nzc)= two*ap3c(nzc)
    elseif(BcOption(zp_dir)==BC_FreeSlip) then
      ap3cForCN(nzc)= zero
    endif
    
    ! FEXP 0: full explicit
    ! PIMP 1: partial implicit, only use C-N in y-dir 
    ! FIMP 2: full implicit, use C-N in all 3 dirs.
    if( (BcOption(xm_dir)==BC_PERIOD .and. BcOption(xp_dir)/=BC_PERIOD) .or.  (BcOption(ym_dir)==BC_PERIOD .and. BcOption(yp_dir)/=BC_PERIOD) .or. &
        (BcOption(zm_dir)==BC_PERIOD .and. BcOption(zp_dir)/=BC_PERIOD)   ) then
      call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Periodic Bc Wrong ")
    endif
    if( BcOption(xm_dir)==BC_OutFlow .or. BcOption(ym_dir)==BC_OutFlow .or. BcOption(zm_dir)==BC_OutFlow &
    .or.BcOption(zp_dir)==BC_OutFlow) call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","OutFlow Is ONLY Supported in xp-dir and yp-dir")

    ! Note here, if use full implicit time scheme, the potential periodic bc should be in x-dir first,a nd then z-dir, then y-dir
    if(IsImplicit==2) then
      if(BcOption(xm_dir) /=BC_PERIOD .and. (BcOption(ym_dir)==BC_PERIOD .or. BcOption(zm_dir)==BC_PERIOD) ) then
        call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Bc type wrong 1 ")
      endif
      if(BcOption(zm_dir) /=BC_PERIOD .and. BcOption(ym_dir) ==BC_PERIOD ) then
        call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","Bc type wrong 2 ")
      endif
    endif

    if(LES_type==0) then
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
        if(BcOption(ym_dir)==BC_PERIOD) then
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
        if(BcOption(xm_dir) ==BC_PERIOD) then ! There can be several periodic Bcs
          if(BcOption(zm_dir) ==BC_PERIOD) then
            if(BcOption(ym_dir) ==BC_PERIOD) then !
              clcU1Hat=> clcU1Hat_FIMP_000
              clcU2Hat=> clcU2Hat_FIMP_000 
              clcU3Hat=> clcU3Hat_FIMP_000    
            else
              clcU1Hat=> clcU1Hat_FIMP_010
              clcU2Hat=> clcU2Hat_FIMP_010
              clcU3Hat=> clcU3Hat_FIMP_010
            endif
          else
            clcU1Hat=> clcU1Hat_FIMP_011
            clcU2Hat=> clcU2Hat_FIMP_011
            clcU3Hat=> clcU3Hat_FIMP_011
          endif
        else                                  ! NO periodic Bc exist
          clcU1Hat=> clcU1Hat_FIMP_111
          clcU2Hat=> clcU2Hat_FIMP_111
          clcU3Hat=> clcU3Hat_FIMP_111
        endif
        clcPrSrc=> clcPrSrc_FIMP
        PressureUpdate => PressureUpdate_FIMP
      END SELECT
    else
      if(BcOption(xp_dir)==BC_OutFlow .or. BcOption(yp_dir)==BC_OutFlow) then
        call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","OutFlow+LES Have not finished yet !!!")
      endif
      SELECT CASE(IsImplicit)
      CASE(0)
        clcRhsX => clcRhsX_FEXP_LES
        clcRhsY => clcRhsY_FEXP_LES
        clcRhsZ => clcRhsZ_FEXP_LES
        clcU1Hat=> clcU1Hat_FEXP_LES
        clcU2Hat=> clcU2Hat_FEXP_LES
        clcU3Hat=> clcU3Hat_FEXP_LES
        clcPrSrc=> clcPrSrcOther
        PressureUpdate => PressureUpdate_FEXP_LES
      CASE(1)
        clcRhsX => clcRhsX_PIMP_LES
        clcRhsY => clcRhsY_PIMP_LES
        clcRhsZ => clcRhsZ_PIMP_LES
        if(BcOption(ym_dir)==BC_PERIOD) then
          clcU1Hat=> clcU1Hat_PIMP_LES_0
          clcU2Hat=> clcU2Hat_PIMP_LES_0
          clcU3Hat=> clcU3Hat_PIMP_LES_0
        else
          clcU1Hat=> clcU1Hat_PIMP_LES
          clcU2Hat=> clcU2Hat_PIMP_LES
          clcU3Hat=> clcU3Hat_PIMP_LES
        endif
        clcPrSrc=> clcPrSrcOther
        PressureUpdate => PressureUpdate_PIMP_LES
      CASE(2)
        call MainLog%CheckForError(ErrT_Abort,"InitTimeScheme","FIMP+LES Have not finished yet !!!")
      END SELECT
    endif
  end subroutine InitTimeScheme

  !******************************************************************
  ! clcPrSrcOther
  !******************************************************************  
  subroutine clcPrSrcOther(ux,uy,uz,prsrc,pressure,divmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
    real(RK),intent(out)::divmax
    
    ! locals
    integer::ic,jc,kc,ip,jp,kp,ierror
    real(RK)::sudtal,sucaj,rdiv,divmax1

    divmax1=zero
    sudtal=one/pmAlpha
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
           prsrc(ic,jc,kc)= sudtal * rdiv
         enddo
       enddo
     ENDDO
     call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)
  end subroutine clcPrSrcOther    

  !******************************************************************
  ! clcPrSrc_FIMP
  !******************************************************************  
  subroutine clcPrSrc_FIMP(ux,uy,uz,prsrc,pressure,divmax)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::prsrc
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::pressure
    real(RK),intent(out)::divmax
    
    ! locals
    integer::ic,jc,kc,ip,jp,kp,ierror
    real(RK)::sudtal,sucaj,rdiv,divmax1,xnuhm

    divmax1=zero
    xnuhm= -half*xnu
    sudtal=one/pmAlpha
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
           prsrc(ic,jc,kc)= sudtal * rdiv
           pressure(ic,jc,kc)= pressure(ic,jc,kc)+ xnuhm*rdiv         ! new added for full implicit scheme
         enddo
       enddo
     ENDDO
     call MPI_REDUCE(divmax1,divmax,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,ierror)
  end subroutine clcPrSrc_FIMP
  
  !******************************************************************
  ! FluidVelUpdate
  !******************************************************************
  subroutine FluidVelUpdate(prphiHalo,ux,uy,uz)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in):: prphiHalo
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out):: ux,uy,uz
    
    ! locals
    integer::ic,jc,kc,im,jm,km
    real(RK)::rdxEta2,sucacEta2,rdzEta2,locphi

    rdxEta2=rdx*pmAlpha
    rdzEta2=rdz*pmAlpha
    DO kc=y1start(3),y1end(3)
      km=kc-1
      do jc=y1start(2),y1end(2)
        jm=jc-1
        sucacEta2=rdyc(jc)*pmAlpha
        do ic=y1start(1),y1end(1)
          im=ic-1
          locphi=prphiHalo(ic,jc,kc)
          ux(ic,jc,kc)=ux(ic,jc,kc)-  (locphi-prphiHalo(im,jc,kc))*rdxEta2
          uy(ic,jc,kc)=uy(ic,jc,kc)-  (locphi-prphiHalo(ic,jm,kc))*sucacEta2
          uz(ic,jc,kc)=uz(ic,jc,kc)-  (locphi-prphiHalo(ic,jc,km))*rdzEta2                
        enddo
      enddo
    ENDDO
  end subroutine FluidVelUpdate

end module m_TScheme
