#include "definitions_inc.f90"
module f2_FlowType_TGVortex
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_Decomp2d
  use f2_Parameters
  use f2_MeshAndMetries
  use f2_Variables,only: mb1
  use f2_Tools,only:CalcDissipationRate
  implicit none
  private

  public:: InitVelocity_TG, InitStatVar_TG, clcStat_TG
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitVelocity_TG
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitVelocity_TG(ux,uy,uz); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
  
    ! locals
    integer :: ic,jc,kc
    real(RK):: VelRef,LenRef,xpt,ypt,zpt,xct,yct,zct
      
    VelRef=1.0_RK
    LenRef=1.0_RK
    do kc=y1start(3),y1end(3)
      zpt=real(kc-1,kind=RK)*dz
      zct=zpt+0.5_RK*dz
      do jc=y1start(2),y1end(2)
        ypt=yp(jc)
        yct=yc(jc)
        do ic=y1start(1),y1end(1)
          xpt=real(ic-1,kind=RK)*dx
          xct=xpt+0.5_RK*dx
          ux(ic,jc,kc) =  VelRef*sin(xpt/LenRef)*cos(yct/LenRef)*cos(zct/LenRef)              
          uy(ic,jc,kc) = -VelRef*cos(xct/LenRef)*sin(ypt/LenRef)*cos(zct/LenRef)
          uz(ic,jc,kc) =  0.0_RK !VelRef*cos(xct/LenRef)*cos(yct/LenRef)*sin(zpt/LenRef)
        enddo
      enddo
    enddo    
  end subroutine InitVelocity_TG

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! InitStatVar_TG
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine InitStatVar_TG(); implicit none

    ! locals
    integer:: ierr,iUnit
    character(len=:),allocatable::filename

    if(nrank/=0) return
    filename = strip(ResultsDir_) // "TG_dissp"// int2str(ilast,10)
    open(newunit=iUnit, file=filename,status='replace',form='formatted',IOSTAT=ierr)
    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_TG","Cannot open file: "//filename)
    close(iUnit,IOSTAT=ierr)
  end subroutine InitStatVar_TG

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clcStat_TG
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine clcStat_TG(ux,uy,uz); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz
   
    ! locals
    real(Rk):: sum_dissp,sumr
    integer::ic,jc,kc,iUnit,ierr
    character(len=:),allocatable::filename
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::dissp
 
    sumr= 0.0_RK
    call CalcDissipationRate(ux,uy,uz,dissp)
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          sumr=sumr+ dissp(ic,jc,kc)
        enddo
      enddo
    enddo
    call MPI_REDUCE(sumr,sum_dissp,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if(nrank==0) then
      filename = strip(ResultsDir_) // "TG_dissp" // int2str(ilast,10)
      open(newunit=iUnit, file=filename, status='old',position='append',form='formatted',IOSTAT=ierr )
      IF(ierr/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_TG","Cannot open file: "//filename)
      ELSE
        write(iUnit,'(2ES24.15)')SimTime,xnu*sum_dissp/real(nxc*nyc*nzc,RK)
      ENDIF
      close(iUnit,IOSTAT=ierr)
    endif
  end subroutine clcStat_TG
    
end module f2_FlowType_TGVortex
