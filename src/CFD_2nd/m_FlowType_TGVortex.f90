module m_FlowType_TGVortex
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only: mb1
  use m_Tools,only:CalcDissipationRate
  implicit none
  private

  public:: InitVelocity_TG, Update_uy_ym_TG
  public:: InitStatVar_TG,  clcStat_TG
contains

  !******************************************************************
  ! InitVelocity_TG
  !******************************************************************
  subroutine InitVelocity_TG(ux,uy,uz,Deviation)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1), y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
  
    ! locals
    integer :: ic,jc,kc
    real(RK):: VelRef,LenRef,xpt,ypt,zpt,xct,yct,zct
      
    VelRef=one
    LenRef=one
    do kc=y1start(3),y1end(3)
      zpt=real(kc-1,kind=RK)*dz
      zct=zpt+half*dz
      do jc=y1start(2),y1end(2)
        ypt=yp(jc)
        yct=yc(jc)
        do ic=y1start(1),y1end(1)
          xpt=real(ic-1,kind=RK)*dx
          xct=xpt+half*dx
          ux(ic,jc,kc) =  VelRef*sin(xpt/LenRef)*cos(yct/LenRef)*cos(zct/LenRef)              
          uy(ic,jc,kc) = -VelRef*cos(xct/LenRef)*sin(ypt/LenRef)*cos(zct/LenRef)
          uz(ic,jc,kc) =  zero !VelRef*cos(xct/LenRef)*cos(yct/LenRef)*sin(zpt/LenRef)
        enddo
      enddo
    enddo
    
  end subroutine InitVelocity_TG

  !******************************************************************
  ! Update_uy_ym_TG
  !******************************************************************   
  subroutine Update_uy_ym_TG(uy_ym, duy_ym, TimeNew)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(inout):: uy_ym,duy_ym    
    real(RK),intent(in):: TimeNew
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = zero
     
    duy_ym = uy_ym - duy_ym
    
  end subroutine Update_uy_ym_TG

  !******************************************************************
  ! InitStatVar_TG
  !******************************************************************
  subroutine InitStatVar_TG()
    implicit none

    ! locals
    integer:: ierror,nUnit
    character(len=128)::filename

    if(nrank/=0) return
    write(filename,'(A,I10.10)')trim(ResultsDir)//"TG_dissp",ilast
    open(newunit=nUnit, file=filename,status='replace',form='formatted',IOSTAT=ierror)
    if(ierror/=0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_TG","Cannot open file: "//trim(filename))
    close(nUnit,IOSTAT=ierror)

  end subroutine InitStatVar_TG

  !******************************************************************
  ! clcStat_TG
  !******************************************************************
  subroutine clcStat_TG(ux,uy,uz,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
   
    ! locals
    real(Rk):: sum_dissp,sumr
    character(len=128)::filename
    integer::ic,jc,kc,nUnit,ierror
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3))::dissp
 
    sumr= zero
    call CalcDissipationRate(ux,uy,uz,dissp)
    do kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
          sumr=sumr+ dissp(ic,jc,kc)
        enddo
      enddo
    enddo
    call MPI_REDUCE(sumr,sum_dissp,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)

    if(nrank==0) then
      write(filename,'(A,I10.10)')trim(ResultsDir)//"TG_dissp",ilast
      open(newunit=nUnit, file=filename, status='old',position='append',form='formatted',IOSTAT=ierror )
      IF(ierror/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_TG","Cannot open file: "//trim(filename))
      ELSE
        write(nUnit,'(2ES24.15)')SimTime,xnu*sum_dissp/real(nxc*nyc*nzc,RK)
      ENDIF
      close(nUnit,IOSTAT=ierror)
    endif
  end subroutine clcStat_TG
    
end module m_FlowType_TGVortex
