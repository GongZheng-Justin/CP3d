module ca_BC_and_Halo
  use MPI
  use m_TypeDef
  use m_Parameters
  use Prtcl_Decomp_2d
  use m_MeshAndMetries
  use m_Variables,only:mb1,hi1
  use m_Decomp2d,only:nrank,y1start,y1end,HaloInfo,myProcNghBC,myUpdate_Halo
  implicit none
  private

  public:: SetBC_and_UpdateHalo_VelIBM,Gather_Halo_IBMForce
contains
  !******************************************************************
  ! SetBC_and_UpdateHalo_VelIBM
  !******************************************************************   
  subroutine SetBC_and_UpdateHalo_VelIBM(ux,uy,uz,uy_ym)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(in)::uy_ym
    
    ! locals
    integer::ic,kc
    type(HaloInfo):: hi_ux_interp,hi_uz_interp

   ! yp-dir
    SELECT CASE(BcOption(yp_dir))
    CASE(BC_NoSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,nyp,  kc)= uxBcValue(yp_dir)*2.0_RK -ux(ic, nyc, kc)
          uy(ic,nyp,  kc)= uyBcValue(yp_dir)
          uy(ic,nyp+1,kc)= uyBcValue(yp_dir)*2.0_RK -uy(ic, nyc, kc)
          uz(ic,nyp,  kc)= uzBcValue(yp_dir)*2.0_RK -uz(ic, nyc, kc)
        enddo
      enddo 
    CASE(BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic,nyp,  kc)= uxBcValue(yp_dir)*dyp(nyp) +ux(ic, nyc, kc)
          uy(ic,nyp,  kc)= uyBcValue(yp_dir)
          uy(ic,nyp+1,kc)= uyBcValue(yp_dir)*2.0_RK      -uy(ic, nyc, kc)
          uz(ic,nyp,  kc)= uzBcValue(yp_dir)*dyp(nyp) +uz(ic, nyc, kc)
        enddo
      enddo     
    END SELECT 

    ! ym-dir
    SELECT CASE(BcOption(ym_dir))       
    CASE(BC_NoSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic, 0, kc) = uxBcValue(ym_dir)*2.0_RK-ux(ic, 1, kc)
          uy(ic, 1, kc) = uy_ym(ic,kc)
          uy(ic, 0, kc) = uy_ym(ic,kc)*2.0_RK -uy(ic,2,kc)
          uz(ic, 0, kc) = uzBcValue(ym_dir)*2.0_RK-uz(ic, 1, kc)
        enddo
      enddo
    CASE(BC_FreeSlip)
      do kc=y1start(3),y1end(3)
        do ic=y1start(1),y1end(1)
          ux(ic, 0, kc) = ux(ic, 1, kc)-uxBcValue(ym_dir)*dyp(1)
          uy(ic, 1, kc) = uy_ym(ic,kc)
          uy(ic, 0, kc) = uy_ym(ic,kc)*2.0_RK -uy(ic,2,kc)
          uz(ic, 0, kc) = uz(ic, 1, kc)-uzBcValue(ym_dir)*dyp(1)
        enddo
      enddo    
    END SELECT 

    ! xp-dir
    SELECT CASE(myProcNghBC(y_pencil, 3))
    CASE(BC_NoSlip)
      ux(nxp,  0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)
      ux(nxp+1,0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)*2.0_RK -ux(nxc, 0:nyp,y1start(3):y1end(3))
      uy(nxp,  0:nyp,y1start(3):y1end(3))= uyBcValue(xp_dir)*2.0_RK -uy(nxc, 0:nyp,y1start(3):y1end(3))
      uz(nxp,  0:nyp,y1start(3):y1end(3))= uzBcValue(xp_dir)*2.0_RK -uz(nxc, 0:nyp,y1start(3):y1end(3))
    CASE(BC_FreeSlip)
      ux(nxp,  0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)
      ux(nxp+1,0:nyp,y1start(3):y1end(3))= uxBcValue(xp_dir)*2.0_RK -ux(nxc, 0:nyp,y1start(3):y1end(3))
      uy(nxp,  0:nyp,y1start(3):y1end(3))= uy(nxc, 0:nyp,y1start(3):y1end(3)) +uyBcValue(xp_dir)*dx
      uz(nxp,  0:nyp,y1start(3):y1end(3))= uz(nxc, 0:nyp,y1start(3):y1end(3)) +uzBcValue(xp_dir)*dx        
    END SELECT
    
    ! xm-dir
    SELECT CASE(myProcNghBC(y_pencil, 4))
    CASE(BC_NoSlip)
      ux(1,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)
      ux(0,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)*2.0_RK -ux(2,0:nyp,y1start(3):y1end(3)) 
      uy(0,0:nyp,y1start(3):y1end(3)) = uyBcValue(xm_dir)*2.0_RK -uy(1,0:nyp,y1start(3):y1end(3))
      uz(0,0:nyp,y1start(3):y1end(3)) = uzBcValue(xm_dir)*2.0_RK -uz(1,0:nyp,y1start(3):y1end(3))       
    CASE(BC_FreeSlip)
      ux(1,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)
      ux(0,0:nyp,y1start(3):y1end(3)) = uxBcValue(xm_dir)*2.0_RK -ux(2,0:nyp,y1start(3):y1end(3)) 
      uy(0,0:nyp,y1start(3):y1end(3)) = uy(1,0:nyp,y1start(3):y1end(3)) -uyBcValue(xm_dir)*dx
      uz(0,0:nyp,y1start(3):y1end(3)) = uz(1,0:nyp,y1start(3):y1end(3)) -uzBcValue(xm_dir)*dx         
    END SELECT    
       
    ! zp-dir        
    SELECT CASE(myProcNghBC(y_pencil, 1))
    CASE(BC_NoSlip) 
      ux(y1start(1):y1end(1), 0:nyp, nzp  ) = uxBcValue(zp_dir)*2.0_RK -ux(y1start(1):y1end(1), 0:nyp, nzc) 
      uy(y1start(1):y1end(1), 0:nyp, nzp  ) = uyBcValue(zp_dir)*2.0_RK -uy(y1start(1):y1end(1), 0:nyp, nzc)
      uz(y1start(1):y1end(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
      uz(y1start(1):y1end(1), 0:nyp, nzp+1) = uzBcValue(zp_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, nzc)
    CASE(BC_FreeSlip)
      ux(y1start(1):y1end(1), 0:nyp, nzp  ) = ux(y1start(1):y1end(1), 0:nyp, nzc) +uxBcValue(zp_dir)*dz 
      uy(y1start(1):y1end(1), 0:nyp, nzp  ) = uy(y1start(1):y1end(1), 0:nyp, nzc) +uyBcValue(zp_dir)*dz 
      uz(y1start(1):y1end(1), 0:nyp, nzp  ) = uzBcValue(zp_dir)
      uz(y1start(1):y1end(1), 0:nyp, nzp+1) = uzBcValue(zp_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, nzc)
    END SELECT         
    
    ! zm-dir        
    SELECT CASE(myProcNghBC(y_pencil, 2))
    CASE(BC_NoSlip)
      ux(y1start(1):y1end(1), 0:nyp, 0) = uxBcValue(zm_dir)*2.0_RK -ux(y1start(1):y1end(1), 0:nyp, 1)
      uy(y1start(1):y1end(1), 0:nyp, 0) = uyBcValue(zm_dir)*2.0_RK -uy(y1start(1):y1end(1), 0:nyp, 1)
      uz(y1start(1):y1end(1), 0:nyp, 1) = uzBcValue(zm_dir)
      uz(y1start(1):y1end(1), 0:nyp, 0) = uzBcValue(zm_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, 2)
    CASE(BC_FreeSlip)
      ux(y1start(1):y1end(1), 0:nyp, 0) = ux(y1start(1):y1end(1), 0:nyp, 1) -uxBcValue(zm_dir)*dz 
      uy(y1start(1):y1end(1), 0:nyp, 0) = uy(y1start(1):y1end(1), 0:nyp, 1) -uyBcValue(zm_dir)*dz 
      uz(y1start(1):y1end(1), 0:nyp, 1) = uzBcValue(zm_dir)
      uz(y1start(1):y1end(1), 0:nyp, 0) = uzBcValue(zm_dir)*2.0_RK -uz(y1start(1):y1end(1), 0:nyp, 2)     
    END SELECT     

    ! update halo
    hi_ux_interp%pencil = y_pencil
    hi_ux_interp%xmh=1;  hi_ux_interp%xph=2
    hi_ux_interp%ymh=0;  hi_ux_interp%yph=0
    hi_ux_interp%zmh=1;  hi_ux_interp%zph=1

    hi_uz_interp%pencil = y_pencil
    hi_uz_interp%xmh=1;  hi_uz_interp%xph=1
    hi_uz_interp%ymh=0;  hi_uz_interp%yph=0
    hi_uz_interp%zmh=1;  hi_uz_interp%zph=2

    call myupdate_halo(ux, mb1, hi_ux_interp)
    call myupdate_halo(uy, mb1, hi1)
    call myupdate_halo(uz, mb1, hi_uz_interp)
  
  end subroutine SetBC_and_UpdateHalo_VelIBM

  !******************************************************************
  ! Gather_Halo_IBMForce
  !******************************************************************
  subroutine Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::VolForce_x,VolForce_y,VolForce_z

    ! locals
    integer:: ixs1,ixs2,iys1,iys2,izs1,izs2                ! index for halo sending
    integer,dimension(MPI_STATUS_SIZE,2) :: SRstatus
    integer:: i,j,k,icount,ierror,requests(2),ProcNgh(4)
    real(RK),dimension(:,:),  allocatable:: Mat2HaloS,Mat2HaloR
    real(RK),dimension(:,:,:),allocatable:: Mat3HaloS,Mat3HaloR

    do j=1,4
      if(myProcNghBC(y_pencil,j)<0) then
        ProcNgh(j)=MPI_PROC_NULL
      else
        ProcNgh(j)=myProcNghBC(y_pencil,j)
      endif
    enddo
#ifdef ReserveIBMForce
    do k=mb1%zmm,mb1%zpm
      do i=mb1%xmm,mb1%xpm
        VolForce_x(i,  1,k)=VolForce_x(i,  1,k)+VolForce_x(i,    0,k);  VolForce_x(i,  0,  k)=0.0_RK
        VolForce_y(i,  1,k)=VolForce_y(i,  1,k)+VolForce_y(i,    0,k);  VolForce_y(i,  0,  k)=0.0_RK
        VolForce_z(i,  1,k)=VolForce_z(i,  1,k)+VolForce_z(i,    0,k);  VolForce_z(i,  0,  k)=0.0_RK
        VolForce_x(i,nyc,k)=VolForce_x(i,nyc,k)+VolForce_x(i,  nyp,k);  VolForce_x(i,nyp,  k)=0.0_RK
        VolForce_y(i,nyp,k)=VolForce_y(i,nyp,k)+VolForce_y(i,nyp+1,k);  VolForce_y(i,nyp+1,k)=0.0_RK
        VolForce_z(i,nyc,k)=VolForce_z(i,nyc,k)+VolForce_z(i,  nyp,k);  VolForce_z(i,nyp,  k)=0.0_RK
      enddo
    enddo
#endif
    if(BcOption(ym_dir)==BC_PERIOD) then
      VolForce_x(:,  1,:)=VolForce_x(:,  1,:)+VolForce_x(:,nyp,:)
      VolForce_x(:,nyc,:)=VolForce_x(:,nyc,:)+VolForce_x(:,  0,:)
      VolForce_z(:,  1,:)=VolForce_z(:,  1,:)+VolForce_z(:,nyp,:)
      VolForce_z(:,nyc,:)=VolForce_z(:,nyc,:)+VolForce_z(:,  0,:)
      VolForce_y(:,  1,:)=VolForce_y(:,  1,:)+VolForce_y(:,nyp,:)
      VolForce_y(:,nyc,:)=VolForce_y(:,nyc,:)+VolForce_y(:,  0,:)
      VolForce_y(:,  2,:)=VolForce_y(:,2,:)+VolForce_y(:,nyp+1,:)
    endif
    
    !===================== Step1: receive from xp_dir, and send to xm_dir =====================!
    iys1=1;             iys2=nyp    
    izs1=y1start(3)-1;  izs2=y1end(3)+2
    icount= (iys2-iys1+1)*(izs2-izs1+1)
    allocate(Mat2HaloS(iys1:iys2,izs1:izs2))
    allocate(Mat2HaloR(iys1:iys2,izs1:izs2))

    ixs1=y1start(1)-1;   ixs2=y1start(1)-1
    ! ux ( receive from xp_dir, and send to xm_dir )
    IF(ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= VolForce_x(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= VolForce_x(ixs1,j,k)
        enddo
      enddo
      call MPI_IRECV( Mat2HaloR,icount,real_type,ProcNgh(3),1,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat2HaloS,icount,real_type,ProcNgh(4),1,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(3)>=0) THEN
      i=y1end(1)
      do k=izs1,izs2
        do j=iys1,iys2
          VolForce_x(i,j,k)=VolForce_x(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ! uy ( receive from xp_dir, and send to xm_dir )
    IF(ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= VolForce_y(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= VolForce_y(ixs1,j,k)
        enddo
      enddo
      call MPI_IRECV( Mat2HaloR,icount,real_type,ProcNgh(3),2,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat2HaloS,icount,real_type,ProcNgh(4),2,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(3)>=0) THEN
      i=y1end(1)
      do k=izs1,izs2
        do j=iys1,iys2
          VolForce_y(i,j,k)=VolForce_y(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ! uz ( receive from xp_dir, and send to xm_dir )
    IF(ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= VolForce_z(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= VolForce_z(ixs1,j,k)
        enddo
      enddo
      call MPI_IRECV( Mat2HaloR,icount,real_type,ProcNgh(3),3,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat2HaloS,icount,real_type,ProcNgh(4),3,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(3)>=0) THEN
      i=y1end(1)
      do k=izs1,izs2
        do j=iys1,iys2
          VolForce_z(i,j,k)=VolForce_z(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ixs1=y1end(1)+1;   ixs2=y1end(1)+1
    ! uy ( receive from xm_dir, and send to xp_dir )
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= VolForce_y(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= VolForce_y(ixs1,j,k)
        enddo
      enddo
      call MPI_IRECV( Mat2HaloR,icount,real_type,ProcNgh(4),4,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat2HaloS,icount,real_type,ProcNgh(3),4,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(4)>=0) THEN
      i=y1start(1)
      do k=izs1,izs2
        do j=iys1,iys2
          VolForce_y(i,j,k)=VolForce_y(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ! uz ( receive from xm_dir, and send to xp_dir )
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= VolForce_z(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= VolForce_z(ixs1,j,k)
        enddo
      enddo
      call MPI_IRECV( Mat2HaloR,icount,real_type,ProcNgh(4),5,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat2HaloS,icount,real_type,ProcNgh(3),5,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(4)>=0) THEN
      i=y1start(1)
      do k=izs1,izs2
        do j=iys1,iys2
          VolForce_z(i,j,k)=VolForce_z(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF
    deallocate(Mat2HaloS,Mat2HaloR)

    !===================== Step2: receive from zp_dir, and send to zm_dir =====================!
    ixs1=y1start(1)-1;  ixs2=y1end(1)+2
    iys1=1;             iys2=nyp
    icount= (ixs2-ixs1+1)*(iys2-iys1+1)
    allocate(Mat2HaloS(ixs1:ixs2,iys1:iys2))
    allocate(Mat2HaloR(ixs1:ixs2,iys1:iys2))

    izs1=y1start(3)-1;   izs2=y1start(3)-1
    ! ux ( receive from zp_dir, and send to zm_dir )
    IF(ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= VolForce_x(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= VolForce_x(i,j,izs1)
        enddo
      enddo
      call MPI_IRECV( Mat2HaloR,icount,real_type,ProcNgh(1),6,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat2HaloS,icount,real_type,ProcNgh(2),6,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(1)>=0) THEN
      k= y1end(3)
      do j=iys1,iys2
        do i=ixs1,ixs2
          VolForce_x(i,j,k)= VolForce_x(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    ! uy ( receive from zp_dir, and send to zm_dir )
    IF(ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= VolForce_y(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= VolForce_y(i,j,izs1)
        enddo
      enddo
      call MPI_IRECV( Mat2HaloR,icount,real_type,ProcNgh(1),7,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat2HaloS,icount,real_type,ProcNgh(2),7,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(1)>=0) THEN
      k= y1end(3)
      do j=iys1,iys2
        do i=ixs1,ixs2
          VolForce_y(i,j,k)= VolForce_y(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    ! uz ( receive from zp_dir, and send to zm_dir )
    IF(ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= VolForce_z(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= VolForce_z(i,j,izs1)
        enddo
      enddo
      call MPI_IRECV( Mat2HaloR,icount,real_type,ProcNgh(1),8,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat2HaloS,icount,real_type,ProcNgh(2),8,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(1)>=0) THEN
      k= y1end(3)
      do j=iys1,iys2
        do i=ixs1,ixs2
          VolForce_z(i,j,k)= VolForce_z(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    izs1=y1end(3)+1;   izs2=y1end(3)+1
    ! ux ( receive from zm_dir, and send to zp_dir )
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= VolForce_x(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= VolForce_x(i,j,izs1)
        enddo
      enddo
      call MPI_IRECV( Mat2HaloR,icount,real_type,ProcNgh(2),9,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat2HaloS,icount,real_type,ProcNgh(1),9,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(2)>=0) THEN
      k= y1start(3)
      do j=iys1,iys2
        do i=ixs1,ixs2
          VolForce_x(i,j,k)= VolForce_x(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    ! uy ( receive from zm_dir, and send to zp_dir )
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= VolForce_y(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= VolForce_y(i,j,izs1)
        enddo
      enddo
      call MPI_IRECV( Mat2HaloR,icount,real_type,ProcNgh(2),10,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat2HaloS,icount,real_type,ProcNgh(1),10,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(2)>=0) THEN
      k= y1start(3)
      do j=iys1,iys2
        do i=ixs1,ixs2
          VolForce_y(i,j,k)= VolForce_y(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF
    deallocate(Mat2HaloS,Mat2HaloR)

    !===================== Step3: receive from xm_dir, and send to xp_dir =====================!
    ixs1=y1end(1)+1;    ixs2=y1end(1)+2
    iys1=1;             iys2=nyp    
    izs1=y1start(3)-1;  izs2=y1end(3)+2
    allocate(Mat3HaloS(ixs1:ixs2,            iys1:iys2,izs1:izs2))
    allocate(Mat3HaloR(y1start(1):y1start(1)+1,iys1:iys2,izs1:izs2))   

    ! ux ( receive from xm_dir, and send to xp_dir )
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      Mat3HaloR= VolForce_x(ixs1:ixs2,iys1:iys2, izs1:izs2)
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          do i=ixs1,ixs2
            Mat3HaloS(i,j,k)= VolForce_x(i,j,k)
          enddo
        enddo
      enddo
      icount= 2*(iys2-iys1+1)*(izs2-izs1+1)
      call MPI_IRECV( Mat3HaloR,icount,real_type,ProcNgh(4),11,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat3HaloS,icount,real_type,ProcNgh(3),11,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(4)>=0) THEN
      ixs1=y1start(1);     ixs2=y1start(1)+1
      do k=izs1,izs2
        do j=iys1,iys2
          do i=ixs1,ixs2
            VolForce_x(i,j,k)= VolForce_x(i,j,k)+ Mat3HaloR(i,j,k)
          enddo
        enddo
      enddo
    ENDIF
    deallocate(Mat3HaloS,Mat3HaloR)

    !===================== Step4: receive from zm_dir, and send to zp_dir =====================!
    ixs1=y1start(1)-1;  ixs2=y1end(1)+2
    iys1=1;             iys2=nyp
    izs1=y1end(3)+1;    izs2=y1end(3)+2
    allocate(Mat3HaloS(ixs1:ixs2,iys1:iys2,izs1:izs2))
    allocate(Mat3HaloR(ixs1:ixs2,iys1:iys2,y1start(3):y1start(3)+1))

    ! uz ( receive from zm_dir, and send to zp_dir )
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      Mat3HaloR= VolForce_z(ixs1:ixs2,iys1:iys2, izs1:izs2)
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          do i=ixs1,ixs2
            Mat3HaloS(i,j,k)= VolForce_z(i,j,k)
          enddo
        enddo
      enddo
      icount= (ixs2-ixs1+1)*(iys2-iys1+1)*2
      call MPI_IRECV( Mat3HaloR,icount,real_type,ProcNgh(2),12,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat3HaloS,icount,real_type,ProcNgh(1),12,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(2)>=0) THEN
      izs1=y1start(3);     izs2=y1start(3)+1
      do k=izs1,izs2
        do j=iys1,iys2
          do i=ixs1,ixs2
            VolForce_z(i,j,k)= VolForce_z(i,j,k)+ Mat3HaloR(i,j,k)
          enddo
        enddo
      enddo
    ENDIF
    deallocate(Mat3HaloS,Mat3HaloR)
  end subroutine Gather_Halo_IBMForce

end module ca_BC_and_Halo
