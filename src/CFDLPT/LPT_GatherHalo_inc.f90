  !******************************************************************
  ! Gather_Halo_dist_1
  !******************************************************************
  subroutine Gather_Halo_dist_1()
    implicit none

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

    !===================== Step1: receive from xp_dir, and send to xm_dir =====================!
    iys1=1;             iys2=nyp    
    izs1=y1start(3)-1;   izs2=y1end(3)+1
    icount= (iys2-iys1+1)*(izs2-izs1+1)
    allocate(Mat2HaloS(iys1:iys2,izs1:izs2))
    allocate(Mat2HaloR(iys1:iys2,izs1:izs2))

    ixs1=y1start(1)-1;   ixs2=y1start(1)-1
    ! ux ( receive from xp_dir, and send to xm_dir )
    ! DO nothing

    ! uy ( receive from xp_dir, and send to xm_dir )
    IF(ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= FpForce_y(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= FpForce_y(ixs1,j,k)
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
          FpForce_y(i,j,k)=FpForce_y(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ! uz ( receive from xp_dir, and send to xm_dir )
    IF(ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= FpForce_z(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= FpForce_z(ixs1,j,k)
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
          FpForce_z(i,j,k)=FpForce_z(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ixs1=y1end(1)+1;   ixs2=y1end(1)+1
    ! uy ( receive from xm_dir, and send to xp_dir )
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= FpForce_y(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= FpForce_y(ixs1,j,k)
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
          FpForce_y(i,j,k)=FpForce_y(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ! uz ( receive from xm_dir, and send to xp_dir )
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= FpForce_z(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= FpForce_z(ixs1,j,k)
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
          FpForce_z(i,j,k)=FpForce_z(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF
    deallocate(Mat2HaloS,Mat2HaloR)

    !===================== Step2: receive from zp_dir, and send to zm_dir =====================!
    ixs1=y1start(1)-1;   ixs2=y1end(1)+1
    iys1=1;             iys2=nyp
    icount= (ixs2-ixs1+1)*(iys2-iys1+1)
    allocate(Mat2HaloS(ixs1:ixs2,iys1:iys2))
    allocate(Mat2HaloR(ixs1:ixs2,iys1:iys2))

    izs1=y1start(3)-1;   izs2=y1start(3)-1
    ! ux ( receive from zp_dir, and send to zm_dir )
    IF(ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= FpForce_x(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= FpForce_x(i,j,izs1)
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
          FpForce_x(i,j,k)= FpForce_x(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    ! uy ( receive from zp_dir, and send to zm_dir )
    IF(ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= FpForce_y(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= FpForce_y(i,j,izs1)
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
          FpForce_y(i,j,k)= FpForce_y(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    ! uz ( receive from zp_dir, and send to zm_dir )
    ! Do nothing

    izs1=y1end(3)+1;   izs2=y1end(3)+1
    ! ux ( receive from zm_dir, and send to zp_dir )
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= FpForce_x(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= FpForce_x(i,j,izs1)
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
          FpForce_x(i,j,k)= FpForce_x(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    ! uy ( receive from zm_dir, and send to zp_dir )
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= FpForce_y(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= FpForce_y(i,j,izs1)
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
          FpForce_y(i,j,k)= FpForce_y(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF
    deallocate(Mat2HaloS,Mat2HaloR)

    !===================== Step3: receive from xm_dir, and send to xp_dir =====================!
    ixs1=y1end(1)+1;     ixs2=y1end(1)+1  ! Different to Gather_Halo_dist_2
    iys1=1;             iys2=nyp    
    izs1=y1start(3)-1;   izs2=y1end(3)+1  ! Different to Gather_Halo_dist_2
    allocate(Mat3HaloS(ixs1:ixs2,          iys1:iys2,izs1:izs2))
    allocate(Mat3HaloR(y1start(1):y1start(1),iys1:iys2,izs1:izs2)) ! Different to Gather_Halo_dist_2   

    ! ux ( receive from xm_dir, and send to xp_dir )
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      Mat3HaloR= FpForce_x(ixs1:ixs2,iys1:iys2, izs1:izs2)
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          do i=ixs1,ixs2
            Mat3HaloS(i,j,k)= FpForce_x(i,j,k)
          enddo
        enddo
      enddo
      icount= 1*(iys2-iys1+1)*(izs2-izs1+1) ! Different to Gather_Halo_dist_2
      call MPI_IRECV( Mat3HaloR,icount,real_type,ProcNgh(4),11,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat3HaloS,icount,real_type,ProcNgh(3),11,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(4)>=0) THEN
      ixs1=y1start(1);     ixs2=y1start(1) ! Different to Gather_Halo_dist_2
      do k=izs1,izs2
        do j=iys1,iys2
          do i=ixs1,ixs2
            FpForce_x(i,j,k)= FpForce_x(i,j,k)+ Mat3HaloR(i,j,k)
          enddo
        enddo
      enddo
    ENDIF
    deallocate(Mat3HaloS,Mat3HaloR)

    !===================== Step4: receive from zm_dir, and send to zp_dir =====================!
    ixs1=y1start(1)-1;   ixs2=y1end(1)+1 ! Different to Gather_Halo_dist_2
    iys1=1;             iys2=nyp
    izs1=y1end(3)+1;     izs2=y1end(3)+1 ! Different to Gather_Halo_dist_2
    allocate(Mat3HaloS(ixs1:ixs2,iys1:iys2,izs1:izs2))
    allocate(Mat3HaloR(ixs1:ixs2,iys1:iys2,y1start(3):y1start(3))) ! Different to Gather_Halo_dist_2

    ! uz ( receive from zm_dir, and send to zp_dir )
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      Mat3HaloR= FpForce_z(ixs1:ixs2,iys1:iys2, izs1:izs2)
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          do i=ixs1,ixs2
            Mat3HaloS(i,j,k)= FpForce_z(i,j,k)
          enddo
        enddo
      enddo
      icount= (ixs2-ixs1+1)*(iys2-iys1+1)*1 ! Different to Gather_Halo_dist_2
      call MPI_IRECV( Mat3HaloR,icount,real_type,ProcNgh(2),12,MPI_COMM_WORLD,requests(1),ierror)
      call MPI_ISSEND(Mat3HaloS,icount,real_type,ProcNgh(1),12,MPI_COMM_WORLD,requests(2),ierror)
      call MPI_WAITALL(2,requests,SRstatus,ierror)
    ENDIF
    IF(ProcNgh(2)>=0) THEN
      izs1=y1start(3);     izs2=y1start(3) ! Different to Gather_Halo_dist_2
      do k=izs1,izs2
        do j=iys1,iys2
          do i=ixs1,ixs2
            FpForce_z(i,j,k)= FpForce_z(i,j,k)+ Mat3HaloR(i,j,k)
          enddo
        enddo
      enddo
    ENDIF
    deallocate(Mat3HaloS,Mat3HaloR)
  end subroutine Gather_Halo_dist_1
  
  !******************************************************************
  ! Gather_Halo_dist_2
  !******************************************************************
  subroutine Gather_Halo_dist_2()
    implicit none

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

    !===================== Step1: receive from xp_dir, and send to xm_dir =====================!
    iys1=1;             iys2=nyp    
    izs1=y1start(3)-1;   izs2=y1end(3)+2
    icount= (iys2-iys1+1)*(izs2-izs1+1)
    allocate(Mat2HaloS(iys1:iys2,izs1:izs2))
    allocate(Mat2HaloR(iys1:iys2,izs1:izs2))

    ixs1=y1start(1)-1;   ixs2=y1start(1)-1
    ! ux ( receive from xp_dir, and send to xm_dir )
    IF(ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= FpForce_x(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= FpForce_x(ixs1,j,k)
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
          FpForce_x(i,j,k)=FpForce_x(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ! uy ( receive from xp_dir, and send to xm_dir )
    IF(ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= FpForce_y(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= FpForce_y(ixs1,j,k)
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
          FpForce_y(i,j,k)=FpForce_y(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ! uz ( receive from xp_dir, and send to xm_dir )
    IF(ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(3)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= FpForce_z(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= FpForce_z(ixs1,j,k)
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
          FpForce_z(i,j,k)=FpForce_z(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ixs1=y1end(1)+1;   ixs2=y1end(1)+1
    ! uy ( receive from xm_dir, and send to xp_dir )
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= FpForce_y(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= FpForce_y(ixs1,j,k)
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
          FpForce_y(i,j,k)=FpForce_y(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF

    ! uz ( receive from xm_dir, and send to xp_dir )
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloR(j,k)= FpForce_z(ixs1,j,k)
        enddo
      enddo
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          Mat2HaloS(j,k)= FpForce_z(ixs1,j,k)
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
          FpForce_z(i,j,k)=FpForce_z(i,j,k)+ Mat2HaloR(j,k)
        enddo
      enddo
    ENDIF
    deallocate(Mat2HaloS,Mat2HaloR)

    !===================== Step2: receive from zp_dir, and send to zm_dir =====================!
    ixs1=y1start(1)-1;   ixs2=y1end(1)+2
    iys1=1;             iys2=nyp
    icount= (ixs2-ixs1+1)*(iys2-iys1+1)
    allocate(Mat2HaloS(ixs1:ixs2,iys1:iys2))
    allocate(Mat2HaloR(ixs1:ixs2,iys1:iys2))

    izs1=y1start(3)-1;   izs2=y1start(3)-1
    ! ux ( receive from zp_dir, and send to zm_dir )
    IF(ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= FpForce_x(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= FpForce_x(i,j,izs1)
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
          FpForce_x(i,j,k)= FpForce_x(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    ! uy ( receive from zp_dir, and send to zm_dir )
    IF(ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= FpForce_y(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= FpForce_y(i,j,izs1)
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
          FpForce_y(i,j,k)= FpForce_y(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    ! uz ( receive from zp_dir, and send to zm_dir )
    IF(ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= FpForce_z(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= FpForce_z(i,j,izs1)
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
          FpForce_z(i,j,k)= FpForce_z(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    izs1=y1end(3)+1;   izs2=y1end(3)+1
    ! ux ( receive from zm_dir, and send to zp_dir )
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= FpForce_x(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= FpForce_x(i,j,izs1)
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
          FpForce_x(i,j,k)= FpForce_x(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF

    ! uy ( receive from zm_dir, and send to zp_dir )
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(1)==nrank also !!!
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloR(i,j)= FpForce_y(i,j,izs1)
        enddo
      enddo
    ELSE
      do j=iys1,iys2
        do i=ixs1,ixs2
          Mat2HaloS(i,j)= FpForce_y(i,j,izs1)
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
          FpForce_y(i,j,k)= FpForce_y(i,j,k)+ Mat2HaloR(i,j)
        enddo
      enddo
    ENDIF
    deallocate(Mat2HaloS,Mat2HaloR)

    !===================== Step3: receive from xm_dir, and send to xp_dir =====================!
    ixs1=y1end(1)+1;     ixs2=y1end(1)+2
    iys1=1;             iys2=nyp    
    izs1=y1start(3)-1;   izs2=y1end(3)+2
    allocate(Mat3HaloS(ixs1:ixs2,            iys1:iys2,izs1:izs2))
    allocate(Mat3HaloR(y1start(1):y1start(1)+1,iys1:iys2,izs1:izs2))   

    ! ux ( receive from xm_dir, and send to xp_dir )
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      Mat3HaloR= FpForce_x(ixs1:ixs2,iys1:iys2, izs1:izs2)
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          do i=ixs1,ixs2
            Mat3HaloS(i,j,k)= FpForce_x(i,j,k)
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
            FpForce_x(i,j,k)= FpForce_x(i,j,k)+ Mat3HaloR(i,j,k)
          enddo
        enddo
      enddo
    ENDIF
    deallocate(Mat3HaloS,Mat3HaloR)

    !===================== Step4: receive from zm_dir, and send to zp_dir =====================!
    ixs1=y1start(1)-1;   ixs2=y1end(1)+2
    iys1=1;             iys2=nyp
    izs1=y1end(3)+1;     izs2=y1end(3)+2
    allocate(Mat3HaloS(ixs1:ixs2,iys1:iys2,izs1:izs2))
    allocate(Mat3HaloR(ixs1:ixs2,iys1:iys2,y1start(3):y1start(3)+1))

    ! uz ( receive from zm_dir, and send to zp_dir )
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself, and ProcNgh(4)==nrank also !!!
      Mat3HaloR= FpForce_z(ixs1:ixs2,iys1:iys2, izs1:izs2)
    ELSE
      do k=izs1,izs2
        do j=iys1,iys2
          do i=ixs1,ixs2
            Mat3HaloS(i,j,k)= FpForce_z(i,j,k)
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
            FpForce_z(i,j,k)= FpForce_z(i,j,k)+ Mat3HaloR(i,j,k)
          enddo
        enddo
      enddo
    ENDIF
    deallocate(Mat3HaloS,Mat3HaloR)
  end subroutine Gather_Halo_dist_2
