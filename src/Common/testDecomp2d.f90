program testDecomp2d
  use MPI
  use m_Decomp2d
  implicit none
  integer::ierror
  integer,parameter::RK=8
  integer,parameter::nxc=128,nyc=95,nzc=123
  integer,parameter,dimension(6):: BcOption=[0,0,0,0,0,0]
  
  ! locals
  real(RK)::rTemp
  type(MatBound)::mbInfo
  type(HaloInfo)::hiInfo
  integer::p_row=4, p_col=2
  character(len=10)::RowColStr
  integer::ic,jc,kc,i,j,k,intT,iTV(8)
  real(RK),allocatable,dimension(:,:,:)::Arr1,Arr2,Arr3

  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
  
  intT=command_argument_count()
  if((intT/=0 .and. intT/=2) .and. nrank==0) then
    write(*,*)'command argument wrong!'; stop
  endif  
  if(intT==2) then
    call get_command_argument(1,RowColStr)
    read(RowColStr,*) p_row
    call get_command_argument(2,RowColStr)
    read(RowColStr,*) p_col
  endif
  if(p_row*p_col /= nproc) then
    write(*,*)'p_row*p_col /= nproc, wrong!'; stop
  endif  
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  call decomp_2d_init(nxc,nyc,nzc,nproc,p_row,p_col,2,BcOption)
 
  mbInfo%pencil = 2  
  mbInfo%xme=3;  mbInfo%xpe=2
  mbInfo%yme=1;  mbInfo%ype=3
  mbInfo%zme=2;  mbInfo%zpe=3

  hiInfo%pencil = 2
  hiInfo%xmh=1;  hiInfo%xph=1
  hiInfo%ymh=1;  hiInfo%yph=1
  hiInfo%zmh=1;  hiInfo%zph=1
  call myallocate(Arr1, mbInfo, opt_global=.true.)
  deallocate(Arr1)
  allocate(Arr1(mbInfo%xmm:mbInfo%xpm,mbInfo%ymm:mbInfo%ypm,mbInfo%zmm:mbInfo%zpm))
  do kc=y1start(3),y1end(3)
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        Arr1(ic,jc,kc)=real(ic,RK)+real(jc,RK)+real(kc,RK)
      enddo
    enddo
  enddo
      
  ! updata_halo
  call update_halo(Arr1, mbInfo, hiInfo)
  do k=y1start(3)-hiInfo%zmh,y1end(3)+hiInfo%zph
    kc=k
    if(kc<1)   kc=kc+nzc
    if(kc>nzc) kc=kc-nzc
    do j=y1start(2)-hiInfo%ymh,y1end(2)+hiInfo%yph
      jc=j
      if(jc<1)   jc=jc+nyc
      if(jc>nyc) jc=jc-nyc
      do i=y1start(1)-hiInfo%xmh,y1end(1)+hiInfo%xph
        ic=i
        if(ic<1)   ic=ic+nxc
        if(ic>nxc) ic=ic-nxc
        rTemp=real(ic,RK)+real(jc,RK)+real(kc,RK)
        if(Arr1(i,j,k)/=rTemp) then
          print*,'update_halo Wrong========',i,j,k
          print*,Arr1(i,j,k),rTemp
          stop
        endif
      enddo
    enddo
  enddo
  print*,nrank,'Test updata_halo OK !'
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

  ! sum_halo
  call date_and_time(values=iTV); !iTV=0
  call random_seed(size= ic)
  call random_seed(put = nrank*iTV(8)+[(jc,jc=1,ic)])
  call random_number(Arr1)
  allocate(Arr2(nxc,nyc,nzc)); Arr2=0.0_RK
  allocate(Arr3(nxc,nyc,nzc)); Arr3=0.0_RK
  do k=y1start(3)-hiInfo%zmh,y1end(3)+hiInfo%zph
    kc=k
    if(kc<1)   kc=kc+nzc
    if(kc>nzc) kc=kc-nzc
    do j=y1start(2)-hiInfo%ymh,y1end(2)+hiInfo%yph
      jc=j
      if(jc<1)   jc=jc+nyc
      if(jc>nyc) jc=jc-nyc
      do i=y1start(1)-hiInfo%xmh,y1end(1)+hiInfo%xph
        ic=i
        if(ic<1)   ic=ic+nxc
        if(ic>nxc) ic=ic-nxc
        Arr2(ic,jc,kc)=Arr2(ic,jc,kc)+Arr1(i,j,k)
      enddo
    enddo
  enddo
  call MPI_ALLREDUCE(Arr2,Arr3,nxc*nyc*nzc,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierror)
  
  call sum_halo(Arr1, mbInfo, hiInfo)
  do kc=y1start(3),y1end(3)
    do jc=y1start(2),y1end(2)
      do ic=y1start(1),y1end(1)
        rTemp=abs(Arr3(ic,jc,kc)-Arr1(ic,jc,kc))
        if(rTemp>1.0E-14) then
          print*,'sum_halo Wrong ========',ic,jc,kc
          print*,Arr1(ic,jc,kc),Arr3(ic,jc,kc)
          stop
        endif
      enddo
    enddo
  enddo
  print*,nrank,'Test sum_halo OK !'
  
  call MPI_FINALIZE(ierror)
end program testDecomp2d
