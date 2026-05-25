! Zheng Gong, 2024-05-26
! gfortran -cpp -Wall -ffree-line-length-none -fimplicit-none -O3 AppendRestart.f90 -o AppendRestart
program main
  implicit none
  integer,parameter::RK=8
  integer,parameter::numPrtcl_old=11806
  integer,parameter::numPrtcl_add=70000
  real(RK),parameter::Diameter_new=0.314980262473718E-3_RK*1.01_RK
  character(80),parameter::File_old='RestartForMoveZY_0000000000_'
  character(80),parameter::File_new='RestartForMoveZY_0000000000'
  integer,parameter::exeType=2
  integer,parameter::iTypeAdd=2
  integer,parameter::iUsrMarkAdd=1
  integer,parameter::nTry=100
  integer,parameter::n_Integer=3
  integer,parameter::GPrtcl_list_tsize=2 ! Note here. If PI_Method=AB2, GPrtcl_list_tsize=2; If PI_Method=AB3, GPrtcl_list_tsize=3;
  integer,parameter::GPrtcl_list_rsize=2
  integer,parameter::int_byte=sizeof(1_4)
  integer,parameter::real_byte=sizeof(1.0_8)
  real(RK),parameter::MkPrtclMinpoint(3) = [0.000_RK,  0.010_RK,   0.000_RK]+Diameter_new
  real(RK),parameter::MkPrtclMaxpoint(3) = [0.0680357366943232_RK, 0.0160639933861596_RK, 0.0340178683471616_RK]-Diameter_new
  
  ! locals
  logical:: ExitFlag
  integer(8):: TotalByte_Old, disp, nRemain
  real(RK),dimension(:),allocatable::GPrtcl_nRemain
  integer,dimension(:),allocatable:: GPrtcl_ncv_old, GPrtcl_ncv_add
  integer,dimension(:,:),allocatable::GPrtcl_integer_old, GPrtcl_integer_add
  real(RK):: MinDist2, MkPrtclDomain(3), randt(3), PosTmp(3), DistVec(3), Dist2
  integer:: nReal3,numPrtcl_new, ierror, fido, fidn, np_nCntctTotal(2),k, kt, iTry
  real(RK),dimension(:,:),allocatable::GPrtcl_Pos_old, GPrtcl_Pos_Add, GPrtcl_nReal3_old, GPrtcl_nReal3_add
  
  MinDist2=Diameter_new*Diameter_new
  nReal3=2*(1 +GPrtcl_list_tsize +GPrtcl_list_rsize)
  if(exeType==2) then
    nReal3=nReal3+2
  endif
  nReal3=nReal3-1 ! '-1' corresponds to GPrtcl_Pos
  
  numPrtcl_new = numPrtcl_old+numPrtcl_add
  MkPrtclDomain= MkPrtclMaxpoint-MkPrtclMinpoint
  open(newunit=fido,file=File_old,status='old',form='unformatted',access='stream',action='read',position='append',IOSTAT=ierror)
  inquire(unit=fido, Pos=TotalByte_Old)
  rewind(unit=fido)
  open(newunit=fidn,file=File_new,status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierror)
  
  ! Read numPrtcl_old and nCntctTotal
  read(unit=fido) np_nCntctTotal
  if(np_nCntctTotal(1) /= numPrtcl_old) then
    print*, 'numPrtcl_old Wrong !!!'; stop
  endif
  
  ! Write numPrtcl_new and nCntctTotal
  np_nCntctTotal(1)=numPrtcl_new;
  write(unit=fidn) np_nCntctTotal

  ! GPrtcl_Pos
  allocate( GPrtcl_Pos_old(3,numPrtcl_old) )
  read(unit=fido) GPrtcl_Pos_old
  write(unit=fidn) GPrtcl_Pos_old
  deallocate(GPrtcl_Pos_old)
  allocate( GPrtcl_Pos_add(3,numPrtcl_add) )
  do k=1,numPrtcl_add
    do iTry=1,nTry
      call random_number(randt)
      PosTmp=randt*MkPrtclDomain + MkPrtclMinpoint
      ExitFlag= .true.
      do kt=1, k-1
        DistVec = PosTmp-GPrtcl_Pos_Add(:,kt)
        Dist2 = DistVec(1)*DistVec(1) +DistVec(2)*DistVec(2) +DistVec(3)*DistVec(3)
        if(Dist2<MinDist2) then
          ExitFlag=.false.; exit
        endif
      enddo
      if(ExitFlag) exit
    enddo
    if(.not. ExitFlag) then
      print*, 'ExitFlag wrong, Try again'; stop
    endif
    if(mod(k,1000)==0) write(*, '(A,I8,A,I5)') 'k=',k,' iTry=',iTry
    GPrtcl_Pos_Add(:,k)=PosTmp
  enddo
  write(unit=fidn) GPrtcl_Pos_add
  deallocate(GPrtcl_Pos_add)
  
  ! GPrtcl_integer
  allocate( GPrtcl_integer_old(n_Integer,numPrtcl_old) )
  read(unit=fido) GPrtcl_integer_old
  write(unit=fidn) GPrtcl_integer_old 
  deallocate( GPrtcl_integer_old )
  allocate( GPrtcl_integer_add(n_Integer,numPrtcl_add) )
  do k=1,numPrtcl_add
    kt=k+numPrtcl_old
    GPrtcl_integer_add(1,k)=kt            ! GPrtcl_id
    GPrtcl_integer_add(2,k)=iTypeAdd      ! GPrtcl_pType
    GPrtcl_integer_add(3,k)=iUsrMarkAdd   ! GPrtcl_UsrMark
  enddo
  write(unit=fidn) GPrtcl_integer_add
  deallocate( GPrtcl_integer_add )
  
  ! GPrtcl_nReal3
  allocate( GPrtcl_nReal3_old(3*nReal3,numPrtcl_old) )
  read(unit=fido) GPrtcl_nReal3_old
  write(unit=fidn) GPrtcl_nReal3_old 
  deallocate( GPrtcl_nReal3_old )
  allocate( GPrtcl_nReal3_add(3*nReal3,numPrtcl_add) )
  GPrtcl_nReal3_add = 0.0_RK
  write(unit=fidn) GPrtcl_nReal3_add 
  deallocate( GPrtcl_nReal3_add )
  
  ! GPrtcl_ncv
  allocate( GPrtcl_ncv_old(numPrtcl_old) )
  read(unit=fido) GPrtcl_ncv_old
  write(unit=fidn) GPrtcl_ncv_old
  deallocate(GPrtcl_ncv_old)
  allocate( GPrtcl_ncv_add(numPrtcl_add) )
  GPrtcl_ncv_add = 0
  write(unit=fidn) GPrtcl_ncv_add
  deallocate(GPrtcl_ncv_add)
  
  ! Remainings
  inquire(unit=fido, Pos=disp)
  nRemain= TotalByte_Old - disp
  if(mod(nRemain,real_byte) /=0) then
    print*, 'nRemain Wrong'; stop
  endif 
  if(nRemain>0) then
    nRemain=nRemain/real_byte
    allocate( GPrtcl_nRemain(nRemain) )
    read(unit=fido) GPrtcl_nRemain
    write(unit=fidn) GPrtcl_nRemain
    deallocate(GPrtcl_nRemain)
  endif
  
  close(unit=fido, IOSTAT=ierror)
  close(unit=fidn, IOSTAT=ierror)
end program main
