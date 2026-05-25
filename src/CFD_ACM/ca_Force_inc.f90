#include "definitions_inc.f90"
! This part is used to test the work of Chan-Braun et al., J.Fluid Mech. (2011), vol.684,pp.441-474.
if(mod(itime,ivstats)/=0) cycle
Block
  integer::pid
  real(RK)::pmAlphaC
  real(RK),dimension(3)::Force,Torqu
  real(RK),dimension(8,3)::ForceAndTorqueT1
  
  if(ns==1 .and. GPrtcl_list%nlocal>0) then
    if(allocated(Ave_FpForce))  Ave_FpForce(1:max(GPrtcl_list%nlocal,1))=zero_r3
    if(allocated(Ave_FpTorque)) Ave_FpTorque(1:max(GPrtcl_list%nlocal,1))=zero_r3
  endif
  pmAlphaC= pmAlpha/dt
  do pid=1,GPrtcl_list%nlocal
    Ave_FpForce(pid) =Ave_FpForce(pid) +pmAlphaC*GPrtcl_FpForce(pid)
    Ave_FpTorque(pid)=Ave_FpTorque(pid)+pmAlphaC*GPrtcl_FpTorque(pid)
  enddo
  if(ns/=iadvance) cycle

  nfstime=nfstime+1
  ForceAndTorqueT1=0.0_RK
  do pid=1,GPrtcl_list%nlocal
    Force =[Ave_FpForce(pid)%x, Ave_FpForce(pid)%y, Ave_FpForce(pid)%z ]
    Torqu =[Ave_FpTorque(pid)%x,Ave_FpTorque(pid)%y,Ave_FpTorque(pid)%z]
    ForceAndTorqueT1(1,:)=ForceAndTorqueT1(1,:)+Force
    ForceAndTorqueT1(2,:)=ForceAndTorqueT1(2,:)+Force*Force
    ForceAndTorqueT1(3,:)=ForceAndTorqueT1(3,:)+Force*Force*Force
    ForceAndTorqueT1(4,:)=ForceAndTorqueT1(4,:)+Force*Force*Force*Force
    ForceAndTorqueT1(5,:)=ForceAndTorqueT1(5,:)+Torqu
    ForceAndTorqueT1(6,:)=ForceAndTorqueT1(6,:)+Torqu*Torqu
    ForceAndTorqueT1(7,:)=ForceAndTorqueT1(7,:)+Torqu*Torqu*Torqu
    ForceAndTorqueT1(8,:)=ForceAndTorqueT1(8,:)+Torqu*Torqu*Torqu*Torqu
  enddo
  ForceAndTorque=ForceAndTorque+ForceAndTorqueT1/real(DEM_Opt%numPrtcl,RK)
  if(mod(itime,SaveStat)/=0) cycle
  
  filename = ' '
  call MPI_REDUCE(ForceAndTorque,ForceAndTorqueT1,24,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierr)
  if(nrank==0) then
    filename = strip(DEM_opt%ResultsDir) // "FpForce" // int2str(ilast,10) // '.txt'
    open(newunit=iUnit,file=filename,status='old',form='formatted',position='append',IOSTAT=ierr)
    if(ierr /= 0) then
      call MainLog%CheckForError(ErrT_Pass,"ChannelACM_Iterate: ","Cannot open file: "//filename)
    else
      write(iUnit,'(I7,24ES24.15)') itime,ForceAndTorqueT1/real(nfstime,RK)
      close(iUnit,IOSTAT=ierr)
    endif
  endif
  
  nfstime=0; ForceAndTorque=0.0_RK
EndBlock
