module m_FlowType_Channel
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use m_MeshAndMetries
  use m_Variables,only:mb1
#ifdef CFDACM
  use m_Variables,only:FluidIndicator
#endif
#ifdef CFDLPT_TwoWay
  use m_Variables,only:FpForce_x,FpForce_y,FpForce_z
#endif
  implicit none
  private

  ! statistics variabls
  logical::IsUxConst
  integer:: nfstime
  real(RK):: PrGradsum
  real(RK),allocatable,dimension(:,:):: SumStat 

  public:: InitVelocity_CH,Update_uy_ym_CH,InitStatVar_CH,clcStat_CH
contains

#ifdef CFDLPT_TwoWay
#define NCHASTAT 42
#elif CFDACM
#define NCHASTAT 19
#else
#define NCHASTAT 35
#endif
  !******************************************************************
  ! InitVelocity_CH
  !******************************************************************
  subroutine InitVelocity_CH(ux,uy,uz,Deviation,chFile)
    implicit none
    character(*),intent(in)::chFile
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2), y1start(3):y1end(3)),intent(inout)::Deviation
  
    ! locals
    logical::IsUbulkGlobal
    real(RK)::ubulk,ubulkTemp,ybulk1,ybulk2
    integer::ii,ierror,ic,jc,kc,m1,m2,nUnit
    real(RK)::retau_guass,utau_guass,height,rem,wx,wz,xlxPlus,zlzPlus
    real(RK)::xplus,yplus,zplus,yct,ybar,xpt,zpt,uxmean,ratiot,uzmean(nyc),uzmeanR(nyc)
    NAMELIST/ubulk_Param/IsUxConst,IsUbulkGlobal,ubulk,ybulk1,ybulk2

    open(newunit=nUnit, file=chFile, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitVelocity_CH","Cannot open file: "//trim(chFile))
    endif
    read(nUnit,nml=ubulk_Param)
    close(nUnit,IOSTAT=ierror)
    height=ybulk2-ybulk1
    ubulkTemp=ubulk
    if(IsUbulkGlobal) then
      ubulkTemp=ubulk*yly/height
    endif
    if(FlowType==FT_CH) then
      height=half*height
    endif
    rem=ubulkTemp*height/xnu
    ux=zero; uy=zero; uz=zero
    if(abs(ubulk)<1.0D-12)return
    
    retau_guass = 0.1538_RK*rem**0.887741_RK
    utau_guass  = retau_guass*xnu/height
    if(nrank==0)print*,'************** retau_gauss=',retau_guass
    if(nrank==0)print*,'************** utau_gauss= ',utau_guass

    call system_clock(count=ierror); !ierror=0
    call random_seed(size = ii)
    call random_seed(put = ierror+63946*(/(ic-1,ic=1,ii)/))
    call random_number(Deviation)
    Deviation= 0.2_RK* Deviation + 0.9_RK ! [0.8, 1.2]

    !modulation of the random noise + initial velocity profile
    uxmean=zero; uzmean=zero
    wx=twopi/500.0_RK; wz=twopi/200.0_RK
    xlxPlus=xlx*utau_guass/xnu;   zlzPlus=zlz*utau_guass/xnu;
    m1=floor(xlxPlus*wx/twopi)+1; wx=real(m1,RK)*twopi/xlxPlus
    m2=floor(zlzPlus*wz/twopi)+1; wz=real(m2,RK)*twopi/zlzPlus
    do jc=y1start(2),y1end(2)
      yct = height-abs(height-(yc(jc)-ybulk1))
      if(yc(jc)<ybulk1 .or. yc(jc)>=ybulk2)yct=zero
      ybar= yct/height; yplus=utau_guass*yct/xnu
      do kc=y1start(3),y1end(3)
        zpt  =real(kc-1,kind=RK)*dz+dz*half
        zplus=utau_guass*zpt/xnu
        do ic=y1start(1),y1end(1)
          xpt  =real(ic-1,kind=RK)*dx+dx*half
          xplus=utau_guass*xpt/xnu
          !ux(ic,jc,kc) = 0.0052_RK*ubulkTemp*yplus*exp(-yplus*yplus/1800.0_RK)*cos(wz*zplus)*Deviation(ic,jc,kc) ! original expression
          !uz(ic,jc,kc) = 0.0050_RK*ubulkTemp*yplus*exp(-yplus*yplus/1800.0_RK)*sin(wx*xplus)*Deviation(ic,jc,kc) ! original expression
          ux(ic,jc,kc) = ubulkTemp*ybar*exp(-4.5_RK*ybar*ybar)*cos(wz*zplus)*Deviation(ic,jc,kc)
          uz(ic,jc,kc) = ubulkTemp*ybar*exp(-4.5_RK*ybar*ybar)*sin(wx*xplus)*Deviation(ic,jc,kc)
          ux(ic,jc,kc) = ux(ic,jc,kc)+ three*ubulkTemp*(ybar-half*ybar*ybar)
          uxmean     = uxmean    + ux(ic,jc,kc)*dyp(jc)
          uzmean(jc) = uzmean(jc)+ uz(ic,jc,kc)
        enddo
      enddo
      uzmean(jc)=uzmean(jc)/real(nxc*nzc,RK)
    enddo
    call MPI_ALLREDUCE(uxmean,ratiot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    ratiot=ratiot/real(nxc*nzc,RK)/(ybulk2-ybulk1)/ubulkTemp
    ux=ux*ratiot
    call MPI_ALLREDUCE(uzmean,uzmeanR,nyc,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
    do jc=y1start(2),y1end(2)
      uz(:,jc,:)=uz(:,jc,:)-uzmeanR(jc)
    enddo
    Deviation=zero
  end subroutine InitVelocity_CH

  !******************************************************************
  ! Update_uy_ym_CH
  !******************************************************************   
  subroutine Update_uy_ym_CH(uy_ym, duy_ym, TimeNew)
    implicit none
    real(RK),dimension(y1start(1):y1end(1),y1start(3):y1end(3)),intent(inout):: duy_ym,uy_ym    
    real(RK),intent(in):: TimeNew
  
    duy_ym = uy_ym
     
    ! update uy_ym here
    uy_ym = uyBcValue(ym_dir)
    duy_ym = uy_ym - duy_ym
  end subroutine Update_uy_ym_CH

  !******************************************************************
  ! InitStatVar_CH
  !******************************************************************
  subroutine InitStatVar_CH(chFile)
    implicit none
    character(*),intent(in)::chFile

    ! locals
    integer::ierror,nUnit
    logical::IsUbulkGlobal
    character(len=128)::filename
    real(RK)::ubulk,ybulk1,ybulk2
    NAMELIST/ubulk_Param/IsUxConst,IsUbulkGlobal,ubulk,ybulk1,ybulk2

    open(newunit=nUnit, file=chFile, status='old',form='formatted',IOSTAT=ierror )
    if(ierror/=0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","Cannot open file: "//trim(chFile))
    endif
    read(nUnit,nml=ubulk_Param)
    close(nUnit,IOSTAT=ierror)
        
    if(nrank==0) then
      if(mod(saveStat,ivstats)/=0 )    call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","ivstats wrong !!!")
      if(IsUxConst)then
        write(filename,'(A,I10.10)')trim(ResultsDir)//"PrGrad",ilast
        open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
        if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","Cannot open file: "//trim(filename))
        close(nUnit,IOSTAT=ierror)
      endif
    endif
#ifndef CFDACM
    allocate(SumStat(NCHASTAT,nyp),Stat=ierror)
#else
    allocate(SumStat(NCHASTAT,nyc),Stat=ierror)
#endif
    if(ierror /= 0) call MainLog%CheckForError(ErrT_Abort,"InitStatVar_CH","Allocation failed")  
    nfstime=0; SumStat=zero; PrGradsum=zero
  end subroutine InitStatVar_CH

#ifndef CFDACM
  !******************************************************************
  ! clcStat_CH
  !******************************************************************
  subroutine clcStat_CH(ux,uy,uz,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
   
    ! locals
    character(len=128)::filename
    integer::ic,jc,kc,im,jm,km,ip,jp,kp,ierror,nUnit
    real(RK)::uxloc,uyloc,uzloc,prloc,uxCell,uyCell,uzCell,prloc2,inxz,infstime,cac,caj,cacU,dudyU,dvdyM
    real(RK)::dudxx,dudyy,dudzz,dvdxx,dvdyy,dvdzz,dwdxx,dwdyy,dwdzz,SumStatR(NCHASTAT,nyp),SumVec(NCHASTAT),rdxh
    real(RK)::dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,vor_x,vor_y,vor_z,InterpY1,InterpY2,dudzC,dvdzC,dudxC,dvdxU
#ifdef CFDLPT_TwoWay
    real(RK)::uyCellm
#endif
    rdxh=rdx*half
    inxz = one/(real(nxc,RK)*real(nzc,RK))
    DO jc=1,nyc
      jp=jc+1; jm=jc-1;
      InterpY1= half*YinterpCoe(jc); InterpY2=half-InterpY1
      cac=rdyc(jc);cacU=rdyc(jp); caj=rdyp(jc); SumVec=zero
      do kc=y1start(3),y1end(3)
        km=kc-1;kp=kc+1
        do ic=y1start(1),y1end(1)
          im=ic-1;ip=ic+1

          uxloc = ux(ic,jc,kc)
          uyloc = uy(ic,jc,kc)
          uzloc = uz(ic,jc,kc)
          prloc = pressure(ic,jc,kc)
          uxCell= half*(ux(ic,jc,kc)+ux(ip,jc,kc))
          uyCell= half*(uy(ic,jc,kc)+uy(ic,jp,kc))
          uzCell= half*(uz(ic,jc,kc)+uz(ic,jc,kp))
          prloc2= InterpY1*(pressure(im,jm,kc)+pressure(ic,jm,kc))+ InterpY2*(pressure(im,jc,kc)+pressure(ic,jc,kc))
 
          dudx= (ux(ip,jc,kc)-ux(ic,jc,kc))*rdx
          dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
          dudz= (ux(ic,jc,kc)-ux(ic,jc,km))*rdz
          dvdx= (uy(ic,jc,kc)-uy(im,jc,kc))*rdx
          dvdy= (uy(ic,jp,kc)-uy(ic,jc,kc))*caj
          dvdz= (uy(ic,jc,kc)-uy(ic,jc,km))*rdz
          dwdx= (uz(ic,jc,kc)-uz(im,jc,kc))*rdx
          dwdy= (uz(ic,jc,kc)-uz(ic,jm,kc))*cac
          dwdz= (uz(ic,jc,kp)-uz(ic,jc,kc))*rdz

          dudxx= (ux(ip,jc,kc)-two*ux(ic,jc,kc)+ux(im,jc,kc))*rdx2
          dudyy= ap2c(jc)*ux(ic,jp,kc)+ac2c(jc)*ux(ic,jc,kc)+am2c(jc)*ux(ic,jm,kc)
          dudzz= (ux(ic,jc,kp)-two*ux(ic,jc,kc)+ux(ic,jc,km))*rdz2
          dvdxx= (uy(ip,jc,kc)-two*uy(ic,jc,kc)+uy(im,jc,kc))*rdx2
          dvdyy= ap2p(jc)*uy(ic,jp,kc)+ac2p(jc)*uy(ic,jc,kc)+am2p(jc)*uy(ic,jm,kc)
          dvdzz= (uy(ic,jc,kp)-two*uy(ic,jc,kc)+uy(ic,jc,km))*rdz2
          dwdxx= (uz(ip,jc,kc)-two*uz(ic,jc,kc)+uz(im,jc,kc))*rdx2
          dwdyy= ap2c(jc)*uz(ic,jp,kc)+ac2c(jc)*uz(ic,jc,kc)+am2c(jc)*uz(ic,jm,kc)
          dwdzz= (uz(ic,jc,kp)-two*uz(ic,jc,kc)+uz(ic,jc,km))*rdz2
          vor_x= dwdy-dvdz
          vor_y= dudz-dwdx
          vor_z= dvdx-dudy

          dudxC= (ux(ip,jc,kc)-ux(im,jc,kc))*rdxh
          dvdxU= (uy(ic,jp,kc)-uy(im,jp,kc))*rdx
          dudyU= (ux(ic,jp,kc)-ux(ic,jc,kc))*cacU
          dvdyM= (uy(im,jp,kc)-uy(im,jc,kc))*caj
          dudzC= (InterpY1*(ux(ic,jm,kp)-ux(ic,jm,km)) +InterpY2*(ux(ic,jc,kp)-ux(ic,jc,km)))*rdz
          dvdzC= (uy(im,jc,kp)-uy(im,jc,km)+ uy(ic,jc,kp)-uy(ic,jc,km))*rdz*quarter

          SumVec( 1)=SumVec( 1)+ uxloc
          SumVec( 2)=SumVec( 2)+ uyloc                     ! yp
          SumVec( 3)=SumVec( 3)+ uzloc
          SumVec( 4)=SumVec( 4)+ prloc
          SumVec( 5)=SumVec( 5)+ uxloc*uxloc
          SumVec( 6)=SumVec( 6)+ uyloc*uyloc               ! yp  
          SumVec( 7)=SumVec( 7)+ uzloc*uzloc
          SumVec( 8)=SumVec( 8)+ prloc*prloc
          SumVec( 9)=SumVec( 9)+ uxCell*uyCell
          SumVec(10)=SumVec(10)+ uyCell*uzCell
          SumVec(11)=SumVec(11)+ uxCell*uzCell
          SumVec(12)=SumVec(12)+ uxCell*prloc
          SumVec(13)=SumVec(13)+ uyCell*prloc
          SumVec(14)=SumVec(14)+ uzCell*prloc
          SumVec(15)=SumVec(15)+ uxCell*uxCell*uyCell 
          SumVec(16)=SumVec(16)+ uyloc *uyloc *uyloc       ! yp
          SumVec(17)=SumVec(17)+ uzCell*uzCell*uyCell
          SumVec(18)=SumVec(18)+ uxCell*uyCell*uyCell
          SumVec(19)=SumVec(19)+ uxloc*uxloc*uxloc
          SumVec(20)=SumVec(20)+ uzloc*uzloc*uzloc
          SumVec(21)=SumVec(21)+ uxloc*uxloc*uxloc*uxloc
          SumVec(22)=SumVec(22)+ uyloc*uyloc*uyloc*uyloc   ! yp
          SumVec(23)=SumVec(23)+ uzloc*uzloc*uzloc*uzloc
          SumVec(24)=SumVec(24)+ prloc*dudx
          SumVec(25)=SumVec(25)+ prloc*dvdy
          SumVec(26)=SumVec(26)+ prloc*dwdz
          SumVec(27)=SumVec(27)+ prloc2*(dudy+dvdx)        ! yp !
          SumVec(28)=SumVec(28)+ uxloc*(dudxx+dudyy+dudzz)
          SumVec(29)=SumVec(29)+ uyloc*(dvdxx+dvdyy+dvdzz) ! yp
          SumVec(30)=SumVec(30)+ uzloc*(dwdxx+dwdyy+dwdzz)
          SumVec(31)=SumVec(31)+ dudxC*(dvdxU+dvdx)*half+ (dudyU+dudy)*(dvdyM+dvdy)*quarter
          SumVec(32)=SumVec(32)+ dudzC*dvdzC          ! yp
          SumVec(33)=SumVec(33)+ vor_x*vor_x          ! yp !
          SumVec(34)=SumVec(34)+ vor_y*vor_y
          SumVec(35)=SumVec(35)+ vor_z*vor_z          ! yp !
#ifdef CFDLPT_TwoWay
          SumVec(36)=SumVec(36)+ FpForce_x(ic,jc,kc)
          SumVec(37)=SumVec(37)+ FpForce_y(ic,jc,kc)  ! yp !
          SumVec(38)=SumVec(38)+ FpForce_z(ic,jc,kc)
          SumVec(39)=SumVec(39)+ FpForce_x(ic,jc,kc)*uxloc
          SumVec(40)=SumVec(40)+ FpForce_y(ic,jc,kc)*uyloc ! yp !
          SumVec(41)=SumVec(41)+ FpForce_z(ic,jc,kc)*uzloc
          SumVec(42)=SumVec(42)+ half*(FpForce_y(ic,jc,kc)+FpForce_y(ic,jp,kc))*uxCell
          
          uyCellm=(uy(im,jc,kc)+ uy(im,jp,kc))*half
          SumVec(42)=SumVec(42)+ FpForce_x(ic,jc,kc)*(uyCellm +uyCell)*half
#endif
        enddo
      enddo
      do kc=1,NCHASTAT
        SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
      enddo
    ENDDO

    ! nyp only
    jc=nyp; jm=jc-1; cac=rdyc(jc); SumVec=zero
    InterpY1= half*YinterpCoe(jc); InterpY2=half-InterpY1
    do kc=y1start(3),y1end(3)
      km=kc-1;kp=kc+1
      do ic=y1start(1),y1end(1)
        im=ic-1;ip=ic+1
        prloc2= InterpY1*(pressure(im,jm,kc)+pressure(ic,jm,kc))+ InterpY2*(pressure(im,jc,kc)+pressure(ic,jc,kc))
        dudy= (ux(ic,jc,kc)-ux(ic,jm,kc))*cac
        dwdy= (uz(ic,jc,kc)-uz(ic,jm,kc))*cac
        vor_x=  dwdy
        vor_z= -dudy
        SumVec(27)=SumVec(27)+ prloc2*(dudy+zero)   ! yp !
        SumVec(33)=SumVec(33)+ vor_x*vor_x          ! yp !
        SumVec(35)=SumVec(35)+ vor_z*vor_z          ! yp !
#if defined(CFDLPT_TwoWay)
        SumVec(37)=SumVec(37)+ FpForce_y(ic,jc,kc)              ! yp !
        SumVec(40)=SumVec(40)+ FpForce_y(ic,jc,kc)*uy(ic,jc,kc) ! yp !
#endif
      enddo
    enddo
    do kc=1,NCHASTAT
      SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
    enddo

    ! shear stress and pressure gradient
    PrGradsum   = PrGradsum+ PrGradAver
    if(nrank==0 .and. IsUxConst) then
      write(filename,'(A,I10.10)')trim(ResultsDir)//"PrGrad",ilast
      open(newunit=nUnit,file=filename,status='old',position='append',form='formatted',IOSTAT=ierror)
      if(ierror/=0) then
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      else
        write(nUnit,'(I7,2ES24.15)')itime,SimTime,PrGradAver
      endif
      close(nUnit,IOSTAT=ierror)
    endif
    nfstime= nfstime + 1
    if(mod(itime,SaveStat)/=0) return

    ! Write statistics
    call MPI_REDUCE(SumStat,SumStatR,NCHASTAT*nyp,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      infstime = one/real(nfstime,RK)
      write(filename,"(A,I10.10)") trim(ResultsDir)//'stats',itime
      open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
      IF(ierror/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      ELSE
        write(nUnit,'(a,I7,a,I7,a,I7)')'  The time step range for this fluid statistics is ', &
                                    itime-(nfstime-1)*ivstats, ':', ivstats, ':', itime
        write(nUnit,'(A)')'  '
        if(IsUxConst) then
          write(nUnit,'(A)')'  Constant velocity in x-dir by adding a pressure gradient.'
          write(nUnit,'(A, ES24.15)')'    time averaged pressure gradient is: ',PrGradsum*infstime
        else
          write(nUnit,'(A)')'  Variable velocity in x-dir while adding a constant body force.'
        endif
        write(nUnit,'(A)')'  '
        
        Block 
        character(len=128)::FormatStr
        write(FormatStr,'(A,I3,A)')'(',NCHASTAT,'ES24.15)'
        do jc=1,nyp
          write(nUnit,FormatStr)SumStatR(1:NCHASTAT,jc)*infstime
        enddo
        End block
      ENDIF
      close(nUnit,IOSTAT=ierror)
    endif
    nfstime=0; SumStat=zero; PrGradsum=zero;
  end subroutine clcStat_CH
#else

  !******************************************************************
  ! clcStat_CH
  !******************************************************************
  subroutine clcStat_CH(ux,uy,uz,pressure)
    implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz,pressure
   
    ! locals
    character(len=128)::filename
    integer::ic,jc,kc,ip,jp,kp,ierror,nUnit,Sum_indicator
    real(RK)::uxCell,uyCell,uzCell,prCell,inxz,infstime,sumStatR(NCHASTAT,nyc),SumVec(NCHASTAT)

    inxz=one/(real(nxc,RK)*real(nzc,RK))
    DO jc=y1start(2),y1end(2)
      jp=jc+1
      SumVec=zero; Sum_indicator=0
      do kc=y1start(3),y1end(3)
        kp=kc+1
        do ic=y1start(1),y1end(1)
          ip=ic+1
          uxCell= half*(ux(ic,jc,kc)+ux(ip,jc,kc))
          uyCell= half*(uy(ic,jc,kc)+uy(ic,jp,kc))
          uzCell= half*(uz(ic,jc,kc)+uz(ic,jc,kp))
          prCell= pressure(ic,jc,kc)
          SumVec( 1)=SumVec( 1)+ uxCell
          SumVec( 2)=SumVec( 2)+ uyCell
          SumVec( 3)=SumVec( 3)+ uzCell
          SumVec( 4)=SumVec( 4)+ uxCell*uxCell
          SumVec( 5)=SumVec( 5)+ uyCell*uyCell
          SumVec( 6)=SumVec( 6)+ uzCell*uzCell
          SumVec( 7)=SumVec( 7)+ uxCell*uyCell
          SumVec( 8)=SumVec( 8)+ prCell
          SumVec( 9)=SumVec( 9)+ prCell*prCell
          if(FluidIndicator(ic,jc,kc)=='P')cycle
          SumVec(10)=SumVec(10)+ uxCell
          SumVec(11)=SumVec(11)+ uyCell
          SumVec(12)=SumVec(12)+ uzCell
          SumVec(13)=SumVec(13)+ uxCell*uxCell
          SumVec(14)=SumVec(14)+ uyCell*uyCell
          SumVec(15)=SumVec(15)+ uzCell*uzCell
          SumVec(16)=SumVec(16)+ uxCell*uyCell
          SumVec(17)=SumVec(17)+ prCell
          SumVec(18)=SumVec(18)+ prCell*prCell
          Sum_indicator=Sum_indicator+1
        enddo
      enddo
      SumVec(NCHASTAT)=real(Sum_indicator,RK)
      do kc=1,NCHASTAT
        SumStat(kc,jc)=SumStat(kc,jc)+ SumVec(kc)*inxz
      enddo
    ENDDO

    ! shear stress and pressure graidient
    PrGradsum = PrGradsum +PrGradAver
    if(nrank==0 .and. IsUxConst) then
      write(filename,'(A,I10.10)')trim(ResultsDir)//"PrGrad",ilast
      open(newunit=nUnit,file=filename,status='old',position='append',form='formatted',IOSTAT=ierror)
      if(ierror /= 0) then
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      else
        write(nUnit,'(I7,2ES24.15)')itime,SimTime,PrGradAver
      endif
      close(nUnit,IOSTAT=ierror)
    endif
    nfstime = nfstime + 1
    if(mod(itime,SaveStat)/=0) return

    ! Write statistics
    call MPI_REDUCE(sumStat,sumStatR,NCHASTAT*nyc,real_type,MPI_SUM,0,MPI_COMM_WORLD,ierror)
    if(nrank==0) then
      infstime = one/real(nfstime,RK)
      write(filename,"(A,I10.10)")trim(ResultsDir)//'stats', itime
      open(newunit=nUnit,file=filename,status='replace',form='formatted',IOSTAT=ierror)
      IF(ierror/=0) THEN
        call MainLog%CheckForError(ErrT_Pass,"clcStat_CH","Cannot open file: "//trim(filename))
      ELSE
        write(nUnit,'(A,I7,A,I7,A,I7)')'  The time step range for this fluid statistics is ', &
                                  itime-(nfstime-1)*ivstats, ':', ivstats, ':', itime
        if(IsUxConst) then
          write(nUnit,'(A)'        )'  Constant velocity in x-dir by adding a pressure gradient.'
          write(nUnit,'(A,ES24.15)')'    time averaged pressure gradient is: ',PrGradsum*infstime
        else
          write(nUnit,'(A)')'  Variable velocity in x-dir while adding a constant body force.'
        endif
        write(nUnit,'(A)')'  '
        Block 
        character(len=128)::FormatStr
        write(FormatStr,'(A,I3,A)')'(',NCHASTAT,'ES24.15)'
        do jc=1,nyc
          write(nUnit,FormatStr)SumStatR(1:NCHASTAT,jc)*infstime
        enddo
        End block
      ENDIF
      close(nUnit,IOSTAT=ierror)
    endif
    nfstime=0; SumStat=zero; PrGradsum=zero;
  end subroutine clcStat_CH
#endif

#undef NCHASTAT
end module m_FlowType_Channel
