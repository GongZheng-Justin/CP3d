module m_Interp
  use MPI
  use m_TypeDef
  use m_LogInfo
  use m_Decomp2d
  use m_Parameters
  use m_BC_and_Halo
  use m_MeshAndMetries
  implicit none
  private

  integer, dimension(:),allocatable::indxp,indyp,indzp,indxc,indyc,indzc
  real(RK),dimension(:),allocatable::Ratxp,Ratyp,Ratzp,Ratxc,Ratyc,Ratzc
  type(DECOMP_INFO):: decomp_uxNew,decomp_uyNew,decomp_uzNew,decomp_prNew

  public:: Interp_init,Interp_uvwp
contains
 
  !**********************************************************************
  ! Interp_init
  !**********************************************************************
  subroutine Interp_init()
    implicit none

    ! locals
    integer::i,j,k,ic,jc,kc,js,je,intRed,ierror
    integer::xpSet,xpEnd,xpSize,xcSet,xcEnd,xcSize
    integer::zpSet,zpEnd,zpSize,zcSet,zcEnd,zcSize
    real(RK)::rdxOld,rdzOld,xstCoord,xedCoord,zstCoord,zedCoord

    rdxOld=real(nxcOld,RK);     rdzOld=real(nzcOld,RK)
    xstCoord= xpOld(yStart(1)); xedCoord= xpOld(yend(1)+1)
    zstCoord= zpOld(yStart(3)); zedCoord= zpOld(yend(3)+1)

    ! xp grid
    xpSet=1
    do i=1,nxcNew
      if(xpNew(i)<xstCoord .and. xpNew(i+1)>=xstCoord) then
         xpSet=i+1;exit
      endif
    enddo
    xpEnd=nxcNew
    do i=xpSet,nxcNew
      if(xpNew(i)<xedCoord .and. xpNew(i+1)>=xedCoord) then
         xpEnd=i;exit
      endif
    enddo
    xpSize=xpEnd-xpSet+1
    call MPI_REDUCE(xpSize,intRed,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror);intRed=intRed/p_col
    if(intRed/=nxcNew .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"Interp_init","xpSize wrong !!!") 

    ! xc grid
    xcSet=1
    do i=1,nxcNew
      if(xcNew(i)<xstCoord .and. xcNew(i+1)>=xstCoord) then
         xcSet=i+1;exit
      endif
    enddo
    xcEnd=nxcNew
    do i=xcSet,nxcNew
      if(xcNew(i)<xedCoord .and. xcNew(i+1)>=xedCoord) then
         xcEnd=i;exit
      endif
    enddo
    xcSize=xcEnd-xcSet+1
    call MPI_REDUCE(xcSize,intRed,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror);intRed=intRed/p_col
    if(intRed/=nxcNew .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"Interp_init","xcSize wrong !!!") 

    ! zp grid
    zpSet=1
    do k=1,nzcNew
      if(zpNew(k)<zstCoord .and. zpNew(k+1)>=zstCoord) then
         zpSet=k+1;exit
      endif
    enddo
    zpEnd=nzcNew
    do k=zpSet,nzcNew
      if(zpNew(k)<zedCoord .and. zpNew(k+1)>=zedCoord) then
         zpEnd=k;exit
      endif
    enddo
    zpSize=zpEnd-zpSet+1
    call MPI_REDUCE(zpSize,intRed,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror);intRed=intRed/p_row
    if(intRed/=nzcNew .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"Interp_init","zpSize wrong !!!") 

    ! zc grid
    zcSet=1
    do k=1,nzcNew
      if(zcNew(k)<zstCoord .and. zcNew(k+1)>=zstCoord) then
         zcSet=k+1;exit
      endif
    enddo
    zcEnd=nzcNew
    do k=zcSet,nzcNew
      if(zcNew(k)<zedCoord .and. zcNew(k+1)>=zedCoord) then
         zcEnd=k;exit
      endif
    enddo
    zcSize=zcEnd-zcSet+1
    call MPI_REDUCE(zcSize,intRed,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierror);intRed=intRed/p_row
    if(intRed/=nzcNew .and. nrank==0) call MainLog%CheckForError(ErrT_Abort,"Interp_init","zcSize wrong !!!")

    ! DECOMP_INFO for ux,uy,uz and pr grid respectively.
    decomp_uxNew%yst(1)=xpSet; decomp_uxNew%yst(2)=1;      decomp_uxNew%yst(3)=zcSet
    decomp_uxNew%yen(1)=xpEnd; decomp_uxNew%yen(2)=nycNew; decomp_uxNew%yen(3)=zcEnd
    decomp_uxNew%ysz=decomp_uxNew%yen-decomp_uxNew%yst+1

    decomp_uyNew%yst(1)=xcSet; decomp_uyNew%yst(2)=1;      decomp_uyNew%yst(3)=zcSet
    decomp_uyNew%yen(1)=xcEnd; decomp_uyNew%yen(2)=nycNew; decomp_uyNew%yen(3)=zcEnd
    decomp_uyNew%ysz=decomp_uyNew%yen-decomp_uyNew%yst+1     

    decomp_uzNew%yst(1)=xcSet; decomp_uzNew%yst(2)=1;      decomp_uzNew%yst(3)=zpSet
    decomp_uzNew%yen(1)=xcEnd; decomp_uzNew%yen(2)=nycNew; decomp_uzNew%yen(3)=zpEnd
    decomp_uzNew%ysz=decomp_uzNew%yen-decomp_uzNew%yst+1

    decomp_prNew%yst(1)=xcSet; decomp_prNew%yst(2)=1;      decomp_prNew%yst(3)=zcSet
    decomp_prNew%yen(1)=xcEnd; decomp_prNew%yen(2)=nycNew; decomp_prNew%yen(3)=zcEnd
    decomp_prNew%ysz=decomp_prNew%yen-decomp_prNew%yst+1

    ! interpolation coefficient 
    allocate(indxp(xpSet:xpEnd),Ratxp(xpSet:xpEnd))
    allocate(indxc(xcSet:xcEnd),Ratxc(xcSet:xcEnd))
    allocate(indzp(zpSet:zpEnd),Ratzp(zpSet:zpEnd))
    allocate(indzc(zcSet:zcEnd),Ratzc(zcSet:zcEnd))
    allocate(indyp(1:nycNew),Ratyp(1:nycNew),indyc(1:nycNew),Ratyc(1:nycNew))
    do i=xpSet,xpEnd
       ic= floor(xpNew(i)*rdxOld)+1; ic=min(ic,yend(1)); ic=max(ic,ystart(1))
       indxp(i)=ic; Ratxp(i)= (xpOld(ic+1)-xpNew(i))*rdxOld
    enddo
    do i=xcSet,xcEnd
       ic= floor(xcNew(i)*rdxOld)+1; ic=min(ic,yend(1)); ic=max(ic,ystart(1))
       if(xcNew(i)<xcOld(ic)) ic=ic-1;
       indxc(i)=ic; Ratxc(i)= (xcOld(ic+1)-xcNew(i))*rdxOld
    enddo
    do k=zpSet,zpEnd
       kc= floor(zpNew(k)*rdzOld)+1; kc=min(kc,yend(3)); kc=max(kc,ystart(3))
       indzp(k)=kc; Ratzp(k)= (zpOld(kc+1)-zpNew(k))*rdzOld
    enddo
    do k=zcSet,zcEnd
       kc= floor(zcNew(k)*rdzOld)+1; kc=min(kc,yend(3)); kc=max(kc,ystart(3))
       if(zcNew(k)<zcOld(kc)) kc=kc-1;
       indzc(k)=kc; Ratzc(k)= (zcOld(kc+1)-zcNew(k))*rdzOld
    enddo

    do j=1,nycNew
      ! if pos%y is within [0,yly), jc will be within [1,nyc]
      js=0; je=nypOld+1
      do
        jc=(js+je)/2
        if(je-js==1) exit
        if(ypNew(j)< ypOld(jc)) then
          je =jc
        else
          js =jc
        endif
      enddo
      if(jc>nycOld) call MainLog%CheckForError(ErrT_Abort,"Interp_init","wrong jc 1")
      indyp(j)=jc; Ratyp(j)=(ypOld(jc+1)-ypNew(j))/(ypOld(jc+1)-ypOld(jc))
    enddo

    do j=1,nycNew
      ! if pos%y is within [0,yly), jc will be within [1,nyc]
      js=0; je=nypOld+1
      do
        jc=(js+je)/2
        if(je-js==1) exit
        if(ycNew(j)< ypOld(jc)) then
          je =jc
        else
          js =jc
        endif
      enddo
      if(jc>nycOld) call MainLog%CheckForError(ErrT_Abort,"Interp_init","wrong jc 2")

      if(ycNew(j)<ycOld(jc)) jc=jc-1
      indyc(j)=jc; Ratyc(j)=(ycOld(jc+1)-ycNew(j))/(ycOld(jc+1)-ycOld(jc))
    enddo

  end subroutine Interp_init

  !**********************************************************************
  ! Interp_uvwp
  !**********************************************************************
  subroutine Interp_uvwp(mbOld,ArrOld)
    implicit none
    type(MatBound),intent(in)::mbOld
    real(RK),dimension(mbOld%xmm:mbOld%xpm,mbOld%ymm:mbOld%ypm,mbOld%zmm:mbOld%zpm),intent(inout)::ArrOld
   
    ! locals
    character(256)::chFileOld,chFileNew
    integer(MPI_OFFSET_KIND)::dispOld,dispNew
    real(RK),dimension(:,:,:),allocatable::ArrNew
    real(RK)::prx1,pry1,prz1,prx2,pry2,prz2,MeanValue,MeanValueR
    integer::idx1,idy1,idz1,idx2,idy2,idz2,ic,jc,kc,ierror,fhOld,fhNew
    
    ! open "RestartOld" and "RestartNew"
    write(chFileOld,'(A)')OldRestartName
    write(chFileNew,'(A)')NewRestartName
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFileNew, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fhNew, ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, chFileOld, MPI_MODE_RDONLY,                 MPI_INFO_NULL, fhOld, ierror)
    if(ierror /= 0 .and. nrank==0) then
      call MainLog%CheckForError(ErrT_Abort,"Interp_uvwp","Cannot open RestartOld: "//trim(chFileOld))
    endif
    dispNew=0_MPI_OFFSET_KIND; dispOld=0_MPI_OFFSET_KIND 

    ! ux grid  =============
    allocate(ArrNew(decomp_uxNew%yst(1):decomp_uxNew%yen(1),decomp_uxNew%yst(2):decomp_uxNew%yen(2),decomp_uxNew%yst(3):decomp_uxNew%yen(3)))
    asso_ux: associate( uxOld=>ArrOld, uxNew=>ArrNew); uxOld=zero
    call decomp2d_readVar(fhOld,dispOld,uxOld(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    call SetBC_and_UpdateHalo_ux(uxOld,mbOld)
    do kc=decomp_uxNew%yst(3),decomp_uxNew%yen(3)
      idz1=indzc(kc);  prz1=Ratzc(kc) 
      idz2=idz1+1;     prz2=one-prz1    
      do jc=decomp_uxNew%yst(2),decomp_uxNew%yen(2)
        idy1=indyc(jc);  pry1=Ratyc(jc)
        idy2=idy1+1;     pry2=one-pry1
        do ic=decomp_uxNew%yst(1),decomp_uxNew%yen(1)
          idx1=indxp(ic);  prx1=Ratxp(ic)
          idx2=idx1+1;     prx2=one-prx1
          uxNew(ic,jc,kc)= uxOld(idx1,idy1,idz1)*prx1*pry1*prz1 +uxOld(idx2,idy1,idz1)*prx2*pry1*prz1 &
                          +uxOld(idx1,idy2,idz1)*prx1*pry2*prz1 +uxOld(idx2,idy2,idz1)*prx2*pry2*prz1 &
                          +uxOld(idx1,idy1,idz2)*prx1*pry1*prz2 +uxOld(idx2,idy1,idz2)*prx2*pry1*prz2 &
                          +uxOld(idx1,idy2,idz2)*prx1*pry2*prz2 +uxOld(idx2,idy2,idz2)*prx2*pry2*prz2 
        enddo
      enddo
    enddo
    if(IsUxConst) then
      MeanValue=zero
      do kc=decomp_uxNew%yst(3),decomp_uxNew%yen(3)
        do jc=decomp_uxNew%yst(2),decomp_uxNew%yen(2)
          do ic=decomp_uxNew%yst(1),decomp_uxNew%yen(1)
            MeanValue=MeanValue+uxNew(ic,jc,kc)*(ypNew(jc+1)-ypNew(jc))
          enddo
        enddo
      enddo
      MeanValue=MeanValue/(real(nxcNew,RK))/(real(nzcNew,RK))/yly
      call MPI_ALLREDUCE(MeanValue,MeanValueR,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      uxNew= (ubulk/MeanValueR)*uxNew
    endif
    call decomp2d_writeVar(fhNew,dispNew,uxNew,decomp_uxNew)
    end associate asso_ux
    deallocate(ArrNew)

    ! uy grid  =============
    allocate(ArrNew(decomp_uyNew%yst(1):decomp_uyNew%yen(1),decomp_uyNew%yst(2):decomp_uyNew%yen(2),decomp_uyNew%yst(3):decomp_uyNew%yen(3)))
    asso_uy: associate( uyOld=>ArrOld, uyNew=>ArrNew); uyOld=zero
    call decomp2d_readVar(fhOld,dispOld,uyOld(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    call SetBC_and_UpdateHalo_uy(uyOld,mbOld)
    do kc=decomp_uyNew%yst(3),decomp_uyNew%yen(3)
      idz1=indzc(kc);  prz1=Ratzc(kc) 
      idz2=idz1+1;     prz2=one-prz1    
      do jc=decomp_uyNew%yst(2),decomp_uyNew%yen(2)
        idy1=indyp(jc);  pry1=Ratyp(jc)
        idy2=idy1+1;     pry2=one-pry1
        do ic=decomp_uyNew%yst(1),decomp_uyNew%yen(1)
          idx1=indxc(ic);  prx1=Ratxc(ic)
          idx2=idx1+1;     prx2=one-prx1
          uyNew(ic,jc,kc)= uyOld(idx1,idy1,idz1)*prx1*pry1*prz1 +uyOld(idx2,idy1,idz1)*prx2*pry1*prz1 &
                          +uyOld(idx1,idy2,idz1)*prx1*pry2*prz1 +uyOld(idx2,idy2,idz1)*prx2*pry2*prz1 &
                          +uyOld(idx1,idy1,idz2)*prx1*pry1*prz2 +uyOld(idx2,idy1,idz2)*prx2*pry1*prz2 &
                          +uyOld(idx1,idy2,idz2)*prx1*pry2*prz2 +uyOld(idx2,idy2,idz2)*prx2*pry2*prz2 
        enddo
      enddo
    enddo
    call decomp2d_writeVar(fhNew,dispNew,uyNew,decomp_uyNew)
    end associate asso_uy
    deallocate(ArrNew)

    ! uz grid  =============
    allocate(ArrNew(decomp_uzNew%yst(1):decomp_uzNew%yen(1),decomp_uzNew%yst(2):decomp_uzNew%yen(2),decomp_uzNew%yst(3):decomp_uzNew%yen(3)))
    asso_uz: associate( uzOld=>ArrOld, uzNew=>ArrNew); uzOld=zero
    call decomp2d_readVar(fhOld,dispOld,uzOld(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    call SetBC_and_UpdateHalo_uz(uzOld,mbOld)
    do kc=decomp_uzNew%yst(3),decomp_uzNew%yen(3)
      idz1=indzp(kc);  prz1=Ratzp(kc) 
      idz2=idz1+1;     prz2=one-prz1    
      do jc=decomp_uzNew%yst(2),decomp_uzNew%yen(2)
        idy1=indyc(jc);  pry1=Ratyc(jc)
        idy2=idy1+1;     pry2=one-pry1
        do ic=decomp_uzNew%yst(1),decomp_uzNew%yen(1)
          idx1=indxc(ic);  prx1=Ratxc(ic)
          idx2=idx1+1;     prx2=one-prx1
          uzNew(ic,jc,kc)= uzOld(idx1,idy1,idz1)*prx1*pry1*prz1 +uzOld(idx2,idy1,idz1)*prx2*pry1*prz1 &
                          +uzOld(idx1,idy2,idz1)*prx1*pry2*prz1 +uzOld(idx2,idy2,idz1)*prx2*pry2*prz1 &
                          +uzOld(idx1,idy1,idz2)*prx1*pry1*prz2 +uzOld(idx2,idy1,idz2)*prx2*pry1*prz2 &
                          +uzOld(idx1,idy2,idz2)*prx1*pry2*prz2 +uzOld(idx2,idy2,idz2)*prx2*pry2*prz2 
        enddo
      enddo
    enddo
    call decomp2d_writeVar(fhNew,dispNew,uzNew,decomp_uzNew)
    end associate asso_uz
    deallocate(ArrNew)

    ! pr grid  =============
    allocate(ArrNew(decomp_prNew%yst(1):decomp_prNew%yen(1),decomp_prNew%yst(2):decomp_prNew%yen(2),decomp_prNew%yst(3):decomp_prNew%yen(3)))
    asso_pr: associate( prOld=>ArrOld, prNew=>ArrNew); prOld=zero
    call decomp2d_readVar(fhOld,dispOld,prOld(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    call SetBC_and_UpdateHalo_pr(prOld,mbOld)
    do kc=decomp_prNew%yst(3),decomp_prNew%yen(3)
      idz1=indzc(kc);  prz1=Ratzc(kc) 
      idz2=idz1+1;     prz2=one-prz1    
      do jc=decomp_prNew%yst(2),decomp_prNew%yen(2)
        idy1=indyc(jc);  pry1=Ratyc(jc)
        idy2=idy1+1;     pry2=one-pry1
        do ic=decomp_prNew%yst(1),decomp_prNew%yen(1)
          idx1=indxc(ic);  prx1=Ratxc(ic)
          idx2=idx1+1;     prx2=one-prx1
          prNew(ic,jc,kc)= prOld(idx1,idy1,idz1)*prx1*pry1*prz1 +prOld(idx2,idy1,idz1)*prx2*pry1*prz1 &
                          +prOld(idx1,idy2,idz1)*prx1*pry2*prz1 +prOld(idx2,idy2,idz1)*prx2*pry2*prz1 &
                          +prOld(idx1,idy1,idz2)*prx1*pry1*prz2 +prOld(idx2,idy1,idz2)*prx2*pry1*prz2 &
                          +prOld(idx1,idy2,idz2)*prx1*pry2*prz2 +prOld(idx2,idy2,idz2)*prx2*pry2*prz2 
        enddo
      enddo
    enddo
    call decomp2d_writeVar(fhNew,dispNew,prNew,decomp_prNew)
    end associate asso_pr
    deallocate(ArrNew)

#ifdef ScalarFlow
    ! saclar grid  =============
    allocate(ArrNew(decomp_prNew%yst(1):decomp_prNew%yen(1),decomp_prNew%yst(2):decomp_prNew%yen(2),decomp_prNew%yst(3):decomp_prNew%yen(3)))
    asso_scalar: associate( scalarOld=>ArrOld, scalarNew=>ArrNew); scalarOld=zero
    call decomp2d_readVar(fhOld,dispOld,scalarOld(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
    call SetBC_and_UpdateHalo_scalar(scalarOld,mbOld)
    do kc=decomp_prNew%yst(3),decomp_prNew%yen(3)
      idz1=indzc(kc);  prz1=Ratzc(kc) 
      idz2=idz1+1;     prz2=one-prz1    
      do jc=decomp_prNew%yst(2),decomp_prNew%yen(2)
        idy1=indyc(jc);  pry1=Ratyc(jc)
        idy2=idy1+1;     pry2=one-pry1
        do ic=decomp_prNew%yst(1),decomp_prNew%yen(1)
          idx1=indxc(ic);  prx1=Ratxc(ic)
          idx2=idx1+1;     prx2=one-prx1
          scalarNew(ic,jc,kc)= scalarOld(idx1,idy1,idz1)*prx1*pry1*prz1 +scalarOld(idx2,idy1,idz1)*prx2*pry1*prz1 &
                              +scalarOld(idx1,idy2,idz1)*prx1*pry2*prz1 +scalarOld(idx2,idy2,idz1)*prx2*pry2*prz1 &
                              +scalarOld(idx1,idy1,idz2)*prx1*pry1*prz2 +scalarOld(idx2,idy1,idz2)*prx2*pry1*prz2 &
                              +scalarOld(idx1,idy2,idz2)*prx1*pry2*prz2 +scalarOld(idx2,idy2,idz2)*prx2*pry2*prz2 
        enddo
      enddo
    enddo
    if(IsScalarConst) then
      MeanValue=zero
      do kc=decomp_prNew%yst(3),decomp_prNew%yen(3)
        do jc=decomp_prNew%yst(2),decomp_prNew%yen(2)
          do ic=decomp_prNew%yst(1),decomp_prNew%yen(1)
            MeanValue=MeanValue+scalarNew(ic,jc,kc)*(ypNew(jc+1)-ypNew(jc))
          enddo
        enddo
      enddo
      MeanValue=MeanValue/(real(nxcNew,RK))/(real(nzcNew,RK))/yly
      call MPI_ALLREDUCE(MeanValue,MeanValueR,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierror)
      scalarNew= (ScalarMean/MeanValueR)*scalarNew
    endif
    call decomp2d_writeVar(fhNew,dispNew,scalarNew,decomp_prNew)
    end associate asso_scalar
    deallocate(ArrNew)
#endif

    call MPI_FILE_CLOSE(fhOld,ierror)
    call MPI_FILE_CLOSE(fhNew,ierror)
  end subroutine Interp_uvwp

end module m_Interp
