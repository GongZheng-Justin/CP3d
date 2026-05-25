#include "definitions_inc.f90"
module ca_IBM
  use MPI
  use mc_TypeDef
  use mc_LogInfo
  use mc_EqualSphere
  use mc_Decomp2d,only:nrank,nproc,y1start,y1end,real_type,myProcNghBC
  use f2_Parameters
  use f2_Variables,only:mb1
  use f2_MeshAndMetries,only: dx,dyUniform,dz,rdx,rdyUniform,rdz,xc,yc,zc,dyp
  use ca_BC_and_Halo,only: Gather_Halo_IBMForce, SetBC_and_UpdateHalo_VolForce
  implicit none
  private
#ifdef IBMDistributeLinear
#define nDistribute 1
#else
#define nDistribute 2
#endif
#define THIS_PROC_XP 7
#define THIS_PROC_XM 8
#define THIS_PROC_YP 9
#define THIS_PROC_YM 10
#define THIS_PROC_ZP 11
#define THIS_PROC_ZM 12

  integer:: int_type, real3_type, real4_type
  integer:: int_byte, real3_byte, real4_byte, real_byte
  
  type part_io_size_vec
    integer,dimension(1)::sizes
    integer,dimension(1)::subsizes
    integer,dimension(1)::starts
  end type part_io_size_vec
  type part_io_size_mat
    integer,dimension(2)::sizes
    integer,dimension(2)::subsizes
    integer,dimension(2)::starts
  end type part_io_size_mat
    
  integer,parameter:: IK = 4
  character(len=:),allocatable,public:: SpheresInfoFile_
  integer :: nlocal_Big, numPrtcl_Big
  real(RK)::RetractionRatio_1 =0.0_RK
  real(RK)::RetractionRatio_2 =0.0_RK
  integer,dimension(:),allocatable::GPrtcl_Big_id
  type(real4),dimension(:),allocatable::GPrtcl_Big_PosR
  type(real3),dimension(:),allocatable::GPrtcl_Big_FpForce
  type(real3),dimension(:),allocatable::GPrtcl_Big_FpTorque
  integer,dimension(:,:),allocatable:: GPrtcl_Big_ind_Sml
  
  integer :: nlocal_Sml, numPrtcl_Sml
  type(real4),dimension(:),allocatable::GPrtcl_Sml_PosR
  
  ! Ghost particles
  integer::nGhost_Big,nGhost_Sml,msend
  integer,dimension(:),allocatable::GhostP_Big_id
  integer,dimension(:),allocatable::GhostP_Big_id_global
  integer,dimension(:),allocatable::GhostP_Big_nSmall
  integer,dimension(:),allocatable::GhostP_Big_Direction
  type(real4),dimension(:),allocatable::GhostP_Big_PosR
  type(real3),dimension(:),allocatable::GhostP_Big_FpForce
  type(real3),dimension(:),allocatable::GhostP_Big_FpTorque
  integer,dimension(:,:),allocatable:: GhostP_Big_ind_Sml
  type(real4),dimension(:),allocatable::GhostP_Sml_PosR
  integer,dimension(:),allocatable::sendlist
  
  real(RK):: CellRadius,CellVolumeIBM,dxhalf,dyhalf,dzhalf,SMALL
  real(RK):: xstCoord,xedCoord,ystCoord,yedCoord,zstCoord,zedCoord

  ! Immersed boundary points parts
  integer:: nIBP         ! number of Immersed Boundary Point in local processor
  integer:: mIBP         ! the possible maxium # of Immersed Boundary Point in local processor
  real(RK),dimension(:),allocatable::    IBP_VolRatio
  integer,dimension(:),allocatable::     IBP_idlocal
  integer(kind=2),dimension(:,:),allocatable::IBP_indxyz
  type(real3),dimension(:),allocatable:: IBP_Pos
  type(real3),dimension(:),allocatable:: IBP_Vel
  type(real3),dimension(:),allocatable:: IBP_Force
  type(real3),dimension(:),allocatable:: IBP_ForceAmplify

  ABSTRACT INTERFACE
    function deltaFunction_(ratio_in) result(delta)
      use mc_TypeDef,only:RK; implicit none
      real(RK),intent(in)::ratio_in
      real(RK):: delta
    end function deltaFunction_
  END INTERFACE
  procedure(deltaFunction_),pointer::deltaFunction
  procedure(),pointer::AdditionalForceIBM

  integer:: IBMForceScheme = 0
  integer:: nForcingExtra  = 2
  real(RK)::ForceAmplifyCoe= 1.0_RK

  ! useful interfaces
  interface Prtcl_dump
    module procedure Prtcl_dump_int_vector,  Prtcl_dump_int_matrix
    module procedure Prtcl_dump_real_vector, Prtcl_dump_real3_vector
  end interface Prtcl_dump
    
  ! public variables and functions/subroutines
  public:: nForcingExtra, Prtcl_Init_Visu, Prtcl_Dump_Visu
  public:: Update_FluidIndicator, IntegrateFluidPrtclForce, PrepareIBM_IbpForce
  public:: Init_IBM, updateRhsIBM, prepareIBM_interp, Clc_NoSlipErr, AdditionalForceIBM
contains
#include "pIBM_IBM_common_inc.f90"
#include "Prtcl_Dump_MPI_inc.f90"

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PrepareIBM_interp
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PrepareIBM_interp(VolForce_x,VolForce_y,VolForce_z); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::VolForce_x,VolForce_y,VolForce_z
    
    ! locals
    type(real3):: Cposition,LagrangeP,PosDiff
    real(RK),allocatable,dimension(:,:)::LagrangePoints
    real(RK)::rTmp,dxyz,maxR,max_dxyz,min_xz,IbpVolRatio
    real(RK)::xBigO,yBigO,zBigO, xBig,yBig,zBig, xSml,ySml,zSml, dTmp
    integer:: j,k,iUnit,ierr,ind,ind_old,indt, ierrTmp,nPartition,nIbpSum1,nIbpSum2,nIbpTmp,ids_b, ids_e, np_Big, np_Sml
    integer:: i,ic,jc,kc,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    integer,allocatable,dimension(:)::GhostP_Big_id_global_sort, GhostP_Big_sort_id
            
    dxyz= CellVolumeIBM**(0.33333333333333333333_RK)
    
    !NOTE (Gong Zheng, 2020/07/01):
    ! [1] Here ONLY the IBP points within the Processor physical Domains will be considered.
    !       The IBP points outsides will be skipped.
    ! [2] ONLY the uniform meshes near y-dir are used.
    nIBP= 0
    open(newunit=iUnit, file=SpheresInfoFile_,form='formatted',status='old', action='read', iostat=ierr)
    if(ierr /=0) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","Cannot open file:"//strip(SpheresInfoFile_))
    
    ! Check consistency 1st
    if(nrank==0) then
      xBigO=-1.0E+29_RK; yBigO=-1.0E+29_RK; zBigO=-1.0E+29_RK; ind_old=-314159265
      DO
        read(unit=iUnit, fmt=*, iostat=ierr) ind, xBig,yBig,zBig, xSml,ySml,zSml, dTmp
        if(ierr /= 0) exit
        if(ind ==ind_old .and. (xBig /= xBigO .or. yBig /= yBigO .or. zBig /= zBigO)) then
          call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","Big coordinates wrong")
        endif
        ind_old=ind; xBigO=xBig; yBigO=yBig; zBigO=zBig;
        rTmp=sqrt((xBig-xSml)**2 +(yBig-ySml)**2 +(zBig-zSml)**2)
        if(rTmp > 0.5_RK*xlx .or. rTmp > 0.5_RK*yly .or. rTmp > 0.5_RK*zlz) then
          call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","So big distance between ensemble sphere and small sphere")        
        endif
      ENDDO
      rewind(unit=iUnit, iostat=ierr)
    endif
   
    ! Count nlocal_Big
    nlocal_Big = 0
    maxR = -1.0_RK
    ind_old = -314159265
    DO
      read(unit=iUnit, fmt=*, iostat=ierr) ind, xBig,yBig,zBig, xSml,ySml,zSml, dTmp
      if(ierr /= 0) exit
      rTmp=sqrt((xBig-xSml)**2 +(yBig-ySml)**2 +(zBig-zSml)**2) +dTmp*0.5_RK
      if(maxR < rTmp) maxR = rTmp        
      if(ind ==ind_old) cycle
      ind_old = ind
      if(xBig < xstCoord .or. xBig >=xedCoord .or. yBig < ystCoord .or. &
         yBig >=yedCoord .or. zBig < zstCoord .or. zBig >=zedCoord) cycle
      nlocal_Big=nlocal_Big+1          
    ENDDO
    max_dxyz=max(max(dx,dyUniform),dz)
    min_xz=min(xedCoord-xstCoord,zedCoord-zstCoord)
    if(min_xz<=maxR+ 2.0_RK*max_dxyz) then
      call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","so big Diameter")
    endif
    
    if(nlocal_Big>0) then
      ierr=0
      allocate(GPrtcl_Big_id(nlocal_Big), stat=ierrTmp);        ierr = ierr + ierrTmp
      allocate(GPrtcl_Big_PosR(nlocal_Big), stat=ierrTmp);      ierr = ierr + ierrTmp
      allocate(GPrtcl_Big_FpForce(nlocal_Big), stat=ierrTmp);   ierr = ierr + ierrTmp
      allocate(GPrtcl_Big_FpTorque(nlocal_Big), stat=ierrTmp);  ierr = ierr + ierrTmp
      allocate(GPrtcl_Big_ind_Sml(2,nlocal_Big), stat=ierrTmp); ierr = ierr + ierrTmp
      if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","Allocation failed 1")
      GPrtcl_Big_FpForce = zero_r3
      GPrtcl_Big_FpTorque= zero_r3
    endif

    nlocal_Big =0; k=0; np_Sml=0; ids_e =0
    ind_old = -314159265;  np_Big=0; indt=ind_old
    rewind(unit=iUnit, iostat=ierr)
    DO
      read(unit=iUnit, fmt=*, iostat=ierr) ind, xBig,yBig,zBig, xSml,ySml,zSml, dTmp
      if(ierr /= 0) exit
      np_Sml = np_Sml+1

      if(ind /= indt) then
        np_Big=np_Big+1; indt=ind
      endif
      if(xBig < xstCoord .or. xBig >=xedCoord .or. yBig < ystCoord .or. &
         yBig >=yedCoord .or. zBig < zstCoord .or. zBig >=zedCoord) cycle
         
      if(ind /= ind_old) then
        ids_e = 1
        nlocal_Big = nlocal_Big + 1
      else
        ids_e = ids_e +1
      endif
      ind_old = ind

      if(nlocal_Big <0 .or. np_Big<0) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","nlocal_Big/np_Big wrong")
      GPrtcl_Big_id(nlocal_Big)   = np_Big
      GPrtcl_Big_PosR(nlocal_Big) = real3(xBig,yBig,zBig)        
      GPrtcl_Big_ind_Sml(2, nlocal_Big) = ids_e
    ENDDO
    call MPI_ALLREDUCE(nlocal_Big, numPrtcl_Big, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if(np_Big /= numPrtcl_Big .and. nrank==0)  call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","numPrtcl_Big wrong")
    
    ids_b =1
    do k=1,nlocal_Big
      GPrtcl_Big_ind_Sml(1, k)= ids_b
      ids_b = ids_b+GPrtcl_Big_ind_Sml(2, k)
      ids_e = ids_b-1
      GPrtcl_Big_ind_Sml(2, k)= ids_e
    enddo
    if(nlocal_Big>0) then
      nlocal_Sml= GPrtcl_Big_ind_Sml(2,nlocal_Big)
    else
      nlocal_Sml=0
    endif   
    call MPI_ALLREDUCE(nlocal_Sml, numPrtcl_Sml, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if(np_Sml /= numPrtcl_Sml .and. nrank==0)  call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","numPrtcl_Sml wrong")
        
    ! Allocate small particle properties
    if(nlocal_Sml >0) then
      ierr=0
      allocate(GPrtcl_Sml_PosR(nlocal_Sml), stat=ierrTmp); ierr = ierr + ierrTmp
      if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","Allocation failed 2")
    endif
    k=0
    rewind(unit=iUnit, iostat=ierr)
    DO
      read(unit=iUnit, fmt=*, iostat=ierr) ind, xBig,yBig,zBig, xSml,ySml,zSml, dTmp
      if(ierr /= 0) exit
      if(xBig < xstCoord .or. xBig >=xedCoord .or. yBig < ystCoord .or. &
         yBig >=yedCoord .or. zBig < zstCoord .or. zBig >=zedCoord) cycle
      k=k+1
      GPrtcl_Sml_PosR(k) = real4(xSml,ySml,zSml, 0.5_RK*dTmp)      
    ENDDO
    if(k /= nlocal_Sml) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","nlocal_Sml wrong")
    
    ! Compute the big diameters and check consistency 2nd
    block
    real(RK)::Coord_low(3), Coord_high(3)
    DO i=1,nlocal_Big
      maxR = -1.0_RK
      Coord_low = +1.00E+99
      Coord_high= -1.00E+99
      xBig = GPrtcl_Big_PosR(i)%x
      yBig = GPrtcl_Big_PosR(i)%y
      zBig = GPrtcl_Big_PosR(i)%z
      do k=GPrtcl_Big_ind_Sml(1,i), GPrtcl_Big_ind_Sml(2,i)
        xSml=GPrtcl_Sml_PosR(k)%x; Coord_low(1)=min(Coord_low(1),xSml); Coord_high(1)=max(Coord_high(1),xSml)
        ySml=GPrtcl_Sml_PosR(k)%y; Coord_low(2)=min(Coord_low(2),ySml); Coord_high(2)=max(Coord_high(2),ySml)
        zSml=GPrtcl_Sml_PosR(k)%z; Coord_low(3)=min(Coord_low(3),zSml); Coord_high(3)=max(Coord_high(3),zSml)        
        rTmp=sqrt((xBig-xSml)**2 +(yBig-ySml)**2 +(zBig-zSml)**2) +GPrtcl_Sml_PosR(k)%w
        if(maxR < rTmp) maxR = rTmp
      enddo
      GPrtcl_Big_PosR(i)%w = maxR + SMALL
      if(GPrtcl_Big_ind_Sml(1,i) == GPrtcl_Big_ind_Sml(2,i)) then
        if(xBig /= xSml .or. yBig /= ySml .or. zBig /= zSml) then
          call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","big/small coordinates not consistent")
        endif
      else
        if(xBig < Coord_low(1)  .or. xBig > Coord_high(1) .or. yBig < Coord_low(2) .or. &
           yBig > Coord_high(2) .or. zBig < Coord_low(3)  .or. zBig > Coord_high(3)) then
          call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","Big Sphere center wrong") 
        endif
      endif
    ENDDO
    endblock
        
    call PComm_IBM_forward()
    write(*, '(A, I4, A, I7, A, I8, A, I5, A)') 'rank',nrank,' contains',nlocal_Big,' multispheres,',nlocal_Sml, &
                                                ' small spheres, and',nGhost_Big,' ghost multispheres'
                                                    
    ! Compute properties for ghost small particles
    block
    integer::pid, pjd, intTmp
    if(nGhost_Big>0) then
      allocate(GhostP_Big_id_global_sort(nGhost_Big), GhostP_Big_sort_id(nGhost_Big), stat=ierr)
      if(ierr /=0) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","Allocation failed 3")
      
      do pid=1, nGhost_Big
        GhostP_Big_id_global_sort(pid) = GhostP_Big_id_global(pid)
        GhostP_Big_sort_id(pid) = pid
      enddo
      
      do pid=1, nGhost_Big-1
        do pjd=pid+1, nGhost_Big
          if(GhostP_Big_id_Global_sort(pid) <= GhostP_Big_id_Global_sort(pjd)) cycle
          
          intTmp = GhostP_Big_id_Global_sort(pjd)
          GhostP_Big_id_Global_sort(pjd) = GhostP_Big_id_Global_sort(pid)
          GhostP_Big_id_Global_sort(pid) = intTmp
          
          intTmp = GhostP_Big_sort_id(pjd)
          GhostP_Big_sort_id(pjd) = GhostP_Big_sort_id(pid)
          GhostP_Big_sort_id(pid) = intTmp      
        enddo
      enddo
    endif
    end block
    
    ids_e= 0                                  
    ind_old = -314159265;  np_Big=0;  indt=ind_old
    rewind(unit=iUnit, iostat=ierr)
    Block
    logical::FindGhostGlobal
    integer::id_Ghost, nGhost_Big_Tmp, iRepeat, kt
    id_Ghost = 0
    nGhost_Big_Tmp = 0
    FindGhostGlobal= .false.
    if(nGhost_Big > 0) then
      DO
        read(unit=iUnit, fmt=*, iostat=ierr) ind, xBig,yBig,zBig, xSml,ySml,zSml, dTmp
        if(ierr /= 0) exit

        if(ind /= indt) then
          np_Big=np_Big+1; indt=ind
          FindGhostGlobal = .false.
          do k=id_Ghost+1, nGhost_Big
            if(GhostP_Big_id_global_sort(k) /= np_Big) cycle
            iRepeat = 1 ! Note here 
            id_Ghost= k
            do kt= k+1, nGhost_Big
              if(GhostP_Big_id_global_sort(kt) /= np_Big) exit
              iRepeat = iRepeat + 1
            enddo
            FindGhostGlobal= .true.
            exit
          enddo
        endif
        if(.not. FindGhostGlobal) cycle
        
        if(ind /= ind_old) then
          ids_e = 1
          nGhost_Big_Tmp = nGhost_Big_Tmp + iRepeat
        else
          ids_e = ids_e + 1
        endif
        ind_old = ind
        
        if(nGhost_Big_Tmp <0 .or. np_Big<0) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","nGhost_Big_Tmp/np_Big wrong")
        GhostP_Big_ind_Sml(1, id_Ghost:nGhost_Big_Tmp)=ids_e        
      ENDDO
    endif
    if(nGhost_Big_Tmp /= nGhost_Big) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","nGhost_Big_Tmp wrong")
    End Block 
    
    do k=1,nGhost_Big
      kc = GhostP_Big_sort_id(k)
      GhostP_Big_ind_Sml(2, kc) = GhostP_Big_ind_Sml(1, k)
    enddo
    
    ids_b =1
    do k=1,nGhost_Big
      GhostP_Big_ind_Sml(1, k)= ids_b
      ids_b = ids_b+GhostP_Big_ind_Sml(2, k)
      ids_e = ids_b-1
      GhostP_Big_ind_Sml(2, k)= ids_e
    enddo
    if(nGhost_Big>0) then
      nGhost_Sml= GhostP_Big_ind_Sml(2,nGhost_Big)
    else
      nGhost_Sml= 0
    endif
    Block
    integer::nGhost_Sml_Tmp
    nGhost_Sml_Tmp=0
    do k=1,nGhost_Big
      nGhost_Sml_Tmp = nGhost_Sml_Tmp + GhostP_Big_nSmall(k)
    enddo
    if(nGhost_Sml /= nGhost_Sml_Tmp .and. nrank==0)  call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","nGhost_Sml_Tmp wrong")
    end block
            
    ! Fill ghost small particles' position
    rewind(unit=iUnit, iostat=ierr)
    Block
    logical::FindGhostGlobal
    real(RK)::xDiff, yDiff, zDiff
    integer::id_Ghost, nGhost_Sml_Tmp, iAppendix, kt, kid, nBase
    nGhost_Sml_Tmp=0
    ind_old = -314159265;  np_Big=0; indt = ind_old
    id_Ghost = 0;
    FindGhostGlobal= .false.    
    if(nGhost_Sml > 0) then
      ierr=0
      allocate(GhostP_Sml_PosR(nGhost_Sml), stat=ierrTmp); ierr = ierr + ierrTmp
      if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","Allocation failed 4")
      DO
        read(unit=iUnit, fmt=*, iostat=ierr) ind, xBig,yBig,zBig, xSml,ySml,zSml, dTmp
        if(ierr /= 0) exit

        if(ind /= indt) then
          np_Big=np_Big+1; indt=ind
          FindGhostGlobal = .false.
          do k=id_Ghost+1, nGhost_Big
            if (GhostP_Big_id_global_sort(k) /= np_Big) cycle
            iAppendix = 0  ! Note here 
            id_Ghost= k
            do kt= k+1, nGhost_Big
              if (GhostP_Big_id_global_sort(kt) /= np_Big) exit
              iAppendix = iAppendix + 1
            enddo
            FindGhostGlobal = .true.            
            exit
          enddo
        endif
        if (.not. FindGhostGlobal) cycle
        
        if (ind /= ind_old) then
          ids_e = 1
        else
          ids_e = ids_e + 1
        endif
        ind_old = ind
                
        kid = GhostP_Big_sort_id(id_Ghost)
        kt = GhostP_Big_nSmall(kid)        
        do k = 0, iAppendix
          kc = k + id_Ghost
          kid = GhostP_Big_sort_id(kc)
          nBase = GhostP_Big_ind_Sml(1, kid)-1

          if (kt /= GhostP_Big_nSmall(kid) ) &
            call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","GhostP_Big_nSmall wrong")
                     
          xDiff = GhostP_Big_PosR(kid)%x - xBig
          yDiff = GhostP_Big_PosR(kid)%y - yBig
          zDiff = GhostP_Big_PosR(kid)%z - zBig
          if(abs(abs(xDiff/xlx) -1.0_RK)> 1.0E-10_RK .and. abs(xDiff/xlx)> 1.0E-10_RK) &
            call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","xDiff wrong")
          if(abs(abs(yDiff/yly) -1.0_RK)> 1.0E-10_RK .and. abs(yDiff/yly)> 1.0E-10_RK) &
            call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","yDiff wrong")
          if(abs(abs(zDiff/zlz) -1.0_RK)> 1.0E-10_RK .and. abs(zDiff/zlz)> 1.0E-10_RK) &
            call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","zDiff wrong")                                   
          nGhost_Sml_Tmp = nGhost_Sml_Tmp + 1
                    
          kid = ids_e + nBase
          GhostP_Sml_PosR(kid) = real4(xSml+xDiff,ySml+yDiff,zSml+zDiff, 0.5_RK*dTmp)
        enddo
      ENDDO
    endif
    if(nGhost_Sml_Tmp /= nGhost_Sml) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","nGhost_Sml wrong")
    End Block
    close(unit=iUnit, iostat=ierr)
    if(allocated(GhostP_Big_nSmall)) deallocate(GhostP_Big_nSmall)
    if(allocated(GhostP_Big_id_global_sort)) deallocate(GhostP_Big_id_global_sort)
    if(allocated(GhostP_Big_sort_id)) deallocate(GhostP_Big_sort_id)
        
    nIBP= 0
    nIbpTmp = 0
    Block
    real(RK):: rTmp2
    type(real4)::PosR_Tmp
    logical :: Is_Within_Other_Small
    
    DO ind = 1, nlocal_Big
      ids_b= GPrtcl_Big_ind_Sml(1, ind)
      ids_e= GPrtcl_Big_ind_Sml(2, ind)
      do i = ids_b, ids_e
        Cposition = GPrtcl_Sml_PosR(i)
        rTmp = GPrtcl_Sml_PosR(i)%w -RetractionRatio_1*dxyz
        nPartition= int(4.0_RK*Pi*(rTmp/dxyz)**2, 4)
        if(nPartition < 1) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","nPartition wrong -1")
        
        allocate(LagrangePoints(nPartition, 3))
        call eq_sphere(LagrangePoints)
        IbpVolRatio = Pi*dxyz*(4.0_RK*rTmp*rTmp+ dxyz*dxyz/3.0_RK)/real(nPartition,RK)*(rdx*rdyUniform*rdz)
        do j=1,nPartition
          PosDiff%x = LagrangePoints(j,1)*rTmp
          PosDiff%y = LagrangePoints(j,2)*rTmp
          PosDiff%z = LagrangePoints(j,3)*rTmp      
          LagrangeP = Cposition+ PosDiff
          
          Is_Within_Other_Small = .false.
          do k=ids_b, ids_e
            if(i==k) cycle
            PosR_Tmp  = GPrtcl_Sml_PosR(k)
            PosDiff%x = PosR_Tmp%x - LagrangeP%x
            PosDiff%y = PosR_Tmp%y - LagrangeP%y
            PosDiff%z = PosR_Tmp%z - LagrangeP%z
            rTmp2 = PosR_Tmp%w -RetractionRatio_2*dxyz
            if(PosDiff%x*PosDiff%x +PosDiff%y*PosDiff%y +PosDiff%z*PosDiff%z < rTmp2*rTmp2) then
              Is_Within_Other_Small = .true.
              exit
            endif
          enddo
          if(Is_Within_Other_Small) cycle
          nIbpTmp = nIbpTmp + 1
          
          ic= floor(LagrangeP%x*rdx)+1
          jc= floor(LagrangeP%y*rdyUniform)+1
          kc= floor(LagrangeP%z*rdz)+1        
          if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 .or. &
             (jc-y1start(2))*(jc-y1end(2))>0) cycle          
          nIBP= nIBP+1
          
          if(nIBP> mIBP) call Reallocate_IbpVar()
          IBP_idlocal(nIBP)= ind   !!! note here
          IBP_Pos(nIBP)= LagrangeP
          IBP_Vel(nIBP)= zero_r3
          IBP_VolRatio(nIBP) = IbpVolRatio
#define   clc_Point_indxyz_IBPMove
#include "pIBM_clc_Point_indxyz_inc.f90"
#undef    clc_Point_indxyz_IBPMove
        enddo
        deallocate(LagrangePoints)
      enddo
    ENDDO
    
    DO ind = 1,nGhost_Big
      ids_b= GhostP_Big_ind_Sml(1, ind)
      ids_e= GhostP_Big_ind_Sml(2, ind)
      do i = ids_b, ids_e
        Cposition = GhostP_Sml_PosR(i)
        rTmp = GhostP_Sml_PosR(i)%w -RetractionRatio_1*dxyz
        nPartition= int(4.0_RK*Pi*(rTmp/dxyz)**2, 4)
        if(nPartition < 1) call MainLog%CheckForError(ErrT_Abort,"PrepareIBM_interp","nPartition wrong -2")
        
        allocate(LagrangePoints(nPartition, 3))
        call eq_sphere(LagrangePoints)      
        IbpVolRatio = Pi*dxyz*(4.0_RK*rTmp*rTmp+ dxyz*dxyz/3.0_RK)/real(nPartition,RK)*(rdx*rdyUniform*rdz)
        do j=1,nPartition
          PosDiff%x = LagrangePoints(j,1)*rTmp
          PosDiff%y = LagrangePoints(j,2)*rTmp
          PosDiff%z = LagrangePoints(j,3)*rTmp
          LagrangeP = Cposition+ PosDiff
          
          Is_Within_Other_Small = .false.
          do k=ids_b, ids_e
            if(i==k) cycle
            PosR_Tmp  = GhostP_Sml_PosR(k)
            PosDiff%x = PosR_Tmp%x - LagrangeP%x
            PosDiff%y = PosR_Tmp%y - LagrangeP%y
            PosDiff%z = PosR_Tmp%z - LagrangeP%z
            rTmp2 = PosR_Tmp%w -RetractionRatio_2*dxyz
            if(PosDiff%x*PosDiff%x +PosDiff%y*PosDiff%y +PosDiff%z*PosDiff%z <rTmp2*rTmp2) then
              Is_Within_Other_Small = .true.
              exit
            endif
          enddo
          if(Is_Within_Other_Small) cycle
          
          ic= floor(LagrangeP%x*rdx)+1
          jc= floor(LagrangeP%y*rdyUniform)+1
          kc= floor(LagrangeP%z*rdz)+1        
          if((ic-y1start(1))*(ic-y1end(1))>0 .or. (kc-y1start(3))*(kc-y1end(3))>0 .or. &
             (jc-y1start(2))*(jc-y1end(2))>0) cycle          
          nIBP= nIBP+1
          if(nIBP> mIBP) call Reallocate_IbpVar()
          IBP_idlocal(nIBP)= ind+nlocal_Big   !!! note here
          IBP_Pos(nIBP)= LagrangeP
          IBP_Vel(nIBP)= zero_r3
          IBP_VolRatio(nIBP) = IbpVolRatio
#define   clc_Point_indxyz_IBPMove
#include "pIBM_clc_Point_indxyz_inc.f90"
#undef    clc_Point_indxyz_IBPMove                                 
        enddo
        deallocate(LagrangePoints)
      enddo
    ENDDO
    End Block

    call MPI_REDUCE(nIBP,nIbpSum1,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(nIbpTmp,nIbpSum2,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(nrank==0) print*, 'Total number of IBP is',nIbpSum1,nIbpSum2
    
    do i=1,nIBP
      IBP_Force(i) =zero_r3
    enddo
    IF(IBMForceScheme==1) call PrepareIBM_ForceAmplify(VolForce_x,VolForce_y,VolForce_z)
    
  end subroutine PrepareIBM_interp
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Init_Prctl_MPI_TYPE
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Init_Prctl_MPI_TYPE(); implicit none
    integer::ierr
    integer,dimension(4)::disp,blocklen,blocktype
  
    ! integer
    int_type = MPI_INTEGER
    call MPI_TYPE_SIZE(int_type,   int_byte,   ierr)

    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,  real_byte,  ierr)

    ! real3 type
    blocklen(1:3)=1
    blocktype(1:3)=real_type
    disp(1)=0
    disp(2)=disp(1)+real_byte
    disp(3)=disp(2)+real_byte
    call MPI_TYPE_STRUCT(3,blocklen(1:3),disp(1:3),blocktype(1:3),real3_type,ierr)
    call MPI_TYPE_COMMIT(real3_type,ierr)
    call MPI_TYPE_SIZE(real3_type, real3_byte, ierr)

    ! real4 type
    blocklen(1:4)=1
    blocktype(1:4)=real_type
    disp(1)=0
    disp(2)=disp(1)+real_byte
    disp(3)=disp(2)+real_byte
    disp(4)=disp(3)+real_byte
    call MPI_TYPE_STRUCT(4,blocklen(1:4),disp(1:4),blocktype(1:4),real4_type,ierr)
    call MPI_TYPE_COMMIT(real4_type,ierr)
    call MPI_TYPE_SIZE(real4_type, real4_byte, ierr)
  end subroutine Init_Prctl_MPI_TYPE
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_Init_Visu
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Prtcl_Init_Visu(); implicit none
    integer(kind=MPI_OFFSET_KIND)::disp
    integer:: iUnit, ierr, nfld, ifld, ptime, dims, iprec
    character(len=:),allocatable::XdmfFile
    
    if(nrank/=0) return
    XdmfFile = strip(ResultsDir_)//"PartVisuFor"//strip(RunName_)//".xmf"
    open(newunit=iUnit, file=XdmfFile,status='replace',form='formatted',IOSTAT=ierr)
    if(ierr /= 0) call MainLog%CheckForError(ErrT_Abort,"Prtcl_Init_Visu","Cannot open file:  "//XdmfFile)
    write(iUnit,'(A)') '<?xml version="1.0" ?>'
    write(iUnit,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
    write(iUnit,'(A)') '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
    write(iUnit,'(A)') '<Domain>'
    
    nfld = (ilast - ifirst +1)/SaveVisu + 1
    write(iUnit,'(A)')'  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
    write(iUnit,'(A)')'    <Time TimeType="List">'
    write(iUnit,'(A,I6,A)')'      <DataItem Format="XML" NumberType="Int" Dimensions="',nfld,'">' 
    write(iUnit,'(A)',advance='no')'        '
    do ifld =1, nfld
      if(mod(ifld,10)==0) then
        write(iUnit,'(I10)') ((ifld-1)*SaveVisu + ifirst-1)
        if(ifld < nfld) write(iUnit,'(A)',advance='no') '        '
      else
        write(iUnit,'(I9)',advance='no') ((ifld-1)*SaveVisu + ifirst-1)
      endif
    enddo    
    if(mod(nfld,10) /=0) write(iUnit,*)' '
    write(iUnit,'(A)') '      </DataItem>'
    write(iUnit,'(A)') '    </Time>'

    do ifld = 1,nfld
      disp = 0_MPI_OFFSET_KIND
      ptime = (ifld-1)*SaveVisu + ifirst-1
      XdmfFile = "PartVisuFor"//strip(RunName_)
      dims=3; iprec=RK
      write(iUnit,'(A,I10.10,A)') '    <Grid Name="T',ptime,'" GridType="Uniform">'
      write(iUnit,'(A,I9,A)') '      <Topology TopologyType="Polyvertex" NodesPerElement="',numPrtcl_Sml,'"/>'
      write(iUnit,'(A)') '      <Geometry GeometryType="'//"XYZ"//'">'
      write(iUnit,'(A,I1,A,I2,I9,A,I15,A)')  '        <DataItem Format="Binary"' // &
          ' DataType="Float" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,numPrtcl_Sml,'" Seek="',disp,'">'
      disp = disp+numPrtcl_Sml*dims*iprec
      write(iUnit,'(A,I10.10)') '          ' // XdmfFile, ptime
      write(iUnit,'(A)') '        </DataItem>'
      write(iUnit,'(A)') '      </Geometry>'

      ! id
      dims=1; iprec=IK
      call Write_XDMF_One(iUnit,dims,iprec,numPrtcl_Sml,ptime,XdmfFile,"ID","Scalar","Int",disp)
      
      ! diameter
      dims=1; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,numPrtcl_Sml,ptime,XdmfFile,"Diameter","Scalar","Float",disp)

      ! FpForce      
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,numPrtcl_Sml,ptime,XdmfFile,"FpForce","Vector","Float",disp)

      ! FpTorque    
      dims=3; iprec=RK
      call Write_XDMF_One(iUnit,dims,iprec,numPrtcl_Sml,ptime,XdmfFile,"FpToruqe","Vector","Float",disp)
      
      write(iUnit,'(A)')'    </Grid>'
    enddo
                    
    ! XDMF/XMF Tail
    open(newunit=iUnit, file=XdmfFile,status='old',position='append',form='formatted',IOSTAT=ierr)
    write(iUnit,'(A)') '  </Grid>'
    write(iUnit,'(A)') '</Domain>'
    write(iUnit,'(A)') '</Xdmf>'
    close(iUnit, IOSTAT=ierr)
  end subroutine Prtcl_Init_Visu

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Write_XDMF_One
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Write_XDMF_One(iUnit,dims,iprec,np,itime,XdmfFile,chName,chAttribute,chDataType,disp); implicit none
    integer,intent(in)::iUnit,dims,iprec,np,itime
    character(len=*),intent(in)::XdmfFile,chName,chAttribute,chDataType
    integer(kind=MPI_OFFSET_KIND),intent(inout)::disp

    write(iUnit,'(A)') '      <Attribute Type="'//strip(chAttribute)//'" Center="Node" Name="'//strip(chName)//'">'
    write(iUnit,'(3A,I1,A,I2,I9,A,I15,A)')  '        <DataItem Format="Binary"' // &
          ' DataType="',chDataType,'" Precision="',iprec,'" Endian="Native"' // &
          ' Dimensions="',dims,np,'" Seek="',disp,'">'
    disp = disp+np*dims*iprec
    write(iUnit,'(A,I10.10)') '          ' // XdmfFile, itime
    write(iUnit,'(A)') '        </DataItem>'
    write(iUnit,'(A)') '      </Attribute>'
  end subroutine Write_XDMF_One
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Prtcl_Dump_Visu
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Prtcl_Dump_Visu(ptime); implicit none
    integer,intent(in)::ptime
    
    ! locals
    type(part_io_size_vec)::pvsize
    integer(kind=MPI_OFFSET_KIND)::disp
    character(len=:),allocatable::chFile
    integer,allocatable,dimension(:)::intVec
    real(RK),allocatable,dimension(:)::realVec
    type(real3),allocatable,dimension(:)::real3Vec
    integer:: bgn_ind,ierr,fh, color,key,Prtcl_WORLD,i,iBig,k,ids_b,ids_e
    
    ! update the bgn_ind
    bgn_ind= clc_bgn_ind(nlocal_Sml)

    ! Create and empty file
    chFile = strip(ResultsDir_)//"PartVisuFor"//strip(RunName_)//int2str(ptime,10)
    if(nrank==0) then
      open(newunit=fh,file=chFile,status='replace',form='unformatted',access='stream',action='write',IOSTAT=ierr)
      close(fh,IOSTAT=ierr)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
    ! create the Prtcl_GROUP
    color = 1; key=nrank
    if(nlocal_Sml <= 0) color=2
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,Prtcl_WORLD,ierr)
    if(color==2) return

    ! begin to dump
    call MPI_FILE_OPEN(Prtcl_WORLD, chFile, MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
    call MPI_BARRIER(Prtcl_WORLD,ierr)
    disp = 0_MPI_OFFSET_KIND
    pvsize%sizes(1)    = numPrtcl_Sml
    pvsize%subsizes(1) = nlocal_Sml
    pvsize%starts(1)   = bgn_ind

    ! position
    allocate(real3Vec(nlocal_Sml))
    do i=1,nlocal_Sml
      real3Vec(i)=GPrtcl_Sml_PosR(i)
    enddo
    call Prtcl_dump(fh,disp, real3Vec,  pvsize)
    deallocate(real3Vec)
    
    ! id
    k=0
    allocate(intVec(nlocal_Sml))
    do iBig = 1, nlocal_Big
      ids_b= GPrtcl_Big_ind_Sml(1, iBig)
      ids_e= GPrtcl_Big_ind_Sml(2, iBig) 
      do i = ids_b, ids_e
        k =k + 1
        intVec(k) = GPrtcl_Big_id(iBig)
      enddo
    enddo
    call Prtcl_dump(fh,disp, intVec,  pvsize)
    deallocate(intVec)
    
    ! diameter
    allocate(realVec(nlocal_Sml))
    do i=1,nlocal_Sml
      realVec(i)= 2.0_RK*GPrtcl_Sml_PosR(i)%w
    enddo   
    call Prtcl_dump(fh,disp, realVec,  pvsize)
    deallocate(realVec)

    ! FpForce
    k=0
    allocate(real3Vec(nlocal_Sml))
    do iBig = 1, nlocal_Big
      ids_b= GPrtcl_Big_ind_Sml(1, iBig)
      ids_e= GPrtcl_Big_ind_Sml(2, iBig) 
      do i = ids_b, ids_e
        k =k + 1
        real3Vec(k) = GPrtcl_Big_FpForce(iBig)
      enddo
    enddo
    call Prtcl_dump(fh,disp, real3Vec,  pvsize)

    ! FpTorque
    k=0
    do iBig = 1, nlocal_Big
      ids_b= GPrtcl_Big_ind_Sml(1, iBig)
      ids_e= GPrtcl_Big_ind_Sml(2, iBig)
      do i = ids_b, ids_e
        k =k + 1
        real3Vec(k) = GPrtcl_Big_FpTorque(iBig)
      enddo
    enddo
    call Prtcl_dump(fh,disp, real3Vec,  pvsize)    
    deallocate(real3Vec)
        
    call MPI_FILE_CLOSE(fh, ierr) 
    call MPI_COMM_FREE( Prtcl_WORLD, ierr)    
  end subroutine Prtcl_Dump_Visu
    
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Init_IBM
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Init_IBM(chFile); implicit none
    character(len=*),intent(in)::chFile
    
    ! locals
    character(512)::SpheresInfoFile
    integer:: ierr,ierrTmp,iUnit, deltaFunction_Scheme
    Namelist/IBMForce_Options/ IBMForceScheme, ForceAmplifyCoe, nforcingextra, &
      deltaFunction_Scheme, RetractionRatio_1, RetractionRatio_2, SpheresInfoFile

    call Init_Prctl_MPI_TYPE()
    open(newunit=iUnit, file=chFile, status='old', form='formatted', IOSTAT=ierr)
    if(ierr /= 0) call MainLog%CheckForError(ErrT_Abort,"Init_IBM","Cannot open file:"//strip(chFile))
    read(iUnit, nml=IBMForce_Options)
    if(nrank==0) write(MainLog%iUnit, nml=IBMForce_Options)
    close(unit=iUnit, IOSTAT=ierr)
    SpheresInfoFile_ = strip(SpheresInfoFile)
    
    ! check integer(kind=2) is enough or not.
    if(nrank==0 .and. (nxc>huge(int(0,2))-2 .or. nyc>huge(int(0,2))-2 .or. nzc>huge(int(0,2))-2)) then
      call MainLog%CheckForError(ErrT_Abort,"Init_IBM","kind=2 is not enough for indxyz and IBPFix_indxyz")
    endif

#ifndef IBMDistributeLinear    
    if(deltaFunction_Scheme ==0) then
      deltaFunction => deltaFunction_Roma
    else
      deltaFunction => deltaFunction_Yang
    endif     
#endif

    select case(IBMForceScheme)
    case(0)  ! Kempe(2012,JCP), Gsell(2021,JCP)
      AdditionalForceIBM => AdditionalForceIBM_0
    case(1)  ! Zhao(2021,JCP), Cheylan(2023,JCP)
      AdditionalForceIBM => AdditionalForceIBM_1
    end select
        
    xstCoord= real(y1start(1)-1, kind=RK)*dx
    xedCoord= real(y1end(1),     kind=RK)*dx
    ystCoord= 0.0_RK
    yedCoord= yly
    zstCoord= real(y1start(3)-1, kind=RK)*dz
    zedCoord= real(y1end(3),     kind=RK)*dz
    SMALL= 1.0E-9_RK*dx
    
    !=========================== Immersed boundary points parts ===========================!
    CellRadius= 0.5_RK*sqrt(dx*dx+dyUniform*dyUniform+dz*dz)
    CellVolumeIBM= dx* dyUniform* dz
    dxhalf=dx*0.5_RK; dyhalf=dyUniform*0.5_RK; dzhalf=dz*0.5_RK
    mIBP = 1000

    allocate(IBP_VolRatio(mIBP),Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(IBP_idlocal(mIBP), Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(IBP_indxyz(6,mIBP),Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(IBP_Pos(mIBP),     Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(IBP_Vel(mIBP),     Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    allocate(IBP_Force(mIBP),   Stat=ierrTmp); ierr =ierr +abs(ierrTmp)
    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"Init_IBM","Allocation failed 2")
    IBP_idlocal=0; IBP_indxyz=0
    IBP_Pos=zero_r3; IBP_Vel=zero_r3; IBP_Force=zero_r3
  end subroutine Init_IBM

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PrepareIBM_IbpForce
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PrepareIBM_IbpForce( ); implicit none
    integer::i
    do i=1,nIBP
      IBP_Force(i) =zero_r3
    enddo  
  end subroutine PrepareIBM_IbpForce
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Update_FluidIndicator
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Update_FluidIndicator(FluidIndicator); implicit none
    character,dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::FluidIndicator

    ! locals
    type(real3)::Cposition,CenterDiff
    real(RK)::Radius,CenterDistY,CenterDistZ,CenterDist
    integer::i,j,k,pid,DomainIx1,DomainIx2,DomainIy1,DomainIy2,DomainIz1,DomainIz2
    
    FluidIndicator= 'F'

    DO pid=1,nlocal_Sml
      Cposition= GPrtcl_Sml_PosR(pid)
      Radius   = GPrtcl_Sml_PosR(pid)%w
      DomainIx1= floor((Cposition%x-Radius)*rdx )+1;        DomainIx1=max(y1start(1),DomainIx1)
      DomainIx2= ceiling((Cposition%x+Radius)*rdx );        DomainIx2=min(y1end(1),  DomainIx2)
      DomainIy1= floor((Cposition%y-Radius)*rdyUniform )+1; DomainIy1=max(y1start(2),DomainIy1)
      DomainIy2= ceiling((Cposition%y+Radius)*rdyUniform ); DomainIy2=min(y1end(2),  DomainIy2)
      DomainIz1= floor((Cposition%z-Radius)*rdz )+1;        DomainIz1=max(y1start(3),DomainIz1)
      DomainIz2= ceiling((Cposition%z+Radius)*rdz );        DomainIz2=min(y1end(3),  DomainIz2)
      do k=DomainIz1,DomainIz2
        CenterDiff%z= real(k,RK)*dz-dzhalf- Cposition%z
        CenterDistZ=CenterDiff%z*CenterDiff%z
        do j=DomainIy1,DomainIy2
          CenterDiff%y= real(j,RK)*dyUniform-dyhalf- Cposition%y
          CenterDistY=CenterDistZ+CenterDiff%y*CenterDiff%y
          do i=DomainIx1,DomainIx2
            CenterDiff%x= real(i,RK)*dx-dxhalf- Cposition%x
            CenterDist = CenterDiff%x*CenterDiff%x +CenterDistY
            if(CenterDist > Radius*Radius) cycle
            FluidIndicator(i,j,k)='P'
          enddo
        enddo
      enddo
    ENDDO
    DO pid=1,nGhost_Sml
      Cposition= GhostP_Sml_PosR(pid)
      Radius   = GhostP_Sml_PosR(pid)%w
      DomainIx1= floor((Cposition%x-Radius)*rdx )+1;        DomainIx1=max(y1start(1),DomainIx1)
      DomainIx2= ceiling((Cposition%x+Radius)*rdx );        DomainIx2=min(y1end(1),  DomainIx2)
      DomainIy1= floor((Cposition%y-Radius)*rdyUniform )+1; DomainIy1=max(y1start(2),DomainIy1)
      DomainIy2= ceiling((Cposition%y+Radius)*rdyUniform ); DomainIy2=min(y1end(2),  DomainIy2)
      DomainIz1= floor((Cposition%z-Radius)*rdz )+1;        DomainIz1=max(y1start(3),DomainIz1)
      DomainIz2= ceiling((Cposition%z+Radius)*rdz );        DomainIz2=min(y1end(3),  DomainIz2)
      do k=DomainIz1,DomainIz2
        CenterDiff%z= real(k,RK)*dz-dzhalf- Cposition%z
        CenterDistZ=CenterDiff%z*CenterDiff%z
        do j=DomainIy1,DomainIy2
          CenterDiff%y= real(j,RK)*dyUniform-dyhalf- Cposition%y
          CenterDistY=CenterDistZ+CenterDiff%y*CenterDiff%y
          do i=DomainIx1,DomainIx2
            CenterDiff%x= real(i,RK)*dx-dxhalf- Cposition%x
            CenterDist = CenterDiff%x*CenterDiff%x +CenterDistY
            if(CenterDist > Radius*Radius) cycle
            FluidIndicator(i,j,k)='P'
          enddo
        enddo
      enddo
    ENDDO
  end subroutine Update_FluidIndicator

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! IntegrateFluidPrtclForce
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine IntegrateFluidPrtclForce(); implicit none

    ! locals
    real(RK)::DenAlpha
    type(real3)::FpVolume,Cposition,PosDiff
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer::i,pid,idlocal,idGhost,nsend,nrecv,nsend2,nrecv2,ierr
    integer::ProcNgh(4),GIBM_size_backward,request(4),SRstatus(MPI_STATUS_SIZE)

    ! (idlocal 1) +(FpForce 3) +(FpTorque 3) = 7
    GIBM_size_backward= 7
    do i=1,4
      if(myProcNghBC(y_pencil,i)<0) then
        ProcNgh(i)=MPI_PROC_NULL
      else
        ProcNgh(i)=myProcNghBC(y_pencil,i)
      endif
    enddo
    
    GPrtcl_Big_FpForce = zero_r3
    GPrtcl_Big_FpTorque= zero_r3
    GhostP_Big_FpForce = zero_r3
    GhostP_Big_FpTorque= zero_r3
    DO pid=1,nIBP
      idlocal  = IBP_idlocal(pid)
      if(idlocal> nlocal_Big) then
        idGhost= idlocal- nlocal_Big
        FpVolume= IBP_Force(pid)*IBP_VolRatio(pid)*CellVolumeIBM
        GhostP_Big_FpForce(idGhost)=  GhostP_Big_FpForce(idGhost)+ FpVolume

        Cposition = GhostP_Big_PosR(idGhost)
        PosDiff= IBP_Pos(pid)- Cposition
        GhostP_Big_FpTorque(idGhost)= GhostP_Big_FpTorque(idGhost)+ (PosDiff .cross. FpVolume)
      else
        FpVolume= IBP_Force(pid)*IBP_VolRatio(pid)*CellVolumeIBM
        GPrtcl_Big_FpForce(idlocal) = GPrtcl_Big_FpForce(idlocal)+ FpVolume

        Cposition = GPrtcl_Big_PosR(idlocal)
        PosDiff= IBP_Pos(pid)- Cposition
        GPrtcl_Big_FpTorque(idlocal)= GPrtcl_Big_FpTorque(idlocal)+ (PosDiff .cross. FpVolume)
      endif
    ENDDO

#ifdef TestSumFpForce
    block
    real(RK)::Fp1(6),Fp_sum(6)
    Fp1 = 0.0_RK
    do pid=1,nGhost_Big
      Fp1(1) = Fp1(1) + GhostP_Big_FpForce(pid)%x
      Fp1(2) = Fp1(2) + GhostP_Big_FpForce(pid)%y
      Fp1(3) = Fp1(3) + GhostP_Big_FpForce(pid)%z    
      Fp1(4) = Fp1(4) + GhostP_Big_FpTorque(pid)%x
      Fp1(5) = Fp1(5) + GhostP_Big_FpTorque(pid)%y
      Fp1(6) = Fp1(6) + GhostP_Big_FpTorque(pid)%z
    enddo
    do pid=1,nlocal_Big
      Fp1(1) = Fp1(1) + GPrtcl_Big_FpForce(pid)%x
      Fp1(2) = Fp1(2) + GPrtcl_Big_FpForce(pid)%y
      Fp1(3) = Fp1(3) + GPrtcl_Big_FpForce(pid)%z  
      Fp1(4) = Fp1(4) + GPrtcl_Big_FpTorque(pid)%x
      Fp1(5) = Fp1(5) + GPrtcl_Big_FpTorque(pid)%y
      Fp1(6) = Fp1(6) + GPrtcl_Big_FpTorque(pid)%z
    enddo
    call MPI_REDUCE(Fp1, Fp_sum, 6, real_type, MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(nrank==0) print*, '============Fp_sum-01',Fp_sum
    end block
#endif

    ! Step5: send to zp_axis, and receive from zm_dir
    nsend=0; nrecv=0
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
      DO i=1, nGhost_Big
        if(GhostP_Big_Direction(i)/= THIS_PROC_ZM) cycle
        idlocal= GhostP_Big_id(i)
        if(idlocal> nlocal_Big) then
          idGhost= idlocal- nlocal_Big
          GhostP_Big_FpForce(idGhost) = GhostP_Big_FpForce(idGhost) +GhostP_Big_FpForce(i)
          GhostP_Big_FpTorque(idGhost)= GhostP_Big_FpTorque(idGhost)+GhostP_Big_FpTorque(i)
        else
          GPrtcl_Big_FpForce(idlocal) = GPrtcl_Big_FpForce(idlocal) +GhostP_Big_FpForce(i)
          GPrtcl_Big_FpTorque(idlocal)= GPrtcl_Big_FpTorque(idlocal)+GhostP_Big_FpTorque(i)
        endif
      ENDDO
    ELSEIF(ProcNgh(1) /= MPI_PROC_NULL) THEN
      DO i=1, nGhost_Big
        if(GhostP_Big_Direction(i) /= zm_dir) cycle
        nsend=nsend+1
        if(nsend > msend) call reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, MPI_INTEGER, ProcNgh(1), 1, &
                      nrecv, 1, MPI_INTEGER, ProcNgh(2), 1, MPI_COMM_WORLD,SRstatus,ierr)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,ProcNgh(2),2,MPI_COMM_WORLD,request(1),ierr)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,ProcNgh(1),2,MPI_COMM_WORLD,ierr)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(1),SRstatus,ierr)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send)) deallocate(buf_send) 
    if(allocated(buf_recv)) deallocate(buf_recv)
        
    ! Step4: send to zm_axis, and receive from zp_dir
    nsend=0; nrecv=0
    IF(ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself. 
      DO i=1, nGhost_Big
        if(GhostP_Big_Direction(i)/= THIS_PROC_ZP) cycle
        idlocal= GhostP_Big_id(i)
        if(idlocal> nlocal_Big) then
          idGhost= idlocal- nlocal_Big
          GhostP_Big_FpForce(idGhost) = GhostP_Big_FpForce(idGhost) +GhostP_Big_FpForce(i)
          GhostP_Big_FpTorque(idGhost)= GhostP_Big_FpTorque(idGhost)+GhostP_Big_FpTorque(i)
        else
          GPrtcl_Big_FpForce(idlocal) = GPrtcl_Big_FpForce(idlocal) +GhostP_Big_FpForce(i)
          GPrtcl_Big_FpTorque(idlocal)= GPrtcl_Big_FpTorque(idlocal)+GhostP_Big_FpTorque(i)
        endif
      ENDDO
    ELSEIF(ProcNgh(2) /= MPI_PROC_NULL) THEN
      DO i=1, nGhost_Big
        if(GhostP_Big_Direction(i) /= zp_dir) cycle
        nsend=nsend+1     
        if(nsend > msend) call reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, MPI_INTEGER, ProcNgh(2), 3, &
                      nrecv, 1, MPI_INTEGER, ProcNgh(1), 3, MPI_COMM_WORLD,SRstatus,ierr)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,ProcNgh(1),4,MPI_COMM_WORLD,request(2),ierr)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,ProcNgh(2),4,MPI_COMM_WORLD,ierr)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(2),SRstatus,ierr)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send)) deallocate(buf_send) 
    if(allocated(buf_recv)) deallocate(buf_recv)

    ! Step3: send to xp_axis, and receive from xm_dir
    nsend=0; nrecv=0
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.    
      DO i=1, nGhost_Big
        if(GhostP_Big_Direction(i)/= THIS_PROC_XM) cycle
        idlocal= GhostP_Big_id(i)
        if(idlocal> nlocal_Big) then
          idGhost= idlocal- nlocal_Big
          GhostP_Big_FpForce(idGhost) = GhostP_Big_FpForce(idGhost) +GhostP_Big_FpForce(i)
          GhostP_Big_FpTorque(idGhost)= GhostP_Big_FpTorque(idGhost)+GhostP_Big_FpTorque(i)
        else
          GPrtcl_Big_FpForce(idlocal) = GPrtcl_Big_FpForce(idlocal) +GhostP_Big_FpForce(i)
          GPrtcl_Big_FpTorque(idlocal)= GPrtcl_Big_FpTorque(idlocal)+GhostP_Big_FpTorque(i)
        endif
      ENDDO
    ELSEIF(ProcNgh(3) /= MPI_PROC_NULL) THEN
      DO i=1, nGhost_Big
        if(GhostP_Big_Direction(i) /= xm_dir) cycle
        nsend=nsend+1     
        if(nsend > msend) call reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, MPI_INTEGER, ProcNgh(3), 5, &
                      nrecv, 1, MPI_INTEGER, ProcNgh(4), 5, MPI_COMM_WORLD,SRstatus,ierr)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,ProcNgh(4),6,MPI_COMM_WORLD,request(3),ierr)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,ProcNgh(3),6,MPI_COMM_WORLD,ierr)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(3),SRstatus,ierr)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send)) deallocate(buf_send) 
    if(allocated(buf_recv)) deallocate(buf_recv)

    ! Step2: send to xm_axis, and receive from xp_dir
    nsend=0; nrecv=0
    IF(ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.    
      DO i=1, nGhost_Big
        if(GhostP_Big_Direction(i)/= THIS_PROC_XP) cycle
        idlocal= GhostP_Big_id(i)
        if(idlocal> nlocal_Big) then
          idGhost= idlocal- nlocal_Big
          GhostP_Big_FpForce(idGhost) = GhostP_Big_FpForce(idGhost) +GhostP_Big_FpForce(i)
          GhostP_Big_FpTorque(idGhost)= GhostP_Big_FpTorque(idGhost)+GhostP_Big_FpTorque(i)        
        else
          GPrtcl_Big_FpForce(idlocal) = GPrtcl_Big_FpForce(idlocal) +GhostP_Big_FpForce(i)
          GPrtcl_Big_FpTorque(idlocal)= GPrtcl_Big_FpTorque(idlocal)+GhostP_Big_FpTorque(i)
        endif
      ENDDO
    ELSEIF(ProcNgh(4) /= MPI_PROC_NULL) THEN
      DO i=1, nGhost_Big
        if(GhostP_Big_Direction(i) /= xp_dir) cycle
        nsend=nsend+1
        if(nsend > msend) call reallocate_sendlist(nsend)
        sendlist(nsend)=i   
      ENDDO
    ENDIF
    call MPI_SENDRECV(nsend, 1, MPI_INTEGER, ProcNgh(4), 7, &
                      nrecv, 1, MPI_INTEGER, ProcNgh(3), 7, MPI_COMM_WORLD,SRstatus,ierr)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_backward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,ProcNgh(3),8,MPI_COMM_WORLD,request(4),ierr)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_backward
      allocate(buf_send(nsend2))
      call pack_IBM_Fpforce(buf_send,nsend)
      call MPI_SEND(buf_send,nsend2,real_type,ProcNgh(4),8,MPI_COMM_WORLD,ierr)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(4),SRstatus,ierr)
      call unpack_IBM_Fpforce(buf_recv,nrecv)
    endif
    if(allocated(buf_send)) deallocate(buf_send)
    if(allocated(buf_recv)) deallocate(buf_recv)

    ! Step1: Handle y-dir
    IF(BcOption(yp_dir)==BC_PERIOD) THEN
      DO i=1, nGhost_Big
        if(GhostP_Big_Direction(i) /= THIS_PROC_YP .and. GhostP_Big_Direction(i) /= THIS_PROC_YM) cycle
        idlocal= GhostP_Big_id(i)
        GPrtcl_Big_FpForce(idlocal) = GPrtcl_Big_FpForce(idlocal) +GhostP_Big_FpForce(i)
        GPrtcl_Big_FpTorque(idlocal)= GPrtcl_Big_FpTorque(idlocal)+GhostP_Big_FpTorque(i)
      ENDDO      
    ENDIF  

#ifdef TestSumFpForce
    block
    real(RK)::Fp1(6),Fp_sum(6)
    Fp1 = 0.0_RK
    do pid=1,nlocal_Big
      Fp1(1) = Fp1(1) + GPrtcl_Big_FpForce(pid)%x
      Fp1(2) = Fp1(2) + GPrtcl_Big_FpForce(pid)%y
      Fp1(3) = Fp1(3) + GPrtcl_Big_FpForce(pid)%z  
      Fp1(4) = Fp1(4) + GPrtcl_Big_FpTorque(pid)%x
      Fp1(5) = Fp1(5) + GPrtcl_Big_FpTorque(pid)%y
      Fp1(6) = Fp1(6) + GPrtcl_Big_FpTorque(pid)%z
    enddo
    call MPI_REDUCE(Fp1, Fp_sum, 6, real_type, MPI_SUM,0,MPI_COMM_WORLD,ierr)
    if(nrank==0) print*, '============Fp_sum-02',Fp_sum
    end block
#endif
    
    ! Finally, get the final GPrtcl_Big_FpForce
    DenAlpha=-FluidDensity/pmAlpha
    do i=1, nlocal_Big
      GPrtcl_Big_FpForce(i) = DenAlpha*GPrtcl_Big_FpForce(i)
      GPrtcl_Big_FpTorque(i)= DenAlpha*GPrtcl_Big_FpTorque(i)
    enddo
  end subroutine IntegrateFluidPrtclForce

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! pack_IBM_FpForce
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine pack_IBM_FpForce(buf_send,nsend); implicit none
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsend

    ! locals
    integer::i,pid,m

    m=1
    do i=1,nsend
      pid= sendlist(i)
      buf_send(m)= real(GhostP_Big_id(pid),RK); m=m+1  ! 01
      buf_send(m)= GhostP_Big_FpForce(pid)%x;   m=m+1  ! 02
      buf_send(m)= GhostP_Big_FpForce(pid)%y;   m=m+1  ! 03
      buf_send(m)= GhostP_Big_FpForce(pid)%z;   m=m+1  ! 04
      buf_send(m)= GhostP_Big_FpTorque(pid)%x;  m=m+1  ! 05
      buf_send(m)= GhostP_Big_FpTorque(pid)%y;  m=m+1  ! 06
      buf_send(m)= GhostP_Big_FpTorque(pid)%z;  m=m+1  ! 07
    enddo
  end subroutine pack_IBM_FpForce
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! unpack_IBM_FpForce
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine unpack_IBM_FpForce(buf_recv,nrecv); implicit none
    real(RK),dimension(:),intent(in):: buf_recv
    integer,intent(in)::nrecv

    !locals
    type(real3)::FpForce,FpTorque
    integer::i,pid,m,gid

    m=1
    do i=1,nrecv
      pid= nint(buf_recv(m)); m=m+1  ! 01
      FpForce%x =buf_recv(m); m=m+1  ! 02
      FpForce%y =buf_recv(m); m=m+1  ! 03
      FpForce%z =buf_recv(m); m=m+1  ! 04
      FpTorque%x=buf_recv(m); m=m+1  ! 05
      FpTorque%y=buf_recv(m); m=m+1  ! 06
      FpTorque%z=buf_recv(m); m=m+1  ! 07
      if(pid>nlocal_Big) then
        gid= pid-nlocal_Big
        GhostP_Big_FpForce(gid) = GhostP_Big_FpForce(gid) +FpForce
        GhostP_Big_FpTorque(gid)= GhostP_Big_FpTorque(gid)+FpTorque
      else
        GPrtcl_Big_FpForce(pid) = GPrtcl_Big_FpForce(pid) +FpForce
        GPrtcl_Big_FpTorque(pid)= GPrtcl_Big_FpTorque(pid)+FpTorque
      endif   
    enddo
  end subroutine unpack_IBM_FpForce
      
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PComm_IBM_forward
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PComm_IBM_forward(); implicit none
    ! locals
    real(RK)::px,pxt,py,pyt,pz,pzt
    real(RK),dimension(:),allocatable::buf_send,buf_recv
    integer::nsend,nsend2,nrecv,nrecv2,nsendg,ng,ngp,ngpp,ProcNgh(4),ierrTmp
    integer::i,ierr,request(4),SRstatus(MPI_STATUS_SIZE),mGhostIBM,GIBM_size_forward
        
    ierr = 0
    mSend =100
    mGhostIBM = 100 
    allocate(GhostP_Big_id(mGhostIBM), stat=ierrTmp);         ierr = ierr + ierrTmp
    allocate(GhostP_Big_id_Global(mGhostIBM), stat=ierrTmp);  ierr = ierr + ierrTmp
    allocate(GhostP_Big_nSmall(mGhostIBM), stat=ierrTmp);     ierr = ierr + ierrTmp  
    allocate(GhostP_Big_Direction(mGhostIBM), stat=ierrTmp);  ierr = ierr + ierrTmp
    allocate(GhostP_Big_PosR(mGhostIBM), stat=ierrTmp);       ierr = ierr + ierrTmp    
    allocate(GhostP_Big_FpForce(mGhostIBM), stat=ierrTmp);    ierr = ierr + ierrTmp
    allocate(GhostP_Big_FpTorque(mGhostIBM),stat=ierrTmp);    ierr = ierr + ierrTmp
    allocate(GhostP_Big_ind_Sml(2,mGhostIBM),stat=ierrTmp);   ierr = ierr + ierrTmp
    allocate(sendlist(mSend), stat=ierrTmp);                  ierr = ierr + ierrTmp
    if(ierr /= 0) call MainLog%CheckForError(ErrT_Abort,"PComm_IBM_forward", "Allocation failed-1")
    
    ! (idlocal 1) +(idglobal 1) +(nSmall 1) +(direction 1) +(PosR 4)=8
    GIBM_size_forward=8
    do i=1,4
      if(myProcNghBC(y_pencil,i)<0) then
        ProcNgh(i)=MPI_PROC_NULL
      else
        ProcNgh(i)=myProcNghBC(y_pencil,i)
      endif
    enddo
    
    ! This part is similar to what is done in subroutine PC_Comm_For_Cntct of file "sp_Comm.f90" (only y_pencil)
    ng=0  
    ! Step1: Handle y-dir
    IF(BcOption(yp_dir)==BC_PERIOD) THEN
      do i=1,nlocal_Big
        py = GPrtcl_Big_PosR(i)%y
        pyt= py +GPrtcl_Big_PosR(i)%w
        if(pyt +SMALL >yedCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
          GhostP_Big_id(ng)       = i
          GhostP_Big_id_Global(ng)= GPrtcl_Big_id(i)
          GhostP_Big_nSmall(ng)   = GPrtcl_Big_ind_Sml(2,i)-GPrtcl_Big_ind_Sml(1,i)+1
          GhostP_Big_PosR(ng)     = GPrtcl_Big_PosR(i)
          GhostP_Big_Direction(ng)= THIS_PROC_YP
          GhostP_Big_PosR(ng)%y   = py -yly            
        endif
        pyt= py -GPrtcl_Big_PosR(i)%w
        if(pyt -SMALL <ystCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
          GhostP_Big_id(ng)       = i
          GhostP_Big_id_Global(ng)= GPrtcl_Big_id(i)
          GhostP_Big_nSmall(ng)   = GPrtcl_Big_ind_Sml(2,i)-GPrtcl_Big_ind_Sml(1,i)+1
          GhostP_Big_PosR(ng)     = GPrtcl_Big_PosR(i)
          GhostP_Big_Direction(ng)= THIS_PROC_YM
          GhostP_Big_PosR(ng)%y   = py +yly
        endif        
      enddo
    ENDIF

    ! Step2: send to xp_axis, and receive from xm_dir
    nsend=0; nrecv=0; ngp=ng; ngpp=ng
    IF(ProcNgh(3)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,ngpp    ! consider the previous ghost particle firstly
        px = GhostP_Big_PosR(i)%x
        pxt= px +GhostP_Big_PosR(i)%w
        if(pxt +SMALL >xedCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
          GhostP_Big_id(ng)       = nlocal_Big+ i
          GhostP_Big_id_Global(ng)= GhostP_Big_id_Global(i)
          GhostP_Big_nSmall(ng)   = GhostP_Big_nSmall(i)
          GhostP_Big_PosR(ng)     = GhostP_Big_PosR(i)
          GhostP_Big_Direction(ng)= THIS_PROC_XP
          GhostP_Big_PosR(ng)%x   = px -xlx
        endif
      enddo

      do i=1,nlocal_Big
        px = GPrtcl_Big_PosR(i)%x
        pxt= px +GPrtcl_Big_PosR(i)%w
        if(pxt +SMALL >xedCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
          GhostP_Big_id(ng)       = i
          GhostP_Big_id_Global(ng)= GPrtcl_Big_id(i)
          GhostP_Big_nSmall(ng)   = GPrtcl_Big_ind_Sml(2,i)-GPrtcl_Big_ind_Sml(1,i)+1
          GhostP_Big_PosR(ng)     = GPrtcl_Big_PosR(i)
          GhostP_Big_Direction(ng)= THIS_PROC_XP
          GhostP_Big_PosR(ng)%x   = px -xlx
        endif
      enddo  
    ELSEIF(ProcNgh(3) /= MPI_PROC_NULL) THEN
      do i=1,ngpp    ! consider the previous ghost particle firstly
        pxt= GhostP_Big_PosR(i)%x +GhostP_Big_PosR(i)%w
        if(pxt +SMALL >xedCoord) then
          nsend=nsend+1
          if(nsend > msend) call reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
      
      nsendg=nsend
      do i=1,nlocal_Big
        pxt= GPrtcl_Big_PosR(i)%x +GPrtcl_Big_PosR(i)%w
        if(pxt +SMALL >xedCoord) then
          nsend=nsend+1
          if(nsend > msend) call reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, MPI_INTEGER, ProcNgh(3), 1, &
                      nrecv, 1, MPI_INTEGER, ProcNgh(4), 1, MPI_COMM_WORLD,SRstatus,ierr)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,ProcNgh(4),2,MPI_COMM_WORLD,request(1),ierr)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward(buf_send,nsendg,nsend,xp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,ProcNgh(3),2,MPI_COMM_WORLD,ierr)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(1),SRstatus,ierr)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
      call unpack_IBM_forward(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send)) deallocate(buf_send)
    if(allocated(buf_recv)) deallocate(buf_recv)
    
    ! Step3: send to xm_axis, and receive from xp_dir
    nsend=0; nrecv=0; ngp=ng; !ngpp=ng   
    IF(ProcNgh(4)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,ngpp    ! consider the previous ghost particle firstly
        px = GhostP_Big_PosR(i)%x
        pxt= px -GhostP_Big_PosR(i)%w
        if(pxt -SMALL <xstCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
          GhostP_Big_id(ng)       = nlocal_Big +i
          GhostP_Big_id_Global(ng)= GhostP_Big_id_Global(i)
          GhostP_Big_nSmall(ng)   = GhostP_Big_nSmall(i)
          GhostP_Big_PosR(ng)     = GhostP_Big_PosR(i)
          GhostP_Big_Direction(ng)= THIS_PROC_XM
          GhostP_Big_PosR(ng)%x   = px +xlx
        endif
      enddo

      do i=1,nlocal_Big
        px = GPrtcl_Big_PosR(i)%x
        pxt= px- GPrtcl_Big_PosR(i)%w
        if(pxt -SMALL <xstCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
          GhostP_Big_id(ng)       = i
          GhostP_Big_id_Global(ng)= GPrtcl_Big_id(i)
          GhostP_Big_nSmall(ng)   = GPrtcl_Big_ind_Sml(2,i)-GPrtcl_Big_ind_Sml(1,i)+1
          GhostP_Big_PosR(ng)     = GPrtcl_Big_PosR(i)
          GhostP_Big_Direction(ng)= THIS_PROC_XM
          GhostP_Big_PosR(ng)%x   = px +xlx
        endif
      enddo 
    ELSEIF(ProcNgh(4) /= MPI_PROC_NULL) THEN
      do i=1,ngpp    ! consider the previous ghost particle firstly
        pxt= GhostP_Big_PosR(i)%x -GhostP_Big_PosR(i)%w
        if(pxt -SMALL <xstCoord) then
          nsend=nsend+1
          if(nsend > msend) call reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
      nsendg=nsend
      
      do i=1,nlocal_Big
        pxt= GPrtcl_Big_PosR(i)%x -GPrtcl_Big_PosR(i)%w
        if(pxt -SMALL <xstCoord) then
          nsend=nsend+1
          if(nsend > msend) call reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, MPI_INTEGER, ProcNgh(4), 3, &
                      nrecv, 1, MPI_INTEGER, ProcNgh(3), 3, MPI_COMM_WORLD,SRstatus,ierr)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,ProcNgh(3),4,MPI_COMM_WORLD,request(2),ierr)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward(buf_send,nsendg,nsend,xm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,ProcNgh(4),4,MPI_COMM_WORLD,ierr)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(2),SRstatus,ierr)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
      call unpack_IBM_forward(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step4: send to zp_axis, and receive from zm_dir
    nsend=0; nrecv=0; ngp=ng; ngpp=ng
    IF(ProcNgh(1)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,ngpp    ! consider the previous ghost particle firstly
        pz = GhostP_Big_PosR(i)%z 
        pzt= pz+ GhostP_Big_PosR(i)%w
        if(pzt+SMALL>zedCoord ) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
          GhostP_Big_id(ng)       = nlocal_Big+ i
          GhostP_Big_id_Global(ng)= GhostP_Big_id_Global(i)
          GhostP_Big_nSmall(ng)   = GhostP_Big_nSmall(i)
          GhostP_Big_PosR(ng)     = GhostP_Big_PosR(i)
          GhostP_Big_Direction(ng)= THIS_PROC_ZP
          GhostP_Big_PosR(ng)%z   = pz -zlz
        endif
      enddo

      do i=1,nlocal_Big
        pz = GPrtcl_Big_PosR(i)%z 
        pzt= pz+ GPrtcl_Big_PosR(i)%w
        if(pzt+SMALL>zedCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
          GhostP_Big_id(ng)       = i
          GhostP_Big_id_Global(ng)= GPrtcl_Big_id(i)
          GhostP_Big_nSmall(ng)   = GPrtcl_Big_ind_Sml(2,i)-GPrtcl_Big_ind_Sml(1,i)+1
          GhostP_Big_PosR(ng)     = GPrtcl_Big_PosR(i)
          GhostP_Big_Direction(ng)= THIS_PROC_ZP
          GhostP_Big_PosR(ng)%z   = pz -zlz
        endif
      enddo  
    ELSEIF(ProcNgh(1) /= MPI_PROC_NULL) THEN
      do i=1,ngpp    ! consider the previous ghost particle firstly
        pzt= GhostP_Big_PosR(i)%z+ GhostP_Big_PosR(i)%w
        if(pzt+SMALL>zedCoord) then
          nsend=nsend+1
          if(nsend > msend) call reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
      nsendg=nsend

      do i=1,nlocal_Big
        pzt=GPrtcl_Big_PosR(i)%z+ GPrtcl_Big_PosR(i)%w
        if(pzt+SMALL>zedCoord) then
          nsend=nsend+1
          if(nsend > msend) call reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, MPI_INTEGER, ProcNgh(1), 5, &
                      nrecv, 1, MPI_INTEGER, ProcNgh(2), 5, MPI_COMM_WORLD,SRstatus,ierr)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,ProcNgh(2),6,MPI_COMM_WORLD,request(3),ierr)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward(buf_send,nsendg,nsend,zp_dir)
      call MPI_SEND(buf_send,nsend2,real_type,ProcNgh(1),6,MPI_COMM_WORLD,ierr)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(3),SRstatus,ierr)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
      call unpack_IBM_forward(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send))deallocate(buf_send) 
    if(allocated(buf_recv))deallocate(buf_recv)

    ! Step5: send to zm_axis, and receive from zp_dir
    nsend=0; nrecv=0; ngp=ng; !ngpp=ng
    IF(ProcNgh(2)==nrank) THEN  ! neighbour is nrank itself.
      do i=1,ngpp    ! consider the previous ghost particle firstly
        pz = GhostP_Big_PosR(i)%z
        pzt= pz- GhostP_Big_PosR(i)%w
        if(pzt-SMALL<zstCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
          GhostP_Big_id(ng)       = nlocal_Big+ i
          GhostP_Big_id_Global(ng)= GhostP_Big_id_Global(i)
          GhostP_Big_nSmall(ng)   = GhostP_Big_nSmall(i)
          GhostP_Big_PosR(ng)     = GhostP_Big_PosR(i)
          GhostP_Big_Direction(ng)= THIS_PROC_ZM
          GhostP_Big_PosR(ng)%z   = pz +zlz
        endif
      enddo

      do i=1,nlocal_Big
        pz = GPrtcl_Big_PosR(i)%z 
        pzt= pz- GPrtcl_Big_PosR(i)%w
        if(pzt-SMALL<zstCoord) then
          ng=ng+1
          if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
          GhostP_Big_id(ng)       = i
          GhostP_Big_id_Global(ng)= GPrtcl_Big_id(i)
          GhostP_Big_nSmall(ng)   = GPrtcl_Big_ind_Sml(2,i)-GPrtcl_Big_ind_Sml(1,i)+1
          GhostP_Big_PosR(ng)     = GPrtcl_Big_PosR(i)         
          GhostP_Big_Direction(ng)= THIS_PROC_ZM
          GhostP_Big_PosR(ng)%z   = pz +zlz
        endif
      enddo  
    ELSEIF(ProcNgh(2) /= MPI_PROC_NULL) then
      do i=1,ngpp    ! consider the previous ghost particle firstly
        pz = GhostP_Big_PosR(i)%z
        pzt= pz- GhostP_Big_PosR(i)%w
        if(pzt-SMALL<zstCoord) then
          nsend=nsend+1
          if(nsend > msend) call reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
      nsendg=nsend

      do i=1,nlocal_Big
        pzt= GPrtcl_Big_PosR(i)%z- GPrtcl_Big_PosR(i)%w
        if(pzt-SMALL<zstCoord) then
          nsend=nsend+1
          if(nsend > msend) call reallocate_sendlist(nsend)
          sendlist(nsend)=i
        endif
      enddo
    ENDIF
    call MPI_SENDRECV(nsend, 1, MPI_INTEGER, ProcNgh(2), 7, &
                      nrecv, 1, MPI_INTEGER, ProcNgh(1), 7, MPI_COMM_WORLD,SRstatus,ierr)
    if(nrecv>0) then
      nrecv2=nrecv*GIBM_size_forward
      allocate(buf_recv(nrecv2))
      call MPI_IRECV(buf_recv,nrecv2,real_type,ProcNgh(1),8,MPI_COMM_WORLD,request(4),ierr)
    endif
    if(nsend>0) then
      nsend2=nsend*GIBM_size_forward
      allocate(buf_send(nsend2))
      call pack_IBM_forward(buf_send,nsendg,nsend,zm_dir)
      call MPI_SEND(buf_send,nsend2,real_type,ProcNgh(2),8,MPI_COMM_WORLD,ierr)
    endif
    if(nrecv>0) then
      call MPI_WAIT(request(4),SRstatus,ierr)
      ng=ng+nrecv
      if(ng>mGhostIBM) call Reallocate_ghost_for_IBM(ng, mGhostIBM)
      call unpack_IBM_forward(buf_recv,ngp+1,ng)
    endif
    if(allocated(buf_send)) deallocate(buf_send) 
    if(allocated(buf_recv)) deallocate(buf_recv)
    
    nGhost_Big= ng
    Block
    integer,dimension(:),allocatable:: IntVec
    integer,dimension(:,:),allocatable:: IntMat
    type(real3),dimension(:),allocatable::Real3Vec
    type(real4),dimension(:),allocatable::Real4Vec
    if(nGhost_Big>0) then
      ierr=0
      ! ======= integer vector part =======
      call move_alloc(GhostP_Big_id, IntVec)
      allocate(GhostP_Big_id(nGhost_Big),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
      GhostP_Big_id=IntVec(1:nGhost_Big)

      call move_alloc(GhostP_Big_id_Global, IntVec)
      allocate(GhostP_Big_id_Global(nGhost_Big),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
      GhostP_Big_id_Global=IntVec(1:nGhost_Big)

      call move_alloc(GhostP_Big_nSmall, IntVec)
      allocate(GhostP_Big_nSmall(nGhost_Big),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
      GhostP_Big_nSmall=IntVec(1:nGhost_Big)
            
      call move_alloc(GhostP_Big_Direction, IntVec)
      allocate(GhostP_Big_Direction(nGhost_Big),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
      GhostP_Big_Direction=IntVec(1:nGhost_Big)            
      deallocate(IntVec)
      
      ! ======= integer matrix part =======
      call move_alloc(GhostP_big_ind_Sml, IntMat)
      allocate(GhostP_big_ind_Sml(2,nGhost_Big),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
      GhostP_big_ind_Sml = IntMat(:,1:nGhost_Big)
      deallocate(IntMat)
      
      ! ======= real3 vector part =======
      call move_alloc(GhostP_Big_FpForce, Real3Vec)
      allocate(GhostP_Big_FpForce(nGhost_Big),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
      GhostP_Big_FpForce=Real3Vec(1:nGhost_Big)

      call move_alloc(GhostP_Big_FpTorque, Real3Vec)
      allocate(GhostP_Big_FpTorque(nGhost_Big),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
      GhostP_Big_FpTorque=Real3Vec(1:nGhost_Big)
      deallocate(Real3Vec)      

      ! ======= real4 vector part =======
      call move_alloc(GhostP_Big_PosR, Real4Vec)
      allocate(GhostP_Big_PosR(nGhost_Big),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
      GhostP_Big_PosR=Real4Vec(1:nGhost_Big)
      deallocate(Real4Vec)            
      if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"PComm_IBM_forward", "Allocation failed-2")

    else
      if(allocated(GhostP_Big_id))        deallocate(GhostP_Big_id)
      if(allocated(GhostP_Big_id_Global)) deallocate(GhostP_Big_id_Global)
      if(allocated(GhostP_Big_nSmall))    deallocate(GhostP_Big_nSmall)
      if(allocated(GhostP_Big_Direction)) deallocate(GhostP_Big_Direction)
      if(allocated(GhostP_Big_FpForce))   deallocate(GhostP_Big_FpForce)
      if(allocated(GhostP_Big_FpTorque))  deallocate(GhostP_Big_FpTorque)
      if(allocated(GhostP_Big_PosR))      deallocate(GhostP_Big_PosR)
      if(allocated(GhostP_big_ind_Sml))   deallocate(GhostP_big_ind_Sml)
    endif
    End Block
  end subroutine PComm_IBM_forward  

  !**********************************************************************
  ! Reallocate_ghost_for_IBM
  !**********************************************************************
  subroutine Reallocate_ghost_for_IBM(ng,mGhost)
    implicit none
    integer,intent(in)::ng
    integer,intent(inout)::mGhost
    
    ! locals
    integer::sizep,sizen,ierrTmp,ierr
    integer,dimension(:),allocatable:: IntVec
    integer,dimension(:,:),allocatable:: IntMat
    type(real3),dimension(:),allocatable::Real3Vec
    type(real4),dimension(:),allocatable::Real4Vec
    
    ierr = 0
    sizep= mGhost
    sizen= int(1.2_RK*real(sizep,kind=RK))
    sizen= min(sizen,numPrtcl_Big)
    sizen= max(sizen,ng+1)
    mGhost=sizen

    ! ======= integer vector part =======
    call move_alloc(GhostP_Big_id,IntVec)
    allocate(GhostP_Big_id(sizen),stat=ierrTmp);ierr=ierr+abs(ierrTmp)
    GhostP_Big_id(1:sizep)=IntVec

    call move_alloc(GhostP_Big_id_Global,IntVec)
    allocate(GhostP_Big_id_Global(sizen),stat=ierrTmp);ierr=ierr+abs(ierrTmp)
    GhostP_Big_id_Global(1:sizep)=IntVec

    call move_alloc(GhostP_Big_nSmall,IntVec)
    allocate(GhostP_Big_nSmall(sizen),stat=ierrTmp);ierr=ierr+abs(ierrTmp)
    GhostP_Big_nSmall(1:sizep)=IntVec
        
    call move_alloc(GhostP_Big_Direction,IntVec)
    allocate(GhostP_Big_Direction(sizen),stat=ierrTmp);ierr=ierr+abs(ierrTmp)
    GhostP_Big_Direction(1:sizep)=IntVec
    deallocate(IntVec)

    ! ======= integer matrix part =======
    call move_alloc(GhostP_big_ind_Sml, IntMat)
    allocate(GhostP_big_ind_Sml(2,sizen),Stat=ierrTmp); ierr=ierr+abs(ierrTmp)
    GhostP_big_ind_Sml(1:2,1:sizep) = IntMat
    deallocate(IntMat)
      
    ! ======= real3 vercor part =======
    call move_alloc(GhostP_Big_FpForce,Real3Vec)
    allocate(GhostP_Big_FpForce(sizen),stat=ierrTmp);ierr=ierr+abs(ierrTmp)
 
    call move_alloc(GhostP_Big_FpTorque,Real3Vec)
    allocate(GhostP_Big_FpTorque(sizen),stat=ierrTmp);ierr=ierr+abs(ierrTmp)
    deallocate(Real3Vec)
        
    ! ======= real4 vercor part =======
    call move_alloc(GhostP_Big_PosR,Real4Vec)
    allocate(GhostP_Big_PosR(sizen),stat=ierrTmp);ierr=ierr+abs(ierrTmp)
    GhostP_Big_PosR(1:sizep)=Real4Vec
    deallocate(Real4Vec)
    
    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"Reallocate_ghost_for_IBM", "Allocation failed")    
  end subroutine Reallocate_ghost_for_IBM

  !**********************************************************************
  ! Reallocate_sendlist
  !**********************************************************************  
  subroutine Reallocate_sendlist(ns)
    implicit none
    integer,intent(in)::ns

    ! locals
    integer::sizep,sizen,ierr
    integer,dimension(:),allocatable:: IntVec

    sizep=msend
    sizen=int(1.2_RK*real(sizep,kind=RK))
    sizen=min(sizen,numPrtcl_Big)
    sizen=max(sizen,ns+1)
    msend=sizen
   
    call move_alloc(sendlist,IntVec)
    allocate(sendlist(sizen),stat=ierr)
    sendlist(1:sizep)=IntVec
    deallocate(IntVec)

    if(ierr/=0) call MainLog%CheckForError(ErrT_Abort,"Reallocate_sendlist", "Allocation failed") 
  end subroutine Reallocate_sendlist

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! pack_IBM_forward
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine pack_IBM_forward(buf_send,nsendg,nsend,dir); implicit none
    real(RK),dimension(:),intent(out)::buf_send
    integer,intent(in)::nsendg,nsend,dir

    ! locals
    integer::i,id,m, coord1,coord2
    real(RK)::pdx,pdy,pdz,dx_pbc(6),dy_pbc(6),dz_pbc(6)
    
    coord1 = int(nrank /p_col)
    coord2 = mod(nrank, p_col)
    dx_pbc=0.0_RK; dy_pbc=0.0_RK; dz_pbc=0.0_RK
    if(BcOption(xp_dir)==BC_PERIOD) then
      if(coord1==0)       dx_pbc(xm_dir)= xlx
      if(coord1==p_row-1) dx_pbc(xp_dir)=-xlx 
    endif
    if(BcOption(yp_dir)==BC_PERIOD) then
      dy_pbc(ym_dir)= yly
      dy_pbc(yp_dir)=-yly
    endif
    if(BcOption(zp_dir)==BC_PERIOD) then
      if(coord2==0)       dz_pbc(zm_dir)= zlz
      if(coord2==p_col-1) dz_pbc(zp_dir)=-zlz
    endif
    
    m=1
    pdx=dx_pbc(dir)
    pdy=dy_pbc(dir)
    pdz=dz_pbc(dir)
    do i=1,nsendg
      id=sendlist(i)
      buf_send(m)=real(nlocal_Big+ id,RK);       m=m+1 ! 01
      buf_send(m)=GhostP_Big_id_Global(id);      m=m+1 ! 02
      buf_send(m)=GhostP_Big_nSmall(id);         m=m+1 ! 03
      buf_send(m)=real(dir,RK);                  m=m+1 ! 04
      buf_send(m)=GhostP_Big_PosR(id)%x+pdx;     m=m+1 ! 05
      buf_send(m)=GhostP_Big_PosR(id)%y+pdy;     m=m+1 ! 06
      buf_send(m)=GhostP_Big_PosR(id)%z+pdz;     m=m+1 ! 07
      buf_send(m)=GhostP_Big_PosR(id)%w;         m=m+1 ! 08
    enddo
    do i=nsendg+1,nsend
      id=sendlist(i)
      buf_send(m)=real(id,RK);                   m=m+1 ! 01
      buf_send(m)=GPrtcl_Big_id(id);             m=m+1 ! 02
      buf_send(m)=GPrtcl_Big_ind_Sml(2,id)-GPrtcl_Big_ind_Sml(1,id)+1; m=m+1 ! 03
      buf_send(m)=real(dir,RK);                  m=m+1 ! 04
      buf_send(m)=GPrtcl_Big_PosR(id)%x+pdx;     m=m+1 ! 05
      buf_send(m)=GPrtcl_Big_PosR(id)%y+pdy;     m=m+1 ! 06
      buf_send(m)=GPrtcl_Big_PosR(id)%z+pdz;     m=m+1 ! 07
      buf_send(m)=GPrtcl_Big_PosR(id)%w;         m=m+1 ! 08     
    enddo
  end subroutine pack_IBM_forward

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! unpack_IBM_forward
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine unpack_IBM_forward(buf_recv,n1,n2); implicit none
    real(RK),dimension(:),intent(in)::buf_recv
    integer,intent(in)::n1,n2

    ! locals
    integer::i,m
   
    m=1
    do i=n1,n2
      GhostP_Big_id(i)       =nint(buf_recv(m)); m=m+1 ! 01
      GhostP_Big_id_Global(i)=nint(buf_recv(m)); m=m+1 ! 02
      GhostP_Big_nSmall(i)   =nint(buf_recv(m)); m=m+1 ! 03
      GhostP_Big_Direction(i)=nint(buf_recv(m)); m=m+1 ! 04
      GhostP_Big_PosR(i)%x   =buf_recv(m);       m=m+1 ! 05
      GhostP_Big_PosR(i)%y   =buf_recv(m);       m=m+1 ! 06
      GhostP_Big_PosR(i)%z   =buf_recv(m);       m=m+1 ! 07
      GhostP_Big_PosR(i)%w   =buf_recv(m);       m=m+1 ! 08 
    enddo
  end subroutine unpack_IBM_forward
        
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! PrepareIBM_ForceAmplify
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine PrepareIBM_ForceAmplify(VolForce_x,VolForce_y,VolForce_z); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::VolForce_x,VolForce_y,VolForce_z

    ! locals
    type(real3):: LagrangeP
    real(RK),dimension(0:nDistribute)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir,VolRatio
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    if(allocated(IBP_ForceAmplify)) deallocate(IBP_ForceAmplify)
    if(nIBP >0) allocate(IBP_ForceAmplify(nIBP))
    VolForce_x=0.0_RK; VolForce_y=0.0_RK; VolForce_z=0.0_RK
    
    ! step1: spread unit Lagrangian force
    DO pid=1,nIBP
      VolRatio   = IBP_VolRatio(pid)
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "pIBM_clc_InterpolateCoe_inc.f90"

      !=====   spreading   =====!
#define spreadIBM_force_inc_unity
#include "pIBM_spreadIBM_force_inc.f90"
#undef  spreadIBM_force_inc_unity
    ENDDO
    call Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)
    call SetBC_and_UpdateHalo_VolForce(VolForce_x,VolForce_y,VolForce_z)
    
    ! step 2: Interpolate unify force to Lagrangian point
    DO pid=1,nIBP
      VolRatio   = IBP_VolRatio(pid)
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "pIBM_clc_InterpolateCoe_inc.f90"

      !===== interpolation =====!
#define spreadIBM_force_inc_unity
#include "pIBM_clc_Point_Interpolation_inc.f90"
#undef  spreadIBM_force_inc_unity

      !=====    forcing    =====!
      IBP_ForceAmplify(pid)= real3(1.0_RK/SumXDir,1.0_RK/SumYDir,1.0_RK/SumZDir)
    ENDDO
  end subroutine PrepareIBM_ForceAmplify
        
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Clc_NoSlipErr
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Clc_NoSlipErr(ux,uy,uz); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::ux,uy,uz

    ! locals
    type(real3):: LagrangeP
    real(RK),dimension(0:nDistribute)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir,VolRatio
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp
    real(RK):: NoSlipErr(4), ErrTmp(3), ErrSum(4)
    
    ErrTmp=-1.0_RK; NoSlipErr=-1.0_RK; ErrSum=0.0_RK
    DO pid=1,nIBP
      VolRatio   = IBP_VolRatio(pid)
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "pIBM_clc_InterpolateCoe_inc.f90"

      !===== interpolation =====!
#include "pIBM_clc_Point_Interpolation_inc.f90"

      !=====    forcing    =====!
      ErrSum(1) = ErrSum(1) +abs(IBP_Vel(pid)%x -SumXDir)
      ErrSum(2) = ErrSum(2) +abs(IBP_Vel(pid)%y -SumYDir)
      ErrSum(3) = ErrSum(3) +abs(IBP_Vel(pid)%z -SumZDir)
      ErrSum(4) = ErrSum(4) +1.0_RK
      ErrTmp(1) = max(ErrTmp(1), abs(IBP_Vel(pid)%x -SumXDir))
      ErrTmp(2) = max(ErrTmp(2), abs(IBP_Vel(pid)%y -SumYDir))
      ErrTmp(3) = max(ErrTmp(3), abs(IBP_Vel(pid)%z -SumZDir))
    ENDDO
    call MPI_REDUCE(ErrTmp, NoSlipErr, 3, real_type, MPI_MAX, 0, MPI_COMM_WORLD, i)
    if(nrank==0) print*, itime,'^^^^^^^^^^^^ Max No-Slip Error = ', NoSlipErr(1:3)
    call MPI_REDUCE(ErrSum, NoSlipErr, 4, real_type, MPI_SUM, 0, MPI_COMM_WORLD, i)
    if(nrank==0) print*, itime,'^^^^^^^^^^^^ Ave No-Slip Error = ', NoSlipErr(1:3)/NoSlipErr(4)
  end subroutine Clc_NoSlipErr

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! updateRhsIBM
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine updateRhsIBM(uxStar,uyStar,uzStar,ux,uy,uz,RhsX,RhsY,RhsZ); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in)::uxStar,uyStar,uzStar
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(in) ::ux,uy,uz
    real(RK),dimension(y1start(1):y1end(1),y1start(2):y1end(2),y1start(3):y1end(3)),intent(out)::RhsX,RhsY
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(out)::RhsZ
    
    ! locals
    integer:: ic,jc,kc

    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        do ic=y1start(1),y1end(1)
           RhsX(ic,jc,kc)= uxStar(ic,jc,kc) -ux(ic,jc,kc)
           RhsY(ic,jc,kc)= uyStar(ic,jc,kc) -uy(ic,jc,kc)
           RhsZ(ic,jc,kc)= uzStar(ic,jc,kc) -uz(ic,jc,kc)        
        enddo
      enddo
    ENDDO  
  end subroutine updateRhsIBM

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! NormVolForce_x
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine NormVolForce_x(VolForce_x); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::VolForce_x
    
    ! locals
    integer::ic,jc,kc,ierr
    real(RK)::SumVolForceX,SumVolForceXTot,ForcedCoe
    
    SumVolForceX=0.0_RK
    DO kc=y1start(3),y1end(3)
      do jc=y1start(2),y1end(2)
        ForcedCoe=dyp(jc)
        do ic=y1start(1),y1end(1)
          SumVolForceX=SumVolForceX+ForcedCoe*VolForce_x(ic,jc,kc)
        enddo
      enddo
    ENDDO
    call MPI_ALLREDUCE(SumVolForceX,SumVolForceXTot,1,real_type,MPI_SUM,MPI_COMM_WORLD,ierr)
    SumVolForceX = -SumVolForceXTot/(real(nxc,RK)*real(nzc,RK))/yly
    VolForce_x =VolForce_x+SumVolForceX
    
    PrGradData(3)=PrGradData(3) +SumVolForceX/pmAlpha
    PrGradData(1)=PrGradData(1) +SumVolForceX/dt
  end subroutine NormVolForce_x
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! AdditionalForceIBM_0
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine AdditionalForceIBM_0(ux,uy,uz,VolForce_x,VolForce_y,VolForce_z); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,VolForce_x,VolForce_y,VolForce_z

    ! locals
    type(real3):: LagrangeP,IbpForce
    real(RK),dimension(0:nDistribute)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir,VolRatio
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    VolForce_x=0.0_RK; VolForce_y=0.0_RK; VolForce_z=0.0_RK
    DO pid=1,nIBP
      VolRatio   = IBP_VolRatio(pid)
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "pIBM_clc_InterpolateCoe_inc.f90"

      !===== interpolation =====!
#include "pIBM_clc_Point_Interpolation_inc.f90"

      !=====    forcing    =====!
      IbpForce= (IBP_Vel(pid) -real3(SumXDir,SumYDir,SumZDir))*ForceAmplifyCoe
      IBP_Force(pid)= IBP_Force(pid)+IbpForce

      !=====   spreading   =====!
#include "pIBM_spreadIBM_force_inc.f90"
    ENDDO
    call Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)
    IF(IsUxConst) call NormVolForce_x(VolForce_x)
    
    ! Velocity corrcetion
    DO k=y1start(3),y1end(3)
      do j=y1start(2),y1end(2)
        do i=y1start(1),y1end(1)
           ux(i,j,k)= ux(i,j,k)+ VolForce_x(i,j,k)
           uy(i,j,k)= uy(i,j,k)+ VolForce_y(i,j,k)
           uz(i,j,k)= uz(i,j,k)+ VolForce_z(i,j,k)        
        enddo
      enddo
    ENDDO
  end subroutine AdditionalForceIBM_0

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! AdditionalForceIBM_1
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine AdditionalForceIBM_1(ux,uy,uz,VolForce_x,VolForce_y,VolForce_z); implicit none
    real(RK),dimension(mb1%xmm:mb1%xpm,mb1%ymm:mb1%ypm,mb1%zmm:mb1%zpm),intent(inout)::ux,uy,uz,VolForce_x,VolForce_y,VolForce_z

    ! locals
    type(real3):: LagrangeP,IbpForce
    real(RK),dimension(0:nDistribute)::RatioXc,RatioYc,RatioZc,RatioXp,RatioYp,RatioZp
    real(RK):: prx,pry,prz,prxc,pryc,przc,prxp,pryp,przp,SumXDir,SumYDir,SumZDir,VolRatio
    integer::i,j,k,id,jd,kd,pid,idxp_interp,idyp_interp,idzp_interp,idxc_interp,idyc_interp,idzc_interp

    VolForce_x=0.0_RK; VolForce_y=0.0_RK; VolForce_z=0.0_RK
    DO pid=1,nIBP
      VolRatio   = IBP_VolRatio(pid)
      LagrangeP  = IBP_Pos(pid)
      idxc_interp= IBP_indxyz(1,pid)
      idxp_interp= IBP_indxyz(2,pid)
      idyc_interp= IBP_indxyz(3,pid)
      idyp_interp= IBP_indxyz(4,pid)
      idzc_interp= IBP_indxyz(5,pid)
      idzp_interp= IBP_indxyz(6,pid)
#include "pIBM_clc_InterpolateCoe_inc.f90"

      !===== interpolation =====!
#include "pIBM_clc_Point_Interpolation_inc.f90"

      !=====    forcing    =====!
      IbpForce= (IBP_Vel(pid) -real3(SumXDir,SumYDir,SumZDir))*IBP_ForceAmplify(pid)
      IBP_Force(pid)= IBP_Force(pid)+IbpForce

      !=====   spreading   =====!
#include "pIBM_spreadIBM_force_inc.f90"
    ENDDO
    call Gather_Halo_IBMForce(VolForce_x,VolForce_y,VolForce_z)
    IF(IsUxConst) call NormVolForce_x(VolForce_x)
    
    ! Velocity corrcetion
    DO k=y1start(3),y1end(3)
      do j=y1start(2),y1end(2)
        do i=y1start(1),y1end(1)
           ux(i,j,k)= ux(i,j,k)+ VolForce_x(i,j,k)
           uy(i,j,k)= uy(i,j,k)+ VolForce_y(i,j,k)
           uz(i,j,k)= uz(i,j,k)+ VolForce_z(i,j,k)        
        enddo
      enddo
    ENDDO
  end subroutine AdditionalForceIBM_1

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! Reallocate_IbpVar
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine Reallocate_IbpVar(); implicit none

    ! locals
    integer:: sizep,sizen,ierrTemp, ierr
    integer,dimension(:),allocatable:: IntVec
    real(RK),dimension(:),allocatable::RealVec
    type(real3),dimension(:),allocatable:: Real3Vec
    integer(kind=2),dimension(:,:),allocatable::IntArr   
 
    ierr=0
    sizep= mIBP
    sizen= int(1.1_RK*real(sizep,kind=RK))
    sizen= max(sizen, nIBP+1)
    mIBP = sizen    

    ! ======= integer vector part =======
    call move_alloc(IBP_idlocal,IntVec)
    allocate(IBP_idlocal(sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)
    IBP_idlocal(1:sizep)=IntVec
    deallocate(IntVec)
    
    ! ======= real vector part =======
    call move_alloc(IBP_VolRatio,RealVec)
    allocate(IBP_VolRatio(sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)
    IBP_VolRatio(1:sizep)=RealVec
    deallocate(RealVec)

    ! ======= integer matrix part =======
    call move_alloc(IBP_indxyz,IntArr)
    allocate(IBP_indxyz(6,sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)
    IBP_indxyz(1:6,1:sizep)=IntArr
    deallocate(IntArr)

    ! ======= real3 vercor part =======
    call move_alloc(IBP_Pos,Real3Vec)
    allocate(IBP_Pos(sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)
    IBP_Pos(1:sizep)=Real3Vec

    call move_alloc(IBP_Vel,Real3Vec)
    allocate(IBP_Vel(sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)
    IBP_Vel(1:sizep)=Real3Vec
    deallocate(Real3Vec)

    deallocate(IBP_Force)
    allocate(IBP_Force(sizen),stat=ierrTemp);ierr=ierr+abs(ierrTemp)

    if(ierr/=0) then
      call MainLog%CheckForError(ErrT_Abort," Reallocate_IbpVar"," Reallocate wrong!")
      call MainLog%OutInfo("The present processor  is :"//strip(num2str(nrank)),3)
    endif    
    !call MainLog%CheckForError(ErrT_Pass,"Reallocate_IbpVar","Need to reallocate IBP variables")
    !call MainLog%OutInfo("The present processor  is :"//strip(num2str(nrank)),3)
    !call MainLog%OutInfo("Previous matirx length is :"//strip(num2str(sizep)),3)
    !call MainLog%OutInfo("Updated  matirx length is :"//strip(num2str(sizen)),3)
  end subroutine Reallocate_IbpVar

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! clc_bgn_ind
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  function clc_bgn_ind(nlocal) result(bgn_ind); implicit none
    integer,intent(in)::nlocal
    integer::bgn_ind

    ! locals
    integer::end_ind,ierr,SRstatus(MPI_STATUS_SIZE)
  
    bgn_ind=0
    if(nproc<=1) return
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    IF(nrank==0)THEN
      end_ind=bgn_ind+nlocal
      call MPI_SEND(end_ind, 1,MPI_INTEGER, nrank+1,0,MPI_COMM_WORLD,ierr)
    ELSEIF(nrank /= nproc-1) THEN
      call MPI_RECV(bgn_ind, 1,MPI_INTEGER, nrank-1,0,MPI_COMM_WORLD,SRstatus,ierr)
      end_ind=bgn_ind+nlocal
      call MPI_SEND(end_ind, 1,MPI_INTEGER, nrank+1,0,MPI_COMM_WORLD,ierr)
    ELSE
      call MPI_RECV(bgn_ind, 1,MPI_INTEGER, nrank-1,0,MPI_COMM_WORLD,SRstatus,ierr)
    ENDIF
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end function clc_bgn_ind
  
#undef nDistribute
#undef THIS_PROC_XP
#undef THIS_PROC_XM
#undef THIS_PROC_YP
#undef THIS_PROC_YM
#undef THIS_PROC_ZP
#undef THIS_PROC_ZM

end module ca_IBM
