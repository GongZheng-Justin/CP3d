module Prtcl_Property
  use m_TypeDef
  use m_LogInfo
  use Prtcl_Parameters
#ifdef CFDDEM
  use m_Decomp2d,only: nrank
  use m_Parameters,only: FluidDensity
#elif  CFDACM
  use Prtcl_EqualSphere
  use m_Decomp2d,only: nrank
  use m_Parameters,only: FluidDensity,xnu
  use m_MeshAndMetries,only: dx,dyUniform,dz
#else
  use Prtcl_decomp_2d,only: nrank
#endif
  implicit none
  private
    
  type PureProperty
    real(RK):: Radius = 0.005_RK
    real(RK):: Density= 2500.0_RK
    real(RK):: PoissonRatio= 0.25_RK
    real(RK):: YoungsModulus=5.0E6_RK
    real(RK):: Mass
    real(RK):: Volume
    real(RK):: Inertia
#ifdef CFDDEM
    real(RK):: MassInFluid
    real(RK):: MassOfFluid
#endif
  end type PureProperty

  type,public:: BinaryProperty
    real(RK):: RadEff       ! effective Radiusius
    real(RK):: MassEff      ! effective Mass

    real(RK):: StiffnessCoe_n
    real(RK):: StiffnessCoe_t        
    real(RK):: DampingCoe_n = 0.0_RK
    real(RK):: DampingCoe_t = 0.0_RK
    real(RK):: RestitutionCoe_n = 0.95_RK  ! Normal Resitution Coefficient
    real(RK):: FrictionCoe_s = 0.80_RK     ! Coefficient of static  friction
    real(RK):: FrictionCoe_k = 0.15_RK     ! Coefficient of kinetic friction
#ifdef CFDACM
    real(RK):: Vel_Crit                    ! Turn off fluid forces for large St collsions
    real(RK):: Kn_Grav                
#endif
  end type BinaryProperty
    
  type PhysicalProperty
    integer,allocatable,dimension(:) :: nPrtcl_in_Bin
    integer,allocatable,dimension(:) :: CS_Hrchl_level ! level for NBS-Munjiza-Hierarchy Contact Search Method  
    type(pureProperty),allocatable,dimension(:) :: Prtcl_PureProp
    type(pureProperty),allocatable,dimension(:) :: Wall_PureProp
    type(BinaryProperty),allocatable,dimension(:,:):: Prtcl_BnryProp    
    type(BinaryProperty),allocatable,dimension(:,:):: PrtclWall_BnryProp         
  contains
    procedure:: InitPrtclProperty
    procedure:: InitWallProperty
  end type PhysicalProperty
  type(PhysicalProperty),public::DEMProperty
  
#ifdef CFDACM
  type IBMProperty
    integer::  nPartition
    real(RK):: IBPVolume
    real(RK):: MassIBL
    real(RK):: InertiaIBL
    real(RK):: MassInFluid
    real(RK):: MassEff
    real(RK):: InertiaEff
  end type IBMProperty
  type(IBMProperty),dimension(:),allocatable,public:: PrtclIBMProp
  type(real3),allocatable,dimension(:,:),public::SpherePartitionCoord
  real(RK),allocatable,dimension(:,:),public:: dlub_pp,LubCoe_pp
  real(RK),allocatable,dimension(:),  public:: dlub_pw,LubCoe_pw
#endif

contains
  !*******************************************************
  ! initializing the size distribution with property 
  !*******************************************************
  subroutine InitPrtclProperty(this,chFile)
    implicit none
    class(PhysicalProperty)::this
    character(*),intent(in)::chFile
        
    ! locals
    real(RK)::FrictionCoe_s_PP,FrictionCoe_k_PP,RestitutionCoe_n_PP
    real(RK),dimension(:),allocatable:: Bin_Divided,Density,Diameter,YoungsModulus_P,PoissonRatio_P
    namelist/ParticlePhysicalProperty/Bin_Divided, Density, Diameter,YoungsModulus_P,PoissonRatio_P, &
                                      FrictionCoe_s_PP,FrictionCoe_k_PP,RestitutionCoe_n_PP
    integer:: i,j,iTV(8),nPType,nUnitFile,ierror,sum_prtcl,bin_pnum,prdiff,bin_id
    real(RK):: sum_divided,rtemp,Radius
    type(PureProperty)::pari,parj
    type(BinaryProperty)::Bnry
#ifdef CFDACM
    integer::nSpSize
    real(RK)::dxyz
    integer,dimension(:),allocatable::  nPartition
    real(RK),dimension(:),allocatable:: RetractionRatio
    real(RK),dimension(:,:),allocatable:: PointTemp
    namelist/ParticleIBMProperty/ nPartition, RetractionRatio
#endif
        
    nPType  = DEM_Opt%numPrtcl_Type
    allocate( Bin_Divided(nPType))
    allocate( Density(nPType))
    allocate( Diameter(nPType))
    allocate( YoungsModulus_P(nPType))
    allocate( PoissonRatio_P(nPType))
        
    allocate( this%CS_Hrchl_level(nPType))
    allocate( this%nPrtcl_in_Bin(nPType))
    allocate( this%Prtcl_PureProp(nPType))
    allocate( this%Prtcl_BnryProp(nPType,nPType)) 
           
    open(newunit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=ierror)
    if(ierror /= 0 ) call DEMLogInfo%CheckForError(ErrT_Abort,"InitPrtclProperty","Cannot open file:"//trim(chFile))
    read(nUnitFile, nml=ParticlePhysicalProperty)
    if(nrank==0)write(DEMLogInfo%nUnit, nml=ParticlePhysicalProperty)
#ifdef CFDACM
    dxyz= (dx*dyUniform*dz)**(0.333333333333333333_RK)
    allocate( nPartition(nPType))
    allocate( RetractionRatio(nPType))
    allocate( PrtclIBMProp(nPType))
    allocate(dlub_pp(nPType,nPType))
    allocate(LubCoe_pp(nPType,nPType))
    rewind(nUnitFile)
    read(nUnitFile, nml=ParticleIBMProperty)
    if(nrank==0)write(DEMLogInfo%nUnit, nml=ParticleIBMProperty)
#endif
    close(nUnitFile,IOSTAT=ierror)
        
    ! calculate this%nPrtcl_in_Bin
    call date_and_time(values=iTV); !iTV=0
    call random_seed(size= i)
    call random_seed(put = iTV(7)*iTV(8)+[(j,j=1,i)])
    sum_divided=0.0_RK
    do i=1,nPType
      sum_divided = sum_divided + Bin_Divided(i)
    enddo
    sum_prtcl = 0
    do i=1,nPType
      bin_pnum = int(DEM_Opt%numPrtcl*Bin_Divided(i)/sum_divided)
      this%nPrtcl_in_Bin(i)= bin_pnum
      sum_prtcl = sum_prtcl + bin_pnum
    enddo
    prdiff = DEM_Opt%numPrtcl - sum_prtcl
    if( prdiff > 0 ) then
      do i=1, prdiff 
        call random_number(rtemp)
        bin_id = int(rtemp*nPType) + 1
        this%nPrtcl_in_Bin(bin_id) = this%nPrtcl_in_Bin(bin_id)  + 1
      enddo
    endif
        
    ! calculate particle properties
    do i = 1, nPType
      Radius= 0.5_RK*Diameter(i)
      this%Prtcl_PureProp(i)%Density= Density(i)
      this%Prtcl_PureProp(i)%Radius = Radius
      this%Prtcl_PureProp(i)%YoungsModulus = YoungsModulus_P(i)
      this%Prtcl_PureProp(i)%PoissonRatio = PoissonRatio_P(i)

      this%Prtcl_PureProp(i)%Volume = 1.333333333333333333_RK*Pi*Radius**3
      this%Prtcl_PureProp(i)%Mass = Density(i)*this%Prtcl_PureProp(i)%Volume
      this%Prtcl_PureProp(i)%Inertia = 0.4_RK*this%Prtcl_PureProp(i)%Mass*Radius**2
#ifdef CFDDEM
      this%Prtcl_PureProp(i)%MassInFluid= (Density(i)-FluidDensity)*this%Prtcl_PureProp(i)%Volume
      this%Prtcl_PureProp(i)%MassOfFluid= FluidDensity*this%Prtcl_PureProp(i)%Volume
#endif
#ifdef CFDACM
      rtemp= Diameter(i)*0.5_RK -RetractionRatio(i)*dxyz
      PrtclIBMProp(i)%nPartition = nPartition(i)
      PrtclIBMProp(i)%IBPVolume= PI*dxyz*(4.0_RK*rtemp*rtemp+ dxyz*dxyz/3.0_RK)/real(nPartition(i),RK)
      PrtclIBMProp(i)%MassIBL  = FluidDensity* PI*dxyz*(4.0_RK*rtemp*rtemp+ dxyz*dxyz/3.0_RK)
      PrtclIBMProp(i)%InertiaIBL=0.6666666666666666667_RK*(rtemp**2)*PrtclIBMProp(i)%MassIBL
      PrtclIBMProp(i)%MassInFluid= (Density(i)-FluidDensity)*this%Prtcl_PureProp(i)%Volume

      select case(IBM_Scheme)
      case(0)   ! 0: Explicit,Uhlmann(2005,JCP)
        PrtclIBMProp(i)%MassEff   = (1.0_RK-FluidDensity/Density(i))*this%Prtcl_PureProp(i)%Mass
        PrtclIBMProp(i)%InertiaEff= (1.0_RK-FluidDensity/Density(i))*this%Prtcl_PureProp(i)%Inertia

      case(1)   ! 1: Explicit,Kempe(2012,JCP)
        PrtclIBMProp(i)%MassEff   = this%Prtcl_PureProp(i)%Mass
        PrtclIBMProp(i)%InertiaEff= this%Prtcl_PureProp(i)%Inertia

      case(2)   ! 2: Semi-implicit,Tschisgale(2017,JCP)
        PrtclIBMProp(i)%MassEff   = this%Prtcl_PureProp(i)%Mass    + PrtclIBMProp(i)%MassIBL
        PrtclIBMProp(i)%InertiaEff= this%Prtcl_PureProp(i)%Inertia + PrtclIBMProp(i)%InertiaIBL

      case default
        call DEMLogInfo%CheckForError(ErrT_Abort,"InitPrtclProperty","IBM scheme wrong!!!")
      end select
#endif
    enddo
        
    do j=1,nPType
      parj = this%Prtcl_PureProp(j)
      do i=1,nPType
        pari = this%Prtcl_PureProp(i)       
        Bnry= clc_BnryPrtcl_Prop(pari,parj,FrictionCoe_s_PP,FrictionCoe_k_PP,RestitutionCoe_n_PP,.false.)
        this%Prtcl_BnryProp(i,j)= Bnry  
#ifdef CFDACM
        dlub_pp(i,j)= Lub_ratio* dxyz
        LubCoe_pp(i,j)= Klub_pp*xnu*FluidDensity*Bnry%RadEff*Bnry%RadEff/dlub_pp(i,j)
#endif
      enddo
    enddo

#ifdef CFDACM
    ! Sphere partition informations
    nSpSize=maxval(nPartition)
    allocate(SpherePartitionCoord(nSpSize,nPType))
    SpherePartitionCoord=zero_r3
    do i=1,nPType
      allocate(PointTemp(nPartition(i),3))
      call eq_Sphere(PointTemp)
      do j=1,nPartition(i)
        rtemp= Diameter(i)*0.5_RK -RetractionRatio(i)*dxyz
        SpherePartitionCoord(j,i)=real3(PointTemp(j,1),PointTemp(j,2),PointTemp(j,3))*rtemp
      enddo
      deallocate(PointTemp)
    enddo
#endif
  end subroutine InitPrtclProperty

  !*****************************************************************************    
  ! setting the physical property of wall 
  !*****************************************************************************
  subroutine InitWallProperty(this,chFile)
    implicit none
    class(PhysicalProperty)::this
    character(*),intent(in)::chFile
        
    ! locals
    type(PureProperty)pari,wall
    integer::i,j,nPType,nWType,nUnitFile,ierror
    real(RK)::FrictionCoe_s_PW,FrictionCoe_k_PW,RestitutionCoe_n_PW
    real(RK),dimension(:),allocatable:: YoungsModulus_W,PoissonRatio_W
    namelist /WallPhysicalProperty/YoungsModulus_W,PoissonRatio_W,FrictionCoe_s_PW,FrictionCoe_k_PW,RestitutionCoe_n_PW
            
    nPType= DEM_Opt%numPrtcl_Type
    nWType= DEM_Opt%numWall_type
    allocate(YoungsModulus_W(nWType))
    allocate(PoissonRatio_W(nWType))
    allocate(this%Wall_PureProp(nWType))
    allocate(this%PrtclWall_BnryProp(nPType,nWType))
#ifdef CFDACM 
    allocate(dlub_pw(nPType))
    allocate(LubCoe_pw(nPType))
#endif
   
    open(newunit=nUnitFile, file=chFile, status='old', form='formatted', IOSTAT=ierror)
    if(ierror /= 0 .and. nrank==0) then
      call DEMLogInfo%CheckForError(ErrT_Abort,"InitWallProperty: " ,"Cannot open file:"//trim(chFile))
    endif
    read(nUnitFile, nml=WallPhysicalProperty)
    if(nrank==0)write(DEMLogInfo%nUnit, nml=WallPhysicalProperty)
    close(nUnitFile)        
        
    do i=1,nWType
       this%Wall_PureProp(i)%Density = 1.0E50_RK
       this%Wall_PureProp(i)%Radius  = 1.0E50_RK
       this%Wall_PureProp(i)%YoungsModulus= YoungsModulus_W(i)
       this%Wall_PureProp(i)%PoissonRatio = PoissonRatio_W(i)
           
       this%Wall_PureProp(i)%Mass    = 1.0E50_RK
       this%Wall_PureProp(i)%Volume  = 1.0E50_RK
       this%Wall_PureProp(i)%Inertia = 1.0E50_RK
    enddo

    do j = 1, nWType
      wall = this%Wall_PureProp(j)
      do i = 1, nPType   
        pari = this%Prtcl_PureProp(i)        
        this%PrtclWall_BnryProp(i,j)=clc_BnryPrtcl_Prop(pari,wall,FrictionCoe_s_PW,FrictionCoe_k_PW,RestitutionCoe_n_PW,.true.)
      enddo
#ifdef CFDACM
      dlub_pw(i)  = Lub_ratio*(dx*dyUniform*dz)**(0.333333333333333333_RK)
      LubCoe_pw(i)= Klub_pw*xnu*FluidDensity*(this%Prtcl_PureProp(i)%Radius)**2/dlub_pw(i)
#endif
    enddo       
  end subroutine InitWallProperty

  !*****************************************************************************    
  ! Calculating the binary contact properties
  !*****************************************************************************
  function clc_BnryPrtcl_Prop(pari,parj,FrictionCoe_s,FrictionCoe_k,RestitutionCoe_n,iswall) result(Bnry)
    implicit none
    class(PureProperty),intent(in):: pari, parj
    real(RK),intent(in)::FrictionCoe_s,FrictionCoe_k,RestitutionCoe_n
    logical,intent(in) ::iswall

    ! locals
    type(BinaryProperty):: Bnry
    real(RK)::kappa,pri,prj,ShearModi,ShearModj,YoungsModEff,ShearModEff,Eta,K_hertz
#ifdef CFDACM
    real(RK)::DensityEff,Vel_Crit,Kn_Grav,GravityMag,GRAV_OVERLAP,Lamda,Beta,rA,rB,rC,AlphaTau2,Tau_c0
#endif
    
    if(.not.iswall) then
      Bnry%RadEff = (pari%Radius*parj%Radius)/(pari%Radius +parj%Radius)
      Bnry%MassEff= (pari%Mass  *parj%Mass  )/(pari%Mass   +parj%Mass )
    else
      Bnry%RadEff = pari%Radius
      Bnry%MassEff= pari%Mass         
    endif
    Bnry%RestitutionCoe_n = RestitutionCoe_n
    Bnry%FrictionCoe_s = FrictionCoe_s
    Bnry%FrictionCoe_k = FrictionCoe_k

    pri=pari%PoissonRatio
    prj=parj%PoissonRatio
    ShearModi=pari%YoungsModulus/(2.0_RK*(1.0_RK+pri))
    ShearModj=parj%YoungsModulus/(2.0_RK*(1.0_RK+prj))
    YoungsModEff=1.0_RK/((1.0_RK-pri*pri)/pari%YoungsModulus + (1.0_RK-prj*prj)/parj%YoungsModulus) ! 2.47, p31
    ShearModEff =1.0_RK/((2.0_RK- pri)/ShearModi + (2.0_RK-prj)/ShearModj)                          ! 2.71, p43
    kappa=((1.0_RK-pri)/ShearModi+(1.0_RK-prj)/ShearModj)/((1.0_RK-0.5_RK*pri)/ShearModi+(1.0_RK-0.5_RK*prj)/ShearModj)
    kappa=abs(kappa)
    Eta=log(RestitutionCoe_n); Eta=Eta*Eta
    if(DEM_Opt%CF_Type == DEM_LSD) then                                             
      Bnry%StiffnessCoe_n = 1.2024_RK*(sqrt(Bnry%MassEff)*(YoungsModEff**2)*Bnry%RadEff)**(0.4_RK)        ! 2.46, p31
      Bnry%DampingCoe_n   =-2.0_RK*log(RestitutionCoe_n)*sqrt(Bnry%MassEff*Bnry%StiffnessCoe_n)/sqrt(PI*PI+Eta) ! 2.44, p31
      Bnry%StiffnessCoe_t = Bnry%StiffnessCoe_n*kappa                                                     ! 2.52, p34
      Bnry%DampingCoe_t   = Bnry%DampingCoe_n*sqrt(kappa)
    elseif(DEM_Opt%CF_Type == DEM_nLin) then
      K_hertz =1.333333333333333333_RK*YoungsModEff*sqrt(Bnry%RadEff)               ! 2.62, P39
      Bnry%StiffnessCoe_n= K_hertz
      Bnry%DampingCoe_n  = -2.2664_RK*log(RestitutionCoe_n)*sqrt(Bnry%MassEff*K_hertz)/sqrt(Eta+10.1354_RK) ! 2.66, p40
      Bnry%StiffnessCoe_t= 5.3333333333333_RK*ShearModEff*sqrt(Bnry%RadEff) ! 2.72, P44
      Bnry%DampingCoe_t  = 0.0_RK ! No equation is considered for tangential damping yet
    endif
#ifdef CFDACM
    if(DEM_Opt%CF_Type == ACM_LSD) then      ! Costa et al./Physics Review E 92,053012 (2015)
      Bnry%StiffnessCoe_n= Bnry%MassEff*(PI*PI+Eta)
      Bnry%DampingCoe_n  =-2.0_RK*Bnry%MassEff*log(RestitutionCoe_n)
      Bnry%StiffnessCoe_t= Bnry%StiffnessCoe_n*kappa
      Bnry%DampingCoe_t  = Bnry%DampingCoe_n*sqrt(kappa)
    elseif(DEM_Opt%CF_Type == ACM_nLin) then ! E. Biegert et al./Journal of Computational Physics 340(2017): 105-127
      rA=0.716_RK; rB=0.830_RK; 
      rC=0.744_RK; Tau_c0=3.218_RK
      AlphaTau2=1.111_RK*1.111_RK +3.218_RK*3.218_RK
      Lamda= (sqrt(0.25_RK*rC*rC*Eta*Eta+AlphaTau2*Eta)-0.5_RK*rC*Eta)/AlphaTau2
      Beta= Tau_c0/sqrt(abs(1.0_RK-rA*Lamda-rB*Lamda*Lamda))
      Bnry%StiffnessCoe_n= Bnry%MassEff*(Beta**(2.5_RK))
      Bnry%DampingCoe_n  = 2.0_RK*Bnry%MassEff*Lamda*Beta
      Bnry%StiffnessCoe_t= Bnry%MassEff*(PI*PI+Eta)*kappa
      Bnry%DampingCoe_t  =-2.0_RK*Bnry%MassEff*log(RestitutionCoe_n)*sqrt(kappa)
    endif
    DensityEff= (pari%Density * parj%Density)/(pari%Density + parj%Density)
    Vel_Crit= 4.5_RK*St_Crit*FluidDensity*xnu/(Bnry%RadEff*DensityEff)
    if(IsDryColl) then
      Bnry%Vel_Crit= Vel_Crit
    else
      Bnry%Vel_Crit= 1.0E100_RK ! If IsDryColl=F, we set a very huge Vel_Crit
    endif

    GRAV_OVERLAP=0.001_RK ! E. Biegert et al./Journal of Computational Physics 340(2017): 105-127
    GravityMag=norm(DEM_Opt%gravity)
    if(DEM_Opt%CF_Type == ACM_LSD) then
      kn_Grav= pari%Mass*GravityMag/(GRAV_OVERLAP*pari%Radius)
      if(.not.iswall) kn_Grav= max(parj%Mass*GravityMag/(GRAV_OVERLAP*parj%Radius),kn_Grav)
    elseif(DEM_Opt%CF_Type == ACM_nLin) then
      kn_Grav= pari%Mass*GravityMag/(GRAV_OVERLAP*pari%Radius)**(1.5_RK)
      if(.not.iswall) kn_Grav= max(parj%Mass*GravityMag/(GRAV_OVERLAP*parj%Radius)**(1.5_RK),kn_Grav)
    endif
    Bnry%kn_Grav=kn_Grav
#endif
  end function clc_BnryPrtcl_Prop
end module Prtcl_Property
