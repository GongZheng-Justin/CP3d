module ATP_decomp_2d
  use MPI
  use m_TypeDef
  use ATP_Parameters
  use m_MeshAndMetries,only:dx,dz
  use m_Parameters,only:p_row,p_col
  use m_Decomp2d,only:nrank,y1start,y1end,myProcNghBC
  implicit none
  private

  integer,public:: int_type,real_type,real3_type
  integer,public:: int_byte,real_byte,real3_byte
  TYPE Prtcl_DECOMP_INFO

    ! define neighboring blocks
    ! second dimension 8 neighbour processors:
    !        1:4, 4 edge neighbours; 5:8, 4 cornor neighbours; 0, current processor(if any)
    integer::prow
    integer::pcol
    integer::coord1
    integer::coord2
    integer,dimension(4):: ProcNgh
    real(RK)::xSt,ySt,zSt !min domain of the current and neighbor processors 
    real(RK)::xEd,yEd,zEd !max domain of the current and neighbor processors
  contains
    procedure:: Init_DECOMP => PDI_Init_DECOMP
  end type Prtcl_DECOMP_INFO
  TYPE(Prtcl_DECOMP_INFO), public :: ATP_decomp

contains

  !**********************************************************************
  ! PDI_Init_DECOMP
  !**********************************************************************
  subroutine PDI_Init_DECOMP(this)
    implicit none
    class(Prtcl_DECOMP_INFO)::this

    ! locals
    integer::i

    this%prow= p_row
    this%pcol= p_col
    call Init_Prctl_MPI_TYPE()
    
    ! ------------------------- neigbor information begins------------------------
    ! 
    ! 2D domain decomposition method in ATP_MPI(from left to right:  x1-pencil, y1-pencil, z1-pencil)
    ! 
    !   If we have 6 processors = 3 row * 2 col
    ! 
    !     the arrangement of the subdomains(nrank) is as follow:
    !       y               x               x
    !       |        4 5    |        4 5    |        4 5   
    !       |        2 3    |        2 3    |        2 3
    !       |_ _ _z  0 1    |_ _ _z  0 1    |_ _ _y  0 1
    !
    !     the arrangement of the coord1 is as follow:
    !       y               x               x
    !       |        2 2    |        2 2    |        2 2
    !       |        1 1    |        1 1    |        1 1
    !       |_ _ _z  0 0    |_ _ _z  0 0    |_ _ _y  0 0
    ! 
    !     the arrangement of the coord2 is as follow:
    !       y               x               x
    !       |        0 1    |        0 1    |        0 1
    !       |        0 1    |        0 1    |        0 1
    !       |_ _ _ z 0 1    |_ _ _z  0 1    |_ _ _y  0 1
    ! 
    !     neighbor index:
    ! 
    !       y               x               x
    !       |        6 3 5  |        6 3 5  |        6 3 5
    !       |        2 0 1  |        2 0 1  |        2 0 1
    !       |_ _ _ z 7 4 8  |_ _ _z  7 4 8  |_ _ _y  7 4 8
    ! 
    !        Here 0 means the center subdomain, and 1-8 stands for the relative location of the eight neighbors

    this%coord1 = int ( nrank / p_col)
    this%coord2 = mod ( nrank,  p_col)
    DO i=1,4
      if(myProcNghBC(2,i)<0) then
        this%ProcNgh(i)= MPI_PROC_NULL
      else
        this%ProcNgh(i)= myProcNghBC(2,i)
      endif
    ENDDO
    this%xSt= real(y1start(1)-1,kind=RK)*dx
    this%xEd= real(y1end(1),    kind=RK)*dx
    this%zSt= real(y1start(3)-1,kind=RK)*dz
    this%zEd= real(y1end(3),    kind=RK)*dz
    this%ySt= ATP_Opt%SimDomain_min%y
    this%yEd= ATP_Opt%SimDomain_max%y
  end subroutine PDI_Init_DECOMP

  !**********************************************************************
  ! Init_Prctl_MPI_TYPE
  !**********************************************************************
  subroutine Init_Prctl_MPI_TYPE()
    implicit none
    integer::ierror
    integer,dimension(3)::disp,blocklen,blocktype
  
    ! integer
    int_type = MPI_INTEGER
    call MPI_TYPE_SIZE(int_type,int_byte,ierror)

    ! real
    if(RK==4) then
      real_type = MPI_REAL
    else
      real_type = MPI_DOUBLE_PRECISION
    endif
    call MPI_TYPE_SIZE(real_type,  real_byte,  ierror)

    ! real3 type
    blocklen(1:3)=1
    blocktype(1:3)=real_type
    disp(1)=0
    disp(2)=disp(1)+real_byte
    disp(3)=disp(2)+real_byte
    call MPI_TYPE_STRUCT(3,blocklen(1:3),disp(1:3),blocktype(1:3),real3_type,ierror)
    call MPI_TYPE_COMMIT(real3_type,ierror)
    call MPI_TYPE_SIZE(real3_type, real3_byte, ierror)
  end subroutine Init_Prctl_MPI_TYPE
  
end module ATP_decomp_2d
