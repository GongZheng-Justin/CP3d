module PPD_TypeDef
  implicit none

  integer,parameter::RKP=4
  integer,parameter::RK =kind(0.0D0)
  real(RK),parameter::Pi = 3.141592653589793238462643383279502884_RK    ! Pi constant

  type DumpPrtclVar
    integer:: itime
    integer:: id
    integer:: pType
    integer:: CntctFlag
    real(RKP):: Pos(3)
    real(RKP):: linVel(3)
    real(RKP):: FpForce(3)
    real(RKP):: CntctForce(3)
    real(RKP):: RotVel(3)
    real(RKP):: FpTorque(3)
    real(RKP):: CntctTorque(3)
  end type DumpPrtclVar

  type DumpPrtclVarOut
    integer:: itime
    integer:: pType
    integer:: CntctFlag
    real(RKP):: Pos(3)
    real(RKP):: linVel(3)
    real(RKP):: FpForce(3)
    real(RKP):: CntctForce(3)
    real(RKP):: RotVel(3)
    real(RKP):: FpTorque(3)
    real(RKP):: CntctTorque(3)
  end type DumpPrtclVarOut

  interface assignment(=)
    module procedure DumpPrtcl_Var_Out
  end interface

contains

  !******************************************************************
  ! DumpPrtcl_Var_Out
  !******************************************************************
  subroutine DumpPrtcl_Var_Out(Dump2,Dump1)
    implicit none
    type(DumpPrtclVar),intent(in)::Dump1
    type(DumpPrtclVarOut),intent(out)::Dump2

    Dump2%itime        = Dump1%itime
    Dump2%pType        = Dump1%pType
    Dump2%CntctFlag    = Dump1%CntctFlag
    Dump2%Pos          = Dump1%Pos
    Dump2%linVel       = Dump1%linVel
    Dump2%FpForce      = Dump1%FpForce
    Dump2%CntctForce   = Dump1%CntctForce
    Dump2%RotVel       = Dump1%RotVel
    Dump2%FpTorque     = Dump1%FpTorque
    Dump2%CntctTorque  = Dump1%CntctTorque
  end subroutine 

end module PPD_TypeDef
