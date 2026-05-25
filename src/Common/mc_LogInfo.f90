#include "definitions_inc.f90"
module mc_LogInfo
  implicit none
  private
  integer,parameter,public:: ErrT_NoError = 0
  integer,parameter,public:: ErrT_Abort   = 1
  integer,parameter,public:: ErrT_Pass    = 2
    
  integer::g_ET_LastReportedError
  character(len=:),allocatable::g_Ch_LastReportedError
  character(len=*),dimension(10),parameter::bullet= [character(len=2)::  &
    ">>", ">", "++", "+", "--", "-", "**", "*", "==", "="]
  type LogType
    integer:: iUnit=100
    integer:: rprt_lvl_file=2
    integer:: rprt_lvl_cmdw=3
    character(len=:),allocatable::chLogFileName
  contains
    procedure:: InitLog   => LI_InitLog
    procedure:: OpenFile  => LI_OpenFile
    procedure:: CloseFile => LI_CloseFile
    procedure:: CreateFile=> LI_CreateFile
    procedure:: LI_OutInfo
    procedure:: LI_OutInfo2
    generic:: OutInfo         => LI_OutInfo, LI_OutInfo2
    procedure:: CheckForError => LI_CheckForError
  end type LogType
  public::LogType
contains

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! LI_InitLog
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine LI_InitLog(this,Dir_Res,RunName,file_lvl,cmdw_lvl); implicit none
    class(LogType)::this
    integer,intent(in)::file_lvl,cmdw_lvl
    character(len=*),intent(in)::Dir_Res,RunName
    
    this%rprt_lvl_file = file_lvl
    this%rprt_lvl_cmdw = cmdw_lvl
    this%chLogFileName = trim(adjustl(Dir_Res))//trim(adjustl(RunName))//".log"
  end subroutine LI_InitLog

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! LI_CreateFile
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine LI_CreateFile(this,RunName); implicit none
    class(LogType)::this
    character(len=*),intent(in)::RunName

    ! locals
    integer::ierr
    open(newunit=this%iUnit, file=this%chLogFileName, status='replace', IOSTAT=ierr)
    write(this%iUnit,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
    write(this%iUnit,*) "Log file for run: "//trim(adjustl(RunName))
    write(this%iUnit,*) "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
    write(this%iUnit,*)
    close(this%iUnit,IOSTAT=ierr)
  end subroutine LI_CreateFile
  
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! LI_OpenFile
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine LI_OpenFile(this); implicit none
    class(LogType)::this

    ! locals
    integer::ierr
    open(newunit=this%iUnit,file=this%chLogFileName,status='old',IOSTAT=ierr)
  end subroutine LI_OpenFile

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! LI_CloseFile
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine LI_CloseFile(this); implicit none
    class(LogType)::this

    ! locals
    integer::ierr
    close(this%iUnit,IOSTAT=ierr)
  end subroutine LI_CloseFile

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! LI_OutInfo
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine LI_OutInfo(this, chInfo, lvl, no_bull); implicit none
    class(LogType)::this
    integer,intent(in)::lvl
    character(len=*),intent(in)::chInfo
    logical,optional,intent(in)::no_bull

    ! locals
    logical::l_no
    character(len=3):: lvl_str
    character(len=:),allocatable::ch100,ch101
    
    l_no = .false.
    if(present(no_bull))l_no=no_bull
    write(lvl_str,'(I3)') 2*(lvl-1)
    
    if( lvl > 1 ) then   
      ch100 = "("//lvl_str//"x, A2, x, A)" 
      ch101 = "("//lvl_str//"x, A)"
    else
      ch100 = "(x, A2, x, A)"
      ch101 = "(x, A)"
    endif
    
    if( lvl <= this%rprt_lvl_file ) then
      if(lvl == 1) write(this%iUnit,*)
      if(l_no)then
        write(this%iUnit,ch101) trim(adjustl(chInfo)) 
      else
        write(this%iUnit,ch100) bullet(lvl), trim(adjustl(chInfo))   
      endif
    endif
    
    if( lvl <= this%rprt_lvl_cmdw )then
      if(lvl == 1 )write(*,*)
      if(l_no)then
        write(*,ch101) trim(adjustl(chInfo)) 
      else
        write(*,ch100) bullet(lvl), trim(adjustl(chInfo)) 
      endif
    endif
  end subroutine LI_OutInfo

  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! LI_OutInfo2
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine LI_OutInfo2(this, chInfo, chInfo2, lvl, no_bull); implicit none
    logical::l_no
    class(LogType)::this
    integer,intent(in)::lvl
    logical,optional,intent(in)::no_bull
    character(len=*),intent(in)::chInfo,chInfo2
      
    l_no = .false.
    if( present(no_bull) ) l_no = no_bull
    call this%OutInfo(chInfo, lvl, l_no)
    call this%OutInfo(chInfo2,lvl, l_no)
  end subroutine LI_OutInfo2
            
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  ! LI_CheckForError
  !  checking the input error message (Err_type) and creating a  
  !  message in command  window and logfile, then stopping the 
  !  execution of the program if necessary
  !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^80^^^^^^^^90
  subroutine LI_CheckForError(this, Err_type, chMethod, chMessage); implicit none
    class(LogType)::this
    integer,intent(in)::Err_type
    character(len=*),intent(in)::chMethod,chMessage
    
    select case(Err_type)
    case(ErrT_NoError)
      ! no error occurred, the program will continue its normal execution
      g_ET_LastReportedError = Err_type
      g_ch_LastReportederror = "No error occurred in: "//trim(adjustl(chMethod))//". Message: "//trim(adjustl(chMessage))
      return
        
    case(ErrT_Abort)
      ! a severe error occurred and program should be aborted 
      call this%OutInfo( "A severe error occurred in program", 1 )
      call this%OutInfo( "A Error occur in: "//trim(adjustl(chMethod)), "Error message is: "//trim(adjustl(chMessage)) , 2)
      g_ET_LastReportedError = Err_type
      g_ch_LastReportederror = "A severe error occurred in: "//trim(adjustl(chMethod))//". Message: "//trim(adjustl(chMessage))
      stop
        
    case(ErrT_Pass)
      ! a warning occurred in the program, a message will appear on the screen and 
      ! a message will be sent to log file but program continues running        
      call this%OutInfo( "A warning is reported in program", 1 )
      call this%OutInfo( "A warning occur in: "//trim(adjustl(chMethod)), "Warning message is: "//trim(adjustl(chMessage)), 2 )
      g_ET_LastReportedError = Err_type
      g_ch_LastReportederror = "A warning occurred in: "//trim(adjustl(chMethod))//". Message: "//trim(adjustl(chMessage))
      return
    end select
  end subroutine LI_CheckForError
  
end module mc_LogInfo
