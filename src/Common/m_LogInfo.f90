module m_LogInfo
  implicit none
  private
  integer,parameter,public:: ErrT_NoError = 0
  integer,parameter,public:: ErrT_Abort   = 1
  integer,parameter,public:: ErrT_Pass    = 2
    
  integer::g_ET_LastReportedError
  character(512)::g_Ch_LastReportedError
  character(2),dimension(10),parameter::bullet= (/">>"," >","++"," +", &
                                        "--"," -","**"," *","==" ," ="/)   
  type LogType
    integer:: nUnit=100
    integer:: rprt_lvl_file=2
    integer:: rprt_lvl_cmdw=3
    character(128)::chLogFileName = "file.log"
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

  !*****************************************************************
  ! LI_InitLog
  !*****************************************************************
  subroutine LI_InitLog(this,Dir_Res,RunName,file_lvl,cmdw_lvl)
    implicit none
    class(LogType)::this
    integer,intent(in)::file_lvl,cmdw_lvl
    character(*),intent(in)::Dir_Res,RunName
    
    this%rprt_lvl_file = file_lvl
    this%rprt_lvl_cmdw = cmdw_lvl
    write(this%chLogFileName,"(A)") trim(Dir_Res)//trim(RunName)//".log"
  end subroutine LI_InitLog

  !*****************************************************************
  ! LI_CreateFile
  !*****************************************************************
  subroutine LI_CreateFile(this,RunName)
    implicit none
    class(LogType)::this
    character(*),intent(in)::RunName

    ! locals
    integer::ierror
    open(newunit=this%nUnit, file=this%chLogFileName, status='replace', IOSTAT=ierror)
    write(this%nUnit,*) "**************************************************************************"
    write(this%nUnit,*) "Log file for run: "//trim(RunName)
    write(this%nUnit,*) "**************************************************************************"
    write(this%nUnit,*)
    close(this%nUnit,IOSTAT=ierror)
  end subroutine LI_CreateFile
  
  !*****************************************************************
  ! LI_OpenFile
  !*****************************************************************
  subroutine LI_OpenFile(this)
    implicit none
    class(LogType)::this

    ! locals
    integer::ierror
    open(newunit=this%nUnit,file=this%chLogFileName,status='old',IOSTAT=ierror)
  end subroutine LI_OpenFile

  !*****************************************************************
  ! LI_CloseFile
  !*****************************************************************
  subroutine LI_CloseFile(this)
    implicit none
    class(LogType)::this

    ! locals
    integer::ierror
    close(this%nUnit,IOSTAT=ierror)
  end subroutine LI_CloseFile

  !*****************************************************************
  ! LI_OutInfo
  !*****************************************************************
  subroutine LI_OutInfo( this, chInfo, lvl , no_bull )
    implicit none
    class(LogType)::this
    integer,intent(in)::lvl
    character(*),intent(in)::chInfo
    logical,optional,intent(in)::no_bull

    ! locals
    logical::l_no
    character(64)::ch100,ch101
    
    l_no = .false.
    if( present(no_bull))l_no=no_bull
    if( lvl > 1 ) then    
      write(ch100,"(A,I3,A)")"(",2*(lvl-1), "x , A2, x ,A)" 
      write(ch101,"(A,I3,A)")"(",2*(lvl-1), "x ,A)"
    else
      write(ch100,"(A)") "(x , A2, x ,A)"
      write(ch101,"(A)") "(x ,A)"
    endif
    
    if( lvl <= this%rprt_lvl_file )then
      if(lvl == 1 )write(this%nUnit,*)
      if(l_no)then
        write(this%nUnit,ch101) trim(chInfo)    
      else
        write(this%nUnit,ch100) bullet(lvl), trim(chInfo)    
      endif
    endif
    
    if( lvl <= this%rprt_lvl_cmdw )then
      if(lvl == 1 )write(*,*)
      if(l_no)then
        write(*,ch101) trim(chInfo)
      else
        write(*,ch100) bullet(lvl), trim(chInfo)
      endif
    endif
  end subroutine LI_OutInfo

  !*****************************************************************
  ! LI_OutInfo2
  !*****************************************************************
  subroutine LI_OutInfo2( this, chInfo , chInfo2, lvl , no_bull )
    implicit none
    logical::l_no
    class(LogType)::this
    integer,intent(in)::lvl
    logical,optional,intent(in)::no_bull
    character(*),intent(in)::chInfo,chInfo2
      
    l_no = .false.
    if( present(no_bull) ) l_no = no_bull
    call this%OutInfo(chInfo, lvl, l_no)
    call this%OutInfo(chInfo2,lvl, l_no)
  end subroutine LI_OutInfo2
            
  !*****************************************************************
  ! LI_CheckForError
  !  checking the input error message (Err_type) and creating a  
  !  message in command  window and logfile, then stopping the 
  !  execution of the program if necessary
  !*****************************************************************
  subroutine LI_CheckForError( this, Err_type, chMethod, chMessage )
    implicit none
    class(LogType)::this
    integer,intent(in)::Err_type
    character(*),intent(in)::chMethod,chMessage
    
    select case(Err_type)
    case(ErrT_NoError)
      ! no error occurred, the program will continue its normal execution
      g_ET_LastReportedError = Err_type
      write(g_ch_LastReportederror,"(A)") "No error occurred in: "//trim(chMethod)//". Message: "//trim(chMessage)
      return
        
    case(ErrT_Abort)
      ! a severe error occurred and program should be aborted 
      call this%OutInfo( "A severe error occurred in program", 1 )
      call this%OutInfo( "A Error occur in: "//trim(chMethod), "Error message is: "//trim(chMessage) , 2 )
      g_ET_LastReportedError = Err_type
      write(g_ch_LastReportederror,"(A)") "A severe error occurred in: "//trim(chMethod)//". Message: "//trim(chMessage)
      stop
        
    case(ErrT_Pass)
      ! a warning occurred in the program, a message will appear on the screen and 
      ! a message will be sent to log file but program continues running        
      call this%OutInfo( "A warning is reported in program", 1 )
      call this%OutInfo( "A warning occur in: "//trim(chMethod), "Warning message is: "//trim(chMessage) , 2 )
      g_ET_LastReportedError = Err_type
      write(g_ch_LastReportederror,"(A)") "A warning occurred in: "//trim(chMethod)//". Message: "//trim(chMessage)
      return
    end select
  end subroutine LI_CheckForError
  
end module m_LogInfo
