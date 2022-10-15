program GeomToVTK
  implicit none
  integer,parameter::RK = KIND(0.0D0)
  integer i, j, npoint,  nw, nUnitFile1,  nUnitFile2,myistat,face
  character(256) :: chFile1, chFile2, chline, chMsg,  Res_Dir,RunName
  type real3
    real(RK)::x
    real(RK)::y
    real(RK)::z
  end type real3
  type(real3)::point
    
    
  ! ============   change the following four lines if necessary ============
  Res_Dir = "../DEM/Results/"
  RunName="CollideWall"
  chFile1 = trim(Res_Dir)//"WallsFor"//trim(RunName)//".backup"
  chFile2 = trim(Res_Dir)//"WallsFor"//trim(RunName)//".vtk"
    
  nUnitFile1 = 142
  open( unit =  nUnitFile1 , file = chFile1, status='old',form='formatted',IOSTAT=myistat )
  if( myistat /=0 ) then
    STOP "opening backup Geometry file failed !!!"    
  endif

  npoint=0
  do 
    npoint = npoint + 1
    read( nUnitFile1 , "(A)" ,IOSTAT = myistat ) chline
    if ( myistat < 0)  exit
  enddo
  npoint = npoint - 1
   
  if( npoint == 0 .or. mod(npoint,4) /= 0  ) then
    STOP " number of point wrong!!!  "    
  endif
  rewind(nUnitFile1)
  nw=npoint/4
    
  nUnitFile2 = 143
  open( unit =  nUnitFile2 , file = chFile2, status='replace',form='formatted',IOSTAT=myistat )
  if( myistat /=0 ) then
    STOP "opening vtk file  failed !!!"    
  endif    
  chMsg = "Geometry for run name: " // trim(RunName)   
  write( nUnitFile2, "(A)")"# vtk DataFile Version 2.0"
  write( nUnitFile2, "(A)") trim( chMsg )
  write( nUnitFile2, "(A)") "ASCII"        
  write( nUnitFile2, "(A)" )"DATASET POLYDATA"
  write( nUnitFile2, "(A8,I20,A)" )"POINTS  ", npoint , "  float"    
  do i = 1, npoint
    read(nUnitFile1, *) point
    write(nUnitFile2, *) point   
  enddo    
  write( nUnitFile2, "(A8, I20, A1, I20)") "POLYGONS ", nw , " " , 5*nw    
  do i =1, nw
    face=4;             write( nUnitFile2, *) face
    face=(i-1)*4;       write( nUnitFile2, *) face
    face=(i-1)*4+1;     write( nUnitFile2, *) face
    face=(i-1)*4+2;     write( nUnitFile2, *) face
    face=(i-1)*4+3;     write( nUnitFile2, *) face
  enddo    
     
  close(nUnitFile1)
  close(nUnitFile2)
    
  print*,'GeomToVTK finished sucessfully ! '
    
end program GeomToVTK
