    ! Write 1-D Energy Spectra
    IF(clcSpectra1D) THEN
      infstime= 1.0_RK/real(nSpectime,RK)

      ! Spectra in x-dir
      allocate(EnergySpecXR(nxhp,x1size(2),NEnergySpec1D));EnergySpecXR=0.0_RK
      call MPI_REDUCE(EnergySpecX,EnergySpecXR,NEnergySpec1D*x1size(2)*nxhp,real_type,MPI_SUM,0,DECOMP_2D_COMM_ROW,ierror)
      call MPI_COMM_RANK(DECOMP_2D_COMM_ROW,nrankX,ierror)
      if(nrankX==0) then
        EnergySpecXR=EnergySpecXR*infstime
        disp_inc=int(mytype_bytes,8)*int(nyc,8)*int(nxhp,8)
        disp=int(mytype_bytes,8)*int(x1start(2)-1,8)*int(nxhp,8)
        
        write(filename,"(A,A,I10.10)") trim(ResultsDir),'SpecX',itime
        call MPI_FILE_OPEN(DECOMP_2D_COMM_COL,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,nUnit,ierror)
        call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
        call MPI_BARRIER(DECOMP_2D_COMM_COL,ierror)
        do kc=1,NEnergySpec1D
          call MPI_FILE_WRITE_AT_ALL(nUnit,disp,EnergySpecXR(:,:,kc),x1size(2)*nxhp,real_type,MPI_STATUS_IGNORE,ierror)
          disp=disp+disp_inc
        enddo
        call MPI_FILE_CLOSE(nUnit,ierror)
      endif
      deallocate(EnergySpecXR)
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)

      ! Spectra in z-dir
      allocate(EnergySpecZR(nzhp,z1size(2),NEnergySpec1D));EnergySpecZR=0.0_RK
      call MPI_REDUCE(EnergySpecZ,EnergySpecZR,NEnergySpec1D*z1size(2)*nzhp,real_type,MPI_SUM,0,DECOMP_2D_COMM_COL,ierror)
      call MPI_COMM_RANK(DECOMP_2D_COMM_COL,nrankZ,ierror)
      if(nrankZ==0) then
        EnergySpecZR=EnergySpecZR*infstime
        disp_inc=int(mytype_bytes,8)*int(nyc,8)*int(nzhp,8)
        disp=int(mytype_bytes,8)*int(z1start(2)-1,8)*int(nzhp,8)
        
        write(filename,"(A,A,I10.10)") trim(ResultsDir),'SpecZ',itime
        call MPI_FILE_OPEN(DECOMP_2D_COMM_ROW,filename,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL,nUnit,ierror)
        call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
        call MPI_BARRIER(DECOMP_2D_COMM_ROW,ierror)
        do kc=1,NEnergySpec1D
          call MPI_FILE_WRITE_AT_ALL(nUnit,disp,EnergySpecZR(:,:,kc),z1size(2)*nzhp,real_type,MPI_STATUS_IGNORE,ierror)
          disp=disp+disp_inc
        enddo
        call MPI_FILE_CLOSE(nUnit,ierror)
      endif
      deallocate(EnergySpecZR)
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    ENDIF

    ! Write 2-D Energy Spectra
    IF(clcSpectra2D) THEN
      infstime= 1.0_RK/real(nSpectime,RK)
      Block
      real(RK),dimension(:,:),allocatable::Ratio2D
      allocate(Ratio2D(decomp_xhzh%y2sz(1),decomp_xhzh%y2sz(3)))
      Ratio2D=4.0_RK
      if(decomp_xhzh%y2st(1)==1)    Ratio2D(1,:)=2.0_RK
      if(decomp_xhzh%y2st(3)==1)    Ratio2D(:,1)=2.0_RK
      if(decomp_xhzh%y2en(1)==nxhp) Ratio2D(decomp_xhzh%y2sz(1),:)=2.0_RK
      if(decomp_xhzh%y2en(3)==nzhp) Ratio2D(:,decomp_xhzh%y2sz(3))=2.0_RK
      if(decomp_xhzh%y2st(1)==1    .and. decomp_xhzh%y2st(3)==1   )  then
        Ratio2D(1,1)=1.0_RK
      endif
      if(decomp_xhzh%y2st(1)==1    .and. decomp_xhzh%y2en(3)==nzhp)  then
        Ratio2D(1,decomp_xhzh%y2sz(3))=1.0_RK
      endif
      if(decomp_xhzh%y2en(1)==nxhp .and. decomp_xhzh%y2en(3)==nzhp)  then
        Ratio2D(decomp_xhzh%y2sz(1),decomp_xhzh%y2sz(3))=1.0_RK
      endif
      if(decomp_xhzh%y2en(1)==nxhp .and. decomp_xhzh%y2st(3)==1   )  then
        Ratio2D(decomp_xhzh%y2sz(1),1)=1.0_RK
      endif
      Ratio2D=infstime*Ratio2D
      do kc=1,decomp_xhzh%y2sz(3)
        do ic=1,decomp_xhzh%y2sz(1)
          EnergySpec2D(ic,:,kc,:)=Ratio2D(ic,kc)*EnergySpec2D(ic,:,kc,:)
        enddo
      enddo   
      deallocate(Ratio2D)
      end Block
 
      Block
      integer(MPI_OFFSET_KIND)::disp
      integer::data_type,data_byte,newtype,nySpec2D
      integer,dimension(3)::sizes,subsizes,starts
#ifdef SAVE_SINGLE_Spec2D
      real(4),allocatable,dimension(:,:,:,:)::EnergySpecOut
#endif      
      if(FlowType==FT_CH) then
        nySpec2D=nyc/2
      else
        nySpec2D=nyc
      endif
#ifdef SAVE_SINGLE_Spec2D
      data_type= MPI_REAL
      allocate(EnergySpecOut(decomp_xhzh%y2sz(1),nySpec2D,decomp_xhzh%y2sz(3),NEnergySpec2D))
      do kt=1,NEnergySpec2D
        do kc=1,decomp_xhzh%y2sz(3)
          do jc=1,nySpec2D
            do ic=1,decomp_xhzh%y2sz(1)
              EnergySpecOut(ic,jc,kc,kt)=real(EnergySpec2D(ic,jc,kc,kt),4)
            enddo
          enddo
        enddo
      enddo
#else
      data_type= MPI_DOUBLE_PRECISION
#endif
      call MPI_TYPE_SIZE(data_type,data_byte,ierror)
      sizes(1)= nxhp
      sizes(2)= nySpec2D
      sizes(3)= nzhp 
      subsizes(1)=decomp_xhzh%y2sz(1)
      subsizes(2)=nySpec2D
      subsizes(3)=decomp_xhzh%y2sz(3)
      starts(1)=decomp_xhzh%y2st(1)-1
      starts(2)=0
      starts(3)=decomp_xhzh%y2st(3)-1
      call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,data_type,newtype,ierror)
      call MPI_TYPE_COMMIT(newtype,ierror)
      
      disp = 0_MPI_OFFSET_KIND
      write(filename,"(A,A,I10.10)") trim(ResultsDir),'Spec2D',itime
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, nUnit, ierror)
      call MPI_FILE_SET_SIZE(nUnit,0_MPI_OFFSET_KIND,ierror)  ! guarantee overwriting
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      do jc=1,NEnergySpec2D
        call MPI_FILE_SET_VIEW(nUnit,disp,data_type,newtype,'native',MPI_INFO_NULL,ierror)

#ifdef SAVE_SINGLE_Spec2D
        call MPI_FILE_WRITE_ALL(nUnit,EnergySpecOut(:,:,:,jc),subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
#else
        call MPI_FILE_WRITE_ALL(nUnit,EnergySpec2D(:,:,:,jc), subsizes(1)*subsizes(2)*subsizes(3),data_type,MPI_STATUS_IGNORE,ierror)
#endif
        disp=disp+ int(sizes(1),8)*int(sizes(2),8)*int(sizes(3),8)*int(data_byte,8) ! Update displacement
      enddo
      call MPI_FILE_CLOSE(nUnit,ierror)
      call MPI_TYPE_FREE(newtype,ierror)
      End Block
    ENDIF
