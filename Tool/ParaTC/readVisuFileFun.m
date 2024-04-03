function dataOut=readVisuFileFun(FileStr,real_prec,nxc,nyc,nzc,nxDomain,nyDomain,nzDomain)
  if(strcmp(real_prec,'real*8'))
    real_byte=8;
  elseif(strcmp(real_prec,'real*4'))
    real_byte=4;  
  else
    error('readVisuFile: real_prec Wrong');  
  end
  nxRead=nxDomain(1):nxDomain(2):nxDomain(3);
  nyRead=nyDomain(1):nyDomain(2):nyDomain(3);
  nzRead=nzDomain(1):nzDomain(2):nzDomain(3);
  
  fid=fopen(FileStr,'r');
  fseek(fid,0,'eof');
  totalByte1=nxc*nyc*nzc*real_byte;
  totalByte2=ftell(fid);
  if(totalByte1~=totalByte2)
    error('readVisuFile: File Byte Wrong');
  end
  dataOut=zeros(length(nxRead),length(nyRead),length(nzRead));
  for k=1:length(nzRead)
    kt=nzRead(k);
    offset=nxc*nyc*(kt-1)*real_byte;
    fseek(fid,offset,'bof');
    dataTemp=fread(fid,[nxc,nyc],real_prec);
    dataTemp=dataTemp(nxRead,nyRead);
    dataOut(:,:,k)=dataTemp;
  end
  fclose(fid);
end