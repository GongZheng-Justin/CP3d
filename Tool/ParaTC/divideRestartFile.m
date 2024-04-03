clear;clc;
FileDir = './';         % Directory
ReadFileName='RestartForOCT0550L_020000150000';
nxc=9216;                % Grid number in x-dir 
nyc=144;                 % Grid number in y-dir
nzc=1400;                % Grid number in z-dir
real_prec='real*8';
varNumber=4;

%% Nothing need to be changed below
if(strcmp(real_prec,'real*8'))
  real_byte=8;
elseif(strcmp(real_prec,'real*4'))
  real_byte=4;  
else
  error('readVisuFile: real_prec Wrong');  
end

FileInVar=[FileDir,ReadFileName];
fidIn=fopen(FileInVar,'r');
fseek(fidIn,0,'eof');
SingleByte=nxc*nyc*nzc*real_byte;
totalByte1=varNumber*SingleByte;
totalByte2=ftell(fidIn);
if(totalByte1~=totalByte2)
  error('readRestaartFile: File Byte Wrong');
end

readNumMax=9999999;
for k=1:varNumber
  FileOutVar=sprintf('%s%s%2.2d',FileInVar,'_',k)
  
  nLeft=nxc*nyc*nzc;
  offWrite=0;
  offRead=(k-1)*SingleByte;
  fidOut=fopen(FileOutVar,'w');
  while(nLeft>0)
    readNum=min(readNumMax,nLeft);
    nLeft=nLeft-readNum
    
    fseek(fidIn,offRead,'bof');
    readVec=fread(fidIn,readNum,real_prec);
    offRead=offRead+readNum*real_byte;
    
    fseek(fidOut,offWrite,'bof');
    fwrite(fidOut,readVec,real_prec);
    offWrite=offWrite+readNum*real_byte;
  end
  fclose(fidOut);
end
fclose(fidIn);