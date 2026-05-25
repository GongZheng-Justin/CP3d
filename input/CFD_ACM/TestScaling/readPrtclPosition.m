clc;clear;
FileStr='./PartVisuForCFDACMScaling_DEM';

real_prec='real*8';
if(strcmp(real_prec,'real*8'))
  real_byte=8;
elseif(strcmp(real_prec,'real*4'))
  real_byte=4;  
else
  error('readVisuFile: real_prec Wrong');  
end

nPrtcl=1050624;
fid=fopen(FileStr,'r');
Position=fread(fid,[3,nPrtcl],real_prec);
fclose(fid);

DataOut=zeros(5,nPrtcl);
DataOut(1:3,:)=Position;
DataOut(4,:)=0.001;
DataOut(5,:)=1.0;

fid=fopen('SpheresCoord.dat','w');
fwrite(fid,DataOut,real_prec);
fclose(fid);