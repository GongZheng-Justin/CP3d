%
% Zheng Gong, 2024-12-03
%
clc; clear;

fileName='SpheresCoord.Phi=0.5';

%=============================
size_one_p=40;
fid=fopen(fileName,'rb');
fseek(fid,0,'eof');
disp=ftell(fid);
if(mod(disp, size_one_p) ~=0 )
    error('disp wrong');
end
fseek(fid,0,'bof');
np=disp/size_one_p;

dataRead=fread(fid,[5,np],'real*8');
fclose(fid);

fid=fopen([fileName,'.txt'],'wt');

myFormat = '%24.15E  %24.15E  %24.15E  %24.15E  %d';
for k=1:np-1
    fprintf(fid, [myFormat,'\n'], dataRead(:,k));
end
k=np;
fprintf(fid, myFormat, dataRead(:,k));
fclose(fid);
