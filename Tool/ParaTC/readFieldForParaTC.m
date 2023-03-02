clear;clc;
nxc=1536;                % Grid number in x-dir 
nyc=384;                 % Grid number in y-dir
nzc=1024;                % Grid number in z-dir
nxDomain=[1, 1,nxc];      % Here [a:b:c] means reading from ix=a to ix=c, with step=b. 
nyDomain=[1,1,nyc];      % Here [a:b:c] means reading from iy=a to iy=c, with step=b. 
nzDomain=[nzc,1,nzc];      % Here [a:b:c] means reading from iz=a to iz=c, with step=b. 
real_prec='real*8';      % Binary precision
FileName='RestartForCha550_040000056000';

%% Read flow field, and then store in 'FieldOut'
FieldOut=readVisuFileFun2(FileName,real_prec,nxc,nyc,nzc,nxDomain,nyDomain,nzDomain);

fid=fopen('Ri000','w');
fwrite(fid,FieldOut/0.035983895051567,'real*4');
fclose(fid);