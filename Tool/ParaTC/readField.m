clear;clc;
FileDir = './';          % Directory 
FindVar='uxFlua_';       % Variable name, ux/uy/uz/pr
FieldIndex=1;            % Field index. e.g. '1' means the first output field.
nxc=1680;                % Grid number in x-dir 
nyc=144;                 % Grid number in y-dir
nzc=1050;                % Grid number in z-dir
nxDomain=[1,1,nxc];      % Here [a:b:c] means reading from ix=a to ix=c, with step=b. 
nyDomain=[1,1,nyc];      % Here [a:b:c] means reading from iy=a to iy=c, with step=b. 
nzDomain=[1,1,nzc];      % Here [a:b:c] means reading from iz=a to iz=c, with step=b. 
real_prec='real*4';      % Binary precision
uTau=0.041958;
%% Read flow field, and then store in 'FieldOut'
dir_output=dir(fullfile(FileDir,['*',FindVar,'*']));
file_names={dir_output.name};
file_num=length(file_names);
if(FieldIndex>file_num)
  error('FieldIndex>file_num, Wrong !');
end
FileName=[FileDir, cell2mat(file_names(FieldIndex))];
FieldOut=readVisuFileFun(FileName,real_prec,nxc,nyc,nzc,nxDomain,nyDomain,nzDomain);

for k=1:nyc
  k,uxRms=sqrt(mean(mean(FieldOut(:,k,:).^2)))
  FieldOut(:,k,:)=FieldOut(:,k,:)/uxRms;
end
%FieldOut=FieldOut/uTau;

fid=fopen('uxRms','w');
fwrite(fid,FieldOut,'real*4');
fclose(fid);