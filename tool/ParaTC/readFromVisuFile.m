%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to read data from visu files                          %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2021-12-17                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
clear;clc;
FileDir = '../CFD/Results/';         % Directory 
FindVar='uy_';           % Variable name, ux/uy/uz/pr
FieldIndex=1;            % Field index. e.g. '1' means the first output field.
nxc=1680;                % Grid number in x-dir 
nyc=144;                 % Grid number in y-dir
nzc=1050;                % Grid number in z-dir
nxDomain=[1,1,nxc];      % Here [a:b:c] means reading from ix=a to ix=c, with step=b. 
nyDomain=[1,1,nyc];      % Here [a:b:c] means reading from iy=a to iy=c, with step=b. 
nzDomain=[1,1,nzc];      % Here [a:b:c] means reading from iz=a to iz=c, with step=b. 
real_prec='real*8';      % Binary precision
MeshName='../CFD/Results/yMeshForCha550_04.txt'; % yMesh file

%% Read flow field, and then store in 'FieldOut'
dir_output=dir(fullfile(FileDir,['*',FindVar,'*']));
file_names={dir_output.name};
file_num=length(file_names);
if(FieldIndex>file_num)
  error('FieldIndex>file_num, Wrong !');
end
FileName=[FileDir, cell2mat(file_names(FieldIndex))];
FieldOut=readVisuFileFun(FileName,real_prec,nxc,nyc,nzc,nxDomain,nyDomain,nzDomain);

%% yOut
nyRead=nyDomain(1):nyDomain(2):nyDomain(3);
yOut=zeros(length(nyRead),1);
yp=load(MeshName);yp=yp(:,2);
if(strcmp(FindVar,'uy_')) 
  for k=1:length(nyRead)
    kt=nyRead(k);
    yOut(k)=yp(kt);
  end
else
  for k=1:length(nyRead)
    kt=nyRead(k);
    yOut(k)=0.5*(yp(kt)+yp(kt+1));
  end 
end