%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to plot 2D-Spectra picture                            %
% Note:                                                                   %
%   You need to run readSpec2D.m firstly.                                 %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2022-01-07                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&% 
clc;clear;

xlx=8*pi;
zlz=3*pi;
nxc=1680;
nzc=1050;
nyc=144;
IsOpenChannel=1;
real_prec='real*4';
height=1;
FileDir='../../StatOut/OC0550/';
VarStr='U';
nyIndex=40;
jForLCS=12;
utau=0.036245;

IsSetCaxis=0;
MaxCaxis=0.07;
PicYscale=0.75;
CBarXScale=0.6;
  
%% ========== Normally no need to change anything below ==========
if(IsOpenChannel==1) 
  nySpec2D=nyc;  
else
  nySpec2D=nyc/2;
end
nxh=nxc/2; nxhp=nxh+1;
nzh=nzc/2; nzhp=nzh+1;
if(strcmp(real_prec,'real*4')==1) 
  real_byte=4;
elseif(strcmp(real_prec,'real*8')==1)
  real_byte=8;
else
  error('readSpec2D: real_prec wrong');
end

% Numerator
FileStr=[FileDir,'Spec2D_L',VarStr];
fid=fopen(FileStr,'r');
fseek(fid, 0, 'eof');
totalbyte1=nxhp*nzhp*nySpec2D*real_byte;
totalbyte2=ftell(fid);
if(totalbyte1~=totalbyte2) 
  error('readSpec2D: File Byte Wrong');
end
fseek(fid,0,'bof');
SpecData=fread(fid,nySpec2D*nxhp*nzhp,real_prec);
SpecData=reshape(SpecData,[nxhp,nzhp,nySpec2D]);
SpecData=SpecData(:,:,nyIndex)/utau^2;
SpecDataX=sum(SpecData,2);
SpecDataZ=sum(SpecData,1);
SpecDataZ=SpecDataZ';
fclose(fid);