%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to calculate the statistic results for 1-way Particle-%
%   laden Closed-channel flows, calculated by CP3d.                       %
% Channel3d can be freely downloaded from :                               %
%   https://github.com/GongZheng-Justin/CP3d                              %
%                                                                         %
% There are 8 input parameters below                                      %
%   * uTau:   Firction velocity                                           %
%   * height: half-channel height of closed-channel                       %
%   * xnu:    Fluid kinematic viscosity                                   %
%   * iTimeSet:      Starting time for statistics calculation             %
%   * iTimeEnd:      Ending time for statistics calculation               %
%   * dir_statIn:    The folder to store the original/raw statistic data  %
%   * dir_statOut:   The folder to dump the final statistic results       %
%   * file_filter:   filter for the statistical file                      %
%                                                                         %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2022-07-01                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

clear;clc;

uTau=0.11775;
height=0.02;
xnu=1.570E-5;
iTimeSet=100000; %6000;
iTimeEnd=200000;

dir_statIn = '../../LPT/Results4th/';
dir_statOut='../../StatOut/LPTOneWay4th/';
file_filter = 'pstats*_set0003';

%% ========== Normally no need to change anything below ==========
if ( exist(dir_statOut,'dir') ==false )
  mkdir(dir_statOut);
  fprintf( '%s\n\n',['crete directory: ',dir_statOut,'  sucessfully'] );
end

%% Calculate averaged data
dir_output=dir(fullfile(dir_statIn,file_filter));
file_names={dir_output.name};
file_num=length(file_names);

for k=1:file_num
  datapatht=cell2mat(file_names(k));
  if(str2double(datapatht(7:16))>iTimeSet); 
    datapath = [dir_statIn,datapatht];
    break;
  end
end

fid=fopen(datapath,'r');
line1 = 0;

while (feof(fid)==0)
  str=strtrim( fgets(fid) );
  dlt=sscanf(str,'%f');
  line1=line1 +1;
  if(isempty(dlt)==0 && length(dlt)>2)
    break;
  end
end
real_num=length(dlt);
frewind(fid);
line2 = 0;

while (feof(fid)==0)
  fgets(fid);
  line2= line2 +1;
end
fclose(fid);
file_len=line2-line1+1;

file_ave=0;
prgrad=zeros(file_ave,1);
dataE=zeros(file_len,real_num);
data_emerget=zeros(file_len,real_num);
for k=1:file_num
  datapath=cell2mat(file_names(k));
  if(str2double(datapath(7:16))<=iTimeSet || str2double(datapath(7:16))>iTimeEnd);continue;end
  file_ave=file_ave+1;  
  datapath = [dir_statIn, cell2mat(file_names(k)) ];
  fid=fopen(datapath,'r');
  for kt=1:line1-1
    str=strtrim( fgets(fid) );
    prgradt=sscanf(str,'time averaged pressure gradient is:%f');
    if(isempty(prgradt)==0)
      prgrad(file_ave)=prgradt;
    end      
  end    
  
  for kt=line1:line2
    idl=kt-line1+1;
    str=strtrim( fgets(fid) );
    dlt=sscanf(str,'%f');
    data_emerget(idl,:)=dlt';
  end
  fclose(fid);
  dataE=dataE+ data_emerget;
  disp( ['read:   ',datapath,'  sucessfully'] );
end
dataE=dataE/file_ave;

file_lenh=file_len/2;
dataPrf=zeros(file_lenh,8);
dataPrf(1:file_lenh,1)=dataE(1:file_lenh,1)/height;
dataPrf(1:file_lenh,2)=dataE(1:file_lenh,1)*uTau/xnu;
dataPrf(1:file_lenh,3)=0.5*(dataE(1:file_lenh,3)+dataE(end:-1:file_lenh+1,3))/uTau;

vectemp=sqrt(dataE(:,6)-dataE(:,3).*dataE(:,3));
dataPrf(1:file_lenh,4)=0.5*(vectemp(1:file_lenh)+vectemp(end:-1:file_lenh+1))/uTau;

vectemp=sqrt(dataE(:,7)-dataE(:,4).*dataE(:,4));
dataPrf(1:file_lenh,5)=0.5*(vectemp(1:file_lenh)+vectemp(end:-1:file_lenh+1))/uTau;

vectemp=sqrt(dataE(:,8)-dataE(:,5).*dataE(:,5));
dataPrf(1:file_lenh,6)=0.5*(vectemp(1:file_lenh)+vectemp(end:-1:file_lenh+1))/uTau;

vectemp=(dataE(:,9)-dataE(:,3).*dataE(:,4));
dataPrf(1:file_lenh,7)=0.5*(vectemp(1:file_lenh)-vectemp(end:-1:file_lenh+1))/(uTau*uTau);

dataPrf(1:file_lenh,8)=0.5*(dataE(1:file_lenh,10)+dataE(end:-1:file_lenh+1,10))/uTau;
dataPrf(1:file_lenh,8)=dataPrf(1:file_lenh,8)-dataPrf(1:file_lenh,3);

myformat=[repmat('%24.15E',1,8),'\n'];
FileName=['PrtclProfile_',file_filter(end-6:end),'.txt'];
fid=fopen([dir_statOut,FileName],'wt');
strout=' ybar yplus up upRms vpRms wpRms up''vp'' uf-up';
fprintf(fid,'%s\n',strout);
for k=1:file_lenh
  fprintf(fid,myformat,dataPrf(k,:));
end
fclose(fid);
