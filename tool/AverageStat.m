%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to calculate the statistic results for wall-bounded   %
%   turbulent Closed-channel flows, calculated by ParaTC.                 %
% Channel3d can be freely downloaded from :                               %
%   https://github.com/GongZheng-Justin/ParaTC                            %
%                                                                         %
% There are 13 input parameters below                                     %
%   * xlx: Domain length in x-dir                                         %
%   * zlz: Domain lenght in z-dir                                         %
%   * nxc: Grid number in x-dir (nxc=nxp-1)                               %
%   * nzc: Gird number in z-dir (nzc=nzp-1)                               %
%   * xnu: Fluid kinematic viscosity                                      %
%   * iTimeSet:      Starting time for statistics calculation             %
%   * IsUxConst:     Does the mean streamwise velocity keep constant?     %
%   * BodyForceX:    If IsUxConst=0, use BodyForceX to calculate u_tau    %
%   * nEnergySpec1D: Number of 1D energy spectra                          %
%   * jForLCS:       Reference j-index for Linear coherent structure      %
%   * dir_statIn:    The folder to store the original/raw statistic data  %
%   * dir_statOut:   The folder to dump the final statistic results       %
%   * yMesh_str:     Ymesh file name                                      %
%                                                                         %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2022-01-09                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
clc;clear;

dir_statIn = './File/';
iTimeSet=0;

%% Calculate averaged data
dir_output=dir(fullfile(dir_statIn, 'vGrad0*') );
file_names={dir_output.name};
file_num=length(file_names);

for k=1:file_num
  datapatht=cell2mat(file_names(k));
  if(str2double(datapatht(6:15))>iTimeSet); 
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
  if(str2double(datapath(6:15))<=iTimeSet);continue;end
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

fid=fopen('AveStat.txt','w');
myformat=[repmat('%24.15E',1,real_num),'\n'];
for k=1:file_len
  fprintf(fid,myformat,dataE(k,:))
end
fclose(fid);
