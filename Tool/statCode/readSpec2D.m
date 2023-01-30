%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to read and calculate 2D-Spectra data                 %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2022-01-07                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&% 
clc;clear;

nxc=9216;
nzc=1400;
nyc=144;
IsOpenChannel=1;
NEnergySpec2D=9;
real_prec='real*4';

iTimeSet=0;
VarStr={'UU','VV','WW','PP','UV','LU','LV','LW','LP','CC','LC'};

dir_statIn = '../../CFD/ResultsOC0550L/';
dir_statOut='../../StatOut/OC0550L/';

%% ========== Normally no need to change anything below ==========
if(IsOpenChannel==1) 
  nySpec2D=nyc;  
else
  nySpec2D=nyc/2;
end

if(strcmp(real_prec,'real*4')==1) 
  real_byte=4;
elseif(strcmp(real_prec,'real*8')==1)
  real_byte=8;
else
  error('readSpec2D: real_prec wrong');
end
if ( exist(dir_statOut,'dir') ==false )
  mkdir(dir_statOut);
  fprintf( '%s\n\n',['crete directory: ',dir_statOut,'  sucessfully'] );
end
nxh=nxc/2; nxhp=nxh+1;
nzh=nzc/2; nzhp=nzh+1;
dir_output=dir(fullfile(dir_statIn,'Spec2D*'));
file_names={dir_output.name};
file_num=length(file_names);
if(file_num>0)
  for k=1:file_num
    datapath = cell2mat(file_names(k));
    if(str2double(datapath(7:16))<=iTimeSet);continue;end;
    datapath = [dir_statIn,cell2mat(file_names(k))];
    fid=fopen(datapath,'r');
    fseek(fid, 0, 'eof');
    totalbyte1=NEnergySpec2D*nySpec2D*nxhp*nzhp*real_byte;
    totalbyte2=ftell(fid);
    if(totalbyte1~=totalbyte2) 
      error('readSpec2D: File Byte Wrong');
    end
    fclose(fid);      
  end
  for NE=1:NEnergySpec2D
    file_ave=0;
    Spectra2D=zeros(nxhp*nySpec2D*nzhp,1);
    offset=(NE-1)*nySpec2D*nxhp*nzhp*real_byte;
    for k=1:file_num
      datapath = cell2mat(file_names(k));
      if(str2double(datapath(7:16))<=iTimeSet);continue;end;
      datapath = [dir_statIn,cell2mat(file_names(k))];
      file_ave=file_ave+1;
      fid=fopen(datapath,'r');
      fseek(fid, offset, 'bof');
      SpecData=fread(fid,nySpec2D*nxhp*nzhp,real_prec);
      Spectra2D=Spectra2D+SpecData;
      fclose(fid);
      disp( [cell2mat(VarStr(NE)),' read:   ',datapath,'  sucessfully'] );
    end
    Spectra2D=Spectra2D/file_ave;
    writename=[dir_statOut,'Spec2D_',cell2mat(VarStr(NE))];
    fid=fopen(writename,'w');
    fseek(fid,0,'bof');
    SpecData=reshape(Spectra2D,[nxhp,nySpec2D,nzhp]);
    Spectra2D=zeros(nxhp,nzhp,nySpec2D);
    for k=1:nySpec2D
      Spectra2D(:,:,k)=SpecData(:,k,:);
    end
    fwrite(fid,Spectra2D,real_prec);
    fclose(fid);
  end
end