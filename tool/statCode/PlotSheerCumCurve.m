%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to plot 1D-Spectra picture                            %
% Note:                                                                   %
%   You should run 'readStat_CH.m'/'readStat_HC.m' firstly                %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2022-01-02                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
  height=1;
  PicYscale=0.75;
  CBarXScale=0.6;
  
  VarStr='uv';DirStr='z';
  SpecDataDir='../../StatOut/OC0550/';
  SpecFileDir=[SpecDataDir,'spec',DirStr,'_',VarStr,'.txt'];
  
  SpecDataT=load(SpecFileDir);
  nk=size(SpecDataT,1)-3;
  ny=size(SpecDataT,2)-1;
  Retau=SpecDataT(1,1);
  yplus=SpecDataT(1,2:end);
  waveNumber=SpecDataT(3:end-1,1);
  lamda=2*pi./waveNumber;

  SpecData=zeros(nk,1);
  ybar=yplus/Retau;
  RatioY=-3.0*(1-ybar);
  ybarP=zeros(ny+1,1);
  ybarP(1)=0.0;
  ybarP(end)=1;
  for m=1:ny-1
    ybarP(m+1)=2*ybar(m)-ybarP(m);
  end
  for m=1:ny
    ybar(m)=ybarP(m+1)-ybarP(m); 
  end
  RatioY=RatioY.*ybar;
  SpecDataT=SpecDataT(3:end-1,2:end);
  dWave=1.0/(waveNumber(2)-waveNumber(1));
  for k=1:nk
    Ratio=waveNumber(k)*dWave;
    for m=1:ny
      SpecData(k)=SpecData(k)+SpecDataT(k,m)*Ratio*RatioY(m);
    end
  end
  x=lamda/height; y=SpecData;
  semilogx(x,y)