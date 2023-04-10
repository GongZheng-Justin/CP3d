%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to plot 1D-Spectra data from Lee & Moser (2015)       %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2021-12-17                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
  clc; clear;
  FileStr='./LM_Channel_0180_1d_energy_spectra.h5';
  yly=h5read(FileStr,'/Ly');height=0.5*yly;
  
  VarStr='uv';DirStr='z';
  utau=h5read(FileStr,'/u_tau');
  SpecData=h5read(FileStr,['/E',VarStr,'_k',DirStr]);
  SpecData=SpecData';SpecData=SpecData(2:end,2:end);
  Ratio=1.0/utau^2;
  SpecData=Ratio*SpecData;
  waveNumber=h5read(FileStr,['/k',DirStr]);
  
  lamda=2*pi./waveNumber;
  yplus=h5read(FileStr,'/Y_plus');
  yplus=yplus';yplus=yplus(2:end);
  Retau=h5read(FileStr,'/Re_tau');
  ybar=yplus/Retau;