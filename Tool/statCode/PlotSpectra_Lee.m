%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to plot 1D-Spectra data from Lee & Moser (2015)       %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2021-12-17                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
function PlotSpectra_Lee
  global PicYscale CBarXScale
  PicYscale=0.75;
  CBarXScale=0.6;
  SpecFileDir='../../LM_Channel_0550_1d_energy_spectra.h5';
  yly=h5read(SpecFileDir,'/Ly');height=0.5*yly;
  
  VarStr='uu';DirStr='x';
  [SpecData,yplus,lamda,Retau]=GetSpecData(SpecFileDir,VarStr,DirStr);
  x=lamda/height*Retau; y=yplus;
  PlotSpec(x,y,SpecData,VarStr,DirStr,'inner',3);
  x=lamda/height; y=yplus/Retau;
  PlotSpec(x,y,SpecData,VarStr,DirStr,'outer',4);
end

function [SpecData,yplus,lamda,Retau]=GetSpecData(FileStr,VarStr,DirStr)
  utau=h5read(FileStr,'/u_tau');
  SpecData=h5read(FileStr,['/E',VarStr,'_k',DirStr]);
  SpecData=SpecData';SpecData=SpecData(2:end,2:end);
  Ratio=1.0/utau^2;
  SpecData=Ratio*SpecData;
  
  waveNumber=h5read(FileStr,['/k',DirStr]);
  waveNumber=waveNumber(2:end);
  nk=size(SpecData,1);
  ny=size(SpecData,2);
  
  for k=1:nk
    Ratio=waveNumber(k);
    for m=1:ny
      SpecData(k,m)=SpecData(k,m)*Ratio;
    end
  end
  lamda=2*pi./waveNumber;
  yplus=h5read(FileStr,'/Y_plus');
  yplus=yplus';yplus=yplus(2:end);
  Retau=h5read(FileStr,'/Re_tau');
end

function PlotSpec(x,y,SpecData,VarStr,DirStr,PicType,nFig)
  global PicYscale CBarXScale
  figure(nFig);hold on;
  [X,Y]=meshgrid(x,y);
  if(strcmp(VarStr,'uv')==0) 
    %pcolor(X,Y, SpecData');
    contourf(X,Y, SpecData',10,'k','LineWidth',0.1);
  else
    %pcolor(X,Y,-SpecData');  
    contourf(X,Y,-SpecData',10,'k','LineWidth',0.1);
  end
  shading interp;box on;
  if(strcmp(PicType,'inner')==1)
    set(gca,'XScale','log');set(gca,'YScale','log');
  else
    set(gca,'XScale','log');set(gca,'YScale','lin');  
  end
  set(gca,'FontSize',12,'FontName','Times New Roman');
  set(gca,'linewidth',0.7)
  set(gca,'TickDir','out','TickLength',[0.015,0.015]);

  set(gca,'XLim',[min(x),max(x)]);
  set(gca,'xMinorTick','on');
  if(strcmp(PicType,'inner')==1)
    xlabel(['\fontsize{16}\it\lambda_',DirStr,'^+']);
  else
    xlabel(['\fontsize{16}\it\lambda_',DirStr,'/h']);  
  end

  if(strcmp(PicType,'inner')==1)
    set(gca,'YLim',[1,max(y)]);
    ylabel('\fontsize{16}\ity^+')
    set(gca,'YMinorTick','on');
  else
    set(gca,'YLim',[0,1]);
    ylabel('\fontsize{16}\ity/h')
    set(gca,'YMinorTick','off');
  end

  ch=colorbar;
  set(ch,'FontSize',12,'FontName','Times New Roman');
  set(ch,'linewidth',0.6);
  if(strcmp(VarStr,'uv')==0)
    title(ch,['\fontsize{12}\itk_',DirStr,'E_{',VarStr,'}/u_{\tau}^2']);
  else
    title(ch,['\fontsize{12}-\itk_',DirStr,'E_{',VarStr,'}/u_{\tau}^2']);  
  end
  set(ch,'TickDir','in','TickLength',0.02);
  color=[255 255 255;185 255 255;155 255 255;
          45 255 255; 35 255 220;135 255 120;
         236 255  19;255 201   0;255 134   0;
         255  67   0;255  0    0]/255;
  colormap(color); %get current colormap
  caxis([0,1.5]);set(ch,'YLim',[0,1.5]);
  if(strcmp(VarStr,'uv')==1)
   %caxis([0,-min(min(SpecData))]);   
   set(ch,'YLim',[0,-min(min(SpecData))]);
  end
  
  axpos=get(gca,'Position');
  axpos(4)=PicYscale*axpos(4);
  ch.Position(3)=CBarXScale*ch.Position(3);
  ch.Position(4)=axpos(4);
  set(gca,'Position',axpos);
  hold off;
end