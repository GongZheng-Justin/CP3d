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
function PlotSheerCum
  global PicYscale CBarXScale
  height=1;
  SpecDataDir='../../StatOut/CC0550L/';
  PicYscale=0.75;
  CBarXScale=0.6;
  
  VarStr='uv';DirStr='x';
  SpecFileDir=[SpecDataDir,'spec',DirStr,'_',VarStr,'.txt'];
  [SpecData,yplus,lamda,Retau]=GetSpecData(SpecFileDir);
  x=lamda/height*Retau; y=yplus;
  PlotSpec(x,y,SpecData,VarStr,DirStr,'inner',1);
  x=lamda/height; y=yplus/Retau;
  PlotSpec(x,y,SpecData,VarStr,DirStr,'outer',2);
end

function [SpecData,yplus,lamda,Retau]=GetSpecData(SpecFileDir)
  SpecData=load(SpecFileDir);
  nk=size(SpecData,1)-3;
  ny=size(SpecData,2)-1;
  Retau=SpecData(1,1);
  yplus=SpecData(1,2:end);
  waveNumber=SpecData(3:end-1,1);
  lamda=2*pi./waveNumber;

  ybar=yplus/Retau;
  RatioY=3.0*ybar.*(1-ybar);
  SpecData=SpecData(3:end-1,2:end);
  dWave=1.0/(waveNumber(2)-waveNumber(1));
  for k=1:nk
    Ratio=waveNumber(k)*dWave;
    for m=1:ny
      SpecData(k,m)=SpecData(k,m)*Ratio*RatioY(m);
    end
  end
end

function PlotSpec(x,y,SpecData,VarStr,DirStr,PicType,nFig)
  global PicYscale CBarXScale
  figure(nFig);hold on;
  [X,Y]=meshgrid(x,y);
  
  %SpecDataP=SpecData;
  SpecDataP=filter2d(SpecData,1);
  %pcolor(X,Y,-SpecData');
  contourf(X,Y,-SpecDataP',20,'k','LineWidth',0.1,'LineColor','none');
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
    set(gca,'YLim',[5,max(y)]);
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
    title(ch,['\fontsize{12}\itk_',DirStr,'E_{',VarStr,'}^+']);
  else
    title(ch,['\fontsize{12}-\itk_',DirStr,'E_{',VarStr,'}^+']);  
  end
  set(ch,'TickDir','in','TickLength',0.02);
color=  [ 255 255 255;
    217 217	217;
189	189	189;
150	150	150;
99	99	99;
219	218	236;
188	189	220;
157	154	199;
117	107	176;
198	233	191;
159	218	152;
115	196	117;
48	165	87;
198	219	240;
157	202	223;
107	175	214;
50	130	189;
252	208	161;
252	175	107;
253	140	60;
228	85	15;
]/255;
  colormap(color); %get current colormap
  %caxis([0,0.5]);set(ch,'YLim',[0,0.5]);
  if(strcmp(VarStr,'uv')==1)
   set(ch,'YLim',[0,-min(min(SpecDataP))]);
  end
  
  axpos=get(gca,'Position');
  axpos(4)=PicYscale*axpos(4);
  ch.Position(3)=CBarXScale*ch.Position(3);
  ch.Position(4)=axpos(4);
  set(gca,'Position',axpos);
  hold off;
end