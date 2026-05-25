function PlotLCS
  global PicYscale CBarXScale RefLine YminInner LCSFileName CValue
  SpecDataDir='../../StatOut/CC0550L/';
  height=1;
  PicYscale=0.75;
  CBarXScale=0.6;
  YminInner=4;
  RefLine=[10^3 10^4 10^2 10^3];
  CValue=[0.01,0.02,0.05,0.1:0.1:1.0];
  %CValue=[0.005 0.5];
  
  VarStr='uu';DirStr='x';
  LCSFileName=['LCS_',VarStr,'_',DirStr];
  LCSFileDir=[SpecDataDir,LCSFileName,'.txt'];
  SpecData=load(LCSFileDir);
  Retau=SpecData(1,1);
  yplus=SpecData(1,2:end);
  waveNumber=SpecData(3:end-1,1);
  lamda=2*pi./waveNumber;
  SpecData=SpecData(3:end-1,2:end);
  x=lamda/height*Retau; y=yplus;
  PlotSpecLCS(x,y,SpecData,DirStr,'inner',1);
  x=lamda/height; y=yplus/Retau;
  PlotSpecLCS(x,y,SpecData,DirStr,'outer',2);
end

function PlotSpecLCS(x,y,SpecData,DirStr,PicType,nFig)
  global PicYscale CBarXScale RefLine YminInner LCSFileName CValue
  figure(nFig);hold on;
  [X,Y]=meshgrid(x,y);
  
  SpecDataP=SpecData;
  %SpecDataP=filter2d(SpecData,1);
  %pcolor(X,Y, SpecData');
  contourf(X,Y, SpecDataP',CValue,...
      'LineColor',[0.5,0.5,0.5],'ShowText','off');
  shading interp;box on;
  if(strcmp(PicType,'inner')==1)
    set(gca,'XScale','log');set(gca,'YScale','log');
    plot(RefLine(1:2),RefLine(3:4),'k','linewidth',1.5);
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
    set(gca,'YLim',[YminInner,max(y)]);
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
  title(ch,'\fontsize{12}\it\gamma_{uu}^2');
  %title(ch,['\fontsize{12}\it\theta']);
  set(ch,'TickDir','in','TickLength',0.02);
  set(ch,'YLim',[0,1]);
    color=[255 255 255;
         185 255 255; 
         155 255 255;
          45 255 255;
          35 255 220;
         135 255 120;
         236 255  19;
         255 201   0;
         255 134   0;
         255  67   0;
         255  0    0]/255;
  colormap(color); %get current colormap
  %caxis([0 35]);
  
  axpos=get(gca,'Position');
  axpos(4)=PicYscale*axpos(4);
  ch.Position(3)=CBarXScale*ch.Position(3);
  ch.Position(4)=axpos(4);
  set(gca,'Position',axpos);
  print(gcf,'-dtiff','-r240',['../../',LCSFileName,PicType,'.tiff'])
  hold off;
end