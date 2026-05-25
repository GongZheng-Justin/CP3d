%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%
% This file is used to plot Lumley Triangle picture                       %
% Note:                                                                   %
%   You should run 'readStat_CH.m'/'readStat_HC.m' firstly                %
% Author:                                                                 %
%   Zheng Gong, Department of Hydraulic Engineering, Tsinghua University  %
% E-mail:                                                                 %
%   gongzheng_justin@outlook.com                                          %
% Last modification date:                                                 %
%   2021-12-17                                                            %
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%

clc;clear;
LumleyTri_str='../../StatOut/LumleyTri.txt';
LumleyPicName='../../LumleyTriangle';

%png jpeg, tiff, bmp,  gif, emf  
FormatType=6;
resolution=120;
yPlusTick=0:30:180;
titleStr='Re_{\tau}\rm=180, Closed-channel';

%% normally nothing need to be changed below
rstr=['-r',num2str(resolution)];
if(FormatType==1) 
  format='png'; 
  subfix='png';
elseif(FormatType==2)
  format='jpeg'; 
  subfix='jpeg'; 
elseif(FormatType==3)
  format='tiff'; 
  subfix='tif'; 
elseif(FormatType==4)
  format='bitmap'; 
  subfix='bmp'; 
elseif(FormatType==5)
  format='gif';
  subfix='gif'; 
elseif(FormatType==6)
  format='meta'; %jpeg, tiff, bitmap, gif, meta
  subfix='emf';  
else
  error('wrong pictuer type');     
end

%% Read Lumley triangle data
fid=fopen(LumleyTri_str,'r');
nyc = 0;
fgets(fid);
while (feof(fid)==0)
  nyc=nyc +1;
  str=strtrim( fgets(fid) );
  dlt=sscanf(str,'%f');
  if(isempty(dlt)==1)
    break;
  end
end
frewind(fid);
dataE=zeros(nyc,4);

fgets(fid);
for k=1:nyc
  str=strtrim( fgets(fid) );
  dlt=sscanf(str,'%f');
  dataE(k,:)=dlt;
end
fclose(fid);

%% Plot figure
figure(1);hold on;
set(gca,'FontSize',14,'FontName','Times New Roman');

Xi=(0:0.001:1)*0.5-1/6;
Eta=sqrt(1/27+2*Xi.^3);
plot(Xi,Eta,'k-','LineWidth',0.8);
Xi=[0,1/3];
Eta=[0,1/3];
plot(Xi,Eta,'k-','LineWidth',0.8);
Xi=[0,-1/6];
Eta=[0,1/6];
plot(Xi,Eta,'k-','LineWidth',0.8);

xlabel('\fontsize{16}\xi','FontWeight','bold')
ylabel('\fontsize{16}\eta','FontWeight','bold')
title(['\fontsize{16}\it',titleStr,]);

shading interp;box on;
set(gca,'DataAspectRatio',[1.2 1 1],'linewidth',0.7)
set(gca,'TickDir','out','TickLength',[0.015,0.015]);
set(gca,'XLim',[-1/6,1/3],'XTick',(-1:2)/6);
set(gca,'XTickLabel',{'-1/6','0','1/6','1/3'});
set(gca,'YLim',[0,1/3],'YTick',(0:2)/6);
set(gca,'YTickLabel',{'0','1/6','1/3'});
set(gca,'xMinorTick','on','YMinorTick','on');
s=scatter(dataE(:,3),dataE(:,4),25,dataE(:,2));
s.LineWidth=1.0;

colormap('jet');ch=colorbar;
set(ch,'FontSize',14,'FontName','Times New Roman');
set(ch,'linewidth',0.5)
set(ch,'YLim',[0,max(dataE(:,2))]);
set(ch,'YTick',yPlusTick);
title(ch,'\ity^+');
set(ch,'TickDir','out','TickLength',0.01);
axpos=get(gca,'Position');
ch.Position(3)=0.6*ch.Position(3);
set(gca,'Position',axpos);
print(gcf, ['-d',format], rstr,[LumleyPicName,'.',subfix])
hold off;