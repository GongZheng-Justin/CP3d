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

xlx=150.796447372310;
zlz=12.5663706143592;
nxc=9216;
nzc=1400;
nyc=144;
IsOpenChannel=1;
real_prec='real*4';
height=1;
FileDir='../../StatOut/OC0550L/';
VarStr='U';
jForLCS=12;
CValue=[0.1 0.1];
CountourText='off';

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
nyIndex=1:nySpec2D;
nxh=nxc/2; nxhp=nxh+1;
nzh=nzc/2; nzhp=nzh+1;
if(strcmp(real_prec,'real*4')==1) 
  real_byte=4;
elseif(strcmp(real_prec,'real*8')==1)
  real_byte=8;
else
  error('readSpec2D: real_prec wrong');
end

for kny=nyIndex
% Denominator
FileStr=[FileDir,'Spec2D_',VarStr,VarStr];
fid=fopen(FileStr,'r');
fseek(fid, 0, 'eof');
totalbyte1=nxhp*nzhp*nySpec2D*real_byte;
totalbyte2=ftell(fid);
if(totalbyte1~=totalbyte2) 
  error('readSpec2D: File Byte Wrong');
end
offset=(kny-1)*nxhp*nzhp*real_byte;
fseek(fid,offset,'bof');
SpecData=fread(fid,nxhp*nzhp,real_prec);
SpecData=reshape(SpecData,[nxhp,nzhp]);
offset=(jForLCS-1)*nxhp*nzhp*real_byte;
fseek(fid,offset,'bof');
SpecDataD=fread(fid,nxhp*nzhp,real_prec);
SpecDataD=reshape(SpecDataD,[nxhp,nzhp]);
SpecDataD=SpecDataD.*SpecData;
SpecDataD=SpecDataD(2:end-1,2:end-1);
fclose(fid);

% Numerator
FileStr=[FileDir,'Spec2D_L',VarStr];
fid=fopen(FileStr,'r');
fseek(fid, 0, 'eof');
totalbyte1=nxhp*nzhp*nySpec2D*real_byte;
totalbyte2=ftell(fid);
if(totalbyte1~=totalbyte2) 
  error('readSpec2D: File Byte Wrong');
end
offset=(kny-1)*nxhp*nzhp*real_byte;
fseek(fid,offset,'bof');
SpecData=fread(fid,nxhp*nzhp,real_prec);
SpecData=reshape(SpecData,[nxhp,nzhp]);
SpecData=SpecData(2:end-1,2:end-1);
SpecData=SpecData.^2;
SpecData=SpecData./SpecDataD;
fclose(fid);

waveNumberX=(1:nxh-1)*2*pi/xlx;
lamdaX=2*pi./waveNumberX';
waveNumberZ=(1:nzh-1)*2*pi/zlz;
lamdaZ=2*pi./waveNumberZ';
[X,Y]=meshgrid(lamdaX,lamdaZ);

hold on
%pcolor(X,Y, SpecData');
%contourf(X,Y, SpecData',10,'k','LineWidth',0.1);
SpecData=filter2d(SpecData,3);
Mc=contour(X,Y,SpecData',CValue,'LineColor','k');
shading interp;box on;
set(gca,'XScale','log');set(gca,'YScale','log');
set(gca,'FontSize',12,'FontName','Times New Roman');
set(gca,'linewidth',0.7)
set(gca,'TickDir','out','TickLength',[0.015,0.015]);
set(gca,'xMinorTick','on');
set(gca,'yMinorTick','on');
xlabel('\fontsize{16}\it\lambda_x/h');
ylabel('\fontsize{16}\it\lambda_z/h');

ch=colorbar;
set(ch,'FontSize',12,'FontName','Times New Roman');
set(ch,'linewidth',0.6);
title(ch,['\fontsize{12}','\it\gamma']);
set(ch,'TickDir','in','TickLength',0.02);
color=[255 255 255;185 255 255;155 255 255;
        45 255 255; 35 255 220;135 255 120;
       236 255  19;255 201   0;255 134   0;
       255  67   0;255  0    0]/255;
colormap(jet); %get current colormap
if(IsSetCaxis==1)
 set(ch,'YLim',[0,MaxCaxis]);caxis([0,MaxCaxis]);
end

axpos=get(gca,'Position');
axpos(4)=PicYscale*axpos(4);
ch.Position(3)=CBarXScale*ch.Position(3);
ch.Position(4)=axpos(4);
set(gca,'Position',axpos);
hold off

    writename=sprintf('%s%4.4d','../../lcs2d/Mc',kny);
    fid=fopen(writename,'w');
    fwrite(fid,Mc,'real*8');
    fclose(fid);
    fprintf('kny=%d\n',kny)
end