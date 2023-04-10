clc;clear;
format long

%% Input file (Need to change if necessary !!!)
xlx=1;
zlz=1;
nxc=128;
nzc=128;
FileName='../CFD/Restart/RestartForLidCavity0000001500';
MeshName='../CFD/Results/yMeshForLidCavity.txt';

%% ==== Nothing needed to change below ====
int_prec='integer*4';
int_byte=4;
real_prec='real*8';
real_byte=8;
fontname='Times New Roman';

%% mesh
yp=load(MeshName);yp=yp(:,2);yly=yp(end);
yc=0.5*(yp(1:end-1)+yp(2:end));
nyc=length(yc);

nxp=nxc+1;  nyp=nyc+1; nzp=nzc+1;
dx=xlx/nxc; dy=yly/nyc; dz=zlz/nzc;
xp=zeros(nxp,1);
for kt=1:nxp
  xp(kt)=dx*(kt-1);
end
xc=0.5*(xp(1:end-1)+xp(2:end));

%% read binary files
varsize=[nxc,nyc]; offset=0;
fid=fopen(FileName,'r');

% ux
ux=zeros(nxp,nyc,nzc);
for kt=1:nzc
  fseek(fid, offset, 'bof');
  ux(1:nxc,1:nyc,kt)=fread(fid, varsize, real_prec);
  offset=offset+nxc*nyc*real_byte;   
end

% uy
uy=zeros(nxc,nyp,nzc);
for kt=1:nzc
  fseek(fid, offset, 'bof');
  uy(1:nxc,1:nyc,kt)=fread(fid, varsize, real_prec);
  offset=offset+nxc*nyc*real_byte;   
end

% uz
uz=zeros(nxc,nyc,nzp);
for kt=1:nzc
  fseek(fid, offset, 'bof');
  uz(1:nxc,1:nyc,kt)=fread(fid, varsize, real_prec);
  offset=offset+nxc*nyc*real_byte;   
end

% pr
pr=zeros(nxc,nyc,nzc);
for kt=1:nzc
  fseek(fid, offset, 'bof');
  pr(1:nxc,1:nyc,kt)=fread(fid, varsize, real_prec);
  offset=offset+nxc*nyc*real_byte;   
end
fclose(fid);

%% ux
[x,y]=meshgrid(xp,yc);x=x';y=y';
figure(1);
h=pcolor(x,y,ux(:,:,nzc/2));
shading interp;box off;
h=xlabel('x'); set(h,'Fontsize',12)
h=ylabel('y'); set(h,'Fontsize',12)
set(gca,'FontSize',12,'FontName',fontname);
set(gca,'DataAspectRatio',[1 1 1],'linewidth',1.2)
set(gca,'TickDir','out','TickLength',[0.005,0.005]);
set(gca,'XLim',[0,xlx],'XTick',0:0.2:1);
set(gca,'xTickLabel',num2str(get(gca,'xTick')','%.1f'));
set(gca,'YLim',[0,yly],'YTick',0:0.2:1);
set(gca,'yTickLabel',num2str(get(gca,'yTick')','%.1f'));

caxis([-0.2,1.0]);colormap(jet);%colormap(parula);
ch=colorbar;title(ch,'u'); set(ch,'Fontsize',12)
set(ch,'YTick',-0.2:0.2:1.0);
set(ch,'TickDir','out');
set(ch,'yTickLabel',num2str(get(ch,'yTick')','%.1f'));
print(gcf, '-dpng','-r360','veloc_x.png');

format short
