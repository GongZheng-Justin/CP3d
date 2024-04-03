clc;clear;
format long

xlx=0.01;
yly=0.05;
zlz=0.01;
Diam=1.66666666666666666666667E-3;
Radius=0.5*Diam;
CoordFileOutName='SpheresCoord.dat';

%================================================================================
real_prec='real*8';
real_byte=8;

%PosR%x, PosR%y, PosR%z, Diameter, pType
fid=fopen(CoordFileOutName,'w');

% Particle 1
Prtcl_PosD=zeros(4,1);
Prtcl_PosD(1)= 2.97/600;
Prtcl_PosD(2)=18.96/600;
Prtcl_PosD(3)= 2.97/600;
Prtcl_PosD(4)=Diam;
fwrite(fid,[Prtcl_PosD;1],real_prec);

% Particle 2
Prtcl_PosD=zeros(4,1);
Prtcl_PosD(1)= 3.03/600;
Prtcl_PosD(2)= 21.0/600;
Prtcl_PosD(3)= 3.03/600;
Prtcl_PosD(4)=Diam;
fwrite(fid,[Prtcl_PosD;1],real_prec);

fclose(fid);
format short
