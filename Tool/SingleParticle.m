clc;clear;
format long

Diam=1/6;
xlx=2;
yly=20;
zlz=2;
CoordFileOutName='SpheresCoord.dat';

%================================================================================
real_prec='real*8';
real_byte=8;

%PosR%x, PosR%y, PosR%z, Diameter, pType
Prtcl_PosD=zeros(4,1);
Prtcl_PosD(1)=xlx*0.5;
Prtcl_PosD(2)=19.5;
Prtcl_PosD(3)=zlz*0.5;
Prtcl_PosD(4)=Diam;

fid=fopen(CoordFileOutName,'w');
fwrite(fid,[Prtcl_PosD;1],real_prec);
fclose(fid);
format short
