clc;clear;
format long

%%
ySet=0;
yEnd=0.036;
CoordFileInName ='./DEM/Results/PartVisuForLaminarSettling0000200000';
OutPrefix='SpheresCoord';

OutName=[OutPrefix,'.dat'];
CoordFileOutName=['./',OutName];

int_prec='integer*4';
int_byte=4;
real_prec='real*8';
real_byte=8;

%%read binary files
fid1=fopen(CoordFileInName,'r');
fid2=fopen(CoordFileOutName,'w');

fseek(fid1, 0,'eof');
FileSize=ftell(fid1);

%PosR%x, PosR%y, PosR%z, Diameter, pType
SingleSize=real_byte*3+real_byte+int_byte;
nPTotal=FileSize/SingleSize;

%Begin to write
nPOut   =0;
disp_Pos=0;
disp_Diam =nPTotal*real_byte*3;
disp_pType=nPTotal*real_byte*4;
for k=1:nPTotal
  fseek(fid1, disp_Pos, 'bof');
  Prtcl_Pos=fread(fid1, [1,3], real_prec);
  if(Prtcl_Pos(2)>=ySet && Prtcl_Pos(2)<yEnd)
     nPOut=nPOut+1;

     fseek(fid1, disp_Diam, 'bof');
     Prtcl_Diam =fread(fid1, [1,1], real_prec);
     fseek(fid1, disp_pType, 'bof');
     Prtcl_pType=fread(fid1, [1,1], int_prec);
     fwrite(fid2,[Prtcl_Pos,Prtcl_Diam,Prtcl_pType],real_prec);
  end
  disp_Pos=disp_Pos+real_byte*3;
  disp_Diam=disp_Diam+real_byte;
  disp_pType=disp_pType+int_byte; 
end
fclose(fid1);
fclose(fid2);

disp(['nPOut= ',num2str(nPOut)]);
format short
