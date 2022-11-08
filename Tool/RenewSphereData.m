clc;clear;

InputName='SpheresCoord.Case2_';
OutputName='SpheresCoord.Case2';

int_prec='integer*4';
int_byte=4;
real_prec='real*8';
real_byte=8;

fid1=fopen(InputName, 'r');
fid2=fopen(OutputName,'w');

fseek(fid1, 0,'eof');
FileSize=ftell(fid1);
SingleSize=real_byte*3+real_byte+int_byte;
nPTotal=FileSize/SingleSize;

dispS=0;
for k=1:nPTotal
  fseek(fid1, dispS, 'bof');
  PosD=fread(fid1, [1,4], real_prec);
  
  disp_pType=dispS+real_byte*4;
  fseek(fid1, disp_pType, 'bof');
  Prtcl_pType=fread(fid1, [1,1], int_prec);
  
  fwrite(fid2,[PosD,Prtcl_pType],real_prec);
  
  dispS=dispS+SingleSize;
end

fclose(fid1);
fclose(fid2);
