clc;clear;

data_byte=4;
data_prec='real*4';
FileName1='./Spec2D0000000100';
FileName2='./Spec2D0000000100_old';

% Compare file byte
fid1=fopen(FileName1,'r');
fseek(fid1,0,'eof');
Filebyte1=ftell(fid1);
fseek(fid1,0,'bof');

fid2=fopen(FileName2,'r');
fseek(fid2,0,'eof');
Filebyte2=ftell(fid2);
fseek(fid2,0,'bof');
if(Filebyte1 ~= Filebyte2)
  error('CompareTwoBinary: File Byte Wrong');
end
if(mod(Filebyte1,data_byte) ~= 0)
  error('CompareTwoBinary: data_byte Wrong');
end
nData=Filebyte1/data_byte;
DataRead1=fread(fid1, nData, data_prec);
DataRead2=fread(fid2, nData, data_prec);

fclose(fid1);
fclose(fid2);
fprintf('%25.14E \n',sum(abs(DataRead2-DataRead1))/nData)
