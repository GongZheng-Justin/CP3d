clc;clear;
nyc=144;
IsOpenChannel=1;
real_prec='real*8';

if(IsOpenChannel==1)
  nySpec2D=nyc;  
else
  nySpec2D=nyc/2;
end
if(strcmp(real_prec,'real*4')==1) 
  real_byte=4;
elseif(strcmp(real_prec,'real*8')==1)
  real_byte=8;
else
  error('readSpec2D: real_prec wrong');
end

ExactLen=zeros(nySpec2D,3);
MaxLenSave=zeros(nySpec2D,2);
for kny=1:nySpec2D
  FileName=sprintf('%s%4.4d','../../lcs2d/Mc',kny);
  fid=fopen(FileName,'r');
  fseek(fid, 0, 'eof');
  totalbyte=ftell(fid);
  nMc=totalbyte/real_byte/2;
  
  fseek(fid, 0, 'bof');
  Mc=fread(fid,[2,nMc],real_prec);Mc=Mc';
  nLeft=nMc; nStart=1;
  nMaxLen=0; iMaxLen=0;
  while(nLeft>0)
    nLen=Mc(nStart,2);
    if(nLen>nMaxLen)
      nMaxLen=nLen; iMaxLen=nStart;
    end
    nStart=nStart+nLen+1;
    nLeft=nLeft-(nLen+1);
  end
  fclose(fid);
  fprintf('kny=%d nMaxLenT=%d \n',kny,nMaxLen)
  MaxLenSave(kny,1)=nMaxLen;
  MaxLenSave(kny,2)=iMaxLen;
  if(nMaxLen>0)
    [ExactLen(kny,1),iLx]=min(Mc(iMaxLen+1:iMaxLen+nMaxLen,1));
    ExactLen(kny,2)=min(Mc(iMaxLen+1:iMaxLen+nMaxLen,2));
    ExactLen(kny,3)=Mc(iMaxLen+iLx,2);
  end
end

McSave=-99*ones(max(MaxLenSave(:,1))+1,nySpec2D*2);
for kny=1:nySpec2D
  k1=2*kny-1; k2=k1+1;
  nMaxLen=MaxLenSave(kny,1);
  iMaxLen=MaxLenSave(kny,2);
  McSave(1,k1)=kny; McSave(1,k2)=nMaxLen;
  if(nMaxLen<1);continue;end;
  FileName=sprintf('%s%4.4d','../../lcs2d/Mc',kny);
  fid=fopen(FileName,'r');
  fseek(fid, 0, 'eof');
  totalbyte=ftell(fid);
  nMc=totalbyte/real_byte/2;
  fseek(fid, 0, 'bof');
  Mc=fread(fid,[2,nMc],real_prec);Mc=Mc';
  fclose(fid);
  McSave(2:nMaxLen+1,k1:k2)=Mc(iMaxLen+1:iMaxLen+nMaxLen,1:2);
end