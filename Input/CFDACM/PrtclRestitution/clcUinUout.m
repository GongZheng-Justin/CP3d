clc;clear;
tInterval=1/500;
datapath='SphereInfo_Restitution02_01.txt';

%% Get line1 and line2
fid=fopen(datapath,'r');
line1 = 0;
while (feof(fid)==0)
  str=strtrim( fgets(fid) );
  dlt=sscanf(str,'%f');
  line1=line1 +1;
  if(isempty(dlt)==0 && length(dlt)>2)
    break;
  end
end
real_num=length(dlt);
frewind(fid);
line2 = 0;
while (feof(fid)==0)
  fgets(fid);
  line2= line2 +1;
end
fclose(fid);
file_len=line2-line1+1;

%% Read data
dataIn=zeros(file_len,real_num);
fid=fopen(datapath,'r');
for kt=1:line1-1
 str=strtrim( fgets(fid) ); 
end
for kt=line1:line2
    idl=kt-line1+1;
    str=strtrim( fgets(fid) );
    dlt=sscanf(str,'%f');
    dataIn(idl,:)=dlt';    
end
fclose(fid);

%% Calcuate Uin and Uout
dataNeed=zeros(file_len,3);
dataNeed(:,1)=dataIn(:,1);
dataNeed(:,2)=0.5*(dataIn(:,3)-dataIn(:,5))./dataIn(:,5);
dataNeed(:,3)=dataIn(:,7);
for kt=1:file_len-1
    y1=dataNeed(kt,2);
    y2=dataNeed(kt+1,2);
   if(y1>=0 && y2<0) 
       n1=kt; break;
   end
end

for kt=n1:file_len-1
    y1=dataNeed(kt,2);
    y2=dataNeed(kt+1,2);
   if(y1<0 && y2>=0) 
       n2=kt; break;
   end
end

tco1=dataNeed(n1,1)-dataNeed(n1,2)*(dataNeed(n1+1,1)-dataNeed(n1,1))/(dataNeed(n1+1,2)-dataNeed(n1,2));
tco2=dataNeed(n2,1)-dataNeed(n2,2)*(dataNeed(n2+1,1)-dataNeed(n2,1))/(dataNeed(n2+1,2)-dataNeed(n2,2));
% tc1=(tco1+tco2)/2-tInterval; % for Case1 only
% tc2=(tco1+tco2)/2+tInterval; % for Case1 only
tc1=tco1-tInterval;
tc2=tco2+tInterval;

Uin=100;
for kt=1:n1
    t1=dataNeed(kt,1);
    t2=dataNeed(kt+1,1);
    if(t1<=tc1 && t2>tc1)
        Uin=dataNeed(kt,3)+(dataNeed(kt+1,3)-dataNeed(kt,3))/(t2-t1)*(tc1-t1);
        break;
    end
end

Uout=-100;
for kt=n2:file_len-1
    t1=dataNeed(kt,1);
    t2=dataNeed(kt+1,1);
    if(t1<=tc2 && t2>tc2)
        Uout=dataNeed(kt,3)+(dataNeed(kt+1,3)-dataNeed(kt,3))/(t2-t1)*(tc2-t1);
        break;
    end    
end
fprintf('%s %25.15f\n','Uin= ', Uin);
fprintf('%s %25.15f\n','Uout=',Uout);