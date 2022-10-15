clear;clc;
nxc=1680;                % Grid number in x-dir 
nyc=144;                 % Grid number in y-dir
nzc=1050;                % Grid number in z-dir
yc=load('yMesh.txt');
FieldOut=zeros(nxc,nyc,nzc);
for k=1:nyc
 FieldOut(:,k,:)=yc(k);   
end
fid=fopen('ux_ypos','w');
fwrite(fid,FieldOut,'real*4');
fclose(fid);