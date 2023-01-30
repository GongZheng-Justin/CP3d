if(strcmp(VarStr,'uu')) 
  iLCSR=iLCSR1DUU;iLCSI=iLCSI1DUU;iSpeR=iSpec1DUU;iLCSR2=iLCSR1DUU2;
elseif(strcmp(VarStr,'vv')) 
  iLCSR=iLCSR1DVV;iLCSI=iLCSI1DVV;iSpeR=iSpec1DVV;iLCSR2=iLCSR1DVV2;  
elseif(strcmp(VarStr,'ww')) 
  iLCSR=iLCSR1DWW;iLCSI=iLCSI1DWW;iSpeR=iSpec1DWW;iLCSR2=iLCSR1DWW2;    
elseif(strcmp(VarStr,'pp'))
  iLCSR=iLCSR1DPP;iLCSI=iLCSI1DPP;iSpeR=iSpec1DPP;iLCSR2=iLCSR1DPP2;   
end

fid=fopen([dir_statOut,'LCS_',VarStr,'_',DirStr,'.txt'],'wt');
if(strcmp(VarStr,'vv'))
  fprintf(fid,myformat,[Retau,yplus_spec_node]);  
else
  fprintf(fid,myformat,[Retau,yplus_spec_center]);
end
for k=1:nxhp
  Ratio=4.0;
  if(k==1 || k==nxhp);Ratio=1.0;end;
  jt=nyc+1-jForLCS(1); 
  SpecData=zeros(1,nych);
  for jc=1:nych
    jsc=nyc+1-jc;
    rTemp1=Ratio*SpectraX(k,jc, iLCSR)^2/(SpectraX(k,jForLCS(1),iSpeR)*SpectraX(k,jc,iSpeR));     
    rTemp2=Ratio*SpectraX(k,jsc,iLCSR)^2/(SpectraX(k,jt,iSpeR)*SpectraX(k,jsc,iSpeR));
    SpecData(jc)=0.5*(rTemp1+rTemp2);
  end
  fprintf(fid,myformat,[2*pi/xlx*(k-1),SpecData]);
end
fclose(fid);

fid=fopen([dir_statOut,'LCS_',VarStr,'_',DirStr,'2.txt'],'wt');
if(strcmp(VarStr,'vv'))
  fprintf(fid,myformat,[Retau,yplus_spec_node]);  
else
  fprintf(fid,myformat,[Retau,yplus_spec_center]);
end
for k=1:nxhp
  Ratio=4.0;
  if(k==1 || k==nxhp);Ratio=1.0;end;
  jt=nyc+1-jForLCS(2); 
  SpecData=zeros(1,nych);
  for jc=1:nych
    jsc=nyc+1-jc;
    rTemp1=Ratio*SpectraX(k,jc, iLCSR2)^2/(SpectraX(k,jForLCS(2),iSpeR)*SpectraX(k,jc,iSpeR));     
    rTemp2=Ratio*SpectraX(k,jsc,iLCSR2)^2/(SpectraX(k,jt,iSpeR)*SpectraX(k,jsc,iSpeR));
    SpecData(jc)=0.5*(rTemp1+rTemp2);
  end
  fprintf(fid,myformat,[2*pi/xlx*(k-1),SpecData]);
end
fclose(fid);

fid=fopen([dir_statOut,'LCS_',VarStr,'_theta_',DirStr,'.txt'],'wt');
YrDist=zeros(nych,1);
if(strcmp(VarStr,'vv'))
  fprintf(fid,myformat,[Retau,yplus_spec_node]);  
  for jc=1:nych
    YrDist(jc)=abs(yp(jc+1)-yp(jForLCS(1)+1));
  end
else
  fprintf(fid,myformat,[Retau,yplus_spec_center]);
  for jc=1:nych    
    YrDist(jc)=abs(yc(jc)-yc(jForLCS(1)));
  end
end
for k=1:nxhp
  SpecData=zeros(1,nych);
  waveNumber=2*pi/xlx*(k-1);
  for jc=1:nych
    jsc=nyc+1-jc;
    YMinusYr=YrDist(jc);
    rTemp1=-atan(SpectraX(k,jc,iLCSI)/SpectraX(k,jc,iLCSR));
    rTemp1=waveNumber*YMinusYr/rTemp1;
    rTemp1=atan(rTemp1)*180/pi;    
    rTemp2=-atan(SpectraX(k,jsc,iLCSI)/SpectraX(k,jsc,iLCSR));
    rTemp2=waveNumber*YMinusYr/rTemp2;
    rTemp2=atan(rTemp2)*180/pi;
    SpecData(jc)=0.5*(rTemp1+rTemp2);
  end
  fprintf(fid,myformat,[2*pi/xlx*(k-1),SpecData]);
end
fclose(fid);