clc;clear;
height=1;
gammaThreshold=[0.001 0.005 0.01 0.02 0.05 0.06];
VarStr={'uu','vv','ww','pp'};
DirStr={'x','z','x2','z2'};
SpecDataDir='../../StatOut/OC0180/';

nVarStr=length(VarStr);
nDirStr=length(DirStr);
nThreshold=length(gammaThreshold);
for kv=1:nVarStr
  for kd=1:nDirStr
    LCSFileName=['LCS_',VarStr{kv},'_',DirStr{kd}];
    LCSFileDir=[SpecDataDir,LCSFileName,'.txt'];
    SpecData=load(LCSFileDir);
    Retau=SpecData(1,1);
    yplus=SpecData(1,2:end); 
    waveNumber=SpecData(3:end-1,1);
    lamda=2*pi./waveNumber/height;
    SpecData=SpecData(3:end-1,2:end);
    ny=length(yplus); nLamda=length(lamda);
    
    lamdaOut=-ones(ny,nThreshold);
    for kTh=1:nThreshold
      gammaV=gammaThreshold(kTh);
      for ky=1:ny
        nFit=-1;
        for kLam=nLamda:-1:2
          kMinus=kLam-1;
          gamma0=SpecData(kLam,ky);
          gamma1=SpecData(kMinus,ky);
          if(gamma0<=gammaV && gamma1>gammaV) 
            nFit=kLam; break;
          end
        end
        if(nFit<0);continue; end;
        lamda0=lamda(nFit); lamda1=lamda(nFit-1);
        lamdaOut(ky,kTh)=lamda0+(gammaV-gamma0)/(gamma1-gamma0)*(lamda1-lamda0);
      end
    end
    FileOutName=['ExtractLCS_',VarStr{kv},'_',DirStr{kd}];
    FileOutDir=[SpecDataDir,FileOutName,'.txt'];
    myformat=[repmat('%24.15E',1,nThreshold+1),'\n'];
    fid=fopen(FileOutDir,'wt');
    fprintf(fid,myformat,Retau,gammaThreshold);
    for ky=1:ny
       fprintf(fid,myformat,yplus(ky)/Retau,lamdaOut(ky,:));   
    end
    fclose(fid);
    disp( ['read:   ',LCSFileDir,'  sucessfully'] );
  end
end