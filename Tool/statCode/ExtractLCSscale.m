clc;clear;
height=1;
gammaThreshold=[0.005 0.01];
%VarStr={'uu','vv','ww','pp'};
%DirStr={'x','z','x2','z2'};
% VarStr={'uu'};
% DirStr={'x','x2'};
VarStr={'pp'};
DirStr={'x','z'};
SpecDataDir='../StatOut/OC2000_04P_3P/';

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
    SpecData=filter2d(SpecData,1,1);
    ny=length(yplus); nLamda=length(lamda);
    lamdaOut=-ones(ny,nThreshold);

    LCSFileName=['LCS_',VarStr{kv},'_theta_',DirStr{kd}];
    LCSFileDir=[SpecDataDir,LCSFileName,'.txt'];
    ThetaData=load(LCSFileDir);
    ThetaData=ThetaData(3:end-1,2:end);
    %ThetaData=filter2d(ThetaData,1,1);
    thetaOut=-999*ones(ny,nThreshold);
    for kTh=1:nThreshold
      gammaV=gammaThreshold(kTh);
      for ky=1:ny
        nFit=-1;
        for kLam=nLamda:-1:2 % Wavelength: small -> long, for other
        %for kLam=2:1:nLamda % Wavelength: long -> small, for u,x
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

        theta0=ThetaData(kLam,ky);
        theta1=ThetaData(kMinus,ky);
        thetaOut(ky,kTh)=theta0+(theta1-theta0)/(lamda1-lamda0)*(lamdaOut(ky,kTh)-lamda0);
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
    FileOutName=['ExtractTheta_',VarStr{kv},'_',DirStr{kd}];
    FileOutDir=[SpecDataDir,FileOutName,'.txt'];
    myformat=[repmat('%24.15E',1,nThreshold+1),'\n'];
    fid=fopen(FileOutDir,'wt');
    fprintf(fid,myformat,Retau,gammaThreshold);
    for ky=1:ny
       fprintf(fid,myformat,yplus(ky)/Retau,thetaOut(ky,:));   
    end
    fclose(fid);
    disp( ['read:   ',LCSFileDir,'  sucessfully'] );
  end
end