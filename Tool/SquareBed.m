clc;clear;
format long

nLx=10;
nLy=10;
nLz=10;
Diam=1.0E-3;
OutputDir ='./';
OutPutStr1='FixedSpheresCoord.UhlmannLaminar';
OutPutStr2='FixedSpheresCoordXdmf.UhlmannLaminar';
xdmfStr='FixedSpheresCoordUhlmannLaminar.xmf';

CoordFileOutName =[OutputDir,OutPutStr1];
CoordFileOutName2=[OutputDir,OutPutStr2];
FixedXMFName=[OutputDir,xdmfStr];
%================================================================================
int_prec='integer*4';
int_byte=4;
real_prec='real*8';
real_byte=8;

nSumMax=nLx*nLy*nLz;
Prtcl_Pos=zeros(3,nSumMax);

%PosR%x, PosR%y, PosR%z, Diameter, pType
fid=fopen(CoordFileOutName,'w');

nPTotal=0;
for nkx=1:nLx
  xCoordNow=Diam*(nkx-0.5);
  for nky=1:nLy
    yCoordNow=Diam*(nky-0.5);
    for nkz=1:nLz
      nPTotal=nPTotal+1;
      zCoordNow=Diam*(nkz-0.5);
      Prtcl_PosD=[xCoordNow;yCoordNow;zCoordNow;Diam];
      fwrite(fid,[Prtcl_PosD;1],real_prec);
      Prtcl_Pos(:,nPTotal)=[xCoordNow;yCoordNow;zCoordNow];
    end
  end
end
fclose(fid);

fid=fopen(CoordFileOutName2,'w');
fwrite(fid,Prtcl_Pos(:,1:nPTotal),real_prec);
fwrite(fid,Diam*ones(nPTotal,1),real_prec);
fwrite(fid,ones(nPTotal,1),int_prec);
fclose(fid);

% write the xmf file
fidX=fopen(FixedXMFName,'wt');
disp_xmf=0;
fprintf(fidX,'<?xml version="1.0" ?>\n');
fprintf(fidX,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
fprintf(fidX,'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">\n');
fprintf(fidX,'<Domain>\n');
fprintf(fidX,'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n');
fprintf(fidX,'        <Time TimeType="List">\n');
fprintf(fidX,'            <DataItem Format="XML" NumberType="Int" Dimensions="     1">\n');
fprintf(fidX,'                    0            </DataItem>\n');
fprintf(fidX,'        </Time>\n');
fprintf(fidX,'        <Grid Name="T0000000000" GridType="Uniform">\n');
fprintf(fidX,'            <Topology TopologyType="Polyvertex" NodesPerElement="        %d"/>\n',nPTotal);
fprintf(fidX,'            <Geometry GeometryType="XYZ">\n');
fprintf(fidX,'                <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Native" Dimensions=" 3      %d" Seek="   %d">\n',nPTotal,disp_xmf);
fprintf(fidX,['                    ',OutPutStr2,'\n']);
fprintf(fidX,'                </DataItem>\n');
fprintf(fidX,'            </Geometry>\n');
disp_xmf=disp_xmf+ nPTotal*real_byte*3;
fprintf(fidX,'            <Attribute Type="Scalar" Center="Node" Name="Diameter">\n');
fprintf(fidX,'                <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Native" Dimensions=" 1      %d" Seek="   %d">\n',nPTotal,disp_xmf);
fprintf(fidX,['                    ',OutPutStr2,'\n']);
fprintf(fidX,'                </DataItem>\n');
fprintf(fidX,'            </Attribute>\n');
disp_xmf=disp_xmf+ nPTotal*real_byte;
fprintf(fidX,'            <Attribute Type="Scalar" Center="Node" Name="Type">');
fprintf(fidX,'                <DataItem Format="Binary" DataType="int" Precision="4" Endian="Native" Dimensions=" 1      %d" Seek="   %d">\n',nPTotal,disp_xmf);
fprintf(fidX,['                    ',OutPutStr2,'\n']);
fprintf(fidX,'                </DataItem>\n');
fprintf(fidX,'            </Attribute>\n');
fprintf(fidX,'        </Grid>\n');
fprintf(fidX,'    </Grid>\n');
fprintf(fidX,'</Domain>\n');
fprintf(fidX,'</Xdmf>\n');
fclose(fidX);

disp(['nPTotal= ',num2str(nPTotal)]);
format short
