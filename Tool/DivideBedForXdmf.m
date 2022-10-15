clc;clear;
format long

%%
ySet=0;
yEnd=30E-3;
CoordFileInName ='./PartVisuForSettling0000100000';
OutPrefix='SphereXdmf';

OutName=[OutPrefix,'.dat'];
CoordFileOutName=['./',OutName];
OutXmfName=['./',OutPrefix,'.xmf']; 
SourceStr=['                    ',OutPrefix,'.dat\n']; 

int_prec='integer*4';
int_byte=4;
real_prec='real*8';
real_byte=8;

%% read binary files
fid1=fopen(CoordFileInName,'r');
fid2=fopen(CoordFileOutName,'w');

fseek(fid1, 0,'eof');
FileSize=ftell(fid1);

%PosR%x, PosR%y, PosR%z, Diameter, pType
SingleSize=real_byte*3+real_byte+int_byte;
nPTotal=FileSize/SingleSize;

% Determine the particle to output
nPOut   =0;
disp_Pos=0;
for k=1:nPTotal
  fseek(fid1, disp_Pos, 'bof');
  Prtcl_Pos=fread(fid1, [1,3], real_prec);
  if(Prtcl_Pos(2)>=ySet && Prtcl_Pos(2)<yEnd)
     nPOut=nPOut+1;
  end
  disp_Pos=disp_Pos+real_byte*3;
end

%======================== 1
disp_Pos  =0;
for k=1:nPTotal
  fseek(fid1, disp_Pos, 'bof');
  Prtcl_Pos=fread(fid1, [1,3], real_prec);
  if(Prtcl_Pos(2)>=ySet && Prtcl_Pos(2)<yEnd)
     fwrite(fid2,Prtcl_Pos,real_prec);
  end
  disp_Pos=disp_Pos+real_byte*3;
end

%======================== 2
disp_Pos  =0;
disp_Diam =nPTotal*real_byte*3;
for k=1:nPTotal
  fseek(fid1, disp_Pos, 'bof');
  Prtcl_Pos=fread(fid1, [1,3], real_prec);
  if(Prtcl_Pos(2)>=ySet && Prtcl_Pos(2)<yEnd)
     fseek(fid1, disp_Diam, 'bof');
     Prtcl_Diam=fread(fid1, [1,1], real_prec);
     fwrite(fid2,Prtcl_Diam,real_prec);
  end
  disp_Pos=disp_Pos+real_byte*3;
  disp_Diam=disp_Diam+real_byte;
end

%======================== 3
disp_Pos  =0;
disp_pType=nPTotal*real_byte*4;
for k=1:nPTotal
  fseek(fid1, disp_Pos, 'bof');
  Prtcl_Pos=fread(fid1, [1,3], real_prec);  
  if(Prtcl_Pos(2)>=ySet && Prtcl_Pos(2)<yEnd)
     fseek(fid1, disp_pType, 'bof');
     Prtcl_pType=fread(fid1, [1,1], int_prec);
     fwrite(fid2,Prtcl_pType,int_prec);     
  end  
  disp_Pos=disp_Pos+real_byte*3;
  disp_pType=disp_pType+int_byte;  
end
fclose(fid1);
fclose(fid2);

% write the xmf file
fid3=fopen(OutXmfName,'wt');
disp_xmf=0;
fprintf(fid3,'<?xml version="1.0" ?>\n');
fprintf(fid3,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
fprintf(fid3,'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">\n');
fprintf(fid3,'<Domain>\n');
fprintf(fid3,'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n');
fprintf(fid3,'        <Time TimeType="List">\n');
fprintf(fid3,'            <DataItem Format="XML" NumberType="Int" Dimensions="     1">\n');
fprintf(fid3,'                    0            </DataItem>\n');
fprintf(fid3,'        </Time>\n');
fprintf(fid3,'        <Grid Name="T0000000000" GridType="Uniform">\n');
fprintf(fid3,'            <Topology TopologyType="Polyvertex" NodesPerElement="        %d"/>\n',nPOut);
fprintf(fid3,'            <Geometry GeometryType="XYZ">\n');
fprintf(fid3,'                <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Native" Dimensions=" 3      %d" Seek="   %d">\n',nPOut,disp_xmf);
fprintf(fid3,SourceStr);
fprintf(fid3,'                </DataItem>\n');
fprintf(fid3,'            </Geometry>\n');
disp_xmf=disp_xmf+ nPOut*real_byte*3;
fprintf(fid3,'            <Attribute Type="Scalar" Center="Node" Name="Diameter">\n');
fprintf(fid3,'                <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Native" Dimensions=" 1      %d" Seek="   %d">\n',nPOut,disp_xmf);
fprintf(fid3,SourceStr);
fprintf(fid3,'                </DataItem>\n');
fprintf(fid3,'            </Attribute>\n');
disp_xmf=disp_xmf+ nPOut*real_byte;
fprintf(fid3,'            <Attribute Type="Scalar" Center="Node" Name="Type">');
fprintf(fid3,'                <DataItem Format="Binary" DataType="int" Precision="4" Endian="Native" Dimensions=" 1      %d" Seek="   %d">\n',nPOut,disp_xmf);
fprintf(fid3,SourceStr);
fprintf(fid3,'                </DataItem>\n');
fprintf(fid3,'            </Attribute>\n');
fprintf(fid3,'        </Grid>\n');
fprintf(fid3,'    </Grid>\n');
fprintf(fid3,'</Domain>\n');
fprintf(fid3,'</Xdmf>\n');
fclose(fid3);

disp(['nPOut= ',num2str(nPOut)]);
format short
