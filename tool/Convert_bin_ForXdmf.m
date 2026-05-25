clc;clear;
format long

%%
CoordFileInName ='./FixedSpheresCoord.bin';
OutPrefix='FixedSphereXdmf';

OutName=[OutPrefix,'.bin'];
CoordFileOutName=['./',OutName];
OutXmfName=['./',OutPrefix,'.xmf']; 
SourceStr=['                    ',OutPrefix,'.bin\n']; 

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
SingleSize=real_byte*5;
nPTotal=FileSize/SingleSize;

%======================== 1
disp_Pos  =0;
for k=1:nPTotal
  fseek(fid1, disp_Pos, 'bof');
  Prtcl_Pos=fread(fid1, [1,3], real_prec);
  fwrite(fid2, Prtcl_Pos, real_prec);
  disp_Pos=disp_Pos + SingleSize;
end

%======================== 2
disp_Pos  = real_byte*3;
for k=1:nPTotal
  fseek(fid1, disp_Pos, 'bof');
  Prtcl_Diam=fread(fid1, [1,1], real_prec);
  fwrite(fid2,Prtcl_Diam,real_prec);
  disp_Pos=disp_Pos + SingleSize;
end

%======================== 3
disp_Pos  = real_byte*4;
for k=1:nPTotal
  fseek(fid1, disp_Pos, 'bof');
  Prtcl_pType=fread(fid1, [1,1], real_prec);
  fwrite(fid2, floor(Prtcl_pType+0.1), int_prec);
  disp_Pos=disp_Pos + SingleSize;
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
fprintf(fid3,'            <Topology TopologyType="Polyvertex" NodesPerElement="        %d"/>\n',nPTotal);
fprintf(fid3,'            <Geometry GeometryType="XYZ">\n');
fprintf(fid3,'                <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Native" Dimensions=" 3      %d" Seek="   %d">\n',nPTotal,disp_xmf);
fprintf(fid3,SourceStr);
fprintf(fid3,'                </DataItem>\n');
fprintf(fid3,'            </Geometry>\n');
disp_xmf=disp_xmf+ nPTotal*real_byte*3;
fprintf(fid3,'            <Attribute Type="Scalar" Center="Node" Name="Diameter">\n');
fprintf(fid3,'                <DataItem Format="Binary" DataType="Float" Precision="8" Endian="Native" Dimensions=" 1      %d" Seek="   %d">\n',nPTotal,disp_xmf);
fprintf(fid3,SourceStr);
fprintf(fid3,'                </DataItem>\n');
fprintf(fid3,'            </Attribute>\n');
disp_xmf=disp_xmf+ nPTotal*real_byte;
fprintf(fid3,'            <Attribute Type="Scalar" Center="Node" Name="Type">\n');
fprintf(fid3,'                <DataItem Format="Binary" DataType="int" Precision="4" Endian="Native" Dimensions=" 1      %d" Seek="   %d">\n',nPTotal,disp_xmf);
fprintf(fid3,SourceStr);
fprintf(fid3,'                </DataItem>\n');
fprintf(fid3,'            </Attribute>\n');
fprintf(fid3,'        </Grid>\n');
fprintf(fid3,'    </Grid>\n');
fprintf(fid3,'</Domain>\n');
fprintf(fid3,'</Xdmf>\n');
fclose(fid3);

disp(['nPTotal= ',num2str(nPTotal)]);
format short
