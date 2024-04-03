function prepare_xdmf_interface()
  clc;clear;
  int_byte=4;
  real_byte=4;
  nWrite=4000;
  numPrtcl=1000000;
  iSerial=20:20:200000;
  OutXmfName='PrtclInfo.xdmf';
  SourcePrefix='./PrtclArrange_';
  
  % Initialize the XDMF/XDF file
  fid=fopen(OutXmfName,'wt');
  fprintf(fid,'<?xml version="1.0" ?>\n');
  fprintf(fid,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
  fprintf(fid,'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">\n');
  fprintf(fid,'<Domain>\n');

  % Time series
  fprintf(fid,'  <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">\n');
  fprintf(fid,'    <Time TimeType="List">\n');
  fprintf(fid,'      <DataItem Format="XML" NumberType="Int" Dimensions="     %d">\n',length(iSerial));
  for itime=iSerial
    fprintf(fid,'%7d',itime);
    if(mod(itime,300)==0); fprintf(fid,'\n'); end;
  end
  fprintf(fid,'\n');
  fprintf(fid,'      </DataItem>\n');
  fprintf(fid,'    </Time>\n');

  % Information
  disp_xmf=0;
  for itime=iSerial
    iWrite=mod(itime,nWrite);
    if(iWrite==0) 
      iWrite=itime;
    else
      iWrite=ceil(itime/nWrite)*nWrite;
    end
    SourceStr=sprintf('%s%10.10d',SourcePrefix,iWrite);
    fprintf(fid,'    <Grid Name="T%10.10d" GridType="Uniform">\n',itime);
    fprintf(fid,'      <Topology TopologyType="Polyvertex" NodesPerElement="        %d"/>\n',numPrtcl);
    fprintf(fid,'      <Geometry GeometryType="XYZ">\n');
    fprintf(fid,'        <DataItem Format="Binary" DataType="Float" Precision="4" Endian="Native" Dimensions=" 3%9d" Seek="%15d">\n', ...
            numPrtcl,disp_xmf);
    fprintf(fid,'          %s\n',SourceStr); 
    fprintf(fid,'        </DataItem>\n');
    fprintf(fid,'      </Geometry>\n');
    disp_xmf=disp_xmf+ numPrtcl*real_byte*3;
    disp_xmf=Write_XDMF_One(fid,3,real_byte,numPrtcl,SourceStr,'linVec',    'Vector','Float', disp_xmf);
    disp_xmf=Write_XDMF_One(fid,3,real_byte,numPrtcl,SourceStr,'RotVel',    'Vector','Float', disp_xmf);    
    disp_xmf=Write_XDMF_One(fid,3,real_byte,numPrtcl,SourceStr,'FpForce',   'Vector','Float', disp_xmf);
    disp_xmf=Write_XDMF_One(fid,3,real_byte,numPrtcl,SourceStr,'FpTorque',  'Vector','Float', disp_xmf);
    disp_xmf=Write_XDMF_One(fid,3,real_byte,numPrtcl,SourceStr,'CntctForce','Vector','Float', disp_xmf);    
    disp_xmf=Write_XDMF_One(fid,3,real_byte,numPrtcl,SourceStr,'Torque',    'Vector','Float', disp_xmf);
    disp_xmf=Write_XDMF_One(fid,1, int_byte,numPrtcl,SourceStr,'iTime',     'Scalar','Int', disp_xmf);
    disp_xmf=Write_XDMF_One(fid,1, int_byte,numPrtcl,SourceStr,'id',        'Scalar','Int', disp_xmf);
    disp_xmf=Write_XDMF_One(fid,1, int_byte,numPrtcl,SourceStr,'iType',     'Scalar','Int', disp_xmf);
    disp_xmf=Write_XDMF_One(fid,1, int_byte,numPrtcl,SourceStr,'isCntct',   'Scalar','Int', disp_xmf);
    fprintf(fid,'    </Grid>\n');
    if(mod(itime,nWrite)==0); disp_xmf=0; end;
  end
  fprintf(fid,'  </Grid>\n');
  fprintf(fid,'</Domain>\n');
  fprintf(fid,'</Xdmf>\n');
  fclose(fid);
end


function dispOut=Write_XDMF_One(fid,ndim,iprec,nPOut,chFile,chName,chAttribute,chDataType,disp)
  fprintf(fid,'      <Attribute Type="%s" Center="Node" Name="%s">\n',chAttribute,chName);
  fprintf(fid,'%s%s%s%1d%s%s%2d%9d%s%15d%s\n','        <DataItem Format="Binary" DataType="',chDataType, ...
              '" Precision="',iprec,'" Endian="Native"',' Dimensions="',ndim,nPOut,'" Seek="',disp,'">')
  fprintf(fid,'          %s\n',chFile);
  fprintf(fid,'        </DataItem>\n');
  fprintf(fid,'      </Attribute>\n');
  dispOut = disp+nPOut*ndim*iprec;
end

