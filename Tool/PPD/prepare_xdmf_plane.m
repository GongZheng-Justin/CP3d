function prepare_xdmf_plane()
  clc;clear;
  xlx=2;
  yly=0.50223*0.02+0.01;
  zlz=0.1;
  nxc=22400;
  nzc=1120;
  real_byte=4;
  nWrite=4000;
  numPrtcl=1000000;
  iSerial=20:20:8000;
  OutXmfName='PlaneInfo.xdmf';
  SourcePrefix='./uxcPlane_y02_';
  
  % Initialize the XDMF/XDF file
  fid=fopen(OutXmfName,'wt');
  fprintf(fid,'<?xml version="1.0" ?>\n');
  fprintf(fid,'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n');
  fprintf(fid,'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">\n');
  fprintf(fid,'<Domain>\n');
  
  % Grid
  fprintf(fid,'  <Topology name="TOPO" TopologyType="3DRectMesh" Dimensions="%6d%6d%6d"/>\n',nzc,1,nxc);
  fprintf(fid,'  <Geometry name="GEO" GeometryType="VXVYVZ">\n');
  fprintf(fid,'    <DataItem Format="XML" DataType="Float" Precision="4" Endian="Native" Dimensions="%d">\n',nxc);
  for ic=1:nxc
    fprintf(fid,'%14.7f',(ic-1)/nxc*xlx);
    if(mod(ic,15)==0); fprintf(fid,'\n'); end;
  end
  fprintf(fid,'    </DataItem>\n');
  fprintf(fid,'    <DataItem Format="XML" DataType="Float" Precision="4" Endian="Native" Dimensions="%d">\n',1);
  fprintf(fid,'%14.7f\n',yly);
  fprintf(fid,'    </DataItem>\n');  
  fprintf(fid,'    <DataItem Format="XML" DataType="Float" Precision="4" Endian="Native" Dimensions="%d">\n',nzc);
  for kc=1:nzc
    fprintf(fid,'%14.7f',(kc-0.5)/nzc*zlz);
    if(mod(kc,15)==0); fprintf(fid,'\n'); end;
  end
  fprintf(fid,'    </DataItem>\n');
  fprintf(fid,'  </Geometry>\n');
  
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
    fprintf(fid,'      <Topology Reference="/Xdmf/Domain/Topology[1]"/>\n');
    fprintf(fid,'      <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>\n');
    disp_xmf=Write_XDMF_One(fid,real_byte,nxc,1,nzc,SourceStr,'ux',disp_xmf);
    fprintf(fid,'    </Grid>\n');
    if(mod(itime,nWrite)==0); disp_xmf=0; end;
  end
  fprintf(fid,'  </Grid>\n');
  fprintf(fid,'</Domain>\n');
  fprintf(fid,'</Xdmf>\n');
  fclose(fid);
end


function dispOut=Write_XDMF_One(fid,iprec,nxc,nyc,nzc,chFile,chName,disp_in)
  fprintf(fid,'      <Attribute Center="Node" Name="%s">\n',chName);
  fprintf(fid,'%s%1d%s%6d%6d%6d%s%15d%s\n','        <DataItem Format="Binary" DataType="Binary" Precision="',iprec, ...
              '" Endian="Native" Dimensions="',nzc,nyc,nxc,'" Seek="',disp_in,'">')
  fprintf(fid,'          %s\n',chFile);
  fprintf(fid,'        </DataItem>\n');
  fprintf(fid,'      </Attribute>\n');
  dispOut = disp_in+nxc*nyc*nzc*iprec;
end

