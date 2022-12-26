function FilterData=filter1d(VecIn)
  iCoe=10;
  FilterType=1; % 1:BMF; 2:AMF
  FilterData=zeros(size(VecIn));
  
  nLen=length(VecIn);
  [nfd,nfu]=clc_Filter_nd_nu(nLen,FilterType,iCoe);
  for k=1:nLen
    FilterData(k)=mean(VecIn(nfd(k):nfu(k)));   
  end
end