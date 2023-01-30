function FilterData=filter2d(MatIn,FilterFlag,iCoe)
  if(FilterFlag~=1 && FilterFlag~=2 && FilterFlag~=3)
    error('%s %d\n','filter2d: FilterFlag Wrong',FilterFlag);
  end
  %iCoe=1;%3 For LCS-2d
  FilterType=2; % 1:BMF; 2:AMF
  FilterData=zeros(size(MatIn));
  
  nLx=size(MatIn,1); nLy=size(MatIn,2);
  if(FilterFlag~=2)
    [nfd,nfu]=clc_Filter_nd_nu(nLx,FilterType,iCoe);
    for k2=1:nLy
      for k1=1:nLx
        FilterData(k1,k2)=mean(MatIn(nfd(k1):nfu(k1),k2));   
      end
     end
  end
  if(FilterFlag~=1)
    [nfd,nfu]=clc_Filter_nd_nu(nLy,FilterType,iCoe);
    for k2=1:nLy
      nfdk=nfd(k2);nfuk=nfu(k2);
      for k1=1:nLx
        FilterData(k1,k2)=mean(MatIn(k1,nfdk:nfuk));   
      end
     end
  end
end