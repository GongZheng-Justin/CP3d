function [nfd,nfu]=clc_Filter_nd_nu(nLen,FilterType,iCoe)
  nfd=zeros(nLen,1); nfu=zeros(nLen,1);
  nhalf=nLen/2;if(mod(nLen,2)~=0);nhalf=(nLen+1)/2;end;
  if(FilterType==1)      % BMF, Bandwidth Moving Filter
    if(iCoe>99 || iCoe<1)
      error('iCoe Wrong -1');
    end
    pband=iCoe/100;
    for kh=1:nhalf
      km=nLen+1-kh;
      nfd(kh)=round(kh*(1-pband));
      nfu(kh)=round(kh*(1+pband));
      if(km==nhalf);break;end;
      nfd(km)=nLen+1-nfu(kh);
      nfu(km)=nLen+1-nfd(kh);
    end
  elseif(FilterType==2) % AMF, Averaged Moving Filter
    if(iCoe>nhalf-2)
      error('iCoe Wrong -2');
    end
    for k=1:nLen
      kw=min(iCoe,min(k-1,nLen-k));
      nfd(k)=k-kw;
      nfu(k)=k+kw;
    end
  end
end