% Zheng Gong
% 2024-06-10
% Purpose: compute the 3d energy spectrum
function test_energy3
  clear; clc;
  n=16;
  arr_3d_1=randn(n,n,n);
  arr_3d_2=randn(n,n,n);
  Energy1=sum(sum(sum(arr_3d_1.*arr_3d_2)))/n^3;

  %==============
  arr_3d_fft1=fftn(arr_3d_1);
  arr_3d_fft2=fftn(arr_3d_2);
  Energy2=sum(sum(sum(arr_3d_fft1.*conj(arr_3d_fft2))));
  Energy2=real(Energy2)/n^6;
  clear arr_3d_fft1 arr_3d_fft1;

  %==============
  matout1=myfft3d(arr_3d_1);
  matout2=myfft3d(arr_3d_2);
  
  RatioXYZ=2*ones(n,1);
  RatioXYZ(1:2)=1.0;
  RatioXYZ=RatioXYZ/n^2;

  IndXYZ=zeros(n,1);
  for ijk=1:n/2
    it=2*ijk-1;
    IndXYZ(it)=ijk-1;
    IndXYZ(it+1)=ijk-1;
  end
  IndXYZ(2)=n/2;
  
  nnt=ceil(n/2*sqrt(3))+1;
  Energy3=zeros(nnt);
  for iz=1:n
    Ratio3=RatioXYZ(iz);
    iiz=IndXYZ(iz)*IndXYZ(iz);
    for iy=1:n
      Ratio2=RatioXYZ(iy)*Ratio3;
      iiy=IndXYZ(iy)*IndXYZ(iy)+iiz;
      for ix=1:n
        iix=floor(abs(sqrt(IndXYZ(ix)*IndXYZ(ix)+iiy)))+1;
        if(iix<1 || iix>nnt)
          iix
          error('iix Error')
        end
        Ratio1=RatioXYZ(ix)*Ratio2;
        Energy3(iix)=Energy3(iix) +Ratio1*matout1(ix,iy,iz)*matout2(ix,iy,iz);
      end
    end
  end

  fprintf('Err1=%25.16f\n',Energy2/Energy1)
  fprintf('Err2=%25.16f\n',sum(sum(Energy3))/Energy1)
end

function MatOut=myfft3d(MatIn)
  nx=size(MatIn, 1);
  ny=size(MatIn, 2);
  nz=size(MatIn, 3);
  
  MatOut=MatIn;
  VecIn=zeros(nx,1);
  for ik=1:nz
    for ij=1:ny
      for ii=1:nx
        VecIn(ii)=  MatOut(ii, ij, ik);
      end
      VecOut=myfft1d(VecIn);
      for ii=1:nx
        MatOut(ii,ij,ik) = VecOut(ii);
      end
    end
  end
  VecIn=zeros(ny,1);
  for ik=1:nz
    for ii=1:nx
        for ij=1:ny
          VecIn(ij)=  MatOut(ii, ij, ik);  
        end
        VecOut=myfft1d(VecIn);
        for ij=1:ny
            MatOut(ii,ij,ik) = VecOut(ij);
        end        
    end
  end
  VecIn=zeros(nz,1);
  for ij=1:ny
    for ii=1:nx
      for ik=1:nz
        VecIn(ik)=  MatOut(ii, ij, ik);   
      end
      VecOut=myfft1d(VecIn);
      for ik=1:nz
        MatOut(ii,ij,ik) = VecOut(ik);  
      end      
    end
  end
end

function VecOut=myfft1d(VecIn)
  nLen=length(VecIn);
  VecOut=zeros(nLen,1);
  VecOut_imag=fft(VecIn);
  VecOut(1)=real(VecOut_imag(1));
  VecOut(2)=real(VecOut_imag(nLen/2+1));
  for m=2:nLen/2
    m1=2*m-1;
    m2=m1+1;
    VecOut(m1)=real(VecOut_imag(m));
    VecOut(m2)=imag(VecOut_imag(m));
  end
end
