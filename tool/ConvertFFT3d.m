%
% Author: Zheng Gong
% Date: 2024-06-19
%   Test the convertion between real FFT and complex FFT
%
function ConvertFFT3d
clc; clear;
n=128;

IndexFFT=zeros(n,1);
IndexFFT(1)=1;
IndexFFT(2)=n/2+1;
for kk=2:n/2
  k1=kk*2-1;
  k2=k1+1;
  IndexFFT(k1)=kk;
  IndexFFT(k2)=kk;
end
RatioFFT=2*ones(n,1);
RatioFFT(1:2)=1;

ai=rand(n,n,n)+1i*rand(n,n,n);
for kk=[1, n/2+1]
  for jk=[1, n/2+1]
    for ik=[1, n/2+1]
      ai(ik,jk,kk)=real(ai(ik,jk,kk)) +1i*0.0;
    end
  end
end
for kk=1:n/2+1
  kt=n+2-kk;
  if(kk==1); kt=1; end
  for jk=1:n
    jt=n+2-jk;
    if(jk==1); jt=1; end
    for ik=1:n
      it=n+2-ik;
      if(ik==1); it=1; end
      ai(it,jt,kt) = conj(ai(ik,jk,kk));
    end
  end
end
ar=ifftn(ai); clear ai;
fprintf('%s%24.15E\n', 'Real Array, average error: ',sum(sum(sum(abs(imag(ar)))))/n^3);
ar = real(ar);
aFFT_real=myfft3d(ar); % Real FFT
aFFT_Cmplx=fftn(ar); % Complex FFT

%=========== From real FFT to complex FFT ===========%
aFFT_IfR=1i*ones(n/2+1,n,n);
for kk=1:n
  kt1=IndexFFT(kk);
  kt2=n+2-kt1;
  Mod3=mod(kk,2);
  Ratio3=RatioFFT(kk);
  for jk=1:n
    jt1=IndexFFT(jk);
    jt2=n+2-jt1;
    Mod2=mod(jk,2);
    Ratio2=RatioFFT(jk)*Ratio3;
    for ik=1:n
      it1=IndexFFT(ik);
      %it2=n+2-it1;
      Mod1=mod(ik,2);
      Ratio1=RatioFFT(ik)*Ratio2;
      if(Ratio1==1)
        aFFT_IfR(it1, jt1, kt1)=aFFT_real(ik,jk,kk);
      elseif(Ratio1==2)
        if(Ratio3==2)
          if(Mod3==1)
            aFFT_IfR(it1, jt1, kt1)= aFFT_real(ik,jk,kk)+ aFFT_real(ik,jk,kk+1)*1i;
            aFFT_IfR(it1, jt1, kt2)= aFFT_real(ik,jk,kk)- aFFT_real(ik,jk,kk+1)*1i;
          end
        elseif(Ratio2==2)
          if(Mod2==1)
            aFFT_IfR(it1, jt1, kt1)= aFFT_real(ik,jk,kk)+ aFFT_real(ik,jk+1,kk)*1i;
            aFFT_IfR(it1, jt2, kt1)= aFFT_real(ik,jk,kk)- aFFT_real(ik,jk+1,kk)*1i;
          end
        else
          if(Mod1==1)
            aFFT_IfR(it1, jt1, kt1)= aFFT_real(ik,jk,kk)+ aFFT_real(ik+1,jk,kk)*1i;
            %aFFT_IfR(it2, jt1, kt1)= aFFT_real(ik,jk,kk)- aFFT_real(ik+1,jk,kk)*1i;
          end
        end
      elseif(Ratio1==4)
        if(Ratio3==1)
          if(Mod1*Mod2==1)
            x=aFFT_real(ik:ik+1, jk:jk+1, kk);
            aFFT_IfR(it1, jt1, kt1)= (x(1,1)-x(2,2))+(x(1,2)+x(2,1))*1i;
            %aFFT_IfR(it2, jt1, kt1)= (x(1,1)+x(2,2))+(x(1,2)-x(2,1))*1i;
            aFFT_IfR(it1, jt2, kt1)= (x(1,1)+x(2,2))-(x(1,2)-x(2,1))*1i;
            %aFFT_IfR(it2, jt2, kt1)= (x(1,1)-x(2,2))-(x(1,2)+x(2,1))*1i;
          end
        elseif(Ratio2==2)
          if(Mod3*Mod1==1)
            x=aFFT_real(ik:ik+1, jk, kk:kk+1);
            aFFT_IfR(it1, jt1, kt1)= (x(1,1)-x(2,2))+(x(1,2)+x(2,1))*1i;
            %aFFT_IfR(it2, jt1, kt1)= (x(1,1)+x(2,2))+(x(1,2)-x(2,1))*1i;
            aFFT_IfR(it1, jt1, kt2)= (x(1,1)+x(2,2))-(x(1,2)-x(2,1))*1i;
            %aFFT_IfR(it2, jt1, kt2)= (x(1,1)-x(2,2))-(x(1,2)+x(2,1))*1i;
          end
        else
          if(Mod2*Mod3==1)
            x=reshape(aFFT_real(ik, jk:jk+1, kk:kk+1),[2,2]);
            aFFT_IfR(it1, jt1, kt1)= (x(1,1)-x(2,2))+(x(1,2)+x(2,1))*1i;
            aFFT_IfR(it1, jt2, kt1)= (x(1,1)+x(2,2))+(x(1,2)-x(2,1))*1i;
            aFFT_IfR(it1, jt1, kt2)= (x(1,1)+x(2,2))-(x(1,2)-x(2,1))*1i;
            aFFT_IfR(it1, jt2, kt2)= (x(1,1)-x(2,2))-(x(1,2)+x(2,1))*1i;
          end
        end
      else
        if(Mod1*Mod2*Mod3==1)
          x=aFFT_real(ik:ik+1, jk:jk+1, kk:kk+1);
          aFFT_IfR(it1, jt1, kt1)=(x(1,1,1)-x(2,2,1)-x(1,2,2)-x(2,1,2)) +(x(1,1,2)-x(2,2,2)+x(1,2,1)+x(2,1,1))*1i;
          aFFT_IfR(it1, jt2, kt1)=(x(1,1,1)+x(2,2,1)+x(1,2,2)-x(2,1,2)) +(x(1,1,2)+x(2,2,2)-x(1,2,1)+x(2,1,1))*1i;
          %aFFT_IfR(it2, jt1, kt1)=(x(1,1,1)+x(2,2,1)-x(1,2,2)+x(2,1,2)) +(x(1,1,2)+x(2,2,2)+x(1,2,1)-x(2,1,1))*1i;
          %aFFT_IfR(it2, jt2, kt1)=(x(1,1,1)-x(2,2,1)+x(1,2,2)+x(2,1,2)) +(x(1,1,2)-x(2,2,2)-x(1,2,1)-x(2,1,1))*1i;
          aFFT_IfR(it1, jt1, kt2)=(x(1,1,1)-x(2,2,1)+x(1,2,2)+x(2,1,2)) -(x(1,1,2)-x(2,2,2)-x(1,2,1)-x(2,1,1))*1i;
          aFFT_IfR(it1, jt2, kt2)=(x(1,1,1)+x(2,2,1)-x(1,2,2)+x(2,1,2)) -(x(1,1,2)+x(2,2,2)+x(1,2,1)-x(2,1,1))*1i;
          %aFFT_IfR(it2, jt1, kt2)=(x(1,1,1)+x(2,2,1)+x(1,2,2)-x(2,1,2)) -(x(1,1,2)+x(2,2,2)-x(1,2,1)+x(2,1,1))*1i;
          %aFFT_IfR(it2, jt2, kt2)=(x(1,1,1)-x(2,2,1)-x(1,2,2)-x(2,1,2)) -(x(1,1,2)-x(2,2,2)+x(1,2,1)+x(2,1,1))*1i;
        end
      end
    end
  end
end
fprintf('%s%24.15E\n', 'From Real FFT to Complex FFT, average error: ',sum(sum(sum(abs( aFFT_IfR-aFFT_Cmplx(1:n/2+1,:,:) ))))/n^3);

%=========== From Complex FFT to Real FFT ===========%
aFFT_RfI=zeros(n,n,n);
for kk=1:n
  kt1=IndexFFT(kk);
  kt2=n+2-kt1;
  Mod3=mod(kk,2);
  Ratio3=RatioFFT(kk);
  for jk=1:n
    jt1=IndexFFT(jk);
    jt2=n+2-jt1;
    Mod2=mod(jk,2);
    Ratio2=RatioFFT(jk)*Ratio3;
    for ik=1:n
      it1=IndexFFT(ik);
      %it2=n+2-it1;
      Mod1=mod(ik,2);
      Ratio1=RatioFFT(ik)*Ratio2;
      if(Ratio1==1)
        aFFT_RfI(ik,jk,kk)=real(aFFT_Cmplx(it1, jt1, kt1));
      elseif(Ratio1==2)
        if(Ratio3==2)
          if(Mod3==1)
            cmplx1=aFFT_Cmplx(it1, jt1, kt1);
            aFFT_RfI(ik,jk,kk)=real(cmplx1);
            aFFT_RfI(ik,jk,kk+1)=imag(cmplx1);
          end
        elseif(Ratio2==2)
          if(Mod2==1)
            cmplx1=aFFT_Cmplx(it1, jt1, kt1);
            aFFT_RfI(ik,jk,kk)=real(cmplx1);
            aFFT_RfI(ik,jk+1,kk)=imag(cmplx1);
          end
        else
          if(Mod1==1)
            cmplx1=aFFT_Cmplx(it1, jt1, kt1);
            aFFT_RfI(ik,jk,kk)=real(cmplx1);
            aFFT_RfI(ik+1,jk,kk)=imag(cmplx1);
          end
        end
      elseif(Ratio1==4)
        if(Ratio3==1)
          if(Mod1*Mod2==1)
            cmplx1=aFFT_Cmplx(it1, jt1, kt1);
            cmplx2=aFFT_Cmplx(it1, jt2, kt1);
            aFFT_RfI(ik,jk,kk) =0.5*(real(cmplx1)+real(cmplx2)); % x(1,1)
            aFFT_RfI(ik+1,jk,kk)=0.5*(imag(cmplx1)+imag(cmplx2)); % x(2,1)
            aFFT_RfI(ik,jk+1,kk)=0.5*(imag(cmplx1)-imag(cmplx2)); % x(1,2)
            aFFT_RfI(ik+1,jk+1,kk)=0.5*(real(cmplx2)-real(cmplx1)); % x(2,2)
          end
        elseif(Ratio2==2)
          if(Mod3*Mod1==1)
            cmplx1=aFFT_Cmplx(it1, jt1, kt1);
            cmplx2=aFFT_Cmplx(it1, jt1, kt2);
            aFFT_RfI(ik,jk,kk) =0.5*(real(cmplx1)+real(cmplx2)); % x(1,1)
            aFFT_RfI(ik+1,jk,kk)=0.5*(imag(cmplx1)+imag(cmplx2)); % x(2,1)
            aFFT_RfI(ik,jk,kk+1)=0.5*(imag(cmplx1)-imag(cmplx2)); % x(1,2)
            aFFT_RfI(ik+1,jk,kk+1)=0.5*(real(cmplx2)-real(cmplx1)); % x(2,2)
          end
        else
          if(Mod2*Mod3==1)
            cmplx1=aFFT_Cmplx(it1, jt1, kt1);
            cmplx2=aFFT_Cmplx(it1, jt2, kt1);
            aFFT_RfI(ik,jk,kk) =0.5*(real(cmplx1)+real(cmplx2)); % x(1,1)
            aFFT_RfI(ik,jk,kk+1)=0.5*(imag(cmplx1)+imag(cmplx2)); % x(2,1)
            aFFT_RfI(ik,jk+1,kk)=0.5*(imag(cmplx1)-imag(cmplx2)); % x(1,2)
            aFFT_RfI(ik,jk+1,kk+1)=0.5*(real(cmplx2)-real(cmplx1)); % x(2,2)   
          end
        end
      else
        if(Mod1*Mod2*Mod3==1)
          cmplx1=aFFT_Cmplx(it1, jt1, kt1); cr1=real(cmplx1); ci1=imag(cmplx1);
          cmplx2=aFFT_Cmplx(it1, jt2, kt1); cr2=real(cmplx2); ci2=imag(cmplx2);
          cmplx3=aFFT_Cmplx(it1, jt1, kt2); cr3=real(cmplx3); ci3=imag(cmplx3);
          cmplx4=aFFT_Cmplx(it1, jt2, kt2); cr4=real(cmplx4); ci4=imag(cmplx4);
          aFFT_RfI(ik,jk,kk)= 0.25*(cr1+cr3+cr2+cr4); %x(1,1,1)
          aFFT_RfI(ik+1,jk+1,kk)=0.25*(cr2+cr4-cr1-cr3); %x(2,2,1)
          aFFT_RfI(ik+1,jk,kk+1)=0.25*(cr3-cr1-cr2+cr4); %x(2,1,2)
          aFFT_RfI(ik,jk+1,kk+1)=0.25*(cr3-cr1+cr2-cr4); %x(1,2,2)
          aFFT_RfI(ik+1,jk,kk)=0.25*(ci1+ci3+ci2+ci4) ; %x(2,1,1)
          aFFT_RfI(ik,jk+1,kk)=0.25*(ci1+ci3-ci2-ci4) ; %x(1,2,1)
          aFFT_RfI(ik,jk,kk+1)=0.25*(ci2-ci4+ci1-ci3) ; %x(1,1,2)
          aFFT_RfI(ik+1,jk+1,kk+1)=0.25*(ci2-ci4-ci1+ci3) ; %x(2,2,2)
        end
      end
    end
  end
end
fprintf('%s%24.15E\n', 'From Complex FFT to Real FFT, average error: ',sum(sum(sum(abs(aFFT_RfI-aFFT_real))))/n^3);
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
