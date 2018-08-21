function [C,phi,S12,S1,S2,f]=coherencyc_gpu(gamydata,gV1data,params)
[tapers,pad,Fs,fpass,err,trialave]=getparams(params);
[N1,C1]=size(gamydata(:,:,1));
[N2,C2]=size(gV1data(:,:,1));
if N1~=N2; error('inconsistent dimensions'); end;
if C1~=C2; error('inconsistent dimensions'); end;
N=N1;
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass);
tapers=dpsschk(tapers,N,Fs);
gtapers=gpuArray(tapers);
data1=gamydata(:,:,:,ones(size(tapers,2),1));
data1=permute(data1,[1,4,2,3]);
data2=gV1data(:,:,:,ones(size(tapers,2),1));
data2=permute(data2,[1,4,2,3]);
dp1=bsxfun(@times,data1,gtapers);
dp2=bsxfun(@times,data2,gtapers);
J1=fft(dp1,nfft)/Fs;
J2=fft(dp2,nfft)/Fs;
J1=J1(findx,:,:,:);
J2=J2(findx,:,:,:);
S12=squeeze(mean(conj(J1).*J2,2));
S1=squeeze(mean(conj(J1).*J1,2));
S2=squeeze(mean(conj(J2).*J2,2));
if trialave; S12=squeeze(mean(S12,2)); S1=squeeze(mean(S1,2)); S2=squeeze(mean(S2,2)); end;
% C12=S12./sqrt(S1.*S2);
C12=arrayfun(@(x,y,z) x./sqrt(y.*z),S12,S1,S2);
C=abs(C12);
phi=angle(C12);