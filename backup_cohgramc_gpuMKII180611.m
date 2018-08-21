function [C,phi,S12,S1,S2,t,f]=backup_cohgramc_gpuMKII180611(data1,data2,movingwin,params,npairsession)
% Multi-taper time-frequency coherence,cross-spectrum and individual spectra - continuous processes
%
% Usage:
%
% [C,phi,S12,S1,S2,t,f,confC,phistd,Cerr]=cohgramc(data1,data2,movingwin,params)
% Input: 
% Note units have to be consistent. Thus, if movingwin is in seconds, Fs
% has to be in Hz. see chronux.m for more information.
%
%       data1 (in form samples x trials) -- required
%       data2 (in form samples x trials) -- required
%       movingwin (in the form [window winstep] -- required
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms: 
%                    (1) A numeric vector [TW K] where TW is the
%                        time-bandwidth product and K is the number of
%                        tapers to be used (less than or equal to
%                        2TW-1). 
%                    (2) A numeric vector [W T p] where W is the
%                        bandwidth, T is the duration of the data and p 
%                        is an integer such that 2TW-p tapers are used. In
%                        this form there is no default i.e. to specify
%                        the bandwidth, you have to specify T and p as
%                        well. Note that the units of W and T have to be
%                        consistent: if W is in Hz, T must be in seconds
%                        and vice versa. Note that these units must also
%                        be consistent with the units of params.Fs: W can
%                        be in Hz if and only if params.Fs is in Hz.
%                        The default is to use form 1 with TW=3 and K=5
%                     Note that T has to be equal to movingwin(1).
%
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
% Output:
%       C (magnitude of coherency time x frequencies x trials for trialave=0; 
%             time x frequency for trialave=1)
%       phi (phase of coherency time x frequencies x trials for no trial averaging; 
%             time x frequency for trialave=1)
%       S12 (cross spectrum - time x frequencies x trials for no trial averaging; 
%             time x frequency for trialave=1)
%       S1 (spectrum 1 - time x frequencies x trials for no trial averaging; 
%             time x frequency for trialave=1)
%       S2 (spectrum 2 - time x frequencies x trials for no trial averaging; 
%             time x frequency for trialave=1)
%       t (time)
%       f (frequencies)
%       confC (confidence level for C at 1-p %) - only for err(1)>=1
%       phistd - theoretical/jackknife (depending on err(1)=1/err(1)=2) standard deviation for phi
%                Note that phi + 2 phistd and phi - 2 phistd will give 95% confidence
%                bands for phi - only for err(1)>=1 
%       Cerr  (Jackknife error bars for C - use only for Jackknife - err(1)=2)

if nargin < 3; error('Need data1 and data2 and window parameters'); end;
if nargin < 4; params=[];end;

if length(params.tapers)==3 && movingwin(1)~=params.tapers(2);
    error('Duration of data in params.tapers is inconsistent with movingwin(1), modify params.tapers(2) to proceed')
end
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);

[N,Ch]=check_consistency(data1,data2);

Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=max(2^(nextpow2(Nwin)+pad),Nwin);
f=getfgrid(Fs,nfft,fpass); 
Nf=length(f);
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers

winstart=1:Nstep:N-Nwin+1;
nw=length(winstart);

indx=bsxfun(@plus,winstart,transpose(0:Nwin-1));
datawin1=single(data1(indx(:),:,:));
datawin2=single(data2(indx(:),:,:));
datawin1=reshape(datawin1,Nwin,nw,size(data1,2),size(data1,3));
datawin1=permute(datawin1,[1 3 2 4]);
datawin1=reshape(datawin1,Nwin,size(data1,2),[]);
datawin2=reshape(datawin2,Nwin,nw,size(data2,2),size(data2,3));
datawin2=permute(datawin2,[1 3 2 4]);
datawin2=reshape(datawin2,Nwin,size(data2,2),[]);

if nargin<5
npairsession=32;
end
npairmass=size(datawin1,3);
if trialave;
%    C=gpuArray(zeros(nw,Nf,size(data1,3)));
%    S12=gpuArray(zeros(nw,Nf,size(data1,3)));
%    S1=gpuArray(zeros(nw,Nf,size(data1,3)));
%    S2=gpuArray(zeros(nw,Nf,size(data2,3)));
%    phi=gpuArray(zeros(nw,Nf,size(data1,3)));
   C=single(zeros(Nf,npairmass));
   S12=C;
   S1=C;
   S2=C;
   phi=C;
   ngpumcontent=floor(floor(75000000/Nf)/npairsession);
   Ctemp=gpuArray(nan(Nf,ngpumcontent));
   S12temp=gpuArray(nan(Nf,ngpumcontent));
   S1temp=gpuArray(nan(Nf,ngpumcontent));
   S2temp=gpuArray(nan(Nf,ngpumcontent));
   phitemp=gpuArray(nan(Nf,ngpumcontent));
else
   C=single(zeros(Nf,Ch,npairmass));
   S12=C;
   S1=C;
   S2=C;
   phi=C;
   ngpumcontent=floor(floor(75000000/Nf/Ch)/npairsession);
   Ctemp=gpuArray(nan(Nf,Ch,ngpumcontent));
   S12temp=gpuArray(nan(Nf,Ch,ngpumcontent));
   S1temp=gpuArray(nan(Nf,Ch,ngpumcontent));
   S2temp=gpuArray(nan(Nf,Ch,ngpumcontent));
   phitemp=gpuArray(nan(Nf,Ch,ngpumcontent));
end;
startpos=1;
nblockpaircount=ngpumcontent;
for n=1:ceil(npairmass/npairsession)
    if n*npairsession>npairmass
        pairidx=1+npairsession*(n-1):npairmass;        
        tempidx=mod(pairidx-1,ngpumcontent*npairsession)+1;
        nblockpaircount=mod(npairmass-1,ngpumcontent*npairsession)+1;
    else
        pairidx=1+npairsession*(n-1):npairsession*n;
        tempidx=mod(pairidx-1,ngpumcontent*npairsession)+1;                
    end
    gdatawin1=gpuArray(double(datawin1(:,:,pairidx)));
    gdatawin2=gpuArray(double(datawin2(:,:,pairidx)));
    [c,ph,s12,s1,s2,f]=coherencyc_gpu(gdatawin1,gdatawin2,params);
    if trialave
        Ctemp(:,tempidx)=c;
        S12temp(:,tempidx)=s12;
        S1temp(:,tempidx)=s1;
        S2temp(:,tempidx)=s2;
        phitemp(:,tempidx)=ph;
    else
        Ctemp(:,:,tempidx)=c;
        S12temp(:,:,tempidx)=s12;
        S1temp(:,:,tempidx)=s1;
        S2temp(:,:,tempidx)=s2;
        phitemp(:,:,tempidx)=ph;
    end    
    if mod(n,ngpumcontent)==0 && pairidx(end)==npairmass
        if trialave
            C(:,startpos:startpos+nblockpaircount-1)=gather(Ctemp);
            S12(:,startpos:startpos+nblockpaircount-1)=gather(s12temp);
            S1(:,startpos:startpos+nblockpaircount-1)=gather(s1temp);
            S2(:,startpos:startpos+nblockpaircount-1)=gather(s2temp);
            phi(:,startpos:startpos+nblockpaircount-1)=gather(phitemp);
        else
            C(:,:,startpos:startpos+nblockpaircount-1)=gather(Ctemp);
            S12(:,:,startpos:startpos+nblockpaircount-1)=gather(s12temp);
            S1(:,:,startpos:startpos+nblockpaircount-1)=gather(s1temp);
            S2(:,:,startpos:startpos+nblockpaircount-1)=gather(s2temp);
            phi(:,:,startpos:startpos+nblockpaircount-1)=gather(phitemp);
        end
        startpos=startpos+nblockpaircount;
        Ctemp=gpuArray(nan(Nf,Ch,ngpumcontent));
        S12temp=gpuArray(nan(Nf,Ch,ngpumcontent));
        S1temp=gpuArray(nan(Nf,Ch,ngpumcontent));
        S2temp=gpuArray(nan(Nf,Ch,ngpumcontent));
        phitemp=gpuArray(nan(Nf,Ch,ngpumcontent));
    end
end
if trialave
    C=reshape(C,size(C,1),nw,size(data1,3));
    C=permute(C,[2 1 3]);
    S12=reshape(S12,size(S12,1),nw,size(data1,3));
    S12=permute(S12,[2 1 3]);
    S1=reshape(S1,size(S1,1),nw,size(data1,3));
    S1=permute(S1,[2 1 3]);
    S2=reshape(S2,size(S2,1),nw,size(data1,3));
    S2=permute(S2,[2 1 3]);
    phi=reshape(phi,size(phi,1),nw,size(data1,3));
    phi=permute(phi,[2 1 3]);
else
    C=reshape(C,size(C,1),size(C,2),nw,size(data1,3));
    C=permute(C,[3 1 2 4]);
    S12=reshape(S12,size(S12,1),size(S12,2),nw,size(data1,3));
    S12=permute(S12,[3 1 2 4]);
    S1=reshape(S1,size(S1,1),size(S1,2),nw,size(data1,3));
    S1=permute(S1,[3 1 2 4]);
    S2=reshape(S2,size(S2,1),size(S2,2),nw,size(data1,3));
    S2=permute(S2,[3 1 2 4]);
    phi=reshape(phi,size(phi,1),size(phi,2),nw,size(data1,3));
    phi=permute(phi,[3 1 2 4]);
end
winmid=winstart+round(Nwin/2);
t=winmid/Fs;