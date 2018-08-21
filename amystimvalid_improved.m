clear
str='D:\data extracted\180408';
mergeblocknum=2;
load([str '\amyconblockunitpsth.mat']);
LRlist=amyconblockunitpsth.LR;
reverselabel=find(LRlist(1:end-1)~=LRlist(2:end));
reverselabel=reverselabel+1;
startblock=1;
[s1,s2,s3,s4,s5]=size(amyconblockunitpsth.psth);
lfpMKII=[];
for irev=[reverselabel size(amyconblockunitpsth.psth,5)+1];
    for jblock=startblock:irev-mergeblocknum
        temp=permute(amyconblockunitpsth.psth(:,:,:,:,jblock:jblock+mergeblocknum-1),[1 2 5 3 4]);
        temp=reshape(temp,s1,s2*mergeblocknum,s3,s4);
        lfpMKII=cat(5,lfpMKII,temp);
    end
    startblock=irev;
end
reverselabelMKII=reverselabel-mergeblocknum+1;

%% ttest stimvalid
[~,ss2,~,~,ss5]=size(lfpMKII);
psthassemble=reshape(lfpMKII,s1,ss2*s3,s4,ss5);
blwin=1:200;
stimwin=151:350;
bl=squeeze(mean(psthassemble(blwin,:,:,:,:),1));
stimresp=squeeze(mean(psthassemble(stimwin,:,:,:,:),1));
ht=ttest2(bl,stimresp);
ht=squeeze(ht);
%% bootstrap stimvalid
nbootstrap=1000;
realdiff=mean(stimresp-bl,1);
pool=cat(1,bl,stimresp);
nsample=size(pool,1);
tic
for nbs=1:nbootstrap
    randidx=randi(nsample,nsample,1);
    diff(:,:,:,nsample)=mean(pool(randidx(1:floor(nsample/2)),:,:)-pool(randidx(1+floor(nsample/2):end),:,:),1);
end
toc
btstd=std(diff,0,4);
hbt=abs(realdiff)>1.96*btstd;
hbt=reshape(hbt,s4,ss5);
