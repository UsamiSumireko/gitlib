clear
mergeblocklabel=1;
str='D:\data extracted\170505';
% prestr=[str '\_prerev'];
% load([prestr '\V1conblockunitpsth.mat']);
% prerevpsth=V1conblockunitpsth;
% clearvars V1conblockunitpsth
% lr='r';
load([str '\amyconblockunitlfp.mat']);
load([str '\V1conblockunitlfp.mat']);
% load([pre '\V1conblockunitpsth.mat']);
lr=amyconblockunitlfp.LR;
allpuff=[];
allneu=[];
allphip=[];
allphin=[];
multineucnd=[4, 5, 6;
             7, 8, 9;
             10, 11, 12];
if strcmpi(lr,'r')
    cndlist=[1 2 3 4 5 6 4 5 6 4 5 6];
elseif strcmpi(lr,'l')
    cndlist=[4 5 6 1 2 3 1 2 3 1 2 3];
end
params.Fs = 1000; % sampling frequency
params.tapers = [3 5]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 2;
params.err =[2 0.05];
bwin=201:600;
npairsession=32;
nrml=2000;
[~,s2,s3,s4,s5]=size(amyconblockunitlfp.lfp);
[~,~,~,ss4,~]=size(V1conblockunitlfp.lfp);
s1=numel(bwin);
ss1=s1;
if mergeblocklabel
    nmerge=2;
    mb=s5:-nmerge:nmerge;
    mb=fliplr(mb);
    s5origin=s5;
    s5=floor(s5/nmerge);
    temp=reshape(amyconblockunitlfp.spikevalid(:,mb(1)-nmerge+1:s5origin),s4,[],s5);
    temp=permute(temp,[1 3 2]);
    amyvalidcell=prod(temp,3);
    temp=permute(amyconblockunitlfp.lfp(bwin,:,:,:,mb(1)-nmerge+1:s5origin),[1 2 5 3 4]);
    amydata=reshape(temp,s1,[],s5,s3,s4);
    amydata=permute(amydata,[1 2 4 5 3]);
    temp=permute(V1conblockunitlfp.lfp(bwin,:,:,:,mb(1)-nmerge+1:s5origin),[1 2 5 3 4]);
    V1data=reshape(temp,s1,[],s5,s3,ss4);        
    V1data=permute(V1data,[1 2 4 5 3]);
else
    amycellidx=amyconblockunitlfp.spikevalid;
    amydata=amyconblockunitlfp.lfp(bwin,:,:,:,:);
    V1data=V1conblockunitlfp.lfp(bwin,:,:,:,:);
end

for iblock=1:s5
    amyvalidcell=amyconblockunitlfp.spikevalid(:,iblock);
    s4=sum(amyvalidcell);
    amypuff=amydata(:,:,1:3,amyvalidcell,iblock);
    amypuff=reshape(amypuff,s1,[]);
    amyneu=amydata(:,:,4:6,amyvalidcell,iblock);
    amyneu=reshape(amyneu,s1,[]);
    V1puff=V1data(:,:,1:3,:,iblock);
    V1puff=reshape(V1puff,s1,[]);
    V1neu=V1data(:,:,4:6,:,iblock);
    V1neu=reshape(V1neu,s1,[]);
%     amynanidx=isnan(amypuff);
%     V1nanidx=isnan(V1puff);
%     amypuff=amypuff(~amynanidx & ~V1nanidx);
%     V1puff=V1puff(~amynanidx & ~V1nanidx);
    amypuff=rmlinesc(amypuff,params,[],'n',[50 100]);
    V1puff=rmlinesc(V1puff,params,[],'n',[50 100]);
%     amypuff=filt50hz(amypuff,1000);
%     V1puff=filt50hz(V1puff,1000);
    amypuff=reshape(amypuff,s1,[],s4);
    V1puff=reshape(V1puff,s1,[],ss4);
%     amynanidx=isnan(amyneu);
%     V1nanidx=isnan(V1neu);
%     amyneu=amyneu(~amynanidx & ~V1nanidx);
%     V1neu=V1neu(~amynanidx & ~V1nanidx);    
    amyneu=rmlinesc(amyneu,params,[],'n',[50 100]);
    V1neu=rmlinesc(V1neu,params,[],'n',[50 100]);
%     amyneu=filt50hz(amyneu,1000);
%     V1neu=filt50hz(V1neu,1000);
    amyneu=reshape(amyneu,s1,[],s4);
    V1neu=reshape(V1neu,s1,[],ss4);
    
    camypuff=repelem(amypuff,1,1,ss4);
    cV1puff=repmat(V1puff,1,1,s4);
    camyneu=repelem(amyneu,1,1,ss4);
    cV1neu=repmat(V1neu,1,1,s4);
    
    gamypuff=gpuArray(camypuff);
    gamyneu=gpuArray(camyneu);
    gV1puff=gpuArray(cV1puff);
    gV1neu=gpuArray(cV1neu);
    
    for n=1:ceil(ss4*s4/npairsession)
        if n*npairsession>ss4*s4
            pairidx=1+npairsession*(n-1):ss4*s4;
        else
            pairidx=1+npairsession*(n-1):npairsession*n;
        end
%         t_puff=[];
%         t_neu=[];
%         t_phip=[];
%         t_phin=[];
        blockcellidx.acell=repelem(find(amyvalidcell),ss4,1);
        blockcellidx.vcell=repmat(1:ss4,1,s4)';
        [Cpuff,phip,~,~,~,fp]=coherencyc_gpu(gamypuff(:,:,pairidx),gV1puff(:,:,pairidx),params);
        [Cneu,phin,~,~,~,fn]=coherencyc_gpu(gamyneu(:,:,pairidx),gV1neu(:,:,pairidx),params);
        t_puff(:,pairidx)=gather(Cpuff);
        t_neu(:,pairidx)=gather(Cneu);
        t_phip(:,pairidx)=gather(phip);
        t_phin(:,pairidx)=gather(phin);
        t_f=gather(fp);
        coherencepool(iblock).puff=t_puff;
        coherencepool(iblock).neu=t_neu;
        coherencepool(iblock).phip=t_phip;
        coherencepool(iblock).phin=t_phin;
        coherencepool(iblock).f=t_f;
        coherencepool(iblock).blockcellidx=blockcellidx;
    end        
end
for iblock=1:numel(coherencepool)
    figure
    plot(coherencepool(iblock).f,mean(coherencepool(iblock).puff,2),'r');
    hold on
    plot(coherencepool(iblock).f,mean(coherencepool(iblock).neu,2),'g');
end

for iblock=1:numel(coherencepool)
    figure
    celllist=unique(coherencepool(iblock).blockcellidx.acell);
    for jcell=1:numel(celllist)
        subplot(3,8,jcell)
        jcellidx=coherencepool(iblock).blockcellidx.acell==celllist(jcell);
        plot(coherencepool(iblock).f,mean(coherencepool(iblock).puff(:,jcellidx),2),'r');
        hold on
        plot(coherencepool(iblock).f,mean(coherencepool(iblock).neu(:,jcellidx),2),'g');
        xlim([4 100])
        title(num2str(celllist(jcell)))
    end
end

