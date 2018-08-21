% function splitcellMKII(destlist,LRlist,selecidx)
clear
destlist={
    'D:\data extracted\160712';
    'D:\data extracted\160713';
    'D:\data extracted\160714';
    'D:\data extracted\160715';
    'D:\data extracted\160716';
    'D:\data extracted\160719';
    'D:\data extracted\160720';
    'D:\data extracted\160721';
    'D:\data extracted\160722';
    'D:\data extracted\160624';
    'D:\data extracted\160625';
    'D:\data extracted\160627';
    'D:\data extracted\160628';
    'D:\data extracted\170426';
    'D:\data extracted\170427';
    'D:\data extracted\170429';
    'D:\data extracted\170430';
    'D:\data extracted\170502';
    'D:\data extracted\170503';
    'D:\data extracted\170516';
    'D:\data extracted\170518';
    'D:\data extracted\170519';
    'D:\data extracted\170522';
    };
LRlist='rrrrrrrrrllllrrrrrrllll';
minusbaselinelabel=0;
thr=0.6;
dprimepool=[];
%%
normpsthlabel=true;
Splitmarker='_45degwin';
% savedir='I:\forprez_161008\SplitMKIII\';
bin=0.001;
win=[-0.1 0.5];
winidx=201:800;
tuningwin=(0.1-win(1))/bin+(0.030/bin+1:0.070/bin);
smoothlength=0.01/bin;
for nday=1:numel(destlist)
    destination=destlist{nday};
    LR=LRlist(nday);
    lr=lower(LR);
    if strcmp(lr,'l')
        CNDlist=[4 5 6 1 2 3 1 2 3 1 2 3];
    elseif strcmp(lr,'r')
        CNDlist=[1 2 3 4 5 6 4 5 6 4 5 6];
    else
        error('LR input is wrong, L or R');
    end
    % destination='G:\MS\160624';
    % LR='l';
    %%
    load([destination '\V1conblockunitpsth.mat']);
    load([destination '\V1oriblockunitpsth.mat']);
    cmppath = 'D:\SN 4566-001303+Probe MK2.cmp';
    cmpinfo = LoadCmp(cmppath,1,0);
    elec{1,1} =cmpinfo{1,1}.RealElec;
    % elec{1,1} =cmpinfo{1,1}.RealElec(1:12,:);
    lr=lower(LR);
    if strcmp(lr,'l')
        CNDlist=[4 5 6 1 2 3 1 2 3 1 2 3];
    elseif strcmp(lr,'r')
        CNDlist=[1 2 3 4 5 6 4 5 6 4 5 6];
    else
        error('LR input is wrong, L or R');
    end
    %%
    selecidx{nday}=true(size(V1conblockunitpsth.psth,4),1);    
    TrainLenth=(win(2)-win(1))/bin;
    traintemp=[];
    mutemp=[];    
    oripsth=permute(V1oriblockunitpsth.psth(winidx,:,:,selecidx{nday},:),[1 2 3 5 4]);
    for m=1:6
        mut=oripsth(tuningwin,:,m,:,:);
        mut=reshape(mut,size(mut,1),[],size(mut,5));
        mutemp(:,m)=squeeze(mean(mean(mut)));
        traint=oripsth(:,:,m,:,:);
        traint=reshape(traint,size(traint,1),[],size(traint,5));
        traintemp(:,:,m,:)=reshape(traint,size(traint,1),size(traint,2),1,size(traint,3));
    end
    [~,maxresidx]=max(mutemp,[],2);
    [x,y]=meshgrid(maxresidx,1:6);
    maxoriidx=x==y;
    orimaxtrain=traintemp(:,:,maxoriidx);
    %Rpcelltuntrain等 每行是一个细胞，不是一个trial
    traintemp=[];
    mutemp=[];
    conpsth=permute(V1conblockunitpsth.psth(winidx,:,:,selecidx{nday},:),[1 2 3 5 4]);
    for m=1:3
        traint=reshape(conpsth(:,:,m,:,:),size(conpsth,1),[],size(conpsth,5));
        traintemp(:,:,CNDlist(m),:)=reshape(traint,size(traint,1),size(traint,2),1,size(traint,3));
        
        traint=reshape(conpsth(:,:,m+3,:,:),size(conpsth,1),[],size(conpsth,5));
        traintemp(:,:,CNDlist(m+3),:)=reshape(traint,size(traint,1),size(traint,2),1,size(traint,3));
    end
    secdpeakwin=321:360;
    conmaxtrain=traintemp(:,:,maxoriidx);
    consecd=squeeze(mean(conmaxtrain(secdpeakwin,:,:),1));
    orisecd=squeeze(mean(orimaxtrain(secdpeakwin,:,:),1));
    consecdstd=std(consecd);
    orisecdstd=std(orisecd);
    dprime=(mean(consecd)-mean(orisecd))./sqrt((consecdstd.^2+orisecdstd.^2)/2);
    secpeakcellidx=dprime>thr;
    secpeakcellidx=secpeakcellidx';
    save([destlist{nday} '\secpeakcellidx_both.mat'],'secpeakcellidx','-v7.3')
    dprimepool=[dprimepool dprime];
end
figure
histogram(dprimepool)