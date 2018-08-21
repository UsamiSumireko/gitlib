clear
folderlist={
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
allnetcoh=struct('puff',[],'neu',[]);
allpuff=[];
allneu=[];
allpufftare=[];
allneutare=[];
allvalididx=[];
nacell=0;
thr=0.2;
% load(['D:\result\171025 vector sum\validpair.mat']);
load('D:\result\171122 V1cell select\Validpair_V1sigselected.mat');
for iday=1:numel(folderlist)
    netpuff=[];
    netneu=[];
%     str=['D:\result\170828 mergeblock lfpspkcoh\' folderlist{iday}(end-5:end)];
%     str=['D:\result\171122 V1cell select\coherency\' folderlist{iday}(end-5:end)];
    str=['D:\result\171221 posneg amycell coherency\negcells\' folderlist{iday}(end-5:end)];
    if exist([str '\con_lfp_lfp_trialshuffle_coherencepool_g.mat'],'file') && exist([str '\con_lfp_lfp_coherencepool_g.mat'],'file')
    load([str '\con_lfp_lfp_trialshuffle_coherencepool_g.mat'])
    trialshufflepool=coherencepool;
    clearvars coherencepool
    load([str '\con_lfp_lfp_coherencepool_g.mat'])
    taretemp=permute(trialshufflepool.puff,[1 3 2]);
    tare=mean(taretemp,3);
    fidx=coherencepool.f>45 & coherencepool.f<55;
    p50hz=mean(coherencepool.puff(fidx,:));
    valididxp=p50hz<thr;
%     valididxp=validpair{iday};
    allpuff=[allpuff coherencepool.puff(:,valididxp)];
    allpufftare=cat(2,allpufftare,taretemp(:,valididxp,:));
    netpuff=coherencepool.puff(:,valididxp)-tare(:,valididxp);
    allnetcoh.puff=cat(2,allnetcoh.puff,netpuff);
    
    taretemp=permute(trialshufflepool.neu,[1 3 2]);
    
    tare=mean(taretemp,3);
    p50hz=mean(coherencepool.neu(fidx,:));
    valididxn=p50hz<thr;
    allneu=[allneu coherencepool.neu(:,valididxp)];
    allneutare=cat(2,allneutare,taretemp(:,valididxp,:));
    netneu=coherencepool.neu(:,valididxp)-tare(:,valididxp);
    allnetcoh.neu=cat(2,allnetcoh.neu,netneu);
%     aclist=arrayfun(@(x) x.acellidx,coherencepool.blockcellidx);
%     nacell=nacell+numel(unique(aclist));
%     allvalididx=[allvalididx (valididxp & valididxn)];
    end
end
f=coherencepool.f;
puffste=std(allpuff,0,2)/sqrt(size(allpuff,2));
neuste=std(allneu,0,2)/sqrt(size(allneu,2));
netpuffste=std(allnetcoh.puff,0,2)/sqrt(size(allnetcoh.puff,2));
netneuste=std(allnetcoh.puff,0,2)/sqrt(size(allnetcoh.neu,2));
figure
xv=f;
yv=mean(allpuff,2);
patch([xv fliplr(xv) xv(1)],[yv+puffste; flipud(yv-puffste); yv(1)+puffste(1)],[1 0.75 0.75]);
hold on
plot(xv,yv,'r');
hold on
yv=mean(allneu,2);
patch([xv fliplr(xv) xv(1)],[yv+neuste; flipud(yv-neuste); yv(1)+neuste(1)],[0.75 1 0.75]);
hold on
plot(xv,yv,'g');

figure
yv=mean(allnetcoh.puff,2);
patch([xv fliplr(xv) xv(1)],[yv+netpuffste; flipud(yv-netpuffste); yv(1)+netpuffste(1)],[1 0.75 0.75]);
hold on
plot(xv,yv,'r');
hold on
yv=mean(allnetcoh.neu,2);
patch([xv fliplr(xv) xv(1)],[yv+netneuste; flipud(yv-netneuste); yv(1)+netneuste(1)],[0.75 1 0.75]);
hold on
plot(xv,yv,'g');

pufftarestd=std(mean(allpufftare,2),[],3);
neutarestd=std(mean(allneutare,2),[],3);
figure
xv=f;
yv=mean(mean(allpufftare,3),2);
patch([xv fliplr(xv) xv(1)],[yv+pufftarestd; flipud(yv-pufftarestd); yv(1)+pufftarestd(1)],[1 0.75 0.75])
hold on
yv=mean(allpuff,2);
plot(xv,yv,'r')

yv=mean(mean(allneutare,3),2);
patch([xv fliplr(xv) xv(1)],[yv+neutarestd; flipud(yv-neutarestd); yv(1)+neutarestd(1)],[0.75 1 0.75])
hold on
yv=mean(allneu,2);
plot(xv,yv,'g')
