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
for iday=1:numel(folderlist)
    netpuff=[];
    netneu=[];
    str=['D:\result\170828 mergeblock lfpspkcoh\' folderlist{iday}(end-5:end)];
    if exist([str '\ori_lfp_lfp_coherencepool_g.mat'],'file')
    load([str '\ori_lfp_lfp_coherencepool_g.mat'])
    orituningpool=coherencepool;
    clearvars coherencepool
    fidx=orituningpool.f>45 & orituningpool.f<55;
    p50hz=mean(orituningpool.puff(fidx,:));
    n50hz=mean(orituningpool.neu(fidx,:));
    valididxp=p50hz<thr;
    valididxn=n50hz<thr;
    allpufftare=cat(2,allpufftare,orituningpool.puff(:,valididxp));
    allneutare=cat(2,allneutare,orituningpool.neu(:,valididxn));
    end
    if exist([str '\con_lfp_lfp_coherencepool_g.mat'],'file')
        load([str '\con_lfp_lfp_coherencepool_g.mat'])
        fidx=coherencepool.f>45 & coherencepool.f<55;
        p50hz=mean(coherencepool.puff(fidx,:));
        n50hz=mean(coherencepool.neu(fidx,:));
        valididxp=p50hz<thr;
        valididxn=n50hz<thr;
        allpuff=[allpuff coherencepool.puff(:,valididxp)];
        allneu=[allneu coherencepool.neu(:,valididxn)];
    end
end
f=coherencepool.f;
pufftare=mean(allpufftare,2);
neutare=mean(allneutare,2);
allnetcoh.puff=allpuff-repmat(pufftare,1,size(allpuff,2));
allnetcoh.neu=allneu-repmat(neutare,1,size(allneu,2));

puffste=std(allpufftare,0,2)/sqrt(size(allpufftare,2));
neuste=std(allneutare,0,2)/sqrt(size(allneutare,2));
netpuffste=std(allnetcoh.puff,0,2)/sqrt(size(allnetcoh.puff,2));
netneuste=std(allnetcoh.puff,0,2)/sqrt(size(allnetcoh.neu,2));

figure
xv=f;
yv=mean(pufftare,2);
patch([xv fliplr(xv) xv(1)],[yv+puffste; flipud(yv-puffste); yv(1)+puffste(1)],[1 0.75 0.75]);
hold on
plot(xv,yv,'r');
yv=mean(neutare,2);
patch([xv fliplr(xv) xv(1)],[yv+neuste; flipud(yv-neuste); yv(1)+neuste(1)],[0.75 1 0.75]);
hold on
plot(xv,yv,'g');

figure
xv=f;
yv=mean(allnetcoh.puff,2);
patch([xv fliplr(xv) xv(1)],[yv+netpuffste; flipud(yv-netpuffste); yv(1)+netpuffste(1)],[1 0.75 0.75]);
hold on
plot(xv,yv,'r');
hold on
yv=mean(allnetcoh.neu,2);
patch([xv fliplr(xv) xv(1)],[yv+netneuste; flipud(yv-netneuste); yv(1)+netneuste(1)],[0.75 1 0.75]);
hold on
plot(xv,yv,'g');

figure
yv=mean(allnetcoh.puff,2)-mean(allnetcoh.neu,2);
yv=conv(yv,[1 1 1 1 1],'same')/5;
plot(xv,yv,'b')
line([0 100],[0 0],'LineStyle','--')
