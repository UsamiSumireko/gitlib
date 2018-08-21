%   MKII version use high dimension array data structure
clear
folderlist={
%     'D:\data extracted\180111';
%     'D:\data extracted\180112';
%     'D:\data extracted\180117';
%     'D:\data extracted\180118';
%     'D:\data extracted\180119';
%     'D:\data extracted\180120';
%     'D:\data extracted\180124';
%     'D:\data extracted\180125';
%     'D:\data extracted\180126';

%     'D:\data extracted\180207';
%     'D:\data extracted\180208';
%     'D:\data extracted\180209';
%     'D:\data extracted\180210';
%     'D:\data extracted\180212';

%     'D:\data extracted\180227';
%     'D:\data extracted\180228';
%     'D:\data extracted\180301';
%     'D:\data extracted\180302';
    
    'D:\data extracted\180425';
    'D:\data extracted\180426';
    'D:\data extracted\180427';

%     'D:\data extracted\180504';
%     'D:\data extracted\180505';
%     'D:\data extracted\180506';

%     'D:\data extracted\180514';
%     'D:\data extracted\180515';
%     'D:\data extracted\180516';


%     'D:\data extracted\180525';

    };
% LRlist=['ll' 'rrrr' 'lll' 'rrrr' 'rrr' 'lll'];
LRlist=['rrr'];
disklabel=['III'];
RFlabel=1;  %0 contra RF 1 on RF
excude1sttunlabel=0;

bin=0.001;
win=[-0.1 0.5];
winidx=1:600;
blwin=-win(1)/bin+1:(0.1-win(1))/bin;
pkwin=(0.1-win(1))/bin+(0.030/bin+1:0.070/bin);
smoothlength=0.01/bin;
conpsth=cell(4,1);
oripsth=cell(4,1);
allconpsth=cell(4,1);
alloripsth=cell(4,1);
cellconpsth=cell(4,1);
celloripsth=cell(4,1);
mergeconpsth=cell(numel(folderlist),1);
mergeoripsth=cell(numel(folderlist),1);

for iday=1:numel(folderlist)
    iday
    if excude1sttunlabel
        namelist=dir([disklabel(iday) ':\' folderlist{iday}(end-5:end)]);
        numlist=[];
        for i=3:numel(namelist)
            if namelist(i).isdir && strcmp(namelist(i).name(1:9),'DiGratTun')
                numlist=[numlist str2num(namelist(i).name(10:12))];
            end
        end
        blockinterval=numlist(2:end)-numlist(1:end-1);
        poststartnum=find(blockinterval>1)+1;
        assert(numel(poststartnum)==1);        
    end
    
    lr=LRlist(iday);
    if strcmp(lr,'l')
        if RFlabel
            CNDlist=[3 1 1 1 4 2 2 2];
        else
            CNDlist=[2 1 1 1 4 3 3 3];
        end
    elseif strcmp(lr,'r')
        if RFlabel
            CNDlist=[1 3 3 3 2 4 4 4];
        else
            CNDlist=[1 2 2 2 3 4 4 4];            
        end        
    else        
        error('LR input is wrong, L or R');
    end
    load([folderlist{iday} '\V1conblockunitpsth.mat'])
    load([folderlist{iday} '\V1oriblockunitpsth.mat'])
    if RFlabel
        V1oriblockunitpsth.psth=V1oriblockunitpsth.psth(:,:,[1 3 2 4],:,:);
    end
    [s1,s2,s3,s4]=size(V1conblockunitpsth.secstimlabel);
    psthtemp=V1conblockunitpsth.psth(1:600,logical(V1conblockunitpsth.secstimlabel));
%   psthtemp=reshape(psthtemp,800,9,s2,s3,s4);
%   psthtemp=reshape(psthtemp,800,10,s2,s3,s4); % 180227 180228 180301 all 100% secgrating
    psthtemp=reshape(psthtemp,600,sum(V1conblockunitpsth.secstimlabel(:,1,1,1)),s2,s3,s4);
    temp=permute(psthtemp,[1 2 5 3 4]);
    for m=1:4    
    conpsth{m}=reshape(temp(:,:,:,m==CNDlist,:),size(temp,1),[],size(temp,5));
    vidx=conpsth{m}==65535;
    tridx=squeeze(sum(sum(vidx,3),1))==0;
    conpsth{m}=conpsth{m}(:,tridx,:);
    end
    if excude1sttunlabel
        temp=permute(V1oriblockunitpsth.psth(1:600,:,:,:,[2:poststartnum-1 poststartnum+1:end]),[1 2 5 3 4]);
    else
        temp=permute(V1oriblockunitpsth.psth(1:600,:,:,:,:),[1 2 5 3 4]);
    end
    for m=1:4
        oripsth{m}=reshape(temp(:,:,:,m,:),size(temp,1),[],size(temp,5));
        vidx=oripsth{m}==65535;
        tridx=squeeze(sum(sum(vidx,3),1))==0;
        oripsth{m}=oripsth{m}(:,tridx,:);
    end
    for m=1:4
        cellconpsth{m}=cat(2,cellconpsth{m},squeeze(mean(conpsth{m},2)));
        celloripsth{m}=cat(2,celloripsth{m},squeeze(mean(oripsth{m},2)));
    end
    mergeconpsth{iday}=conpsth;
    mergeoripsth{iday}=oripsth;
end
smconpsth=cell(4,1);
smoripsth=cell(4,1);
for m=1:4
    smconpsth{m}=convn(cellconpsth{m},ones(smoothlength,1),'same')/smoothlength;
    smoripsth{m}=convn(celloripsth{m},ones(smoothlength,1),'same')/smoothlength;
    conste{m}=std(smconpsth{m}/bin,0,2)/sqrt(size(smconpsth{m},2));
    conste{m}=conste{m}';
    oriste{m}=std(smoripsth{m}/bin,0,2)/sqrt(size(smoripsth{m},2));
    oriste{m}=oriste{m}';
end

figure
for m=1:4
    subplot(2,2,m)
    vx=(win(1)+bin:bin:win(2))-0.1;
    vy=mean(smconpsth{m}(winidx,:)/bin,2);
    vy=vy';
    patch([vx fliplr(vx) vx(1)],[vy+conste{m}(winidx) fliplr(vy-conste{m}(winidx)) vy(1)+conste{m}(1)],[1 0.8 0.8])
    hold on
    plot(vx,vy,'r')
    vy=mean(smoripsth{m}(winidx,:)/bin,2);
    vy=vy';
    patch([vx fliplr(vx) vx(1)],[vy+oriste{m}(winidx) fliplr(vy-oriste{m}(winidx)) vy(1)+oriste{m}(1)]',[0.8 1 0.8])
    plot(vx,vy,'g')
    title(['RF' num2str(90*(m>2)+45) 'Contra' num2str(90*mod(m+1,2)+45)])
    xlabel('Time (ms)')
    ylabel('FR (Hz)')
end
dp=cell(4,1);
for m=1:4
sd1=var(cellconpsth{m}(1:600,:),0,2);
sd2=var(celloripsth{m}(1:600,:),0,2);
dp{m}=(mean(cellconpsth{m}(1:600,:),2)-mean(celloripsth{m}(1:600,:),2))./sqrt((sd1+sd2)/2);
dp{m}=squeeze(dp{m});
xlabel('Time (ms)')
ylabel('d''')
end

smoothlength=20;
figure
for m=1:4
subplot(2,2,m)
% vx=(win(1)+bin:bin:win(2))-0.1;
vy=convn(dp{m},ones(smoothlength,1),'same')/smoothlength;
% vy=dp45_45;
% yste=std(dp45_45,0,2)/sqrt(size(dp45_45,2));
plot(vx,vy)
% hold on
% patch([vx fliplr(vx) vx(1)],[vy+yste fliplr(vy-yste) vy(1)+yste(1)],[1 0.8 0.8])
end 

%% different cell based dprime
halfmaxidx=cell(4,1);
for m=1:4
    halfmax(m)=max(mean(smoripsth{m},2))/2;
    halfmaxidx{m}=mean(smoripsth{m},2)>halfmax(m);
    peakmeancon=mean(smconpsth{m}(halfmaxidx{m}(1:size(smconpsth{m},1)),:),1);
    peakmeanori=mean(smoripsth{m}(halfmaxidx{m},:),1);
    sd1=var(peakmeancon);
    sd2=var(peakmeanori);
    peakdp(m)=(mean(peakmeancon)-mean(peakmeanori))./sqrt((sd1+sd2)/2);
end
figure
bar([1 2 3 4],peakdp);


%% dprime for every cell
peakdpMKII=cell(4,1);
peakmeanconMKII=cell(4,1);
peakmeanoriMKII=cell(4,1);
for iday=1:numel(folderlist)
    for m=1:4
        peakmeanconMKII{m}=squeeze(mean(mergeconpsth{iday}{m}(halfmaxidx{m}(1:size(smconpsth{m},1)),:,:),1));
        peakmeanoriMKII{m}=squeeze(mean(mergeoripsth{iday}{m}(halfmaxidx{m}(1:size(smoripsth{m},1)),:,:),1));
        sd1=var(peakmeanconMKII{m},0,1);
        sd2=var(peakmeanoriMKII{m},0,1);
        peakdpMKII{m}=[peakdpMKII{m} (mean(peakmeanconMKII{m})-mean(peakmeanoriMKII{m}))./sqrt((sd1+sd2)/2)];
    end
end
for m=1:4
peakdpste(m)=std(peakdpMKII{m})/sqrt(numel(peakdpMKII{m}));
end
figure
yv=cellfun(@mean,peakdpMKII);
errorbar([1 2 3 4],yv,peakdpste);
title(['N= ' num2str(numel(peakmeancon))])
%% dprime time course for every cell
dpMKII=cell(4,1);
for iday=1:numel(folderlist)
    for m=1:4
        conMKII=convn(mergeconpsth{iday}{m}(1:size(smconpsth{m},1),:,:),ones(20,1,1),'same');
        oriMKII=convn(mergeoripsth{iday}{m}(1:size(smconpsth{m},1),:,:),ones(20,1,1),'same');
        sd1=squeeze(var(conMKII,0,2));
        sd2=squeeze(var(oriMKII,0,2));
        dpMKII{m}=cat(2,dpMKII{m},(squeeze(mean(conMKII,2))-squeeze(mean(oriMKII,2)))./sqrt((sd1+sd2)/2));
    end
end
for m=1:4
    dpvalididx(:,m)=sum(isnan(dpMKII{m}),1)==0;
end
dpvalididx=all(dpvalididx,2);
dpste=cell(4,1);
for m=1:4
dpste{m}=std(dpMKII{m}(:,dpvalididx),0,2)/sqrt(size(dpMKII{m}(:,dpvalididx),2));
end

figure
for m=1:4
subplot(2,2,m)
errorbar(-199:380,mean(dpMKII{m}(1:580,dpvalididx),2),dpste{m}(1:580));
title(['RF' num2str(90*(m>2)+45) 'Contra' num2str(90*mod(m+1,2)+45)])
xlabel('Time (ms)')
ylabel('d''')

end
%% Puff vs NEU, merge Ldays and Rdays
conpuffpsth=[];
conneupsth=[];
tunpuffpsth=[];
tunneupsth=[];
for iday=1:numel(folderlist)
    lr=LRlist(iday);
    if strcmp(lr,'l')
        conpuffpsth=cat(2,conpuffpsth,squeeze(mean(cat(2,mergeconpsth{iday}{3},mergeconpsth{iday}{4}),2)));
        conneupsth=cat(2,conneupsth,squeeze(mean(cat(2,mergeconpsth{iday}{1},mergeconpsth{iday}{2}),2)));
        tunpuffpsth=cat(2,tunpuffpsth,squeeze(mean(cat(2,mergeoripsth{iday}{3},mergeoripsth{iday}{4}),2)));
        tunneupsth=cat(2,tunneupsth,squeeze(mean(cat(2,mergeoripsth{iday}{1},mergeoripsth{iday}{2}),2)));
    elseif strcmp(lr,'r')
        conpuffpsth=cat(2,conpuffpsth,squeeze(mean(cat(2,mergeconpsth{iday}{1},mergeconpsth{iday}{2}),2)));
        conneupsth=cat(2,conneupsth,squeeze(mean(cat(2,mergeconpsth{iday}{3},mergeconpsth{iday}{4}),2)));
        tunpuffpsth=cat(2,tunpuffpsth,squeeze(mean(cat(2,mergeoripsth{iday}{1},mergeoripsth{iday}{2}),2)));
        tunneupsth=cat(2,tunneupsth,squeeze(mean(cat(2,mergeoripsth{iday}{3},mergeoripsth{iday}{4}),2)));
    else
        error('LR input is wrong, L or R');
    end
end
figure
subplot(1,2,1)
xv=transpose(-199:400);
yv=mean(conpuffpsth,2);
yste=std(conpuffpsth(1:600,:)/0.001,[],2)/sqrt(size(conpuffpsth,2));
yv=yv(1:600)/0.001;
patch([xv; flipud(xv)],[yv+yste; flipud(yv-yste)],[1 0.7 0.7],'EdgeAlpha',0)
hold on
plot(xv,yv,'r');
yv=mean(tunpuffpsth,2);
yste=std(tunpuffpsth(1:600,:)/0.001,[],2)/sqrt(size(tunpuffpsth,2));
yv=yv(1:600)/0.001;
patch([xv; flipud(xv)],[yv+yste; flipud(yv-yste)],[0.7 1 0.7],'EdgeAlpha',0)
plot(xv,yv,'g');
ylim([0 200])
title('Con Ori PSTH')
xlabel('Time (ms)')
ylabel('FR (Hz)')
subplot(1,2,2)
yv=mean(conneupsth,2);
yste=std(conneupsth(1:600,:)/0.001,[],2)/sqrt(size(conneupsth,2));
yv=yv(1:600)/0.001;
patch([xv; flipud(xv)],[yv+yste; flipud(yv-yste)],[1 0.7 0.7],'EdgeAlpha',0)
hold on
plot(xv,yv,'r');
yv=mean(tunneupsth,2);
yste=std(tunneupsth(1:600,:)/0.001,[],2)/sqrt(size(tunneupsth,2));
yv=yv(1:600)/0.001;
patch([xv; flipud(xv)],[yv+yste; flipud(yv-yste)],[0.7 1 0.7],'EdgeAlpha',0)
plot(xv,yv,'g');
ylim([0 200])
title('Neu Ori PSTH')
xlabel('Time (ms)')
ylabel('FR (Hz)')