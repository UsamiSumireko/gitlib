%   MKII version use high dimension array data structure
clear
folderlist={
    %Monkey S
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
    %Monkey Milk
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

LRlist=['rrrrrrrrr' 'llll' 'rrrrrr' 'llll'];
% LRlist=['rrrrrr' 'llll']; %for Monkey Milk only
bin=0.001;
win=[-0.1 0.5];
winidx=201:800;
blwin=(-0.1-win(1))/bin+1:(0.1-win(1))/bin;
pkwin=(0.1-win(1))/bin+(0.030/bin+1:0.070/bin);
latewin=(0.1-win(1))/bin+(0.12/bin+1:0.16/bin);
smoothlength=0.01/bin;
allconpsth=cell(6,1);
alloripsth=cell(6,1);
cellconpuffpsth=[];
cellconneupsth=[];
celloripuffpsth=[];
cellorineupsth=[];
offpeakcell=[];
mergeh=[];
mergepuffh=[];
mergeneuh=[];
% midconpuffpsth=[];
% midconneupsth=[];
% midoripuffpsth=[];
% midorineupsth=[];
puffdprimepool=[];
neudprimepool=[];
folderlistMKII=cellfun(@(x) ['H:\' x(end-5:end)],folderlist,'UniformOutput',false);
prepostidx=preposttestblockidx(folderlistMKII);
mergehbl=[];
mergehpk=[];
mergeprepsth=[];
mergepostpsth=[];
mergeprepuffpsth=[];
mergepreneupsth=[];
mergepostpuffpsth=[];
mergepostneupsth=[];
for iday=1:numel(folderlist)
%     iday
    if strcmpi(LRlist(iday),'r')
        cndlist=[1 2 3 4 5 6];
    elseif strcmpi(LRlist(iday),'l')
        cndlist=[4 5 6 1 2 3];
    else
        error('wrong LR')
    end
    load([folderlist{iday} '\V1conblockunitpsth.mat'])
    load([folderlist{iday} '\V1oriblockunitpsth.mat'])
    load([folderlist{iday} '\XYvalidid.mat'])
    load([folderlist{iday} '\ORIvalidid.mat'])   
    nanidx=isnan(V1conblockunitpsth.psth(winidx,:,:,:,:));
    connanvalid=squeeze(0==sum(sum(sum(sum(nanidx,1),2),3),4));
    idx65535=65535==V1conblockunitpsth.psth(winidx,:,:,:,:);
    contbinvalid=squeeze(0==sum(sum(sum(sum(idx65535,1),2),3),4));
    nanidx=isnan(V1oriblockunitpsth.psth(winidx,:,:,:,:));
    orinanvalid=squeeze(0==sum(sum(sum(sum(nanidx,1),2),3),4));
    idx65535=65535==V1oriblockunitpsth.psth(winidx,:,:,:,:);
    oritbinvalid=squeeze(0==sum(sum(sum(sum(idx65535,1),2),3),4));
    convalidblockidx=connanvalid & contbinvalid;
    orivalidblockidx=orinanvalid & oritbinvalid;
%     if ~any(convalidblockidx) || ~any(orivalidblockidx)
    if ~(any(convalidblockidx) && any(orivalidblockidx))
        continue % for iday=1:numel(folderlist)
    end
    
    temp=permute(V1conblockunitpsth.psth(winidx,:,:,:,convalidblockidx),[1 2 5 3 4]);    
    conpuffpsth=reshape(temp(:,:,:,1:3,:),size(temp,1),[],size(temp,5));
    midconpuffpsth=reshape(temp(:,:,:,2,:),size(temp,1),[],size(temp,5));  % midpuff/midneu only include 45/135 for d' compute
    conneupsth=reshape(temp(:,:,:,4:12,:),size(temp,1),[],size(temp,5));
    midconneupsth=reshape(temp(:,:,:,5:3:11,:),size(temp,1),[],size(temp,5));  % midpuff/midneu only include 45/135 for d' compute
    temp=permute(V1oriblockunitpsth.psth(winidx,:,:,:,orivalidblockidx),[1 2 5 3 4]);
    oripuffpsth=reshape(temp(:,:,:,cndlist(1:3),:),size(temp,1),[],size(temp,5));
    midoripuffpsth=reshape(temp(:,:,:,cndlist(2),:),size(temp,1),[],size(temp,5)); % midpuff/midneu only include 45/135 for d' compute
    orineupsth=reshape(temp(:,:,:,cndlist(4:6),:),size(temp,1),[],size(temp,5));
    midorineupsth=reshape(temp(:,:,:,cndlist(5),:),size(temp,1),[],size(temp,5)); % midpuff/midneu only include 45/135 for d' compute
%% for compare the response of pretest and posttest
    pretestpsth=V1oriblockunitpsth.psth(winidx,:,:,:,~prepostidx{iday});
    posttestpsth=V1oriblockunitpsth.psth(winidx,:,:,:,prepostidx{iday});
    prepsth=[];
    postpsth=[];
    hbl=[];
    hpk=[];
    for mcnd=1:6
        pretemp=pretestpsth(:,:,mcnd,:,:);
        pretemp=permute(pretemp,[1 2 5 3 4]);
        pretemp=reshape(pretemp,size(pretemp,1),[],size(pretemp,5));
        prebl=squeeze(mean(pretemp(blwin,:,:),1));
        prepk=squeeze(mean(pretemp(pkwin,:,:),1));
        prepsth(:,:,mcnd)=squeeze(mean(pretemp,2));      
        
        posttemp=posttestpsth(:,:,mcnd,:,:);
        posttemp=permute(posttemp,[1 2 5 3 4]);
        posttemp=reshape(posttemp,size(posttemp,1),[],size(posttemp,5));
        postbl=squeeze(mean(posttemp(blwin,:,:),1));
        postpk=squeeze(mean(posttemp(pkwin,:,:),1));
        postpsth(:,:,mcnd)=squeeze(mean(posttemp,2));        
        hbl(mcnd,:)=ttest2(prebl,postbl);
        hpk(mcnd,:)=ttest2(prepk,postpk);
    end
    mergehbl=cat(2,mergehbl,hbl);
    mergehpk=cat(2,mergehpk,hpk);
    mergeprepsth=cat(2,mergeprepsth,prepsth);
    mergepostpsth=cat(2,mergepostpsth,postpsth);
    mergeprepuffpsth=cat(2,mergeprepuffpsth,prepsth(:,:,cndlist(2)));
    mergepreneupsth=cat(2,mergepreneupsth,prepsth(:,:,cndlist(5)));
    mergepostpuffpsth=cat(2,mergepostpuffpsth,postpsth(:,:,cndlist(2)));
    mergepostneupsth=cat(2,mergepostneupsth,postpsth(:,:,cndlist(5)));
%%    
%     offpeak=cat(2,conpuffpsth(631:640,:,:),conneupsth(831:840,:,:));
%     preoff=cat(2,conpuffpsth(501:600,:,:),conneupsth(701:800,:,:));
%     offpeak=squeeze(mean(offpeak));
%     preoff=squeeze(mean(preoff));
%     [th,tp]=ttest(offpeak,preoff);
%     offpeakcell=[offpeakcell th];
%     
%     
    cellconpuffpsth=cat(2,cellconpuffpsth,squeeze(mean(conpuffpsth,2)));
    cellconneupsth=cat(2,cellconneupsth,squeeze(mean(conneupsth,2)));
    celloripuffpsth=cat(2,celloripuffpsth,squeeze(mean(oripuffpsth,2)));
    cellorineupsth=cat(2,cellorineupsth,squeeze(mean(orineupsth,2)));
    %%%%%%%%%%%%%%%%%%%%%%%
    puffdiff=midconpuffpsth(:,:,:)-repmat(mean(midoripuffpsth(:,:,:),2),1,size(midconpuffpsth,2),1);
    neudiff=midconneupsth(:,:,:)-repmat(mean(midorineupsth(:,:,:),2),1,size(midconneupsth,2),1);
    puffpeakdiff=squeeze(mean(puffdiff(pkwin,:,:),1));
    neupeakdiff=squeeze(mean(neudiff(pkwin,:,:),1));
    hpuff=ttest(puffpeakdiff);
    hneu=ttest(neupeakdiff);
    h=ttest2(puffpeakdiff,neupeakdiff);
    V1sigcell_ttest{iday}=h;
    mergepuffh=[mergepuffh hpuff];
    mergeneuh=[mergeneuh hneu];
    mergeh=[mergeh h];
    %%%%%%%%%%%%%%%%%%%%%%%
    %%
    % d-prime for puff ori and neu ori
    % d-prime for every cell, mean plot with error bar.
    conpuffvalididx=squeeze(sum(sum(midconpuffpsth==65535,1),2))==0;
    conneuvalididx=squeeze(sum(sum(midconneupsth==65535,1),2))==0;
    oripuffvalididx=squeeze(sum(sum(midoripuffpsth==65535,1),2))==0;
    orineuvalididx=squeeze(sum(sum(midorineupsth==65535,1),2))==0;
    valididx=conpuffvalididx & conneuvalididx & oripuffvalididx & orineuvalididx;
    sprintf(['Day ' num2str(iday) 'N= ' num2str(sum(valididx))])
    
    smmidconpuffpsth=convn(midconpuffpsth(:,:,valididx),ones(20,1,1),'same')/20;
    conpuffvar=var(smmidconpuffpsth,0,2);
    smmidconneupsth=convn(midconneupsth(:,:,valididx),ones(20,1,1),'same')/20;
    conneuvar=var(smmidconneupsth,0,2);
    smmidoripuffpsth=convn(midoripuffpsth(:,:,valididx),ones(20,1,1),'same')/20;
    oripuffvar=var(smmidoripuffpsth,0,2);
    smmidorineupsth=convn(midorineupsth(:,:,valididx),ones(20,1,1),'smae')/20;
    orineuvar=var(smmidorineupsth,0,2);
    
    puffdprime=(mean(smmidconpuffpsth,2)-mean(smmidoripuffpsth,2))./sqrt((conpuffvar+oripuffvar)/2);
    puffdprime=squeeze(puffdprime);
    notnanidx=sum(isnan(puffdprime),1)==0;
    puffdprime=puffdprime(:,notnanidx);
    neudprime=(mean(smmidconneupsth,2)-mean(smmidorineupsth,2))./sqrt((conneuvar+orineuvar)/2);
    neudprime=squeeze(neudprime);
    notnanidx=sum(isnan(neudprime),1)==0;
    neudprime=neudprime(:,notnanidx);
    puffdprimepool=cat(2,puffdprimepool,puffdprime);
    neudprimepool=cat(2,neudprimepool,neudprime);
end


smconpuffpsth=convn(cellconpuffpsth,ones(smoothlength,1),'same')/smoothlength;
smconneupsth=convn(cellconneupsth,ones(smoothlength,1),'same')/smoothlength;
smoripuffpsth=convn(celloripuffpsth,ones(smoothlength,1),'same')/smoothlength;
smorineupsth=convn(cellorineupsth,ones(smoothlength,1),'same')/smoothlength;

puffste1=std(smconpuffpsth/bin,0,2)/sqrt(size(smconpuffpsth,2));
puffste1=puffste1';
puffste2=std(smoripuffpsth/bin,0,2)/sqrt(size(smoripuffpsth,2));
puffste2=puffste2';
neuste1=std(smconneupsth/bin,0,2)/sqrt(size(smconneupsth,2));
neuste1=neuste1';
neuste2=std(smorineupsth/bin,0,2)/sqrt(size(smorineupsth,2));
neuste2=neuste2';

figure
subplot(2,1,1)
vx=(win(1)+bin:bin:win(2))-0.1;
vx=vx*1000;
vy=mean(smconpuffpsth(:,:)/bin,2);
vy=vy';
patch([vx fliplr(vx) vx(1)],[vy+puffste1 fliplr(vy-puffste1) vy(1)+puffste1(1)],[1 0.8 0.8])
hold on
h1=plot(vx,vy,'r');
vy=mean(smoripuffpsth(:,:)/bin,2);
vy=vy';
patch([vx fliplr(vx) vx(1)],[vy+puffste2 fliplr(vy-puffste2) vy(1)+puffste2(1)],[0.8 1 0.8])
h2=plot(vx,vy,'g');
xlabel('time (ms)');
ylabel('FR (Hz)');
legend([h1,h2],'puff ori','neutral ori');
title(['N= ' num2str(size(cellconpuffpsth,2))])

subplot(2,1,2)
vy=mean(smconneupsth(:,:)/bin,2);
vy=vy';
patch([vx fliplr(vx) vx(1)],[vy+neuste1 fliplr(vy-neuste1) vy(1)+neuste1(1)],[1 0.8 0.8])
hold on
plot(vx,vy,'r')
vy=mean(smorineupsth(:,:)/bin,2);
vy=vy';
patch([vx fliplr(vx) vx(1)],[vy+neuste2 fliplr(vy-neuste2) vy(1)+neuste2(1)],[0.8 1 0.8])
plot(vx,vy,'g')
xlabel('time (ms)');
ylabel('FR (Hz)');

% figure
% vx=-0.1+bin:bin:0.2;
% offsetidx=501:800;
% preoffidx=501:600;
% vy=mean(smconpuffpsth(offsetidx,:)/bin,2);%-mean(mean(smconpuffpsth(preoffidx,:)/bin));
% vy=vy';
% % patch([vx fliplr(vx) vx(1)],[vy+puffste1(offsetidx) fliplr(vy-puffste1(offsetidx)) vy(1)+puffste1(1)],[1 0.8 0.8])
% hold on
% plot(vx,vy,'r')
% vy=mean(smconneupsth(offsetidx,:)/bin,2);%-mean(mean(smconneupsth(preoffidx,:)/bin));
% vy=vy';
% hold on
% % patch([vx fliplr(vx) vx(1)],[vy+neuste1(offsetidx) fliplr(vy-neuste1(offsetidx)) vy(1)+neuste1(1)],[0.8 1 0.8])
% hold on
% plot(vx,vy,'g')
% 
% % csoffpeak=squeeze(mean(smconpuffpsth(631:640,:)/bin));
% % cspreoff=squeeze(mean(smconpuffpsth(501:600,:)/bin));
% % neuoffpeak=squeeze(mean(smconneupsth(631:640,:,:)/bin));
% % a=mean(cellconpuffpsth(231:270,:))-mean(celloripuffpsth(231:270,:));
% % b=mean(cellconneupsth(231:270,:))-mean(cellorineupsth(231:270,:));
% % subplot(1,2,1)
% % histogram(a)
% % subplot(1,2,2)
% % histogram(b)
% % figure
% % histogram(a-b)

%%
puffdprimeste=std(puffdprimepool,0,2)/sqrt(size(puffdprimepool,2));
neudprimeste=std(neudprimepool,0,2)/sqrt(size(neudprimepool,2));

figure
xv=(-199:400)';
yv=mean(puffdprimepool,2);
patch([xv; flipud(xv)],[yv+puffdprimeste; flipud(yv-puffdprimeste)],[1 0.75 0.75]);
hold on
h1=plot(xv,yv,'r');
yv=mean(neudprimepool,2);
patch([xv; flipud(xv)],[yv+neudprimeste; flipud(yv-neudprimeste)], [0.75 1 0.75]);
hold on
h2=plot(xv,yv,'g');
xlabel('time (ms)');
ylabel('d''');
legend([h1,h2],'puff ori','neutral ori');
%%
% statistical test for psth
conpuffpk=mean(smconpuffpsth(pkwin,:));
conneupk=mean(smconneupsth(pkwin,:));
oripuffpk=mean(smoripuffpsth(pkwin,:));
orineupk=mean(smorineupsth(pkwin,:));
pkpuffchangeratio=conpuffpk./oripuffpk;
pkneuchangeratio=conneupk./orineupk;

conpuffbl=mean(smconpuffpsth(blwin,:));
conneubl=mean(smconneupsth(blwin,:));
oripuffbl=mean(smoripuffpsth(blwin,:));
orineubl=mean(smorineupsth(blwin,:));
blpuffchangeratio=conpuffbl./oripuffbl;
blneuchangeratio=conneubl./orineubl;

conpufflate=mean(smconpuffpsth(latewin,:));
conneulate=mean(smconneupsth(latewin,:));
oripufflate=mean(smoripuffpsth(latewin,:));
orineulate=mean(smorineupsth(latewin,:));
latepuffchangeratio=conpufflate./oripufflate;
lateneuchangeratio=conneulate./orineulate;

% ranktest
[ppuffpk_rank,~]=ranksum(conpuffpk',oripuffpk');
[pneupk_rank,~]=ranksum(conneupk',orineupk');
[ppuffbl_rank,~]=ranksum(conpuffbl',oripuffbl');
[pneubl_rank,~]=ranksum(conneubl',orineubl');
[ppufflate_rank,~]=ranksum(conpufflate',oripufflate');
[pneulate_rank,~]=ranksum(conneulate',orineulate');
[penchancepk_rank,~]=ranksum((conpuffpk-oripuffpk)',(conneupk-orineupk)');
[penchancebl_rank,~]=ranksum((conpuffbl-oripuffbl)',(conneubl-orineubl)');
[penchancelate_rank,~]=ranksum((conpufflate-oripufflate)',(conneulate-orineulate)');

%ttest2
[~,ppuffpk_ttest]=ttest2(conpuffpk',oripuffpk');
[~,pneupk_ttest]=ttest2(conneupk',orineupk');
[~,ppuffbl_ttest]=ttest2(conpuffbl',oripuffbl');
[~,pneubl_ttest]=ttest2(conneubl',orineubl');
[~,ppufflate_ttest]=ttest2(conpufflate',oripufflate');
[~,pneulate_ttest]=ttest2(conneulate',orineulate');
[~,penchancepk_ttest]=ttest2((conpuffpk-oripuffpk)',(conneupk-orineupk)');
[~,penchancebl_ttest]=ttest2((conpuffbl-oripuffbl)',(conneubl-orineubl)');
[~,penchancelate_ttest]=ttest2((conpufflate-oripufflate)',(conneulate-orineulate)');
% bootstrap
bootstraptime=1000;
puffpkpool=[conpuffpk oripuffpk];
neupkpool=[conneupk orineupk];
pkratiopool=[pkpuffchangeratio pkneuchangeratio];
pkenhancediffpool=[conpuffpk-oripuffpk conneupk-orineupk];

puffblpool=[conpuffbl oripuffbl];
neublpool=[conneubl orineubl];
blratiopool=[blpuffchangeratio blneuchangeratio];
blenhancediffpool=[conpuffbl-oripuffbl conneubl-orineubl];

pufflatepool=[conpufflate oripufflate];
neulatepool=[conneulate orineulate];
lateratiopool=[latepuffchangeratio lateneuchangeratio];
lateenhancediffpool=[conpufflate-oripufflate conneulate-orineulate];

ncell=size(smconpuffpsth,2);
puffpkdiffsample=nan(bootstraptime,1);
neupkdiffsample=puffpkdiffsample;
puffbldiffsample=puffpkdiffsample;
neubldiffsample=puffpkdiffsample;
pufflatediffsample=puffpkdiffsample;
neulatediffsample=puffpkdiffsample;
pkenhancediffsample=puffpkdiffsample;
blenhancediffsample=puffpkdiffsample;
lateenhancediffsample=puffpkdiffsample;
pkratiodiffsample=puffpkdiffsample;
blratiodiffsample=puffpkdiffsample;
lateratiodiffsample=puffpkdiffsample;

for ibst=1:bootstraptime
%     idx=randperm(2*ncell);
    idx=randi(2*ncell,2*ncell,1);
    temp1=puffpkpool(idx(1:ncell));
    temp2=puffpkpool(idx(1+ncell:2*ncell));    
    puffpkdiffsample(ibst)=mean(temp1-temp2);
    
    temp1=neupkpool(idx(1:ncell));
    temp2=neupkpool(idx(1+ncell:2*ncell));    
    neupkdiffsample(ibst)=mean(temp1-temp2);
    
    temp1=puffblpool(idx(1:ncell));
    temp2=puffblpool(idx(1+ncell:2*ncell));
    puffbldiffsample(ibst)=mean(temp1-temp2);
    
    temp1=neublpool(idx(1:ncell));
    temp2=neublpool(idx(1+ncell:2*ncell));
    neubldiffsample(ibst)=mean(temp1-temp2);
    
    temp1=pufflatepool(idx(1:ncell));
    temp2=pufflatepool(idx(1+ncell:2*ncell));
    pufflatediffsample(ibst)=mean(temp1-temp2);
    
    temp1=neulatepool(idx(1:ncell));
    temp2=neulatepool(idx(1+ncell:2*ncell));
    neulatediffsample(ibst)=mean(temp1-temp2);
    
    temp1=pkenhancediffpool(idx(1:ncell));
    temp2=pkenhancediffpool(idx(1+ncell:2*ncell));
    pkenhancediffsample(ibst)=mean(temp1-temp2);
    
    temp1=blenhancediffpool(idx(1:ncell));
    temp2=blenhancediffpool(idx(1+ncell:2*ncell));
    blenhancediffsample(ibst)=mean(temp1-temp2);
    
    temp1=lateenhancediffpool(idx(1:ncell));
    temp2=lateenhancediffpool(idx(1+ncell:2*ncell));
    lateenhancediffsample(ibst)=mean(temp1-temp2);
    
    temp1=pkratiopool(idx(1:ncell));
    temp2=pkratiopool(idx(1+ncell:2*ncell));
    pkratiodiffsample(ibst)=mean(temp1-temp2);
    
    temp1=blratiopool(idx(1:ncell));
    temp2=blratiopool(idx(1+ncell:2*ncell));
    blratiodiffsample(ibst)=mean(temp1-temp2);
    
    temp1=lateratiopool(idx(1:ncell));
    temp2=lateratiopool(idx(1+ncell:2*ncell));
    lateratiodiffsample(ibst)=mean(temp1-temp2);
end
puffpkdiff=mean(conpuffpk)-mean(oripuffpk);
neupkdiff=mean(conneupk)-mean(orineupk);
puffbldiff=mean(conpuffbl)-mean(oripuffbl);
neubldiff=mean(conneubl)-mean(orineubl);
pufflatediff=mean(conpufflate)-mean(oripufflate);
neulatediff=mean(conneulate)-mean(orineulate);
pkenhancediff=mean(conpuffpk-oripuffpk)-mean(conneupk-orineupk);
blenhancediff=mean(conpuffbl-oripuffbl)-mean(conneubl-orineubl);
lateenhancediff=mean(conpufflate-oripufflate)-mean(conneulate-orineulate);
pkratiodiff=mean(pkpuffchangeratio)-mean(pkneuchangeratio);
blratiodiff=mean(blpuffchangeratio)-mean(blneuchangeratio);
lateratiodiff=mean(latepuffchangeratio)-mean(lateneuchangeratio);

puffpkzsc=(puffpkdiff-mean(puffpkdiffsample))/std(puffpkdiffsample);
neupkzsc=(neupkdiff-mean(neupkdiffsample))/std(neupkdiffsample);
puffblzsc=(puffbldiff-mean(puffbldiffsample))/std(puffbldiffsample);
neublzsc=(neubldiff-mean(neubldiffsample))/std(neubldiffsample);
pufflatezsc=(pufflatediff-mean(pufflatediffsample))/std(pufflatediffsample);
neulatezsc=(neulatediff-mean(neulatediffsample))/std(neulatediffsample);
pkenhancezsc=(pkenhancediff-mean(pkenhancediffsample))/std(pkenhancediffsample);
blenhancezsc=(blenhancediff-mean(blenhancediffsample))/std(blenhancediffsample);
lateenhancezsc=(lateenhancediff-mean(lateenhancediffsample))/std(lateenhancediffsample);
pkratiozsc=(pkratiodiff-mean(pkratiodiffsample))/std(pkratiodiffsample);
blratiozsc=(blratiodiff-mean(blratiodiffsample))/std(blratiodiffsample);
lateratiozsc=(lateratiodiff-mean(lateratiodiffsample))/std(lateratiodiffsample);

meanpuffpkratio=mean(pkpuffchangeratio)
meanneupkratio=mean(pkneuchangeratio)
p_puffpk=1-cdf('Normal',puffpkzsc,0,1)+cdf('Normal',-puffpkzsc,0,1);
p_neupk=1-cdf('Normal',neupkzsc,0,1)+cdf('Normal',-neupkzsc,0,1);
p_puffbl=1-cdf('Normal',puffblzsc,0,1)+cdf('Normal',-puffblzsc,0,1);
p_neubl=1-cdf('Normal',neublzsc,0,1)+cdf('Normal',-neublzsc,0,1);
p_pufflate=1-cdf('Normal',pufflatezsc,0,1)+cdf('Normal',-pufflatezsc,0,1);
p_neulate=1-cdf('Normal',neulatezsc,0,1)+cdf('Normal',-neulatezsc,0,1);
p_enhancepk=1-cdf('Normal',pkenhancezsc,0,1)+cdf('Normal',-pkenhancezsc,0,1);
p_enhancebl=1-cdf('Normal',blenhancezsc,0,1)+cdf('Normal',-blenhancezsc,0,1);
p_enhancelate=1-cdf('Normal',abs(lateenhancezsc),0,1)+cdf('Normal',-abs(lateenhancezsc),0,1);
p_pkratio=1-cdf('Normal',pkratiozsc,0,1)+cdf('Normal',-pkratiozsc,0,1);
p_blratio=1-cdf('Normal',blratiozsc,0,1)+cdf('Normal',-blratiozsc,0,1);
p_lateratio=1-cdf('Normal',abs(lateratiozsc),0,1)+cdf('Normal',-abs(lateratiozsc),0,1);

% fid=fopen('D:\result\Paper Prepare\part 1\psthstat.txt','w');
% fprintf(fid,[' peak puff diff ttest = %6d\r\n peak neu diff ttest = %6d\r\n bl puff diff ttest = %6d\r\n bl neu diff ttest = %6d\r\n' ...
%     ' late puff diff ttest = %6d\r\n late neu diff ttest = %6d\r\n'...
%     ' peak enhance diff ttest = %6d\r\n bl enhance diff ttest = %6d\r\n late enhance diff ttest = %6d\r\n\r\n'],...
%     ppuffpk_ttest,pneupk_ttest,ppuffbl_ttest,pneubl_ttest,ppufflate_ttest,pneulate_ttest,penchancepk_ttest,penchancebl_ttest,penchancelate_ttest)
% fprintf(fid,[' peak puff diff rank = %6d\r\n peak neu diff rank = %6d\r\n bl puff diff rank = %6d\r\n bl neu diff rank = %6d\r\n' ...
%     ' late puff diff rank = %6d\r\n late neu diff rank = %6d\r\n'...
%     ' peak enhance diff rank = %6d\r\n bl enhance diff rank = %6d\r\n late enhance diff rank = %6d\r\n\r\n'],...
%     ppuffpk_rank,pneupk_rank,ppuffbl_rank,pneubl_rank,ppufflate_rank,pneulate_rank,penchancepk_rank,penchancebl_rank,penchancelate_rank)
% fprintf(fid,[' peak puff diff bootstrap = %6d\r\n peak neu diff bootstrap = %6d\r\n bl puff diff bootstrap = %6d\r\n bl neu diff bootstrap = %6d\r\n' ...
%     ' late puff diff bootstrap = %6d\r\n late neu diff bootstrap = %6d\r\n'...
%     ' peak enhance diff bootstrap = %6d\r\n bl enhance diff bootstrap = %6d\r\n late enhance diff bootstrap = %6d\r\n\r\n'],...
%     p_puffpk,p_neupk,p_puffbl,p_neubl,p_pufflate,p_neulate,p_enhancepk,p_enhancebl,p_enhancelate)
% fclose(fid)

%%
% statistical test for d prime
secpkwin=(0.1-win(1))/bin+(0.016/bin+1:0.045/bin);
% secpkwin=316:345;
bootstraptime=1000;
npuffcell=size(puffdprimepool,2);
nneucell=size(neudprimepool,2);

pkdppool=[mean(puffdprimepool(pkwin,:),1) mean(neudprimepool(pkwin,:),1)];
secpkdppool=[mean(puffdprimepool(secpkwin,:),1) mean(neudprimepool(secpkwin,:),1)];
pkdpdiffsample=nan(bootstraptime,1);
secpkdpdiffsample=nan(bootstraptime,1);
for ibst=1:bootstraptime
    idx=randperm(npuffcell+nneucell);
    temp1=pkdppool(idx(1:npuffcell));
    temp2=pkdppool(idx(npuffcell+1:end));
    pkdpdiffsample(ibst)=mean(temp1)-mean(temp2);
    temp1=secpkdppool(idx(1:npuffcell));
    temp2=secpkdppool(idx(npuffcell+1:end));
    secpkdpdiffsample(ibst)=mean(temp1)-mean(temp2);    
end
pkdpdiff=mean(mean(puffdprimepool(pkwin,:)))-mean(mean(neudprimepool(pkwin,:)));
secpkdpdiff=mean(mean(puffdprimepool(secpkwin,:)))-mean(mean(neudprimepool(secpkwin,:)));

pkdpstd=(pkdpdiff-mean(pkdpdiffsample))/std(pkdpdiffsample);
secpkdpstd=(secpkdpdiff-mean(secpkdpdiffsample))/std(secpkdpdiffsample);

p_pkdp=1-cdf('Normal',pkdpstd,0,1)+cdf('Normal',-pkdpstd,0,1)
p_secpkdp=1-cdf('Normal',abs(secpkdpstd),0,1)+cdf('Normal',-abs(secpkdpstd),0,1)

% [hpuffpk,ppuffpk]=ttest(mean(puffdprimepool(pkwin,:),1));
% [hneupk,pneupk]=ttest(mean(neudprimepool(pkwin,:),1));
% [hpuffsecpk,ppuffsecpk]=ttest(mean(puffdprimepool(secpkwin,:),1));
% [hneusecpk,pneusecpk]=ttest(mean(neudprimepool(secpkwin,:),1));
% ppuffpk
% pneupk
% ppuffsecpk
% pneusecpk
%%
figure
xv=(-199:400)';
for mcnd=1:6
    subplot(2,3,mcnd)
    yv1=squeeze(mean(mergeprepsth(:,:,mcnd),2))/bin;
    yste1=std(mergeprepsth(:,:,mcnd),[],2)/sqrt(size(mergeprepsth,2))/bin;
    yv2=squeeze(mean(mergepostpsth(:,:,mcnd),2))/bin;
    yste2=std(mergepostpsth(:,:,mcnd),[],2)/sqrt(size(mergepostpsth,2))/bin;
    patch([xv; flipud(xv)],[yv1+yste1; flipud(yv1-yste1)],[0.75 1 0.75]);
    hold on
    h1=plot(xv,yv1,'g');
    patch([xv; flipud(xv)],[yv2+yste2; flipud(yv2-yste2)],[0.75 0.75 1]);
    h2=plot(xv,yv2,'b');
    if 3==mcnd
    legend([h1,h2],'Pre test','Post test')
    end
    mergeprebl=squeeze(mean(mergeprepsth(blwin,:,mcnd),1));
    mergeprepk=squeeze(mean(mergeprepsth(pkwin,:,mcnd),1));
    mergepostbl=squeeze(mean(mergepostpsth(blwin,:,mcnd),1));
    mergepostpk=squeeze(mean(mergepostpsth(pkwin,:,mcnd),1));
    [hcellbl(mcnd), pcellbl(mcnd)]=ttest2(mergeprebl,mergepostbl);
    [hcellpk(mcnd), pcellpk(mcnd)]=ttest2(mergeprepk,mergepostpk);
end
mergeprepuff=squeeze(mean(mergeprepuffpsth(pkwin,:),1));
mergepostpuff=squeeze(mean(mergepostpuffpsth(pkwin,:),1));
mergepreneu=squeeze(mean(mergepreneupsth(pkwin,:),1));
mergepostneu=squeeze(mean(mergepostneupsth(pkwin,:),1));
[hcellpuff,pcellpuff]=ttest2(mergeprepuff,mergepostpuff);
[hcellneu,pcellneu]=ttest2(mergepreneu,mergepostneu);
