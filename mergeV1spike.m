% old ver, use old data structure, for 5-d data array, use MKII version
clear
folderlist={
%     'I:\160712';
%     'I:\160713';
%     'I:\160714';
%     'I:\160715';
%     'I:\160716';
%     'I:\160719';
%     'I:\160720';
%     'I:\160721';
%     'I:\160722';
    
    'G:\MS\160624';
    'G:\MS\160625';
    'G:\MS\160627';
    'G:\MS\160628';
    
%     'I:\170426'
%     'I:\170427'
%     'I:\170429'
%     'I:\170430'
%     'I:\170502'
%     'I:\170503'
    
    'I:\170516'
    'I:\170518'
    'I:\170519'
    'I:\170522'
    };

LRlist=['rrrrrrrrr' 'llll' 'rrrrrr' 'llll'];
selec='ori';
bin=0.001;
win=[-0.1 0.5];
% for i=1:numel(folderlist)
%     i,
% %     delete([folderlist{i} '\15deg_1msbin_MergetunTrain.mat']);
% %     delete([folderlist{i} '\1msbin_MergeconTrain.mat']);
%     mergespiketrain(folderlist{i},LRlist(i),selec,bin,win);
% end
blwin=-win(1)/bin+1:(0.1-win(1))/bin;
pkwin=(0.1-win(1))/bin+(0.030/bin+1:0.070/bin);
smoothlength=0.01/bin;
allconpsth=cell(6,1);
alloripsth=cell(6,1);
allpuffpsth=cell(2,1);
allneupsth=cell(2,1);
for iday=1:numel(folderlist)
    iday
    if strcmpi(LRlist(iday),'r')
        cndlist=[1 2 3 4 5 6];
    elseif strcmpi(LRlist(iday),'l')
        cndlist=[4 5 6 1 2 3];
    else
        error('wrong LR')
    end
    if exist([folderlist{iday} '\1msbin_rawunitORIPSTH.mat'],'file')
        load([folderlist{iday} '\1msbin_rawunitORIPSTH.mat'])
%     elseif exist([folderlist{iday} '\rawunitORIPSTH.mat'],'file')   
%         load([folderlist{iday} '\rawunitORIPSTH.mat'])
    else
        error('no rawunitORIPSTH')
    end
    
    if exist([folderlist{iday} '\1msbin_rawunitCONPSTH.mat'],'file')
        load([folderlist{iday} '\1msbin_rawunitCONPSTH.mat'])
%     elseif exist([folderlist{iday} '\rawunitCONPSTH.mat'],'file')
%         load([folderlist{iday} '\rawunitCONPSTH.mat'])
    else
        error('no rawunitORIPSTH')
    end
    load([folderlist{iday} '\XYvalidid.mat'])
    load([folderlist{iday} '\ORIvalidid.mat'])
    if strcmpi(selec,'ori')
        validid=ORIvalidid;
    elseif strcmpi(selec,'xy')
        validid=XYvalidid;
    else
        error('wrong validid')
    end   
    for j=1:12
        for k=1:8
            if validid(j,k)
                for m=1:6
                    assert((win(2)-win(1))/bin==size(rawunitCONPSTH{j,k}{m},2))
                    assert((win(2)-win(1))/bin==size(rawunitORIPSTH{j,k}{m},2))
                    allconpsth{m}=[allconpsth{m}; mean(rawunitCONPSTH{j,k}{m})];
                    alloripsth{m}=[alloripsth{m}; mean(rawunitORIPSTH{j,k}{m})];
                end
                psth1=mean([rawunitCONPSTH{j,k}{cndlist(1)}; rawunitCONPSTH{j,k}{cndlist(2)}; rawunitCONPSTH{j,k}{cndlist(3)}]);
                psth2=mean([rawunitORIPSTH{j,k}{cndlist(1)}; rawunitORIPSTH{j,k}{cndlist(2)}; rawunitORIPSTH{j,k}{cndlist(3)}]);
                allpuffpsth{1}=[allpuffpsth{1}; psth1];
                allpuffpsth{2}=[allpuffpsth{2}; psth2];
                psth1=mean([rawunitCONPSTH{j,k}{cndlist(4)}; rawunitCONPSTH{j,k}{cndlist(5)}; rawunitCONPSTH{j,k}{cndlist(6)}]);
                psth2=mean([rawunitORIPSTH{j,k}{cndlist(4)}; rawunitORIPSTH{j,k}{cndlist(5)}; rawunitORIPSTH{j,k}{cndlist(6)}]);
                allneupsth{1}=[allneupsth{1}; psth1];
                allneupsth{2}=[allneupsth{2}; psth2];
            end
        end
    end
end
% allpuffpsth 是con block {1}和 ori block{2}中吹气朝向的psth
% allneupsth 是con block {1}和 ori block{2}中中性朝向的psth
% for m=1:6
%     for n=1:100
%     [h,p,ci,stats]=ttest(allconpsth{m}(:,n),alloripsth{m}(:,n));
%     H(n)=h;P(n)=p;
%     [h2,p2,ci2,stats2]=ttest2(allconpsth{m}(:,n),alloripsth{m}(:,n));
%     H2(n)=h2;P2(n)=p2;
%     end
% end
% for n=1:100
%     [h,p,ci,stats]=ttest(allpuffpsth{1}(:,n),allpuffpsth{2}(:,n));
%     HH(n)=h;PP(n)=p;
%     [h2,p2,ci2,stats2]=ttest2(allpuffpsth{1}(:,n),allpuffpsth{2}(:,n));
%     HH2(n)=h2;PP2(n)=p2;
%     [h,p,ci,stats]=ttest(allneupsth{1}(:,n),allneupsth{2}(:,n));
%     Hh(n)=h;Pp(n)=p;
%     [h2,p2,ci2,stats2]=ttest2(allneupsth{1}(:,n),allneupsth{2}(:,n));
%     Hh2(n)=h2;Pp2(n)=p2;
% end
% figure
% for m=1:6
%     conste=std(allconpsth{m}/bin)/sqrt(size(allconpsth{m},1));
%     oriste=std(alloripsth{m}/bin)/sqrt(size(alloripsth{m},1));
%     subplot(2,3,m)
%     fill([(1:100)*bin fliplr((1:100)*bin)],[mean(allconpsth{m}/bin)+conste fliplr(mean(allconpsth{m}/bin)-conste)],'r')
%     hold on
%     plot((1:100)*bin,mean(allconpsth{m}/bin),'r')
%     plot(find(H==1)*bin,repmat(10,1,sum(H==1)),'ob')
%     plot(find(H2==1)*bin,repmat(20,1,sum(H2==1)),'*b')
%     fill([(1:100)*bin fliplr((1:100)*bin)],[mean(alloripsth{m}/bin)+oriste fliplr(mean(alloripsth{m}/bin)-oriste)],'g')
%     plot((1:100)*bin,mean(alloripsth{m}/bin),'g')
% end

figure

allpuffpsth{1}=conv2(allpuffpsth{1},ones(1,smoothlength),'same')/smoothlength;
allpuffpsth{2}=conv2(allpuffpsth{2},ones(1,smoothlength),'same')/smoothlength;
allneupsth{1}=conv2(allneupsth{1},ones(1,smoothlength),'same')/smoothlength;
allneupsth{2}=conv2(allneupsth{2},ones(1,smoothlength),'same')/smoothlength;
puffste1=std(allpuffpsth{1}/bin)/sqrt(size(allpuffpsth{1},1));
puffste2=std(allpuffpsth{2}/bin)/sqrt(size(allpuffpsth{2},1));
neuste1=std(allneupsth{1}/bin)/sqrt(size(allneupsth{1},1));
neuste2=std(allneupsth{2}/bin)/sqrt(size(allneupsth{2},1));
subplot(2,1,1)
vx=(win(1)+bin:bin:win(2))-0.1;
vy=mean(allpuffpsth{1}/bin);
patch([vx fliplr(vx) vx(1)],[vy+puffste1 fliplr(vy-puffste1) vy(1)+puffste1(1)],'r')
hold on
plot(vx,vy,'r')
vy=mean(allpuffpsth{2}/bin);
patch([vx fliplr(vx) vx(1)]',[vy+puffste2 fliplr(vy-puffste2) vy(1)+puffste2(1)]','g')
plot(vx,vy,'g')
xlim([0 0.3])
ylim([0 120])
% plot(find(HH==1)*bin,repmat(10,1,sum(HH==1)),'ob')
% plot(find(HH2==1)*bin,repmat(20,1,sum(HH2==1)),'*b')

subplot(2,1,2)
vy=mean(allneupsth{1}/bin);
patch([vx fliplr(vx) vx(1)],[vy+neuste1 fliplr(vy-neuste1) vy(1)+neuste1(1)],'r')
hold on
plot(vx,vy,'r')
vy=mean(allneupsth{2}/bin);
patch([vx fliplr(vx) vx(1)],[vy+neuste2 fliplr(vy-neuste2) vy(1)+neuste2(1)],'g')
plot(vx,mean(allneupsth{2}/bin),'g')
xlim([0 0.3])
ylim([0 120])

puffbldiff=mean(mean(allpuffpsth{1}(:,blwin)-allpuffpsth{2}(:,blwin)))/bin;
puffblste=std(mean(allpuffpsth{1}(:,blwin)-allpuffpsth{2}(:,blwin),2))/sqrt(size(allpuffpsth{1},1))/bin;
neubldiff=mean(mean(allneupsth{1}(:,blwin)-allneupsth{2}(:,blwin)))/bin;
neublste=std(mean(allneupsth{1}(:,blwin)-allneupsth{2}(:,blwin),2))/sqrt(size(allneupsth{1},1))/bin;
puffpkdiff=mean(mean(allpuffpsth{1}(:,pkwin)-allpuffpsth{2}(:,pkwin)))/bin;
puffpkste=std(mean(allpuffpsth{1}(:,pkwin)-allpuffpsth{2}(:,pkwin),2))/sqrt(size(allpuffpsth{1},1))/bin;
neupkdiff=mean(mean(allneupsth{1}(:,pkwin)-allneupsth{2}(:,pkwin)))/bin;
neupkste=std(mean(allneupsth{1}(:,pkwin)-allneupsth{2}(:,pkwin),2))/sqrt(size(allneupsth{1},1))/bin;
figure 
bar([1 2 6 7],[puffbldiff puffpkdiff neubldiff neupkdiff]);
hold on
errorbar([1 2 6 7],[puffbldiff puffpkdiff neubldiff neupkdiff],[puffblste puffpkste neublste neupkste],'k*')



% plot(find(Hh==1)*bin,repmat(10,1,sum(Hh==1)),'ob')
% plot(find(Hh2==1)*bin,repmat(20,1,sum(Hh2==1)),'*b')