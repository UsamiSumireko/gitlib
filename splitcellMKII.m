function splitcellMKII(destlist,LRlist,selecidx,plotlabel,savedir,savelabel)
% clear
% destlist=[{'D:\data extracted\160712'}];
% LRlist='r';
% Unitselec=[];
% plotlabel=0;
% savedir='D:\data extracted\160712';
% savelabel=0;
secpeakrmlabel=1;
secpeaklabel=0;
minusbaselinelabel=1;
normpsthlabel=true;
if isempty(selecidx)
    noselecidx=1;
end
% load('D:\result\170721_splitcell\test_prepostpsth\sameori_sele_vaildcell_pre-post_tp(1).mat');
%%
Splitmarker='_45degwin';
% savedir='I:\forprez_161008\SplitMKIII\';
Twidthwin=[0 90];
widthlabel='all';
bin=0.001;
win=[-0.1 0.5];
winidx=201:800; % winidx and win must be compatible
blwin=101:200;
pkwin=231:270;
tuningwin=(0.1-win(1))/bin+(0.030/bin+1:0.070/bin);
smoothlength=0.01/bin;
MergeLpcellStaORI=[];
MergeLpcellStaORIPost=[];
MergeLpcellCON=[];
MergeRpcellStaORI=[];
MergeRpcellStaORIPost=[];
MergeRpcellCON=[];
MergeMpcellStaORI=[];
MergeMpcellStaORIPost=[];
MergeMpcellCON=[];
MergeLpcelltuntrain=cell(6,1);
MergeLpcellcontrain=cell(6,1);
MergeRpcelltuntrain=cell(6,1);
MergeRpcellcontrain=cell(6,1);
MergeMpcelltuntrain=cell(6,1);
MergeMpcellcontrain=cell(6,1);
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
    load([destination '\XYvalidid.mat']);
    load([destination '\ORIvalidid.mat']);
    load([destination '\Oparam15.mat']);
    load([destination '\secpeakcellidx_both.mat']);
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
    if noselecidx
        selecidx{nday}=true(size(V1conblockunitpsth.psth,4),1);
    end
    if secpeakrmlabel
        if secpeaklabel
            selecidx{nday}=logical(selecidx{nday}) & secpeakcellidx;
        else
            selecidx{nday}=logical(selecidx{nday}) & ~secpeakcellidx;
        end
    else
        selecidx{nday}=logical(selecidx{nday});
    end
    Osigma=Oparam15{5}(1:12,:);
    Twidth=Osigma(ORIvalidid(1:12,:) & XYvalidid(1:12,:))*2.3548;
    Twidth=Twidth(selecidx{nday});
    [~,Twidthorder]=sort(Twidth);
    opref=Oparam15{1,4}(1:12,:);
    ORIpref=mod(opref(ORIvalidid(1:12,:) & XYvalidid(1:12,:)),180);
    ORIpref=ORIpref(selecidx{nday});
    
    if strcmp(widthlabel,'thin')
        Twidthcellidx=Twidth<=Twidth(Twidthorder(round(numel(Twidthorder)/2)));
%     elseif  strcmp(widthlabel,'mid')
%         Twidthcellidx= Twidth<=Twidth(Twidthorder(round(numel(Twidthorder)/3*2))) & Twidth>Twidth(Twidthorder(round(numel(Twidthorder)/3)));
    elseif strcmp(widthlabel,'thick')
        Twidthcellidx=Twidth>Twidth(Twidthorder(round(numel(Twidthorder)/2)));
    elseif strcmp(widthlabel,'all')
        Twidthcellidx=true(size(Twidth));
    else
        error('wrong widthlabel');
    end
    
      Lpref=ORIpref>112.5 & ORIpref<157.5 & Twidthcellidx;
      Rpref=ORIpref>22.5 & ORIpref<67.5 & Twidthcellidx;
      Mpref=ORIpref>67.5 & ORIpref<112.5 & Twidthcellidx;

%     Lpref=ORIpref>112.5 & ORIpref<157.5 & Twidth>Twidthwin(1) & Twidth<=Twidthwin(2);
%     Rpref=ORIpref>22.5 & ORIpref<67.5 & Twidth>Twidthwin(1) & Twidth<=Twidthwin(2);
%     Mpref=ORIpref>67.5 & ORIpref<112.5 & Twidth>Twidthwin(1) & Twidth<=Twidthwin(2);

    Lpcelltuntrain=cell(6,1);
    Lpcellcontrain=cell(6,1);
    Rpcelltuntrain=cell(6,1);
    Rpcellcontrain=cell(6,1);
    Mpcelltuntrain=cell(6,1);
    Mpcellcontrain=cell(6,1);
    
    nl=0;
    nr=0;
    nm=0;
    TrainLenth=(win(2)-win(1))/bin;
    traintemp=cell(6,1);
    mutemp=[];
    %Rpcelltuntrain等 每行是一个细胞，不是一个trial
    
    nanidx=isnan(V1conblockunitpsth.psth(winidx,:,:,:,:));
    connanvalid=sum(sum(sum(sum(nanidx))))==0;
    idx65535=65535==V1conblockunitpsth.psth(winidx,:,:,:,:);
    contbinvalid=sum(sum(sum(sum(idx65535))))==0;
    nanidx=isnan(V1oriblockunitpsth.psth(winidx,:,:,:,:));
    orinanvalid=sum(sum(sum(sum(nanidx))))==0;
    idx65535=65535==V1oriblockunitpsth.psth(winidx,:,:,:,:);
    oritbinvalid=sum(sum(sum(sum(idx65535))))==0;
    
    conpsth=permute(V1conblockunitpsth.psth(winidx,:,:,selecidx{nday},connanvalid & contbinvalid),[1 2 3 5 4]);
    oripsth=permute(V1oriblockunitpsth.psth(winidx,:,:,selecidx{nday},orinanvalid & oritbinvalid),[1 2 3 5 4]);
    if isempty(conpsth) || isempty(oripsth)
        continue % for nday=1:numel(destlist)    
    end
    for m=1:3
        mut=reshape(conpsth(tuningwin,:,m,:,:),numel(tuningwin),[],size(conpsth,5));
        if minusbaselinelabel
            bl=reshape(conpsth(blwin,:,m,:,:),100,[],size(conpsth,5));
            mutemp(:,CNDlist(m))=squeeze(mean(mean(mut)-mean(bl)));
        else
            mutemp(:,CNDlist(m))=squeeze(mean(mean(mut)));
        end
        traint=reshape(conpsth(:,:,m,:,:),size(conpsth,1),[],size(conpsth,5));
        traintemp{m}=squeeze(mean(traint,2))';
        Lpcellcontrain{CNDlist(m)}=traintemp{m}(Lpref,:);
        Rpcellcontrain{CNDlist(m)}=traintemp{m}(Rpref,:);
        Mpcellcontrain{CNDlist(m)}=traintemp{m}(Mpref,:);
        
        mut=reshape(conpsth(tuningwin,:,[m+3,m+6,m+9],:,:),numel(tuningwin),[],size(conpsth,5));
        if minusbaselinelabel
            bl=reshape(conpsth(blwin,:,[m+3,m+6,m+9],:,:),100,[],size(conpsth,5));
            mutemp(:,CNDlist(m+3))=squeeze(mean(mean(mut)-mean(bl)));
        else
        mutemp(:,CNDlist(m+3))=squeeze(mean(mean(mut)));
        end
        traint=reshape(conpsth(:,:,[m+3,m+6,m+9],:,:),size(conpsth,1),[],size(conpsth,5));
        traintemp{m+3}=squeeze(mean(traint,2))';
        Lpcellcontrain{CNDlist(m+3)}=traintemp{m+3}(Lpref,:);
        Rpcellcontrain{CNDlist(m+3)}=traintemp{m+3}(Rpref,:);
        Mpcellcontrain{CNDlist(m+3)}=traintemp{m+3}(Mpref,:);
    end
    LpcellCON=mutemp(Lpref,[2 3 4 5 6 1 2]);
    RpcellCON=mutemp(Rpref,[5 6 1 2 3 4 5]);
    MpcellCON=mutemp(Mpref,[1 2 3 4 5 6]);

    for m=1:6
        mut=oripsth(tuningwin,:,m,:,:);
        mut=reshape(mut,size(mut,1),[],size(mut,5));
        if minusbaselinelabel
            bl=oripsth(blwin,:,m,:,:);
            bl=reshape(bl,size(bl,1),[],size(bl,5));
            mutemp(:,m)=squeeze(mean(mean(mut)-mean(bl)));
        else
            mutemp(:,m)=squeeze(mean(mean(mut)));
        end
        traint=oripsth(:,:,m,:,:);
        traint=reshape(traint,size(traint,1),[],size(traint,5));
        traintemp{m}=squeeze(mean(traint,2))';
        Lpcelltuntrain{m}=traintemp{m}(Lpref,:);
        Rpcelltuntrain{m}=traintemp{m}(Rpref,:);
        Mpcelltuntrain{m}=traintemp{m}(Mpref,:);
    end
    LpcellStaORI=mutemp(Lpref,[2 3 4 5 6 1 2]);
    RpcellStaORI=mutemp(Rpref,[5 6 1 2 3 4 5]);
    MpcellStaORI=mutemp(Mpref,[1 2 3 4 5 6]);
   
    MergeLpcellStaORI=[MergeLpcellStaORI; LpcellStaORI];
    MergeLpcellCON=[MergeLpcellCON; LpcellCON];
    MergeRpcellStaORI=[MergeRpcellStaORI; RpcellStaORI];
    MergeRpcellCON=[MergeRpcellCON; RpcellCON];
    MergeMpcellStaORI=[MergeMpcellStaORI; MpcellStaORI];
    MergeMpcellCON=[MergeMpcellCON; MpcellCON];
    
    for m=1:6
        MergeLpcelltuntrain{m,1}=[MergeLpcelltuntrain{m,1}; Lpcelltuntrain{m}];
        MergeLpcellcontrain{m,1}=[MergeLpcellcontrain{m,1}; Lpcellcontrain{m}];
        MergeRpcelltuntrain{m,1}=[MergeRpcelltuntrain{m,1}; Rpcelltuntrain{m}];
        MergeRpcellcontrain{m,1}=[MergeRpcellcontrain{m,1}; Rpcellcontrain{m}];
        MergeMpcelltuntrain{m,1}=[MergeMpcelltuntrain{m,1}; Mpcelltuntrain{m}];
        MergeMpcellcontrain{m,1}=[MergeMpcellcontrain{m,1}; Mpcellcontrain{m}];
    end
    if savelabel
    save([destination '\LpcellStaORI' Splitmarker],'LpcellStaORI','-v7.3');
    save([destination '\Lpcelltuntrain' Splitmarker],'Lpcelltuntrain','-v7.3');
    save([destination '\LpcellCON' Splitmarker],'LpcellCON','-v7.3');
    save([destination '\Lpcellcontrain' Splitmarker],'Lpcellcontrain','-v7.3');
    save([destination '\RpcellStaORI' Splitmarker],'RpcellStaORI','-v7.3');
    save([destination '\Rpcelltuntrain' Splitmarker],'Rpcelltuntrain','-v7.3');
    save([destination '\RpcellCON' Splitmarker],'RpcellCON','-v7.3');
    save([destination '\Rpcellcontrain' Splitmarker],'Rpcellcontrain','-v7.3');
    save([destination '\MpcellStaORI' Splitmarker],'MpcellStaORI','-v7.3');
    save([destination '\Mpcelltuntrain' Splitmarker],'Mpcelltuntrain','-v7.3');
    save([destination '\MpcellCON' Splitmarker],'MpcellCON','-v7.3');
    save([destination '\Mpcellcontrain' Splitmarker],'Mpcellcontrain','-v7.3');
    end
    %%
    if plotlabel
        figure
        subplot(1,3,1)
        plot(-90:30:90,mean(LpcellCON),'r')
        hold on
        plot(-90:30:90,mean(LpcellStaORI),'g')
        subplot(1,3,2)
        plot(-90:30:90,mean(RpcellCON),'r')
        hold on
        plot(-90:30:90,mean(RpcellStaORI),'g')
        subplot(1,3,3)
        plot(15:30:165,mean(MpcellCON),'r')
        hold on
        plot(15:30:165,mean(MpcellStaORI),'g')
        figure
        for i=1:6
            subplot(2,3,i)
            plot(bin:bin:win(2),mean(Lpcelltuntrain{i}),'g')
            hold on
            plot(bin:bin:win(2),mean(Lpcellcontrain{i}),'r')
            hold on
            plot(bin:bin:win(2),mean(Rpcelltuntrain{i}),'b')
            hold on
            plot(bin:bin:win(2),mean(Rpcellcontrain{i}),'m')
            legend('Lpref Tun','Lpref Con','Rpref Tun','Rpref Con')
            title(num2str(i*30-15))
        end
    end
end

figure
% [Lmax,Lmaxidx]=max(MergeLpcellStaORI,[],2);
Lmax=MergeLpcellStaORI(:,4);
Lnormtun=MergeLpcellStaORI./repmat(Lmax,1,size(MergeLpcellStaORI,2));
Lnormtunste=std(Lnormtun)/sqrt(size(Lnormtun,1));
Lnormcon=MergeLpcellCON./repmat(Lmax,1,size(MergeLpcellCON,2));
Lnormconste=std(Lnormcon)/sqrt(size(Lnormcon,1));
% [Rmax,Rmaxidx]=max(MergeRpcellStaORI,[],2);
Rmax=MergeRpcellStaORI(:,4);
Rnormtun=MergeRpcellStaORI./repmat(Rmax,1,size(MergeRpcellStaORI,2));
Rnormtunste=std(Rnormtun)/sqrt(size(Rnormtun,1));
Rnormcon=MergeRpcellCON./repmat(Rmax,1,size(MergeRpcellCON,2));
Rnormconste=std(Rnormcon)/sqrt(size(Rnormcon,1));
[Mmax,Mmaxidx]=max(MergeMpcellStaORI,[],2);
Mnormtun=MergeMpcellStaORI./repmat(Mmax,1,size(MergeMpcellStaORI,2));
Mnormtunste=std(Mnormtun)/sqrt(size(Mnormtun,1));
Mnormcon=MergeMpcellCON./repmat(Mmax,1,size(MergeMpcellCON,2));
Mnormconste=std(Mnormcon)/sqrt(size(Mnormcon,1));
for m=1:7
    [h,p]=ttest(Lnormtun(:,m),Lnormcon(:,m));
    LTh(m)=h;
    LTp(m)=p;
    [h,p]=ttest(Rnormtun(:,m),Rnormcon(:,m));
    RTh(m)=h;
    RTp(m)=p;
end
for m=1:6
    [h,p]=ttest(Mnormtun(:,m),Mnormcon(:,m));
    MTh(m)=h;
    MTp(m)=p;
end
subplot(1,3,1)
plot(-90:30:90,mean(Lnormcon),'r')
hold on
errorbar(-90:30:90,mean(Lnormcon),Lnormconste,'rx')
plot(-90:30:90,mean(Lnormtun),'g')
errorbar(-90:30:90,mean(Lnormtun),Lnormtunste,'gx')
% plot(find(LTh==1)*30-120,repmat(150,1,sum(LTh==1)),'b*')
title(['Lpref cell N=' num2str(size(Lnormtun,1))])
xlabel('ORI (135°=0)')
ylabel('FR')
text
subplot(1,3,2)
plot(-90:30:90,mean(Rnormcon),'r')
hold on
errorbar(-90:30:90,mean(Rnormcon),Rnormconste,'rx')
plot(-90:30:90,mean(Rnormtun),'g')
errorbar(-90:30:90,mean(Rnormtun),Rnormtunste,'gx')
% plot(find(RTh==1)*30-120,repmat(150,1,sum(RTh==1)),'b*')
title(['Rpref cell N=' num2str(size(Rnormcon,1))])
xlabel('ORI (45°=0)')
ylabel('FR')
subplot(1,3,3)
plot(15:30:165,mean(Mnormcon),'r')
hold on
errorbar(15:30:165,mean(Mnormcon),Mnormconste,'rx')
plot(15:30:165,mean(Mnormtun),'g')
errorbar(15:30:165,mean(Mnormtun),Mnormtunste,'gx')
% plot(find(MTh==1)*30-15,repmat(150,1,sum(MTh==1)),'b*')
title(['Mpref cell N=' num2str(size(Mnormtun,1))])
xlabel('ORI ')
ylabel('FR')

% stat=struct();
% for m=1:6
%     for n=1:100
%         [h1,p1]=ttest(MergeLpcelltuntrain{m,1}(:,n),MergeLpcellcontrain{m,1}(:,n));
%         [h2,p2]=ttest(MergeRpcelltuntrain{m,1}(:,n),MergeRpcellcontrain{m,1}(:,n));
%         [h3,p3]=ttest(MergeMpcelltuntrain{m,1}(:,n),MergeMpcellcontrain{m,1}(:,n));
%         stat.Lh(m,n)=h1;
%         stat.Lp(m,n)=p1;
%         stat.Rh(m,n)=h2;
%         stat.Rp(m,n)=p2;
%         stat.Mh(m,n)=h3;
%         stat.Mp(m,n)=p3;
%     end
% end

% figure
% for i=1:6
%     subplot(2,3,i)
%     steLtun=std(MergeLpcelltuntrain{i})/sqrt(size(MergeLpcelltuntrain{i},1));
%     patch([bin:bin:win(2) fliplr(bin:bin:win(2))],[mean(MergeLpcelltuntrain{i})/bin-steLtun/bin fliplr(mean(MergeLpcelltuntrain{i})/bin+steLtun/bin)],[0.8 1 0.8])
%     hold on
%     h1=plot(bin:bin:win(2),mean(MergeLpcelltuntrain{i})/bin,'g');
%     hold on
%     steLcon=std(MergeLpcellcontrain{i})/sqrt(size(MergeLpcellcontrain{i},1));
%     patch([bin:bin:win(2) fliplr(bin:bin:win(2))],[mean(MergeLpcellcontrain{i})/bin-steLcon/bin fliplr(mean(MergeLpcellcontrain{i})/bin+steLcon/bin)],[1 0.8 0.8])
%     hold on
%     h2=plot(bin:bin:win(2),mean(MergeLpcellcontrain{i})/bin,'r');
%     hold on
%     steRtun=std(MergeRpcelltuntrain{i})/sqrt(size(MergeRpcelltuntrain{i},1));
%     patch([bin:bin:win(2) fliplr(bin:bin:win(2))],[mean(MergeRpcelltuntrain{i})/bin-steRtun/bin fliplr(mean(MergeRpcelltuntrain{i})/bin+steRtun/bin)],[0.8 0.8 1])
%     hold on
%     h3=plot(bin:bin:win(2),mean(MergeRpcelltuntrain{i})/bin,'b');
%     hold on
%     steRcon=std(MergeRpcellcontrain{i})/sqrt(size(MergeRpcellcontrain{i},1));
%     patch([bin:bin:win(2) fliplr(bin:bin:win(2))],[mean(MergeRpcellcontrain{i})/bin-steRcon/bin fliplr(mean(MergeRpcellcontrain{i})/bin+steRcon/bin)],[1 0.8 1])
%     hold on
%     h4=plot(bin:bin:win(2),mean(MergeRpcellcontrain{i})/bin,'m');
% %     h5=plot(bin*(find(stat.Lh(i,:)==1)),repmat(40,1,sum(stat.Lh(i,:)==1)),'r*');
% %     h6=plot(bin*(find(stat.Rh(i,:)==1)),repmat(60,1,sum(stat.Rh(i,:)==1)),'m*');
%     legend([h1,h2,h3,h4],'Lpref Tun','Lpref Con','Rpref Tun','Rpref Con')
%     title(num2str(i*30-15))
%     xlabel('Time (ms)')
%     xlabel('FR')
% end
% blwin=101:200;
Lpcellnormbase=MergeLpcelltuntrain{5};
Rpcellnormbase=MergeRpcelltuntrain{2};
for m=1:6
    if normpsthlabel
        MergeLpcellcontrain{m}=normlizepsth(MergeLpcellcontrain{m},blwin,Lpcellnormbase);
        MergeRpcellcontrain{m}=normlizepsth(MergeRpcellcontrain{m},blwin,Rpcellnormbase);
%         MergeMpcellcontrain{m}=normlizepsth(MergeMpcellcontrain{m},blwin,MergeMpcelltuntrain{m});
        MergeLpcelltuntrain{m}=normlizepsth(MergeLpcelltuntrain{m},blwin,Lpcellnormbase);
        MergeRpcelltuntrain{m}=normlizepsth(MergeRpcelltuntrain{m},blwin,Rpcellnormbase);
%         MergeMpcelltuntrain{m}=normlizepsth(MergeMpcelltuntrain{m},blwin,MergeMpcelltuntrain{m});

    end
    MergeLpcelltuntrain{m}=conv2(MergeLpcelltuntrain{m},ones(1,smoothlength),'same')/smoothlength;    
    MergeLpcellcontrain{m}=conv2(MergeLpcellcontrain{m},ones(1,smoothlength),'same')/smoothlength;   
    MergeRpcelltuntrain{m}=conv2(MergeRpcelltuntrain{m},ones(1,smoothlength),'same')/smoothlength;  
    MergeRpcellcontrain{m}=conv2(MergeRpcellcontrain{m},ones(1,smoothlength),'same')/smoothlength;
%     MergeMpcelltuntrain{m}=conv2(MergeMpcelltuntrain{m},ones(1,smoothlength),'same')/smoothlength;
%     MergeMpcellcontrain{m}=conv2(MergeMpcellcontrain{m},ones(1,smoothlength),'same')/smoothlength;
end
if savelabel
% save([savedir 'stat'],'stat','-v7.3');
save([savedir 'MergeLpcellCON'],'MergeLpcellCON','-v7.3');
save([savedir 'MergeLpcellStaORI'],'MergeLpcellStaORI','-v7.3');
save([savedir 'MergeRpcellCON'],'MergeRpcellCON','-v7.3');
save([savedir 'MergeRpcellStaORI'],'MergeRpcellStaORI','-v7.3');
save([savedir 'MergeMpcellCON'],'MergeMpcellCON','-v7.3');
save([savedir 'MergeMpcellStaORI'],'MergeMpcellStaORI','-v7.3');
save([savedir 'MergeLpcelltuntrain'],'MergeLpcelltuntrain','-v7.3');
save([savedir 'MergeLpcellcontrain'],'MergeLpcellcontrain','-v7.3');
save([savedir 'MergeRpcelltuntrain'],'MergeRpcelltuntrain','-v7.3');
save([savedir 'MergeRpcellcontrain'],'MergeRpcellcontrain','-v7.3');
save([savedir 'MergeMpcelltuntrain'],'MergeMpcelltuntrain','-v7.3');
save([savedir 'MergeMpcellcontrain'],'MergeMpcellcontrain','-v7.3');
end

    LcLconpsth=[];
    RcLconpsth=[];
    McLconpsth=[];
    LcLtunpsth=[];
    RcLtunpsth=[];
    McLtunpsth=[];
    
    LcRconpsth=[];
    RcRconpsth=[];
    McRconpsth=[];
    LcRtunpsth=[];
    RcRtunpsth=[];
    McRtunpsth=[];
% for m=1:3
for m=2
    LcLconpsth=[LcLconpsth; MergeLpcellcontrain{m+3}];
    RcLconpsth=[RcLconpsth; MergeRpcellcontrain{m+3}];
    LcLtunpsth=[LcLtunpsth; MergeLpcelltuntrain{m+3}];
    RcLtunpsth=[RcLtunpsth; MergeRpcelltuntrain{m+3}];
    
    LcRconpsth=[LcRconpsth; MergeLpcellcontrain{m}];
    RcRconpsth=[RcRconpsth; MergeRpcellcontrain{m}];
    LcRtunpsth=[LcRtunpsth; MergeLpcelltuntrain{m}];
    RcRtunpsth=[RcRtunpsth; MergeRpcelltuntrain{m}];
end
steLLcon=std(LcLconpsth)/sqrt(size(LcLconpsth,1));
steRLcon=std(RcLconpsth)/sqrt(size(RcLconpsth,1));
steLLtun=std(LcLtunpsth)/sqrt(size(LcLtunpsth,1));
steRLtun=std(RcLtunpsth)/sqrt(size(RcLtunpsth,1));

steLRcon=std(LcRconpsth)/sqrt(size(LcRconpsth,1));
steRRcon=std(RcRconpsth)/sqrt(size(RcRconpsth,1));
steLRtun=std(LcRtunpsth)/sqrt(size(LcRtunpsth,1));
steRRtun=std(RcRtunpsth)/sqrt(size(RcRtunpsth,1));

figure
subplot(1,2,1)
xv=(win(1)+bin:bin:win(2))-0.1; 
xv=xv/bin;
yv=mean(LcLconpsth(:,:));
yste=steLLcon;
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[1 0.8 0.8]);
hold on
h1=plot(xv,yv,'r');

yv=mean(LcLtunpsth(:,:));
yste=steLLtun;
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[0.8 1 0.8]);
hold on
h2=plot(xv,yv,'g');

yv=mean(RcLconpsth(:,:));
yste=steRLcon;
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[1 0.8 1]);
hold on
h3=plot(xv,yv,'m');

yv=mean(RcLtunpsth(:,:));
yste=steRLtun;
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[0.8 0.8 1]);
hold on
h4=plot(xv,yv,'b');
xlim([-100 300])
% ylim([0 180])
h5=legend([h1,h2,h3,h4],'Lori Lpref Con','Lpref tun','Rpref con','Rpref tun');
set(h5,'box','off')

subplot(1,2,2)
xv=(win(1)+bin:bin:win(2))-0.1;
xv=xv/bin;
yv=mean(LcRconpsth(:,:));
yste=steLRcon;
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[1 0.8 0.8]);
hold on
h1=plot(xv,yv,'r');

yv=mean(LcRtunpsth(:,:));
yste=steLRtun;
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[0.8 1 0.8]);
hold on
h2=plot(xv,yv,'g');

yv=mean(RcRconpsth(:,:));
yste=steRRcon;
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[1 0.8 1]);
hold on
h3=plot(xv,yv,'m');

yv=mean(RcRtunpsth(:,:));
yste=steRRtun;
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[0.8 0.8 1]);
hold on
h4=plot(xv,yv,'b');
xlim([-100 300])
% ylim([0 180])
h5=legend([h1,h2,h3,h4],'Rori Lpref Con','Lpref tun','Rpref con','Rpref tun');
set(h5,'box','off');

% d-prime
figure
subplot(1,2,1)
xv=(win(1)+bin:bin:win(2))-0.1; 
xv=xv/bin;
yv1=(mean(LcLconpsth(:,:))-mean(LcLtunpsth(:,:)))./sqrt((std(LcLconpsth(:,:)).^2+std(LcLtunpsth(:,:)).^2)/2);
h1=plot(xv,yv1,'r');

yv2=(mean(RcLconpsth(:,:))-mean(RcLtunpsth(:,:)))./sqrt((std(RcLconpsth(:,:)).^2+std(RcLtunpsth(:,:)).^2)/2);
hold on
h2=plot(xv,yv2,'g');

xlim([-100 300])
% ylim([0 180])
h3=legend([h1,h2],'Lori Lpref d-prime','Rpref d-prime');
set(h3,'box','off')

subplot(1,2,2)
xv=(win(1)+bin:bin:win(2))-0.1; 
xv=xv/bin;
yv1=(mean(LcRconpsth(:,:))-mean(LcRtunpsth(:,:)))./sqrt((std(LcRconpsth(:,:)).^2+std(LcRtunpsth(:,:)).^2)/2);
h1=plot(xv,yv1,'r');

yv2=(mean(RcRconpsth(:,:))-mean(RcRtunpsth(:,:)))./sqrt((std(RcRconpsth(:,:)).^2+std(RcRtunpsth(:,:)).^2)/2);
hold on
h2=plot(xv,yv2,'g');
xlim([-100 300])
% ylim([0 180])
h3=legend([h1,h2],'rori Lpref d-prime','Rpref d-prime');
set(h3,'box','off')

peakwin=231:270;
figure
subplot(1,2,1)
yv=mean(LcLconpsth(:,peakwin),2)./mean(LcLtunpsth(:,peakwin),2);
yv=yv(yv>0 & yv<10);
yste=std(yv)/sqrt(size(yv,1));
errorbar(1,mean(yv),yste);
h1=bar(1,mean(yv),'r');
hold on
yv=mean(RcLconpsth(:,peakwin),2)./mean(RcLtunpsth(:,peakwin),2);
yv=yv(yv>0 & yv<10);
yste=std(yv)/sqrt(size(yv,1));
errorbar(3,mean(yv),yste);
h2=bar(3,mean(yv),'b');
h3=legend([h1,h2],'Lori Lpref','Rpref');
ylim([0.9 1.5]);

subplot(1,2,2)
yv=mean(LcRconpsth(:,peakwin),2)./mean(LcRtunpsth(:,peakwin),2);
yv=yv(yv>0 & yv<10);
yste=std(yv)/sqrt(size(yv,1));
errorbar(1,mean(yv),yste);
h1=bar(1,mean(yv),'r');
hold on
yv=mean(RcRconpsth(:,peakwin),2)./mean(RcRtunpsth(:,peakwin),2);
yv=yv(yv>0 & yv<10);
yste=std(yv)/sqrt(size(yv,1));
errorbar(3,mean(yv),yste);
h2=bar(3,mean(yv),'b');
h3=legend([h1,h2],'Rori Lpref','Rpref');
ylim([0.9 1.5]);
%%
%statistical test for lp and rp cell psth
pkwin=231:270;


LcLdiff=mean(mean(LcLconpsth(:,pkwin)-LcLtunpsth(:,pkwin)));
LcRdiff=mean(mean(LcRconpsth(:,pkwin)-LcRtunpsth(:,pkwin)));
RcLdiff=mean(mean(RcLconpsth(:,pkwin)-RcLtunpsth(:,pkwin)));
RcRdiff=mean(mean(RcRconpsth(:,pkwin)-RcRtunpsth(:,pkwin)));

LcLpool=[mean(LcLconpsth(:,pkwin),2); mean(LcLtunpsth(:,pkwin),2)];
LcRpool=[mean(LcRconpsth(:,pkwin),2); mean(LcRtunpsth(:,pkwin),2)];
RcLpool=[mean(RcLconpsth(:,pkwin),2); mean(RcLtunpsth(:,pkwin),2)];
RcRpool=[mean(RcRconpsth(:,pkwin),2); mean(RcRtunpsth(:,pkwin),2)];

nbootstrap=1000;
LcLdiffpool=nan(nbootstrap,1);
LcRdiffpool=LcLdiffpool;
RcLdiffpool=LcLdiffpool;
RcRdiffpool=LcLdiffpool;
nLcell=size(LcLconpsth,1);
nRcell=size(RcLconpsth,1);
for ibst=1:nbootstrap
    Lidx=randperm(2*nLcell);
    temp1=LcLpool(Lidx(1:nLcell));
    temp2=LcLpool(Lidx(nLcell+1:2*nLcell));
    LcLdiffpool(ibst)=mean(temp1-temp2);
    temp1=LcRpool(Lidx(1:nLcell));
    temp2=LcRpool(Lidx(nLcell+1:2*nLcell));
    LcRdiffpool(ibst)=mean(temp1-temp2);
    
    Ridx=randperm(2*nRcell);
    temp1=RcLpool(Ridx(1:nRcell));
    temp2=RcLpool(Ridx(nRcell+1:2*nRcell));
    RcLdiffpool(ibst)=mean(temp1-temp2);
    temp1=RcRpool(Ridx(1:nRcell));
    temp2=RcRpool(Ridx(nRcell+1:2*nRcell));
    RcRdiffpool(ibst)=mean(temp1-temp2);    
end
LcLstd=(LcLdiff-mean(LcLdiffpool))/std(LcLdiffpool);
LcRstd=(LcRdiff-mean(LcRdiffpool))/std(LcRdiffpool);
RcLstd=(RcLdiff-mean(RcLdiffpool))/std(RcLdiffpool);
RcRstd=(RcRdiff-mean(RcRdiffpool))/std(RcRdiffpool);

p_LcL=1-cdf('Norm',abs(LcLstd),0,1)+cdf('Norm',-abs(LcLstd),0,1)
p_LcR=1-cdf('Norm',abs(LcRstd),0,1)+cdf('Norm',-abs(LcRstd),0,1)
p_RcL=1-cdf('Norm',abs(RcLstd),0,1)+cdf('Norm',-abs(RcLstd),0,1)
p_RcR=1-cdf('Norm',abs(RcRstd),0,1)+cdf('Norm',-abs(RcRstd),0,1)

function nmpsth=normlizepsth(psth,blwin,blpsth)
BL1=repmat(mean(psth(:,blwin),2),1,size(psth,2));
BL2=repmat(mean(blpsth(:,blwin),2),1,size(blpsth,2));
maxbin=repmat(max(blpsth-BL2,[],2),1,size(psth,2));
nmpsth=(psth-BL1)./maxbin;