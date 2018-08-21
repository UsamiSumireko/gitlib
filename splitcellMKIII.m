function splitcellMKIII(destlist,LRlist,selecidx,plotlabel,savedir,savelabel)
% only to get normalized peak=1 splitcell psth

% destlist=[{'D:\data extracted\160712'}];
% LRlist='r';
% Unitselec=[];
% plotlabel=0;
% savedir='D:\data extracted\160712';
% savelabel=0;
secpeakrmlabel=1;
secpeaklabel=0;
minusbaselinelabel=0;
normpsthlabel=true;
% load('D:\result\170721_splitcell\test_prepostpsth\sameori_sele_vaildcell_pre-post_tp(1).mat');
%%
Splitmarker='_45degwin';
% savedir='I:\forprez_161008\SplitMKIII\';
Twidthwin=[0 90];
widthlabel='all';
bin=0.001;
win=[-0.1 0.5];
winidx=1:600;
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
    load([destination '\secpeakcellidx.mat']);
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
    Twidth=Osigma(ORIvalidid(1:12,:))*2.3548;
    Twidth=Twidth(selecidx{nday});
    [~,Twidthorder]=sort(Twidth);
    opref=Oparam15{1,4}(1:12,:);
    ORIpref=mod(opref(ORIvalidid(1:12,:)),180);
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
    conpsth=permute(V1conblockunitpsth.psth(:,:,:,selecidx{nday},:),[1 2 3 5 4]);
    for m=1:3
        mut=reshape(conpsth(tuningwin,:,m,:,:),numel(tuningwin),[],size(conpsth,5));
        if minusbaselinelabel
            bl=reshape(conpsth(101:200,:,m,:,:),100,[],size(conpsth,5));
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
            bl=reshape(conpsth(101:200,:,[m+3,m+6,m+9],:,:),100,[],size(conpsth,5));
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
        mut=permute(V1oriblockunitpsth.psth(tuningwin,:,m,selecidx{nday},:),[1 2 3 5 4]);
        mut=reshape(mut,size(mut,1),[],size(mut,5));
        if minusbaselinelabel
            bl=permute(V1oriblockunitpsth.psth(101:200,:,m,selecidx{nday},:),[1 2 3 5 4]);
            bl=reshape(bl,size(bl,1),[],size(bl,5));
            mutemp(:,m)=squeeze(mean(mean(mut)-mean(bl)));
        else
            mutemp(:,m)=squeeze(mean(mean(mut)));
        end
        traint=permute(V1oriblockunitpsth.psth(:,:,m,selecidx{nday},:),[1 2 3 5 4]);
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
blwin=101:200;
for m=1:6
    MergeLpcelltuntrain{m}=conv2(MergeLpcelltuntrain{m},ones(1,smoothlength),'same')/smoothlength;
    MergeLpcellcontrain{m}=conv2(MergeLpcellcontrain{m},ones(1,smoothlength),'same')/smoothlength;
    MergeRpcelltuntrain{m}=conv2(MergeRpcelltuntrain{m},ones(1,smoothlength),'same')/smoothlength;
    MergeRpcellcontrain{m}=conv2(MergeRpcellcontrain{m},ones(1,smoothlength),'same')/smoothlength;
    %     MergeMpcelltuntrain{m}=conv2(MergeMpcelltuntrain{m},ones(1,smoothlength),'same')/smoothlength;
    %     MergeMpcellcontrain{m}=conv2(MergeMpcellcontrain{m},ones(1,smoothlength),'same')/smoothlength;
end
if normpsthlabel
    blpsthLpc=MergeLpcelltuntrain{5};
    blpsthRpc=MergeRpcelltuntrain{2};
    for m=1:6
        MergeLpcellcontrain{m}=normlizepsth(MergeLpcellcontrain{m},blwin,blpsthLpc);
        MergeRpcellcontrain{m}=normlizepsth(MergeRpcellcontrain{m},blwin,blpsthRpc);
        %         MergeMpcellcontrain{m}=normlizepsth(MergeMpcellcontrain{m},blwin,MergeMpcelltuntrain{m});
        MergeLpcelltuntrain{m}=normlizepsth(MergeLpcelltuntrain{m},blwin,blpsthLpc);
        MergeRpcelltuntrain{m}=normlizepsth(MergeRpcelltuntrain{m},blwin,blpsthRpc);
        %         MergeMpcelltuntrain{m}=normlizepsth(MergeMpcelltuntrain{m},blwin,MergeMpcelltuntrain{m});
    end
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
winidx=1:500;
subplot(1,2,1)
% xv=(win(1)+bin:bin:win(2))-0.1; 
% xv=xv/bin;
xv=-199:300;
yv=mean(LcLconpsth(:,winidx));
yste=steLLcon(winidx);
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[1 0.8 0.8]);
hold on
h1=plot(xv,yv,'r');

yv=mean(LcLtunpsth(:,winidx));
yste=steLLtun(winidx);
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[0.8 1 0.8]);
hold on
h2=plot(xv,yv,'g');

yv=mean(RcLconpsth(:,winidx));
yste=steRLcon(winidx);
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[1 0.8 1]);
hold on
h3=plot(xv,yv,'m');

yv=mean(RcLtunpsth(:,winidx));
yste=steRLtun(winidx);
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[0.8 0.8 1]);
hold on
h4=plot(xv,yv,'b');
xlim([-199 300])
% ylim([0 180])
h5=legend([h1,h2,h3,h4],'Lori Lpref Con','Lpref tun','Rpref con','Rpref tun');
set(h5,'box','off')

subplot(1,2,2)
% xv=(win(1)+bin:bin:win(2))-0.1;
% xv=xv/bin;
xv=-199:300;
yv=mean(LcRconpsth(:,winidx));
yste=steLRcon(winidx);
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[1 0.8 0.8]);
hold on
h1=plot(xv,yv,'r');

yv=mean(LcRtunpsth(:,winidx));
yste=steLRtun(winidx);
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[0.8 1 0.8]);
hold on
h2=plot(xv,yv,'g');

yv=mean(RcRconpsth(:,winidx));
yste=steRRcon(winidx);
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[1 0.8 1]);
hold on
h3=plot(xv,yv,'m');

yv=mean(RcRtunpsth(:,winidx));
yste=steRRtun(winidx);
patch([xv fliplr(xv)],[yv-yste fliplr(yv+yste)],[0.8 0.8 1]);
hold on
h4=plot(xv,yv,'b');
xlim([-199 300])
% ylim([0 180])
h5=legend([h1,h2,h3,h4],'Rori Lpref Con','Lpref tun','Rpref con','Rpref tun');
set(h5,'box','off');


% d-prime
function nmpsth=normlizepsth(psth,blwin,blpsth)
BL1=repmat(mean(psth(:,blwin),2),1,size(psth,2));
BL2=repmat(mean(blpsth(:,blwin),2),1,size(blpsth,2));
% maxbin=repmat(max(blpsth-BL2,[],2),1,size(psth,2));
[maxvalue,maxidx]=max(blpsth-BL2,[],2);
validmaxidx=maxidx>230 & maxidx<270;
maxbin=repmat(maxvalue,1,size(psth,2));
tempnmpsth=(psth-BL1)./maxbin;
tempnmpsth=tempnmpsth(validmaxidx,:);
maxidx=maxidx(validmaxidx);
for n=1:size(tempnmpsth,1)
    nmpsth(n,:)=tempnmpsth(n,maxidx(n)-199:maxidx(n)+300);
end

