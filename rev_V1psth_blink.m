clear
revdaylist={
    'i:\160805'
    'i:\160908'
    'i:\170505'
    };
Nblock=8;
% LRlist='L';
LRlist='LRL';
preframelist=[0.2 0.2 0.1];
blockspkdiff=[];
blockblinkratio=[];
for iday=1:numel(revdaylist)
    if strcmpi(LRlist(iday),'r')
        cndlist=[1 2 3 4 5 6];
    elseif strcmpi(LRlist(iday),'l')
        cndlist=[4 5 6 1 2 3];
    end
    load([revdaylist{iday} '\XYvalidid']);
    load([revdaylist{iday} '\ORIvalidid']);
    load([revdaylist{iday} '_prerev\1msbin_rawunitCONPSTH.mat']);
    load([revdaylist{iday} '\1msbin_rawunitORIPSTH.mat']);
    precspeak=[];
    preneupeak=[];
    oricspeak=[];
    orineupeak=[];
    blinkTdiff=[];
    for j=1:12
        for k=1:8
            if ORIvalidid(j,k)
                for n=1:3
                    temp1=mean(mean(rawunitCONPSTH{j,k}{cndlist(n)}(:,preframelist(iday)/0.001+(31:70))))/0.001;
                    temp2=mean(mean(rawunitCONPSTH{j,k}{cndlist(n+3)}(:,preframelist(iday)/0.001+(31:70))))/0.001;
%                     temp3=mean(mean(rawunitORIPSTH{j,k}{cndlist(n)}(:,preframelist(iday)/0.001+(31:70))))/0.001;
%                     temp4=mean(mean(rawunitORIPSTH{j,k}{cndlist(n+3)}(:,preframelist(iday)/0.001+(31:70))))/0.001;
                    precspeak=[precspeak; temp1];
                    preneupeak=[preneupeak; temp2];
                end
                temp3=mean(mean(rawunitORIPSTH{j,k}{cndlist(2)}(:,preframelist(iday)/0.001+(31:70))))/0.001;
                temp4=mean(mean(rawunitORIPSTH{j,k}{cndlist(5)}(:,preframelist(iday)/0.001+(31:70))))/0.001;
                oricspeak=[oricspeak; temp3];
                orineupeak=[orineupeak; temp4];
            end
        end
    end
    blockcsdiff(iday,1)=mean(precspeak)-mean(oricspeak);
    blockneudiff(iday,1)=mean(preneupeak)-mean(orineupeak);
    blockspkdiff(iday,1)=mean(precspeak)-mean(oricspeak)-mean(preneupeak)+mean(orineupeak);
    
    load([revdaylist{iday} '\V1_unit_block_conpsth.mat']);
    L=cellfun(@length,v1unitblockconpsth);
    v1unitblockconpsth(L==0)=[];
    for nn=1:Nblock
        cspeak=[];
        neupeak=[];
        for l=1:numel(v1unitblockconpsth)
            for m=1:3
                temp1=mean(mean(v1unitblockconpsth{l}{nn}{cndlist(m)}(:,preframelist(iday)/0.001+(31:70))))/0.001;
                temp2=mean(mean(v1unitblockconpsth{l}{nn}{cndlist(m+3)}(:,preframelist(iday)/0.001+(31:70))))/0.001;
                cspeak=[cspeak; temp1];
                neupeak=[neupeak; temp2];
            end
        end
        blockcsdiff(iday,nn+1)=mean(cspeak)-mean(oricspeak);
        blockneudiff(iday,nn+1)=mean(neupeak)-mean(orineupeak);
        blockspkdiff(iday,nn+1)=mean(cspeak)-mean(oricspeak)-mean(neupeak)+mean(orineupeak);
    end
    
    NORICON=0;
    ORICONlist=[];
    namelist=dir([revdaylist{iday} '_prerev']);
    for i=3:numel(namelist)
        if namelist(i).isdir && strcmp(namelist(i).name(1:6),'ORICON')
            NORICON=NORICON+1;
            ORICONlist=[ORICONlist; {namelist(i).name}];
        end
    end
    for n=1:numel(ORICONlist)
        [TraceTimeCS,TraceTimeNEU,PCloseTime,NCloseTime,PNBlink,NNBlink]=BlinkImageMKII([revdaylist{iday} '_prerev\' ORICONlist{n}],0:2,3:11,0);
        preNBlinkCS(n)=sum(sum(PNBlink~=0));
        preNBlinkNEU(n)=sum(sum(NNBlink~=0));
        preTrialNEU(n)=size(NNBlink,1)*size(NNBlink,2);
        preTrialCS(n)=size(PNBlink,1)*size(PNBlink,2);
        nonnanidx= isnan(PCloseTime) | isnan(TraceTimeCS);
        preTBlinkCS(n)=sum(PCloseTime(~nonnanidx))/sum(TraceTimeCS(~nonnanidx));
        nonnanidx= isnan(NCloseTime) | isnan(TraceTimeNEU);
        preTBlinkNEU(n)=sum(NCloseTime(~nonnanidx))/sum(TraceTimeNEU(~nonnanidx));
    end
    blockblinkratio(iday,1)=sum(preNBlinkCS)/sum(preTrialCS)/sum(preNBlinkNEU)*sum(preTrialNEU);
    blockCSblinkT(iday,1)=mean(preTBlinkCS);
    blockNEUblinkT(iday,1)=mean(preTBlinkNEU);
    
    NORICON=0;
    ORICONlist=[];
    namelist=dir(revdaylist{iday});
    for i=3:numel(namelist)
        if namelist(i).isdir && strcmp(namelist(i).name(1:6),'ORICON')
            NORICON=NORICON+1;
            ORICONlist=[ORICONlist; {namelist(i).name}];
        end
    end
    for n=1:Nblock
    [TraceTimeCS,TraceTimeNEU,PCloseTime,NCloseTime,PNBlink,NNBlink]=BlinkImageMKII([revdaylist{iday} '\' ORICONlist{n}],0:2,3:11,0);
    NBlinkCS(n)=sum(sum(PNBlink~=0));
    NBlinkNEU(n)=sum(sum(NNBlink~=0));
    TrialNEU(n)=size(NNBlink,1)*size(NNBlink,2);
    TrialCS(n)=size(PNBlink,1)*size(PNBlink,2); 
    nonnanidx= isnan(PCloseTime) | isnan(TraceTimeCS);
    TBlinkCS(n)=sum(PCloseTime(~nonnanidx))/sum(TraceTimeCS(~nonnanidx));
    nonnanidx= isnan(NCloseTime) | isnan(TraceTimeNEU);
    TBlinkNEU(n)=sum(NCloseTime(~nonnanidx))/sum(TraceTimeNEU(~nonnanidx));
    blockblinkratio(iday,n+1)=NBlinkCS(n)/TrialCS(n)/NBlinkNEU(n)*TrialNEU(n);
    blockCSblinkT(iday,n+1)=TBlinkCS(n);
    blockNEUblinkT(iday,n+1)=TBlinkNEU(n);
    end
end
figure
bar(1:Nblock+1,[mean(blockcsdiff,1)' mean(blockneudiff,1)'])
figure
bar(1:Nblock+1,[mean(blockCSblinkT,1)' mean(blockNEUblinkT,1)'])
figure
yyaxis left
plot(1:Nblock+1,mean(blockblinkratio,1),'b')
line([1 9],[1 1],'color',[0.8 0.8 1],'LineStyle','--');
ylabel('blink rate proportion')
yyaxis right
plot(1:Nblock+1,mean(blockspkdiff,1),'r')
line([1 9],[0 0],'color',[1 0.8 0.8],'LineStyle','--');
ylabel('Firing Rate diff')

% [AX,H1,H2]=plotyy(1:Nblock+1,mean(blockblinkratio,1),1:Nblock+1,mean(blockspkdiff,1));

