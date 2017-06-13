clear
str='I:\170606';
LR='L';
nPuffcnd=3;
if strcmpi(LR,'l')
    CNDlist=[4 5 6 1 2 3 1 2 3 1 2 3];
elseif strcmpi(LR,'r')
    CNDlist=[1 2 3 4 5 6 4 5 6 4 5 6];
else
    error('LR input is wrong, L or R');
end
namelist=dir(str);
ORICONlist=[];
DftORIlist=[];
for i=3:numel(namelist)
    if namelist(i).isdir && strcmp(namelist(i).name(1:6),'ORICON')
        ORICONlist=[ORICONlist; {namelist(i).name}];
%     elseif namelist(i).isdir && strcmp(namelist(i).name(1:6),'DftORI')
%         DftORIlist=[DftORIlist; {namelist(i).name}];
    end
end
if numel(ORICONlist)>12
    ORICONlist=ORICONlist(1:12);
end
NORICON=numel(ORICONlist);

%%
TimeCS=[];
TimeNEU=[];
CloseTCS=[];
CloseTNEU=[];
NCS=[];
NNEU=[];
figure
for n=1:numel(ORICONlist)
    subplot(3,4,n)
    [TraceTimeCS,TraceTimeNEU,PCloseTime,NCloseTime,PNBlink,NNBlink]=BlinkImageMKII([str '\' ORICONlist{n}],0:2,3:11,1);
    TimeCS=[TimeCS; sum(TraceTimeCS(~isnan(TraceTimeCS)))];
    TimeNEU=[TimeNEU; sum(TraceTimeNEU(~isnan(TraceTimeNEU)))];
    CloseTCS=[CloseTCS; sum(PCloseTime(~isnan(TraceTimeCS)))];
    CloseTNEU=[CloseTNEU; sum(NCloseTime(~isnan(TraceTimeNEU)))];
    NCS=[NCS; sum(PNBlink(~isnan(TraceTimeCS)))];
    NNEU=[NNEU; sum(NNBlink(~isnan(TraceTimeNEU)))];
end
sprintf(['CS   ' num2str(sum(CloseTCS)) '/' num2str(sum(TimeCS)) '   ' num2str(sum(CloseTCS)/sum(TimeCS))])
sprintf(['NEU   ' num2str(sum(CloseTNEU)) '/' num2str(sum(TimeNEU)) '   ' num2str(sum(CloseTNEU)/sum(TimeNEU))])
sprintf(['BLINKN   ' num2str(sum(NCS)) '/' num2str(sum(NNEU)/3.0)])

%%
% load([str '\5msbin_MergeconTrain.mat']);
% load([str '\15deg_5msbin_MergetunTrain.mat']);
% load([str '\ORIvalidid']);
% v1blockconpsth=cell(NORICON,1);
% V1oripsth=cell(6,1);
% Nmulti=nan;  % NEU puff 比例
% for a=1:NORICON
%     v1blockconpsth{a}=cell(2*nPuffcnd,1);
% end
% for j=1:12
%     for k=1:8
%         if ORIvalidid(j,k)
%             for m=1:2*nPuffcnd
%                 V1oripsth{m}=[V1oripsth{m}; MergetunTrain{j,k}{12+m*2}.Train];
%             end
%             if isnan(Nmulti)
%                 Nmulti=round((size(MergeconTrain{j,k}{4}.Train,1))/size(MergeconTrain{j,k}{1}.Train,1)-1);
%             end
%             for n=1:NORICON
%                 for m=1:2*nPuffcnd
% %                     ntrial=10+20*(m>3);
%                     ntrial=10+10*Nmulti*(m>3);
%                     v1blockconpsth{n}{CNDlist(m)}=[v1blockconpsth{n}{CNDlist(m)}; MergeconTrain{j,k}{m}.Train((n-1)*ntrial+1:n*ntrial,:)];
%                     % MergeconTrain为完全原始顺序， v1unitblockconpsth从原始顺序变成朝向顺序
%                     % rawunitCONPSTH并没有选细胞，所有记录细胞都保存，后续处理再选细胞
%                 end
%             end
%         end
%     end
% end
% for nn=1:NORICON
%     figure
%     for mm=1:6
%         subplot(2,3,mm)
%         plot(1:size(V1oripsth{mm},2),mean(V1oripsth{mm})/0.005,'g')
%         hold on
%         plot(1:size(V1oripsth{mm},2),mean(v1blockconpsth{nn}{mm})/0.005,'r')
%     end
% end
%%---------------------------------------------------------------------------------------------------------
% [S1merge,S2merge]=V1spec('con');
% [S1merge,S2merge]=V1spec('ori');
% % 
% for i=3:14
%         saveas(i,['C:\Users\Li Zhh\Documents\170604\location transfer\Block' num2str(i-2) 'psth.fig'])
% end