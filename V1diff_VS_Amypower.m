function [allv1diff,allamypower]=V1diff_VS_Amypower()
%V1DIFF_VS_AMYPOWER 此处显示有关此函数的摘要
%   此处显示详细说明
folderlist={
    'I:\160712';
    'I:\160713';
    'I:\160714';
    'I:\160715';
    'I:\160716';
    'I:\160719';
    'I:\160720';
    'I:\160721';
    'I:\160722';
    'G:\MS\160624';
    'G:\MS\160625';
    'G:\MS\160627';
    'G:\MS\160628';
    };
LRlist='rrrrrrrrrllll';
% LRlist='llll';
allv1diff=[];
allamypower=[];

params.Fs = 1000; % sampling frequency
params.tapers = [2 3]; % taper parameters
params.trialave = 0;
params.fpass = [4 100];
params.pad = 0;
params.err =[2 0.05];
params.fqpass=[60 90];
params.V1win=[431 470];
params.Amywin=[601 1000];
rlist=[];
plist=[];
for iday=1:numel(folderlist)
    iday,
    load([folderlist{iday} '\AmylfpMKII1-1000ms.mat']);
    load([folderlist{iday} '\V1_unit_block_diffpsth.mat']);
    l=cellfun(@length,AmylfpMKII);
    AmylfpMKII(l==0)=[];
    BlockTrialPower=cellfun(@(x) blockspec(x,params),AmylfpMKII,'UniformOutput',false);
    MergeBlockTrialPower=cellfun(@cell2mat,BlockTrialPower,'UniformOutput',false);
    TrialPower=cellfun(@mean,MergeBlockTrialPower,'UniformOutput',false);
%     TrialPower=cellfun(@(x) x(2,:),MergeBlockTrialPower,'UniformOutput',false);
    MergeTrialPower=transpose(cell2mat(TrialPower));
    allamypower=[allamypower; MergeTrialPower];
    L=cellfun(@length,v1unitblockdiffpsth);
    v1unitblockdiffpsth(0==L)=[];
    tempBlockTrialDiff=cellfun(@(x) blockdiff(x,params),v1unitblockdiffpsth,'UniformOutput',false);
    BlockTrialDiff=cellfun(@(x) x(l~=0),tempBlockTrialDiff,'UniformOutput',false);
    TrialDiff=cellfun(@(x) cell2mat(x'),BlockTrialDiff,'UniformOutput',false);
    MergeTrialDiff=mean(cell2mat(TrialDiff),2);
    allv1diff=[allv1diff; MergeTrialDiff];
    Blkdiff=mat2cell(MergeTrialDiff,repmat(numel(MergeTrialDiff)/numel(TrialPower),numel(TrialPower),1));
    [bR,bP]=cellfun(@corrcoef,Blkdiff',TrialPower,'UniformOutput',false);
    br=cellfun(@(x) x(1,2),bR);
    bp=cellfun(@(x) x(1,2),bP);
    rlist=[rlist br];
    plist=[plist bp];
end        
figure
plot(allv1diff.*1000,allamypower,'b.');
[R,P]=corrcoef(allv1diff,allamypower);
sprintf(['done,r=' num2str(R(1,2)) ' p=' num2str(P(1,2))])
figure
hist(rlist)
end
function blocktrialpower=blockspec(a,params)
[cellspec,F]=cellfun(@(x) mtspectrumc(x(params.Amywin(1):params.Amywin(2),:),params),a,'UniformOutput',false);
f=F{1};
fidx=f>params.fqpass(1) & f<params.fqpass(2);
% blocktrialpower=cellfun(@(x) mean(x(fidx,:)),cellspec,'UniformOutput',false);    
blocktrialpower=cellfun(@(x) zscore(mean(x(fidx,:))),cellspec,'UniformOutput',false);    
end
function blocktrialdiff=blockdiff(a,params)
mergespk=cellfun(@(x) cell2mat(x'),a,'UniformOutput',false);
blocktrialdiff=cellfun(@(x) mean(x(:,params.V1win(1):params.V1win(2)),2),mergespk,'UniformOutput',false);
end