function [S1merge,S2merge]=V1spec(CONORI)
win=[401 800];
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
%     'G:\MS\160624';
%     'G:\MS\160625';
%     'G:\MS\160627';
%     'G:\MS\160628';
'I:\170603'
% 'I:\170519'
% 'I:\170522'
% 'I:\170523'
% 'I:\170502'
% 'I:\170503'
    };
% LRlist='rrrrrrrrrllll';
% LRlist='rrrrrrrrr';
% LRlist='llll';
LRlist='llllllll';
params.Fs = 1000; % sampling frequency
params.tapers = [2 3]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 0;
params.err =[2 0.05]; 
% params.lr   L==0 R==1;
params.iscon=strcmpi(CONORI,'con'); % ori=0 con=1

S1merge=[];
S2merge=[];
for i=1:numel(folderlist)
    tic
    params.lr=strcmp(LRlist(i),'r');  % L==0 R==1;
    if params.iscon
        load([folderlist{i} '\V1CONblockunitLFP.mat']);
        V1blockunitLFP=V1CONblockunitLFP;
    else
        load([folderlist{i} '\V1ORIblockunitLFP.mat']);
        V1blockunitLFP=V1ORIblockunitLFP;
    end
    spec=cellfun(@(x) merge_block(x,win,params),V1blockunitLFP,'UniformOutput',false);
    s1=cellfun(@(x) cell2mat(x.S1),spec,'UniformOutput',false);
    S1=cell2mat(s1');
    s2=cellfun(@(x) cell2mat(x.S2),spec,'UniformOutput',false);
    S2=cell2mat(s2');
    f=spec{1}.f{1};
    S1merge=[S1merge S1];
    S2merge=[S2merge S2];
    V1_spec=struct('S1',S1,'S2',S2,'f',f);
    save([folderlist{i} '\V1_spec_' CONORI '_' num2str(win(1)-300) '_' num2str(win(2)-300) 'ms.mat'],'V1_spec','-v7.3');
    toc
end
figure
plot(f,mean(S1merge,2),'r')
hold on
plot(f,mean(S2merge,2),'g')
sprintf('done')
function spec=merge_block(a,win,params)
spec=struct;
l=cellfun(@length,a);
a(l==0)=[];
aa=cellfun(@(x) cell2mat(x'),a,'UniformOutput',false);
[spec.S1,spec.S2,spec.f]=cellfun(@(x) cellspec(x,win,params),aa,'UniformOutput',false);
function [S1,S2,f]=cellspec(a,win,params)
if params.iscon
csidx=(1:30)+~params.lr*90;
neuidx=randperm(90,30)+params.lr*30;
else
csidx=randperm(60,30)+~params.lr*60;
neuidx=randperm(60,30)+params.lr*60;
end
rmlfp=rmlinesc(double(a),params,'n',50);
[S1,f]=mtspectrumc(rmlfp(win(1):win(2),csidx),params);
[S2,~]=mtspectrumc(rmlfp(win(1):win(2),neuidx),params);
