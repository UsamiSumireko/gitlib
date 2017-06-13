function [LFPSldCoh]=V1_Amy_lfp_sldwin_coherency(destination,win,shufflelabel,savelabel,LR,CONORI)
% [LFPCoh]=V1_Amy_lfpconherency(destination,win,shufflelabel,savelabel)
% shufflelabel =0;
% savelabel=[1 1 1];%  save Coh Amy V1
% win=[301 1300];
% destination='i:\160720';
% CONORI='CON';
% LR= l r  L==0,R==1
%%
if numel(savelabel)~=1
    error('savelabel should be like [1 0 0]');
end 
tic
shuffleL=shufflelabel;
if strcmpi(CONORI,'con')
    load([destination '\V1CONblockunitLFP.mat'])
    V1blockunitLFP=V1CONblockunitLFP;
    load([destination '\Amy_unit_block_conLFP.mat'])
    amyunitblockLFP=amyunitblockconLFP;
    load([destination '\Blksortedcell.mat'])
    clearvars V1CONblockunitLFP amyunitblockconLFP
elseif strcmpi(CONORI,'ori')
    load([destination '\V1ORIblockunitLFP.mat'])
    V1blockunitLFP=V1ORIblockunitLFP;
    load([destination '\Amy_unit_block_oriLFP.mat'])
    amyunitblockLFP=amyunitblockoriLFP;
    load([destination '\ORIBlksortedcell.mat'])
    Blksortedcell=ORIBlksortedcell;
    clearvars V1ORIblockunitLFP amyunitblockoriLFP ORIBlksortedcell
end
load([destination '\ORIvalidid.mat'])
load([destination '\XYvalidid.mat'])
params.Fs = 1000; % sampling frequency
params.tapers = [2 3]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 0;
params.err =[2 0.05];
if strcmpi(CONORI,'con')
    params.iscon=1;
elseif strcmpi(CONORI,'ori')
    params.iscon=0;
end
if strcmpi(LR,'l')
    params.lr=0;
elseif strcmpi(LR,'r')
    params.lr=1;
end
%% reform AMY lfp
if strcmpi(CONORI,'con')
    str='\AmylfpMKII';
elseif strcmpi(CONORI,'ori')
    str='\AmyorilfpMKII';
end
if ~exist([destination str num2str(win(1)-300) '-' num2str(win(2)-300) 'ms.mat'],'file')
    celldata=cellfun(@(x) merge_block(x,win,params),amyunitblockLFP,'UniformOutput',false);
    unfoldeddata=cell2mat(celldata(:));
    zzz=mat2cell(unfoldeddata,size(unfoldeddata,1),repmat(size(unfoldeddata,2)/numel(Blksortedcell),1,numel(Blksortedcell)));
    ccc=cellfun(@(x) mat2cell(x,repmat(win(2)-win(1)+1,24,1),size(x,2)),zzz,'UniformOutput',false);
    AmylfpMKII=cellfun(@(x,y) x(y(:)),ccc,Blksortedcell,'UniformOutput',false);
else
    load([destination str num2str(win(1)-300) '-' num2str(win(2)-300) 'ms.mat'])
end
%% V1 lfp
if strcmpi(CONORI,'con')
    str='\V1lfp';
elseif strcmpi(CONORI,'ori')
    str='\V1orilfp';
end
if ~exist([destination str num2str(win(1)-300) '-' num2str(win(2)-300) 'ms.mat'],'file')
    ccc=cellfun(@(x) merge_block(x,win,params),V1blockunitLFP,'UniformOutput',false);
    Nv1validcell=sum(sum(0~=cellfun(@length,V1blockunitLFP{1})));
    V1lfp=cellfun(@(x) mat2cell(x,size(x,1),repmat(size(x,2)/Nv1validcell,1,Nv1validcell)),ccc,'UniformOutput',false);
else
    load([destination str num2str(win(1)-300) '-' num2str(win(2)-300) 'ms.mat'])
end
%% compute
LFPSldCoh=cellfun(@(x,y) block_coh(x,y,params,shuffleL),AmylfpMKII',V1lfp,'UniformOutput',false);
if strcmpi(CONORI,'con')
    str='';
elseif strcmpi(CONORI,'ori')
    str='ori';
end
if savelabel(1)
    if shuffleL
    save([destination '\' str 'LFPSldCoh_shuffle' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms'],'LFPSldCoh','-v7.3');    
    else
    save([destination '\' str 'LFPSldCoh' num2str(win(1)-300) '-' num2str(win(2)-300) 'ms'],'LFPSldCoh','-v7.3');
    end
end
toc
%%
function b=merge_block(a,win,params)
id=cellfun('length',a);
a(id==0)=[];
a=a(:);
bb=cellfun(@(x) cell2mat(x'),a,'UniformOutput',false);
bblsc=cellfun(@(x) rmlinesc(double(x),params,'n',50),bb,'UniformOutput',false);
b=cell2mat(bblsc');
b=b(win(1):win(2),:);

function C=block_coh(amy,v1,params,shuffleL)
C=cellfun(@(x) amy_coh(x,v1,params,shuffleL),amy,'UniformOutput',false);

function Coh=amy_coh(amy,v1,params,shuffleL)
Coh=struct;
[c1,phi1,c2,phi2,cerr1,cerr2,t1,t2,f1,f2]=cellfun(@(x) amy_v1_coh(amy,x,params,shuffleL),v1,'UniformOutput',false);
Coh.C1=cell2mat(reshape(c1,1,1,[]));
Coh.Phi1=cell2mat(reshape(phi1,1,1,[]));
Coh.C2=cell2mat(reshape(c2,1,1,[]));
Coh.Phi2=cell2mat(reshape(phi2,1,1,[]));
Coh.Cerr1=cell2mat(reshape(cerr1,1,1,[]));
Coh.Cerr2=cell2mat(reshape(cerr2,1,1,[]));
Coh.f1=f1;
Coh.f2=f2;
Coh.t1=t1;
Coh.t2=t2;

function [c1,phi1,c2,phi2,cerr1,cerr2,t1,t2,f1,f2]=amy_v1_coh(amy,v1,params,shuffleL)
if params.iscon
    conidx=(1:30)+~params.lr*90;
    neuidx=randperm(90,30)+params.lr*30;
else
    conidx=randperm(60,30)+~params.lr*60;
    neuidx=randperm(60,30)+params.lr*60;
end
if shuffleL
    Nshuffle=20;
    randtim=1:Nshuffle;
    randidx=arrayfun(@(x) randperm(30),randtim,'UniformOutput',false);
    [c1c,phi1c,~,~,~,f1c,~,~,cerr1c]=cellfun(@(x) cohgramc(amy(:,conidx),v1(:,x),params),randidx,'UniformOutput',false);
    [c2c,phi2c,~,~,~,f2c,~,~,cerr2c]=cellfun(@(x) cohgramc(amy(:,neuidx),v1(:,neuidx(x)),params),randidx,'UniformOutput',false);
    c1=mean(cell2mat(reshape(c1c,1,1,1,[])),2);
    phi1=mean(cell2mat(reshape(phi1c,1,1,1,[])),2);
    cerr1=mean(cell2mat(reshape(cerr1c,1,1,1,[])),2);
    f1=f1c{1};
    c2=mean(cell2mat(reshape(c2c,1,1,1,[])),2);
    phi2=mean(cell2mat(reshape(phi2c,1,1,1,[])),2);
    cerr2=mean(cell2mat(reshape(cerr2c,1,1,1,[])),2);
    f2=f2c{1};
else
    [c1,phi1,~,~,~,t1,f1,~,~,cerr1]=cohgramc(amy(:,conidx),v1(:,conidx),[0.2 0.01],params);
    [c2,phi2,~,~,~,t2,f2,~,~,cerr2]=cohgramc(amy(:,neuidx),v1(:,neuidx),[0.2 0.01],params);
end
% 
% C1poll=[];
% C2poll=[];
% for i=1:numel(LFPSldCoh)
%     for j=1:numel(LFPSldCoh{i})
%         C1poll=cat(3,C1poll,LFPSldCoh{i}{j}.C1);
%         C2poll=cat(3,C2poll,LFPSldCoh{i}{j}.C2);
%     end
% end
% figure
% subplot(1,2,1)
% plot_matrix(squeeze(mean(C1poll,3)),0.1:0.01:0.3,1:24,'n');
% subplot(1,2,2)
% plot_matrix(squeeze(mean(C2poll,3)),0.1:0.01:0.3,1:24,'n');
% 
% figure
% plot_matrix(squeeze(mean(C1poll,3))./squeeze(mean(C2poll,3)),0.1:0.01:0.3,1:24,'n');