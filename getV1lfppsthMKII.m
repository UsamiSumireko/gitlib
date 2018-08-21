function getV1lfppsthMKII(destination,LR,select,CONORIlabel)
%% for digrating conditioning use /180114
%%
% destination='H:\180112';
% LR='L';
% win=[-0.3 1.1];
% CONORIlabel=[1 1]  [CON ORI]  
tic
cuelabel=false;
conwin=[-0.1 1.2];
oriwin=[-0.1 0.5];
%%
slashidx=strfind(destination,'\');
slashidx=slashidx(end);
cmppath = 'D:\SN 4566-001277.cmp';
cmpinfo=LoadCmp(cmppath,1,0);
elec=cmpinfo{1,1}.RealElec(1:12,:);
load([destination '\XYvalidid']);
load([destination '\ORIvalidid']);
selecti=lower(select);
fs=1000;
nsample=30;
switch selecti
    case 'xy'
        valididx=XYvalidid;
    case 'ori'
        valididx=ORIvalidid;
    case 'both'
        valididx=ORIvalidid & XYvalidid;
    case 'all'
        valididx=ones(12,8);
    otherwise
        msgbox('wrong select idx')
end
valididx=valididx(1:12,:);
validelec=elec(valididx);
lr=lower(LR);
if strcmp(lr,'l')
    CNDlist=[2 1 1 1 4 3 3 3];
elseif strcmp(lr,'r')
    CNDlist=[1 2 2 2 3 4 4 4];
else
    error('LR input is wrong, L or R');
end
% 1 RF45 contra45
% 2 RF45 contra135
% 3 FR135 contra45
% 4 RF135 contra135
%%
namelist=dir(destination);
NCON=0;
CONlist=[];
NORI=0;
ORIlist=[];
for i=3:numel(namelist)
    if namelist(i).isdir && strcmpi(namelist(i).name(1:9),'DiGratCon')
        if ~xor(strcmpi(namelist(i).name(11),'c'),cuelabel)
            NCON=NCON+1;
            CONlist=[CONlist; {namelist(i).name}];
        end
    elseif namelist(i).isdir && strcmpi(namelist(i).name(1:9),'DiGratTun')
        NORI=NORI+1;
        ORIlist=[ORIlist; {namelist(i).name}];
    end
end

%%
if CONORIlabel(1)
    if ~exist(['D:\data extracted\' destination(slashidx+1:end) '\V1conblockunitlfp.mat'],'file')
        V1conblockunitlfp=struct('lfp',[],'timelabel',[],'LR',LR);
        V1conblockunitlfp.lfp=nan(floor((conwin(2)-conwin(1))*fs*30/nsample),10,8,sum(sum(valididx)),NCON);
        for n=1:NCON
            n,
            for k=1:sum(sum(valididx))
                sprintf(repmat('.',1,k))
                load([destination '\' CONlist{n} '\ExpMonitor.mat']);
                validtrialid=1==ExpMonitor.RespCode;
%                 vlaidstimid=ExpMonitor.StimCnd(validtrialid);
                for mstim=0:max(ExpMonitor.StimCnd)
                    mstimidx=mstim==ExpMonitor.StimCnd;
                    secstimlabel(:,mstim+1,k,n)=ExpMonitor.ISecPosStim(validtrialid & mstimidx);
                end
                lfpinfo=TruncateLFP([destination '\' CONlist{n} '.ns5'],validelec(k),conwin,'MKII');
                [splitlfp]=SplitInfo(lfpinfo, [1 0 1 0]);
                if isempty(V1conblockunitlfp.timelabel)
                    temp=downsample(lfpinfo.TimeLabel,nsample);
                    V1conblockunitlfp.timelabel=temp(2:end);
                end
                for m=1:numel(splitlfp)
                    dslfp=downsample(splitlfp{m}{1}.LFP',nsample);
                    V1conblockunitlfp.lfp(:,:,m,k,n)=dslfp(2:end,:);
                end
            end
        end
        V1conblockunitlfp.secstimlabel=secstimlabel;
        if ~exist(['D:\data extracted\' destination(slashidx+1:end)],'dir')
            mkdir(['D:\data extracted\' destination(slashidx+1:end)])
        end
        save(['D:\data extracted\' destination(slashidx+1:end) '\V1conblockunitlfp.mat'],'V1conblockunitlfp','-v7.3');
%     else
%         load(['D:\data extracted\' destination(end-5:end) '\V1conblockunitlfp_' num2str(nsample) '.mat']);
    end
end
sprintf('CON done')
%%
if CONORIlabel(2)
    if ~exist(['D:\data extracted\' destination(slashidx+1:end) '\V1oriblockunitlfp.mat'],'file')
        V1oriblockunitlfp=struct('lfp',[],'timelabel',[],'LR',LR);
        V1oriblockunitlfp.lfp=nan(floor((oriwin(2)-oriwin(1))*fs*30/nsample),20,4,sum(sum(valididx)),NORI);
        for n=1:NORI
            n,
            for k=1:sum(sum(valididx))
                sprintf(repmat('.',1,k))
                lfpinfo=TruncateLFP([destination '\' ORIlist{n} '.ns5'],validelec(k),oriwin,'MKII');
                [splitlfp]=SplitInfo(lfpinfo, [1 0 1 0]);
                if isempty(V1oriblockunitlfp.timelabel)
                    temp=downsample(lfpinfo.TimeLabel,nsample);
                    V1oriblockunitlfp.timelabel=temp(2:end);
                end
                for m=1:4
                    dslfp=downsample(splitlfp{m}{1}.LFP',nsample);
                    V1oriblockunitlfp.lfp(:,:,m,k,n)=dslfp(2:end,:);
                end
            end
        end
        if ~exist(['D:\data extracted\' destination(slashidx+1:end)],'dir')
            mkdir(['D:\data extracted\' destination(slashidx+1:end)])
        end
        save(['D:\data extracted\' destination(slashidx+1:end) '\V1oriblockunitlfp.mat'],'V1oriblockunitlfp','-v7.3');
%     else
%         load(['D:\data extracted\' destination(end-5:end)  '\V1oriblockunitlfp_' num2str(nsample) '.mat']);
    end
end
sprintf('ORI done')
toc
%%
tic
V1conblockunitpsth=struct('psth',[],'LR',LR);
V1conblockunitpsth.psth=nan(floor((conwin(2)-conwin(1))*fs),10,8,sum(sum(valididx)),NCON);
V1oriblockunitpsth=struct('psth',[],'LR',LR);
V1oriblockunitpsth.psth=nan(floor((oriwin(2)-oriwin(1))*fs),20,4,sum(sum(valididx)),NORI);
if exist('secstimlabel','var')
    V1conblockunitpsth.secstimlabel=secstimlabel;
else
    for n=1:NCON
        for k=1:sum(sum(valididx))
            load([destination '\' CONlist{n} '\ExpMonitor.mat']);
            validtrialid=1==ExpMonitor.RespCode;
            %                 vlaidstimid=ExpMonitor.StimCnd(validtrialid);
            for mstim=0:max(ExpMonitor.StimCnd)
                mstimidx=mstim==ExpMonitor.StimCnd;
                secstimlabel(:,mstim+1,k,n)=ExpMonitor.ISecPosStim(validtrialid & mstimidx);
            end
        end
    end
    V1conblockunitpsth.secstimlabel=secstimlabel;
end
if ~exist(['D:\data extracted\' destination(slashidx+1:end) '\V1conblockunitpsth.mat'],'file')
    for n=1:NCON
        for k=1:sum(sum(valididx))
            spikeinfo=BinSpike([destination '\' CONlist{n} '.nev'],validelec(k),'t',conwin,0.001);
            splitspk=SplitInfo(spikeinfo,[1 0 1 0]);
            for m=1:numel(splitspk)
                V1conblockunitpsth.psth(:,:,m,k,n)=splitspk{m}{1}.Train';
            end
        end
    end
    save(['D:\data extracted\' destination(slashidx+1:end) '\V1conblockunitpsth.mat'],'V1conblockunitpsth','-v7.3')
end
if ~exist(['D:\data extracted\' destination(slashidx+1:end) '\V1oriblockunitpsth.mat'],'file')
    for n=1:NORI
        for k=1:sum(sum(valididx))
            spikeinfo=BinSpike([destination '\' ORIlist{n} '.nev'],validelec(k),'t',oriwin,0.001);
            splitspk=SplitInfo(spikeinfo,[1 0 1 0]);
            for m=1:4
                V1oriblockunitpsth.psth(:,:,m,k,n)=splitspk{m}{1}.Train';
            end
        end
    end
    save(['D:\data extracted\' destination(slashidx+1:end) '\V1oriblockunitpsth.mat'],'V1oriblockunitpsth','-v7.3')
end
%%
toc
