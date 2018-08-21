function getV1lfppsth(destination,LR,select,CONORIlabel,filtlabel)
tic
% destination='I:\160720';
% LR='R';
% win=[-0.3 1.1];
% CONORIlabel=[1 1]  [CON ORI]  
conwin=[-0.3 1.1];
oriwin=[-0.3 0.5];
load('D:\data extracted\LFPMedium.mat');
%%
slashidx=strfind(destination,'\');
slashidx=slashidx(end);
% cmppath = 'D:\SN 4566-001303+Probe MK2.cmp';
cmppath =  'D:\SN 4566-001277.cmp';
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
    CNDlist=[4 5 6 1 2 3 1 2 3 1 2 3];
elseif strcmp(lr,'r')
    CNDlist=[1 2 3 4 5 6 4 5 6 4 5 6];
else
    error('LR input is wrong, L or R');
end
%%
namelist=dir(destination);
NORICON=0;
ORICONlist=[];
NDftORI=0;
DftORIlist=[];
for i=3:numel(namelist)
    if strcmpi(namelist(i).name(1:6),'ORICON') && namelist(i).isdir
        NORICON=NORICON+1;
        ORICONlist=[ORICONlist; {namelist(i).name}];
    elseif strcmpi(namelist(i).name(1:6),'DftORI') && namelist(i).isdir
        NDftORI=NDftORI+1;
        DftORIlist=[DftORIlist; {namelist(i).name}];
    end
end

%%
if CONORIlabel(1)
    if ~exist(['D:\data extracted\' destination(slashidx+1:end) '\V1conblockunitlfp.mat'],'file')
        V1conblockunitlfp=struct('lfp',[],'timelabel',[],'LR',LR);
        V1conblockunitlfp.lfp=nan(floor((conwin(2)-conwin(1))*fs*30/nsample),10,12,sum(sum(valididx)),NORICON);
        for n=1:NORICON
            n,
            %%
%             load([destination '\' ORICONlist{n} '\ExpMonitor.mat']);
%             validtrial=1==ExpMonitor.RespCode;
%             Tstart=ExpMonitor.StartT(validtrial);
%             if max(Tstart)>=conwin(1)
%                 trialnum=sum(0==ExpMonitor.StimCnd(validtrial));
%                 cndnum=numel(unique(ExpMonitor.StimCnd));
%                 V1conblockunitlfp.lfp(:,:,:,:,n)=nan(floor((conwin(2)-conwin(1))*fs*30/nsample),trialnum,cndnum,sum(sum(valididx)));
%                 continue % for n=1:NORICON
%             end
            %%
            for k=1:sum(sum(valididx))
                sprintf(repmat('.',1,k))
                lfpinfo=TruncateLFP([destination '\' ORICONlist{n} '.ns5'],validelec(k),conwin,'MKII');
                [splitlfp]=SplitInfo(lfpinfo, [1 0 1 0]);
                if isempty(V1conblockunitlfp.timelabel)
                    temp=downsample(lfpinfo.TimeLabel,nsample);
                    V1conblockunitlfp.timelabel=temp(2:end);
                end
                for m=1:numel(splitlfp)
                    lfptemp=double(splitlfp{m}{1}.LFP');
                    if filtlabel
                        lfptemp=filtfilt(SOS,G,lfptemp);
                    end
                    dslfp=downsample(lfptemp,nsample);
                    V1conblockunitlfp.lfp(:,:,m,k,n)=dslfp(2:end,:);
                end
            end
        end
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
        V1oriblockunitlfp.lfp=nan(floor((oriwin(2)-oriwin(1))*fs*30/nsample),20,6,sum(sum(valididx)),NDftORI);
        for n=1:NDftORI
            n,
            %%
%             load([destination '\' DftORIlist{n} '\ExpMonitor.mat']);
%             validtrial=1==ExpMonitor.RespCode;
%             Tstart=ExpMonitor.StartT(validtrial);
%             if max(Tstart)>=oriwin(1)
%                 trialnum=sum(0==ExpMonitor.StimCnd(validtrial));
%                 V1oriblockunitlfp.lfp(:,:,:,:,n)=nan(floor((conwin(2)-conwin(1))*fs*30/nsample),trialnum,6,sum(sum(valididx)));
%                 continue % for n=1:NDftORI
%             end
            %%
            for k=1:sum(sum(valididx))
                sprintf(repmat('.',1,k))
                lfpinfo=TruncateLFP([destination '\' DftORIlist{n} '.ns5'],validelec(k),oriwin,'MKII');
                [splitlfp]=SplitInfo(lfpinfo, [1 0 1 0]);
                if isempty(V1oriblockunitlfp.timelabel)
                    temp=downsample(lfpinfo.TimeLabel,nsample);
                    V1oriblockunitlfp.timelabel=temp(2:end);
                end
                for m=1:6
                    lfptemp=double([splitlfp{m*2}{1}.LFP' splitlfp{m*2+12}{1}.LFP']);
                    if filtlabel
                        lfptemp=filtfilt(SOS,G,lfptemp);
                    end
                    dslfp=downsample(lfptemp,nsample);
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
% V1conblockunitpsth.psth=nan(floor((conwin(2)-conwin(1))*fs),10,12,sum(sum(valididx)),NORICON);
V1oriblockunitpsth=struct('psth',[],'LR',LR);
% V1oriblockunitpsth.psth=nan(floor((oriwin(2)-oriwin(1))*fs),20,6,sum(sum(valididx)),NDftORI);
if CONORIlabel(1)
if ~exist(['D:\data extracted\' destination(slashidx+1:end) '\V1conblockunitpsth.mat'],'file')
    for n=1:NORICON
        %%
%         load([destination '\' ORICONlist{n} '\ExpMonitor.mat']);
%         validtrial=1==ExpMonitor.RespCode;
%         Tstart=ExpMonitor.StartT(validtrial);
%         if max(Tstart)>=-0.2
%             trialnum=sum(0==ExpMonitor.StimCnd(validtrial));
%             cndnum=numel(unique(ExpMonitor.StimCnd));
%             V1conblockunitpsth.psth(:,:,:,:,n)=nan(floor((conwin(2)-conwin(1))*fs*30/nsample),trialnum,cndnum,sum(sum(valididx)));
%             continue % for n=1:NORICON
%         end
        %%
        for k=1:sum(sum(valididx))
            spikeinfo=BinSpike([destination '\' ORICONlist{n} '.nev'],validelec(k),'t',conwin,0.001);
            splitspk=SplitInfo(spikeinfo,[1 0 1 0]);
            for m=1:numel(splitspk)
                V1conblockunitpsth.psth(:,:,m,k,n)=splitspk{m}{1}.Train';
            end
        end
    end
    save(['D:\data extracted\' destination(slashidx+1:end) '\V1conblockunitpsth.mat'],'V1conblockunitpsth','-v7.3')
end
end
if CONORIlabel(2)
if ~exist(['D:\data extracted\' destination(slashidx+1:end) '\V1oriblockunitpsth.mat'],'file')
    for n=1:NDftORI
        %%
%         load([destination '\' DftORIlist{n} '\ExpMonitor.mat']);
%         validtrial=1==ExpMonitor.RespCode;
%         Tstart=ExpMonitor.StartT(validtrial);
%         if max(Tstart)>=-0.2
%             trialnum=sum(0==ExpMonitor.StimCnd(validtrial));
%             V1oriblockunitpsth.psth(:,:,:,:,n)=nan(floor((conwin(2)-conwin(1))*fs*30/nsample),trialnum,6,sum(sum(valididx)));
%             continue % for n=1:NDftORI
%         end
        %%
        for k=1:sum(sum(valididx))
            spikeinfo=BinSpike([destination '\' DftORIlist{n} '.nev'],validelec(k),'t',oriwin,0.001);
            splitspk=SplitInfo(spikeinfo,[1 0 1 0]);
            for m=1:6
                V1oriblockunitpsth.psth(:,:,m,k,n)=[splitspk{m*2}{1}.Train' splitspk{m*2+12}{1}.Train'];
            end
        end
    end
    save(['D:\data extracted\' destination(slashidx+1:end) '\V1oriblockunitpsth.mat'],'V1oriblockunitpsth','-v7.3')
end
end
%%
toc
