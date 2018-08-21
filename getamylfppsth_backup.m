function getamylfppsth(destination,LR,filtlabel)
% destination='I:\160720';
% LR='R';
tic
nPuffcnd=3;
cmppath =  'D:\SN 4566-001277.cmp';
% cmppath = 'I:\SN 4566-001303+Probe MK2.cmp';
cmpinfo = LoadCmp(cmppath,1,0);
elec=cmpinfo{1,1}.RealElec(13:15,:);
elec=reshape(elec',[],1);
load('D:\data extracted\LFPMedium.mat');
conwin=[-0.1 1.1];
oriwin=[-0.1 0.5];

fs=1000;
slashidx=strfind(destination,'\');
slashidx=slashidx(end);
baselinewin=[-0.1 0.1];
respwin=[0.15 0.35];
nbootstrap=1000;
nsample=30;
% cnd=(ExpMonitor.StimID>=3)*3+mod(ExpMonitor.StimID,3)+1;
NDftORI=0;
NORICON=0;
namelist=dir(destination);
DftORIlist=[];
ORICONlist=[];
for i=3:numel(namelist)
    if namelist(i).isdir && strcmp(namelist(i).name(1:6),'DftORI')
        NDftORI=NDftORI+1;
        DftORIlist=[DftORIlist; {namelist(i).name}];
    elseif namelist(i).isdir && strcmp(namelist(i).name(1:6),'ORICON')
        NORICON=NORICON+1;
        ORICONlist=[ORICONlist; {namelist(i).name}];
    end
end
%%
strr='AmyConLFP'
if ~exist(['D:\data extracted\' destination(slashidx+1:end) '\amyconblockunitlfp_' num2str(nsample) '.mat'],'file')
    amyconblockunitlfp=struct('timelabel',[],'lfp',[]);
    amyconblockunitlfp.lfp=nan(floor((conwin(2)-conwin(1))*fs*30/nsample),10,12,24,NORICON);
    amyconblockunitlfp.LR=LR;
    
    for n=1:numel(ORICONlist)
        n,
        for k=1:24
            sprintf(repmat('.',1,k))
            path=[destination '\' ORICONlist{n}];
            lfpinfo=TruncateLFP([path '.ns5'],elec(k),conwin,'MKII');
            [splitlfp]=SplitInfo(lfpinfo, [1 0 1 0]);
            if isempty(amyconblockunitlfp.timelabel)
                temp=downsample(lfpinfo.TimeLabel,30);
                amyconblockunitlfp.timelabel=temp(2:end);
            end
            for m=1:numel(splitlfp)
                lfptemp=double(splitlfp{m}{1}.LFP');
                if filtlabel
                    lfptemp=filtfilt(SOS,G,lfptemp);
                end
                dslfp=downsample(lfptemp,nsample);
                amyconblockunitlfp.lfp(:,:,m,k,n)=dslfp(2:end,:);
            end
        end
    end
    if ~exist(['D:\data extracted\' destination(slashidx+1:end)],'dir')
        mkdir(['D:\data extracted\' destination(slashidx+1:end)])
    end
    %     save(['D:\data extracted\' destination(end-5:end) '\amyconblockunitlfp_' num2str(nsample) '.mat'],'amyconblockunitlfp','-v7.3');
else
    load(['D:\data extracted\' destination(slashidx+1:end) '\amyconblockunitlfp_' num2str(nsample) '.mat']);
end
%%
strr='AmyOriLFP'
if ~exist(['D:\data extracted\' destination(slashidx+1:end)  '\amyoriblockunitlfp_' num2str(nsample) '.mat'],'file')
    amyoriblockunitlfp=struct('timelabel',[],'lfp',[]);
    amyoriblockunitlfp.lfp=nan(floor((oriwin(2)-oriwin(1))*fs*30/nsample),20,6,24,NDftORI);
    for n=1:numel(DftORIlist)
        n,
        for k=1:24
            sprintf(repmat('.',1,k))
            path=[destination '\' DftORIlist{n}];
            lfpinfo=TruncateLFP([path '.ns5'],elec(k),oriwin,'MKII');
            [splitlfp]=SplitInfo(lfpinfo, [1 0 1 0]);
            if isempty(amyoriblockunitlfp.timelabel)
                temp=downsample(lfpinfo.TimeLabel,30);
                amyoriblockunitlfp.timelabel=temp(2:end);
            end
            for m=1:2*nPuffcnd
                lfptemp=double([splitlfp{2*m}{1}.LFP' splitlfp{2*m+12}{1}.LFP']);
                if filtlabel
                    lfptemp=filtfilt(SOS,G,lfptemp);
                end
                dslfp=downsample(lfptemp,nsample);
                %                 dslfp=[downsample(splitlfp{2*m}{1,1}.LFP',nsample) downsample(splitlfp{2*m+12}{1,1}.LFP',nsample)];
                amyoriblockunitlfp.lfp(:,:,m,k,n)=dslfp(2:end,:);
            end
        end
    end
    if ~exist(['D:\data extracted\' destination(slashidx+1:end)],'dir')
        mkdir(['D:\data extracted\' destination(slashidx+1:end)])
    end
    %     save(['D:\data extracted\' destination(end-5:end)  '\amyoriblockunitlfp_' num2str(nsample) '.mat'],'amyoriblockunitlfp','-v7.3');
else
    load(['D:\data extracted\' destination(slashidx+1:end)  '\amyoriblockunitlfp_' num2str(nsample) '.mat']);
end
toc
%%
tic
strr='AmyCONPSTH'
if ~exist(['D:\data extracted\' destination(slashidx+1:end) '\amyconblockunitpsth.mat'],'file')
    amyconblockunitpsth=struct('psth',[],'LR',LR,'spikevalid',zeros(24,NORICON),'stimvalid',zeros(24,NORICON));
    amyconblockunitpsth.psth=nan(floor((conwin(2)-conwin(1))/0.001),10,12,24,NORICON);
    resptidx=floor((respwin(1)-conwin(1))*fs)+1:floor((respwin(2)-conwin(1))*fs);
    bltidx=floor((baselinewin(1)-conwin(1))*fs)+1:floor((baselinewin(2)-conwin(1))*fs);
    for n=1:numel(ORICONlist)
        for k=1:24
            path=[destination '\' ORICONlist{n} ];
            if exist([path '\Elec' num2str(elec(k)) 'Quality.mat'],'file')
                load([path '\Elec' num2str(elec(k)) 'Quality.mat'])
                [~,snridx]=max(Quality.SNR);
                amyconblockunitpsth.spikevalid(k,n)=Quality.SNR(snridx)>1;
                amyconblockunitlfp.spikevalid(k,n)=Quality.SNR(snridx)>1;
                spikeinfo=BinSpike([path '.nev'],elec(k),char(96+snridx),conwin,0.001);
            else
                spikeinfo=BinSpike([path '.nev'],elec(k),'t',conwin,0.001);
            end
            [splitspike]=SplitInfo(spikeinfo,[1 0 1 0]);
            mergepsth=[];
            diff=nan(nbootstrap,1);
            for m=1:numel(splitspike)
                amyconblockunitpsth.psth(:,:,m,k,n)=splitspike{m}{1}.Train';
                mergepsth=[mergepsth splitspike{m}{1}.Train'];
            end
            pool=[mean(mergepsth(bltidx,:)) mean(mergepsth(resptidx,:))];
            for x=1:nbootstrap
                ridx=randperm(2*size(mergepsth,2));
                diff(x)=mean(pool(ridx(1:size(mergepsth,2))))-mean(pool(ridx(1+size(mergepsth,2):end)));
            end
            diffstd=std(diff);%/sqrt(size(mergepsth,2));
            Rdiff=mean(mean(mergepsth(resptidx,:)))-mean(mean(mergepsth(bltidx,:)));
            amyconblockunitpsth.stimvalid(k,n)=Rdiff>1.96*diffstd;
            amyconblockunitlfp.stimvalid(k,n)=Rdiff>1.96*diffstd;
            
            ht=ttest(mean(mergepsth(bltidx,:)),mean(mergepsth(resptidx,:)));
            if isnan(ht)
                amyconblockunitpsth.stimvalid_ttest(k,n)=false;
                amyconblockunitlfp.stimvalid_ttest(k,n)=false;
            else
                amyconblockunitpsth.stimvalid_ttest(k,n)=logical(ht);
                amyconblockunitlfp.stimvalid_ttest(k,n)=logical(ht);
            end
        end
    end
    save(['D:\data extracted\' destination(slashidx+1:end) '\amyconblockunitpsth.mat'],'amyconblockunitpsth','-v7.3')
end

strr='AmyORIPSTH'
if ~exist(['D:\data extracted\' destination(slashidx+1:end) '\amyoriblockunitpsth.mat'],'file')
    amyoriblockunitpsth=struct('psth',[],'LR',LR,'spikevalid',zeros(24,NDftORI),'stimvalid',zeros(24,NDftORI));
    amyoriblockunitpsth.psth=nan(floor((oriwin(2)-oriwin(1))/0.001),20,6,24,NDftORI);
    resptidx=floor((respwin(1)-oriwin(1))*fs)+1:floor((respwin(2)-oriwin(1))*fs);
    bltidx=floor((baselinewin(1)-oriwin(1))*fs)+1:floor((baselinewin(2)-oriwin(1))*fs);    
    for n=1:NDftORI
        for k=1:24
            path=[destination '\' DftORIlist{n} ];
            if exist([path '\Elec' num2str(elec(k)) 'Quality.mat'],'file')
                load([path '\Elec' num2str(elec(k)) 'Quality.mat'])
                [~,snridx]=max(Quality.SNR);
                amyoriblockunitpsth.spikevalid(k,n)=Quality.SNR(snridx)>1;
                amyoriblockunitlfp.spikevalid(k,n)=Quality.SNR(snridx)>1;
                spikeinfo=BinSpike([path '.nev'],elec(k),char(96+snridx),oriwin,0.001);
            else
                spikeinfo=BinSpike([path '.nev'],elec(k),'t',oriwin,0.001);
            end
            [splitspike]=SplitInfo(spikeinfo,[1 0 1 0]);
            mergepsth=[];
            mergecs=[];
            mergeneu=[];
            diff=nan(nbootstrap,1);
            for m=1:6
                amyoriblockunitpsth.psth(:,:,m,k,n)=[splitspike{m*2}{1}.Train' splitspike{12+m*2}{1}.Train'];
                mergepsth=[mergepsth splitspike{m*2}{1}.Train' splitspike{12+m*2}{1}.Train'];
            end
            
            pool=[mean(mergepsth(bltidx,:)) mean(mergepsth(resptidx,:))];
            for x=1:nbootstrap
                ridx=randperm(2*size(mergepsth,2));
                diff(x)=mean(pool(ridx(1:size(mergepsth,2))))-mean(pool(ridx(1+size(mergepsth,2):end)));
            end
            diffstd=std(diff);%/sqrt(size(mergepsth,2));
            Rdiff=mean(mean(mergepsth(resptidx,:)))-mean(mean(mergepsth(bltidx,:)));
            amyoriblockunitpsth.stimvalid(k,n)=Rdiff>1.96*diffstd;
            amyoriblockunitlfp.stimvalid(k,n)=Rdiff>1.96*diffstd;
            ht=ttest(mean(mergepsth(bltidx,:)),mean(mergepsth(resptidx,:)));
            if isnan(ht)
                amyoriblockunitpsth.stimvalid(k,n)=false;
                amyoriblockunitlfp.stimvalid(k,n)=false;
            else
                amyoriblockunitpsth.stimvalid(k,n)=logical(ht);
                amyoriblockunitlfp.stimvalid(k,n)=logical(ht);
            end
        end
    end
    save(['D:\data extracted\' destination(slashidx+1:end) '\amyoriblockunitpsth.mat'],'amyoriblockunitpsth','-v7.3')
end
save(['D:\data extracted\' destination(slashidx+1:end) '\amyconblockunitlfp_' num2str(nsample) '.mat'],'amyconblockunitlfp','-v7.3');
save(['D:\data extracted\' destination(slashidx+1:end)  '\amyoriblockunitlfp_' num2str(nsample) '.mat'],'amyoriblockunitlfp','-v7.3');
toc
