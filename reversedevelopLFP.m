function reversedevelopLFP(str,plotlabel)
%%
% plotlabel=0;
% str='H:\180408';
%%
filtlabel=1;
load([str '\XYvalidid.mat']);
load([str '\ORIvalidid.mat']);
load('D:\data extracted\LFPMedium.mat');
valididx=XYvalidid(1:12,:) & ORIvalidid(1:12,:);
cmppath =  'D:\SN 4566-001277.cmp';
cmpinfo=LoadCmp(cmppath,1,0);
elec=cmpinfo{1,1}.RealElec(1:12,:);
validelec=elec(valididx);
puffcnd=1:3;
neucnd=4:12;
nsample=30;
powerband=[{[8 26]}; {[26 54]}; {[54 80]};];
params.Fs = 1000; % sampling frequency
params.tapers = [3 5]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 2;
params.err =[2 0.05]; 
params.sldwinlable=0;
movingwin=[0.2 0.005];
params.movingwin=movingwin;

slashidx=strfind(str,'\');
slashidx=slashidx(end);
namelist=dir(str);
ORICONnum=[];
ORICONlist=[];
ORICONlr=[];
DftORIlist=[];
concount=0;
oricount=0;
oriwin=[-0.3 0.500];
conwin=[-0.3 1.100];
for i=3:numel(namelist)
    if strcmp(namelist(i).name(1:6),'ORICON') && namelist(i).isdir
        concount=concount+1;
        ORICONnum(concount)=str2num(namelist(i).name(end-2:end));
        ORICONlist=[ORICONlist; {namelist(i).name}];
        ORICONlr(concount)=strcmpi(namelist(i).name(7),'l');% R=0 L=1
    elseif strcmp(namelist(i).name(1:6),'DftORI') && namelist(i).isdir
        oricount=oricount+1;
        DftORIlist=[DftORIlist; {namelist(i).name}];
    end
end
[~,Bindex]=sort(ORICONnum);
temp=ORICONlr(Bindex);
reverselabel=find((temp(1:end-1)-temp(2:end))~=0);
reverselabel=reverselabel+1;
%%
if ~exist(['D:\data extracted\' str(slashidx+1:end) '\V1conblockunitlfp.mat'],'file')
    V1conblockunitlfp=struct('lfp',[],'LRlist',[],'reverselabel',[],'timelabel',[]);
    for i=1:numel(Bindex)
        i,
        load([str '\' ORICONlist{Bindex(i)} '\LFPBasicMKII.mat']);
        LFPBasic.Fs=30000;
        save([str '\' ORICONlist{Bindex(i)} '\LFPBasicMKII.mat'],'LFPBasic','-v7.3');
        for j=1:numel(validelec)
                sprintf(repmat('.',1,j))
                lfpinfo=TruncateLFP([str '\' ORICONlist{Bindex(i)} '.ns5'],validelec(j),conwin,'MKII');
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
                    V1conblockunitlfp.lfp(:,:,m,i,j)=dslfp(2:end,:);
                end
        end
    end
    V1conblockunitlfp.LRlist=ORICONlr(Bindex);
    V1conblockunitlfp.reverselabel=reverselabel;
    if ~exist(['D:\data extracted\' str(slashidx+1:end)],'dir')
        mkdir(['D:\data extracted\' str(slashidx+1:end)])
    end
    V1conblockunitlfp.lfp=permute(V1conblockunitlfp.lfp,[ 1 2 3 5 4]);
    save(['D:\data extracted\' str(slashidx+1:end) '\V1conblockunitlfp.mat'],'V1conblockunitlfp','-v7.3')
    V1conblockunitlfp.lfp=permute(V1conblockunitlfp.lfp,[ 1 2 3 5 4]);
else
    load(['D:\data extracted\' str(slashidx+1:end) '\V1conblockunitlfp.mat']);
end
if ~exist(['D:\data extracted\' str(slashidx+1:end) '\V1oriblockunitlfp.mat'],'file')
    V1oriblockunitlfp=struct('lfp',[],'timelabel',[]);
    for i=1:numel(DftORIlist)
        i,
        for j=1:numel(validelec)
            sprintf(repmat('.',1,j))
            lfpinfo=TruncateLFP([str '\' DftORIlist{i} '.ns5'],validelec(j),oriwin,'MKII');
            [splitlfp]=SplitInfo(lfpinfo,[1 0 1 0]);
            if isempty(V1oriblockunitlfp.timelabel)
                temp=downsample(lfpinfo.TimeLabel,nsample);
                V1oriblockunitlfp.timelabel=temp(2:end);
            end
            for m=1:2
                lfptemp=double([splitlfp{6*m-2}{1}.LFP' splitlfp{6*m+10}{1}.LFP']);
                if filtlabel
                    lfptemp=filtfilt(SOS,G,lfptemp);
                end
                dslfp=downsample(lfptemp,nsample);
                V1oriblockunitlfp.lfp(:,:,m,i,j)=dslfp(2:end,:);
            end
        end
    end
    V1oriblockunitlfp.lfp=permute(V1oriblockunitlfp.lfp,[1 2 3 5 4]);
    save(['D:\data extracted\' str(slashidx+1:end) '\V1oriblockunitlfp.mat'],'V1oriblockunitlfp','-v7.3')
    V1oriblockunitlfp.lfp=permute(V1oriblockunitlfp.lfp,[1 2 3 5 4]);
else
    load(['D:\data extracted\' str(slashidx+1:end) '\V1oriblockunitlfp.mat']);
end
%%
if plotlabel
    ntrial=60;
    prerevconS1=[];
    prerevconS2=[];
    conS1=[];
    conS2=[];
    oriS1=[];
    oriS2=[];
    [sc1,sc2,sc3,sc4,sc5]=size(V1conblockunitlfp.lfp);
    [so1,so2,so3,so4,so5]=size(V1oriblockunitlfp.lfp);
    winidx=401:680;
    Llfp=numel(winidx);
    
    for nrev=0:numel(V1conblockunitlfp.reverselabel)
        if 0==nrev
            blockidx=1:V1conblockunitlfp.reverselabel(1)-1;
            prerevconpufflfp=V1conblockunitlfp.lfp(winidx,:,puffcnd,blockidx,:);
            temp=reshape(prerevconpufflfp,Llfp,[]);
            temp=rmlinesc(temp,params,[],'n',[50 100]);
            prerevconpufflfp=reshape(temp,Llfp,[],sc5);
            
            prerevconneulfp=V1conblockunitlfp.lfp(winidx,:,neucnd,blockidx,:);
            temp=reshape(prerevconneulfp,Llfp,[]);
            temp=rmlinesc(temp,params,[],'n',[50 100]);
            prerevconneulfp=reshape(temp,Llfp,[],sc5);
            prerevLR=int8(V1conblockunitlfp.LRlist(blockidx(1)));
            
            oripufflfp=V1oriblockunitlfp.lfp(winidx,:,prerevLR+1,:,:);
            temp=reshape(oripufflfp,Llfp,[]);
            temp=rmlinesc(temp,params,[],'n',[50 100]);
            oripufflfp=reshape(temp,Llfp,[],so5);
            
            orineulfp=V1oriblockunitlfp.lfp(winidx,:,2-prerevLR,:,:);
            temp=reshape(orineulfp,Llfp,[]);
            temp=rmlinesc(temp,params,[],'n',[50 100]);
            orineulfp=reshape(temp,Llfp,[],so5);
            
            for ncell=1:sc5
                [S,f]=mtspectrumc(prerevconpufflfp(:,:,ncell),params);
                prerevconS1(:,ncell)=S;
                [S,~]=mtspectrumc(prerevconneulfp(:,:,ncell),params);
                prerevconS2(:,ncell)=S;
                [S,~]=mtspectrumc(oripufflfp(:,:,ncell),params);
                oriS1(:,ncell)=S;
                [S,~]=mtspectrumc(orineulfp(:,:,ncell),params);
                oriS2(:,ncell)=S;
            end
            notnanidx=0==sum(isnan(prerevconS1),1);
            prerevconS1=prerevconS1(:,notnanidx);
            prepuffchange=prerevconS1./oriS1(:,notnanidx);
            
            notnanidx=0==sum(isnan(prerevconS2),1);
            prerevconS2=prerevconS2(:,notnanidx);
            preneuchange=prerevconS2./oriS2(:,notnanidx);
            
            figure
            subplot(1,2,1)
            xv=f;
            yv=squeeze(mean(prerevconS1,2));
            yv=log(yv);
            plot(xv,yv,'r');
            hold on
            yv=squeeze(mean(prerevconS2,2));
            yv=log(yv);
            plot(xv,yv,'g');
            title('Prerev Con')
            subplot(1,2,2)
            yv=squeeze(mean(oriS1,2));
            yv=log(yv);
            plot(xv,yv,'r');
            hold on
            yv=squeeze(mean(oriS2,2));
            yv=log(yv);
            plot(xv,yv,'g');
            title('Tun Blocks')
            
            figure
            plot(f,mean(prepuffchange,2)-mean(preneuchange,2))
        else
            if numel(V1conblockunitlfp.reverselabel)==nrev
                blockidx=V1conblockunitlfp.reverselabel(nrev):numel(ORICONlist);
            else
                blockidx=V1conblockunitlfp.reverselabel(nrev):V1conblockunitlfp.reverselabel(nrev+1)-1;
            end
            pufflfp{nrev}=reshape(V1conblockunitlfp.lfp(winidx,:,puffcnd,blockidx,:),Llfp,[],sc5);
            pufflfp{nrev}=reshape(pufflfp{nrev},Llfp,ntrial/4,[],sc5);
            neulfp{nrev}=reshape(V1conblockunitlfp.lfp(winidx,:,neucnd,blockidx,:),Llfp,[],sc5);
            neulfp{nrev}=reshape(neulfp{nrev},Llfp,ntrial/4*3,[],sc5);
            LR=int8(V1conblockunitlfp.LRlist(blockidx(1)));
            conS1=[];
            conS2=[];
            neuchange=[];
            puffchange=[];
            revpowerchangeidx=[];
            for nsession=1:size(pufflfp{nrev},3)
                if 0==mod(nsession-1,4)
                    figure
                end
                subplot(2,2,mod(nsession-1,4)+1)
                for mcell=1:sc5
                    temp=squeeze(pufflfp{nrev}(:,:,nsession,mcell));
                    [S,f]=mtspectrumc(temp,params);
                    conS1(:,mcell,nsession)=S;
                    temp=squeeze(neulfp{nrev}(:,:,nsession,mcell));
                    [S,~]=mtspectrumc(temp,params);
                    conS2(:,mcell,nsession)=S;
                end
                xv=f;
                notnanidx=0==sum(isnan(conS1(:,:,nsession)),1);
                yv=squeeze(mean(conS1(:,notnanidx,nsession),2));
                for mband=1:numel(powerband)
                    pidx=f<=powerband{mband}(2) & f>powerband{mband}(1);
                    puffchange{mband,nsession}=conS1(pidx,notnanidx,nsession)./oriS2(pidx,notnanidx);
                end
                yv=log(yv);
                plot(xv,yv,'r')
                hold on
                notnanidx=0==sum(isnan(conS2(:,:,nsession)),1);
                yv=squeeze(mean(conS2(:,notnanidx,nsession),2));
                for mband=1:numel(powerband)
                    pidx=f<=powerband{mband}(2) & f>powerband{mband}(1);
                    neuchange{mband,nsession}=conS2(pidx,notnanidx,nsession)./oriS1(pidx,notnanidx);
                end
                yv=log(yv);
                plot(xv,yv,'g')
                title(['Con Block ' num2str(nsession)]);
            end
            for nsession=1:size(pufflfp{nrev},3)
                for mband=1:numel(powerband)
                    revpowerchangeidx(nsession,mband)=mean(mean(puffchange{mband,nsession}-neuchange{mband,nsession}));
                end
            end
            figure
            plot(1:nsession,revpowerchangeidx);
        end
    end
end
