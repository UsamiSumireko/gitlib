function reversedevelop(str,plotlabel)
%%
% str='H:\180408';
% plotlabel=0;
%%
load([str '\XYvalidid.mat']);
load([str '\ORIvalidid.mat']);
valididx=XYvalidid(1:12,:) & ORIvalidid(1:12,:);
cmppath =  'D:\SN 4566-001277.cmp';
cmpinfo=LoadCmp(cmppath,1,0);
elec=cmpinfo{1,1}.RealElec(1:12,:);
validelec=elec(valididx);
puffcnd=1:3;
neucnd=4:12;
winidx=201:780;
peakwin=231:270;
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
conwin=[-0.3 0.500];
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

if ~exist(['D:\data extracted\' str(slashidx+1:end) '\V1conblockunitpsth.mat'],'file')
    V1conblockunitpsth=struct('psth',[],'LRlist',[],'reverselabel',[]);
    for i=1:numel(Bindex)
        for j=1:numel(validelec)
            spikeinfo=BinSpike([str '\' ORICONlist{Bindex(i)} '.nev'],validelec(j),'t',conwin,0.001);
            splitspk=SplitInfo(spikeinfo,[1 0 1 0]);
            for m=1:numel(splitspk)
                V1conblockunitpsth.psth(:,:,m,i,j)=splitspk{m}{1}.Train';
            end
        end
    end
    V1conblockunitpsth.LRlist=ORICONlr(Bindex);
    V1conblockunitpsth.reverselabel=reverselabel;
    if ~exist(['D:\data extracted\' str(slashidx+1:end)],'dir')
        mkdir(['D:\data extracted\' str(slashidx+1:end)])
    end
    V1conblockunitpsth.psth=permute(V1conblockunitpsth.psth,[1 2 3 5 4]);
    save(['D:\data extracted\' str(slashidx+1:end) '\V1conblockunitpsth.mat'],'V1conblockunitpsth','-v7.3')
    V1conblockunitpsth.psth=permute(V1conblockunitpsth.psth,[1 2 3 5 4]);
else
    load(['D:\data extracted\' str(slashidx+1:end) '\V1conblockunitpsth.mat']);
    V1conblockunitpsth.psth=permute(V1conblockunitpsth.psth,[1 2 3 5 4]);
end
if ~exist(['D:\data extracted\' str(slashidx+1:end) '\V1oriblockunitpsth.mat'],'file')
    V1oriblockunitpsth=struct('psth',[]);
    for i=1:numel(DftORIlist)
        for j=1:numel(validelec)
            spikeinfo=BinSpike([str '\' DftORIlist{i} '.nev'],validelec(j),'t',oriwin,0.001);
            splitspk=SplitInfo(spikeinfo,[1 0 1 0]);
            for m=1:2
                V1oriblockunitpsth.psth(:,:,m,i,j)=[splitspk{6*m-2}{1}.Train' splitspk{6*m+10}{1}.Train'];
            end
        end
    end
    V1oriblockunitpsth.psth=permute(V1oriblockunitpsth.psth,[1 2 3 5 4]);
    save(['D:\data extracted\' str(slashidx+1:end) '\V1oriblockunitpsth.mat'],'V1oriblockunitpsth','-v7.3')
    V1oriblockunitpsth.psth=permute(V1oriblockunitpsth.psth,[1 2 3 5 4]);
else
    load(['D:\data extracted\' str(slashidx+1:end) '\V1oriblockunitpsth.mat']);
    V1oriblockunitpsth.psth=permute(V1oriblockunitpsth.psth,[1 2 3 5 4]);
end
%%
if plotlabel
    ntrial=60;
    [~,sc2,sc3,sc4,sc5]=size(V1conblockunitpsth.psth);
    [~,so2,so3,so4,so5]=size(V1oriblockunitpsth.psth);
    sc1=numel(winidx);
    so1=sc1;
    for nrev=0:numel(V1conblockunitpsth.reverselabel)
        if 0==nrev
            blockidx=1:V1conblockunitpsth.reverselabel(1)-1;
            prerevconpuffpsth=V1conblockunitpsth.psth(winidx,:,puffcnd,blockidx,:);
            prerevconneupsth=V1conblockunitpsth.psth(winidx,:,neucnd,blockidx,:);
            prerevLR=int8(V1conblockunitpsth.LRlist(blockidx(1)));
            oripuffpsth=V1oriblockunitpsth.psth(winidx,:,prerevLR+1,:,:);
            orineupsth=V1oriblockunitpsth.psth(winidx,:,2-prerevLR,:,:);
            
            figure
            subplot(1,2,1)
            xv=(1:sc1)';
            temp=reshape(prerevconpuffpsth(:,:,:,:,:),sc1,[],sc5);
            temp=squeeze(mean(temp,2));
            yste=std(temp,[],2)/sqrt(sc5);
            yste=yste(xv)/0.001;
            prerevconpuffpeak=mean(temp(peakwin,:),1);
            yv=mean(temp,2);
            yv=yv(xv)/0.001;
            patch([xv; flipud(xv)],[yv+yste; flipud(yv-yste)],[1 0.75 0.75]);
            hold on
            plot(xv,yv,'r')
            temp=reshape(oripuffpsth(:,:,:,:,:),so1,[],so5);
            temp=squeeze(mean(temp,2));
            yste=std(temp,[],2)/sqrt(so5);
            prerevoripuffpeak=mean(temp(peakwin,:),1);
            yste=yste(xv)/0.001;
            yv=mean(temp,2);
            yv=yv(xv)/0.001;
            patch([xv; flipud(xv)],[yv+yste; flipud(yv-yste)],[0.75 1 0.75]);
            plot(xv,yv,'g')
            title('puff ori')
            xlabel('time(ms)')
            ylabel('FR (Hz)')
            
            subplot(1,2,2)
            temp=reshape(prerevconneupsth(:,:,:,:,:),sc1,[],sc5);
            temp=squeeze(mean(temp,2));
            yste=std(temp,[],2)/sqrt(sc5);
            yste=yste(xv)/0.001;
            prerevconneupeak=mean(temp(peakwin,:),1);
            yv=mean(temp,2);
            yv=yv(xv)/0.001;
            patch([xv; flipud(xv)],[yv+yste; flipud(yv-yste)],[1 0.75 0.75]);
            hold on
            plot(xv,yv,'r')
            temp=reshape(orineupsth(:,:,:,:,:),so1,[],so5);
            temp=squeeze(mean(temp,2));
            yste=std(temp,[],2)/sqrt(so5);
            prerevorineupeak=mean(temp(peakwin,:),1);
            yste=yste(xv)/0.001;
            yv=mean(temp,2);
            yv=yv(xv)/0.001;
            patch([xv; flipud(xv)],[yv+yste; flipud(yv-yste)],[0.75 1 0.75]);
            plot(xv,yv,'g')
            title('neu ori')
            xlabel('time(ms)')
            ylabel('FR (Hz)')
            
            figure
            bar([2 4],[mean(prerevconpuffpeak,2)-mean(prerevoripuffpeak,2) mean(prerevconneupeak,2)-mean(prerevorineupeak,2)]/0.001);
            ylabel('peak FR diff (Hz)')
        else
            if numel(V1conblockunitpsth.reverselabel)==nrev
                blockidx=V1conblockunitpsth.reverselabel(nrev):numel(ORICONlist);
            else
                blockidx=V1conblockunitpsth.reverselabel(nrev):V1conblockunitpsth.reverselabel(nrev+1)-1;
            end
            conpuffpsth{nrev}=reshape(V1conblockunitpsth.psth(winidx,:,puffcnd,blockidx,:),sc1,[],sc5);
            conpuffpsth{nrev}=reshape(conpuffpsth{nrev},sc1,ntrial/4,[],sc5);
            conneupsth{nrev}=reshape(V1conblockunitpsth.psth(winidx,:,neucnd,blockidx,:),sc1,[],sc5);
            conneupsth{nrev}=reshape(conneupsth{nrev},sc1,ntrial/4*3,[],sc5);
            LR=int8(V1conblockunitpsth.LRlist(blockidx(1)));
            oripuffpsth=V1oriblockunitpsth.psth(winidx,:,LR+1,:,:);
            orineupsth=V1oriblockunitpsth.psth(winidx,:,2-LR,:,:);
            temp=reshape(oripuffpsth(xv,:,:,:,:),so1,[],so5);
            temp=squeeze(mean(temp,2));
            oripuffyste=std(temp,[],2)/sqrt(so5);
            oripuffyv=mean(temp,2);
            oripuffpk=mean(temp(peakwin,:),1);
            temp=reshape(orineupsth(xv,:,:,:,:),so1,[],so5);
            temp=squeeze(mean(temp,2));
            orineuyste=std(temp,[],2)/sqrt(so5);
            orineuyv=mean(temp,2);
            orineupk=mean(temp(peakwin,:),1);
            
            nsession=size(conpuffpsth{nrev},3);
            conneupk=[];
            conpuffpk=[];
            for k=1:nsession
                if 0==mod(k-1,4)
                    figure
                end
                subplot(2,4,mod(k-1,4)+1)
                xv=(1:sc1)';
                temp=mean(conpuffpsth{nrev}(:,:,k,:),2);
                yv=squeeze(temp);
                yste=std(yv,[],2)/sqrt(sc5);
                yste=yste(xv);
                conpuffpk(:,k)=squeeze(mean(yv(peakwin,:),1));
                yv=mean(yv,2);
                yv=yv(xv);
                patch([xv; flipud(xv)],[oripuffyv+oripuffyste; flipud(oripuffyv-oripuffyste)],[0.75 1 0.75])
                hold on
                plot(xv,oripuffyv,'g')
                patch([xv; flipud(xv)],[yv+yste;flipud(yv-yste)],[1 0.75 0.75])
                plot(xv,yv,'r')
                title('puff ori')
                
                subplot(2,4,mod(k-1,4)+5)
                temp=mean(conneupsth{nrev}(:,:,k,:),2);
                yv=squeeze(temp);
                yste=std(yv,[],2)/sqrt(sc5);
                yste=yste(xv);
                conneupk(:,k)=squeeze(mean(yv(peakwin,:),1));
                yv=mean(yv,2);
                yv=yv(xv);
                patch([xv; flipud(xv)],[orineuyv+orineuyste; flipud(orineuyv-orineuyste)],[0.75 1 0.75])
                hold on
                plot(xv,orineuyv,'g')
                patch([xv; flipud(xv)],[yv+yste;flipud(yv-yste)],[1 0.75 0.75])
                plot(xv,yv,'r')
                title('neu ori')
            end
            figure
            yv=mean(conpuffpk,1)-repmat(mean(oripuffpk),1,nsession)-mean(conneupk,1)+repmat(mean(orineupk),1,nsession);
            plot(1:nsession,yv/0.001,'bo')
            xlabel('block number')
            ylabel('peak FR diff')
        end
    end
end