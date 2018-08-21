%   MKII version use high dimension array data structure
clear
folderlist={
%     'D:\data extracted\180111';
%     'D:\data extracted\180112';
%     'D:\data extracted\180117';
%     'D:\data extracted\180118';
%     'D:\data extracted\180119';
%     'D:\data extracted\180120';
%     'D:\data extracted\180124';
%     'D:\data extracted\180125';
%     'D:\data extracted\180126';

%     'D:\data extracted\180227';
%     'D:\data extracted\180228';
%     'D:\data extracted\180301';
%     'D:\data extracted\180302';

    'D:\data extracted\180425';
    'D:\data extracted\180426';
    'D:\data extracted\180427';

%     'D:\data extracted\180504';
%     'D:\data extracted\180505';
%     'D:\data extracted\180506';
    };

selecti='both';

% LRlist=['ll' 'rrrr' 'llll' 'rrrr' 'rrr'];
LRlist=['rrr'];
RFlabel=1;  %0 contra RF 1 on RF
excude1sttunlabel=0;
ploteachcelllabel=0;

preferwinr=22.5;
bin=0.001;
win=[-0.1 0.5];
winidx=1:600;
blwin=-win(1)/bin+1:(0.1-win(1))/bin;
pkwin=(0.1-win(1))/bin+(0.030/bin+1:0.070/bin);
smoothlength=0.01/bin;
conpsth=cell(4,1);
oripsth=cell(4,1);
allconpsth=cell(4,1);
alloripsth=cell(4,1);
Lpcellconpsth=cell(4,1);
Rpcellconpsth=cell(4,1);
Mpcellconpsth=cell(4,1);
Lpcelloripsth=cell(4,1);
Rpcelloripsth=cell(4,1);
Mpcelloripsth=cell(4,1);
mergeLpconpsth=cell(numel(folderlist),1);
mergeRpconpsth=cell(numel(folderlist),1);
mergeMpconpsth=cell(numel(folderlist),1);
mergeLporipsth=cell(numel(folderlist),1);
mergeRporipsth=cell(numel(folderlist),1);
mergeMporipsth=cell(numel(folderlist),1);

for iday=1:numel(folderlist)
    iday
    if excude1sttunlabel
        namelist=dir(['H:\' folderlist{iday}(end-5:end)]);
        numlist=[];
        for i=3:numel(namelist)
            if namelist(i).isdir && strcmp(namelist(i).name(1:9),'DiGratTun')
                numlist=[numlist str2num(namelist(i).name(10:12))];
            end
        end
        blockinterval=numlist(2:end)-numlist(1:end-1);
        poststartnum=find(blockinterval>1)+1;
        assert(numel(poststartnum)==1);        
    end
    
    lr=LRlist(iday);
    if strcmp(lr,'l')
        if RFlabel
            CNDlist=[3 1 1 1 4 2 2 2];
        else
            CNDlist=[2 1 1 1 4 3 3 3];
        end
    elseif strcmp(lr,'r')
        if RFlabel
            CNDlist=[1 3 3 3 2 4 4 4];
        else
            CNDlist=[1 2 2 2 3 4 4 4];            
        end            
    else        
        error('LR input is wrong, L or R');
    end
    load([folderlist{iday} '\V1conblockunitpsth.mat'])
    load([folderlist{iday} '\V1oriblockunitpsth.mat'])
    if RFlabel
        V1oriblockunitpsth.psth=V1oriblockunitpsth.psth(:,:,[1 3 2 4],:,:);
    end
    load([folderlist{iday} '\ORIvalidid.mat'])
    load([folderlist{iday} '\XYvalidid.mat'])
    load([folderlist{iday} '\Oparam15.mat'])    
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
    perefori=mod(Oparam15{4}(1:12,:),180);
    Lprefid=perefori>135-preferwinr & perefori<135+preferwinr;
    Rprefid=perefori>45-preferwinr & perefori<45+preferwinr;
    Vprefid=perefori>90-preferwinr & perefori<90+preferwinr;
    
    [s1,s2,s3,s4]=size(V1conblockunitpsth.secstimlabel);
    psthtemp=V1conblockunitpsth.psth(1:600,logical(V1conblockunitpsth.secstimlabel));
%     psthtemp=reshape(psthtemp,800,9,s2,s3,s4);
    psthtemp=reshape(psthtemp,600,sum(V1conblockunitpsth.secstimlabel(:,1,1,1)),s2,s3,s4);
    temp=permute(psthtemp,[1 2 5 3 4]);
    for m=1:4
    conpsth{m}=reshape(temp(:,:,:,m==CNDlist,:),size(temp,1),[],size(temp,5));
    vidx=conpsth{m}==65535;
    tridx=squeeze(sum(sum(vidx,3),1))==0;
    conpsth{m}=conpsth{m}(:,tridx,:);
    end
    if excude1sttunlabel
        temp=permute(V1oriblockunitpsth.psth(1:600,:,:,:,[2:poststartnum-1 poststartnum+1:end]),[1 2 5 3 4]);
    else
        temp=permute(V1oriblockunitpsth.psth(1:600,:,:,:,:),[1 2 5 3 4]);
    end
    for m=1:4
        oripsth{m}=reshape(temp(:,:,:,m,:),size(temp,1),[],size(temp,5));
        vidx=oripsth{m}==65535;
        tridx=squeeze(sum(sum(vidx,3),1))==0;
        oripsth{m}=oripsth{m}(:,tridx,:);
    end
    Lprefcell=find(Lprefid(valididx));
    Rprefcell=find(Rprefid(valididx));
    Mprefcell=find(Vprefid(valididx));
    for m=1:4
        Lpcellconpsth{m}=cat(2,Lpcellconpsth{m},squeeze(mean(conpsth{m}(:,:,Lprefcell),2)));
        Rpcellconpsth{m}=cat(2,Rpcellconpsth{m},squeeze(mean(conpsth{m}(:,:,Rprefcell),2)));
        Mpcellconpsth{m}=cat(2,Mpcellconpsth{m},squeeze(mean(conpsth{m}(:,:,Mprefcell),2)));
        
        Lpcelloripsth{m}=cat(2,Lpcelloripsth{m},squeeze(mean(oripsth{m}(:,:,Lprefcell),2)));
        Rpcelloripsth{m}=cat(2,Rpcelloripsth{m},squeeze(mean(oripsth{m}(:,:,Rprefcell),2)));
        Mpcelloripsth{m}=cat(2,Mpcelloripsth{m},squeeze(mean(oripsth{m}(:,:,Mprefcell),2)));
        
        mergeLpconpsth{iday}{m}=conpsth{m}(:,:,Lprefcell);
        mergeRpconpsth{iday}{m}=conpsth{m}(:,:,Rprefcell);
        mergeMpconpsth{iday}{m}=conpsth{m}(:,:,Mprefcell);
        
        mergeLporipsth{iday}{m}=oripsth{m}(:,:,Lprefcell);
        mergeRporipsth{iday}{m}=oripsth{m}(:,:,Rprefcell);
        mergeMporipsth{iday}{m}=oripsth{m}(:,:,Mprefcell);
    end
    if ploteachcelllabel
        figure
        for ncell=1:numel(Lprefcell)
            subplot(6,6,ncell)
            xv=-199:400;
            yv=mean([Lpcellconpsth{3}(1:600,ncell) Lpcellconpsth{4}(1:600,ncell)],2);
            plot(xv,yv,'r')
            hold on
            yv=mean([Lpcelloripsth{3}(1:600,ncell) Lpcelloripsth{4}(1:600,ncell)],2);
            plot(xv,yv,'g')
        end
        figure
        for ncell=1:numel(Rprefcell)
            subplot(6,6,ncell)
            xv=-199:400;
            yv=mean([Rpcellconpsth{1}(1:600,ncell) Rpcellconpsth{2}(1:600,ncell)],2);
            plot(xv,yv,'r')
            hold on
            yv=mean([Rpcelloripsth{1}(1:600,ncell) Rpcelloripsth{2}(1:600,ncell)],2);
            plot(xv,yv,'g')
        end
    end
    
end
smLpconpsth=cell(4,1);
smRpconpsth=cell(4,1);
smMpconpsth=cell(4,1);

smLporipsth=cell(4,1);
smRporipsth=cell(4,1);
smMporipsth=cell(4,1);
for m=1:4
    smLpconpsth{m}=convn(Lpcellconpsth{m},ones(smoothlength,1),'same')/smoothlength;
    smRpconpsth{m}=convn(Rpcellconpsth{m},ones(smoothlength,1),'same')/smoothlength;
    smMpconpsth{m}=convn(Mpcellconpsth{m},ones(smoothlength,1),'same')/smoothlength;
    smLporipsth{m}=convn(Lpcelloripsth{m},ones(smoothlength,1),'same')/smoothlength;
    smRporipsth{m}=convn(Rpcelloripsth{m},ones(smoothlength,1),'same')/smoothlength;
    smMporipsth{m}=convn(Mpcelloripsth{m},ones(smoothlength,1),'same')/smoothlength;
    
    Lpconste{m}=std(smLpconpsth{m}/bin,0,2)/sqrt(size(smLpconpsth{m},2));
    Lpconste{m}=Lpconste{m}';
    Rpconste{m}=std(smRpconpsth{m}/bin,0,2)/sqrt(size(smRpconpsth{m},2));
    Rpconste{m}=Rpconste{m}';
    Mpconste{m}=std(smMpconpsth{m}/bin,0,2)/sqrt(size(smMpconpsth{m},2));
    Mpconste{m}=Mpconste{m}';
    
    Lporiste{m}=std(smLporipsth{m}/bin,0,2)/sqrt(size(smLporipsth{m},2));
    Lporiste{m}=Lporiste{m}';
    Rporiste{m}=std(smRporipsth{m}/bin,0,2)/sqrt(size(smRporipsth{m},2));
    Rporiste{m}=Rporiste{m}';
    Mporiste{m}=std(smMporipsth{m}/bin,0,2)/sqrt(size(smMporipsth{m},2));
    Mporiste{m}=Mporiste{m}';
end

figure
for m=1:4
    subplot(2,2,m)
    vx=(win(1)+bin:bin:win(2))-0.1;
    vy=mean(smLpconpsth{m}(winidx,:)/bin,2);
    vy=vy';
    patch([vx fliplr(vx) vx(1)],[vy+Lpconste{m}(winidx) fliplr(vy-Lpconste{m}(winidx)) vy(1)+Lpconste{m}(1)],[1 0.8 0.8])
    hold on
    plot(vx,vy,'r')
    vy=mean(smLporipsth{m}(winidx,:)/bin,2);
    vy=vy';
    patch([vx fliplr(vx) vx(1)],[vy+Lporiste{m}(winidx) fliplr(vy-Lporiste{m}(winidx)) vy(1)+Lporiste{m}(1)]',[0.8 1 0.8])
    plot(vx,vy,'g')
    title(['Lp RF' num2str(90*(m>2)+45) 'Contra' num2str(90*mod(m+1,2)+45)])
    xlabel('Time (ms)')
    ylabel('FR (Hz)')
end
title(['Lp RF' num2str(90*(m>2)+45) 'Contra' num2str(90*mod(m+1,2)+45) 'N=' num2str(size(Lpcellconpsth{1},2))])

figure
for m=1:4
    subplot(2,2,m)
    vx=(win(1)+bin:bin:win(2))-0.1;
    vy=mean(smRpconpsth{m}(winidx,:)/bin,2);
    vy=vy';
    patch([vx fliplr(vx) vx(1)],[vy+Rpconste{m}(winidx) fliplr(vy-Rpconste{m}(winidx)) vy(1)+Rpconste{m}(1)],[1 0.8 0.8])
    hold on
    plot(vx,vy,'r')
    vy=mean(smRporipsth{m}(winidx,:)/bin,2);
    vy=vy';
    patch([vx fliplr(vx) vx(1)],[vy+Rporiste{m}(winidx) fliplr(vy-Rporiste{m}(winidx)) vy(1)+Rporiste{m}(1)]',[0.8 1 0.8])
    plot(vx,vy,'g')
    title(['Rp RF' num2str(90*(m>2)+45) 'Contra' num2str(90*mod(m+1,2)+45)])
    xlabel('Time (ms)')
    ylabel('FR (Hz)')
end
title(['Rp RF' num2str(90*(m>2)+45) 'Contra' num2str(90*mod(m+1,2)+45) 'N=' num2str(size(Rpcellconpsth{1},2))])
figure
for m=1:4
    subplot(2,2,m)
    vx=(win(1)+bin:bin:win(2))-0.1;
    vy=mean(smMpconpsth{m}(winidx,:)/bin,2);
    vy=vy';
    patch([vx fliplr(vx) vx(1)],[vy+Mpconste{m}(winidx) fliplr(vy-Mpconste{m}(winidx)) vy(1)+Mpconste{m}(1)],[1 0.8 0.8])
    hold on
    plot(vx,vy,'r')
    vy=mean(smMporipsth{m}(winidx,:)/bin,2);
    vy=vy';
    patch([vx fliplr(vx) vx(1)],[vy+Mporiste{m}(winidx) fliplr(vy-Mporiste{m}(winidx)) vy(1)+Mporiste{m}(1)]',[0.8 1 0.8])
    plot(vx,vy,'g')
    title(['Mp RF' num2str(90*(m>2)+45) 'Contra' num2str(90*mod(m+1,2)+45)])
    xlabel('Time (ms)')
    ylabel('FR (Hz)')
end
title(['Mp RF' num2str(90*(m>2)+45) 'Contra' num2str(90*mod(m+1,2)+45) 'N=' num2str(size(Mpcellconpsth{1},2))])