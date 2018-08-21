function [S1merge,S2merge]=merge2blockamyspec(CONORI,sldwinlabel,cellselec)
% win=[401 800];
twin=201:800;
blwin=1:200;
stimwin=201:600;
movingwin=[0.2 0.010];
rmnoisycell=0;
rmbl=1;

% CONORI= 'con' 'roi'
% S2 neu条件做了trial数匹配，90中选了30 
% win=[901 1300];
folderlist={
    'D:\data extracted\160712';
    'D:\data extracted\160713';
    'D:\data extracted\160714';
    'D:\data extracted\160715';
    'D:\data extracted\160716';
    'D:\data extracted\160719';
    'D:\data extracted\160720';
    'D:\data extracted\160721';
    'D:\data extracted\160722';
    'D:\data extracted\160624';
    'D:\data extracted\160625';
    'D:\data extracted\160627';
    'D:\data extracted\160628';
    'D:\data extracted\170426';
    'D:\data extracted\170427';
    'D:\data extracted\170429';
    'D:\data extracted\170430';
    'D:\data extracted\170502';
    'D:\data extracted\170503';
    'D:\data extracted\170516';
    'D:\data extracted\170518';
    'D:\data extracted\170519';
    'D:\data extracted\170522';
    };
LRlist=['rrrrrrrrr' 'llll' 'rrrrrr' 'llll'];

params.Fs = 1000; % sampling frequency
params.tapers = [3 5]; % taper parameters
params.trialave = 0;
params.fpass = [4 100];
params.pad = 2;
params.err =[2 0.05]; 
params.movingwin=movingwin;
% params.lr   L==0 R==1;
params.iscon=strcmpi(CONORI,'con');  % ori=0 con=1

S1merge=[];
S2merge=[];
proportionlist=[];
load('D:\result\170913\noiseyacell.mat');
for iday=1:numel(folderlist)
    tic
    sprintf(['day ' num2str(iday)])
    S1daymerge=[];
    S2daymerge=[];
    params.lr=strcmp(LRlist(iday),'r');  % L==0 R==1;
    if strcmpi(params.lr,'r')
        cndlist=[1 2 3 4 5 6];
    elseif strcmpi(params.lr,'l')
        cndlist=[4 5 6 1 2 3];
    end
    if params.iscon
        load([folderlist{iday} '\amyconblockunitlfp.mat']);
        [~,s2,s3,s4,s5]=size(amyconblockunitlfp.lfp);
        s1=numel(twin);
        m2b=floor(s5/2);
        amyrawlfp=permute(amyconblockunitlfp.lfp(twin,:,:,:,1:2*m2b),[1 2 5 3 4]);
        amyrawlfp=reshape(amyrawlfp,s1,[],m2b,s3,s4);
        amyrawlfp=permute(amyrawlfp,[1,2,4,5,3]);
        if  cellselec
            load([folderlist{iday} '\amyconblockunitpsth.mat']);
            amyrawspk=permute(amyconblockunitpsth.psth(twin,:,:,:,1:2*m2b),[1,2,5,3,4]);
            amyrawspk=reshape(amyrawspk,s1,[],m2b,s3,s4);
            amyrawspk=permute(amyrawspk,[1,2,4,5,3]);
            [~,sa1,sa2,sa3,sa4]=size(amyrawspk);
            %
            amyblfr=squeeze(mean(amyrawspk(blwin,:,:,:,:)));
            amyblfr=reshape(amyblfr,sa1*sa2,[]);
            amystimfr=squeeze(mean(amyrawspk(stimwin,:,:,:,:)));
            amypufffr=amystimfr(:,1:3,:,:);
            amyneufr=amystimfr(:,4:end,:,:);
            amystimfr=reshape(amystimfr,sa1*sa2,[]);
            [t1h,t1p]=ttest(amyblfr,amystimfr);
            amystimvalid=reshape(t1h==1,sa3,sa4);
            
            amypufffr=reshape(amypufffr,[],sa3*sa4);
            amyneufr=reshape(amyneufr,[],sa3*sa4);
            [t2h,t2p]=ttest2(amypufffr,amyneufr);
            posneg=mean(amypufffr,1)<mean(amyneufr,1);
            valencevalid=reshape(t2h==1 & posneg,sa3,sa4);
            amyspkvalid=amyconblockunitlfp.spikevalid(:,1:2:m2b*2-1) & amyconblockunitlfp.spikevalid(:,2:2:m2b*2);
%             amycellidx=amystimvalid & amyspkvalid & valencevalid;
            amycellidx=amystimvalid & valencevalid;
            if rmnoisycell
                if exist(['D:\result\170828 mergeblock lfpspkcoh\' folderlist{iday}(end-5:end) '\con_lfp_lfp_coherencepool_g.mat'],'file')
                    load(['D:\result\170828 mergeblock lfpspkcoh\' folderlist{iday}(end-5:end) '\con_lfp_lfp_coherencepool_g.mat']);
                    if ~isempty(noisyacell{iday}) && sum(noisyacell{iday})>0
                        acellid=arrayfun(@(x) x.acellidx,coherencepool.blockcellidx(noisyacell{iday}));
                        blockid=arrayfun(@(x) x.blockidx,coherencepool.blockcellidx(noisyacell{iday}));
                        x=1:sa3;
                        y=unique(acellid);
                        [X,Y]=meshgrid(x,y);
                        noisycell=bsxfun(@eq,X,Y);
                        noisycell=sum(noisycell,1);
                        x=1:sa4;
                        y=unique(blockid);
                        [X,Y]=meshgrid(x,y);
                        noisyblock=bsxfun(@eq,X,Y);
                        noisyblock=sum(noisyblock,1);
                        [xx,yy]=meshgrid(noisyblock,noisycell);
                        zz=logical(xx.*yy);
                        amycellidx(zz)=0;
                    end
                end
            end
        else
            amycellidx=ones(24,m2b);
        end                
        amypufflfp=amyrawlfp(:,:,1:3,amycellidx);
        amyneulfp=amyrawlfp(:,:,4:6,amycellidx);
    else
        load([folderlist{iday} '\amyoriblockunitlfp.mat']);
        load([folderlist{iday} '\amyoriblockunitpsth.mat']);
        amyrawspk=amyoriblockunitpsth.psth(twin,:,:,:,:);
        
        [~,sp1,sp2,sp3,sp4]=size(amyrawdata);
        amyblfr=mean(amyrawspk(blwin,:,:,:,:),1);
        amyblfr=reshape(amyblfr,sp1,sp2,sp3,sp4);
        amyblfr=reshape(amyblfr,sp1*sp2,[]);
        amystimfr=mean(amyrawspk(stimwin,:,:,:,:),1);
        amystimfr=reshape(amystimfr,sp1,sp2,sp3,sp4);
        amystimfr=reshape(amystimfr,sp1*sp2,[]);
        [t1h,t1p]=ttest(amyblfr,amystimfr);
        amystimvalid=reshape(t1h==1,sp3,sp4);
        valencevalid=ones(24,s5);
        amyspkvalid=amyoriblockunitlfp.spikevalid;
        amycellidx=amystimvalid & amyspkvalid & valencevalid;
        
        amypufflfp=amyoriblockunitlfp.lfp(twin,:,cndlist(1:3),amycellidx);
        amyneulfp=amyoriblockunitlfp.lfp(twin,:,cndlist(4:6),amycellidx);
    end
    if ~isempty(amypufflfp)
    [~,sz2,sz3,sz4,sz5]=size(amypufflfp);
    szbl1=numel(blwin);
    temp=reshape(amypufflfp(blwin,:,:,:),szbl1,[]);
    temp=rmlinesc(temp,params,[],'n',[50 100]);
    amypufflfp_bl=reshape(temp,szbl1,[],sz4,sz5);
    temp=reshape(amyneulfp(blwin,:,:,:),szbl1,[]);
    temp=rmlinesc(temp,params,[],'n',[50 100]);
    amyneulfp_bl=reshape(temp,szbl1,[],sz4,sz5);
    szstim1=numel(stimwin);
    temp=reshape(amypufflfp(stimwin,:,:,:),szstim1,[]);
    temp=rmlinesc(temp,params,[],'n',[50 100]);
    amypufflfp_stim=reshape(temp,szstim1,[],sz4,sz5);
    temp=reshape(amyneulfp(stimwin,:,:,:),szstim1,[]);
    temp=rmlinesc(temp,params,[],'n',[50 100]);
    amyneulfp_stim=reshape(temp,szstim1,[],sz4,sz5);

    for nblock=1:sz5
        for ncell=1:sz4
            Sstim1=[];
            Sstim2=[];
            if sldwinlabel
                [Sstim1,t,f]=mtspecgramc(amypufflfp_stim(:,:,ncell,nblock),params.movingwin,params);
                [Sstim2,~,~]=mtspecgramc(amyneulfp_stim(:,:,ncell,nblock),params.movingwin,params);
                if rmbl
                    [Sbl1,~,~]=mtspecgramc(amypufflfp_bl(:,:,ncell,nblock),params.movingwin,params);
                    [Sbl2,~,~]=mtspecgramc(amyneulfp_bl(:,:,ncell,nblock),params.movingwin,params);
                    if ndims(Sbl1)==2
                        temp=Sbl1;
                        Sbl1=[];
                        Sbl1(1,:,:)=temp;
                        temp=Sbl2;
                        Sbl2=[];
                        Sbl2(1,:,:)=temp;
                    end
                    Stemp=Sstim1./repmat(Sbl1,size(Sstim1,1),1,1);
                    S1daymerge=cat(3,S1daymerge,mean(Stemp,3));
                    Stemp=Sstim2./repmat(Sbl2,size(Sstim2,1),1,1);
                    S2daymerge=cat(3,S2daymerge,mean(Stemp,3));
                else
                    S1daymerge=cat(3,S1daymerge,mean(Sstim1,3));
                    S2daymerge=cat(3,S2daymerge,mean(Sstim2,3));
                end
            else
                [Sstim1,f]=mtspectrumc(amypufflfp_bl(:,:,ncell,nblock),params);
                [Sstim2,~]=mtspectrumc(amyneulfp_bl(:,:,ncell,nblock),params);
                if rmbl
                    [Sbl1,~]=mtspectrumc(amypufflfp_bl(:,:,ncell,nblock),params);
                    [Sbl2,~]=mtspectrumc(amyneulfp_bl(:,:,ncell,nblock),params);
                    if ndims(Sbl1)==1
                        temp=Sbl1;
                        Sbl1=[];
                        Sbl1(1,:)=temp;
                        temp=Sbl2;
                        Sbl2=[];
                        Sbl2(1,:)=temp;
                    end
                    Stemp=Sstim1./repmat(Sbl1,size(Sstim1,1),1);
                    S1daymerge=cat(2,S1daymerge,mean(Stemp,2));
                    Stemp=Sstim2./repmat(Sbl2,size(Sstim2,1),1);
                    S2daymerge=cat(2,S2daymerge,mean(Stemp,2));
                else
                    S1daymerge=cat(2,S1daymerge,mean(Sstim1,2));
                    S2daymerge=cat(2,S2daymerge,mean(Sstim2,2));
                end               
            end
        end
    end
    %     figure;
    %     for ncell=1:sl5*sl4
    %         if ncell<24
    %             subplot(4,6,ncell)
    %             semilogy(f,S1daymerge(:,ncell),'r')
    %             hold on
    %             semilogy(f,S2daymerge(:,ncell),'g')
    %             title(num2str(ncell))
    %         end
    %     end
    %     title(['day ' num2str(iday)])
    if sldwinlabel
        S1merge=cat(3,S1merge,S1daymerge);
        S2merge=cat(3,S2merge,S2daymerge);
    else
        S1merge=cat(3,S1merge,S1daymerge);
        S2merge=cat(3,S2merge,S2daymerge);
    end
    end
    if sldwinlabel
        sprintf(['Day ' num2str(iday) ' Ncell= ' num2str(size(S1daymerge,3))])
    else
        sprintf(['Day ' num2str(iday) ' Ncell= ' num2str(size(S1daymerge,2))])
    end
end
if sldwinlabel    
    figure
    subplot(1,2,1)
    if rmbl
    plot_matrix(mean(S1merge,3),t,f,'n');
    else
    plot_matrix(mean(S1merge,3),t,f,'l');    
    end
    title('puff ori')
    subplot(1,2,2)
    if rmbl
    plot_matrix(mean(S2merge,3),t,f,'n');
    else
    plot_matrix(mean(S2merge,3),t,f,'l');
    end
    title('neu ori')
    
%     subplot(1,3,3)
%     plot_matrix(mean(S1merge,3)./mean(S2merge,3),t,f,'l');
else
    if rmbl
    figure
    plot(f,mean(S1merge,2),'r')
    hold on
    plot(f,mean(S2merge,2),'g')
    else
    figure
    semilogy(f,mean(S1merge,2),'r')
    hold on
    semilogy(f,mean(S2merge,2),'g')
    end
end
sprintf('done')
