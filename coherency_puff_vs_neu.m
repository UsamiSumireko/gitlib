function coherency_puff_vs_neu(folderlist,str,twin,CONORI,LfpSpk,rmlinenoiselabel,sldwinlabel,trialshufflelabel)
% folderlist={'D:\data extracted\160720'};
% str='D:\result\170713coherency\';
% twin=[0.1 0.5]-(-0.1);
gpulabel=true;
% CONORI='con';
% LfpSpk=[{'lfp'} {'lfp'}]; % LfpSpk(1) amygdala  LfpSpk(2) V1   'lfp' or 'psth'
% rmlinenoiselabel=true;
% sldwinlabel=true;
% trialshufflelabel=true;
resamplelabel=0;
blockmergelabel=1;
%%
assert(strcmp(CONORI,'con') || strcmp(CONORI,'ori'));
assert((strcmp(LfpSpk{1},'lfp') || strcmp(LfpSpk{1},'psth')) && (strcmp(LfpSpk{2},'lfp') || strcmp(LfpSpk{2},'psth'))); 
multineucnd=[4, 5, 6;
             7, 8, 9;
             10, 11, 12];
params.Fs = 1000; % sampling frequency
params.tapers = [2 3]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 0;
params.err =[2 0.05];
bwin=(twin(1)*params.Fs+1):twin(2)*params.Fs;
npairsession=200;

for iday=1:numel(folderlist)
    iday,
    coherencepool=struct('puff',[],'neu',[]);
    fileexistlabel=exist([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} 'data.mat'],'file') && ...
                   exist([str folderlist{iday}(end-5:end) '\V1' CONORI LfpSpk{2} 'data.mat'],'file');
    if ~fileexistlabel           
    load([folderlist{iday} '\amy' CONORI 'blockunit' LfpSpk{1} '.mat']);
    load([folderlist{iday} '\V1' CONORI 'blockunit' LfpSpk{2} '.mat']);
    amyblockunit.data=eval(['amy' CONORI 'blockunit' LfpSpk{1} '.' LfpSpk{1}]);
    amyblockunit.spikevalid=eval(['amy' CONORI 'blockunit' LfpSpk{1} '.spikevalid']);
    V1blockunit.data=eval(['V1' CONORI 'blockunit' LfpSpk{2} '.' LfpSpk{2}]);
%     V1blockunit.data=eval(['V1' CONORI 'blockunit' LfpSpk{2} '.' LfpSpk{2}]);
    acount=0;
    vcount=0;
    amylfpidx=struct();
    V1lfpidx=struct();
    blocknum=size(amyblockunit.data,5);
    for nblock=1:blocknum
        amycellidx=find(amyblockunit.spikevalid(:,nblock));
        for acell=1:numel(amycellidx)
            acount=acount+1;
            data1=[];
            data2=[];
            for mstim=1:3
                temp1=squeeze(amyblockunit.data(:,:,mstim,acell,nblock));
                temp2=reshape(amyblockunit.data(:,:,multineucnd(:,mstim),acell,nblock),size(amyblockunit.data,1),[]);
                data1=[data1 temp1];           
                data2=[data2 temp2];
            end
            amydata(acount).puff=data1;
            amydata(acount).neu=data2;
            amyidx.block(acount)=nblock;
            amyidx.cell(acount)=amycellidx(acell);
        end
        for vcell=1:size(V1blockunit.data,4)
            vcount=vcount+1;
            data1=[];
            data2=[];
            for mstim=1:3
                temp1=squeeze(V1blockunit.data(:,:,mstim,vcell,nblock));
                temp2=reshape(V1blockunit.data(:,:,multineucnd(:,mstim),vcell,nblock),size(V1blockunit.data,1),[]);
                data1=[data1 temp1];           
                data2=[data2 temp2];
            end
            V1data(vcount).puff=data1;
            V1data(vcount).neu=data2;
            V1idx.block(vcount)=nblock;
            V1idx.cell(vcount)=vcell;
        end
    end
%     eval(['amy' CONORI LfpSpk{1} 'data=amydata;']);
%     eval(['amy' CONORI LfpSpk{1} 'idx=amyidx;']);
%     eval(['V1' CONORI LfpSpk{2} 'data=V1data;']);
%     eval(['V1' CONORI LfpSpk{2} 'idx=V1idx;']);
    if ~exist([str folderlist{iday}(end-5:end)],'dir')
        mkdir([str folderlist{iday}(end-5:end)])
    end
    save([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} 'data.mat'],'amydata','-v7.3');
    save([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} 'idx.mat'],'amyidx','-v7.3');
    save([str folderlist{iday}(end-5:end) '\V1' CONORI LfpSpk{2} 'data.mat'],'V1data','-v7.3');
    save([str folderlist{iday}(end-5:end) '\V1' CONORI LfpSpk{2} 'idx.mat'],'V1idx','-v7.3');
    else
        load([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} 'data.mat'])
        load([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} 'idx.mat'])
        load([str folderlist{iday}(end-5:end) '\V1' CONORI LfpSpk{2} 'data.mat'])
        load([str folderlist{iday}(end-5:end) '\V1' CONORI LfpSpk{2} 'idx.mat'])
    end
    tic
    mixcount=0;
    for nblock=1:numel(unique(amyidx.block))
        ablockidx=find(amyidx.block==nblock);
        vblockidx=find(V1idx.block==nblock);
        for acell=1:numel(ablockidx)
            for vcell=1:numel(vblockidx)
                mixcount=mixcount+1;
                nanidx1=isnan(amydata(ablockidx(acell)).puff) | isnan(V1data(vblockidx(vcell)).puff);
                nanidx2=isnan(amydata(ablockidx(acell)).neu) | isnan(V1data(vblockidx(vcell)).neu);
                nidx1=sum(0==sum(nanidx1));
                assert(nidx1==30);
                nidx2=sum(0==sum(nanidx2));
                nannumidx2=find(0==sum(nanidx2));
                selecidx=nannumidx2(randperm(nidx2,nidx1));                

                amypuff(:,:,mixcount)=amydata(ablockidx(acell)).puff(bwin,0==sum(nanidx1));
                V1puff(:,:,mixcount)=V1data(vblockidx(vcell)).puff(bwin,0==sum(nanidx1));
                amyneu(:,:,mixcount)=amydata(ablockidx(acell)).neu(bwin,selecidx);
                V1neu(:,:,mixcount)=V1data(vblockidx(vcell)).neu(bwin,selecidx);
                blockcellidx(mixcount).blockidx=nblock;
                blockcellidx(mixcount).acellidx=acell;
                blockcellidx(mixcount).vcellidx=vcell;
                trialidx(mixcount).idx1=(0==sum(nanidx1));
                trialidx(mixcount).idx2=selecidx;
            end
        end
    end
    toc 
    tic
    if strcmp(LfpSpk{1},'lfp') && rmlinenoiselabel
        if ~exist([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} '_rmln_puff.mat'],'file') ||...
                ~exist([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} '_rmln_neu.mat'],'file') || resamplelabel
            tamypuff=reshape(amypuff,size(amypuff,1),[]);
            tamyneu=reshape(amyneu,size(amyneu,1),[]);
            camypuff=reshape(rmlinesc(tamypuff,params,[],'n',50),size(amypuff));
            toc
            camyneu=reshape(rmlinesc(tamyneu,params,[],'n',50),size(amyneu));
            toc
            clear tamypuff tamyneu
            save([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} '_rmln_puff.mat'],'camypuff','-v7.3');
            save([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} '_rmln_neu.mat'],'camyneu','-v7.3');
        else
            load([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} '_rmln_puff.mat']);
            load([str folderlist{iday}(end-5:end) '\amy' CONORI LfpSpk{1} '_rmln_neu.mat']);
        end
    else
        camypuff=amypuff;
        camyneu=amyneu;
    end
    if strcmp(LfpSpk{1},'lfp') && rmlinenoiselabel
        if ~exist([str folderlist{iday}(end-5:end) '\v1' CONORI LfpSpk{1} '_rmln_puff.mat'],'file') ||...
                ~exist([str folderlist{iday}(end-5:end) '\v1' CONORI LfpSpk{1} '_rmln_neu.mat'],'file') || resamplelabel
            g_tV1puff=reshape(V1puff,size(V1puff,1),[]);
            g_tV1neu=reshape(V1neu,size(V1neu,1),[]);
            cV1puff=reshape(rmlinesc(g_tV1puff,params,[],'n',50),size(V1puff));
            toc
            cV1neu=reshape(rmlinesc(g_tV1neu,params,[],'n',50),size(V1neu));
            toc
            clear g_tV1puff g_tV1neu
            save([str folderlist{iday}(end-5:end) '\v1' CONORI LfpSpk{1} '_rmln_puff.mat'],'cV1puff','-v7.3');
            save([str folderlist{iday}(end-5:end) '\v1' CONORI LfpSpk{1} '_rmln_neu.mat'],'cV1neu','-v7.3');
        else
            load([str folderlist{iday}(end-5:end) '\v1' CONORI LfpSpk{1} '_rmln_puff.mat']);
            load([str folderlist{iday}(end-5:end) '\v1' CONORI LfpSpk{1} '_rmln_neu.mat']);
        end
    else
        cV1puff=V1puff;
        cV1neu=V1neu;
    end
    
    savestr='_';
    if trialshufflelabel
        savestr=[savestr 'trialshuffle_'];
        if sldwinlabel
            t_puff=[];
            t_neu=[];
            for mshuffle=1:size(camypuff,2)-1
                sprintf(repmat('*',1,mshuffle))
                shuffleidx=[mshuffle+1:size(camypuff,2) 1:mshuffle];
%                 camypuff=gather(camypuff);
%                 cV1puff=gather(cV1puff);
                [Cpuff,phip,~,~,~,t,fp]=cohgramc_gpuMKII(camypuff,cV1puff(:,shuffleidx,:),[0.2 0.01],params);
%                 camyneu=gather(gamyneu);
%                 cV1neu=gather(cV1neu);
                [Cneu,phin,~,~,~,~,fn]=cohgramc_gpuMKII(camyneu,cV1neu(:,shuffleidx,:),[0.2 0.01],params);
                t_puff(:,:,:,mshuffle)=gather(Cpuff);
                t_neu(:,:,:,mshuffle)=gather(Cneu);
            end
            t_puff=permute(t_puff,[1 2 4 3]);
            t_neu=permute(t_neu,[1 2 4 3]);
            t_t=gather(t);
            t_f=gather(fp);
%             t_neuidx(:,n)=trialidx(n).idx2;
%             t_puffidx(:,n)=trialidx(n).idx1;
        else
            gamypuff=gpuArray(camypuff);
            gamyneu=gpuArray(camyneu);
            gV1puff=gpuArray(cV1puff);
            gV1neu=gpuArray(cV1neu);
            t_puff=[];
            t_neu=[];
            for mshuffle=1:size(camypuff,2)-1
                shuffleidx=[mshuffle+1:size(camypuff,2) 1:mshuffle];
%                 for n=1:mixcount
                for n=1:ceil(mixcount/npairsession)
                    if n*npairsession>mixcount
                        pairidx=1+npairsession*(n-1):mixcount;
                    else
                        pairidx=1+npairsession*(n-1):npairsession*n;
                    end
                    [Cpuff,phip,~,~,~,fp]=coherencyc_gpu(gamypuff(:,:,pairidx),gV1puff(:,shuffleidx,pairidx),params);
                    [Cneu,phin,~,~,~,fn]=coherencyc_gpu(gamyneu(:,:,pairidx),gV1neu(:,shuffleidx,pairidx),params);
                    
                    t_puff(:,mshuffle,pairidx)=gather(Cpuff);
                    t_neu(:,mshuffle,pairidx)=gather(Cneu);
%                     t_puffphi=gather(phip);
%                     t_neuphi=gather(phin);
                end
                t_f=gather(fp);
%                 t_neuidx(:,pairidx)=trialidx(pairidx).idx2;
%                 t_puffidx(:,pairidx)=trialidx(pairidx).idx1;
            end
        end
    else
        if sldwinlabel
%             cV1puff=gather(cV1puff);
            [Cpuff,phip,~,~,~,t,fp]=cohgramc_gpuMKII(camypuff,cV1puff,[0.2 0.01],params);
%             camyneu=gather(gamyneu);
%             cV1neu=gather(cV1neu);
            [Cneu,phin,~,~,~,~,fn]=cohgramc_gpuMKII(camyneu,cV1neu,[0.2 0.01],params);
            t_puff=gather(Cpuff);
%             t_puffphi=gather(phip);
            t_neu=gather(Cneu);
%             t_neuphi=gather(phin);
            t_t=gather(t);
            t_f=gather(fp);
%             t_neuidx(:,n)=trialidx(n).idx2;
%             t_puffidx(:,n)=trialidx(n).idx1;
        else
            gamypuff=gpuArray(camypuff);
            gamyneu=gpuArray(camyneu);
            gV1puff=gpuArray(cV1puff);
            gV1neu=gpuArray(cV1neu);
            t_puff=[];
            t_neu=[];
            for n=1:ceil(mixcount/npairsession)
                if n*npairsession>mixcount
                    pairidx=1+npairsession*(n-1):mixcount;
                else
                    pairidx=1+npairsession*(n-1):npairsession*n;
                end
                [Cpuff,phip,~,~,~,fp]=coherencyc_gpu(gamypuff(:,:,pairidx),gV1puff(:,:,pairidx),params);
                [Cneu,phin,~,~,~,fn]=coherencyc_gpu(gamyneu(:,:,pairidx),gV1neu(:,:,pairidx),params);
                t_puff(:,pairidx)=gather(Cpuff);
                t_neu(:,pairidx)=gather(Cneu);
%                 t_puffphi=gather(phip);
%                 t_neuphi=gather(phin);
                t_f=gather(fp);
%                 t_neuidx(:,pairidx)=trialidx(pairidx).idx2;
%                 t_puffidx(:,pairidx)=trialidx(pairidx).idx1;
            end
        end
    end
    coherencepool.puff=t_puff;
    coherencepool.neu=t_neu;
%     coherencepool.puffphi=t_puffphi;
%     coherencepool.neuphi=t_neuphi;
    coherencepool.f=t_f;
%     coherencepool.neuidx=t_neuidx;
%     coherencepool.puffidx=t_puffidx;
    coherencepool.blockcellidx=blockcellidx;
    coherencepool.trialidx=trialidx;
    if sldwinlabel
        savestr=[savestr 'sldwin_'];
        coherencepool.t=t_t;
    end
    if ~exist([str folderlist{iday}(end-5:end)],'dir')
        mkdir([str folderlist{iday}(end-5:end)])
    end
    if gpulabel
        save([str folderlist{iday}(end-5:end) '\' savestr 'coherencepool_g.mat'],'coherencepool','-v7.3');
    else
        save([str folderlist{iday}(end-5:end) '\coherencepool.mat'],'coherencepool','-v7.3');
    end
    toc
    clear gV1puff gV1neu gamypuff gamyneu
end