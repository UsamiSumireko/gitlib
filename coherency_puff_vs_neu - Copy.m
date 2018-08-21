clear
folderlist={'D:\data extracted\160720'};
str='D:\result\170713coherency\';
twin=[0.1 0.5]-0.1;
gpulabel=true;
CONORI='con';
LfpSpk=[{'lfp'} {'lfp'}]; % LfpSpk(1) amygdala  LfpSpk(2) V1   'lfp' or 'psth'
rmlinenoiselabel=true;
sldwinlabel=true;
trialshufflelabel=true;

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
for iday=1:numel(folderlist)
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
        g_tamypuff=reshape(amypuff,size(amypuff,1),[]);
        g_tamyneu=reshape(amyneu,size(amyneu,1),[]);
        gamypuff=gpuArray(reshape(rmlinesc(g_tamypuff,params,[],'n',50),size(amypuff)));
        toc
        gamyneu=gpuArray(reshape(rmlinesc(g_tamyneu,params,[],'n',50),size(amyneu)));
        toc
        clear g_tamypuff g_tamyneu
    else
        gamypuff=gpuArray(amypuff);
        gamyneu=gpuArray(amyneu);
    end
    if strcmp(LfpSpk{1},'lfp') && rmlinenoiselabel
        g_tV1puff=reshape(V1puff,size(V1puff,1),[]);
        g_tV1neu=reshape(V1neu,size(V1neu,1),[]);
        gV1puff=gpuArray(reshape(rmlinesc(g_tV1puff,params,[],'n',50),size(V1puff)));
        toc
        gV1neu=gpuArray(reshape(rmlinesc(g_tV1neu,params,[],'n',50),size(V1neu)));
        toc
        clear g_tV1puff g_tV1neu
    else
        gV1puff=gpuArray(V1puff);
        gV1neu=gpuArray(V1neu);
    end
    savestr='_';
    if trialshufflelabel
        savestr=[savestr 'trialshuffle_'];
        for n=1:mixcount
            if sldwinlabel
                for mshuffle=1:size(gamypuff,2)-1
                    shuffleidx=[mshuffle+1:size(gamypuff,2) 1:mshuffle];
                    [Cpuff,phip,~,~,~,t,fp]=cohgramc_gpu(gamypuff(:,:,n),gV1puff(:,shuffleidx,n),[0.2 0.01],params);
                    [Cneu,phin,~,~,~,~,fn]=cohgramc_gpu(gamyneu(:,:,n),gV1neu(:,shuffleidx,n),[0.2 0.01],params);
                    t_puff(:,:,mshuffle,n)=gather(Cpuff);
                    t_neu(:,:,mshuffle,n)=gather(Cneu);
                end
                t_t(:,n)=gather(t);
                t_f(:,n)=gather(fp);
                t_neuidx(:,n)=trialidx(n).idx2;
                t_puffidx(:,n)=trialidx(n).idx1;
            else
                parfor mshuffle=1:size(gamypuff,2)-1
                    shuffleidx=[mshuffle+1:size(gamypuff,2) 1:mshuffle];
                    [Cpuff,phip,~,~,~,fp]=coherencyc(gamypuff(:,:,n),gV1puff(:,shuffleidx,n),params);
                    [Cneu,phin,~,~,~,fn]=coherencyc(gamyneu(:,:,n),gV1neu(:,shuffleidx,n),params);
                    t_puff(:,mshuffle,n)=gather(Cpuff);
                    t_neu(:,mshuffle,n)=gather(Cneu);
                end
                t_f(:,n)=gather(fp);
                t_neuidx(:,n)=trialidx(n).idx2;
                t_puffidx(:,n)=trialidx(n).idx1;
            end
        end
    else
        for n=1:mixcount
            if sldwinlabel
                [Cpuff,phip,~,~,~,t,fp]=cohgramc_gpu(gamypuff(:,:,n),gV1puff(:,:,n),[0.2 0.01],params);
                [Cneu,phin,~,~,~,~,fn]=cohgramc_gpu(gamyneu(:,:,n),gV1neu(:,:,n),[0.2 0.01],params);
                t_puff(:,:,n)=gather(Cpuff);
                t_neu(:,:,n)=gather(Cneu);
                t_t(:,n)=gather(t);
                t_f(:,n)=gather(fp);
                t_neuidx(:,n)=trialidx(n).idx2;
                t_puffidx(:,n)=trialidx(n).idx1;
            else
                [Cpuff,phip,~,~,~,fp]=coherencyc(gamypuff(:,:,n),gV1puff(:,:,n),params);
                [Cneu,phin,~,~,~,fn]=coherencyc(gamyneu(:,:,n),gV1neu(:,:,n),params);
                t_puff(:,n)=gather(Cpuff);
                t_neu(:,n)=gather(Cneu);
                t_f(:,n)=gather(fp);
                t_neuidx(:,n)=trialidx(n).idx2;
                t_puffidx(:,n)=trialidx(n).idx1;
            end
        end
    end
    coherencepool.puff=t_puff;
    coherencepool.neu=t_neu;
    coherencepool.f=t_f;
    coherencepool.neuidx=t_neuidx;
    coherencepool.puffidx=t_puffidx;
    coherencepool.blockcellidx=blockcellidx;
    coherencepool.trialidx=trialidx;
    if sldwinlabel
        savestr=[savestr 'sldwin_'];
        coherencepool.t=t_t;
    end
    if gpulabel
        save([str folderlist{iday}(end-5:end) '\' savestr 'coherencepool_g.mat'],'coherencepool','-v7.3');
    else
        save([str folderlist{iday}(end-5:end) '\coherencepool.mat'],'coherencepool','-v7.3');
    end
end
toc