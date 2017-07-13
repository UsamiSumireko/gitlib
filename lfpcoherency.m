clear
folderlist={'D:\data extracted\160720'};
str='D:\result\170713coherency\';
gpulabel=true;
multineucnd=[4, 5, 6;
             7, 8, 9
             10, 11, 12];
         
params.Fs = 1000; % sampling frequency
params.tapers = [2 3]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 0;
params.err =[2 0.05];
for iday=1:numel(folderlist)
    coherencepool=struct('puff',[],'neu',[]);
    fileexistlabel=exist([str folderlist{iday}(end-5:end) 'amylfpdata.mat'],'file') && ...
                   exist([str folderlist{iday}(end-5:end) 'V1lfpdata.mat'],'file');
    if ~fileexistlabel           
    load([folderlist{iday} '\amyconblockunitlfp.mat']);
    load([folderlist{iday} '\V1conblockunitlfp.mat']);
    acount=0;
    vcount=0;
    amylfpidx=struct();
    V1lfpidx=struct();
    blocknum=size(amyconblockunitlfp.lfp,5);
    for nblock=1:blocknum
        amycellidx=find(amyconblockunitlfp.spikevalid(:,nblock));
        for acell=1:numel(amycellidx)
            acount=acount+1;
            data1=[];
            data2=[];
            for mstim=1:3
                temp1=squeeze(amyconblockunitlfp.lfp(:,:,mstim,acell,nblock));
                temp2=reshape(amyconblockunitlfp.lfp(:,:,multineucnd(:,mstim),acell,nblock),size(amyconblockunitlfp.lfp,1),[]);
                data1=[data1 temp1];           
                data2=[data2 temp2];
            end
            amylfpdata(acount).puff=data1;
            amylfpdata(acount).neu=data2;
            amylfpidx.block(acount)=nblock;
            amylfpidx.cell(acount)=amycellidx(acell);
        end
        for vcell=1:size(V1conblockunitlfp.lfp,4)
            vcount=vcount+1;
            data1=[];
            data2=[];
            for mstim=1:3
                temp1=squeeze(V1conblockunitlfp.lfp(:,:,mstim,vcell,nblock));
                temp2=reshape(V1conblockunitlfp.lfp(:,:,multineucnd(:,mstim),vcell,nblock),size(V1conblockunitlfp.lfp,1),[]);
                data1=[data1 temp1];           
                data2=[data2 temp2];
            end
            V1lfpdata(vcount).puff=data1;
            V1lfpdata(vcount).neu=data2;
            V1lfpidx.block(vcount)=nblock;
            V1lfpidx.cell(vcount)=vcell;
        end
    end
    save([str folderlist{iday}(end-5:end) '\amylfpdata.mat'],'amylfpdata','-v7.3');
    save([str folderlist{iday}(end-5:end) '\amylfpidx.mat'],'amylfpidx','-v7.3');
    save([str folderlist{iday}(end-5:end) '\V1lfpdata.mat'],'V1lfpdata','-v7.3');
    save([str folderlist{iday}(end-5:end) '\V1lfpidx.mat'],'V1lfpidx','-v7.3');
    else
        load([str folderlist{iday}(end-5:end) '\amylfpdata.mat'])
        load([str folderlist{iday}(end-5:end) '\amylfpidx.mat'])
        load([str folderlist{iday}(end-5:end) '\V1lfpdata.mat'])
        load([str folderlist{iday}(end-5:end) '\V1lfpidx.mat'])
    end
    tic
    amylfppuff=[];
    V1lfppuff=[];
    amylfpneu=[];
    V1lfpneu=[];
    blockcellidx=[];
    trialidx=[];
    mixcount=0;
    for nblock=1:blocknum
        ablockidx=find(amylfpidx.block==nblock);
        vblockidx=find(V1lfpidx.block==nblock);
        for acell=1:numel(ablockidx)
            for vcell=1:numel(vblockidx)
                mixcount=mixcount+1;
%                 mixcount=vcell+(acell-1)*numel(vblockidx)+(nblock-1)*numel(ablockidx)*numel(vblockidx);
                nanidx1=isnan(amylfpdata(ablockidx(acell)).puff) | isnan(V1lfpdata(vblockidx(vcell)).puff);
                nanidx2=isnan(amylfpdata(ablockidx(acell)).neu) | isnan(V1lfpdata(vblockidx(vcell)).neu);
                nidx1=sum(0==sum(nanidx1));
                assert(nidx1==30);
                nidx2=sum(0==sum(nanidx2));
                nannumidx2=find(0==sum(nanidx2));
                selecidx=nannumidx2(randperm(nidx2,nidx1));
%                 t_amylfppuff(1:1200,1:30,vcell,acell)=amylfpdata(ablockidx(acell)).puff(:,0==sum(nanidx1));
%                 t_V1lfppuff(1:1200,1:30,vcell,acell)=V1lfpdata(vblockidx(vcell)).puff(:,0==sum(nanidx1));
%                 t_amylfpneu(1:1200,1:30,vcell,acell)=amylfpdata(ablockidx(acell)).neu(:,selecidx);
%                 t_V1lfpneu(1:1200,1:30,vcell,acell)=V1lfpdata(vblockidx(vcell)).neu(:,selecidx);
%                 t_blockcellidx(vcell,acell).blockidx=nblock;
%                 t_blockcellidx(vcell,acell).acellidx=acell;
%                 t_blockcellidx(vcell,acell).vcellidx=vcell;
%                 t_trialidx(vcell,acell).idx1=(0==sum(nanidx1));
%                 t_trialidx(vcell,acell).idx2=selecidx;
                amylfppuff(:,:,mixcount)=amylfpdata(ablockidx(acell)).puff(:,0==sum(nanidx1));
                V1lfppuff(:,:,mixcount)=V1lfpdata(vblockidx(vcell)).puff(:,0==sum(nanidx1));
                amylfpneu(:,:,mixcount)=amylfpdata(ablockidx(acell)).neu(:,selecidx);
                V1lfpneu(:,:,mixcount)=V1lfpdata(vblockidx(vcell)).neu(:,selecidx);
                blockcellidx(mixcount).blockidx=nblock;
                blockcellidx(mixcount).acellidx=acell;
                blockcellidx(mixcount).vcellidx=vcell;
                trialidx(mixcount).idx1=(0==sum(nanidx1));
                trialidx(mixcount).idx2=selecidx;
            end
        end
%         amylfppuff=cat(3,amylfppuff,reshape(t_amylfppuff,size(t_amylfppuff,1),size(t_amylfppuff,2),[]));
%         V1lfppuff=cat(3,V1lfppuff,reshape(t_V1lfppuff,size(t_V1lfppuff,1),size(t_V1lfppuff,2),[]));
%         amylfpneu=cat(3,amylfpneu,reshape(t_amylfpneu,size(t_amylfpneu,1),size(t_amylfpneu,2),[]));
%         V1lfpneu=cat(3,V1lfpneu,reshape(t_V1lfpneu,size(t_V1lfpneu,1),size(t_V1lfpneu,2),[]));
%         blockcellidx=cat(1,blockcellidx,t_blockcellidx(:));
%         trialidx=cat(1,trialidx,t_trialidx(:));
    end
    toc 
    tic
    gamylfppuff=gpuArray(amylfppuff);
    gV1lfppuff=gpuArray(V1lfppuff);
    gamylfpneu=gpuArray(amylfpneu);
    gV1lfpneu=gpuArray(V1lfpneu); 
    
    [Cpuff,phip,~,~,~,fp]=coherencyc(gamylfppuff(:,:,1),gV1lfppuff(:,:,1),params);
    t_puff=nan(numel(Cpuff),mixcount);
    t_neu=nan(numel(Cpuff),mixcount);
    t_f=nan(numel(fp),mixcount);
    t_neuidx=nan(30,mixcount);
    t_puffidx=nan(30,mixcount);
    parfor n=1:mixcount
        [Cpuff,phip,~,~,~,fp]=coherencyc(gamylfppuff(:,:,n),gV1lfppuff(:,:,n),params);
        [Cneu,phin,~,~,~,fn]=coherencyc(gamylfpneu(:,:,n),gV1lfpneu(:,:,n),params);
        t_puff(:,n)=gather(Cpuff);
        t_neu(:,n)=gather(Cneu);
        t_f(:,n)=gather(fp);
        t_neuidx(:,n)=gather(selecidx);
        t_puffidx(:,n)=gather(find(0==sum(nanidx1)));
    end
    coherencepool.puff=t_puff;
    coherencepool.neu=t_neu;
    coherencepool.f=t_f;
    coherencepool.neuidx=t_neuidx;
    coherencepool.puffidx=t_puffidx;
    coherencepool.blockcellidx=blockcellidx;
    coherencepool.trialidx=trialidx;
    if gpulabel
        save([str folderlist{iday}(end-5:end) '\coherencepool_g.mat'],'coherencepool','-v7.3');
    else
        save([str folderlist{iday}(end-5:end) '\coherencepool.mat'],'coherencepool','-v7.3');
    end
    
end
toc