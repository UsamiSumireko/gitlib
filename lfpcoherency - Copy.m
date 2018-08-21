clear
folderlist={'D:\data extracted\160720'};
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
tic
coherencepool=struct('puff',[],'neu',[]);
for iday=1:numel(folderlist)
    fileexistlabel=exist(['D:\result\170713coherency\' folderlist{iday}(end-5:end) 'amylfpdata.mat'],'file') && ...
                   exist(['D:\result\170713coherency\' folderlist{iday}(end-5:end) 'V1lfpdata.mat'],'file');
    if ~fileexistlabel           
    load([folderlist{iday} '\amyconblockunitlfp.mat']);
    load([folderlist{iday} '\V1conblockunitlfp.mat']);
    for nblock=1:size(amyconblockunitlfp.lfp,5)
        amycelidx=find(amyconblockunitlfp.spikevalid(:,nblock));
        for acell=1:numel(amycelidx)
            data1=[];
            data2=[];
            for mstim=1:3
                temp1=squeeze(amyconblockunitlfp.lfp(:,:,mstim,acell,nblock));
                temp2=reshape(amyconblockunitlfp.lfp(:,:,multineucnd(:,mstim),acell,nblock),size(amyconblockunitlfp.lfp,1),[]);
                data1=[data1 temp1];           
                data2=[data2 temp2];
            end
            amylfpdata(acell,nblock).puff=data1;
            amylfpdata(acell,nblock).neu=data2;
        end
        for vcell=1:size(V1conblockunitlfp.lfp,4)
            data1=[];
            data2=[];
            for mstim=1:3
                temp1=squeeze(V1conblockunitlfp.lfp(:,:,mstim,vcell,nblock));
                temp2=reshape(V1conblockunitlfp.lfp(:,:,multineucnd(:,mstim),vcell,nblock),size(V1conblockunitlfp.lfp,1),[]);
                data1=[data1 temp1];           
                data2=[data2 temp2];
            end
            V1lfpdata(vcell,nblock).puff=data1;
            V1lfpdata(vcell,nblock).neu=data2;
        end
    end
    save(['D:\result\170713coherency\' folderlist{iday}(end-5:end) '\amylfpdata.mat'],'amylfpdata','-v7.3');
    save(['D:\result\170713coherency\' folderlist{iday}(end-5:end) '\V1lfpdata.mat'],'V1lfpdata','-v7.3');
    else
        load(['D:\result\170713coherency\' folderlist{iday}(end-5:end) '\amylfpdata.mat'])
        load(['D:\result\170713coherency\' folderlist{iday}(end-5:end) '\V1lfpdata.mat'])
    end
    
    for nblock=1:size(amylfpdata,2)
        for acell=1:size(amylfpdata,1)
            if ~isempty(amylfpdata(acell,nblock).puff)
                for vcell=1:size(V1lfpdata,1)
                    if ~isempty(V1lfpdata(vcell,nblock).puff)
                        nanidx1=isnan(amylfpdata(acell,nblock).puff) | isnan(V1lfpdata(vcell,nblock).puff);
                        nanidx2=isnan(amylfpdata(acell,nblock).neu) | isnan(V1lfpdata(vcell,nblock).neu);
                        nidx1=sum(0==sum(nanidx1));
                        nidx2=sum(0==sum(nanidx2));
                        nannumidx2=find(0==sum(nanidx2));
                        selecidx=nannumidx2(randperm(nidx2,nidx1));
                        if gpulabel
                            gamylfppuff=gpuArray(amylfpdata(acell,nblock).puff(:,0==sum(nanidx1)));
                            gV1lfppuff=gpuArray(V1lfpdata(vcell,nblock).puff(:,0==sum(nanidx1)));
                            gamylfpneu=gpuArray(amylfpdata(acell,nblock).neu(:,selecidx));
                            gV1lfpneu=gpuArray(V1lfpdata(vcell,nblock).neu(:,selecidx));
                            
                            [Cpuff,phip,~,~,~,fp]=coherencyc(gamylfppuff,gV1lfppuff,params);
                            [Cneu,phin,~,~,~,fn]=coherencyc(gamylfpneu,gV1lfpneu,params);

                            t_puff(:,vcell,acell,nblock)=gather(Cpuff);
                            t_neu(:,vcell,acell,nblock)=gather(Cneu);
                            t_f(:,vcell,acell,nblock)=gather(fp);
                            t_neuidx(:,vcell,acell,nblock)=gather(selecidx);
                            t_puffidx(:,vcell,acell,nblock)=gather(find(0==sum(nanidx1)));
%                             coherencepool.puff(:,vcell,acell,nblock)=gather(Cpuff);
%                             coherencepool.neu(:,vcell,acell,nblock)=gather(Cneu);
%                             coherencepool.f(:,vcell,acell,nblock)=gather(fp);
%                             coherencepool.neuidx(:,vcell,acell,nblock)=gather(selecidx);
%                             coherencepool.puffidx(:,vcell,acell,nblock)=gather(find(0==sum(nanidx1)));
                        else
                            namylfppuff=amylfpdata(acell,nblock).puff(:,0==sum(nanidx1));
                            nV1lfppuff=V1lfpdata(vcell,nblock).puff(:,0==sum(nanidx1));
                            namylfpneu=amylfpdata(acell,nblock).neu(:,selecidx);
                            nV1lfpneu=V1lfpdata(vcell,nblock).neu(:,selecidx);
                            
                            [Cpuff,phip,~,~,~,fp]=coherencyc(namylfppuff,nV1lfppuff,params);
                            [Cneu,phin,~,~,~,fn]=coherencyc(namylfpneu,nV1lfpneu,params);
                            
                            t_puff(:,vcell,acell,nblock)=gather(Cpuff);
                            t_neu(:,vcell,acell,nblock)=gather(Cneu);
                            t_f(:,vcell,acell,nblock)=gather(fp);
                            t_neuidx(:,vcell,acell,nblock)=gather(selecidx);
                            t_puffidx(:,vcell,acell,nblock)=gather(find(0==sum(nanidx1)));
%                             coherencepool.puff(:,vcell,acell,nblock)=Cpuff;
%                             coherencepool.neu(:,vcell,acell,nblock)=Cneu;
%                             coherencepool.f(:,vcell,acell,nblock)=fp;
%                             coherencepool.neuidx(:,vcell,acell,nblock)=selecidx;
%                             coherencepool.puffidx(:,vcell,acell,nblock)=find(0==sum(nanidx1));
                        end
                    end
                end
            end
        end
    end
    coherencepool.puff=t_puff;
    coherencepool.neu=t_neu;
    coherencepool.f=t_f;
    coherencepool.neuidx=t_neuidx;
    coherencepool.puffidx=t_puffidx;
    if gpulabel
        save(['D:\result\170713coherency\' folderlist{iday}(end-5:end) '\coherencepool_g.mat'],'coherencepool','-v7.3');
    else
        save(['D:\result\170713coherency\' folderlist{iday}(end-5:end) '\coherencepool.mat'],'coherencepool','-v7.3');
    end
    
end
toc