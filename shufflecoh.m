function [t_puff,t_neu,t_phip,t_phin,t_f]=shufflecoh(camypuff,cV1puff,camyneu,cV1neu,params)
% camypuff=apuff;
% camyneu=aneu;
% cV1puff=vpuff;
% cV1neu=vneu;
npairsession=30;
mixcount=size(camypuff,3);

t_puff=single([]);
t_neu=single([]);
t_phip=single([]);
t_phin=single([]);
if ~params.sldwinlabel
    gamypuff=gpuArray(camypuff);
    gamyneu=gpuArray(camyneu);
    gV1puff=gpuArray(cV1puff);
    gV1neu=gpuArray(cV1neu);

    for mshuffle=1:size(camypuff,2)-1
        sprintf(repmat('*',1,mshuffle))
        shuffleidx=[mshuffle+1:size(camypuff,2) 1:mshuffle];
        for n=1:ceil(mixcount/npairsession)
            if n*npairsession>mixcount
                pairidx=1+npairsession*(n-1):mixcount;
            else
                pairidx=1+npairsession*(n-1):npairsession*n;
            end
            [Cpuff,phip,~,~,~,fp]=coherencyc_gpu(gamypuff(:,:,pairidx),gV1puff(:,shuffleidx,pairidx),params);
            [Cneu,phin,~,~,~,fn]=coherencyc_gpu(gamyneu(:,:,pairidx),gV1neu(:,shuffleidx,pairidx),params);
            t_puff(:,pairidx,mshuffle)=gather(Cpuff);
            t_neu(:,pairidx,mshuffle)=gather(Cneu);
            t_phip(:,pairidx,mshuffle)=gather(phip);
            t_phin(:,pairidx,mshuffle)=gather(phin);
        end
    end
    t_puff=permute(t_puff,[1 3 2]);
    t_neu=permute(t_neu,[1 3 2]);
    t_phip=permute(t_phip,[1 3 2]);
    t_phin=permute(t_phin,[1 3 2]);
    t_f=gather(fp);    
else % sliding window
    sldwin=params.sldwin;
%     if sldwin(1)<=0.256
%         params.pad = 3;
%     end
    for mshuffle=1:size(camypuff,2)-1
        sprintf(repmat('*',1,mshuffle))
        shuffleidx=[mshuffle+1:size(camypuff,2) 1:mshuffle];
        [Cpuff,phip,~,~,~,t,fp]=cohgramc_gpuMKII(camypuff,cV1puff(:,shuffleidx,:),sldwin,params);
        [Cneu,phin,~,~,~,~,~]=cohgramc_gpuMKII(camyneu,cV1neu(:,shuffleidx,:),sldwin,params);
        t_puff(:,:,:,mshuffle)=gather(Cpuff);
        t_neu(:,:,:,mshuffle)=gather(Cneu);
        t_phip(:,:,:,mshuffle)=gather(phip);
        t_phin(:,:,:,mshuffle)=gather(phin);
    end
    t_puff=permute(t_puff,[1 2 4 3]);
    t_neu=permute(t_neu,[1 2 4 3]);
    t_phip=permute(t_phip,[1 2 4 3]);
    t_phin=permute(t_phin,[1 2 4 3]);
    t_f.t=gather(t);
    t_f.f=gather(fp);
end