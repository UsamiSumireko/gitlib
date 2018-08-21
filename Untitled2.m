clear
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
%     'D:\data extracted\160624';
%     'D:\data extracted\160625';
%     'D:\data extracted\160627';
%     'D:\data extracted\160628';
    'D:\data extracted\170426';
    'D:\data extracted\170427';
    'D:\data extracted\170429';
    'D:\data extracted\170430';
    'D:\data extracted\170502';
    'D:\data extracted\170503';
%     'D:\data extracted\170516';
%     'D:\data extracted\170518';
%     'D:\data extracted\170519';
%     'D:\data extracted\170522';
    };
count=0;
resptidx=201:600;
nbs=1000;
for nday=1:numel(folderlist)
    nday
%     if ~exist([folderlist{nday} '\n2b_valencevalid.mat'],'file')
        load([folderlist{nday} '\amyconblockunitpsth.mat'])
        [s1,s2,s3,s4,s5]=size(amyconblockunitpsth.psth);
        n2b=floor(s5/2);
        bsdiff=nan(s4,n2b,nbs);
        
        conpsthMKII=permute(amyconblockunitpsth.psth(:,:,:,:,1:2*n2b),[1 2 5 3 4]);
        conpsthMKII=reshape(conpsthMKII,s1,[],n2b,s3,s4);
        conpsthMKII=permute(conpsthMKII,[1 2 4 5 3]);
        mergecs=reshape(conpsthMKII(:,:,1:3,:,:),s1,[],s4,n2b);
        mergecs=squeeze(mean(mergecs(resptidx,:,:,:),1));
        mergeneu=reshape(conpsthMKII(:,:,4:12,:,:),s1,[],s4,n2b);
        mergeneu=squeeze(mean(mergeneu(resptidx,:,:,:),1));
        pool=cat(1,mergecs,mergeneu);
        nalltrial=size(pool,1);
        ncstrial=size(mergecs,1);
        nneutrial=size(mergeneu,1);
        realdiff=squeeze(mean(mergecs,1))-squeeze(mean(mergeneu,1));
        for i=1:nbs
            randidx=randperm(nalltrial);
            sample1=pool(randidx(1:ncstrial),:,:);
            sample2=pool(randidx(ncstrial+1:end),:,:);
            bsdiff(:,:,i)=squeeze(mean(sample1,1))-squeeze(mean(sample2,1));
        end
        diffstd=std(bsdiff,0,3);%/sqrt(ncstrial);
        valencevalid=abs(realdiff)>1.96*diffstd;
        valencevalid=reshape(valencevalid,s4,n2b);
        count=count+sum(sum(valencevalid));
        save([folderlist{nday} '\n2b_Svalencevalid.mat'],'valencevalid','-v7.3')
%     else
%         load([folderlist{nday} '\n2b_valencevalid.mat']);
%     end
end
count

% coherency_puff_vs_neu(folderlist,'D:\result\170719coherency\',[0.2 0.6],'con',[{'lfp'} {'lfp'}],1,0,0)
% coherency_puff_vs_neu(folderlist,'D:\result\170719coherency\',[0 0.6],'con',[{'lfp'} {'lfp'}],1,1,0)

% folderlist=dir('D:\result\170719coherency');
% puff=[];
% neu=[];
% for i=3:numel(folderlist)
%     load(['D:\result\170719coherency\' folderlist(i).name '\_coherencepool_g.mat'])
%     puff=cat(2,puff,coherencepool.puff);
%     neu=cat(2,neu,coherencepool.neu);
% end
% % t=coherencepool.t;
% f=coherencepool.f;
% p=mean(puff,2);
% n=mean(neu,2);
% figure
% plot(f,p,'r')
% hold on
% plot(f,n,'g')
% figure
% plot_matrix(mean(coherencepool.puff,3)./mean(coherencepool.neu,3),t,f)


% gamydata=gpuArray(camypuff);
% gV1data=gpuArray(cV1puff);
% npair=size(cV1puff,3);
% tic
% for n=1:ceil(npair/200)
%     if n*200>npair
%         pairidx=1+200*(n-1):npair;
%     else
%         pairidx=1+200*(n-1):200*n;
%     end
%     [C,phi,S12,S1,S2,f]=coherencyc_gpu(gamydata(:,:,pairidx),gV1data(:,:,pairidx),params);
%     t_CI(:,pairidx)=C;
% end
% toc
% clearvars -except gamydata gV1data params camypuff cV1puff npair
% tic
% for n=1:npair
%     [Cpuff,phip,~,~,~,fp]=coherencyc(gamydata(:,:,n),gV1data(:,:,n),params);
%     t_CII(:,n)=Cpuff;
% end
% toc
% C1=gather(C);
% C2=gather(t_C);
% a=(C1==C2);
% sum(sum(~a))
% tic
% [C,phi,S12,S1,S2,t,f]=cohgramc_gpuMKII(amypuff,V1puff,[0.2 0.010],params);
% toc
% res=(C==coherencepool.puff);



