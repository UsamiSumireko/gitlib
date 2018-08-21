function [S1merge,S2merge,f]=V1spec(CONORI,sldwinlabel)
win=[201 600];
movingwin=[0.2 0.005];
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

load('D:\result\170913\noiseyacell.mat');
load('D:\result\170913\noiseyvcell.mat');

params.Fs = 1000; % sampling frequency
params.tapers = [3 5]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 3;
params.err =[2 0.05]; 
params.sldwinlable=sldwinlabel;
params.movingwin=movingwin;
% params.lr   L==0 R==1;
params.iscon=strcmpi(CONORI,'con'); % ori=0 con=1
S1merge=[];
S2merge=[];
proportionlist=[];
mergenoisyidx=[];
% ncount=0;
for iday=1:numel(folderlist)
    S1daymerge=[];
    S2daymerge=[];
    tic
    params.lr=strcmpi(LRlist(iday),'r');  % L==0 R==1;
    if params.lr
        cndlist=[1 2 3 4 5 6];
    else
        cndlist=[4 5 6 1 2 3];
    end
    if params.iscon
        load([folderlist{iday} '\V1conblockunitlfp.mat']);
        V1pufflfp=V1conblockunitlfp.lfp(win(1):win(2),:,1:3,:,:);
        V1neulfp=V1conblockunitlfp.lfp(win(1):win(2),:,4:6,:,:);
    else
        load([folderlist{iday} '\V1oriblockunitlfp.mat']);
        V1pufflfp=V1oriblockunitlfp.lfp(win(1):win(2),:,cndlist(1:3),:,:);
        V1neulfp=V1oriblockunitlfp.lfp(win(1):win(2),:,cndlist(4:6),:,:);
    end
    [s1,s2,s3,s4,s5]=size(V1pufflfp);
    V1pufflfp=permute(V1pufflfp,[1 2 3 5 4]);
    temp=reshape(V1pufflfp,s1,[]);
%     temp=transpose(filt50hz(temp',params.Fs));
    temp=rmlinesc(temp,params,[],'n',[50 100]);
    V1pufflfp=reshape(temp,s1,[],s4);
    [s1,s2,s3,s4,s5]=size(V1neulfp);
    V1neulfp=permute(V1neulfp,[1 2 3 5 4]);
    temp=reshape(V1neulfp,s1,[]);
%     temp=transpose(filt50hz(temp',params.Fs));
    temp=rmlinesc(temp,params,[],'n',[50 100]);
    V1neulfp=reshape(temp,s1,[],s4);
    if params.sldwinlable
        if movingwin(1)<=0.256
            params.pad = 3;
        end
        for ncell=1:s4
            [S1,t,f]=mtspecgramc(V1pufflfp(:,:,ncell),params.movingwin,params);
            [S2,~,~]=mtspecgramc(V1neulfp(:,:,ncell),params.movingwin,params);
            S1daymerge(:,:,ncell)=S1;
            S2daymerge(:,:,ncell)=S2;
            S1merge=cat(3,S1merge,S1);
            S2merge=cat(3,S2merge,S2);
        end
    else
        for ncell=1:s4
            [S1,f]=mtspectrumc(V1pufflfp(:,:,ncell),params);
            [S2,~]=mtspectrumc(V1neulfp(:,:,ncell),params);
            S1daymerge(:,ncell)=S1;
            S2daymerge(:,ncell)=S2;
            S1merge=cat(2,S1merge,S1);
            S2merge=cat(2,S2merge,S2);
        end
    end 
    x=1:s4;
    y=noisyvcell{iday};
    if isempty(y)
        y=0;
    end
    [X,Y]=meshgrid(x,y);
    vnoisyidx=bsxfun(@eq,X,Y);
    vnoisyidx=sum(vnoisyidx,1);
    mergenoisyidx=[mergenoisyidx vnoisyidx];
    toc
    if params.sldwinlable
    else
%         h1=figure;
%         for ncell=1:s4
%             subplot(5,10,ncell)
%             plot(f,S1daymerge(:,ncell),'r')
%             hold on
%             plot(f,S2daymerge(:,ncell),'g')
%             title(num2str(ncell))
%         end
%         saveas(gcf,['D:\result\170803resamplelfp\' num2str(iday) '.fig'])
%         close(gcf)
%         fidx=f>45 &f<75;
%         pp=log(mean(S1daymerge(fidx,:),1)./mean(S2daymerge(fidx,:),1));
%         proportionlist=[proportionlist pp];
%         figure
%         histogram(proportionlist)
%         saveas(gcf,['D:\result\170803resamplelfp\' num2str(iday) 'MKII.fig'])
%         close(gcf)
    end
end

if params.sldwinlable
    figure
    subplot(1,3,1)
    plot_matrix(mean(S1merge(:,:,~mergenoisyidx),3),t,f,'l')
    colormap jet
    subplot(1,3,2)
    plot_matrix(mean(S2merge(:,:,~mergenoisyidx),3),t,f,'l')
    colormap jet
    subplot(1,3,3)
    plot_matrix(mean(S1merge(:,:,~mergenoisyidx),3)./mean(S2merge(:,:,~mergenoisyidx),3),t,f,'l')
    colormap jet
else
%     figure
%     plot(f,log(mean(S1merge(:,~mergenoisyidx),2)),'r')
%     hold on
%     plot(f,log(mean(S2merge(:,~mergenoisyidx),2)),'g')
    puffste=std(log(S1merge),0,2)/sqrt(size(S1merge,2));
    neuste=std(log(S2merge),0,2)/sqrt(size(S1merge,2));
    figure
    xv=f';
    yv=log(mean(S1merge,2));
    patch([xv; flipud(xv)],[yv-puffste; flipud(yv+puffste)],[1 0.8 0.8]);
    hold on
    plot(xv,yv,'r')
    yv=log(mean(S2merge,2));
    patch([xv; flipud(xv)],[yv-neuste; flipud(yv+puffste)],[0.8 1 0.8]);
    plot(xv,yv,'g')
end
% save(['D:\result\170913\S1merge_' CONORI '_' num2str(win(1)-300) '_' num2str(win(2)-300) 'ms.mat'],'S1merge','-v7.3');
% save(['D:\result\170913\S2merge_' CONORI '_' num2str(win(1)-300) '_' num2str(win(2)-300) 'ms.mat'],'S2merge','-v7.3');
sprintf('done')

% fidx=f>45 & f<55;
% S1sus=mean(S1merge(fidx,:));
% S1valididx=S1sus<20;
% S2sus=mean(S2merge(fidx,:));
% S2valididx=S2sus<20;
% save(['D:\result\170913\S1validid.mat'],'S1valididx','-v7.3');
% save(['D:\result\170913\S2validid.mat'],'S2valididx','-v7.3');

% tidx=t>0.15 & t<0.2;
% fidx=f>60 & f<80;
% S1sus=squeeze(mean(mean(S1merge(tidx,fidx,:))));
% S2sus=squeeze(mean(mean(S2merge(tidx,fidx,:))));
% figure
% subplot(1,3,1)
% histogram(S1sus)
% subplot(1,3,2)
% histogram(S2sus)
% subplot(1,3,3)
% histogram(S1sus-S2sus)

function Y = filt50hz(X,sf)  % X's format: trials x bins
    Wo = 50/(sf/2);  BW = Wo/35;
    [b,a] = iirnotch(Wo,BW);  
    siz = size(X);
    Y = zeros(siz);
    for i=1:1:siz(1)
        Y(i,:) = filter(b,a,X(i,:));
%         Y(i,:) = lowp(Y(i,:),100,120,0.1,30,sf);
    end