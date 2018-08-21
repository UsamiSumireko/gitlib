function [S1merge,S2merge]=amyspec(CONORI,sortlabel)
win=[401 800];
% win=201:600;
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
params.trialave = 1;
params.fpass = [4 100];
params.pad = 3;
params.err =[2 0.05]; 
% params.lr   L==0 R==1;
params.iscon=strcmpi(CONORI,'con');  % ori=0 con=1

blockcellidx=struct('day',[],'block',[],'cell',[]);
S1merge=[];
S2merge=[];
proportionlist=[];
for iday=1:numel(folderlist)
    iday
    S1daymerge=[];
    S2daymerge=[];
    params.lr=strcmp(LRlist(iday),'r');  % L==0 R==1;
    if params.lr
        cndlist=[1 2 3 4 5 6];
    else
        cndlist=[4 5 6 1 2 3];
    end
    if params.iscon
        load([folderlist{iday} '\n2b_Svalencevalid.mat']);
        load([folderlist{iday} '\amyconblockunitlfp.mat']);
        if sortlabel
            valididx=logical(amyconblockunitlfp.spikevalid);% & valencevalid;
        else
            valididx=true(size(amyconblockunitlfp.spikevalid));
        end
        amypufflfp=amyconblockunitlfp.lfp(win,:,1:3,valididx);
        amyneulfp=amyconblockunitlfp.lfp(win,:,4:6,valididx);
    else

        load([folderlist{iday} '\amyoriblockunitlfp.mat']);
        if sortlabel
            valididx=logical(amyoriblockunitlfp.spikevalid);% & valencevalid;
        else
            valididx=true(size(amyconblockunitlfp.spikevalid));
        end
        amypufflfp=amyoriblockunitlfp.lfp(win,:,cndlist(1:3),valididx);
        amyneulfp=amyoriblockunitlfp.lfp(win,:,cndlist(4:6),valididx);
    end
    [s1,s2,s3,s4,s5]=size(amypufflfp);
    amypufflfp=permute(amypufflfp,[1 2 3 5 4]);
    temp=reshape(amypufflfp,s1,[]);
    temp=rmlinesc(temp,params,[],'n',[50 100]);
    amypufflfp=reshape(temp,s1,[],s4);
    amyneulfp=permute(amyneulfp,[1 2 3 5 4]);
    temp=reshape(amyneulfp,s1,[]);
    temp=rmlinesc(temp,params,[],'n',[50 100]);
    amyneulfp=reshape(temp,s1,[],s4);
    for ncell=1:s4
        [S1,f]=mtspectrumc(amypufflfp(:,:,ncell),params);
        [S2,~]=mtspectrumc(amyneulfp(:,:,ncell),params);
        S1daymerge(:,ncell)=S1;
        S2daymerge(:,ncell)=S2;
        S1merge=cat(2,S1merge,S1);
        S2merge=cat(2,S2merge,S2);
%         blockcellidx.cell=[blockcellidx.cell; ncell];
%         blockcellidx.block=[blockcellidx.block; ncell];
    end
%     h1=figure;
%     for ncell=1:24
%         subplot(3,8,ncell)
%         plot(f,S1daymerge(:,ncell),'r')
%         hold on
%         plot(f,S2daymerge(:,ncell),'g')
%         title(num2str(ncell))
%     end
%     saveas(gcf,['D:\result\170803resamplelfp\amy\' num2str(iday) '.fig'])
%     close(gcf)
%     fidx=f>45 &f<75;
%     pp=log(mean(S1daymerge(fidx,:),1)./mean(S2daymerge(fidx,:),1));
%     proportionlist=[proportionlist pp];
%     figure
%     histogram(proportionlist)
%     saveas(gcf,['D:\result\170803resamplelfp\amy\' num2str(iday) 'MKII.fig'])
%     close(gcf)
end
figure
yv=mean(S1merge,2);
yste=std(S1merge,[],2)/sqrt(size(S1merge,2));
patch([f'; flipud(f')],[yv+yste; flipud(yv-yste)],[1 0.75 0.75]);
hold on
plot(f,yv,'r')

yv=mean(S2merge,2);
yste=std(S2merge,[],2)/sqrt(size(S2merge,2));
patch([f'; flipud(f')],[yv+yste; flipud(yv-yste)],[0.75 1 0.75]);
hold on
plot(f,yv,'g')
xlabel('f (Hz)')
ylabel('power')

figure
semilogy(f,mean(S1merge,2),'r')
hold on
semilogy(f,mean(S2merge,2),'g')
xlabel('f (Hz)')
ylabel('log power')
sprintf('done')

    