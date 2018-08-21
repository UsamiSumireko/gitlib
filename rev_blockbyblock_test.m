clear
str='D:\data extracted\160908';
lr='l';

%%
params.Fs = 1000; % sampling frequency
params.tapers = [3 5]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 2;
params.err =[2 0.05]; 
prestr=[str '_prerev'];
load([prestr '\V1conblockunitpsth.mat']);
load([prestr '\V1conblockunitlfp.mat']);
load([prestr '\amyconblockunitlfp.mat']);
prerevpsth=V1conblockunitpsth;
prerevvlfp=V1conblockunitlfp;
prerevalfp=amyconblockunitlfp;
clearvars V1conblockunitpsth
clearvars V1conblockunitlfp

load([str '\V1conblockunitpsth.mat']);
load([str '\V1conblockunitlfp.mat']);
load([str '\amyconblockunitlfp.mat']);

lfpwin=301:700;

[as1,as2,as3,as4,as5]=size(prerevalfp.lfp(lfpwin,:,:,:,:));
aprepufflfp=permute(prerevalfp.lfp(lfpwin,:,1:3,:,:),[1 2 3 5 4]);
temp=reshape(aprepufflfp,as1,[]);
temp=rmlinesc(temp,params,[],'n',[50 100]);
aprepufflfp=reshape(temp,as1,[],as4);
apreneulfp=permute(prerevalfp.lfp(lfpwin,:,4:9,:,:),[1 2 3 5 4]);
temp=reshape(apreneulfp,as1,[]);
temp=rmlinesc(temp,params,[],'n',[50 100]);
apreneulfp=reshape(temp,as1,[],as4);
for ncell=1:as4
    [S1,f]=mtspectrumc(aprepufflfp(:,:,ncell),params);
    [S2,~]=mtspectrumc(apreneulfp(:,:,ncell),params);
    aS1prepower(:,ncell)=S1;
    aS2prepower(:,ncell)=S2;
end

% [s1,s2,s3,s4,s5]=size(prerevvlfp.lfp(lfpwin,:,:,:,:));
% prpuffelfp=permute(prerevvlfp.lfp(lfpwin,:,1:3,:,:),[1 2 3 5 4]);
% temp=reshape(prpuffelfp,s1,[]);
% temp=rmlinesc(temp,params,[],'n',[50 100]);
% prpuffelfp=reshape(temp,s1,[],s4);
% preneulfp=permute(prerevvlfp.lfp(lfpwin,:,4:9,:,:),[1 2 3 5 4]);
% temp=reshape(preneulfp,s1,[]);
% temp=rmlinesc(temp,params,[],'n',[50 100]);
% preneulfp=reshape(temp,s1,[],s4);
% for ncell=1:s4
%     [S1,f]=mtspectrumc(prpuffelfp(:,:,ncell),params);
%     [S2,~]=mtspectrumc(preneulfp(:,:,ncell),params);
%     S1prepower(:,ncell)=S1;
%     S2prepower(:,ncell)=S2;
% end

%%
% psthwin=1:600;
% peakwin=331:370;
% [s1,s2,s3,s4,s5]=size(V1conblockunitpsth.psth(psthwin,:,:,:,:));
% [ps1,ps2,ps3,ps4,ps5]=size(prerevpsth.psth(psthwin,:,:,:,:));
% prepuff=prerevpsth.psth(psthwin,:,1:3,:,:);
% preneu=prerevpsth.psth(psthwin,:,4:end,:,:);
% prepuff=reshape(prepuff,ps1,[]);
% preneu=reshape(preneu,ps1,[]);
% validtrialidx=sum(isnan(preneu))==0;
% preneu=preneu(:,validtrialidx);
% 
% peakdiff(1)=mean(mean(preneu(peakwin,:)))-mean(mean(prepuff(peakwin,:)));
% figure
% plot(-99:500,mean(prepuff,2),'g')
% hold on
% plot(-99:500,mean(preneu,2),'r')

%%
h1=figure;
% h2=figure;
for nblock=1:6
    [as1,as2,as3,as4,~]=size(amyconblockunitlfp.lfp(lfpwin,:,:,:,nblock));    
    ablockpufflfp=permute(amyconblockunitlfp.lfp(lfpwin,:,1:3,:,nblock),[1 2 3 5 4]);
    temp=reshape(ablockpufflfp,as1,[]);
    temp=rmlinesc(temp,params,[],'n',[50 100]);
    ablockpufflfp=reshape(temp,as1,[],as4);
    
    ablockneulfp=permute(amyconblockunitlfp.lfp(lfpwin,:,4:12,:,nblock),[1 2 3 5 4]);
    temp=reshape(ablockneulfp,as1,[]);
    temp=rmlinesc(temp,params,[],'n',[50 100]);
    ablockneulfp=reshape(temp,as1,[],as4);
    for ncell=1:as4
        [S1,f]=mtspectrumc(ablockpufflfp(:,:,ncell),params);
        [S2,~]=mtspectrumc(ablockneulfp(:,:,ncell),params);
        aS1blockpower(:,ncell,nblock)=S1;
        aS2blockpower(:,ncell,nblock)=S2;
    end
    figure(h1);
    subplot(2,6,nblock)
    semilogy(f,mean(aS1blockpower(:,:,nblock),2),'r')
    hold on
    semilogy(f,mean(aS2prepower,2),'g')
    ylabel('power(log)')
    xlabel('F hz')
%     figure(h2);
    subplot(2,6,nblock+6)
    semilogy(f,mean(aS2blockpower(:,:,nblock),2),'m')
    hold on
    semilogy(f,mean(aS1prepower,2),'b')
    ylabel('power(log)')
    xlabel('F hz')
end

%%
% h1=figure;
% % h2=figure;
% for nblock=1:6
%     [s1,s2,s3,s4,~]=size(V1conblockunitlfp.lfp(lfpwin,:,:,:,nblock));    
%     blockpufflfp=permute(V1conblockunitlfp.lfp(lfpwin,:,1:3,:,nblock),[1 2 3 5 4]);
%     temp=reshape(blockpufflfp,s1,[]);
%     temp=rmlinesc(temp,params,[],'n',[50 100]);
%     blockpufflfp=reshape(temp,s1,[],s4);
%     
%     blockneulfp=permute(V1conblockunitlfp.lfp(lfpwin,:,4:12,:,nblock),[1 2 3 5 4]);
%     temp=reshape(blockneulfp,s1,[]);
%     temp=rmlinesc(temp,params,[],'n',[50 100]);
%     blockneulfp=reshape(temp,s1,[],s4);
%     for ncell=1:s4
%         [S1,f]=mtspectrumc(blockpufflfp(:,:,ncell),params);
%         [S2,~]=mtspectrumc(blockneulfp(:,:,ncell),params);
%         S1blockpower(:,ncell,nblock)=S1;
%         S2blockpower(:,ncell,nblock)=S2;
%     end
%     figure(h1);
%     subplot(2,6,nblock)
%     semilogy(f,mean(S1blockpower(:,:,nblock),2),'r')
%     hold on
%     semilogy(f,mean(S2prepower,2),'g')
%     ylabel('power(log)')
%     xlabel('F hz')
% %     figure(h2);
%     subplot(2,6,nblock+6)
%     semilogy(f,mean(S2blockpower(:,:,nblock),2),'m')
%     hold on
%     semilogy(f,mean(S1prepower,2),'b')
%     ylabel('power(log)')
%     xlabel('F hz')
% end
%%
% figure
% for nblock=1:s5
%     subplot(floor(sqrt(s5))+1,floor(sqrt(s5))+1,nblock);
%     blockpuff=V1conblockunitpsth.psth(psthwin,:,1:3,:,nblock);
%     blockneu=V1conblockunitpsth.psth(psthwin,:,4:end,:,nblock);
%     blockpuff=reshape(blockpuff,s1,[]);
%     blockneu=reshape(blockneu,s1,[]);
%     plot(-99:500,mean(blockpuff,2),'r');
%     hold on
%     plot(-99:500,mean(blockneu,2),'g');
%     validtrialidx=sum(isnan(blockneu))==0;
%     blockneu=blockneu(:,validtrialidx);
%     peakdiff(1+nblock)=mean(mean(blockpuff(peakwin,:)))-mean(mean(blockneu(peakwin,:)));
% end
% figure
% bar(peakdiff);


% str='D:\result\170905 blockbyblock\170505\';
% for i=9:16
%     saveas(i,[str num2str(i-8) '.fig']);
% end