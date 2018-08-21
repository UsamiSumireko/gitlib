clear
folderlist={
%     'I:\160712';
%     'I:\160713';
%     'I:\160714';
%     'I:\160715';
%     'I:\160716';
%     'I:\160719';
%     'I:\160720';
%     'I:\160721';
%     'I:\160722';
    
    'I:\170426'
    'I:\170427'
    'I:\170429'
    'I:\170430'
    'I:\170502'
    'I:\170503'
};
stimC1=[];
stimC2=[];
shuffleC1=[];
shuffleC2=[];
for idays=1:numel(folderlist)
    shufflelabel=cell(6,1);
    str=folderlist{idays};
    load([folderlist{idays} '\LFPCoh_shuffle101-500ms.mat']);
    shufflecoh=LFPCoh;
    load([folderlist{idays} '\LFPCoh101-500ms.mat']);
    for j=1:numel(LFPCoh)
        for k=1:numel(LFPCoh{j})
            stimC1=[stimC1 LFPCoh{j}{k}.C1];
            stimC2=[stimC2 LFPCoh{j}{k}.C2];
            shuffleC1=[shuffleC1 squeeze(mean(shufflecoh{j}{k}.C1,3))];
            shuffleC2=[shuffleC2 squeeze(mean(shufflecoh{j}{k}.C2,3))];
        end
    end
end
stimf=LFPCoh{1}{1}.f1{1};
figure
plot(stimf,mean(stimC1,2),'r')
hold on
plot(stimf,mean(stimC2,2),'g')
title('stim time Coh')

figure
plot(stimf,mean(stimC1,2),'r')
hold on
plot(stimf,mean(stimC2,2),'g')
hold on
plot(stimf,mean(shuffleC1,2),'m')
hold on
plot(stimf,mean(shuffleC2,2),'c')
title('stim time Coh')
%%
% i=1;j=2;k=3;
% f=shufflecoh{1}{1}.f1{1};
% c1=LFPCoh{i}{j}.C1(:,k);
% c2=LFPCoh{i}{j}.C2(:,k);
% shuffle1=squeeze(mean(shufflecoh{i}{j}.C1(:,k,:),3));
% ste=std(shufflecoh{i}{j}.C1(:,k,:),0,3)/sqrt(size(shufflecoh{i}{j}.C1,3));
% figure
% fill([f fliplr(f)],[shuffle1-ste;flipud(shuffle1+ste)],[0.8 0.8 1],'EdgeAlpha',0);
% hold on
% plot(f,shuffle1,'-b');
% plot(f,c2,'-g');
% plot(f,c1,'-r');