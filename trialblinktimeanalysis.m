clear
folderlist={
%     'D:\data extracted\180111';
%     'D:\data extracted\180112';
%     'D:\data extracted\180117';
%     'D:\data extracted\180118';
%     'D:\data extracted\180119';
%     'D:\data extracted\180120';
%     'D:\data extracted\180212';

%     'D:\data extracted\180425';
%     'D:\data extracted\180426';
%     'D:\data extracted\180427';
    
%     'D:\data extracted\180504';
%     'D:\data extracted\180505';
%     'D:\data extracted\180506';
    
%     'D:\data extracted\180514';
%     'D:\data extracted\180516';
    'D:\data extracted\180525';
    };
% LRlist=['rrrr'];
LRlist=['r'];
RFlabel=0;

blinkpsth=cell(4,1);
nonblinkpsth=cell(4,1);
blinktime=cell(4,1);
blinkN=cell(4,1);
for iday=1:numel(folderlist)
    if RFlabel
        lr=LRlist(iday);
        if strcmp(lr,'l')
            CNDlist=[3 1 1 1 4 2 2 2];
        elseif strcmp(lr,'r')
            CNDlist=[1 3 3 3 2 4 4 4];
        else
            error('LR input is wrong, L or R');
        end
    else
        lr=LRlist(iday);
        if strcmp(lr,'l')
            CNDlist=[2 1 1 1 4 3 3 3];
        elseif strcmp(lr,'r')
            CNDlist=[1 2 2 2 3 4 4 4];
        else
            error('LR input is wrong, L or R');
        end
    end
    load([folderlist{iday} '\eyeblinkdata.mat']);
    load([folderlist{iday} '\V1conblockunitpsth.mat'])
    load([folderlist{iday} '\V1oriblockunitpsth.mat'])
    nanidx=isnan(eyeblinkdata.eyeclosetime);
    shorttracetimeidx=eyeblinkdata.tracetime<1000;
    blinkidx=eyeblinkdata.Nblink>0;
    psth=permute(V1conblockunitpsth.psth,[1 4 2 3 5]);
    for m=1:4
        cndidx=CNDlist==m;
        idx1=~shorttracetimeidx & blinkidx;
        idx1=idx1(:,cndidx,:);
        idx2=~shorttracetimeidx & ~blinkidx;
        idx2=idx2(:,cndidx,:);
        psthm=psth(:,:,:,cndidx,:);
        psthmb=psthm(:,:,idx1);
        psthmb=reshape(psthmb,size(psthmb,1),[]);
        psthmn=psthm(:,:,idx2);
        psthmn=reshape(psthmn,size(psthmn,1),[]);
        blinkpsth{m}=cat(2,blinkpsth{m},psthmb);
        nonblinkpsth{m}=cat(2,nonblinkpsth{m},psthmn);
        temp=eyeblinkdata.eyeclosetime(:,cndidx,:);
        blinktime{m}=[blinktime{m}; temp(~shorttracetimeidx(:,cndidx,:))];
        temp=eyeblinkdata.Nblink(:,cndidx,:);
        blinkN{m}=[blinkN{m}; temp(~shorttracetimeidx(:,cndidx,:))];        
    end
end
figure
for m=1:4
xv=-199:400;
valididx=squeeze(sum(blinkpsth{m}(1:600,:)==65535,1)==0);
yv=mean(blinkpsth{m}(1:600,valididx),2);
subplot(2,2,m)
plot(xv,yv,'r')
hold on
valididx=squeeze(sum(nonblinkpsth{m}(1:600,:)==65535,1)==0);
yv=mean(nonblinkpsth{m}(1:600,valididx),2);
plot(xv,yv,'g')
title(['RF' num2str(90*(m>2)+45) 'Contra' num2str(90*mod(m+1,2)+45)])
legend('Blinked trial','Unblinked trial')
end
figure
for m=1:4
    yv1(m)=mean(blinktime{m})/2;
    yv2(m)=mean(blinkN{m});
end
subplot(1,2,1)
bar(2:2:8,yv1);
ax=gca;
ax.XTick=2:2:8;
ax.XTickLabels={'R45 C45','R45 C135','R135 C45','R135 C135'};
title('Eye close period')
ylabel('Time (ms)')

subplot(1,2,2)
bar(2:2:8,yv2);
ax=gca;
ax.XTick=2:2:8;
ax.XTickLabels={'R45 C45','R45 C135','R135 C45','R135 C135'};
title('Blink time')
ylabel('N')



for i=1:3
    subplot(1,3,i)
    ylim([0 0.35])
end