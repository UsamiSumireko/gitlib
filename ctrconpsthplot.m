clear
savelabel=1;
str='H:\180316';
lr='r';
win=[-0.100 0.500];
binwidth=0.001;
% ctrvalue=[1 2 4 8 16 32 64 100];
ctrvalue=[];
puffcnd=[1:3];
neucnd=[4:6];
conblocktrialnum=10;
oriblocktrialnum=20;

cmppath =  'D:\SN 4566-001277.cmp';
cmpinfo = LoadCmp(cmppath,1,0);
elec = cmpinfo{1,1}.RealElec(1:12,:);
unit = lower(char(zeros(size(elec,1),size(elec,2))+84));  

if strcmpi(lr,'l')
    CNDlist=[2 2 2 1 1 1 1 1 1 1 1 1]; % 1=45 2=135
%     CNDlist=[2 1 1 1 2 1 1 1]; 
elseif strcmp(lr,'r')
    CNDlist=[1 1 1 2 2 2 2 2 2 2 2 2];
%     CNDlist=[1 2 2 2 1 2 2 2]; 
else
    error('LR input is wrong, L or R');
end

load([str '\ORIvalidid.mat'])
load([str '\XYvalidid.mat'])
valididx=ORIvalidid & XYvalidid;
valididx=valididx(1:12,:);
celllist=elec(valididx);

namelist=dir(str);
orinamelist=[];
orictrlist=[];
connamelist=[];
conctrlist=[];

if ~isempty(ctrvalue)
    for i=3:numel(namelist)
        if namelist(i).isdir
            name=namelist(i).name;
            if strcmp(name(1:6),'DftORI')
                anchor=strfind(name,'_');
                ctrnumidx=anchor(2)+1:anchor(3)-1;
                ctr=str2num(name(ctrnumidx));
                orinamelist=[orinamelist; {name}];
                orictrlist=[orictrlist; ctr];
            elseif strcmp(name(1:6),'ORICON')
                anchor=strfind(name,'_');
                ctrnumidx=anchor(1)+1:anchor(2)-1;
                ctr=str2num(name(ctrnumidx));
                connamelist=[connamelist; {name}];
                conctrlist=[conctrlist; ctr];
            end
        end
    end
    V1conctrunitpsth=[];%ones(floor((win(2)-win(1))/binwidth),conblocktrialnum,numel(CNDlist),8,numel(celllist))*2^16-1;
    V1orictrunitpsth=[];%ones(floor((win(2)-win(1))/binwidth),oriblocktrialnum,2,8,numel(celllist))*2^16-1;
    for jcrt=1:numel(ctrvalue)
        for ncell=1:numel(celllist)
            idx=find(orictrlist==ctrvalue(jcrt));
            for kblock=1:numel(idx)
                oripath=[str '\' orinamelist{idx(kblock)} '.nev'];
                spkstruct=BinSpike(oripath,celllist(ncell),'t',win,binwidth);
                splitspk=SplitInfo(spkstruct,[1 0 1 0]);
%                 for m=1:2
                    [s1,s2]=size(splitspk{CNDlist(puffcnd(1))*2}{1}.Train);
                    V1orictrunitpsth(1:s2,(kblock-1)*s1+1:kblock*s1,1,jcrt,ncell)=transpose(splitspk{CNDlist(puffcnd(1))*2}{1}.Train);
                    [s1,s2]=size(splitspk{CNDlist(neucnd(1))*2}{1}.Train);
                    V1orictrunitpsth(1:s2,(kblock-1)*s1+1:kblock*s1,2,jcrt,ncell)=transpose(splitspk{CNDlist(neucnd(1))*2}{1}.Train);
%                 end
            end
            idx=find(conctrlist==ctrvalue(jcrt));
            for kblock=1:numel(idx)
                conpath=[str '\' connamelist{idx(kblock)} '.nev'];
                spkstruct=BinSpike(conpath,celllist(ncell),'t',win,binwidth);
                splitspk=SplitInfo(spkstruct,[1 0 1 0]);
                for m=1:numel(splitspk)
                    [s1,s2]=size(splitspk{m}{1}.Train);
                    V1conctrunitpsth(1:s2,(kblock-1)*s1+1:kblock*s1,m,jcrt,ncell)=transpose(splitspk{m}{1}.Train);
                end
            end
        end
    end
    if savelabel
        if ~exist(['D:\data extracted\' str(end-5:end)],'dir')
            mkdir(['D:\data extracted\' str(end-5:end)]);
        end
        save(['D:\data extracted\' str(end-5:end) '\V1conctrunitpsth'],'V1conctrunitpsth','-v7.3')
        save(['D:\data extracted\' str(end-5:end) '\V1orictrunitpsth'],'V1orictrunitpsth','-v7.3')
    end
    figure
    for x=1:8
        subplot(2,8,x)
        xv=transpose(-199:400);
        temp=reshape(V1conctrunitpsth(:,:,puffcnd,x,:),size(V1conctrunitpsth,1),[],size(V1conctrunitpsth,5));
        trialidx=0==sum(sum(temp==2^16-1,1),3);
        yv1=mean(mean(temp(:,trialidx,:),3),2);
        yv1=conv(yv1,ones(10,1),'same')/10/binwidth;
        yste1=std(squeeze(mean(temp(:,trialidx,:)/binwidth,2)),[],2)/sqrt(size(temp,3));
        patch([xv; flipud(xv)],[yv1+yste1; flipud(yv1-yste1)],[1 0.75 0.75],'EdgeAlpha',0)
        hold on
        plot(xv,yv1,'r')
        
        temp=reshape(V1orictrunitpsth(:,:,1,x,:),size(V1orictrunitpsth,1),[],size(V1orictrunitpsth,5));
        trialidx=0==sum(sum(temp==2^16-1,1),3);
        yv2=mean(mean(temp(:,trialidx,:),3),2);
        yv2=conv(yv2,ones(10,1),'same')/10/binwidth;
        yste2=std(squeeze(mean(temp(:,trialidx,:)/binwidth,2)),[],2)/sqrt(size(temp,3));
        patch([xv; flipud(xv)],[yv2+yste2; flipud(yv2-yste2)],[0.75 1 0.75],'EdgeAlpha',0)
        hold on
        plot(xv,yv2,'g')
        xlim([-199 400])
        ylim([0 300])
        title(['ctr' num2str(ctrvalue(x)) 'PuffOri PSTH'])
        xlabel('Time (ms)')
        ylabel('FR (Hz)')
        
        if(x<=2)
            hpkidx=1:600>250 & 1:600<401;
        else
        bl=mean(yv2(1:200));
        pk=max(yv2);
        rawhpkidx=yv2>(bl+pk)/2;
        L=bwlabel(rawhpkidx);
        hpkidx=L==1;
        end
        puffpeakFRdiff(x)=mean(yv1(hpkidx)-yv2(hpkidx));
        puffpeakFRdiffidx(x)=mean(yv1(hpkidx)-yv2(hpkidx))/(mean(yv1(hpkidx)-yv2(hpkidx))+mean(yv1(hpkidx)+yv2(hpkidx)));
        
        
        subplot(2,8,x+8)
        temp=reshape(V1conctrunitpsth(:,:,neucnd,x,:),size(V1conctrunitpsth,1),[],size(V1conctrunitpsth,5));
        trialidx=0==sum(sum(temp==2^16-1,1),3);
        yv1=mean(mean(temp(:,trialidx,:),3),2);
        yv1=conv(yv1,ones(10,1),'same')/10/binwidth;
        yste1=std(squeeze(mean(temp(:,trialidx,:)/binwidth,2)),[],2)/sqrt(size(temp,3));
        patch([xv; flipud(xv)],[yv1+yste1; flipud(yv1-yste1)],[1 0.75 0.75],'EdgeAlpha',0)
        hold on
        plot(xv,yv1,'r')
        
        temp=reshape(V1orictrunitpsth(:,:,2,x,:),size(V1orictrunitpsth,1),[],size(V1orictrunitpsth,5));
        trialidx=0==sum(sum(temp==2^16-1,1),3);
        yv2=mean(mean(temp(:,trialidx,:),3),2);
        yv2=conv(yv2,ones(10,1),'same')/10/binwidth;
        yste2=std(squeeze(mean(temp(:,trialidx,:)/binwidth,2)),[],2)/sqrt(size(temp,3));
        patch([xv; flipud(xv)],[yv2+yste2; flipud(yv2-yste2)],[0.75 1 0.75],'EdgeAlpha',0)
        hold on
        plot(xv,yv2,'g')
        xlim([-199 400])
        ylim([0 300])
        title(['ctr' num2str(ctrvalue(x)) 'NeuOri PSTH'])
        xlabel('Time (ms)')
        ylabel('FR (Hz)')
        
        if(x<=2)
            hpkidx=1:600>250 & 1:600<401;
        else
            bl=mean(yv2(1:200));
            pk=max(yv2);
            rawhpkidx=yv2>(bl+pk)/2;
            L=bwlabel(rawhpkidx);
            hpkidx=L==1;
        end
        neupeakFRdiff(x)=mean(yv1(hpkidx)-yv2(hpkidx));
        neupeakFRdiffidx(x)=mean(yv1(hpkidx)-yv2(hpkidx))/(mean(yv1(hpkidx)-yv2(hpkidx))+mean(yv1(hpkidx)+yv2(hpkidx)));
    end
    figure
    subplot(1,2,1)
    plot(1:8,puffpeakFRdiff-neupeakFRdiff)
    xlabel('contrast level')
    ylabel('conditioning effect diff')
    subplot(1,2,2)
    plot(1:8,puffpeakFRdiffidx-neupeakFRdiffidx)
    xlabel('contrast level')
    ylabel('conditioning effect idx diff')
    
else %ctrvalue is empty
    for i=3:numel(namelist)
        if namelist(i).isdir
            name=namelist(i).name;
            if strcmp(name(1:7),'DftORI0')
                orinamelist=[orinamelist; {name}];
            elseif strcmp(name(1:7),'ORICON0')
                connamelist=[connamelist; {name}];
            end            
        end
    end
    V1conunitpsth=[];%ones(floor((win(2)-win(1))/binwidth),conblocktrialnum*numel(connamelist),8,numel(celllist))*2^16-1;
    V1oriunitpsth=[];%ones(floor((win(2)-win(1))/binwidth),oriblocktrialnum*numel(orinamelist),2,numel(celllist))*2^16-1;
    for ncell=1:numel(celllist)
        for iblock=1:numel(orinamelist)            
            oripath=[str '\' orinamelist{iblock} '.nev'];
            spkstruct=BinSpike(oripath,celllist(ncell),'t',win,binwidth);
            splitspk=SplitInfo(spkstruct,[1 0 1 0]);
%             for m=1:2
                [s1,s2]=size(splitspk{CNDlist(puffcnd(1))*2}{1}.Train);
                V1oriunitpsth(1:s2,(iblock-1)*s1+1:iblock*s1,1,ncell)=transpose(splitspk{CNDlist(puffcnd(1))*2}{1}.Train);
%                 V1oriunitpsth(1:s2,(iblock-1)*2*s1+1:iblock*2*s1,1,ncell)=transpose([splitspk{CNDlist(puffcnd(1))*6-2}{1}.Train; splitspk{CNDlist(puffcnd(1))*6+10}{1}.Train;]);
                [s1,s2]=size(splitspk{CNDlist(neucnd(1))*2}{1}.Train);
                V1oriunitpsth(1:s2,(iblock-1)*s1+1:iblock*s1,2,ncell)=transpose(splitspk{CNDlist(neucnd(1))*2}{1}.Train);
%                 V1oriunitpsth(1:s2,(iblock-1)*2*s1+1:iblock*2*s1,2,ncell)=transpose([splitspk{CNDlist(neucnd(1))*6-2}{1}.Train; splitspk{CNDlist(neucnd(1))*6+10}{1}.Train;]);

                %             end
        end
        for iblock=1:numel(connamelist)
            conpath=[str '\' connamelist{iblock} '.nev'];
            spkstruct=BinSpike(conpath,celllist(ncell),'t',win,binwidth);
            splitspk=SplitInfo(spkstruct,[1 0 1 0]);
            for m=1:numel(splitspk)
                [s1,s2]=size(splitspk{m}{1}.Train);
                V1conunitpsth(1:s2,(iblock-1)*s1+1:iblock*s1,m,ncell)=transpose(splitspk{m}{1}.Train);
            end
        end
    end
figure
subplot(2,1,1)
xv=transpose(-199:400);
temp=reshape(V1conunitpsth(:,:,puffcnd,:),size(V1conunitpsth,1),[],size(V1conunitpsth,4));
trialidx=0==sum(sum(temp==2^16-1,1),3);
yv1=mean(mean(temp(:,trialidx,:),3),2);
yv1=conv(yv1,ones(10,1),'same')/10/binwidth;
yste1=std(squeeze(mean(temp(:,trialidx,:)/binwidth,2)),[],2)/sqrt(size(temp,3));
patch([xv; flipud(xv)],[yv1+yste1; flipud(yv1-yste1)],[1 0.75 0.75],'EdgeAlpha',0)
hold on
plot(xv,yv1,'r')

temp=reshape(V1oriunitpsth(:,:,1,:),size(V1oriunitpsth,1),[],size(V1oriunitpsth,4));
trialidx=0==sum(sum(temp==2^16-1,1),3);
yv2=mean(mean(temp(:,trialidx,:),3),2);
yv2=conv(yv2,ones(10,1),'same')/10/binwidth;
yste2=std(squeeze(mean(temp(:,trialidx,:)/binwidth,2)),[],2)/sqrt(size(temp,3));
patch([xv; flipud(xv)],[yv2+yste2; flipud(yv2-yste2)],[0.75 1 0.75],'EdgeAlpha',0)
hold on
plot(xv,yv2,'g')
title('puff ori')
xlabel('Time (ms)')
ylabel('FR (Hz)')

subplot(2,1,2)
xv=transpose(-199:400);
temp=reshape(V1conunitpsth(:,:,neucnd,:),size(V1conunitpsth,1),[],size(V1conunitpsth,4));
trialidx=0==sum(sum(temp==2^16-1,1),3);
yv1=mean(mean(temp(:,trialidx,:),3),2);
yv1=conv(yv1,ones(10,1),'same')/10/binwidth;
yste1=std(squeeze(mean(temp(:,trialidx,:)/binwidth,2)),[],2)/sqrt(size(temp,3));
patch([xv; flipud(xv)],[yv1+yste1; flipud(yv1-yste1)],[1 0.75 0.75],'EdgeAlpha',0)
hold on
plot(xv,yv1,'r')

temp=reshape(V1oriunitpsth(:,:,2,:),size(V1oriunitpsth,1),[],size(V1oriunitpsth,4));
trialidx=0==sum(sum(temp==2^16-1,1),3);
yv2=mean(mean(temp(:,trialidx,:),3),2);
yv2=conv(yv2,ones(10,1),'same')/10/binwidth;
yste2=std(squeeze(mean(temp(:,trialidx,:)/binwidth,2)),[],2)/sqrt(size(temp,3));
patch([xv; flipud(xv)],[yv2+yste2; flipud(yv2-yste2)],[0.75 1 0.75],'EdgeAlpha',0)
hold on
plot(xv,yv2,'g')
title('neu ori')
xlabel('Time (ms)')
ylabel('FR (Hz)')

end

% str='H:\180311';
% namelist=dir(str);
% for i=3:numel(namelist)
%     if strcmp(namelist(i).name(1:6),'DftORI') && ~namelist(i).isdir
% %         if 3>numel(strfind(namelist(i).name,'_'))
% %             newname=[namelist(i).name(1:10) '_' namelist(i).name(11:21)];
% %             movefile([str '\' namelist(i).name],[str '\' newname])
% %         end
%             newname=[namelist(i).name(1:10) namelist(i).name(19:end)];
%             movefile([str '\' namelist(i).name],[str '\' newname])
%     end    
% end