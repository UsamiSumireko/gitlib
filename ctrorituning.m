clear
%%
str='H:\180304\FlexTuning';
load([str '\XYvalidid.mat'])
selecti='xy';
valididx=XYvalidid;
valididx=valididx(1:12,:);
strRootPath = {
    [str '\FlexTuning002'];
    [str '\FlexTuning003'];
    [str '\FlexTuning004'];
    };
cmppath =  'D:\SN 4566-001277.cmp';
cmpinfo = LoadCmp(cmppath,1,0);
elec = cmpinfo{1,1}.RealElec(1:12,:);
unit = lower(char(zeros(size(elec,1),size(elec,2))+84));  %'t'
%o
win = [0.200 1.000];
binwidth=0.001;
currunit='t';
celllist=elec(valididx);
ctroripsth=[];
for nblock=1:3
    path = strcat(strRootPath{nblock},'.nev');
    psth=[];
    for icell=1:numel(celllist)
        spkstruct=BinSpike(path,celllist(icell),currunit,win,binwidth);
        splitspk=SplitInfo(spkstruct,[1 0 1 0]);
        for m=1:numel(splitspk)
            psth(:,:,m,icell)=transpose(splitspk{1,m}{1,1}.Train);
        end
    end
    ctroripsth=cat(2,ctroripsth,psth);
end
%%
%             idx=sum(splitspk{1,m}{1,1}.Train==(2^16-1), 2)==0;
winidx=31:71;
mu45=[];
mu135=[];
for acontrast=0:7
    cndidx=[acontrast*24+4 acontrast*24+16];
    psth45=reshape(ctroripsth(:,:,cndidx,:),size(ctroripsth,1),size(ctroripsth,2)*2,size(ctroripsth,4));
    valididx=sum(sum(psth45==2^16-1,3),1)==0;
    psth45=psth45(:,valididx,:);
    cndidx=[acontrast*24+10 acontrast*24+22];
    psth135=reshape(ctroripsth(:,:,cndidx,:),size(ctroripsth,1),size(ctroripsth,2)*2,size(ctroripsth,4));
    valididx=sum(sum(psth135==2^16-1,3),1)==0;
    psth135=psth135(:,valididx,:);
    mu45(acontrast+1)=mean(mean(mean(psth45(winidx,:,:),3),2),1);
    mu135(acontrast+1)=mean(mean(mean(psth135(winidx,:,:),3),2),1);
end
figure
subplot(1,2,1)
plot(1:8,mu45/binwidth,'g')
title('45deg contrast tuning')
subplot(1,2,2)
plot(1:8,mu135/binwidth,'b')
title('135deg contrast tuning')
