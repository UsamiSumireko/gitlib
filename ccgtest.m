% load('G:\V1&V2Contour\DataAnalysis\CCG\train1.mat');
% load('G:\V1&V2Contour\DataAnalysis\CCG\train2.mat');
clear
iNormMethod=3;
iShufMethod=2;
nShuffle=1000;
iJitterWin=25;
% [fRawCCG, fShuffleCCG, fShuffleCi] = CCG(Train1, Train2, iNormMethod, iShufMethod, nShuffle, iJitterWin);

%%
twin=1:600;
str='D:\data extracted\160720';
load([str '\amyconblockunitpsth.mat']);
load([str '\V1conblockunitpsth.mat']);
amydata=amyconblockunitpsth.psth(twin,:,:,:,:);
V1data=V1conblockunitpsth.psth(twin,:,:,:,:);
[s1,s2,s3,s4,s5]=size(amydata);
[~,~,~,ss4,~]=size(V1data);
cscnd=1:3;
neucnd=4:12;
PairNum=0;
ccgav=struct('ccg',[],'blockcellidx',[]);
for iblock=1:s5
    sprintf(['Nblock=' num2str(iblock)])
    amyvalid=amyconblockunitpsth.stimvalid & amyconblockunitpsth.spikevalid;
    amydatab=amydata(:,:,:,logical(amyvalid(:,iblock)),iblock);
    validcellnum=sum(amyvalid(:,iblock));
    validcelllist=find(amyvalid(:,iblock));
    amydatab=reshape(amydatab,s1,s2,s3,validcellnum);
    v1datab=V1data(:,:,:,:,iblock);
    v1datab=reshape(v1datab,s1,s2,s3,ss4);
    for acell=1:validcellnum
        sprintf(num2str(acell))
        amydatac=amydatab(:,:,cscnd,acell);
        amydatac=reshape(amydatac,s1,[]);        
        for vcell=1:ss4
            PairNum=PairNum+1;
            v1datac=v1datab(:,:,cscnd,vcell);
            v1datac=reshape(v1datac,s1,[]);
            [fRawCCG,fShuffleCCG,fShuffleCi] = CCG(amydatac',v1datac',iNormMethod,iShufMethod,nShuffle,iJitterWin);
            ccgav.ccg(PairNum,:)=mean(fRawCCG);%-fShuffleCCG;
            ccgav.shuffleccg(PairNum,:)=fShuffleCCG;
            ccgav.blockcellidx.block(PairNum)=iblock;
            ccgav.blockcellidx.acell(PairNum)=validcelllist(acell);
            ccgav.blockcellidx.vcell(PairNum)=vcell;
        end
    end
end
%%
HalfSmoothWindow=10;
HSmoothTrain=1:HalfSmoothWindow;
SmoothTrain=[HSmoothTrain,fliplr(HSmoothTrain(1:end-1))];
SmoothTrain=SmoothTrain./sum(SmoothTrain);
figure;
% clr=colormap(hsv(5));
% idx=ccgav.blockcellidx.acell==1;
rawy=mean(ccgav.ccg,1)'-mean(ccgav.shuffleccg,1)';
y=conv(rawy,SmoothTrain,'same');
plot(-twin(end)+1:twin(end)-1,y);
% ylim([-5*10^(-3) 5*10^(-3)]);
% title(['v1(' num2str(i) ')&v2(' num2str(j) ')']);