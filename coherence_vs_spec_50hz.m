clear
CONORI='con';
spklfp=[{'lfp'},{'lfp'}];%  (1)amy spk/lfp (2)V1 spk/lfp
savelabel=1;
sldwinlabel=0;
trialshufflelabel=1;
bwin=201:600;
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
str='D:\result\170828 mergeblock lfpspkcoh\';
multineucnd=[4, 5, 6;
             7, 8, 9;
             10, 11, 12];
LRlist=['rrrrrrrrr' 'llll' 'rrrrrr' 'llll'];         

params.Fs = 1000; % sampling frequency
params.tapers = [3 5]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 2;
params.err =[2 0.05];
npairsession=32;
nrml=2000;
% mergeperi50=[];
thr=0.23;
mergeSa=[];
mergeSv=[];
mergeC=[];
mergepuffcoh=[];
mergeneucoh=[];
mergevalididx=[];
noisyacell=[];
noisyvcell=[];
for iday=1:numel(folderlist)
    path1=[folderlist{iday} '\' 'amy' CONORI 'blockunitlfp.mat'];
    path2=[folderlist{iday} '\' 'V1' CONORI 'blockunitlfp.mat'];
    path3=[str folderlist{iday}(end-5:end) '\' CONORI '_' spklfp{1} '_' spklfp{2} '_coherencepool_g.mat'];
    if exist(path3,'file')
        load(path1)
        load(path2)
        load(path3)
        peri50fidx=coherencepool.f>45 & coherencepool.f<55;
        peri50=mean(coherencepool.puff(peri50fidx,:));
%         mergeperi50=[mergeperi50, peri50];
        noiseidx=find(peri50>thr); 
        noisyacell{iday}=zeros(size(peri50));
        noisyvcell{iday}=[];
        if ~isempty(noiseidx)            
            blocklist=arrayfun(@(x) x.blockidx,coherencepool.blockcellidx);
            acelllist=arrayfun(@(x) x.acellidx,coherencepool.blockcellidx);
            vcelllist=arrayfun(@(x) x.vcellidx,coherencepool.blockcellidx);
            noisyacell{iday}=noiseidx;
            noisyvcell{iday}=unique(vcelllist(noiseidx));            
            for npair=1:numel(noiseidx)
                blockidx=blocklist(noiseidx(npair));
                acellidx=acelllist(noiseidx(npair));
                vcellidx=vcelllist(noiseidx(npair));
                [s1,s2,s3,s4,s5]=size(amyconblockunitlfp.lfp);
                m2b=floor(s5/2);
                temp=permute(amyconblockunitlfp.lfp(:,:,:,:,1:2*m2b),[1,2,5,3,4]);
                temp=reshape(temp,s1,[],m2b,s3,s4);
                rawamylfp=permute(temp,[1,2,4,5,3]);                
                rawamylfp=rawamylfp(:,:,1:3,acellidx,blockidx);
                rawamylfp=squeeze(rawamylfp);
                amylfp=reshape(rawamylfp,s1,[]);
                amylfp=rmlinesc(amylfp,params,[],'n',[50 100]);
                [Sa,f]=mtspectrumc(amylfp,params);
                fidx= f>45 & f<55;
                mergeSa=[mergeSa; mean(Sa(fidx))];
                [~,~,~,s4,~]=size(V1conblockunitlfp.lfp);
                temp=permute(V1conblockunitlfp.lfp(:,:,:,:,1:2*m2b),[1,2,5,3,4]);
                temp=reshape(temp,s1,[],m2b,s3,s4);
                rawv1lfp=permute(temp,[1,2,4,5,3]);
                rawv1lfp=rawv1lfp(:,:,1:3,vcellidx,blockidx);
                rawv1lfp=squeeze(rawv1lfp);
                v1lfp=reshape(rawv1lfp,s1,[]);
                v1lfp=rmlinesc(v1lfp,params,[],'n',[50 100]);
                [Sv,f]=mtspectrumc(v1lfp,params);
                fidx= f>45 & f<55;
                mergeSv=[mergeSv; mean(Sv(fidx))];
                mergeC=[mergeC; peri50(noiseidx(npair))];                
            end
        end
        valididx=peri50<thr;
        mergevalididx=[mergevalididx valididx];
        mergepuffcoh=[mergepuffcoh coherencepool.puff];
        mergeneucoh=[mergeneucoh coherencepool.neu];
    end
end
figure
subplot(1,2,1)
plot(mergeC,mergeSa,'o');
subplot(1,2,2)
plot(mergeC,mergeSv,'o');
figure
scatter3(mergeSa,mergeSv,mergeC)
figure
plot(coherencepool.f,mean(mergepuffcoh(:,logical(mergevalididx)),2),'r')
hold on
plot(coherencepool.f,mean(mergeneucoh(:,logical(mergevalididx)),2),'g')
% figure
% histogram(mergeperi50)