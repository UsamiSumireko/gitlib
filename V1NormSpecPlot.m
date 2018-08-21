clear
folderlist={
    'I:\160712';
    'I:\160713';
    'I:\160714';
    'I:\160715';
    'I:\160716';
    'I:\160719';
    'I:\160720';
    'I:\160721';
    'I:\160722';
    
    'G:\MS\160624';
    'G:\MS\160625';
    'G:\MS\160627';
    'G:\MS\160628';
    
    'I:\170426'
    'I:\170427'
    'I:\170429'
    'I:\170430'
    'I:\170502'
    'I:\170503'
    
    'I:\170518'
    'I:\170519'
    'I:\170522'
    'I:\170523'
    };
allnormS1=[];
allnormS2=[];
allconS1=[];
allconS2=[];
alloriS1=[];
alloriS2=[];
for iday=1:numel(folderlist)
    NORI=0;
    NCON=0;
    namelist=dir(folderlist{iday});
    for n=3:numel(namelist)
        if namelist(n).isdir && strcmp(namelist(n).name(1:6),'DftORI')
           NORI=NORI+1;
        end
        if namelist(n).isdir && strcmp(namelist(n).name(1:6),'ORICON')
           NCON=NCON+1;
        end
    end
    load([folderlist{iday} '\V1_spec_ori_101_500ms.mat']);
    nS=size(V1_spec.S1,2);
    ncell=nS/NORI;
    if ncell~=floor(ncell)
        error('ncell wrong')
    end
    nf=numel(V1_spec.f);
    orispecS1=reshape(V1_spec.S1,nf,ncell,NORI);
    orispecS1=squeeze(mean(orispecS1,3));
    orispecS2=reshape(V1_spec.S2,nf,ncell,NORI);
    orispecS2=squeeze(mean(orispecS2,3));
    orispecf=V1_spec.f;
    alloriS1=[alloriS1 orispecS1];
    alloriS2=[alloriS2 orispecS2];
    clearvars V1_spec
    load([folderlist{iday} '\V1_spec_con_101_500ms.mat']);
    nf=numel(V1_spec.f);
    conspecS1=reshape(V1_spec.S1,nf,ncell,NCON);
    conspecS1=squeeze(mean(conspecS1,3));
    conspecS2=reshape(V1_spec.S2,nf,ncell,NCON);
    conspecS2=squeeze(mean(conspecS2,3));
    normspecf=V1_spec.f;
    allconS1=[allconS1 conspecS1];
    allconS2=[allconS2 conspecS2];
    normspecS1=conspecS1./orispecS1;
    normspecS2=conspecS2./orispecS2;
    allnormS1=[allnormS1 normspecS1];
    allnormS2=[allnormS2 normspecS2];
end
normS1ste=std(allnormS1,0,2)/sqrt(size(allnormS1,2));
normS2ste=std(allnormS2,0,2)/sqrt(size(allnormS2,2));
figure
errorbar(normspecf,mean(allnormS1,2),normS1ste,'r')
hold on
errorbar(normspecf,mean(allnormS2,2),normS2ste,'g')

conS1ste=std(allconS1,0,2)/sqrt(size(allconS1,2));
conS2ste=std(allconS2,0,2)/sqrt(size(allconS2,2));
figure
errorbar(normspecf,mean(allconS1,2),conS1ste,'r')
hold on
errorbar(normspecf,mean(allconS2,2),conS2ste,'g')

oriS1ste=std(alloriS1,0,2)/sqrt(size(alloriS1,2));
oriS2ste=std(alloriS2,0,2)/sqrt(size(alloriS2,2));
figure
errorbar(normspecf,mean(alloriS1,2),oriS1ste,'r')
hold on
errorbar(normspecf,mean(alloriS2,2),oriS2ste,'g')

ratioste=std(allnormS1./allnormS2,0,2)/sqrt(size(alloriS2,2));
figure
errorbar(normspecf,mean(allnormS1./allnormS2,2),ratioste,'r')