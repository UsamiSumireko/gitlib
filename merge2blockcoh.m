clear
CONORI='con';
spklfp=[{'lfp'},{'lfp'}];%  (1)amy spk/lfp (2)V1 spk/lfp
savelabel=0;
sldwinlabel=0;
trialshufflelabel=0;
bwin=401:800;
sldwin=[0.2 0.005];
LooseCriterionLabel=0; % Monkey MS avaliable amy cell too few, loose the cell selection criterion to p<0.10 at ttest
V1selectlabel=0; % V1 elec stable in 1 day, pretest consistant with post test
v1secpklabel=-1; % -1 all cell 0 cell without second peak 1 cell with second peak
meansubstract='intertrial';% intertrial intratial both
% meansubstract='';
%%%%%%%%%
if V1selectlabel
    load('D:\result\171122 V1cell select\V1sigcell_ttest.mat');
end
%%%%%%%%
folderlist={
%     'D:\data extracted\160712';
%     'D:\data extracted\160713';
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
allpuff=[];
allneu=[];
allphip=[];
allphin=[];
str='D:\result\180606 granger and coherency ultra\resting state coh';
% str='D:\result\171122 V1cell select\coherency';
% str='D:\result\170828 mergeblock lfpspkcoh\';
multineucnd=[4, 5, 6;
             7, 8, 9;
             10, 11, 12];
LRlist=['rrrrrrrrr' 'llll' 'rrrrrr' 'llll'];         

params.Fs = 1000; % sampling frequency
params.tapers = [2 3]; % taper parameters
params.trialave = 1;
params.fpass = [4 100];
params.pad = 2;
params.err =[2 0.05];
npairsession=32;
nrml=2000;
tic
for iday=1:numel(folderlist)
    coherencepool=struct('puff',[],'neu',[]);
    load([folderlist{iday} '\secpeakcellidx_both.mat']);
    %% load data and amy stim valence select
    if strcmpi(CONORI,'con') 
        if strcmpi(spklfp{1},'lfp')
            load([folderlist{iday} '\amyconblockunitlfp.mat']);
            load([folderlist{iday} '\amyconblockunitpsth.mat']);
            amyblockunitdata.data=amyconblockunitlfp.lfp;
            amyblockunitdata.spikevalid=amyconblockunitpsth.spikevalid;
        elseif strcmpi(spklfp{1},'spk')
            load([folderlist{iday} '\amyconblockunitpsth.mat']);
            amyblockunitdata.data=amyconblockunitpsth.psth;
            amyblockunitdata.spikevalid=amyconblockunitpsth.spikevalid;
        end
        if strcmpi(spklfp{2},'lfp')
            load([folderlist{iday} '\V1conblockunitlfp.mat']);
            if V1selectlabel
                V1selectidx=logical(V1sigcell_ttest{iday});
            else
                V1selectidx=true(size(V1conblockunitlfp.lfp,4),1);
            end
            if v1secpklabel<0
                secpklabel=true(size(V1conblockunitlfp.lfp,4),1);
            elseif v1secpklabel==0
                secpklabel=~secpeakcellidx;
            elseif v1secpklabel==1
                secpklabel=secpeakcellidx;
            else
                warning('secpklabel may be wrong')
            end
            V1blockunitdata.data=V1conblockunitlfp.lfp(:,:,:,V1selectidx & secpklabel,:);
%             V1blockunitdata.spikevalid=V1conblockunitlfp.spikevalid;
        elseif strcmpi(spklfp{2},'spk')
            load([folderlist{iday} '\V1conblockunitpsth.mat']);
            if V1selectlabel
                V1selectidx=logical(V1sigcell_ttest{iday});
            else
                V1selectidx=true(size(V1conblockunitpsth.psth,4),1);
            end
            if v1secpklabel<0
                secpklabel=true(size(V1conblockunitpsth.psth,4),1);
            elseif v1secpklabel==0
                secpklabel=~secpeakcellidx;
            elseif v1secpklabel==1
                secpklabel=secpeakcellidx;
            else
                warning('secpklabel may be wrong')
            end
            V1blockunitdata.data=V1conblockunitpsth.psth(:,:,:,V1selectidx & secpklabel,:);
%             V1blockunitdata.spikevalid=V1conblockunitpsth.spikevalid;
        end
        
        clearvars V1conblockunitlfp V1conblockunitlfp 
        [s1,s2,s3,s4,s5]=size(amyblockunitdata.data);
        n2b=floor(s5/2);     
        [ss1,ss2,ss3,ss4,ss5]=size(V1blockunitdata.data);
        amyrawdata=permute(amyblockunitdata.data(:,:,:,:,1:2*n2b),[1 2 5 3 4]);
        amyrawdata=reshape(amyrawdata,s1,[],n2b,s3,s4);
        amyrawdata=permute(amyrawdata,[1 2 4 5 3]);
        v1rawdata=permute(V1blockunitdata.data(:,:,:,:,1:2*n2b),[1 2 5 3 4]);
        v1rawdata=reshape(v1rawdata,ss1,[],n2b,ss3,ss4);
        v1rawdata=permute(v1rawdata,[1 2 4 5 3]);
        if strcmpi(spklfp{1},'spk')
            amyrawspk=amyrawdata;
        elseif strcmpi(spklfp{1},'lfp')
            amyrawspk=permute(amyconblockunitpsth.psth(:,:,:,:,1:2*n2b),[1 2 5 3 4]);
%             amyrawspk=permute(amyconblockunitpsth.psth(201:end,:,:,:,1:2*n2b),[1 2 5 3 4]);
            amyrawspk=reshape(amyrawspk,s1,[],n2b,s3,s4);
            amyrawspk=permute(amyrawspk,[1 2 4 5 3]);
        end
        [~,sss1,sss2,sss3,sss4]=size(amyrawspk);        
        amyblfr=squeeze(mean(amyrawspk(201:400,:,:,:,:)));
        amyblfr=reshape(amyblfr,sss1*sss2,[]);        
        amystimfr=squeeze(mean(amyrawspk(401:800,:,:,:,:)));
        amypufffr=amystimfr(:,1:3,:,:);
        amyneufr=amystimfr(:,4:end,:,:);
        amystimfr=reshape(amystimfr,sss1*sss2,[]);
        [t1h,t1p]=ttest(amyblfr,amystimfr);
        if LooseCriterionLabel
        amystimvalid=reshape(t1p<0.1,sss3,sss4);
        else
        amystimvalid=reshape(t1h==1,sss3,sss4);
        end
        
        amypufffr=reshape(amypufffr,[],sss3*sss4);
        amyneufr=reshape(amyneufr,[],sss3*sss4);
        [t2h,t2p]=ttest2(amypufffr,amyneufr,'Tail','left');
        if LooseCriterionLabel
        valencevalid=reshape(t2p<0.1,sss3,sss4);
        else
        valencevalid=reshape(t2h==1,sss3,sss4);
        end
        
%         amystimvalid=amyblockunitlfp.stimvalid(:,1:2:n2b*2-1) & amyblockunitlfp.stimvalid(:,2:2:n2b*2);
        amyspkvalid=amyblockunitdata.spikevalid(:,1:2:n2b*2-1) & amyblockunitdata.spikevalid(:,2:2:n2b*2);
        amycellidx=amystimvalid & amyspkvalid & valencevalid;
%         amycellidx=amyspkvalid & valencevalid;
        
    elseif strcmpi(CONORI,'ori')
        lr=LRlist(iday);
        if strcmpi(lr,'r')
            cndlist=[1 2 3 4 5 6];
        elseif strcmpi(lr,'l')
            cndlist=[4 5 6 1 2 3];
        end
        if strcmpi(spklfp{1},'lfp')
            load([folderlist{iday} '\amyoriblockunitlfp.mat']);
            load([folderlist{iday} '\amyoriblockunitpsth.mat']);
            amyblockunitdata.data=amyoriblockunitlfp.lfp;
            amyblockunitdata.spikevalid=amyoriblockunitpsth.spikevalid;
        elseif strcmpi(spklfp{1},'spk')
            load([folderlist{iday} '\amyoriblockunitpsth.mat']);
            amyblockunitdata.data=amyoriblockunitpsth.psth;
            amyblockunitdata.spikevalid=amyoriblockunitpsth.spikevalid;
        end        
        if strcmpi(spklfp{2},'lfp')
            load([folderlist{iday} '\V1oriblockunitlfp.mat']);
            if V1selectlabel
                V1selectidx=logical(V1sigcell_ttest{iday});
            else
                V1selectidx=true(size(V1oriblockunitlfp.lfp,4),1);
            end
            if v1secpklabel<0
                secpklabel=true(size(V1oriblockunitlfp.lfp,4),1);
            elseif v1secpklabel==0
                secpklabel=~secpeakcellidx;
            elseif v1secpklabel==1
                secpklabel=secpeakcellidx;
            else
                warning('secpklabel may be wrong')
            end
            V1blockunitdata.data=V1oriblockunitlfp.lfp(:,:,:,V1selectidx & secpklabel,:);
%             V1blockunitdata.spikevalid=V1oriblockunitlfp.spikevalid;
        elseif strcmpi(spklfp{2},'spk')
            load([folderlist{iday} '\V1oriblockunitpsth.mat']);
            if V1selectlabel
                V1selectidx=logical(V1sigcell_ttest{iday});
            else
                V1selectidx=true(size(V1oriblockunitpsth.psth,4),1);
            end
            if v1secpklabel<0
                secpklabel=true(size(V1oriblockunitpsth.psth,4),1);
            elseif v1secpklabel==0
                secpklabel=~secpeakcellidx;
            elseif v1secpklabel==1
                secpklabel=secpeakcellidx;
            else
                warning('secpklabel may be wrong')
            end
            V1blockunitdata.data=V1oriblockunitpsth.psth(:,:,:,V1selectidx & secpklabel,:);
%             V1blockunitdata.spikevalid=V1oriblockunitpsth.spikevalid;
        end                        
        clearvars V1oriblockunitlfp V1oriblockunitlfp
%         if size(amyblockunitlfp.lfp,5)>4
%             amyblockunitlfp.lfp=amyblockunitlfp.lfp(:,:,:,:,2:3);
%             amyblockunitlfp.stimvalid=amyblockunitlfp.stimvalid(:,2:3);
%             amyblockunitlfp.spikevalid=amyblockunitlfp.spikevalid(:,2:3);
%             V1blockunitlfp.lfp=V1blockunitlfp.lfp(:,:,:,:,2:3);
%         else
%             amyblockunitlfp.lfp=amyblockunitlfp.lfp(:,:,:,:,1:2);
%             amyblockunitlfp.stimvalid=amyblockunitlfp.stimvalid(:,1:2);
%             amyblockunitlfp.spikevalid=amyblockunitlfp.spikevalid(:,1:2);
%             V1blockunitlfp.lfp=V1blockunitlfp.lfp(:,:,:,:,1:2);
%         end
        [s1,s2,s3,s4,s5]=size(amyblockunitdata.data);
        amyrawdata=amyblockunitdata.data;
        v1rawdata=V1blockunitdata.data;
        
        amyrawspk=amyoriblockunitpsth.psth;
        
        [~,sss1,sss2,sss3,sss4]=size(amyrawdata);        
        amyblfr=squeeze(mean(amyrawspk(201:400,:,:,:,:)));
        amyblfr=reshape(amyblfr,sss1*sss2,[]);
        amystimfr=squeeze(mean(amyrawspk(401:800,:,:,:,:)));
        amystimfr=reshape(amystimfr,sss1*sss2,[]);
        [t1h,t1p]=ttest(amyblfr,amystimfr);
        amystimvalid=reshape(t1h==1,sss3,sss4);
        valencevalid=ones(24,s5);
%         amystimvalid=amyblockunitlfp.stimvalid;
        amyspkvalid=amyblockunitdata.spikevalid;
%         amycellidx=amystimvalid & amyspkvalid & valencevalid;
        amycellidx= amyspkvalid & valencevalid;
    end
    %% prepare data format and remove line noise
    acount=0;
    vcount=0;
    amydata=[];
    V1data=[];
    amyidx=struct('block',[],'cell',[]);
    V1idx=struct('block',[],'cell',[]);
    for nblock=1:size(amycellidx,2)
        amyblockcellidx=find(amycellidx(:,nblock));
        validblockcellidx=[];
        if ~isempty(amyblockcellidx)
            amylfptemp=amyrawdata(bwin,:,:,amyblockcellidx,nblock);
            nanidx=isnan(amylfptemp);
            nanvalid=0==squeeze(sum(sum(sum(nanidx,1),2),3));
            idx65535=65535==amylfptemp;
            tbinvalid=0==squeeze(sum(sum(sum(idx65535,1),2),3));
            validblockcellidx=find(nanvalid & tbinvalid);
            if ~all(nanvalid & tbinvalid)
                sprintf([num2str(iday) ' ' num2str(nblock) ' '])
            end
        end
        for acell=validblockcellidx'
            acount=acount+1;
            if strcmpi(CONORI,'con')
                data1=amylfptemp(:,:,1:3,acell);
                data2=amylfptemp(:,:,multineucnd(:),acell);
                switch meansubstract
                    case 'intertrial'
                        mu=mean(data1,2);
                        data1=data1-repmat(mu,1,size(data1,2),1);
                        mu=mean(data2,2);
                        data2=data2-repmat(mu,1,size(data2,2),1);
                    case 'intratrial'
                        mu=mean(data1,1);
                        data1=data1-repmat(mu,size(data1,1),1,1);
                        mu=mean(data2,1);
                        data2=data2-repmat(mu,size(data2,1),1,1);
                    case 'both'
                        mu=mean(data1,1);
                        data1=data1-repmat(mu,size(data1,1),1,1);
                        mu=mean(data1,2);
                        data1=data1-repmat(mu,1,size(data1,2),1);
                        mu=mean(data2,1);
                        data2=data2-repmat(mu,size(data2,1),1,1);
                        mu=mean(data2,2);
                        data2=data2-repmat(mu,1,size(data2,2),1);
                    otherwise
                end               
                data1=reshape(data1,size(data1,1),[]);                
                data2=reshape(data2,size(data2,1),[]);             
            elseif strcmpi(CONORI,'ori')
                data1=amylfptemp(:,:,cndlist(1:3),acell);                
                data2=amylfptemp(:,:,cndlist(4:6),acell);
                switch meansubstract
                    case 'intertrial'
                        mu=mean(data1,2);
                        data1=data1-repmat(mu,1,size(data1,2),1);
                        mu=mean(data2,2);
                        data2=data2-repmat(mu,1,size(data2,2),1);
                    case 'intratrial'
                        mu=mean(data1,1);
                        data1=data1-repmat(mu,size(data1,1),1,1);
                        mu=mean(data2,1);
                        data2=data2-repmat(mu,size(data2,1),1,1);
                    case 'both'
                        mu=mean(data1,1);
                        data1=data1-repmat(mu,size(data1,1),1,1);
                        mu=mean(data1,2);
                        data1=data1-repmat(mu,1,size(data1,2),1);
                        mu=mean(data2,1);
                        data2=data2-repmat(mu,size(data2,1),1,1);
                        mu=mean(data2,2);
                        data2=data2-repmat(mu,1,size(data2,2),1);
                    otherwise
                end               
                data1=reshape(data1,size(data1,1),[]);                
                data2=reshape(data2,size(data2,1),[]);
            end
            if strcmpi(spklfp{1},'lfp')
                amydata(acount).puff=rmlinesc(data1,params,[],'n',[50 100]);
                amydata(acount).neu=rmlinesc(data2,params,[],'n',[50 100]);
            elseif strcmpi(spklfp{1},'spk')
                amydata(acount).puff=data1;
                amydata(acount).neu=data2;
            end
            amyidx.block(acount)=nblock;
            amyidx.cell(acount)=amyblockcellidx(acell);
        end
        
        v1lfptemp=v1rawdata(bwin,:,:,:,nblock);
        nanidx=isnan(v1lfptemp);
        nanvalid=0==squeeze(sum(sum(sum(nanidx,1),2),3));
        idx65535=65535==v1lfptemp;
        tbinvalid=0==squeeze(sum(sum(sum(idx65535,1),2),3));
        v1blockcellidx=find(nanvalid & tbinvalid);
        for vcell=v1blockcellidx'                    
            vcount=vcount+1;
            if strcmpi(CONORI,'con')
                data1=v1lfptemp(:,:,1:3,vcell);
                data2=v1lfptemp(:,:,multineucnd(:),vcell);
                switch meansubstract
                    case 'intertrial'
                        mu=mean(data1,2);
                        data1=data1-repmat(mu,1,size(data1,2),1);
                        mu=mean(data2,2);
                        data2=data2-repmat(mu,1,size(data2,2),1);
                    case 'intratrial'
                        mu=mean(data1,1);
                        data1=data1-repmat(mu,size(data1,1),1,1);
                        mu=mean(data2,1);
                        data2=data2-repmat(mu,size(data2,1),1,1);
                    case 'both'
                        mu=mean(data1,1);
                        data1=data1-repmat(mu,size(data1,1),1,1);
                        mu=mean(data1,2);
                        data1=data1-repmat(mu,1,size(data1,2),1);
                        mu=mean(data2,1);
                        data2=data2-repmat(mu,size(data2,1),1,1);
                        mu=mean(data2,2);
                        data2=data2-repmat(mu,1,size(data2,2),1);
                    otherwise
                end               
                data1=reshape(data1,size(data1,1),[]);                
                data2=reshape(data2,size(data2,1),[]);
            elseif strcmpi(CONORI,'ori')
                data1=squeeze(v1lfptemp(:,:,cndlist(1:3),vcell));
                data2=squeeze(v1lfptemp(:,:,cndlist(4:6),vcell));
                switch meansubstract
                    case 'intertrial'
                        mu=mean(data1,2);
                        data1=data1-repmat(mu,1,size(data1,2),1);
                        mu=mean(data2,2);
                        data2=data2-repmat(mu,1,size(data2,2),1);
                    case 'intratrial'
                        mu=mean(data1,1);
                        data1=data1-repmat(mu,size(data1,1),1,1);
                        mu=mean(data2,1);
                        data2=data2-repmat(mu,size(data2,1),1,1);
                    case 'both'
                        mu=mean(data1,1);
                        data1=data1-repmat(mu,size(data1,1),1,1);
                        mu=mean(data1,2);
                        data1=data1-repmat(mu,1,size(data1,2),1);
                        mu=mean(data2,1);
                        data2=data2-repmat(mu,size(data2,1),1,1);
                        mu=mean(data2,2);
                        data2=data2-repmat(mu,1,size(data2,2),1);
                    otherwise
                end               
                data1=reshape(data1,size(data1,1),[]);                    
                data2=reshape(data2,size(data2,1),[]);    
            end
            if strcmpi(spklfp{2},'lfp')
                V1data(vcount).puff=rmlinesc(data1,params,[],'n',[50 100]);
                V1data(vcount).neu=rmlinesc(data2,params,[],'n',[50 100]);
            elseif strcmpi(spklfp{2},'spk')
                V1data(vcount).puff=data1;
                V1data(vcount).neu=data2;
            end
            V1idx.block(vcount)=nblock;
            V1idx.cell(vcount)=vcell;
        end
    end

    mixcount=0;
    amypuff=[];
    amyneu=[];
    V1puff=[];
    V1neu=[];
    blockcellidx=[];
    trialidx=[];
    validblock=unique(amyidx.block);
    for nblock=1:numel(validblock)
        ablockidx=find(amyidx.block==validblock(nblock));
        vblockidx=find(V1idx.block==validblock(nblock));
        for acell=1:numel(ablockidx)
            for vcell=1:numel(vblockidx)
                mixcount=mixcount+1;
                nanidx1=isnan(amydata(ablockidx(acell)).puff) | isnan(V1data(vblockidx(vcell)).puff);
                nanidx2=isnan(amydata(ablockidx(acell)).neu) | isnan(V1data(vblockidx(vcell)).neu);
                nidx1=sum(0==sum(nanidx1));
%                 assert(nidx1==30);
                nidx2=sum(0==sum(nanidx2));
                nannumidx2=find(0==sum(nanidx2));
                selecidx=nannumidx2(randperm(nidx2,nidx1));
                amypuff(:,:,mixcount)=amydata(ablockidx(acell)).puff;
                V1puff(:,:,mixcount)=V1data(vblockidx(vcell)).puff;
                amyneu(:,:,mixcount)=amydata(ablockidx(acell)).neu(:,selecidx);
                V1neu(:,:,mixcount)=V1data(vblockidx(vcell)).neu(:,selecidx);
                blockcellidx(mixcount).blockidx=nblock;
                blockcellidx(mixcount).acellidx=amyidx.cell(ablockidx(acell));
                blockcellidx(mixcount).vcellidx=V1idx.cell(vblockidx(vcell));
                trialidx(mixcount).idx1=(0==sum(nanidx1));
                trialidx(mixcount).idx2=selecidx;
            end
        end
    end
    %% compute coherency
    if ~isempty(amypuff)
        savestr=[CONORI '_' spklfp{1} '_' spklfp{2} '_'];
        tamypuff=reshape(amypuff,size(amypuff,1),[]);
        tamyneu=reshape(amyneu,size(amyneu,1),[]);
        tV1puff=reshape(V1puff,size(V1puff,1),[]);
        tV1neu=reshape(V1neu,size(V1neu,1),[]);
%         for nloop=1:ceil(size(tamypuff,2)/nrml)
%             if nloop*nrml<size(tamypuff,2)
%                 loopidx=(nloop-1)*nrml+1:nloop*nrml;
%                 camypuff(:,loopidx)=rmlinesc(tamypuff(:,loopidx),params,[],'n',[50 100]);
%                 camyneu(:,loopidx)=rmlinesc(tamyneu(:,loopidx),params,[],'n',[50 100]);
%                 cV1puff(:,loopidx)=rmlinesc(tV1puff(:,loopidx),params,[],'n',[50 100]);
%                 cV1neu(:,loopidx)=rmlinesc(tV1neu(:,loopidx),params,[],'n',[50 100]);
%             else
%                 loopidx=(nloop-1)*nrml+1:size(tamypuff,2);
%                 camypuff(:,loopidx)=rmlinesc(tamypuff(:,loopidx),params,[],'n',[50 100]);
%                 camyneu(:,loopidx)=rmlinesc(tamyneu(:,loopidx),params,[],'n',[50 100]);
%                 cV1puff(:,loopidx)=rmlinesc(tV1puff(:,loopidx),params,[],'n',[50 100]);
%                 cV1neu(:,loopidx)=rmlinesc(tV1neu(:,loopidx),params,[],'n',[50 100]);
%             end        
%         end
        camypuff=reshape(tamypuff,size(amypuff));
        camyneu=reshape(tamyneu,size(amyneu));
        cV1puff=reshape(tV1puff,size(amypuff));
        cV1neu=reshape(tV1neu,size(amyneu));        
        
        if ~trialshufflelabel % not shuffle
            if sldwinlabel 
                if sldwin(1)<=0.256
                    params.pad = 3;
                end
                [Cpuff,phip,~,~,~,t,fp]=cohgramc_gpuMKII(camypuff,cV1puff,sldwin,params,64);
                [Cneu,phin,~,~,~,~,fn]=cohgramc_gpuMKII(camyneu,cV1neu,sldwin,params,64);
                t_puff=gather(Cpuff);
                t_neu=gather(Cneu);
                t_phip=gather(phip);
                t_phin=gather(phin);
                t_t=gather(t);
                t_f=gather(fp);
                allphip=cat(3,allphip,t_phip);
                allphin=cat(3,allphin,t_phin);
                allpuff=cat(3,allpuff,t_puff);
                allneu=cat(3,allneu,t_neu);
            else % not sliding window
                gamypuff=gpuArray(camypuff);
                gamyneu=gpuArray(camyneu);
                gV1puff=gpuArray(cV1puff);
                gV1neu=gpuArray(cV1neu);
                t_puff=[];
                t_neu=[];
                t_phip=[];
                t_phin=[];
                for n=1:ceil(mixcount/npairsession)
                    if n*npairsession>mixcount
                        pairidx=1+npairsession*(n-1):mixcount;
                    else
                        pairidx=1+npairsession*(n-1):npairsession*n;
                    end
                    [Cpuff,phip,~,~,~,fp]=coherencyc_gpu(gamypuff(:,:,pairidx),gV1puff(:,:,pairidx),params);
                    [Cneu,phin,~,~,~,fn]=coherencyc_gpu(gamyneu(:,:,pairidx),gV1neu(:,:,pairidx),params);
                    t_puff(:,pairidx)=gather(Cpuff);
                    t_neu(:,pairidx)=gather(Cneu);
                    t_phip(:,pairidx)=gather(phip);
                    t_phin(:,pairidx)=gather(phin);
                    t_f=gather(fp);
                end
                allpuff=cat(2,allpuff,t_puff);
                allneu=cat(2,allneu,t_neu);
                allphip=[allphip t_phip];
                allphin=[allphin t_phin];
            end
        else % trial shuffle
            savestr=[savestr 'trialshuffle_'];
            if sldwinlabel
                if sldwin(1)<=0.256
                    params.pad = 3;
                end
                t_puff=[];
                t_neu=[];
                t_phip=[];
                t_phin=[];
                for mshuffle=1:size(camypuff,2)-1
                    sprintf(repmat('*',1,mshuffle))
                    shuffleidx=[mshuffle+1:size(camypuff,2) 1:mshuffle];
                    [Cpuff,phip,~,~,~,t,fp]=cohgramc_gpuMKII(camypuff,cV1puff(:,shuffleidx,:),sldwin,params);
                    [Cneu,phin,~,~,~,~,fn]=cohgramc_gpuMKII(camyneu,cV1neu(:,shuffleidx,:),sldwin,params);
                    t_puff(:,:,:,mshuffle)=gather(Cpuff);
                    t_neu(:,:,:,mshuffle)=gather(Cneu);
                    t_phip(:,:,:,mshuffle)=gather(phip);
                    t_phin(:,:,:,mshuffle)=gather(phin);
                end
                t_puff=permute(t_puff,[1 2 4 3]);
                t_neu=permute(t_neu,[1 2 4 3]);
                t_phip=permute(t_phip,[1 2 4 3]);
                t_phin=permute(t_phin,[1 2 4 3]);
                t_t=gather(t);
                t_f=gather(fp);
                allphip=cat(4,allphip,t_phip);
                allphin=cat(4,allphin,t_phin);
                allpuff=cat(4,allpuff,t_puff);
                allneu=cat(4,allneu,t_neu);
            else
                gamypuff=gpuArray(camypuff);
                gamyneu=gpuArray(camyneu);
                gV1puff=gpuArray(cV1puff);
                gV1neu=gpuArray(cV1neu);
                t_puff=[];
                t_neu=[];
                t_phip=[];
                t_phin=[];
                for mshuffle=1:size(camypuff,2)-1
                    sprintf(repmat('*',1,mshuffle))
                    shuffleidx=[mshuffle+1:size(camypuff,2) 1:mshuffle];
                    for n=1:ceil(mixcount/npairsession)
                        if n*npairsession>mixcount
                            pairidx=1+npairsession*(n-1):mixcount;
                        else
                            pairidx=1+npairsession*(n-1):npairsession*n;
                        end
                        [Cpuff,phip,~,~,~,fp]=coherencyc_gpu(gamypuff(:,:,pairidx),gV1puff(:,shuffleidx,pairidx),params);
                        [Cneu,phin,~,~,~,fn]=coherencyc_gpu(gamyneu(:,:,pairidx),gV1neu(:,shuffleidx,pairidx),params);
                        t_puff(:,pairidx,mshuffle)=gather(Cpuff);
                        t_neu(:,pairidx,mshuffle)=gather(Cneu);
                        t_phip(:,pairidx,mshuffle)=gather(phip);
                        t_phin(:,pairidx,mshuffle)=gather(phin);
                    end
                end
                t_puff=permute(t_puff,[1 3 2]);
                t_neu=permute(t_neu,[1 3 2]);
                t_phip=permute(t_phip,[1 3 2]);
                t_phin=permute(t_phin,[1 3 2]);
                t_f=gather(fp);
                allpuff=cat(3,allpuff,t_puff);
                allneu=cat(3,allneu,t_neu);
                allphip=cat(3,allphip,t_phip);
                allphin=cat(3,allphin,t_phin);                
            end            
        end        
        coherencepool.puff=t_puff;
        coherencepool.neu=t_neu;
        coherencepool.phip=t_phip;
        coherencepool.phin=t_phin;
        coherencepool.f=t_f;
        coherencepool.blockcellidx=blockcellidx;
        coherencepool.trialidx=trialidx;
        if sldwinlabel
            savestr=[savestr 'sldwin_'];
            coherencepool.t=t_t;
        end
        if ~exist([str '\' folderlist{iday}(end-5:end)],'dir')
            mkdir([str '\' folderlist{iday}(end-5:end)])
        end
        if savelabel
            savestr=[savestr num2str(bwin(1)-400) '-' num2str(bwin(2)-400) 'ms_'];
            save([str '\' folderlist{iday}(end-5:end) '\' savestr 'coherencepool_g.mat'],'coherencepool','-v7.3');
        end
        clear gV1puff gV1neu gamypuff gamyneu
    end
    toc
    sprintf(['Day_' num2str(iday) 'acell = ' num2str(mixcount/numel(v1blockcellidx))])
end
if sldwinlabel
    if trialshufflelabel
        figure
        plot_matrix(mean(mean(allpuff,4),3),t_t,t_f,'n');
        title('puff ori')
        figure
        plot_matrix(mean(mean(allneu,4),3),t_t,t_f,'n');
        title('neu ori')
        figure
        plot_matrix(mean(mean(allpuff,4),3)./mean(mean(allneu,4),3),t_t,t_f,'n');
    else
        figure
        plot_matrix(mean(allpuff,3),t_t,t_f,'n');
        title('puff ori')
        figure
        plot_matrix(mean(allneu,3),t_t,t_f,'n');
        title('neu ori')
        figure
        plot_matrix(mean(allpuff,3)./mean(allneu,3),t_t,t_f,'n');
    end
else
    if trialshufflelabel
        figure
        plot(t_f,mean(mean(allpuff,3),2),'r');
        hold on
        plot(t_f,mean(mean(allneu,3),2),'g');
    else
    figure
    plot(t_f,mean(allpuff,2),'r');
    hold on
    plot(t_f,mean(allneu,2),'g');
    end
end

% figure
% puffste=std(allpuff,0,2)/sqrt(size(allpuff,2));
% neuste=std(allneu,0,2)/sqrt(size(allneu,2));
% xv=t_f;
% yv=mean(allpuff,2);
% patch([xv,fliplr(xv)],[yv-puffste;flipud(yv+puffste)],[1 0.8 0.8]);
% hold on
% plot(xv,yv,'r');
% yv=mean(allneu,2);
% patch([xv,fliplr(xv)],[yv-neuste;flipud(yv+neuste)],[0.8 1 0.8]);
% hold on
% plot(xv,yv,'g');