clear
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
str='D:\result\170726 mergeblock coh';
sldwinlabel=0;
featurestr='Smergeblkcoherency.mat';
%%

% RP=[];
% RN=[];
% LP=[];
% LN=[];
% MP=[];
% MN=[];
% for iday=1:numel(folderlist)
%     if exist([str '\' folderlist{iday}(end-5:end) '\' featurestr],'file')
%         load([folderlist{iday} '\Oparam15.mat']);
%         load([str '\' folderlist{iday}(end-5:end) '\' featurestr])
% 
%         load([folderlist{iday} '\n2b_Svalencevalid.mat']);
%         load([folderlist{iday} '\ORIvalidid.mat']);
%         opref=Oparam15{1,4}(1:12,:);
%         ORIpref=mod(opref(ORIvalidid(1:12,:)),180);
%         Lpref=find(ORIpref>112.5 & ORIpref<157.5);
%         Rpref=find(ORIpref>22.5 & ORIpref<67.5);
%         Mpref=find(ORIpref>67.5 & ORIpref<112.5);
%         a=coherencepool.blockcellidx;
%         avc=arrayfun(@(x) x.vcellidx,a);
%         [x,y]=meshgrid(avc,Lpref);
%         Lidx=logical(sum(x==y));
%         [x,y]=meshgrid(avc,Rpref);
%         Ridx=logical(sum(x==y));
%         [x,y]=meshgrid(avc,Mpref);
%         Midx=logical(sum(x==y));
%         if sldwinlabel
%             LP=cat(3,LP,coherencepool.puff(:,:,Lidx));
%             RP=cat(3,RP,coherencepool.puff(:,:,Ridx));
%             MP=cat(3,MP,coherencepool.puff(:,:,Midx));
%             LN=cat(3,LN,coherencepool.neu(:,:,Lidx));
%             RN=cat(3,RN,coherencepool.neu(:,:,Ridx));
%             MN=cat(3,MN,coherencepool.neu(:,:,Midx));
%             t=coherencepool.t;
%         else
%             LP=[LP coherencepool.puff(:,Lidx)];
%             RP=[RP coherencepool.puff(:,Ridx)];
%             MP=[MP coherencepool.puff(:,Midx)];
%             LN=[LN coherencepool.neu(:,Lidx)];
%             RN=[RN coherencepool.neu(:,Ridx)];
%             MN=[MN coherencepool.neu(:,Midx)];
%         end
%         f=coherencepool.f;
%     end
% end
% if sldwinlabel
%     figure
%     plot_matrix(mean(LP,3)./mean(LN,3),t,f,'l');
%     title(['N=' num2str(size(LP,3))])
%     figure
%     plot_matrix(mean(RP,3)./mean(RN,3),t,f,'l');
%     title(['N=' num2str(size(RP,3))])
%     figure
%     plot_matrix(mean(MP,3)./mean(MN,3),t,f,'l');
%     title(['N=' num2str(size(MP,3))])
% else
% figure
% subplot(1,3,1)
% plot(f,mean(LP,2),'r');
% hold on
% plot(f,mean(LN,2),'g');
% title(['N=' num2str(size(LP,2))])
% subplot(1,3,2)
% plot(f,mean(RP,2),'r');
% hold on
% plot(f,mean(RN,2),'g');
% title(['N=' num2str(size(RP,2))])
% subplot(1,3,3)
% plot(f,mean(MP,2),'r');
% hold on
% plot(f,mean(MN,2),'g');
% title(['N=' num2str(size(MP,2))])
% end
%%
peakwin=231:270;
multineucnd=[4, 5, 6;
             7, 8, 9;
             10, 11, 12];
v1peakP=[];
v1peakN=[];
v1peakphip=[];
v1peakphin=[];
for iday=1:numel(folderlist)
    if strcmpi(LRlist(iday),'r')
        cndlist=[1 2 3 4 5 6];
    elseif strcmpi(LRlist(iday),'l')
        cndlist=[4 5 6 1 2 3];
    else
        error('wrong LR')
    end
    if exist([str '\' folderlist{iday}(end-5:end) '\' featurestr],'file')
        load([str '\' folderlist{iday}(end-5:end) '\' featurestr])
        load([folderlist{iday} '\V1conblockunitpsth.mat']);
        [s1,s2,s3,s4,s5]=size(V1conblockunitpsth.psth);
        load([folderlist{iday} '\V1oriblockunitpsth.mat']);
        temp=permute(V1conblockunitpsth.psth,[1 2 3 5 4]);
        V1conpuffpeak=reshape(mean(temp(peakwin,:,1:3,:,:)),[],s4);
        V1conneupeak=reshape(mean(temp(peakwin,:,4:12,:,:)),[],s4);
        temp=permute(V1oriblockunitpsth.psth,[1 2 3 5 4]);
        V1oripuffpeak=reshape(mean(temp(peakwin,:,cndlist(1:3),:,:)),[],s4);
        V1orineupeak=reshape(mean(temp(peakwin,:,cndlist(4:6),:,:)),[],s4);
        V1oripuffpeak=mean(V1oripuffpeak);
        V1orineupeak=mean(V1orineupeak);
        V1conpuffpeak=bsxfun(@minus,V1conpuffpeak,V1oripuffpeak);
        V1conneupeak=bsxfun(@minus,V1conneupeak,V1orineupeak);
        [ht,pt]=ttest2(V1conpuffpeak,V1conneupeak);
        vht=find(ht)';
        a=coherencepool.blockcellidx;
        avc=arrayfun(@(x) x.vcellidx,a);
        [x,y]=meshgrid(avc,vht);
        v1peakidx=logical(sum(x==y));        
        if sldwinlabel
            v1peakP=cat(3,v1peakP,coherencepool.puff(:,:,v1peakidx));
            v1peakN=cat(3,v1peakN,coherencepool.neu(:,:,v1peakidx));
            v1peakphip=cat(3,v1peakphip,coherencepool.phip(:,:,v1peakidx));
            v1peakphin=cat(3,v1peakphin,coherencepool.phin(:,:,v1peakidx));
        else
            v1peakP=[v1peakP coherencepool.puff(:,v1peakidx)];
            v1peakN=[v1peakN coherencepool.neu(:,v1peakidx)];
            v1peakphip=[v1peakphip coherencepool.phip(:,v1peakidx)];
            v1peakphin=[v1peakphin coherencepool.phin(:,v1peakidx)];
        end
        f=coherencepool.f;
    end
end

if sldwinlabel
    t=coherencepool.t;
    figure
    plot_matrix(mean(v1peakP,3),t,f,'l');
    figure
    plot_matrix(mean(v1peakN,3),t,f,'l');
    figure
    plot_matrix(mean(v1peakP,3)./mean(v1peakN,3),t,f,'l');
    title(['N=' num2str(size(v1peakP,3))])
else
    figure
    plot(f,mean(v1peakP,2),'r');
    hold on
    plot(f,mean(v1peakN,2),'g');
    title(['N=' num2str(size(v1peakP,2))])
end
v1peakPr=bsxfun(@(x,y) x.*cos(y),v1peakP,v1peakphip);
v1peakPv=bsxfun(@(x,y) x.*sin(y),v1peakP,v1peakphip);
v1peakPm=bsxfun(@(x,y) sqrt(x.^2+y.^2),v1peakPr,v1peakPv);
v1peakNr=bsxfun(@(x,y) x.*cos(y),v1peakN,v1peakphin);
v1peakNv=bsxfun(@(x,y) x.*sin(y),v1peakN,v1peakphin);
v1peakNm=bsxfun(@(x,y) sqrt(x.^2+y.^2),v1peakNr,v1peakNv);
if sldwinlabel
    t=coherencepool.t;
    figure
    plot_matrix(mean(v1peakPm,3),t,f,'l');
    figure
    plot_matrix(mean(v1peakNm,3),t,f,'l');
    figure
    plot_matrix(mean(v1peakPm,3)./mean(v1peakNm,3),t,f,'l');
    title(['N=' num2str(size(v1peakP,3))])
else
    figure
    plot(f,mean(v1peakPm,2),'r');
    hold on
    plot(f,mean(v1peakNm,2),'g');
    title(['N=' num2str(size(v1peakPm,2))])
end