% function BlinkCmpMKII(destination)
% % 
%
% % strlist={
% %     'I:\160712\ORICON008'
% %     'I:\160712\ORICON009'
% %     'I:\160712\ORICON010'
% %     'I:\160712\ORICON011'
% %     'I:\160712\ORICON012'
% %     'I:\160712\ORICON013'
% %     }
% destination='I:\160712';
% CSplusID=0:2;
% CSminID=3:14;
% Twin=[0 1000];
% fixtim=0.5;
% Fwin=[0.1 0.2];
%%
clear
eyeblinkdata=struct('eyeclosetime',[],'Nblink',[]);
destination='I:\180525';
slashidx=strfind(destination,'\');
slashidx=slashidx(end);
THR=-30000;
NORICON=0;
namelist=dir(destination);
strlist=[];
for i=3:numel(namelist)
    if namelist(i).isdir && strcmpi(namelist(i).name(1:9),'DiGratCon')
        if ~strcmp(namelist(i).name(11),'C')
            NORICON=NORICON+1;
            strlist=[strlist; {namelist(i).name}];
        end
    end
end
for n=1:numel(strlist)
    load([destination '\' strlist{n} '\ExpMonitor.mat']);
    if ~exist([destination '\' strlist{n} '\Elec140LFP.mat'],'file');    
        LoadLFPMKII([destination '\' strlist{n} '.ns3'],140,2,false);  % LoadLFP检查ns*文件从4:9，LoadLFPMKII从0:9
    end
    load([destination '\' strlist{n} '\Elec140LFP.mat']);
    StartTime=ExpMonitor.StartT;
    ValidTrial=ExpMonitor.RespCode==1;
    TTLT=ExpMonitor.TTLT;
    % if 2TTL in 10ms, repeat recording for the same TTL
    TTLinterval=TTLT(:,2:end)-TTLT(:,1:end-1);
    wrongTTL=TTLinterval<0.008;
    for itrial=1:size(ExpMonitor.TTLT,1)
        for jTTL=8:-1:1
            if wrongTTL(itrial,jTTL)
                TTLT(itrial,jTTL+1:10-1)=TTLT(itrial,jTTL+2:10);
            end
        end
    end
%     pagecycstartlabel=(TTLT(:,1:10)-repmat(TTLT(:,1),1,10))>0.008;
    tracestartlabel=(TTLT(:,1:10)-repmat(TTLT(:,1),1,10))>0.300;
    for m=1:8
        stimidx=ExpMonitor.StimID==m-1;
        stimvidx=find(stimidx & ValidTrial);
        for i=1:numel(stimvidx)
%             pagecycstartTTL=find(pagecycstartlabel(stimvidx(i)));
%             pagecycstartTTL=pagecycstartTTL(1);
            tracestartTTL=find(tracestartlabel(stimvidx(i),:));
            tracestartTTL=tracestartTTL(1);
            pagecyct=TTLT(stimvidx(i),tracestartTTL)-TTLT(stimvidx(i),1);
            tracet=TTLT(stimvidx(i),tracestartTTL+1)-TTLT(stimvidx(i),tracestartTTL);
            if isnan(tracet) || isnan(pagecyct)
               NBlink(i,m,n)=nan;
            else
            pipulsize=LFP{stimvidx(i)}(floor(2000*(-StartTime(stimvidx(i))+pagecyct))+1:...
                floor(2000*(-StartTime(stimvidx(i))+pagecyct+tracet)));
            eyeclosetime(i,m,n)=sum(pipulsize<THR);
            tracetime(i,m,n)=numel(pipulsize);
            [~,NBlink(i,m,n)]=bwlabel(pipulsize<THR, 4);
            end
        end
    end
end
eyeblinkdata.eyeclosetime=eyeclosetime;
eyeblinkdata.tracetime=tracetime;
eyeblinkdata.Nblink=NBlink;
if ~exist(['D:\data extracted\' destination(slashidx+1:end)],'dir')
    mkdir(['D:\data extracted\' destination(slashidx+1:end)])
end
save(['D:\data extracted\' destination(slashidx+1:end) '\eyeblinkdata.mat'],'eyeblinkdata','-v7.3')