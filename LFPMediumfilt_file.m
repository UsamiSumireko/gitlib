folderlist={
    'D:\data extracted\160712';
%     'D:\data extracted\160713';
%     'D:\data extracted\160714';
%     'D:\data extracted\160715';
%     'D:\data extracted\160716';
%     'D:\data extracted\160719';
%     'D:\data extracted\160720';
%     'D:\data extracted\160721';
%     'D:\data extracted\160722';
%     'D:\data extracted\160624';
%     'D:\data extracted\160625';
%     'D:\data extracted\160627';
%     'D:\data extracted\160628';
%     'D:\data extracted\170426';
%     'D:\data extracted\170427';
%     'D:\data extracted\170429';
%     'D:\data extracted\170430';
%     'D:\data extracted\170502';
%     'D:\data extracted\170503';
%     'D:\data extracted\170516';
%     'D:\data extracted\170518';
%     'D:\data extracted\170519';
%     'D:\data extracted\170522';
    };
load('D:\data extracted\LFPMedium.mat');
for iday=1:numel(folderlist)
    if ~exist([folderlist{iday} '\amyconblockunitlfp_mfilted.mat'],'file')
        load([folderlist{iday} '\amyconblockunitlfp.mat'])
        [s1,s2,s3,s4,s5]=size(amyconblockunitlfp.lfp);
        lfptemp=amyconblockunitlfp.lfp;
        lfptemp=reshape(lfptemp,s1,s2*s3*s4*s5);
        lfptemp=filtfilt(SOS,G,lfptemp);
        amyconblockunitlfp.lfp=reshape(lfptemp,s1,s2,s3,s4,s5);
        save([folderlist{iday} '\amyconblockunitlfp_mfilted.mat'],'amyconblockunitlfp','-v7.3');
    end
    if ~exist([folderlist{iday} '\V1conblockunitlfp_mfilted.mat'],'file')
        load([folderlist{iday} '\V1conblockunitlfp.mat'])
        [s1,s2,s3,s4,s5]=size(V1conblockunitlfp.lfp);
        lfptemp=V1conblockunitlfp.lfp;
        lfptemp=reshape(lfptemp,s1,s2*s3*s4*s5);
        lfptemp=filtfilt(SOS,G,lfptemp);
        V1conblockunitlfp.lfp=reshape(lfptemp,s1,s2,s3,s4,s5);
        save([folderlist{iday} '\V1conblockunitlfp_mfilted.mat'],'V1conblockunitlfp','-v7.3');
    end    
end

