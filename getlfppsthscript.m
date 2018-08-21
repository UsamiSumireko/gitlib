clear
folderlist={
%     'I:\160712';
%     'I:\160713';
%     'I:\160714';
%     'I:\160715';
%     'I:\160716';
%     'I:\160719';
%     'I:\160720';
%     'I:\160721';
%     'I:\160722';
%     
%     'G:\MS\160624';
%     'G:\MS\160625';
%     'G:\MS\160627';
%     'G:\MS\160628';
%     
%     'I:\170426';
%     'I:\170427';
%     'I:\170429';
%     'I:\170430';
%     'I:\170502';
%     'I:\170503';
%  
%     'I:\170516';
%     'I:\170518';
%     'I:\170519';
%     'I:\170522';
    };
% LRlist=['rrrrrrrrr' 'llll' 'rrrrrr' 'llll'];
for n=1:numel(folderlist)
    sprintf(repmat('*',1,n))
    getV1lfppsth(folderlist{n},LRlist(n),'ori',[1,1]);
end
% for n=1:numel(folderlist)
%     getV1lfppsth(folderlist{n},LRlist(n),'ori',[1,1]);
% end
% folderlist=dir('D:\data extracted');
% for iday=3:numel(folderlist)
%     sprintf(repmat('*',1,iday))
%     load(['D:\data extracted\' folderlist(iday).name '\amyconblockunitlfp.mat']);
%     load(['D:\data extracted\' folderlist(iday).name '\amyconblockunitpsth.mat']);
%     load(['D:\data extracted\' folderlist(iday).name '\amyoriblockunitlfp.mat']);
%     load(['D:\data extracted\' folderlist(iday).name '\amyoriblockunitpsth.mat']);
%     amyconblockunitlfp.spikevalid=amyconblockunitpsth.spikevalid;
%     amyoriblockunitlfp.spikevalid=amyoriblockunitpsth.spikevalid;
%     amyconblockunitlfp.stimvalid=amyconblockunitpsth.stimvalid;
%     amyoriblockunitlfp.stimvalid=amyoriblockunitpsth.stimvalid;
%     save(['D:\data extracted\' folderlist(iday).name '\amyconblockunitlfp.mat'],'amyconblockunitlfp','-v7.3');
%     save(['D:\data extracted\' folderlist(iday).name '\amyoriblockunitlfp.mat'],'amyoriblockunitlfp','-v7.3');
%     clearvars  amyconblockunitlfp amyconblockunitpsth amyoriblockunitlfp amyoriblockunitpsth
% end