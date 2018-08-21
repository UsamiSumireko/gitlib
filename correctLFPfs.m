str='I:\180421';
namelist=dir(str);
for i=3:numel(namelist)
    if namelist(i).isdir && ~strcmp(namelist(i).name,'OFSnev')
        load([str '\' namelist(i).name '\LFPBasicMKII.mat']);
        if LFPBasic.Fs~=30000
            LFPBasic.Fs=30000;
            save([str '\' namelist(i).name '\LFPBasicMKII.mat'],'LFPBasic','-v7.3');
        end
    end
end