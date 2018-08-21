function prepostidx=preposttestblockidx(folderlist)
% prepostidx 0=pretest 1=posttest
for iday=1:numel(folderlist)
    namelist=dir(folderlist{iday});
    CONlist=[];
    ORIlist=[];
    oriblocknum=[];
    for j=3:numel(namelist)
        if namelist(j).isdir && strcmpi(namelist(j).name(1:6),'ORICON')
           CONlist=[CONlist; {namelist(j).name}];           
        elseif namelist(j).isdir && strcmpi(namelist(j).name(1:6),'DftORI')
           ORIlist=[ORIlist; {namelist(j).name}]; 
           oriblocknum=[oriblocknum; str2double(namelist(j).name(end-2:end))];           
        end        
    end
    prepostidx{iday}=oriblocknum>str2double(CONlist{1}(end-2:end));
end
end
