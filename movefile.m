clear
remotedir='X:\LiZhh\amysorting';
localdir='I:\';
rdaylist=dir(remotedir);
for i=3:numel(rdaylist)
    if rdaylist(i).isdir && ~strcmp(rdaylist(i).name,'170515')
        i
        day=rdaylist(i).name;
        rnamelist=dir([remotedir '\' day ]);
        for j=3:numel(rnamelist)
            if rnamelist(j).isdir && (strcmp(rnamelist(j).name(1:6),'DftORI') || strcmp(rnamelist(j).name(1:6),'ORICON') )
                j
                matlist=dir([remotedir '\' day '\' rnamelist(j).name]);
                cd([remotedir '\' day '\' rnamelist(j).name]);
                for k=97:127
                    if exist(['Elec' num2str(k) 'Cluster.png'],'file')
                        copyfile(['Elec' num2str(k) 'Cluster.png'],[localdir '\' day '\' rnamelist(j).name '\Elec' num2str(k) 'Cluster.png'])
                    end
                    if exist(['Elec' num2str(k) 'Quality.mat'],'file')
                        copyfile(['Elec' num2str(k) 'Quality.mat'],[localdir '\' day '\' rnamelist(j).name '\Elec' num2str(k) 'Quality.mat'])
                    end
                    if exist(['Elec' num2str(k) 'Unit.mat'],'file')
                        copyfile(['Elec' num2str(k) 'Unit.mat'],[localdir '\' day '\' rnamelist(j).name '\Elec' num2str(k) 'Unit.mat'])
                    end
                end
            end
        end
    end
end