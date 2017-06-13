function MergeV1LFP(destination,LR,win,select,CONORIlabel)
tic
% destination='I:\160720';
% LR='R';
% win=[-0.3 1.1];
% CONORIlabel=[1 1]  [CON ORI]  
oriwin=[-0.3 0.5];
%%
cmppath = 'G:\SN 4566-001303+Probe MK2.cmp';
cmpinfo=LoadCmp(cmppath,1,0);
elec=cmpinfo{1,1}.RealElec;
load([destination '\XYvalidid']);
load([destination '\ORIvalidid']);
selecti=lower(select);
switch selecti
    case 'xy'
        valididx=XYvalidid;
    case 'ori'
        valididx=ORIvalidid;
    case 'all'
        valididx=ones(12,8);
    otherwise
        msgbox('wrong select idx')
end
lr=lower(LR);
if strcmp(lr,'l')
    CNDlist=[4 5 6 1 2 3 1 2 3 1 2 3];
elseif strcmp(lr,'r')
    CNDlist=[1 2 3 4 5 6 4 5 6 4 5 6];
else
    error('LR input is wrong, L or R');
end
%%
namelist=dir(destination);
NORICON=0;
ORICONlist=[];
NDftORI=0;
DftORIlist=[];
for i=3:numel(namelist)
    if strcmp(namelist(i).name(1:6),'ORICON') && namelist(i).isdir 
        NORICON=NORICON+1;
        ORICONlist=[ORICONlist; {namelist(i).name}];
    elseif strcmp(namelist(i).name(1:6),'DftORI') && namelist(i).isdir 
        NDftORI=NDftORI+1;
        DftORIlist=[DftORIlist; {namelist(i).name}];
    end
end

%%
if CONORIlabel(1)
    V1CONblockunitLFP=cell(numel(ORICONlist),1);
    for i=1:NORICON
        V1CONblockunitLFP{i}=cell(12,8);
        i,
        for j=1:12
            for k=1:8
                sprintf(repmat('.',1,(j-1)*8+k))
                if valididx(j,k)
                    V1CONblockunitLFP{i}{j,k}=cell(6,1);
                    lfpinfo=TruncateLFP([destination '\' ORICONlist{i} '.ns5'],elec(j,k),win,'MKII');
                    [splitlfp]=SplitInfo(lfpinfo, [1 0 1 0]);
                    for m=1:numel(splitlfp)
                        dslfp=downsample(splitlfp{m}{1}.LFP',30);
                        V1CONblockunitLFP{i}{j,k}{CNDlist(mod((m-1),3)+1+3*(m>3))}=[V1CONblockunitLFP{i}{j,k}{CNDlist(mod((m-1),3)+1+3*(m>3))} dslfp];
                    end
                end
            end
        end
    end
    save([destination '\V1CONblockunitLFP.mat'],'V1CONblockunitLFP','-v7.3');
end
sprintf('CON done')
%%
if CONORIlabel(2)
    V1ORIblockunitLFP=cell(NDftORI,1);
    for i=1:NDftORI
        V1ORIblockunitLFP{i}=cell(12,8);
        i,
        for j=1:12
            for k=1:8
                sprintf(repmat('.',1,(j-1)*8+k))
                if valididx(j,k)
                    V1ORIblockunitLFP{i}{j,k}=cell(6,1);
                    lfpinfo=TruncateLFP([destination '\' DftORIlist{i} '.ns5'],elec(j,k),oriwin,'MKII');
                    [splitlfp]=SplitInfo(lfpinfo, [1 0 1 0]);
                    for m=2:2:numel(splitlfp)
                        dslfp=downsample(splitlfp{m}{1}.LFP',30);
                        V1ORIblockunitLFP{i}{j,k}{mod(m-1,12)/2+1/2}=[V1ORIblockunitLFP{i}{j,k}{mod(m-1,12)/2+1/2} dslfp];
                    end
                end
            end
        end
    end
    save([destination '\V1ORIblockunitLFP.mat'],'V1ORIblockunitLFP','-v7.3');
end
sprintf('ORI done')