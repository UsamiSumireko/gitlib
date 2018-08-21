clear
str='h:\180413';
puffcnd=0:2;
neucnd=3:11;
slashidx=strfind(str,'\');
slashidx=slashidx(end);

namelist=dir(str);
ORICONnum=[];
ORICONlist=[];
ORICONlr=[];
DftORIlist=[];
concount=0;
oricount=0;
for i=3:numel(namelist)
    if strcmp(namelist(i).name(1:6),'ORICON') && namelist(i).isdir
        concount=concount+1;
        ORICONnum(concount)=str2num(namelist(i).name(end-2:end));
        ORICONlist=[ORICONlist; {namelist(i).name}];
        ORICONlr(concount)=strcmpi(namelist(i).name(7),'l');% R=0 L=1
    elseif strcmp(namelist(i).name(1:6),'DftORI') && namelist(i).isdir
        oricount=oricount+1;
        DftORIlist=[DftORIlist; {namelist(i).name}];
    end
end
[~,Bindex]=sort(ORICONnum);
temp=ORICONlr(Bindex);
reverselabel=find((temp(1:end-1)-temp(2:end))~=0);
reverselabel=reverselabel+1;


for nrev=0:numel(reverselabel)
    if 0==nrev
         blockidx=1:reverselabel(1)-1;
         prerevLR=int8(ORICONlr(blockidx(1)));
         pretraceCS=[];
         pretraceNEU=[];
         precloseCS=[];
         precloseNEU=[];         
         for iblock=blockidx
            path=[str '\' ORICONlist{Bindex(iblock)}];
            [TraceTimeCS,TraceTimeNEU,PCloseTime,NCloseTime,PNBlink,NNBlink]=BlinkImageMKII(path,puffcnd,neucnd,0);
            pretraceCS=[pretraceCS TraceTimeCS'];
            pretraceNEU=[pretraceNEU TraceTimeNEU'];
            precloseCS=[precloseCS PCloseTime'];
            precloseNEU=[precloseNEU NCloseTime'];
         end
         figure
         bar([2 4],[mean(mean(precloseCS)) mean(mean(precloseNEU))]);
         xlabel('puff/neu trials')
         ylabel('mean eye close time')
         figure
         xv=1:size(precloseCS,2);
         yv=[mean(precloseCS,1)' mean(precloseNEU,1)'];
         bar(xv,yv)
         xlabel('block num')
         ylabel('mean eye close time')
    else 
        if numel(reverselabel)==nrev
            blockidx=reverselabel(nrev):numel(ORICONlist);
        else
            blockidx=reverselabel(nrev):reverselabel(nrev+1)-1;
        end
        traceCS=[];
        traceNEU=[];
        closeCS=[];
        closeNEU=[];
        count=0;
        for iblock=blockidx
            count=count+1;
            path=[str '\' ORICONlist{Bindex(iblock)}];
            [TraceTimeCS,TraceTimeNEU,PCloseTime,NCloseTime,PNBlink,NNBlink]=BlinkImageMKII(path,puffcnd,neucnd,0);
            traceCS{nrev}(:,count)=TraceTimeCS';
            traceNEU{nrev}(:,count)=TraceTimeNEU';
            closeCS{nrev}(:,count)=PCloseTime';
            closeNEU{nrev}(:,count)=NCloseTime';
        end
        xv=2:2:2*count;
        yv1=mean(closeCS{nrev},1);
        yv2=mean(closeNEU{nrev},1);
        yv=[yv1' yv2'];
        figure
        bar(xv,yv);
        xlabel('block number')
        ylabel('mean eye close time')
    end
end

