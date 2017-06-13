function [V1diff,amyspecdiff]=diff_amyspec()
%%
% run after execute amyspec
%%
folderlist={
    'I:\160712';
    'I:\160713';
    'I:\160714';
    'I:\160715';
    'I:\160716';
    'I:\160719';
    'I:\160720';
    'I:\160721';
    'I:\160722';
    'G:\MS\160624';
    'G:\MS\160625';
    'G:\MS\160627';
    'G:\MS\160628';
    };
LRlist='rrrrrrrrrllll';
% LRlist='llll';
fqpass=[60 90];
win=[430 470];
V1diff=[];
amyspecdiff=[];
figure
for idays=1:numel(folderlist)
    lr=strcmp(LRlist(idays),'r');  % L==0 R==1;
    load([folderlist{idays} '\amy_spec_601_1000ms.mat']);
    load([folderlist{idays} '\Blksortedcell.mat']);
    load([folderlist{idays} '\V1_unit_block_diffpsth.mat']);
    fidx=amy_spec.f>fqpass(1) & amy_spec.f<fqpass(2);
    l=cellfun(@length,v1unitblockdiffpsth);
    v1unitblockdiffpsth(l==0)=[];
    celldiff=cellfun(@(x) merge_cell(x,win,lr),v1unitblockdiffpsth,'UniformOutput',false);
    alldiff=cell2mat(celldiff);
    mergediff=mean(alldiff,2);  %[CSdiff_block1; NEUdiff_block1; CSdiff_block2; NEUdiff_block2; CSdiff_block3; NEUdiff_block3....]
    V1diff=[V1diff; mergediff];
    
    clockcellnum=cellfun(@(x) sum(sum(x)),Blksortedcell');
    edidx=cumsum(clockcellnum);
    stidx=[0; edidx(1:end-1)]+1;
    blockdiff=arrayfun(@(x,y) genspecdiffvector(x,y,amy_spec,fidx),stidx,edidx,'UniformOutput',false);
    specdiff=cell2mat(blockdiff);
    amyspecdiff=[amyspecdiff; specdiff];
    subplot(3,5,idays)
    plot(mergediff.*1000,specdiff,'or')% ,'color',[1-double(idays)/numel(folderlist) double(idays)/numel(folderlist) double(idays)/numel(folderlist)]
    [R,P]=corrcoef(mergediff(~isnan(specdiff)),specdiff(~isnan(specdiff)));
    title([ 'R=' num2str(R(1,2)) '  P=' num2str(P(1,2))]);
    hold on
end
figure
plot(V1diff.*1000,amyspecdiff,'or')
sprintf('done')
function celldiff=merge_cell(unit,win,lr)
blockdiff=cellfun(@(x) merge_block(x,win,lr),unit,'UniformOutput',false);
celldiff=cell2mat(blockdiff');

function blockdiff=merge_block(block,win,lr)
cnddiff=cellfun(@(x) mean(mean(x(:,win(1):win(2)))),block);
blockdiff=[mean(cnddiff((1:3)+~lr*3)); mean(cnddiff((4:6)-~lr*3))];

function specdiff=genspecdiffvector(a,b,amy_spec,fidx)
specdiff=[mean(sum(amy_spec.S1(fidx,a:b)),2); mean(sum(amy_spec.S2(fidx,a:b)),2)];







