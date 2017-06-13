function diff_svm(mthd)
%%
% diff_svm(mthd)
% mthd= 'spk' or 'lfp'
% run after execute amysvm_lfp or amysvm_spk
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
    };
win=[430 470];
score1all=[];
score2all=[];
V1diff1=[];
V1diff2=[];
crct1=[];
crct2=[];
right=0;
allc=0;
%%
for n=1:numel(folderlist)
    str=folderlist{n};
    load([str '\V1_unit_block_diffpsth.mat'])
    l=cellfun(@length,v1unitblockdiffpsth);
    v1unitblockdiffpsth(l==0)=[];
    V1diff=cellfun(@(x) merge_cell(x,win),v1unitblockdiffpsth,'UniformOutput',false);
    v1spkdiff=cell2mat(V1diff);
    Z=mean(v1spkdiff,2);
    z=[];
%%
    if strcmpi(mthd,'spk')
        [X,Y,ranidx,SVMModel,score1,score2,pfm]=amysvm_spk(str);
        crct1=[crct1; Y(ranidx)'==(score1(:,2)>0)];
        crct2=[crct2; Y(~ranidx)'==(score2(:,2)>0)];
        score1all=[score1all; score1(:,2)];
        score2all=[score2all; score2(:,2)];
        V1diff1=[V1diff1; Z(ranidx)];
        V1diff2=[V1diff2; Z(~ranidx)];
        right=right+pfm;
        allc=allc+sum(~ranidx);
        sprintf([num2str(right) '/' num2str(allc)])
%%
    elseif strcmpi(mthd,'lfp')
        load([str '\Blksortedcell']);
        blkcell=cellfun(@(x) sum(sum(x)),Blksortedcell);
        [X,Y,ranidx,SVMModel,score1,score2,pfm]=amysvm_lfp(str);
        for i=1:6:numel(Z)
            z1=mean(Z(i:i+2));
            z2=mean(Z(i+3:i+5));
            z=[z; repmat([z1;z2],blkcell(floor(i/6)+1),1)];
        end
        crct1=[crct1; Y(ranidx)'==(score1(:,2)>0)];
        crct2=[crct2; Y(~ranidx)'==(score2(:,2)>0)];
        score1all=[score1all; score1(:,2)];
        score2all=[score2all; score2(:,2)];
        V1diff1=[V1diff1; z(ranidx)];
        V1diff2=[V1diff2; z(~ranidx)];
        right=right+pfm;
        allc=allc+sum(~ranidx);
        sprintf([num2str(right) '/' num2str(allc)])
    end
end
figure
subplot(1,2,1)
plot(V1diff1(logical(crct1)),score1all(logical(crct1)),'b.')
hold on
plot(V1diff1(~logical(crct1)),score1all(~logical(crct1)),'r.')
subplot(1,2,2)
plot(V1diff2(logical(crct2)),score2all(logical(crct2)),'b.')
hold on
plot(V1diff2(~logical(crct2)),score2all(~logical(crct2)),'r.')
sprintf('done')
% figure
% subplot(1,2,1)
% plot(Z(ranidx),score1(:,2),'.')
% subplot(1,2,2)
% plot(Z(~ranidx),score2(:,2),'.')
%%
function a=merge_cell(unit,win)
aa=cellfun(@(x) merge_block(x,win),unit,'UniformOutput',false);
a=cell2mat(aa);
a=a';
function a=merge_block(block,win)
a=cellfun(@(x) mean(mean(x(:,win(1):win(2)))),block);






