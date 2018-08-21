function p=bootstraptest(sample1,sample2,nbootstrap)
%%
% p=bootstraptest(sample1,sample2,nbootstrap)
% bootstrap test with put back
% resample from sample 1 and sample 2, shuold be column vector
%%
if nargin<3
    nbootstrap=1000;
end
if 1==size(sample1,1)
    sample1=transpose(sample1);
end
if 1==size(sample2,1)
    sample2=transpose(sample2);
end

samplecount1=size(sample1,1);
samplepool=cat(1,sample1,sample2);
allsamplecount=size(samplepool,1);
randidx=randi(allsamplecount,allsamplecount,nbootstrap);
randsample=samplepool(randidx);
samplediff=mean(randsample(1:samplecount1,:),1)-mean(randsample(1+samplecount1:end,:),1);
samplediffstd=std(samplediff);
realdiff=mean(sample1,1)-mean(sample2,1);
zs=abs(realdiff-mean(samplediff))/samplediffstd;
p=1-cdf('Normal',zs,0,1)+cdf('Normal',-zs,0,1);    
end