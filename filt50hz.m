function Y = filt50hz(X,sf)  % X's format: trials x bins
    nfilt=3000;
    Wo = 50/(sf/2);  BW = Wo/10;
    [b,a] = iirnotch(Wo,BW);  
    Wo2 = 100/(sf/2);  BW2 = Wo2/35;
    [b2,a2] = iirnotch(Wo2,BW2);  
    siz = size(X);
    Y = zeros(siz);
    for i=1:ceil(siz(2)/nfilt)
        if i*nfilt<siz(2)
            loopidx=(i-1)*nfilt+1:i*nfilt;
        else
            loopidx=(i-1)*nfilt+1:siz(2);
        end
        temp=filter(b,a,X(:,loopidx));
        Y(:,loopidx)=filter(b2,a2,temp);
    end
end