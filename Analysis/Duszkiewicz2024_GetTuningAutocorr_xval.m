
   
load(fullfile('Analysis','HDtuning_xval_moveEp'));

hAllA = hAll1;
hAllB = hAll2;
hAllA_sh = hAll1_sh;
hAllB_sh = hAll2_sh;

[totBins,totC,totTh] = size(hAllA);


 %% calculate TC crosscorrs

ac = nan(size(hAllA));
for nTh = 1:totTh 
    h1 = zscore(hAllA(:,:,nTh));
    h2 = zscore(hAllB(:,:,nTh));
    for n = 0:size(ac,1)-1
        hs = circshift(h2,n);
        ac(n+1,:,nTh) = diag(h1' * hs);
    end
end

acAll = ac / (size(ac,1)-1);

ac_sh = nan(size(hAllA));
for nTh = 1:totTh 
    h1 = zscore(hAllA_sh(:,:,nTh));
    h2 = zscore(hAllB_sh(:,:,nTh));
    for n = 0:size(ac_sh,1)-1
        hs = circshift(h2,n);
        ac_sh(n+1,:,nTh) = diag(h1' * hs);
    end
end

acAll_sh = ac_sh / (size(ac_sh,1)-1);

step = 360./totBins;
lags = 0:step:360;
lags = lags(1:end-1);


 %% Save crosscorrs and FFT
 
 SaveAnalysis(pwd,'TCautocorrs_xval_moveEp',{acAll; acAll_sh; lags},{'ac'; 'ac_sh'; 'lags'});