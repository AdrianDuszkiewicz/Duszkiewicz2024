function PoSubPaper_GetTCxcorrs

load(fullfile('Analysis','HdTuning_moveEp'));

smoothTC = 3;  % smoothing of tuning curves (3 used in Duszkiewicz et al., 2024)

hAll = hAll(:,:,smoothTC+1);

if size(hAll,1) > 360    
    hAll = hAll(1:end-1,:);
end

[bins,totC] = size(hAll); 
allPairs = nchoosek(1:totC,2); 
   
% calculate TC xcorrs and get pairwise angular differences

totPairs = size(allPairs,1);
ccgTC = nan(bins,totPairs);

for nP = 1:totPairs    
    ccg = TCcrosscorr(hAll(:,allPairs(nP,1)),hAll(:,allPairs(nP,2)));
    ccgTC(:,nP) = ccg; 
end

SaveAnalysis(pwd,'TCxcorrs',{allPairs, ccgTC},{'pairsTC', 'ccgTC'});


