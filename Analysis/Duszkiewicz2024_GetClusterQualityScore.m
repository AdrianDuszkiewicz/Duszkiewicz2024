

load(fullfile('Analysis','AutoCorrs'));

midBin = floor(length(bins)./2);
mxAc = max(acorrs,[],1);
midAc = acorrs(midBin,:);
cqs = [midAc./mxAc]';

SaveAnalysis(pwd,'ClusterQualityScore',{cqs},{'cqs'});


