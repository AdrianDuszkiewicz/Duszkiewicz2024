

%% Parameters

binsize = 1;
nbins = 100;


%% Load data

load(fullfile('Data','SpikeData'));


%% Compute AutoCorrs

totC = length(S);
allAC = nan(nbins+1,totC);


for n = 1:totC                    
    [ac, bins] = CrossCorr(Range(S{n},'ms'),Range(S{n},'ms'),binsize,nbins); 
    ac(bins == 0) = 0;
     rx = length(Range(S{n}));
     ac = ac*rx*binsize/1000;
    allAC(:,n) = ac;   
end
              
SaveAnalysis(pwd,'AutoCorrs',{allAC; bins},{'acorrs'; 'bins'});
    