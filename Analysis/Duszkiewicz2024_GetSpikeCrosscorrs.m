
%% parameters and data
   
binsize = 0.2;
nbins = 250;

load(fullfile('Data','SpikeData'),'S');

%% Compute CrossCorrs

fprintf('Computing spike train CCGs \n');

totC = length(S);
allCCG = nan(nbins+1,nchoosek(totC,2)*2);
allPairs = [];

noCCG = 1;

for nC1 = 1:totC   
    tempPairs = [];   
    for nC2 = 1:totC      
        if nC1 ~= nC2            
            [tempCCG, bins] = CrossCorr(Range(S{nC1},'ms'),Range(S{nC2},'ms'),binsize,nbins);           
            rx = length(Range(S{nC1}));           
            tempCCG = tempCCG*rx*binsize/1000;         
            allCCG(:,noCCG) = tempCCG;           
            tempPairs(end+1,:) = [nC1 nC2];          
            noCCG = noCCG+1;                    
        end
    end        
   allPairs = [allPairs; tempPairs];   
end
  
SaveAnalysis(pwd,'SynCrossCorrs',{allPairs, allCCG, bins},{'pairs', 'ccg', 'bins'});



