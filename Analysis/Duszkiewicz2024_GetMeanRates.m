

getREM = 0;

load(fullfile('Data','BehavEpochs'));
load(fullfile('Data','SpikeData'),'S');
[~, foldername, ~] = fileparts(pwd);

if getREM == 1    
    load(fullfile('Sleep',[foldername '.SleepState.states.mat'] ));
    % get REM episodes
    rem = SleepState.ints.REMstate;
    remDur = sum(rem(:,2) - rem(:,1));
    epRem = intervalSet(rem(:,1),rem(:,2));
    rateREM = Rate(S,epRem); 
else
    rateREM = nan(length(S),1);
end
% calculate mean wake FR

% rate in square 
epS = wake1Ep;
rateS = Rate(S,epS);

% rate in triangle
epT = wake2Ep;
rateT = Rate(S,epT);

SaveAnalysis(pwd,'MeanFR',{rateS, rateT, rateREM},{'rateS', 'rateT', 'rateREM'});
%SaveAnalysis(pwd,'MeanFR',{rateS},{'rateS'});    