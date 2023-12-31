


smoothTC = 3; % smoothing of tuning curves (3 used in Duszkiewicz et al., 2024)


%% load data
load(fullfile('Analysis','MeanFR'));
load(fullfile('Data','WaveformFeatures'));
load(fullfile('Analysis','HdTuning_moveEp'));
%load(fullfile('Analysis','CellDepth'));

goodSpk = ~isnan(tr2pk);

%frate = rateO; % for Opto dataset

frate = rateS;
gd = frate > 0.5 & goodSpk == 1;
hdInfo = hdInfo(:,smoothTC+1);

%gd(cellDep > 0) = 0; % For the neuropixel rec (A8303)
%% define groups

% FS cells
fs = frate >= 10 & gd == 1 & tr2pk <= 0.35;

% Excitatory cells
ex = gd == 1 & tr2pk > 0.35 & frate <= 10;

% HD cells
hd = ex == 1 & hdInfo >= 0.2;

% Non-HD cells
nhd = ex == 1 & hdInfo < 0.2;

%% 
SaveAnalysis(pwd,'CellTypes',{gd; hd; fs; ex; nhd},{'gd'; 'hd'; 'fs'; 'ex'; 'nhd'});
    
