 function CueRotPaper_GetCueEpochs

clear all

%% Load data

    led1 = load(fullfile('Data','LED1_TS.csv'));
    led2 = load(fullfile('Data','LED2_TS.csv'));


%% Get cue change times

    % get timestamps of cue change
    led1 = led1(:,1);
    led2 = led2(:,1);
    ix = find(diff(led1) > 150 & diff(led1) < 250);
    cue1_ts = led1(ix+1);
    ix = find(diff(led2) > 150 & diff(led2) < 250);
    cue2_ts = led2([1; ix+1]);
    [tsCue,sIx] = sort([cue1_ts; cue2_ts]);
    cueID = [ones(length(cue1_ts),1); ones(length(cue2_ts),1)*2]; % 1's and 2's for cue ID
    cueID = cueID(sIx); % sort exactly like cue ts
    cueEp = [tsCue [tsCue(2:end); tsCue(end)+200]]; % vector with starts and ends of each cue epoch

    
%% Save data

SaveAnalysis(pwd,'CueEpochs',{cueEp; cueID},{'cueEp'; 'cueID'});
