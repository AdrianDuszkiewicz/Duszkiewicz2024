

%% Parameters and data

% parameters    
nbins = 360;
smoothTh = 3; 
velTh = 2; % velocity threshold
mergeTh = 0; % threshold for merging run epochs (0 = no merging)
lagVec = [-0.5:0.01:0.5]; 

%load data
    load(fullfile('Data','BehavEpochs'));
    load(fullfile('Data','SpikeData'),'S');
    load(fullfile('Data','Angle'));
    load(fullfile('Data','Velocity'));

 
%% Calculate tuning curves for all cells in first environment

totLags = length(lagVec); 
totC = length(S);    
hAll = nan(nbins,totC,totLags);
hdInfo = nan(totLags,totC);
    
%restrict to epoch (ADN opto)

%     ep = wakeOptoEp;
%     wake2Ep = intervalSet([],[]);
%     totC = length(S);   
%     angR = Restrict(ang,ep); % make sure only tracked times are used
%     angrange = Range(angR);
%     ep = intervalSet(angrange(1),angrange(end));    
%     optoEp = csvread(fullfile('Analysis','Opto_TS.csv'));
%     epStart = optoEp(find(diff(optoEp(:,2)) > 10),2)+4; 
%     epEnd = optoEp(find(diff(optoEp(:,1)) > 10) + 1,1);
%     ep = [[Start(ep) optoEp(1,1)]; [epStart epEnd]; [optoEp(end,2) End(ep)]]; %  this this the epoch in between light stim bursts
%     ep = intervalSet(ep(:,1),ep(:,2));
%     wake1Ep = ep;
    
%restrict to epoch (ADN opto - only A1106-180619)

%     ep = wakeOptoEp;
%     wake2Ep = intervalSet([],[]);
%     totC = length(S);   
%     angR = Restrict(ang,ep); % make sure only tracked times are used
%     angrange = Range(angR);
%     ep = intervalSet(angrange(1),angrange(end));    
%     optoEp = csvread(fullfile('Analysis','Opto_TS.csv'));
%     ep = [[Start(ep) optoEp(1,1)]; [optoEp(end,2) End(ep)]]; %  this this the epoch in between light stim bursts
%     ep = intervalSet(ep(:,1),ep(:,2));
%     wake1Ep = ep;

% restrict to epoch (Not opto)

    if ~exist('wake1Ep','var')   
        wake1Ep = wakeEp;
    end
    
    ep = wake1Ep;  
    angR = Restrict(ang,ep);
    angrange = Range(angR);  % make sure only tracked times are used
    ep = intervalSet(angrange(1),angrange(end));
    
% Restrict ep to times with movement 
    epVel = thresholdIntervals(vel,velTh,'Direction','Above');
    if mergeTh > 0
        epVel = mergeCloseIntervals(epVel, mergeTh); 
    end
    fprintf(['Duration of the whole 1st epoch: ', num2str(tot_length(ep)), ' seconds\n'])
    ep = intersect(ep,epVel);
    fprintf(['Duration of 1st movement epoch: ', num2str(tot_length(ep)), ' seconds\n'])
    
% calculate tuning curves and HD info for each cell   
    rng = Range(ang);    
    for nL = 1:totLags
        t = rng - lagVec(nL); %Negative time slide means past
        angS = tsd(t,Data(ang));       
        for nC = 1:totC
            [h b hdi] = HeadDirectionField(S{nC},angS,ep,nbins,smoothTh);
            hAll(:,nC,nL) = h(1:end-1);
            hdInfo(nL,nC) = hdi;
        end
    end
                          
b = b(1:end-1);

%% Save analysis       
        
SaveAnalysis(pwd,'HdTuning_AnticipCorr',{hAll; lagVec; hdInfo; b},{'hAll'; 'lagVec'; 'hdInfo_lag';'b'});
          