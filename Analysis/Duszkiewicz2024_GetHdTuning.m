


%% Parameters and data

% parameters    
nbins = 360;
smoothTh = 0:12; 
velTh = 2; % velocity threshold (2 = default)
mergeTh = 0; % threshold for merging run epochs (0 = no merging)

%load data
    load(fullfile('Data','BehavEpochs'));
    load(fullfile('Data','SpikeData'),'S');
    load(fullfile('Data','Angle'));
    load(fullfile('Data','Velocity'));
        
totC = length(S);  
totTh = length(smoothTh);    
%% Calculate tuning curves for all cells in first environment

hAll = nan(nbins,totC,totTh);
hdInfo = nan(totC,totTh);
hAll_sh = nan(nbins,totC,totTh);
hdInfo_sh = nan(totC,totTh);
    
%restrict to epoch (only A1106-180619 and A1106-180621)

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
    
% restrict to epoch (Cue Rotation dataset)
%     epTemp = csvread(fullfile('Analysis','LED1_TS.csv'));
%     wakeEp = intervalSet(epTemp(2,1),epTemp(2,2));

% restrict to epoch (All the rest)

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

% Real tuning curves
    angR = Restrict(angR,ep); % restrict angle to movement epoch
    B = 2*pi*(0:1:nbins-1)'/nbins;
    occ = hist(mod(Data(angR),2*pi),B); % occupancy in movement epoch

    for nT = 1:totTh
        for nC = 1:totC  
            [h,b] = HeadDirectionField(S{nC},angR,ep,nbins,smoothTh(nT)); % using restricted angle on purpose
            h = h(1:end-1); 
            hAll(:,nC,nT) = h;
            hdInfo(nC,nT) = SpatialInfo(h,occ);
        end
    end
    
% Shuffled tuning curves (reversed angle)
    angR_rev = tsd(Range(angR),flipud(Data(angR))); % flip the times of restricted angle (so occupancy is the same as in real)   
    B = 2*pi*(0:1:nbins-1)'/nbins;
    occ_sh = hist(mod(Data(angR_rev),2*pi),B); % occupancy in movement epoch, should be the same as real occ

    for nT = 1:totTh
        for nC = 1:totC     
            [h,b] = HeadDirectionField(S{nC},angR_rev,ep,nbins,smoothTh(nT)); % epoch doesn't change so ok to use the same ang restriction
            h = h(1:end-1); 
            hAll_sh(:,nC,nT) = h;
            hdInfo_sh(nC,nT) = SpatialInfo(h,occ_sh);
        end     
    end
    
%% Calculate tuning curves of all cells in the second environment (if present)

hAllT = nan(nbins,totC,totTh);
hAllT_sh = nan(nbins,totC,totTh);
hdInfoT = nan(totC,totTh);
hdInfoT_sh = nan(totC,totTh);
occT = nan(size(occ));
occT_sh = nan(size(occ_sh));


if ~exist('wake2Ep','var') 
    wake2Ep = intervalSet([],[]);
end

if isempty(Start(wake2Ep)) % for processing of sessions without the triangle
    tri = zeros(totC,1);
else
    tri = ones(totC,1);
   
% restrict to epoch (PoSub)
    ep = wake2Ep; 
    angR = Restrict(ang,ep); 
    angrange = Range(angR); % make sure only tracked times are used
    ep = intervalSet(angrange(1),angrange(end));
    
% Restrict ep to times with movement 
    velR = Restrict(vel, ep);
    ep = thresholdIntervals(velR,velTh,'Direction','Above');
    if mergeTh > 0
        ep = mergeCloseIntervals(ep, mergeTh); 
    end
    

   fprintf(['Duration of 2nd movement epoch: ', num2str(tot_length(ep)), ' seconds\n'])

% Real tuning curves    
    angR = Restrict(angR,ep); % restrict angle to movement epoch
    B = 2*pi*(0:1:nbins-1)'/nbins;
    occT = hist(mod(Data(angR),2*pi),B); % occupancy in movement epoch

    for nT = 1:totTh
        for nC = 1:totC  
            [h,b] = HeadDirectionField(S{nC},angR,ep,nbins,smoothTh(nT));
            h = h(1:end-1); 
            hAllT(:,nC,nT) = h;
            hdInfoT(nC,nT) = SpatialInfo(h,occT);
        end
    end
    
% Shuffled tuning curves (reversed angle)       
    angR_rev = tsd(Range(angR),flipud(Data(angR)));
    B = 2*pi*(0:1:nbins-1)'/nbins;
    occT_sh = hist(mod(Data(angR_rev),2*pi),B); % occupancy (should be the same as real occ)

    for nT = 1:totTh
        for nC = 1:totC     
            [h,b] = HeadDirectionField(S{nC},angR_rev,ep,nbins,smoothTh(nT));
            h = h(1:end-1); 
            hAllT_sh(:,nC,nT) = h;
            hdInfoT_sh(nC,nT) = SpatialInfo(h,occT_sh);
        end     
    end
end    

b = b(1:end-1);
%% Save analysis       
        
SaveAnalysis(pwd,'HdTuning_moveEp',{hAll; hAll_sh; occ; occ_sh; hdInfo; hdInfo_sh; hAllT; hAllT_sh; occT; occT_sh; hdInfoT; hdInfoT_sh; smoothTh; tri; b},{'hAll'; 'hAll_sh'; 'occ'; 'occ_sh'; 'hdInfo'; 'hdInfo_sh'; 'hAllT'; 'hAllT_sh'; 'occT'; 'occT_sh'; 'hdInfoT'; 'hdInfoT_sh'; 'smoothTh'; 'tri'; 'b'});
          