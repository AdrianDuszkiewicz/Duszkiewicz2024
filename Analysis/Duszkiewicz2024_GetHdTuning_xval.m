


%% Load data

    load(fullfile('Data','BehavEpochs'));
    load(fullfile('Data','SpikeData'),'S');
    load(fullfile('Data','Angle'));
    load(fullfile('Data','Velocity'));    
    
        
%     epTemp = csvread(fullfile('Analysis','LED1_TS.csv'));
%     wakeEp = intervalSet(epTemp(2,1),epTemp(2,2));

% parameters    
nbins = 360;
smoothTh = 0:12;
velTh = 2; % velocity threshold
mergeTh = 0; % threshold for merging run epochs (0 = no merging)

totC = length(S);  
totTh = length(smoothTh);    

%restrict to epoch (PoSub and ADN)

    if ~exist('wake1Ep','var')   
        wake1Ep = wakeEp;
    end
    
    ep = wake1Ep;   
    angR = Restrict(ang,ep); % make sure only tracked times are used
    angrange = Range(angR);
    ep = intervalSet(angrange(1),angrange(end)); 
    


%restrict to epoch (only A1106-180619 and A1106-180621)

%     ep = wakeOptoEp;
%     totC = length(S);   
%     angR = Restrict(ang,ep); % make sure only tracked times are used
%     angrange = Range(angR);
%     ep = intervalSet(angrange(1),angrange(end));    
%     optoEp = csvread(fullfile('Analysis','Opto_TS.csv'));
%     ep = [[Start(ep) optoEp(1,1)]; [optoEp(end,2) End(ep)]]; %  this this the epoch in between light stim bursts
%     ep = intervalSet(ep(:,1),ep(:,2));
%     fprintf(['Duration of the whole 1st epoch: ', num2str(tot_length(ep)), ' seconds\n'])

% Restrict ep to times with movement 

    epVel = thresholdIntervals(vel,velTh,'Direction','Above');   
    if mergeTh > 0
        epVel = mergeCloseIntervals(epVel, mergeTh); 
    end  
    ep = intersect(ep,epVel);
    fprintf(['Duration of 1st movement epoch: ', num2str(tot_length(ep)), ' seconds\n'])
    eps = regIntervals(ep,2);  

%% Calculate tuning curves of all cells 

    totTh = length(smoothTh);
    hAll1 = nan(nbins,totC,totTh);
    hAll2 = nan(nbins,totC,totTh);
    hdInfo1 = nan(totC,totTh);
    hdInfo2 = nan(totC,totTh);
    
    B = 2*pi*(0:1:nbins-1)'/nbins;
    ang1 = Restrict(angR,eps{1}); % Restrict angle to first half
    ang2 = Restrict(angR,eps{2}); % Restrict angle to 2nd half
    occ1 = hist(mod(Data(ang1),2*pi),B); % occupancy
    occ2 = hist(mod(Data(ang2),2*pi),B); % occupancy    

    for nT = 1:totTh
        for nC = 1:totC  
            [h,b] = HeadDirectionField(S{nC},ang1,eps{1},nbins,smoothTh(nT));
            h = h(1:end-1); 
            hAll1(:,nC,nT) = h;
            hdInfo1(nC,nT) = SpatialInfo(h,occ1);
        end
    end
    
    for nT = 1:totTh
        for nC = 1:totC  
            [h,b] = HeadDirectionField(S{nC},ang2,eps{2},nbins,smoothTh(nT));
            h = h(1:end-1); 
            hAll2(:,nC,nT) = h;
            hdInfo2(nC,nT) = SpatialInfo(h,occ2);
        end
    end
    
%% Calculate shuffled tuning curves by reversing the angle
   
    hAll1_sh = nan(nbins,totC,totTh);
    hAll2_sh = nan(nbins,totC,totTh);
    hdInfo1_sh = nan(totC,totTh);
    hdInfo2_sh = nan(totC,totTh);
        
    ang1_rev = tsd(Range(ang1),flipud(Data(ang1)));
    ang2_rev = tsd(Range(ang2),flipud(Data(ang2)));
    B = 2*pi*(0:1:nbins-1)'/nbins;
    occ1_sh = hist(mod(Data(ang1_rev),2*pi),B); % occupancy
    occ2_sh = hist(mod(Data(ang2_rev),2*pi),B); % occupancy 

    for nT = 1:totTh
        for nC = 1:totC     
            [h,b] = HeadDirectionField(S{nC},ang1_rev,eps{1},nbins,smoothTh(nT));
            h = h(1:end-1); 
            hAll1_sh(:,nC,nT) = h;
            hdInfo1_sh(nC,nT) = SpatialInfo(h,occ1_sh);
        end     
    end
    
    for nT = 1:totTh
        for nC = 1:totC     
            [h,b] = HeadDirectionField(S{nC},ang2_rev,eps{2},nbins,smoothTh(nT));
            h = h(1:end-1); 
            hAll2_sh(:,nC,nT) = h;
            hdInfo2_sh(nC,nT) = SpatialInfo(h,occ2_sh);
        end     
    end
        
%% Save analysis       
        
SaveAnalysis(pwd,'HdTuning_xval_moveEp',{hAll1; hAll2; b; occ1; occ2; occ1_sh; occ2_sh; hAll1_sh; hAll2_sh; hdInfo1; hdInfo2; hdInfo1_sh; hdInfo2_sh; smoothTh},{'hAll1'; 'hAll2'; 'b'; 'occ1'; 'occ2'; 'occ1_sh'; 'occ2_sh'; 'hAll1_sh'; 'hAll2_sh'; 'hdInfo1'; 'hdInfo2'; 'hdInfo1_sh';  'hdInfo2_sh';  'smoothTh'});
          