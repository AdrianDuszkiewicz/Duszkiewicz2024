
%function CueRotPaper_GetSigmoid


% Copyright Adrian Duszkiewicz and Eliott Owczarek 2020-2021
% Sigmoid fitting code incorporates elements of:
% R P (2021). sigm_fit (https://www.mathworks.com/matlabcentral/fileexchange/42641-sigm_fit), MATLAB Central File Exchange. Retrieved June 9, 2021.

clear all

%% Load data

    load(fullfile('Analysis','BayesianDecoding_coarse')); 
    load(fullfile('Analysis','BayesianDecoding_fine')); 
    load(fullfile('Data','CueEpochs'));
    load(fullfile('Data','Ahv')); 

    tsCue = cueEp(:,1);
    totEps = length(tsCue);

%% Fit sigmoids into transitions (coarse)

   % sigmoid function parameters for nlinfit
   sigmoidFun = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4))); % create function handle for a sigmoid function

    cueWin = [50 200]; % before and after cue rotation

    binBayes = round(median(diff(Range(errBayes))),6); % get bin size of Bayesian decoding
    sigVal = []; % sigmoid values
    sigTS = []; % sigmoid timestamps
    sigErr = []; % decoder error in sigmoid TS

    sigDiff = [];
    sigDir = zeros(totEps,1); 

    for nEp = 1:totEps
        int = intervalSet(cueEp(nEp,1)-cueWin(1), cueEp(nEp,1)+cueWin(2)); % make intervalSet to fit sigmoid into
        err = Restrict(errBayes,int); % restrict error to the above  
        x = Range(err); % x values for estimation
        y = round(Data(err),6);
        
        % Fit sigmoid       
        beta0 = [quantile(y,0.05) quantile(y,0.95) NaN 1]; % initial parameters for possible sigmoid (nlinfit needs some starter values)                
        % third beta0 parameter is the value of x at 0.5 quantile of y
         if sum(y == quantile(y,0.5)) == 0
            beta0(3) = x(y==quantile(y(2:end),0.5)); % sometimes it fails to find the y value equal to this quantile, no idea why, but (2:end) fixes it
         else
            beta0(3) = x(y==quantile(y,0.5));
         end
        
        beta = nlinfit(x,y,sigmoidFun,beta0); % fit sigmoid using function 'modelfit' and initial paramaters 'beta0'
        yval = sigmoidFun(beta,x); % evaluate function with new parameters 'beta' over the timestamps 
        
        % TO DO: catch warnings and pad with NANs
        sigVal = [sigVal yval]; % values of sigmoid
        sigTS = [sigTS x]; % timestamps of sigmoid
        sigErr = [sigErr y];

        %evaluate the sigmoid 
        ev = abs(angdiff([y' ; yval']))';
        sigDiff = [sigDiff ev]; 

        %find direction of sigmoid   
        if yval(1) < yval(end)
            sigDir(nEp) = 1;
        elseif yval(1) > yval(end)
            sigDir(nEp) = -1;
        end    

    end

    meanSigDiff = rad2deg(nanmean(sigDiff,1))'; % mean difference between sigmoid fit and actual data

%% Determine valid transitions 
 
    valEp = 50; % how many seconds at the end of epoch to use for validation 
    lateEps = [tsCue(1) - 200; tsCue; tsCue(end)+200]; % add first and last timepoint to have something to compare to
    meanErr = nan(totEps+1,1); 
    varErr = nan(totEps+1,1);
    epAhv = nan(totEps+1,1);

    for nEp = 1:totEps+1    
        ep = intervalSet(lateEps(nEp)+200-valEp, lateEps(nEp+1)); 
        err = Restrict(errBayes,ep);
        err = Data(err); 
        err(isnan(err)) = []; % remove NANs, circstat functions don't take them
        meanErr(nEp) = circ_mean(err); % mean error to determine if shift happened   
        varErr(nEp) = circ_var(err); % variance of error to determine stability
        a = Restrict(ahv,ep);
        epAhv(nEp) = mean(abs(Data(a)));   
    end

    errChange = rad2deg(abs(angdiff(meanErr))); % how much did error change from previous epoch 
    varErr = varErr(1:end-1); % we don't need the last value, it's after the last cue change
    epAhv = epAhv(1:end-1); % we don't need the last value, it's after the last cue change
    meanErr = rad2deg(meanErr); 

%% Exclusion criteria 

    errTh = 45; % threshold for rotation
    varTh = 10; % threshold for error variance
    ahvTh = 10; % threshold for AHV
    sigTh = 15; % threshold for mean sigmoid difference

    % indices of included and excluded epochs
    isExcl = abs(diff(meanErr)) < errTh | rad2deg(varErr) > varTh | epAhv < ahvTh; % conditions
    ixExcl = find(isExcl); % excluded
    ixIncl = find(~isExcl); % included


%% Find beginning, mid and end of sigmoids (coarse)

    % start of sigmoid based on whether higher than initial error (in std)
    % unreliable, try calculating percentiles (say 10 and 90% amplitude)
    
    thStart = 0.01; % set the start point of sigmoid
    thEnd = 1 - thStart; 
    
    sigStart = nan(totEps,1); 
    sigMid = nan(totEps,1);
    sigEnd = nan(totEps,1); 
    
    for nEp = 1:totEps          
        sig = sigVal(:,nEp);
        sig = rescale(sig); % rescale sigmoid to 0:1
        if sigDir(nEp) == 1
            [~,ix] = min(abs(thStart - sig)); % find the point in sigmoid closest to start threshold
            sigStart(nEp) = sigTS(ix,nEp); 
            [~,ix] = min(abs(thEnd - sig)); % find the point in sigmoid closest to end threshold
            sigEnd(nEp) = sigTS(ix,nEp);  
        elseif sigDir(nEp) == -1 % reverse is true if the sigmoid is reversed
            [~,ix] = min(abs(thEnd - sig)); 
            sigStart(nEp) = sigTS(ix,nEp); 
            [~,ix] = min(abs(thStart - sig)); 
            sigEnd(nEp) = sigTS(ix,nEp);       
        end
            [~,ix] = min(abs(0.5 - sig)); % midpoint same for both directions
            sigMid(nEp) = sigTS(ix,nEp);          
    end
 
sigTimes = [sigStart sigMid sigEnd];
  
% STD-based method   
%     noStd = 2; 
%     sigStart = nan(totEps,1); 
%     for nEp = 1:totEps   
%         ix1 = find(sigTS(:,nEp) < tsCue(nEp)); % find sigmoid timestamps before cue change
%         errStd = nanstd(sigErr(ix1,nEp)); % not ideal because angular value
% 
%         if sigDir(nEp) == 1
%            th = nanmean(sigErr(ix1,nEp))+errStd*noStd;
%            ix2 = min(find(sigVal(:,nEp) > th)); 
%         elseif sigDir(nEp) == -1
%            th = nanmean(sigErr(ix1,nEp))-errStd*noStd;
%            ix2 = min(find(sigVal(:,nEp) < th));        
%         end
%         if ~isempty(ix2)
%             sigStart(nEp) = sigTS(ix2,nEp); 
%         end
%     end

%% Fit sigmoids (fine res)

   % sigmoid function parameters for nlinfit
   sigmoidFun = @(param,xval) param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4))); % create function handle for a sigmoid function

    cueWin = [5 5]; % before and after sigmoid

    sigFine = {}; % (1) timestamps (2) sigmoid values (3) error values


    meanSigDiffF = [];

    for nEp = 1:totEps
        int = intervalSet(sigStart(nEp)-cueWin(1), sigEnd(nEp)+cueWin(2)); % make intervalSet to fit sigmoid into
        err = Restrict(errBayesFine,int); % restrict error to the above       
        x = Range(err); % x values for estimation
        y = Data(err);
        
        % Fit sigmoid       
        beta0 = [quantile(y,0.05) quantile(y,0.95) NaN 1]; % initial parameters for possible sigmoid (nlinfit needs some starter values)                 
        % third beta0 parameter is the value of x at 0.5 quantile of y
         if sum(y == quantile(y,0.5)) == 0
            beta0(3) = x(y==quantile(y(2:end),0.5)); % sometimes it fails to find the y value equal to this quantile, no idea why, but (2:end) fixes it
         else
            beta0(3) = x(y==quantile(y,0.5));
         end
        
        beta = nlinfit(x,y,sigmoidFun,beta0); % fit sigmoid using function 'modelfit' and initial paramaters 'beta0'
        yval = sigmoidFun(beta,x); % evaluate function with new parameters 'beta' over the timestamps 
        
        %evaluate the sigmoid 
        ev = abs(angdiff([y' ; yval']))';
        meanSigDiffF = [meanSigDiffF; rad2deg(nanmean(ev))]; 
        
        % TO DO: catch warnings and pad with NANs
        sigFine{nEp} = [x yval y]; % values of sigmoid
  

    end
    
%% Find beginning, mid and end of sigmoids (fine)

    % start of sigmoid based on whether higher than initial error (in std)
    % unreliable, try calculating percentiles (say 10 and 90% amplitude)
    
    thStart = 0.2; % set the start point of sigmoid
    thEnd = 1 - thStart; 
    
    sigStartF = nan(totEps,1); 
    sigMidF = nan(totEps,1);
    sigEndF = nan(totEps,1); 
    
    for nEp = 1:totEps          
        d = sigFine{nEp};
        sig = d(:,2);
        sig = rescale(sig); % rescale sigmoid to 0:1
        if sigDir(nEp) == 1
            [~,ix] = min(abs(thStart - sig)); % find the point in sigmoid closest to start threshold
            sigStartF(nEp) = d(ix,1); 
            [~,ix] = min(abs(thEnd - sig)); % find the point in sigmoid closest to end threshold
            sigEndF(nEp) = d(ix,1); 
        elseif sigDir(nEp) == -1 % reverse is true if the sigmoid is reversed
            [~,ix] = min(abs(thEnd - sig)); 
            sigStartF(nEp) = d(ix,1); 
            [~,ix] = min(abs(thStart - sig)); 
            sigEndF(nEp) = d(ix,1);       
        end
            [~,ix] = min(abs(0.5 - sig)); % midpoint same for both directions
            sigMidF(nEp) = d(ix,1);          
    end
 
% calculate duration, if start < tsCue then take tsCue

 remapDur = nan(16,1); 
    for nEp = 1:totEps
        if ~isExcl(nEp)
            %if sigStartF(nEp) > tsCue(nEp)
                remapDur(nEp) = sigEndF(nEp) - sigStartF(nEp);
%             else
%                 remapDur(nEp) = sigEndF(nEp) - tsCue(nEp); 
            %end
        end
    end
    
sigTimesF = [sigStartF sigMidF sigEndF];  


%% Save analysis

%SaveAnalysis(pwd,'Sigmoids',{sigVal; sigTS; sigErr; isExcl; sigTimes},{'sigVal'; 'sigTS'; 'sigErr'; 'isExcl'; 'sigTimes'});

%return

%%% Sanity checks %%% 

%% Figure 1: Plot decoding error descriptives

figure (1), clf
set(gcf,'Color','w')
tiledlayout(3,2); 

ax = sgtitle(['Excluded cue epochs: ' num2str(ixExcl')]);
ax.FontSize = 18;

nexttile
    smErr = movmean(Data(errBayes),100,'omitnan'); 
    vbl = rad2deg(smErr);
    xVal = Range(errBayes); 
    ax = plot(xVal,vbl,'b'); 
    hold on
    for nEp = 1:totEps
        plot([tsCue(nEp) tsCue(nEp)],[-120 120],'--r')
    end
    
    ax.LineWidth = 1;
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    xlim([min(xVal) max(xVal)])
   
    xlabel('Time (s)')
    ylabel('Decoding error (deg)')

nexttile
    vbl = errChange;
    bins = 1:totEps; 
    ax = bar(bins,vbl); 
    ax.FaceColor = [0.5 0.5 0.5];
    ax.EdgeColor = 'none';
    hold on 
    plot([0 17], [errTh errTh],'r')
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.XTick = 1:totEps;
    ylim([0 max(vbl)*1.2])
    ylabel('Mean error change (deg)')
    xlabel('No. cue epoch')

nexttile
    vbl = rad2deg(Data(errBayes));
    xVal = Range(errBayes); 
    ax = plot(xVal,vbl,'b'); 
    ax.LineWidth = 1;
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    xlim([min(xVal) max(xVal)])
    xlabel('Time (s)')
    ylabel('Decoding error (deg)')
    
nexttile
    vbl = rad2deg(varErr);
    bins = 1:totEps; 
    ax = bar(bins,vbl); 
    ax.FaceColor = [0.5 0.5 0.5];
    ax.EdgeColor = 'none';
    hold on 
    plot([0 17], [varTh varTh],'r')
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.XTick = 1:totEps;
    ylim([0 max(vbl)*1.2])
    ylabel('Circular variance (deg)')
    xlabel('No. cue epoch')
    
 nexttile
    vbl = abs(Data(Restrict(ahv,errBayes)));
    vbl = movmean(vbl,10);
    xVal = Range(errBayes);
    ax = plot(xVal,vbl,'b'); 
    ax.LineWidth = 1;
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    xlim([min(xVal) max(xVal)])
    xlabel('Time (s)')
    ylabel('Ahv (deg/s)')
 
nexttile
    vbl = epAhv;
    bins = 1:totEps; 
    ax = bar(bins,vbl); 
    ax.FaceColor = [0.5 0.5 0.5];
    ax.EdgeColor = 'none';
    hold on 
    plot([0 17], [ahvTh ahvTh],'r')
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.XTick = 1:totEps;
    ylim([0 max(vbl)*1.2])
    ylabel('Mean ahv (deg/s)')
    xlabel('No. cue epoch')
    
%% Figure 2: Plot all sigmoids (coarse res)

figure (2), clf
set(gcf,'Color','w')
tiledlayout(4,4)
    
for nEp = 1:totEps
    nexttile
    hold on
    err = rad2deg(sigErr(:,nEp));
    plot(sigTS(:,nEp),err,'b')
    ax = plot(sigTS(:,nEp), rad2deg(sigVal(:,nEp)),'r');
    ax.LineWidth = 1.5;
    ax = plot([tsCue(nEp) tsCue(nEp)],[min(err) max(err)],'--k');
    ax.LineWidth = 1.5;
    ax = plot([sigStart(nEp) sigStart(nEp)],[min(err) max(err)],'--m');
    ax.LineWidth = 1.5;
    ax = plot([sigEnd(nEp) sigEnd(nEp)],[min(err) max(err)],'--g');
    ax.LineWidth = 1.5;
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.YTick = [-180:30:180];
    ylim([min(err) max(err)])
    xlim([min(sigTS(:,nEp)) max(sigTS(:,nEp))])
    if sigDir(nEp) == -1
       ax.YDir = 'reverse';
    end
end


%% Figure 3: Plot all sigmoids (fine res)

figure (3), clf
set(gcf,'Color','w')
tiledlayout(4,4)
    
for nEp = 1:totEps
    nexttile
    sig = sigFine{nEp};
    hold on
    plot(sig(:,1),rad2deg(sig(:,3)),'b')
    ax = plot(sig(:,1), rad2deg(sig(:,2)),'r');
    ax.LineWidth = 1.5;
    ax = plot([tsCue(nEp) tsCue(nEp)],[min(err) max(err)],'--k');
    ax.LineWidth = 1.5;
    ax = plot([sigStartF(nEp) sigStartF(nEp)],[min(err) max(err)],'--m');
    ax.LineWidth = 1.5;
    ax = plot([sigEndF(nEp) sigEndF(nEp)],[min(err) max(err)],'--g');
    ax.LineWidth = 1.5;
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.YTick = [-180:30:180];
    ylim([min(err) max(err)])
    xlim([min(sig(:,1)) max(sig(:,1))])
    if sigDir(nEp) == -1
       ax.YDir = 'reverse';
    end
end


%% Figure 4: Evaluate all sigmoids

figure (4), clf
set(gcf,'Color','w')
tiledlayout(1,3)

nexttile
    vbl = meanSigDiffF;
    bins = 1:totEps; 
    ax = bar(bins,vbl); 
    ax.FaceColor = [0.5 0.5 0.5];
    ax.EdgeColor = 'none';
    hold on 
    plot([0 17], [sigTh sigTh],'r')
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.XTick = 1:totEps;
    ylim([0 max(vbl)*1.2])
    ylabel('Mean sigmoid difference (deg)')
    xlabel('No. cue epoch')
    
nexttile
    remapStart = sigStartF - tsCue;
    vbl = remapStart;
    bins = 1:totEps; 
    ax = bar(bins,vbl); 
    ax.FaceColor = [0.5 0.5 0.5];
    ax.EdgeColor = 'none';
    hold on 
    plot([0 17], [0 0],'r')
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.XTick = 1:totEps;
    %ylim([m max(vbl)*1.2])
    ylabel('Remaping delay (s)')
    xlabel('No. cue epoch')

    
nexttile                
    vbl = remapDur;
    bins = 1:totEps; 
    ax = bar(bins,vbl); 
    ax.FaceColor = [0.5 0.5 0.5];
    ax.EdgeColor = 'none';
    ax = gca; 
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.XTick = 1:totEps;
    %ylim([m max(vbl)*1.2])
    ylabel('Remaping duration (s)')
    xlabel('No. cue epoch')