function Duszkiewicz2024_Figure_GLMRemWake

% This script reporoduces panels from:
% Figure 6, Extended Data Figure 9

% Dependencies: 
% TStoolbox

%TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%% Parameters and variables

clear all 

smoothTC = 3; % smoothing of tuning curves (3 used in Duszkiewicz et al., 2024)

isHd = [];
isFs = [];
isGd = [];
isEx = [];
tcAll = [];
ccAll = [];
pairsTc = [];
pairsGlm = [];
rateWk = [];
rateRem = [];
betaRem = [];
betaWake = [];
warnRem = [];
warnWake = [];
remDur = [];

pairsGlm_rev = [];
betaRem_rev = [];
betaWake_rev = [];
warnRem_rev = [];
warnWake_rev = [];

betaRemP = [];
betaWakeP = [];
warnRemP = [];
warnWakeP = [];
pvalRemP = [];
pvalWakeP = [];

%% Load data
 
dataset = List2Cell('dataset_All_Rem.list');
cellsSoFar = 0;

for ii = 1:length(dataset)
    
    %Load data
    fprintf(['Uploading recording ', dataset{ii},'\n'])
    load(fullfile(dataset{ii},'Data','BehavEpochs'));
    load(fullfile(dataset{ii},'Data','SpikeData'));
    load(fullfile(dataset{ii},'Data','CellTypes'));
    load(fullfile(dataset{ii},'Analysis','HdTuning_moveEp'));
    load(fullfile(dataset{ii},'Analysis','GLM_RemWake'));
    load(fullfile(dataset{ii},'Analysis','GLM_RemWake_Rev'));
    load(fullfile(dataset{ii},'Analysis','TCxcorrs'));
    load(fullfile(dataset{ii},'Analysis','MeanFR'));
    load(fullfile(dataset{ii},'Sleep',[dataset{ii} '.SleepState.states.mat'] ));

    isHd = [isHd; hd];
    isFs = [isFs; fs];
    isGd = [isGd; gd];
    isEx = [isEx; ex];
    
    tcAll = [tcAll, hAll(:,:,smoothTC+1)];
    ccAll = [ccAll, ccgTC];
    pairsTc = [pairsTc; pairsTC+cellsSoFar]; 
    
    pairsGlm = [pairsGlm; pairs+cellsSoFar];
    betaRem = [betaRem beta_rem];
    betaWake = [betaWake beta_wake];
    warnRem = [warnRem warn_rem];
    warnWake = [warnWake warn_wake];   
    
    pairsGlm_rev = [pairsGlm_rev; pairs_rev+cellsSoFar];
    betaRem_rev = [betaRem_rev beta_rem_rev];
    betaWake_rev = [betaWake_rev beta_wake_rev];
    warnRem_rev = [warnRem_rev warn_rem_rev];
    warnWake_rev = [warnWake_rev warn_wake_rev];
       
    binsGlm = bins;
    
    rateWk = [rateWk; rateS];
    rateRem = [rateRem; rateREM]; 
        

    
    % Run GLM to calculate coupling of individual cells to pop activity
    
    Q = MakeQfromS(S,0.1);
    dQ = Data(Q);
    rQ = Range(Q);
    dQm = sum(dQ,2);
    dQ = gaussFilt(dQ,1,0); 
    dQm = gaussFilt(dQm,1,0);
    dQm = zscore(dQm); % regressor needs to be z-scored
    Qm = tsd(rQ,dQm);
    Q = tsd(rQ,dQ);
    
    totC = length(S);
    br = nan(totC,2);
    bw = nan(totC,2);
    wr = zeros(totC,1); % vector for catching glmfit warnings
    ww = zeros(totC,1); % vector for catching glmfit warnings
    pw = nan(totC,2);
    pr = nan(totC,2);

    % calculate population GLMs   
    rem = SleepState.ints.REMstate;
    remEp = intervalSet(rem(:,1),rem(:,2));
    remDur = [remDur; tot_length(remEp)];
    Qm_rem = Data(Restrict(Qm,remEp)); 
    Q_rem = Data(Restrict(Q,remEp));   
    Qm_wake = Data(Restrict(Qm,wake1Ep)); 
    Q_wake = Data(Restrict(Q,wake1Ep));

    for nC = 1:totC
        lastwarn('', ''); % clear last warning
        [b,~,stats] = glmfit(Qm_rem,Q_rem(:,nC),'poisson'); % run GLM
        [warnMsg, ~] = lastwarn(); % catch last warning
        if ~isempty(warnMsg)
            wr(nC) = 1;
        end
        br(nC,:) = b;
        pr(nC,:) = stats.p;
        
        lastwarn('', ''); % clear last warning
        [b,~,stats] = glmfit(Qm_wake,Q_wake(:,nC),'poisson'); % run GLM
        [warnMsg, ~] = lastwarn(); % catch last warning
        if ~isempty(warnMsg)
            ww(nC) = 1;
        end
        bw(nC,:) = b;
        pw(nC,:) = stats.p;
    end
    
    betaRemP = [betaRemP; br];
    betaWakeP = [betaWakeP; bw];
    warnRemP = [warnRemP ; wr];
    warnWakeP = [warnWakeP; ww];
    pvalRemP = [pvalRemP; pr];
    pvalWakeP = [pvalWakeP; pw];
     
    cellsSoFar = cellsSoFar+length(hd);
end

%% define groups

% Good cells 
ixGd = find(isGd & rateRem > 0.5 & rateWk > 0.5);

% HD cells
ixHd = find(isHd & isGd);

% FS cells
ixFs = find(isFs & isGd);

% Ex cells
ixEx = find(isEx & isGd);


%% Calculations

% Get angular differences between cells
    [~,diffAng] = max(ccAll,[],1);
    diffAng(diffAng >=180) = diffAng(diffAng >=180) - 360;
    diffAng = abs(diffAng'); 
    
 % Get 0-lag correlations between cells
    zeroCorr = ccAll(1,:);
    
% find GLMs with warnings 
    [~,colW] = find(warnWake == 1);
    [~,colR] = find(warnRem == 1);
    ixWarn = unique([colW; colR]); 
    
    [~,colW] = find(warnWake_rev == 1);
    [~,colR] = find(warnRem_rev == 1); 
    ixWarn_rev = unique([colW; colR]); 
    ixWarn = unique([ixWarn; ixWarn_rev]);
    
    isWarn = zeros(length(pairsGlm),1);
    isWarn(ixWarn) = 1; % ENABLE WARNING REMOVAL
    
% average two GLMs
    betaRem = (betaRem(:,:,3) + flipud(betaRem_rev(:,:,3)))./2; 
    betaWake = (betaWake(:,:,3) + flipud(betaWake_rev(:,:,3)))./2;   
    mid = round(size(betaWake,1)./2);
    
%%% Figures %%% 

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];
unc = [0.5 0.5 0.5];

fsc_light = [1 0.4 0.4];
fsc_dark = [0.6 0.1 0.1];

hdc_light = [0.7 0.45 1];
hdc_dark = [0.4 0.15 0.8];

wkc = [0.15 0.3 0.5];
rmc = [0.3 0.5 0.2];

fs = 11;
lw = 1;
%% Figure 1: Display HD cell GLM ccg's in WAKE and REM 

figure (1), clf
set(gcf,'color','w');

  
    % find HD-HD cell pairs and sort according to ang diff
    ix = find(ismember(pairsGlm(:,1),ixHd) & ismember(pairsGlm(:,2),ixHd) & ~isWarn);
     [~,sIx] = sort(diffAng(ix));
     %ix = flipud(ix(sIx)); 
     ix = ix(sIx); 
     
    cmap = redblue(256);
    colormap(cmap); 
    sdGauss1 = 2;
    sdGauss2 = 20;
    means = 1;

    
subplot(3,5,[1 2])
hold on
    vbl = betaWake(:,ix); 
    vbl(find(vbl < -1)) = -1;
    vbl(find(vbl > 1)) = 1;
    vbl = rescale(vbl,'InputMin',min(vbl),'InputMax',max(vbl));
    vbl = gaussFilt(vbl,sdGauss1,sdGauss2,means);
    ax = imagesc(flipud(vbl')); 
    ax = gca;
    xlim([0 201])
    ylim([0 size(vbl,2)])
    axis(ax, 'square');
    ax.Visible = 'off';
    ax = axes('Position',ax.Position);
    ax.Color = 'none';
    ax.YDir = 'reverse';
    axis(ax, 'square');
    xlim([-10 10])
    ylim([0 size(vbl,2)])
    ax.YTick = [1 size(vbl,2)];
    xlabel('Time lag (s)')
    ylabel('HD-HD cell pair number')
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YAxis.Exponent = 0;
    box on
    
    ax = title('WAKE');
    ax.FontSize = fs;
    ax.FontWeight = 'normal';

 subplot(3,5,[3 4])
    vbl = betaRem(:,ix); 
    vbl(find(vbl < -1)) = -1;
    vbl(find(vbl > 1)) = 1;
    vbl = rescale(vbl,'InputMin',min(vbl),'InputMax',max(vbl));
    vbl = gaussFilt(vbl,sdGauss1,sdGauss2,means);
    ax = imagesc(vbl');       
    ax = gca;
    xlim([0 201])
    ylim([0 size(vbl,2)])
    axis(ax, 'square');
    ax.Visible = 'off';
    ax = axes('Position',ax.Position);
    ax.Color = 'none';
    ax.YDir = 'reverse';
    axis(ax, 'square');
    xlim([-10 10])
    ylim([0 size(vbl,2)])
    ax.YTick = [];
    xlabel('Time lag (s)')
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YAxis.Exponent = 0;
    box on
    
    ax = title('REM');
    ax.FontSize = fs;
    ax.FontWeight = 'normal';
     
 subplot(3,5,5)   
   vbl = diffAng(ix);   
   ax = plot(vbl,1:length(vbl));   
   ax.LineWidth = 1;
   ax.Color = 'k';
   ax = gca;
   ax.YDir = 'reverse';
   ax.XDir = 'reverse';
   ax.XTick = [0 180];
   ax.YTick = [];
   ylim([0 length(vbl)])
   xlim([0 180])
   xlabel('Ang. diff. (deg)')
   box off
   ax.FontSize = fs;
   ax.LineWidth = lw;   
    
    

    % find FS-FS cell pairs and sort according to ang diff
    ix = find(ismember(pairsGlm(:,1),ixFs) & ismember(pairsGlm(:,2),ixFs) & ~isWarn);
    [~,sIx] = sort(zeroCorr(ix));
    ix = ix(fliplr(sIx)); 

    cmap = redblue(256);
    colormap(cmap); 
    sdGauss1 = 2;
    sdGauss2 = 20;

    
subplot(3,5,[6 7])
hold on
    vbl = betaWake(:,ix); 
    vbl(find(vbl < -1)) = -1;
    vbl(find(vbl > 1)) = 1;
    vbl = rescale(vbl,'InputMin',min(vbl),'InputMax',max(vbl));
    vbl = gaussFilt(vbl,sdGauss1,sdGauss2,means);
    ax = imagesc(flipud(vbl')); 
    ax = gca;
    xlim([0 201])
    ylim([0 size(vbl,2)])
    axis(ax, 'square');
    ax.Visible = 'off';
    ax = axes('Position',ax.Position);
    ax.Color = 'none';
    ax.YDir = 'reverse';
    axis(ax, 'square');
    xlim([-10 10])
    ylim([0 size(vbl,2)])
    ax.YTick = [1 size(vbl,2)];
    xlabel('Time lag (s)')
    ylabel('FS-FS cell pair number')
    ax.FontSize = fs;
    ax.LineWidth = lw;
    box on
    
    ax = title('WAKE');
    ax.FontSize = fs;
    ax.FontWeight = 'normal';

 subplot(3,5,[8 9])
    vbl = betaRem(:,ix); 
    vbl(find(vbl < -1)) = -1;
    vbl(find(vbl > 1)) = 1;
    vbl = rescale(vbl,'InputMin',min(vbl),'InputMax',max(vbl));
    vbl = gaussFilt(vbl,sdGauss1,sdGauss2,means);
    ax = imagesc(vbl');       
    ax = gca;
    xlim([0 201])
    ylim([0 size(vbl,2)])
    axis(ax, 'square');
    ax.Visible = 'off';
    ax = axes('Position',ax.Position);
    ax.Color = 'none';
    ax.YDir = 'reverse';
    axis(ax, 'square');
    xlim([-10 10])
    ylim([0 size(vbl,2)])
    ax.YTick = [];
    xlabel('Time lag (s)')
    ax.FontSize = fs;
    ax.LineWidth = lw;
    box on
    
    ax = title('REM');
    ax.FontSize = fs;
    ax.FontWeight = 'normal';
     
 subplot(3,5,10)   
   vbl = zeroCorr(ix);   
   ax = plot(vbl,1:length(vbl));   
   ax.LineWidth = 1;
   ax.Color = 'k';
   ax = gca;
   ax.YDir = 'reverse';
   ax.XTick = [-1 1];
   ax.YTick = [];
   ylim([0 length(vbl)])
   xlim([-1 1])
   xlabel('Tuning corr. (r)')
   box off
   ax.FontSize = fs;
   ax.LineWidth = lw;   

% find FS-HD cell pairs and sort according to ang diff
    ix1 = find(ismember(pairsGlm(:,1),ixHd) & ismember(pairsGlm(:,2),ixFs) & ~isWarn);
    ix2 = find(ismember(pairsGlm(:,1),ixFs) & ismember(pairsGlm(:,2),ixHd) & ~isWarn);
    ix = [ix1; ix2];
    [~,sIx] = sort(zeroCorr(ix));
     ix = flipud(ix(sIx)); 

    cmap = redblue(256);
    colormap(cmap); 
    sdGauss1 = 2;
    sdGauss2 = 20;

    
subplot(3,5,[11 12])
hold on
    vbl = betaWake(:,ix); 
    vbl(find(vbl < -1)) = -1;
    vbl(find(vbl > 1)) = 1;
    vbl = rescale(vbl,'InputMin',min(vbl),'InputMax',max(vbl));
    vbl = gaussFilt(vbl,sdGauss1,sdGauss2,means);
    ax = imagesc(flipud(vbl')); 
    ax = gca;
    xlim([0 201])
    ylim([0 size(vbl,2)])
    axis(ax, 'square');
    ax.Visible = 'off';
    ax = axes('Position',ax.Position);
    ax.Color = 'none';
    ax.YDir = 'reverse';
    axis(ax, 'square');
    xlim([-10 10])
    ylim([0 size(vbl,2)])
    ax.YTick = [1 size(vbl,2)];
    xlabel('Time lag (s)')
    ylabel('HD-FS cell pair number')
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YAxis.Exponent = 0;
    box on
    
    ax = title('WAKE');
    ax.FontSize = fs;
    ax.FontWeight = 'normal';

 subplot(3,5,[13 14])
    vbl = betaRem(:,ix); 
    vbl(find(vbl < -1)) = -1;
    vbl(find(vbl > 1)) = 1;
    vbl = rescale(vbl,'InputMin',min(vbl),'InputMax',max(vbl));
    vbl = gaussFilt(vbl,sdGauss1,sdGauss2,means);
    ax = imagesc(vbl');    
    ax = gca;
    xlim([0 201])
    ylim([0 size(vbl,2)])
    axis(ax, 'square');
    ax.Visible = 'off';
    ax = axes('Position',ax.Position);
    ax.Color = 'none';
    ax.YDir = 'reverse';
    axis(ax, 'square');
    xlim([-10 10])
    ylim([0 size(vbl,2)])
    ax.YTick = [];
    xlabel('Time lag (s)')
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YAxis.Exponent = 0;
    box on
    
    ax = title('REM');
    ax.FontSize = fs;
    ax.FontWeight = 'normal';
     
 subplot(3,5,15)   
   vbl = zeroCorr(ix);   
   ax = plot(vbl,1:length(vbl));   
   ax.LineWidth = 1;
   ax.Color = 'k';
   ax = gca;
   ax.YDir = 'reverse';
   ax.XTick = [-1 1];
   ax.YTick = [];
   ylim([0 length(vbl)])
   xlim([-1 1])
   xlabel('Tuning corr. (r)')
   box off
   ax.FontSize = fs;
   ax.LineWidth = lw;   
    

%% Figure 2: Correlations (REM vs WAKE)

figure (2), clf
set(gcf,'color','w');
tiledlayout(1,3);

nexttile
hold on
    ix1 = find(ismember(pairsGlm(:,1),ixHd) & ismember(pairsGlm(:,2),ixFs) & ~isWarn);
    ix2 = find(ismember(pairsGlm(:,1),ixFs) & ismember(pairsGlm(:,2),ixHd) & ~isWarn);
    ix = [ix1; ix2];
    mid = round(size(betaWake,1)./2);
    vbl1 = betaWake(mid,ix);       
    vbl2 = betaRem(mid,ix);          
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 1;
    ax.MarkerFaceColor = [0.5 0.5 0.5];
    ax.MarkerEdgeColor = 'none';
    Fit = polyfit(vbl1,vbl2,1);
    xVal = [min(vbl1),max(vbl1)];
    ax = plot(xVal,polyval(Fit,xVal),'-');
    ax.LineWidth = 0.5;
    ax.Color = 'r';
    ax = gca;  
    ax.LineWidth = 1;
    ax.FontSize = 11;
    ax.XTick = [-1 0 1];
    ax.YTick = [-1 0 1];
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylabel('Beta REM')
    xlabel('Beta WAKE')
    [r,p] = corr(vbl1',vbl2','type','Pearson');
    ax.FontSize = fs;
    ax.FontWeight = 'normal';
    xlim([-1 1])
    ylim([-1 1])
    ax = title([r; p]);
    ax.FontWeight = 'normal';
    
      
nexttile
hold on
    ix = find(ismember(pairsGlm(:,1),ixHd) & ismember(pairsGlm(:,2),ixHd) & ~isWarn);
    mid = round(size(betaWake,1)./2);
    vbl1 = betaWake(mid,ix);       
    vbl2 = betaRem(mid,ix);          
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 1;
    ax.MarkerFaceColor = [0.5 0.5 0.5];
    ax.MarkerEdgeColor = 'none';
    Fit = polyfit(vbl1,vbl2,1);
    xVal = [min(vbl1),max(vbl1)];
    ax = plot(xVal,polyval(Fit,xVal),'-');
    ax.LineWidth = 0.5;
    ax.Color = 'r';
    ax = gca;  
    ax.LineWidth = 1;
    ax.FontSize = 11;
    ax.XTick = [-1 0 1];
    ax.YTick = [-1 0 1];
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylabel('Beta REM')
    xlabel('Beta WAKE')
    ax.FontSize = fs;
    ax.FontWeight = 'normal';
    xlim([-1 1])
    ylim([-1 1])
    [r,p] = corr(vbl1',vbl2','type','Pearson');
    ax = title([r; p]);
    ax.FontWeight = 'normal';

nexttile
hold on
    ix = find(ismember(pairsGlm(:,1),ixFs) & ismember(pairsGlm(:,2),ixFs) & ~isWarn);
    mid = round(size(betaWake,1)./2);
    vbl1 = betaWake(mid,ix);       
    vbl2 = betaRem(mid,ix);          
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 1;
    ax.MarkerFaceColor = [0.5 0.5 0.5];
    ax.MarkerEdgeColor = 'none';
    Fit = polyfit(vbl1,vbl2,1);
    xVal = [min(vbl1),max(vbl1)];
    ax = plot(xVal,polyval(Fit,xVal),'-');
    ax.LineWidth = 0.5;
    ax.Color = 'r';
    ax = gca;  
    ax.LineWidth = 1;
    ax.FontSize = 11;
    ax.XTick = [-0.3 0 0.3];
    ax.YTick = [-0.3 0 0.3];
    ax.TickLength = [0.03 0.025];
    ylabel('Beta REM')
    xlabel('Beta WAKE')
    [r,p] = corr(vbl1',vbl2','type','Pearson');
    axis(ax, 'square');
    ax.FontSize = fs;
    ax.FontWeight = 'normal';
    xlim([-0.3 0.3])
    ylim([-0.3 0.3])
    ax = title([r; p]);
    ax.FontWeight = 'normal';
    


%% Figure 3: Correlations (beta vs tuning)


figure (3), clf
set(gcf,'color','w');
tiledlayout(2,2);

%2nd row 
nexttile
hold on
    ix1 = find(ismember(pairsGlm(:,1),ixHd) & ismember(pairsGlm(:,2),ixFs) & ~isWarn);
    ix2 = find(ismember(pairsGlm(:,1),ixFs) & ismember(pairsGlm(:,2),ixHd) & ~isWarn);
    ix = [ix1; ix2];  
    mid = round(size(betaWake,1)./2);
    vbl2 = betaWake(mid,ix);
    vbl1 = zeroCorr(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 1;
    ax.MarkerFaceColor = [0.5 0.5 0.5];
    ax.MarkerEdgeColor = 'none';
    Fit = polyfit(vbl1,vbl2,1);
    xVal = [min(vbl1),max(vbl1)];
    ax = plot(xVal,polyval(Fit,xVal),'-');
    ax.LineWidth = 0.5;
    ax.Color = 'r';
    ax = gca;
    ax.LineWidth = lw;
    ax.FontSize = fs; 
    ax.XTick = [-1 0 1];
    ax.YTick = [-1 0 1];
    ax.TickLength = [0.03 0.025];
    ylabel('Beta (WAKE)')
    xlabel('Tuning corr. (r)')
    xlim([-1 1]);
    ylim([-1 1]);
    axis(ax, 'square');
    [r,p] = corr(vbl1',vbl2','type','Pearson');
    ax = title([r; p]);
    ax.FontWeight = 'normal';
    
 
nexttile
hold on
    ix1 = find(ismember(pairsGlm(:,1),ixHd) & ismember(pairsGlm(:,2),ixFs) & ~isWarn);
    ix2 = find(ismember(pairsGlm(:,1),ixFs) & ismember(pairsGlm(:,2),ixHd) & ~isWarn);
    ix = [ix1; ix2];
    mid = round(size(betaRem,1)./2);
    vbl2 = betaRem(mid,ix);   
    vbl1 = zeroCorr(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 1;
    ax.MarkerFaceColor = [0.5 0.5 0.5];
    ax.MarkerEdgeColor = 'none';
    Fit = polyfit(vbl1,vbl2,1);
    xVal = [min(vbl1),max(vbl1)];
    ax = plot(xVal,polyval(Fit,xVal),'-');
    ax.LineWidth = 0.5;
    ax.Color = 'r';
    ax = gca;  
    ax.LineWidth = lw;
    ax.FontSize = fs;
    ax.XTick = [-1 0 1];
    ax.YTick = [-1 0 1];
    ax.TickLength = [0.03 0.025];
    ylabel('Beta REM')
    xlabel('Tuning corr. (r)')
    xlim([-1 1]);
    ylim([-1 1]);
    axis(ax, 'square');
    [r,p] = corr(vbl1',vbl2','type','Pearson');
    ax = title([r; p]);
    ax.FontWeight = 'normal';
    
%3rd row
nexttile
hold on
    ix = find(ismember(pairsGlm(:,1),ixFs) & ismember(pairsGlm(:,2),ixFs) & ~isWarn);
    mid = round(size(betaWake,1)./2);
    vbl2 = betaWake(mid,ix);
    vbl1 = zeroCorr(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 1;
    ax.MarkerFaceColor = [0.5 0.5 0.5];
    ax.MarkerEdgeColor = 'none';
    Fit = polyfit(vbl1,vbl2,1);
    xVal = [min(vbl1),max(vbl1)];
    ax = plot(xVal,polyval(Fit,xVal),'-');
    ax.LineWidth = 0.5;
    ax.Color = 'r';
    ax = gca;
    ax.LineWidth = lw;
    ax.FontSize = fs; 
    ax.YTick = [-0.3 0 0.3];
    ax.XTick = [-1 0 1];
    ax.TickLength = [0.03 0.025];
    ylabel('Beta (WAKE)')
    xlabel('Tuning corr. (r)')
    ylim([-0.3 0.3]);
    xlim([-1 1]);
    axis(ax, 'square');
    [r,p] = corr(vbl1',vbl2','type','Pearson');
    ax = title([r; p]);
    ax.FontWeight = 'normal';
  
nexttile
hold on
    ix = find(ismember(pairsGlm(:,1),ixFs) & ismember(pairsGlm(:,2),ixFs) & ~isWarn);
    mid = round(size(betaRem,1)./2);
    vbl2 = betaRem(mid,ix);   
    vbl1 = zeroCorr(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 1;
    ax.MarkerFaceColor = [0.5 0.5 0.5];
    ax.MarkerEdgeColor = 'none';
    Fit = polyfit(vbl1,vbl2,1);
    xVal = [min(vbl1),max(vbl1)];
    ax = plot(xVal,polyval(Fit,xVal),'-');
    ax.LineWidth = 0.5;
    ax.Color = 'r';
    ax = gca;  
    ax.LineWidth = lw;
    ax.FontSize = fs;
    ax.YTick = [-0.3 0 0.3];
    ax.XTick = [-1 0 1];
    ax.TickLength = [0.03 0.025];
    ylabel('Beta (REM)')
    xlabel('Tuning corr. (r)')
    axis(ax, 'square');
    ylim([-0.3 0.3]);
    xlim([-1 1]);
    [r,p] = corr(vbl1',vbl2','type','Pearson');
    ax = title([r; p]);
    ax.FontWeight = 'normal';
    
 %% Figure 4: GLM with population firing rate
 
figure (4), clf
set(gcf,'color','w');
tiledlayout(1,4);
 
isWarnP = warnWakeP == 1;
    
nexttile
hold on

    ix = find(isHd & ~isWarnP);
    vbl1 = betaWakeP(ix,2);
    m1 = median(vbl1);
    h1 = histogram(vbl1);
    h1.BinWidth = 0.1;
    h1.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h1.LineWidth = 2;
    h1.EdgeColor = hdc;
    h1.FaceColor = 'none';

    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = gca; 
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    xlim([-1 2]);
    ylim([0 0.2]);
    ylabel(['Prop. HD cells'])
    xlabel(['Beta(pop)'])
    ax = title('Wake');
    ax.FontWeight = 'normal';

    
nexttile
hold on

    ix = find(isHd & ~isWarnP);
    vbl1 = betaRemP(ix,2);
    m1 = median(vbl1);
    h1 = histogram(vbl1);
    h1.BinWidth = 0.1;
    h1.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h1.LineWidth = 2;
    h1.EdgeColor = hdc;
    h1.FaceColor = 'none';

    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = gca; 
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025]; 
    axis(ax, 'square');
    xlim([-1 2]);
    ylim([0 0.2])
    ylabel(['Prop. HD cells'])
    xlabel(['Beta(pop)'])
    ax = title('REM');
    ax.FontWeight = 'normal';
    
nexttile
hold on

    ix = find(isFs & ~isWarnP);
    vbl1 = betaWakeP(ix,2);
    m1 = median(vbl1);
    h1 = histogram(vbl1);
    h1.BinWidth = 0.1;
    h1.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h1.LineWidth = 2;
    h1.EdgeColor = fsc;
    h1.FaceColor = 'none';

    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = fsc;
    ax = gca; 
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025]; 
    axis(ax, 'square');
    xlim([-1 2]);
    ylim([0 0.4])
    ylabel(['Prop. FS cells'])
    xlabel(['Beta(pop)'])
    ax = title('Wake');
    ax.FontWeight = 'normal';
    
nexttile
hold on

    ix = find(isHd & ~isWarnP);
    vbl1 = betaRemP(ix,2);
    m1 = median(vbl1);
    h1 = histogram(vbl1);
    h1.BinWidth = 0.1;
    h1.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h1.LineWidth = 2;
    h1.EdgeColor = fsc;
    h1.FaceColor = 'none';

    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = fsc;
    ax = gca; 
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025]; 
    axis(ax, 'square');
    xlim([-1 2]);
    ylim([0 0.2])
    ylabel(['Prop. FS cells'])
    xlabel(['Beta(pop)'])
    ax = title('REM');
    ax.FontWeight = 'normal';
    


%% Figure 5: Wake/REM firing rates

figure (5), clf
set(gcf,'color','w');
tiledlayout(1,2);

nexttile
hold on

    vbl1 = rateWk(ixHd);
    vbl2 = rateRem(ixHd);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 8;
    ax.MarkerFaceColor = hdc;
    ax.MarkerEdgeColor = 'none';
    Fit = polyfit(vbl1,vbl2,1);
    xVal = [min(vbl1),max(vbl1)];
    ax = plot(xVal,polyval(Fit,xVal),'-');
    ax.LineWidth = 1;
    ax.Color = 'r';
    ax = plot([0 15],[0 15],'--k');
    ax.LineWidth = lw;
    ax = gca;
    ax.LineWidth = lw;
    ax.FontSize = fs; 
    %ax.YTick = [-0.3 0 0.3];
    %ax.XTick = [-1 0 1];
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylabel('Firing ratew REM (Hz)')
    xlabel('Firing rate WAKE (Hz)')
    ylim([0 15]);
    xlim([0 15]);
    [r,p] = corr(vbl1,vbl2,'type','Pearson');
    ax = title([r; p]);
    ax.FontWeight = 'normal';
    
nexttile
hold on
    vbl1 = rateWk(ixFs);
    vbl2 = rateRem(ixFs);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 8;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    Fit = polyfit(vbl1,vbl2,1);
    xVal = [min(vbl1),max(vbl1)];
    ax = plot(xVal,polyval(Fit,xVal),'-');
    ax.LineWidth = 1;
    ax.Color = 'r';
    ax = plot([0 100],[0 100],'--k');
    ax.LineWidth = lw;
    ax = gca;
    ax.LineWidth = lw;
    ax.FontSize = fs; 
    %ax.YTick = [-0.3 0 0.3];
    %ax.XTick = [-1 0 1];
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylabel('Firing ratew REM (Hz)')
    xlabel('Firing rate WAKE (Hz)')
    ylim([0 120]);
    xlim([0 120]);
    [r,p] = corr(vbl1,vbl2,'type','Pearson');
    ax = title([r; p]);
    ax.FontWeight = 'normal';



 %% Figure 6: GLM in Wake and REM - example HD-HD

vbl1 = betaWake;        
vbl1(:,ixWarn) = nan; % this bit is a walkaround to preserve cell indices
vbl2 = betaRem; 
vbl2(:,ixWarn) = nan; % this bit is a walkaround to preserve cell indices
vbl3 = circshift(ccAll,179);
vbl3(:,ixWarn) = nan; % this bit is a walkaround to preserve cell indices

pairIx = [41722 23876];
b = deg2rad(1:360)'; 

figure (6), clf
set(gcf,'color','w');
tiledlayout(length(pairIx),4); 
    
for nP = 1:length(pairIx)
        
    nexttile
        tc1 = tcAll(:,pairsGlm(pairIx(nP),1));
        tc1 = tc1./max(tc1);
        tc2 = tcAll(:,pairsGlm(pairIx(nP),2));
        tc2 = tc2./max(tc2);
        pAx = polarplot(b,tc1,'Color',hdc_light,'LineWidth',2);
        hold on
        pAx = polarplot(b,tc2,'Color',hdc_dark,'LineWidth',2);
            pAx = gca;
            pAx.ThetaDir = 'clockwise';
            pAx.ThetaZeroLocation = 'top';
            pAx.GridColor = 'k';
            thetaticks(0:90:270);  
            pAx.RColorMode = 'manual';
            pAx.RColor = hdc;
            pAx.ThetaColorMode = 'manual';
            pAx.ThetaColor = 'k';
            pAx.GridAlpha = 0.3;
            pAx.ThetaTickLabel = {'';'';'';''};
            rticks([]);    
            pAx.LineWidth = 3;
            
    nexttile
    hold on
    xVal = -179:180; 
    ax = plot(xVal,vbl3(:,pairIx(nP)));
        ax.LineWidth = 1;
        ax.Color = 'k';
        [~,maxIx] = max(vbl3(:,pairIx(nP)));
    ax = plot([xVal(maxIx) xVal(maxIx)],[-1 1],':k');
        ax.LineWidth = 3;
        ax = gca; 
        ax.XTick = [-180 0 180];
        ax.FontSize = fs;
        ax.LineWidth = lw; 
        axis(ax, 'square');
        ax.TickLength = [0.03 0.025];
        ylabel('Corr. coef. (r)')
        xlabel('Offset (deg)')
        xlim([-180 180])
        ylim([-1 1])  
        box off

    nexttile   
    xVal = -10:0.1:10;
    vbl = gaussFilt(vbl1(:,pairIx(nP)),2,0);
    ax = bar(xVal,vbl);
        ax.EdgeColor = wkc;
        ax.FaceColor = wkc;
        ax = gca;
        ax.FontSize = fs;
        ax.LineWidth = lw; 
        axis(ax, 'square');
        ax.TickLength = [0.03 0.025];
        ax.XTick = [-10 -5 0 5 10];
        if ~isnan(vbl1(:,nP))
            ylim([min(vbl1(:,pairIx(nP)))-0.1 max(vbl1(:,pairIx(nP)))+0.1])
        end
        xlabel('Offset ( s)')
        ylabel('Beta (Wake)')
        ylim([-1 1]) 
        box off
    
    nexttile   
    xVal = -10:0.1:10;
    vbl = gaussFilt(vbl2(:,pairIx(nP)),2,0);
    ax = bar(xVal,vbl);
        ax.EdgeColor = rmc;
        ax.FaceColor = rmc;
        ax = gca;
        ax.FontSize = fs;
        ax.LineWidth = lw; 
        axis(ax, 'square');
        ax.TickLength = [0.03 0.025];
        ax.XTick = [-10 -5 0 5 10];
        if ~isnan(vbl2(:,nP))
            ylim([min(vbl2(:,pairIx(nP)))-0.1 max(vbl2(:,pairIx(nP)))+0.1])
        end
        xlabel('Offset (s)')
        axis(ax, 'square');

        ylabel('Beta (REM)')
        ylim([-1 1]) 
        box off
            
end
    
    
 %% Figure 7: GLM in Wake and REM - example FS-FS

vbl1 = betaWake;        
vbl1(:,ixWarn) = nan; % this bit is a walkaround to preserve cell indices
vbl2 = betaRem; 
vbl2(:,ixWarn) = nan; % this bit is a walkaround to preserve cell indices
vbl3 = circshift(ccAll,179);
vbl3(:,ixWarn) = nan; % this bit is a walkaround to preserve cell indices

pairIx = [28005 12128];
b = deg2rad(1:360)'; 

figure (101), clf
set(gcf,'color','w');
tiledlayout(length(pairIx),4); 
    
for nP = 1:length(pairIx)
        
    nexttile
        tc1 = tcAll(:,pairsGlm(pairIx(nP),1));
        tc1 = tc1./max(tc1);
        tc2 = tcAll(:,pairsGlm(pairIx(nP),2));
        tc2 = tc2./max(tc2);
        pAx = polarplot(b,tc1,'Color',fsc_light,'LineWidth',2);
        hold on
        pAx = polarplot(b,tc2,'Color',fsc_dark,'LineWidth',2);
            pAx = gca;
            pAx.ThetaDir = 'clockwise';
            pAx.ThetaZeroLocation = 'top';
            pAx.GridColor = 'k';
            thetaticks(0:90:270);  
            pAx.RColorMode = 'manual';
            pAx.RColor = hdc;
            pAx.ThetaColorMode = 'manual';
            pAx.ThetaColor = 'k';
            pAx.GridAlpha = 0.3;
            pAx.ThetaTickLabel = {'';'';'';''};
            rticks([]);    
            pAx.LineWidth = 3;
            
    nexttile
    hold on
    xVal = -179:180; 
    ax = plot(xVal,vbl3(:,pairIx(nP)));
        ax.LineWidth = 1;
        ax.Color = 'k';
        [~,maxIx] = max(vbl3(:,pairIx(nP)));
    ax = plot([xVal(maxIx) xVal(maxIx)],[-1 1],':k');
        ax.LineWidth = 3;
        ax = gca; 
        ax.XTick = [-180 0 180];
        ax.FontSize = fs;
        ax.LineWidth = lw; 
        axis(ax, 'square');
        ax.TickLength = [0.03 0.025];
        ylabel('Corr. coef. (r)')
        xlabel('Offset (deg)')
        xlim([-180 180])
        ylim([-1 1])  
        box off

    nexttile   
    xVal = -10:0.1:10;
    vbl = gaussFilt(vbl1(:,pairIx(nP)),2,0);
    ax = bar(xVal,vbl);
        ax.EdgeColor = wkc;
        ax.FaceColor = wkc;
        ax = gca;
        ax.FontSize = fs;
        ax.LineWidth = lw; 
        axis(ax, 'square');
        ax.TickLength = [0.03 0.025];
        ax.XTick = [-10 -5 0 5 10];
        if ~isnan(vbl1(:,nP))
            ylim([min(vbl1(:,pairIx(nP)))-0.1 max(vbl1(:,pairIx(nP)))+0.1])
        end
        xlabel('Offset ( s)')
        ylabel('Beta (Wake)')
        ylim([-0.2 0.2]); 
        box off
    
    nexttile   
    xVal = -10:0.1:10;
    vbl = gaussFilt(vbl2(:,pairIx(nP)),2,0);
    ax = bar(xVal,vbl);
        ax.EdgeColor = rmc;
        ax.FaceColor = rmc;
        ax = gca;
        ax.FontSize = fs;
        ax.LineWidth = lw; 
        axis(ax, 'square');
        ax.TickLength = [0.03 0.025];
        ax.XTick = [-10 -5 0 5 10];
        if ~isnan(vbl2(:,nP))
            ylim([min(vbl2(:,pairIx(nP)))-0.1 max(vbl2(:,pairIx(nP)))+0.1])
        end
        xlabel('Offset (s)')
        axis(ax, 'square');

        ylabel('Beta (REM)')
        ylim([-0.2 0.2]);
        box off
            
end
    
 %% Figure 8: GLM in Wake and REM - example HD-FS

vbl1 = betaWake;        
vbl1(:,ixWarn) = nan; % this bit is a walkaround to preserve cell indices
vbl2 = betaRem; 
vbl2(:,ixWarn) = nan; % this bit is a walkaround to preserve cell indices
vbl3 = circshift(ccAll,179);
vbl3(:,ixWarn) = nan; % this bit is a walkaround to preserve cell indices

pairIx = [4760 64624];
b = deg2rad(1:360)'; 

figure (101), clf
set(gcf,'color','w');
tiledlayout(length(pairIx),4); 
    
for nP = 1:length(pairIx)
        
    nexttile
        tc1 = tcAll(:,pairsGlm(pairIx(nP),1));
        tc1 = tc1./max(tc1);
        tc2 = tcAll(:,pairsGlm(pairIx(nP),2));
        tc2 = tc2./max(tc2);
        pAx = polarplot(b,tc1,'Color',hdc,'LineWidth',2);
        hold on
        pAx = polarplot(b,tc2,'Color',fsc,'LineWidth',2);
            pAx = gca;
            pAx.ThetaDir = 'clockwise';
            pAx.ThetaZeroLocation = 'top';
            pAx.GridColor = 'k';
            thetaticks(0:90:270);  
            pAx.RColorMode = 'manual';
            pAx.RColor = hdc;
            pAx.ThetaColorMode = 'manual';
            pAx.ThetaColor = 'k';
            pAx.GridAlpha = 0.3;
            pAx.ThetaTickLabel = {'';'';'';''};
            rticks([]);    
            pAx.LineWidth = 3;
            
    nexttile
    hold on
    xVal = -179:180; 
    ax = plot(xVal,vbl3(:,pairIx(nP)));
        ax.LineWidth = 1;
        ax.Color = 'k';
        [~,maxIx] = max(vbl3(:,pairIx(nP)));
    ax = plot([xVal(maxIx) xVal(maxIx)],[-1 1],':k');
        ax.LineWidth = 3;
        ax = gca; 
        ax.XTick = [-180 0 180];
        ax.FontSize = fs;
        ax.LineWidth = lw; 
        axis(ax, 'square');
        ax.TickLength = [0.03 0.025];
        ylabel('Corr. coef. (r)')
        xlabel('Offset (deg)')
        xlim([-180 180])
        ylim([-1 1]); 
        box off

    nexttile   
    xVal = -10:0.1:10;
    vbl = gaussFilt(vbl1(:,pairIx(nP)),2,0);
    ax = bar(xVal,vbl);
        ax.EdgeColor = wkc;
        ax.FaceColor = wkc;
        ax = gca;
        ax.FontSize = fs;
        ax.LineWidth = lw; 
        axis(ax, 'square');
        ax.TickLength = [0.03 0.025];
        ax.XTick = [-10 -5 0 5 10];
        if ~isnan(vbl1(:,nP))
            ylim([min(vbl1(:,pairIx(nP)))-0.1 max(vbl1(:,pairIx(nP)))+0.1])
        end
        xlabel('Offset ( s)')
        ylabel('Beta (Wake)')
        ylim([-0.5 0.5]); 
        box off
    
    nexttile   
    xVal = -10:0.1:10;
    vbl = gaussFilt(vbl2(:,pairIx(nP)),2,0);
    ax = bar(xVal,vbl);
        ax.EdgeColor = rmc;
        ax.FaceColor = rmc;
        ax = gca;
        ax.FontSize = fs;
        ax.LineWidth = lw; 
        axis(ax, 'square');
        ax.TickLength = [0.03 0.025];
        ax.XTick = [-10 -5 0 5 10];
        if ~isnan(vbl2(:,nP))
            ylim([min(vbl2(:,pairIx(nP)))-0.1 max(vbl2(:,pairIx(nP)))+0.1])
        end
        xlabel('Offset (s)')
        axis(ax, 'square');

        ylabel('Beta (REM)')
        ylim([-0.5 0.5]); 
        box off
            
end
    
