function Duszkiewicz2024_Figure_Optogenetics

% This script reporoduces panels from:
% Figure 5, Extended Data Figures 7 and 8

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

offStart = 4; % start of OFF epoch (from beginning of light pulse)
offEnd = 0; % start of ON epoch (from beginning of light pulse)
smoothTC = 3; % smoothing of tuning curves for half-width (3 used in Duszkiewicz et al., 2024)
smoothFactors = 0; % smoothing of tuning curves for factor analysis (0 used in Duszkiewicz et al., 2024)

isHd = [];
isNhd = [];
isEx = [];
isFs = [];
isGd = [];
group = [];

tcAll_on = [];
tcAll_off = [];
tcAll_base = [];

tcAll60_on = [];
tcAll60_off = [];

frate_on = [];
frate_off = [];
frate_base = [];
allCCG_all = [];
allCCG_ind = [];
hdiOn = [];
hdiOff = [];
hdiBase = [];
rec = [];
isPos = [];
isAdn = [];
id = [];

fratePulseOn = [];
fratePulseOff = [];


%% Load data 
dataset = List2Cell('dataset_Opto.list');

% get animal ID
idRec = [];
for ii = 1:length(dataset)
    d = cell2mat(dataset(ii));
    idRec = [idRec; str2num(d(2:5))];
end

for ii = 1:length(dataset)
    
    fprintf(['Uploading recording ', dataset{ii},'\n'])
    load(fullfile(dataset{ii},'Analysis','BehavEpochs'));
    load(fullfile(dataset{ii},'Analysis','CellTypes'));
    load(fullfile(dataset{ii},'Analysis','SpikeData'));
    load(fullfile(dataset{ii},'Analysis','BrainArea'));
    load(fullfile(dataset{ii},'Analysis','Angle'));
    load(fullfile(dataset{ii},'Analysis','Groups'));
    load(fullfile(dataset{ii},'Analysis','Velocity'));
    optoEp = csvread(fullfile(dataset{ii},'Analysis','Opto_TS.csv'));     
    
    % collect variables
    isHd = [isHd; hd];
    isNhd = [isNhd; nhd];
    isFs = [isFs; fs];
    isGd = [isGd; gd];
    isEx = [isEx; ex];
    group = [group; g];
    isPos = [isPos; pos];
    isAdn = [isAdn; adn];
    id = [id; repmat(idRec(ii),[length(ex),1])];
    
    % get baseline epoch (no opto stim)
    baseEp = wakeEp;
    angrange = Range(Restrict(ang,baseEp));
    baseEp = intervalSet(angrange(1),angrange(end));
    frate_base = [frate_base; Rate(S,baseEp)];
    
    % restrict baseEp to velocity   
    velTh = 2;
    epVel = thresholdIntervals(vel,velTh,'Direction','Above');
    baseEp = intersect(baseEp,epVel);
    
    % calculate epochs  
    onEp    = intervalSet(optoEp(:,1), optoEp(:,2));
    offEp   = intervalSet(optoEp(:,1)-offStart,optoEp(:,1)-offEnd);    
    frate_on = [frate_on; Rate(S,onEp)];
    frate_off = [frate_off; Rate(S,offEp)];

    % calculate firing rates for each individual light pulses
    tstartOn = Start(onEp);
    tendOn = End(onEp);
    tstartOff = Start(offEp);
    tendOff = End(offEp);
    
    frOn = nan(length(tstartOn),length(S));
    frOff = nan(length(tstartOn),length(S));
    
    for nInt = 1:length(tstartOn)
        frOn(nInt,:) = Rate(S,intervalSet(tstartOn(nInt),tendOn(nInt)));
        frOff(nInt,:) = Rate(S,intervalSet(tstartOff(nInt),tendOff(nInt)));
    end
    
    fratePulseOn = [fratePulseOn frOn(1:240,:)];
    fratePulseOff = [fratePulseOff frOff(1:240,:)];
    
    % compute tuning curves and HD info
    
    % get regular tuning curves for width analysis
    for nC = 1:length(S)   
        [h1,b,mu,occ1] = HeadDirectionField(S{nC},ang,onEp,360,smoothTC);
        [h2,b,mu,occ2] = HeadDirectionField(S{nC},ang,offEp,360,smoothTC);    
        [h3,b,mu,occ3] = HeadDirectionField(S{nC},ang,baseEp,360,smoothTC); 
            
        tcAll_on = [tcAll_on h1(1:end-1)];
        tcAll_off = [tcAll_off h2(1:end-1)];  
        tcAll_base = [tcAll_base h3(1:end-1)];  
        
        if sum(h1) > 0 && sum(h2) > 0 && sum(h3) > 0
            hdiOn = [hdiOn; SpatialInfo(h1(1:end-1),occ1)];
            hdiOff = [hdiOff; SpatialInfo(h2(1:end-1),occ2)];
            hdiBase = [hdiBase; SpatialInfo(h3(1:end-1),occ3)];
        else
            hdiOn = [hdiOn; nan];
            hdiOff = [hdiOff; nan];
            hdiBase = [hdiBase; nan];
        end
    end
     
    % get downsampled tuning curves for correlations    
    for nC = 1:length(S)   
        [h1] = HeadDirectionField(S{nC},ang,onEp,60,smoothFactors);
        [h2] = HeadDirectionField(S{nC},ang,offEp,60,smoothFactors);    
        
        tcAll60_on = [tcAll60_on h1(1:end-1)];
        tcAll60_off = [tcAll60_off h2(1:end-1)];   
    end
    
    % calculate xcorr of each cell on opto stimulation (for average)
    tsOpto = ts(optoEp(:,1));
    binsize = 0.025; % binsize in miliseconds
    tot_time = 4; % total time in seconds (before+after)
    nbins_psth = tot_time./binsize; % Total number of bins for ccg 
    bins_all = -nbins_psth*binsize./2:binsize:nbins_psth*binsize./2;

     for nC = 1:length(S)
         c = Restrict(S{nC},wakeOptoEp);
         [ccg] = CrossCorr(tsOpto, c, binsize, nbins_psth);
         ccg = Data(ccg);
         allCCG_all = [allCCG_all, ccg];   
     end

    % calculate xcorr of each cell on opto stimulation (for individual cells)
    tsOpto = ts(optoEp(:,1));
    binsize = 0.05; % binsize in miliseconds
    tot_time = 5; % total time in seconds (before+after)
    nbins_psth = tot_time./binsize; % Total number of bins for ccg 
    bins_ind = -nbins_psth*binsize/2:binsize:nbins_psth*binsize/2;

     for nC = 1:length(S)
         c = Restrict(S{nC},wakeOptoEp);
         [ccg] = CrossCorr(tsOpto, c, binsize, nbins_psth);
         ccg = Data(ccg);
         allCCG_ind = [allCCG_ind, ccg];   
     end

end

totC = length(isFs);

%% define groups
% Good cells 
isGd = isGd == 1 & frate_on > 0.5 & frate_off > 0.5 & frate_base > 0.5;

% HD cells
isHd = isHd == 1 & isGd == 1;
ixHd = find(isHd);

% FS cells
isFs = isFs == 1 & isGd == 1;
ixFs = find(isFs);

% Ex cells
isEx = isEx == 1 & isGd == 1;
ixEx = find(isEx);

% NHD cells
isNhd = isNhd == 1 & isGd == 1;
ixNhd = find(isNhd);

%% Get significantly modulated or correlated cells 

% correlations 
[rTC,pvalTC] =corr(tcAll60_on,tcAll60_off);
rTC = diag(rTC);
pvalTC = diag(pvalTC);

% Good cells for gain analysis 
isGdGain = pvalTC < 0.05;

% select top half of most responsive cells
optoRes = frate_on./frate_off;
mResFs = median(optoRes(find(isFs & group == 1)));
mResPos = median(optoRes(find(isHd & isPos & group ==1)));
mResAdn = median(optoRes(find(isHd & isAdn)));
%% Compute half-width

% compute half-width (tuning curves)

tcAllSide_on = zeros(181,totC);
tcAllSide_off = zeros(181,totC);
tcAllSide_wk = zeros(181,totC);
hwTc_on = nan(totC,1);
hwTc_off = nan(totC,1);
hwTc_base = nan(totC,1);
hwTc_wk = nan(totC,1);

for nC = 1:totC
    
 [~,muMax] = max(tcAll_on(:,nC));  
  muMax = muMax(1);   
  if ~isnan(muMax)    
    hn = circshift(tcAll_on(:,nC),-muMax+1);
    hn = [hn./max(hn); 1];  
    hn = mean([hn(1:181),flipud(hn(181:361))],2);
    tcAllSide_on(:,nC) = hn;     
    hw = find(hn < 0.5,1,'first') * 2;
    if ~isempty(hw)
        hwTc_on(nC) = hw;
    end
  end
  
  [~,muMax] = max(tcAll_off(:,nC));  
  muMax = muMax(1);   
  if ~isnan(muMax)    
    hn = circshift(tcAll_off(:,nC),-muMax+1);
    hn = [hn./max(hn); 1];  
    hn = mean([hn(1:181),flipud(hn(181:361))],2);
    tcAllSide_off(:,nC) = hn;
    hw = find(hn < 0.5,1,'first') * 2;
    if ~isempty(hw)
        hwTc_off(nC) = hw;
    end
  end
  
end

%% Compute additive and multiplicative factors

% normalize
tcAllN_off = tcAll60_off./max(tcAll60_off,[],1);
tcAllN_on = tcAll60_on./max(tcAll60_off,[],1);

factorAdd_opto = nan(totC,1);
pAll_opto = nan(totC,2);

% for multiplication use non-normalized tuning curves 
for nC = 1:totC    
    pAll_opto(nC,:) = polyfit(tcAll60_off(:,nC), tcAll60_on(:,nC),1);   % slope
end
factorMult_opto = pAll_opto(:,1);

% for addition use normalized tuning curves
for nC = 1:totC    
    pAll_opto(nC,:) = polyfit(tcAllN_off(:,nC), tcAllN_on(:,nC),1);   % slope
    factorAdd_opto(nC) = polyval(pAll_opto(nC,:),0); % Y intercept
end

%% process PSTH

psthBase = mean(allCCG_all(1:nbins_psth/2,:),1);
psthProp = allCCG_all ./psthBase; 

      
%%% Figures %%%

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];
unc = [0.5 0.5 0.5];
adc = [1 0.4 1];

fsc_light = [1 0.4 0.4];
fsc_dark = [0.6 0.1 0.1];

hdc_light = [0.7 0.45 1];
hdc_dark = [0.4 0.15 0.8];

adc_light = [1 0.7 1];
adc_dark = [0.8 0.2 0.8];

lw = 1;
fs = 11;

%% Figure 1: Opto responses

figure (1), clf
set(gcf, 'Color', 'w');

%1st row
subplot(2,4,1)
hold on
    vbl = log10(hdiBase);
    ix1 = find(isHd & isAdn);
    ix2 = find(isNhd & isAdn);
    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = adc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    ax = plot([-0.5 -0.5],[0 30],'--');
    ax.LineWidth = 1;
    ax.Color = 'k';
    ax = gca; 
    xlim([-4 1]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. cells'])
    xlabel(['log10(HD info)'])
    
subplot(2,4,2)
hold on

    ix = find(isHd & isAdn & optoRes > mResAdn);
    ccg_mean = mean(psthProp(:,ix),2);
    ccg_sem = sem(psthProp(:,ix),2);
    boundedline(bins_all,ccg_mean,ccg_sem,'cmap',adc)
 
    ix = find(isFs & group == 1 & optoRes > mResFs);
    ccg_mean = mean(psthProp(:,ix),2);
    ccg_sem = sem(psthProp(:,ix),2);
    boundedline(bins_all,ccg_mean,ccg_sem,'cmap',fsc)   
      
    ix = find(isHd & isPos & group == 1 & optoRes > mResPos);
    ccg_mean = mean(psthProp(:,ix),2);
    ccg_sem = sem(psthProp(:,ix),2);
    boundedline(bins_all,ccg_mean,ccg_sem,'cmap',hdc)   
       
    ax = gca;
    ylabel('Firing rate (fold change)')
    xlabel('Time from light onset (s) ')
    ax = gca;
    axis(ax, 'square');
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ylim([0 3]);
    xlim([-1 2]);
 
subplot(2,4,3)
hold on

    ix = find(isHd & isAdn & optoRes > mResAdn);
    ccg_mean = mean(psthProp(:,ix),2);
    ccg_sem = sem(psthProp(:,ix),2);
    boundedline(bins_all,ccg_mean,ccg_sem,'cmap',adc)  
      
    ix = find(isHd & isPos & group == 1 & optoRes > mResPos);
    ccg_mean = mean(psthProp(:,ix),2);
    ccg_sem = sem(psthProp(:,ix),2);
    boundedline(bins_all,ccg_mean,ccg_sem,'cmap',hdc) 
    
    ix = find(isFs & group == 1 & optoRes > mResFs);
    ccg_mean = mean(psthProp(:,ix),2);
    ccg_sem = sem(psthProp(:,ix),2);
    boundedline(bins_all,ccg_mean,ccg_sem,'cmap',fsc) 
       
    ax = gca;
    ylabel('Firing rate (fold change)')
    xlabel('Time from light onset (s) ')
    ax = gca;
    axis(ax, 'square');
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ylim([0 3]);
    xlim([-0.05 0.1]);
    
% 2nd row
subplot(2,4,5)
hold on

    ix = find(isHd & isAdn);
    xVal = 0:180;
    
    for nC = 1:length(ix)
        ax = plot(xVal,tcAllSide_off(:,ix(nC)));
        ax.Color = [0.7 0.7 0.7];
        ax.LineWidth = 0.5;
    end
      
    
    ax = plot(xVal,mean(tcAllSide_off(:,ix),2));
    ax.Color = adc_dark;
    ax.LineWidth = 2;
        
    ax = gca;
    ylabel('Firing rate (fold change)')
    xlabel('Time from light onset (s) ')
    ax = gca;
    axis(ax, 'square');
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ylim([0 1]);
    xlim([0 90]);
    ax.XTick = [0 30 60 90];
    
subplot(2,4,6)
hold on

    ix = find(isHd & isAdn);
    xVal = 0:180;
    
    for nC = 1:length(ix)
        ax = plot(xVal,tcAllSide_on(:,ix(nC)));
        ax.Color = [0.7 0.7 0.7];
        ax.LineWidth = 0.5;
    end
    
    ax = plot(xVal,mean(tcAllSide_off(:,ix),2));
    ax.Color = adc_light;
    ax.LineWidth = 2;
       
    ax = gca;
    ylabel('Firing rate (fold change)')
    xlabel('Time from light onset (s) ')
    ax = gca;
    axis(ax, 'square');
    ax.TickLength = [0.03 0.025];
    ax.XTick = [0 180];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ylim([0 1]);
    xlim([0 90]);
    ax.XTick = [0 30 60 90];
    
subplot(2,4,7)
hold on

    ix = find(isHd & isPos & group == 1);
    xVal = 0:180;
    
    for nC = 1:length(ix)
        ax = plot(xVal,tcAllSide_off(:,ix(nC)));
        ax.Color = [0.7 0.7 0.7];
        ax.LineWidth = 0.5;
    end
      
    
    ax = plot(xVal,mean(tcAllSide_off(:,ix),2));
    ax.Color = hdc_light;
    ax.LineWidth = 2;
        
    ax = gca;
    ylabel('Firing rate (fold change)')
    xlabel('Time from light onset (s) ')
    ax = gca;
    axis(ax, 'square');
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ylim([0 1]);
    xlim([0 90]);
    ax.XTick = [0 30 60 90];
    
subplot(2,4,8)
hold on

    ix = find(isHd & isPos & group == 1);
    xVal = 0:180;
    
    for nC = 1:length(ix)
        ax = plot(xVal,tcAllSide_on(:,ix(nC)));
        ax.Color = [0.7 0.7 0.7];
        ax.LineWidth = 0.5;
    end
    
    ax = plot(xVal,mean(tcAllSide_off(:,ix),2));
    ax.Color = hdc_dark;
    ax.LineWidth = 2;
       
    ax = gca;
    ylabel('Firing rate (fold change)')
    xlabel('Time from light onset (s) ')
    ax = gca;
    axis(ax, 'square');
    ax.TickLength = [0.03 0.025];
    ax.XTick = [0 180];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    xlim([0 90]);
    ax.XTick = [0 30 60 90];
    


%% Figure 2: Opto responses continued

figure (2), clf
set(gcf, 'Color', 'w');
    
subplot(2,3,1)
hold on

    ix1 = find(isHd & isAdn);
    ix2 = find(isNhd & isAdn);

    vbl1 = optoRes(ix1); 
    vbl2 = optoRes(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    %h1.Normalization = 'probability';
    %h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = adc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = adc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    xlim([0 4]);
    ylim([0 30]);
    ax.TickLength = [0.03 0.025];
    ax.XTick = [0 1 2 3 4];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. cells'])
    xlabel(['Fold change in firing rate'])
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
subplot(2,3,2)
hold on

    vbl = optoRes;
    ix1 = find(isHd == 1 & group == 1 & isPos);
    ix2 = find(isHd == 1 & group == 2 & isPos);

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    %h1.Normalization = 'probability';
    %h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = hdc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    xlim([0 3]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. HD cells'])
    xlabel(['Fold change in firing rate'])
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';

 
subplot(2,3,3)
hold on
    vbl = optoRes;
    ix1 = find(isFs == 1 & group == 1 & isPos);
    ix2 = find(isFs == 1 & group == 2 & isPos);

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    %h1.Normalization = 'probability';
    %h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = fsc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = fsc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    xlim([0 3]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. FS cells'])
    xlabel(['Fold change in firing rate'])
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
subplot(2,3,4)
hold on

    ix = find(isHd & isAdn);
    vbl1 = hwTc_off(ix);
    vbl2 = hwTc_on(ix);
    
    ax = plot([1 2],[vbl1,vbl2],'-','Color',adc);
    ax = plot([0.8 1.2],[mean(vbl1),mean(vbl1)],'-','Color',[0 0 0]);
    ax.LineWidth = lw;
    ax = plot([1.8 2.2],[mean(vbl2),mean(vbl2)],'-','Color',[0 0 0]);
    ax.LineWidth = lw;
    
    ax = gca; 
    ax.XTick = [1 2];
    ax.XTickLabel = {'OFF';'ON'};
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    xlim([0.5 2.5]);
    ylim([0 150]);
    ylabel(['Tuning curve width (deg)'])   
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
subplot(2,3,5)
hold on

    ix  = find(isHd & group == 1 & isPos);
    vbl1 = hwTc_off(ix);
    vbl2 = hwTc_on(ix);
    
    ax = plot([1 2],[vbl1,vbl2],'-','Color',hdc);
    ax = plot([0.8 1.2],[mean(vbl1),mean(vbl1)],'-','Color',[0 0 0]);
    ax.LineWidth = lw;
    ax = plot([1.8 2.2],[mean(vbl2),mean(vbl2)],'-','Color',[0 0 0]);
    ax.LineWidth = lw;
    
    ax = gca; 
    ax.XTick = [1 2];
    ax.XTickLabel = {'OFF';'ON'};
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    xlim([0.5 2.5]);
    ylim([0 150]);
    ylabel(['Tuning curve width (deg)'])   
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';

subplot(2,3,6)
hold on
    vbl = rTC;
    ix1 = isFs == 1 & group == 1;
    ix2 = isFs == 1 & group == 2;

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.1;
    h2.BinWidth = 0.1;
    %h1.Normalization = 'probability';
    %h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = fsc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = fsc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    xlim([-1 1]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. FS cells'])
    xlabel(['Correlation (r)'])
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';


%% Figure 3: Additive vs multiplicative factors (opto)

figure (3), clf
set(gcf, 'Color', 'w');
   
% 1st row    
subplot(3,4,1)
hold on
    vbl = factorAdd_opto;
    ix1 = isFs == 1 & group == 1 & isPos & isGdGain;
    ix2 = isFs == 1 & group == 2;

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.1;
    h2.BinWidth = 0.1;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = fsc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = fsc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    xlim([-0.5 1.1]);
    ylim([0 0.6]);
    ax.TickLength = [0.03 0.025];
    %ax.YTick = [0:0.1:0.5];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. FS cells'])
    xlabel(['Additive factor'])
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
subplot(3,4,2)
hold on
    vbl = factorAdd_opto;
    ix1 = isHd == 1 & group == 1 & isPos & isGdGain;
    ix2 = isHd == 1 & group == 2;

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.05;
    h2.BinWidth = 0.05;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = hdc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    xlim([-0.5 1]);
    ylim([0 1]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. HD cells'])
    xlabel(['Additive factor'])
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
   
subplot(3,4,3)
hold on 

    vbl = factorAdd_opto;
    ix1 = isHd == 1 & group == 1 & isAdn & isGdGain;
    ix2 = isHd == 1 & group == 1 & isPos & isGdGain;

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.02;
    h2.BinWidth = 0.02;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = adc;
    h1.FaceColor = 'none';
    h2.EdgeColor = hdc;
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = adc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = gca; 
    xlim([-0.15 0.3]);
    ylim([0 0.6]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. HD cells'])
    xlabel(['Additive factor'])
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
subplot(3,4,4)
hold on
    
    
    ix1 = isHd == 1 & group == 1 & isPos & isGdGain;
    vbl1 = factorMult_opto;
    vbl2 = factorAdd_opto;
    ax = scatter(vbl1(ix1),vbl2(ix1));
    ax.SizeData = 7;
    ax.MarkerFaceColor = hdc;
    ax.MarkerEdgeColor = 'none';
    ax = plot([0 3], [0 0], ':k');
    ax.LineWidth = 1.5;
    ax = plot([1 1], [-1 1], ':k');
    ax.LineWidth = 1.5;
    [r,p] = corr(vbl1,vbl2,'type','Pearson');
    
    ix2 = isFs == 1 & group == 1 & isPos & isGdGain;
    ax = scatter(vbl1(ix2),vbl2(ix2));
    ax.SizeData = 7;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    ax = plot([0 3], [0 0], ':k');
    ax.LineWidth = 1.5;
    ax = plot([1 1], [-1 1], ':k');
    ax.LineWidth = 1.5;
    [r,p] = corr(vbl1,vbl2,'type','Pearson');
    
    ax = gca;  
    ax.LineWidth = lw;
    ax.FontSize = fs;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    xlabel({'Multiplicative factor'})
    ylabel('Additive factor')
    xlim([0 3])
    ylim([-0.5 1])

    [p,h,stats] = ranksum(vbl1(ix1),vbl1(ix2));
    [p,h,stats] = ranksum(vbl2(ix1),vbl2(ix2))
    
% 2nd row
subplot(3,4,5)
hold on
    vbl = factorMult_opto;
    ix1 = isFs == 1 & group == 1 & isPos & isGdGain;
    ix2 = isFs == 1 & group == 2;

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = fsc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = fsc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    xlim([-0.5 3]);
    ylim([0 0.8]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. FS cells'])
    xlabel(['Multiplicative factor'])
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
subplot(3,4,6)
hold on
    vbl = factorMult_opto;
    ix1 = isHd == 1 & group == 1 & isPos & isGdGain;
    ix2 = isHd == 1 & group == 2 & isPos;

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = hdc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    xlim([-0.5 3]);
    ylim([0 0.6]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. HD cells'])
    xlabel(['Multiplicative factor'])
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';

subplot(3,4,7)
hold on
    vbl = factorMult_opto;
    ix1 = isHd == 1 & group == 1 & isAdn & isGdGain;
    ix2 = isHd == 1 & group == 1 & isPos & isGdGain;

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = adc;
    h1.FaceColor = 'none';
    h2.EdgeColor = hdc;
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = adc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = gca; 
    xlim([-0.5 3]);
    ylim([0 0.6]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. HD cells'])
    xlabel(['Multiplicative factor'])
    [p,h,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
subplot(3,4,8)
hold on
        
    ix = isHd == 1 & group == 1 & isPos & isGdGain;
    vbl1 = factorMult_opto(ix);
    vbl2 = factorAdd_opto(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 7;
    ax.MarkerFaceColor = hdc;
    ax.MarkerEdgeColor = 'none';
    ax = plot([0 3], [0 0], ':k');
    ax.LineWidth = 1.5;
    ax = plot([1 1], [-1 1], ':k');
    ax.LineWidth = 1.5;
    [r,p] = corr(vbl1,vbl2,'type','Pearson');

    ix = isHd == 1 & group == 1 & isAdn & isGdGain;
    vbl1 = factorMult_opto(ix);
    vbl2 = factorAdd_opto(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 7;
    ax.MarkerFaceColor = adc;
    ax.MarkerEdgeColor = 'none';
    ax = plot([0 3], [0 0], ':k');
    ax.LineWidth = 1.5;
    ax = plot([1 1], [-1 1], ':k');
    ax.LineWidth = 1.5;
    [r,p] = corr(vbl1,vbl2,'type','Pearson');
    
    ax = gca;  
    ax.LineWidth = lw;
    ax.FontSize = fs;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    xlabel({'Multiplicative factor'})
    ylabel('Additive factor')
    %xlim([0.5 1.5])
    ylim([-0.2 0.2])

%3rd row

subplot(3,4,4)
hold on
    
    
    ix = isFs == 1 & group == 1 & isPos & isGdGain;
    vbl1 = factorMult_opto(ix);
    vbl2 = factorAdd_opto(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 7;
    ax.MarkerFaceColor = hdc;
    ax.MarkerEdgeColor = 'none';
    ax = plot([0 3], [0 0], ':k');
    ax.LineWidth = 1.5;
    ax = plot([1 1], [-1 1], ':k');
    ax.LineWidth = 1.5;
    [r,p] = corr(vbl1,vbl2,'type','Pearson');
    
    ix = isFs == 1 & group == 1 & isPos & isGdGain;
    vbl1 = factorMult_opto(ix);
    vbl2 = factorAdd_opto(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 7;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    ax = plot([0 3], [0 0], ':k');
    ax.LineWidth = 1.5;
    ax = plot([1 1], [-1 1], ':k');
    ax.LineWidth = 1.5;
    [r,p] = corr(vbl1,vbl2,'type','Pearson');
    
    ax = gca;  
    ax.LineWidth = lw;
    ax.FontSize = fs;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    xlabel({'Multiplicative factor'})
    ylabel('Additive factor')
    xlim([0 3])
    ylim([-0.5 1])
    
%% Figure 4: Individual examples ADN-Nhd

figure (4), clf
set(gcf, 'Color', 'w');

ixC = 54;

ix = find(isNhd & isAdn);
totc = length(ixC); 
tiledlayout(totc,5);

for nC = 1:totc
    
    nexttile 
        tc1 = tcAll_on(:,ix(ixC(nC)));  
        tc2 = tcAll_off(:,ix(ixC(nC)));

        ax = polarplot(b,[tc1; tc1(1)]);
        hold on
            ax.Color = 'k';
            ax.LineWidth = 2;

        ax = polarplot(b,[tc2; tc2(1)]);
            ax.Color = unc;
            ax.LineWidth = 2;
            ax = gca;
            ax.ThetaDir = 'clockwise';
            ax.ThetaZeroLocation = 'bottom';
            ax.GridColor = 'k';
            thetaticks(0:90:270);
            ax.RColorMode = 'manual';
            ax.RColor = 'k';
            ax.ThetaColorMode = 'manual';
            ax.ThetaColor = 'k';
            ax.GridAlpha = 0.3;
            ax.ThetaTickLabel = {'';'';'';''};
            rticks([]);  
            ax.LineWidth = 2; 
            
    nexttile   
    hold on   
        tc1 = tcAll_on(:,ix(ixC(nC)));  
        tc2 = tcAll_off(:,ix(ixC(nC))); 
        tc1 = tc1./max(tc2);
        tc2 = tc2./max(tc2);

        xVal = 1:360;
        yVal = [tc1,tc2];
        ax = plot(xVal,yVal(:,1));
        ax.Color = 'k';
        ax.LineWidth = 1;
        ax = plot(xVal,yVal(:,2));
        ax.Color = unc;
        ax.LineWidth = 1;
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.XTick = [0 360];
        xlim([0 360])
        ylim([0 2])
        ax.TickLength = [0.03 0.025];
        xlabel('HD')
        ylabel('Firing rate')
        axis(ax, 'square');  


    nexttile  
    hold on   
        psthBase = mean(allCCG_ind(1:nbins_psth/2,:),1);
        psthProp = allCCG_ind ./psthBase; 
        ccg = psthProp(:,ix(ixC(nC)));
        ax = plot(bins_ind,ccg);
        ax.LineWidth = lw;
        ax.Color = 'k';

        ylabel('Firing rate (fold change)')
        xlabel('Time from light onset (s) ')
        ax = gca;
        axis(ax, 'square');
        ax.FontSize = fs;
        ax.LineWidth = lw;
        ax.TickLength = [0.03 0.025];
        ylim([0 3]);
        xlim([-1 2]);
  

    nexttile  
    hold on    
        tc1 = tcAll_on(:,ix(ixC(nC)));  
        tc2 = tcAll_off(:,ix(ixC(nC)));  
        xVal = -179:180;
        cc = TCcrosscorr(tc1,tc2);
        cc = circshift(cc,179);
        ax = plot(xVal,cc);
        ax.Color = 'k';
        ax.LineWidth = 1;
        [~,maxIx] = max(cc);
        maxIx = maxIx - 179;
        ax = plot([maxIx maxIx],[-1 1],':');
        ax.LineWidth = 2;
        ax.Color = 'k';
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.XTick = [-180 -90 0 90 180];
        ax.XTickLabel = [{'-180'},{''},{'0'},{''},{'180'}];
        xlim([-180 180])
        ax.TickLength = [0.03 0.025];
        xlabel('Offset (deg)')
        ylabel('Correlation (r)')
        axis(ax, 'square');  
    
    nexttile 
    hold on   
        tc1 = tcAll60_on(:,ix(ixC(nC)));  
        tc2 = tcAll60_off(:,ix(ixC(nC))); 
        tc1 = tc1./max(tc2);
        tc2 = tc2./max(tc2);  
        ax = scatter(tc2,tc1,'.');
        ax.SizeData = 50;
        ax.MarkerFaceColor = [0.5 0.5 0.5];
        ax.MarkerEdgeColor = [0.5 0.5 0.5];
        
        Fit = polyfit(tc2,tc1,1);
        xVal = [0 1];
        ax = plot(xVal,polyval(Fit,xVal),'-r');
        ax.LineWidth = lw;

        ax = plot([0 1],[0 1],'--k');
        ax.LineWidth = lw;
               
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.YTick = [0 1 2];
        ax.XTick = [0 1];
        xlim([0 1])
        ylim([0 2])
        ax.TickLength = [0.03 0.025];
        xlabel('Tuning OFF')
        ylabel('Tuning ON')
        axis(ax, 'square'); 
end
    
%% Figure 5: Individual examples ADN-HD

figure (5), clf
set(gcf, 'Color', 'w');

ixC = 37; 

ix = find(isHd & isAdn);
totc = length(ixC);    
tiledlayout(totc,5)

for nC = 1:totc
    
    nexttile   
        tc1 = tcAll_on(:,ix(ixC(nC)));  
        tc2 = tcAll_off(:,ix(ixC(nC)));

        ax = polarplot(b,[tc1; tc1(1)]);
        hold on
            ax.Color = adc_dark;
            ax.LineWidth = 2;

        ax = polarplot(b,[tc2; tc2(1)]);
            ax.Color = adc_light;
            ax.LineWidth = 2;
            ax = gca;
            ax.ThetaDir = 'clockwise';
            ax.ThetaZeroLocation = 'bottom';
            ax.GridColor = [0.7 0.7 0.7];
            thetaticks(0:90:270);
            ax.RColorMode = 'manual';
            ax.RColor = [0.7 0.7 0.7];
            ax.ThetaColorMode = 'manual';
            ax.ThetaColor = [0.7 0.7 0.7];
            ax.GridAlpha = 1;
            ax.ThetaTickLabel = {'';'';'';''};
            rticks([]);  
            ax.LineWidth = 2; 
            
    nexttile  
    hold on    
        tc1 = tcAll_on(:,ix(ixC(nC)));  
        tc2 = tcAll_off(:,ix(ixC(nC))); 
        [~,mIx] = max(tc1);
        tc1 = circshift(tc1,-mIx+180);
        tc2 = circshift(tc2,-mIx+180);
        tc1 = tc1./max(tc2);
        tc2 = tc2./max(tc2);

        xVal = -179:180;
        yVal = [tc1,tc2];
        ax = plot(xVal,yVal(:,1));
        ax.Color = adc_dark;
        ax.LineWidth = 1;
        ax = plot(xVal,yVal(:,2));
        ax.Color = adc_light;
        ax.LineWidth = 1;
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.XTick = [-180 0 180];
        xlim([-180 180])
        ylim([0 2])
        ax.TickLength = [0.03 0.025];
        xlabel('Angle from peak (deg)')
        ylabel('Firing rate')
        axis(ax, 'square');  


    nexttile 
    hold on    
        psthBase = mean(allCCG_ind(1:nbins_psth/2,:),1);
        psthProp = allCCG_ind ./psthBase; 
        ccg = psthProp(:,ix(ixC(nC)));
        ax = plot(bins_ind,ccg);
        ax.LineWidth = lw;
        ax.Color = 'k';

        ylabel('Firing rate (fold change)')
        xlabel('Time from light onset (s) ')
        ax = gca;
        axis(ax, 'square');
        ax.FontSize = fs;
        ax.LineWidth = lw;
        ax.TickLength = [0.03 0.025];
        ylim([0 3]);
        xlim([-1 2]);
  
    nexttile 
    hold on   
        tc1 = tcAll_on(:,ix(ixC(nC)));  
        tc2 = tcAll_off(:,ix(ixC(nC)));  
        xVal = -179:180;
        cc = TCcrosscorr(tc1,tc2);
        cc = circshift(cc,179);
        ax = plot(xVal,cc);
        ax.Color = 'k';
        ax.LineWidth = 1;
        [~,maxIx] = max(cc);
        maxIx = maxIx - 179;
        ax = plot([maxIx maxIx],[-1 1],':');
        ax.LineWidth = 2;
        ax.Color = 'k';
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.XTick = [-180 -90 0 90 180];
        ax.XTickLabel = [{'-180'},{''},{'0'},{''},{'180'}];
        xlim([-180 180])
        ax.TickLength = [0.03 0.025];
        xlabel('Offset (deg)')
        ylabel('Correlation (r)')
        axis(ax, 'square');  
    
    nexttile 
    hold on    
        tc1 = tcAll60_on(:,ix(ixC(nC)));  
        tc2 = tcAll60_off(:,ix(ixC(nC))); 
        tc1 = tc1./max(tc2);
        tc2 = tc2./max(tc2);  
        ax = scatter(tc2,tc1,'.');
        ax.SizeData = 50;
        ax.MarkerFaceColor = [0.5 0.5 0.5];
        ax.MarkerEdgeColor = [0.5 0.5 0.5];
        
        Fit = polyfit(tc2,tc1,1);
        xVal = [0 1];
        ax = plot(xVal,polyval(Fit,xVal),'-r');
        ax.LineWidth = lw;

        ax = plot([0 1],[0 1],'--k');
        ax.LineWidth = lw;
                
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.YTick = [0 1 2];
        ax.XTick = [0 1];
        xlim([0 1])
        ylim([0 2])
        ax.TickLength = [0.03 0.025];
        xlabel('Tuning OFF')
        ylabel('Tuning ON')
        axis(ax, 'square'); 
end       
  
    
%% Figure 6: Individual examples FACTORS (HD OPTO)

figure (6), clf
set(gcf, 'Color', 'w');

ixC = [3 82]; 

totc = length(ixC);  
ix = find(isHd & isPos & group == 1);
tiledlayout(totc,5)

for nC = 1:totc
    
    nexttile  
        tc1 = tcAll_on(:,ix(ixC(nC)));  
        tc2 = tcAll_off(:,ix(ixC(nC)));
        tc1 = gaussFiltAng(tc1,3);
        tc2 = gaussFiltAng(tc2,3);

        ax = polarplot(b,[tc1; tc1(1)]);
        hold on
            ax.Color = hdc_dark;
            ax.LineWidth = 2;

        ax = polarplot(b,[tc2; tc2(1)]);
            ax.Color = hdc_light;
            ax.LineWidth = 2;
            ax = gca;
            ax.ThetaDir = 'clockwise';
            ax.ThetaZeroLocation = 'bottom';
            ax.GridColor = 'k';
            thetaticks(0:90:270);
            ax.RColorMode = 'manual';
            ax.RColor = 'k';
            ax.ThetaColorMode = 'manual';
            ax.ThetaColor = 'k';
            ax.GridAlpha = 0.3;
            ax.ThetaTickLabel = {'';'';'';''};
            rticks([]);  
            ax.LineWidth = 2; 
            
    nexttile 
    hold on   
        tc1 = tcAll_on(:,ix(ixC(nC)));  
        tc2 = tcAll_off(:,ix(ixC(nC))); 
        tc1 = gaussFiltAng(tc1,3);
        tc2 = gaussFiltAng(tc2,3);
        [~,mIx] = max(tc1);
        tc1 = circshift(tc1,-mIx+180);
        tc2 = circshift(tc2,-mIx+180);
        tc1 = tc1./max(tc2);
        tc2 = tc2./max(tc2);

        xVal = -179:180;
        yVal = [tc1,tc2];
        ax = plot(xVal,yVal(:,1));
        ax.Color = hdc_dark;
        ax.LineWidth = 1;
        ax = plot(xVal,yVal(:,2));
        ax.Color = hdc_light;
        ax.LineWidth = 1;
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.XTick = [-180 0 180];
        xlim([-180 180])
        ylim([0 2])
        ax.TickLength = [0.03 0.025];
        xlabel('Angle from peak (deg)')
        ylabel('Firing rate')
        axis(ax, 'square');  

    nexttile  
    hold on   
        psthBase = mean(allCCG_ind(1:nbins_psth/2,:),1);
        psthProp = allCCG_ind ./psthBase; 
        ccg = psthProp(:,ix(ixC(nC)));
        ax = plot(bins_ind,ccg);
        ax.LineWidth = lw;
        ax.Color = 'k';

        ylabel('Firing rate (fold change)')
        xlabel('Time from light onset (s) ')
        ax = gca;
        axis(ax, 'square');
        ax.FontSize = fs;
        ax.LineWidth = lw;
        ax.TickLength = [0.03 0.025];
        ylim([0 3]);
        xlim([-1 2]);
  
    nexttile 
    hold on
        tc1 = tcAll_on(:,ix(ixC(nC)));  
        tc2 = tcAll_off(:,ix(ixC(nC)));  
        xVal = -179:180;
        cc = TCcrosscorr(tc1,tc2);
        cc = circshift(cc,179);
        ax = plot(xVal,cc);
        ax.Color = 'k';
        ax.LineWidth = 1;
        [~,maxIx] = max(cc);
        maxIx = maxIx - 179;
        ax = plot([maxIx maxIx],[-1 1],':');
        ax.LineWidth = 2;
        ax.Color = 'k';
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.XTick = [-180 -90 0 90 180];
        ax.XTickLabel = [{'-180'},{''},{'0'},{''},{'180'}];
        xlim([-180 180])
        ax.TickLength = [0.03 0.025];
        xlabel('Offset (deg)')
        ylabel('Correlation (r)')
        axis(ax, 'square');  
    
    nexttile 
    hold on    
        tc1 = tcAll60_on(:,ix(ixC(nC)));  
        tc2 = tcAll60_off(:,ix(ixC(nC))); 
        tc1 = tc1./max(tc2);
        tc2 = tc2./max(tc2);  
        ax = scatter(tc2,tc1,'.');
        ax.SizeData = 50;
        ax.MarkerFaceColor = [0.5 0.5 0.5];
        ax.MarkerEdgeColor = [0.5 0.5 0.5];
        
         Fit = polyfit(tc2,tc1,1);
         xVal = [0 1];
         ax = plot(xVal,polyval(Fit,xVal),'-r');
         ax.LineWidth = lw;
         
         ax = plot([0 1],[0 1],'--k');
         ax.LineWidth = lw;       
         
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.YTick = [0 1 2];
        ax.XTick = [0 1];
        xlim([0 1])
        ylim([0 2])
        ax.TickLength = [0.03 0.025];
        xlabel('Tuning OFF')
        ylabel('Tuning ON')
        axis(ax, 'square'); 
end       

