
%function Duszkiewicz2024_Figure_SynapticConnections

% This script reporoduces panels from:
% Figure 4, Extended Data Figure 2

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

thSmooth = 3; % smoothing of tuning curves (3 used in Duszkiewicz et al., 2024)

cDep = [];
isHd = [];
isFs = [];
isGd = [];
isEx = [];
cellsSoFar = 0;
ccSynAll = [];
pairsSyn = [];
synLat = [];
synStrZ = [];
totSpk = [];
totRate = [];
tcAll = [];

%% Load data 

dataset = List2Cell('dataset_All.list');

for ii = 1:length(dataset)

    %Load data
    fprintf(['Uploading recording ', dataset{ii},'\n'])  
    load(fullfile(dataset{ii},'Data','CellTypes'));
    load(fullfile(dataset{ii},'Data','CellDepth'));
    load(fullfile(dataset{ii},'Data','SpikeData'));
    load(fullfile(dataset{ii},'Analysis','HdTuning_moveEp'));
    load(fullfile(dataset{ii},'Analysis','SynPairs_P01'));
    load(fullfile(dataset{ii},'Analysis','SynCrossCorrs'));
    
    isHd = [isHd; hd];
    isFs = [isFs; fs];
    isGd = [isGd; gd];
    isEx = [isEx; ex];

    cDep = [cDep; cellDep];
    tcAll = [tcAll hAll(1:360,:,thSmooth+1)];
    ccSynAll = [ccSynAll ccg];
    pairsSyn = [pairsSyn; pairs+cellsSoFar]; 
    cellsSoFar = cellsSoFar + length(hd);
    synLat = [synLat; lat]; 
    synStrZ = [synStrZ; strZ];
    
    % calculate total number of spikes
    
    for nC = 1:length(S)
        ts = Range(S{nC});
        totSpk = [totSpk; length(ts)];
    end
    
    totRate = [totRate; Rate(S)];
end

totC = length(isFs);
cDep = cDep ./ 1000; 

synBins = bins;

%% define groups

isGd(cDep > 0) = 0; % For neuropixel rec

% Good cells 
ixGd = find(isGd == 1 & isGd == 1);

% HD cells
ixHd = find(isHd == 1 & isGd == 1);

% FS cells
ixFs = find(isFs == 1 & isGd == 1);

% Ex cells
ixEx = find(isEx == 1);

%% Calculations
   
% get linear distance between cells    
    diffLinSyn = abs(cDep(pairsSyn(:,1)) - cDep(pairsSyn(:,2)));
    
% find true excitatory synapses 
    latLim = 0.5; % lower bound of synaptic latency
    isSyn = ~isnan(synStrZ(:,2)) & isnan(synStrZ(:,1)) & synLat(:,2) > latLim & ismember(pairsSyn(:,1),ixHd) & ismember(pairsSyn(:,2),ixFs);
    ixSyn = find(isSyn);    
    ixHdFs = find(ismember(pairsSyn(:,1),ixHd) & ismember(pairsSyn(:,2),ixFs));

    fprintf(['Number of excitatory connections: ', num2str(length(ixSyn)),'\n'])
    fprintf(['Average excitatory connection probability: ', num2str(length(ixSyn)./length(ixHdFs)),'\n'])

% values for model
    syn_mean = mean(synStrZ(ixSyn,2));
    syn_sd = sqrt(var(synStrZ(ixSyn,2)));
    sdAsMean = (syn_sd/syn_mean)*100;
    


%% Figures

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];

fs = 11;
lw = 1;

% pick synaptic pairs to display 
ixPairs = [248 1 14];
    
%% Figure 1: Synaptic connectivity HD-FS

figure(1),clf
set(gcf, 'Color', 'w');
tiledlayout(1,4)

nexttile
hold on

    ix1 = ixSyn;
    vbl1 = synLat(ix1,2);
    m1 = median(vbl1);
    h1 = histogram(vbl1);
    h1.BinWidth = 0.2;
    h1.BinEdges = h1.BinEdges -0.1;
    h1.DisplayStyle = 'stairs';
    h1.LineWidth = 2;
    h1.EdgeColor = [0.5 0.5 0.5];
    h1.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    ax.XTick = [0:1:2];
    xlim([0 2]);
    ylabel(['No. FS-HD connections'])
    xlabel(['Latency (ms)'])
    axis(ax,'square');
       
nexttile
hold on

    edges = 0:0.025:0.8;
    bins = edges(1:end-1)+diff(edges)./2;
    
    ix = ixSyn;
    vbl1 = histcounts(diffLinSyn(ix),edges);
    
    ix = ixHdFs;
    vbl2 = histcounts(diffLinSyn(ix),edges);
    
    vbl3 = vbl1./vbl2;
    
    ax = plot(bins,vbl3);   
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    
    ax = gca; 
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    ax.XTick = [0 0.4 0.8];
    xlim([0 0.8]);
    ylabel(['Prop. HD-FS connections'])
    xlabel(['Linear dist. (mm)'])
    axis(ax,'square');
    
 nexttile
 hold on
    ix = ixSyn;
    vbl1 = synStrZ(ix,2);
    xVal = 1:40;
    m1 = mean(vbl1);
    h1 = histogram(vbl1);
    h1.BinWidth = 1;
    %h1.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h1.LineWidth = 2;
    h1.EdgeColor = [0.5 0.5 0.5];
    h1.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    
    ax = gca; 
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    xlim([0 70])
    ylabel(['No. HD-FS connections'])
    xlabel(['Synaptic strength (z-score)'])
    axis(ax,'square');
    
    xVal = synStrZ(ixSyn(ixPairs),2);
    ax = scatter(xVal,[0 0 0]);
    ax.MarkerFaceColor = 'r';
    ax.MarkerEdgeColor = 'none';
    

 nexttile
 hold on
    ix = ixSyn;
    vbl1 = log10(synStrZ(ix,2));
    m1 = mean(vbl1);
    h1 = histogram(vbl1);
    h1.BinWidth = 0.1;
    %h1.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h1.LineWidth = 2;
    h1.EdgeColor = [0.5 0.5 0.5];
    h1.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    ylabel(['No. HD-FS connections'])
    xlabel(['log(synaptic strength)'])
    axis(ax,'square');
        
    xVal = log10(synStrZ(ixSyn(ixPairs),2));
    ax = scatter(xVal,[0 0 0]);
    ax.MarkerFaceColor = 'r';
    ax.MarkerEdgeColor = 'none';
    

%%% Sanity checks %%%

    syn_mean = mean(synStrZ(ixSyn,2));
    syn_sd = sqrt(var(synStrZ(ixSyn,2)));
    
%% Figure 2: Synaptic connectivity HD-FS - examples

figure(2),clf
set(gcf, 'Color', 'w');
tiledlayout(3,5)


pairsSyn(ixSyn(ixPairs),:);

binsize = round(median(diff(synBins)),5) ./ 1000;

bix = synBins >= 0 & synBins <= 2;
alpha = 0.01;
b = deg2rad(1:360);

for nC = 1:length(ixPairs)
       
    nexttile
        pair = pairsSyn(ixSyn(ixPairs(nC)),:);
        tc = tcAll(:,pair(1));  
        ax = polarplot(b,tc);
        hold on
        ax.Color = hdc;
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
        pair = pairsSyn(ixSyn(ixPairs(nC)),:);
        tc = tcAll(:,pair(2));  
        ax = polarplot(b,tc);
        hold on
        ax.Color = fsc;
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
        
        pair = pairsSyn(ixSyn(ixPairs(nC)),:);
        ccg = ccSynAll(:,ixSyn(ixPairs(nC)));
        [~, pred,~] = cch_conv(round(ccg),21);
        
         % find the bound again (should include in the output next time)       
        hiBound = poissinv( 1-alpha, max(pred(bix,:), [], 1 ));
        hiBins = bsxfun( @gt, ccg(bix), hiBound );
        predMin = pred-sqrt(pred)*2;
        predMax = pred+sqrt(pred)*2;
        
        % revert to rate and normalize by postsynaptic
        ccg = ccg /binsize / totSpk(pair(1)) / totRate(pair(2));
        pred = pred /binsize / totSpk(pair(1)) / totRate(pair(2));
        hiBound = hiBound /binsize / totSpk(pair(1)) / totRate(pair(2));
        predMin = predMin /binsize / totSpk(pair(1)) / totRate(pair(2));        
        predMax = predMax /binsize / totSpk(pair(1)) / totRate(pair(2));

        ax = plot([0 0],[0 1.2*max(ccg)],':r');
        ax.LineWidth = 2;

        ax = bar(synBins,ccg);
        ax.BarWidth = 1;
        ax.FaceColor = 'k';
        ax.EdgeColor = 'none';
        ax.LineWidth = 0.5;

        ax = plot(synBins,pred);
        ax.Color = [0.6 0.6 0.6];
        ax.LineWidth = lw;

        ax = plot(synBins,predMin,'--');
        ax.Color = [0.6 0.6 0.6];
        ax.LineWidth = lw;

        ax = plot(synBins,predMax,'--');
        ax.Color = [0.6 0.6 0.6];
        ax.LineWidth = lw;   

        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.TickLength = [0.03 0.025];
        ylabel(['Firing rate (fold increase)'])
        xlabel(['Lag (ms)']) 
        box off
        ylim([0 1.2 * max(ccg)])
        xlim([-10 10])
        axis(ax,'square');
        title(round(log(synStrZ(ixSyn(ixPairs(nC)),2)),1))
    
    nexttile
    hold on
        
        vbl = ccg - pred;

        ax = plot([0 2],[hiBound - max(pred(bix,:)) hiBound - max(pred(bix,:))],':g');
        ax.LineWidth = 2;  
                
        ax = bar(synBins,vbl);
        ax.BaseValue = 0;
        ax.BarWidth = 1;
        ax.FaceColor = 'k';
        ax.EdgeColor = 'none';
        
        
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.TickLength = [0.03 0.025];
        ylabel(['Spike density (Hz)'])
        xlabel(['Lag (ms)']) 
        box off
        ylim([1.2* min(vbl) 1.2 * max(vbl)])
        xlim([-2 4])
        axis(ax,'square');
        
    nexttile
    hold on
        
        ax = plot([0 2],[hiBound hiBound],':g');
        ax.LineWidth = 2;  
                      
        vbl2 = vbl;
        vbl2(vbl2 < hiBound- max(pred(bix,:))) = hiBound- max(pred(bix,:));
        ax = bar(synBins,vbl2,'BaseValue',hiBound- max(pred(bix,:)));
        ax.BaseValue = hiBound- max(pred(bix,:));
        ax.BarWidth = 1;
        ax.FaceColor = 'r';
        ax.EdgeColor = 'none';
        
        ax = gca;
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.TickLength = [0.03 0.025];
        ylabel(['Spike density (Hz)'])
        xlabel(['Lag (ms)']) 
        box off
        xlim([-2 4])
        axis(ax,'square');
        
end
