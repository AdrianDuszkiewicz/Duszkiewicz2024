function Duszkiewicz2024_Figure_IsomapWake

% This script reporoduces panels from:
% Figure 7, Extended Data Figure 10

% Dependencies: 

%TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%% Parameters and variables
clear all 

tcReal = [];
tcWake = [];
tcRv = [];

tcRealSh = [];
tcWakeSh = [];
tcRvSh = [];

whichRec = [];
isHd = [];
isFs = [];
isGd = [];
isEx = [];
totCells = [];

%% Load data
dataset = List2Cell('dataset_All.list');
noRec = length(dataset);
radius = cell(noRec,1);
radius_sh = cell(noRec,1);

for ii = 1:length(dataset)
    
    fprintf(['Uploading recording ', dataset{ii},'\n'])
    load(fullfile(dataset{ii},'Data','WaveformFeatures'));
    load(fullfile(dataset{ii},'Data','CellTypes'));
    load(fullfile(dataset{ii},'Data','SpikeData'));
    
    isHd = [isHd; hd];
    isFs = [isFs; fs];
    isGd = [isGd; gd];
    isEx = [isEx; ex];
    whichRec = [whichRec; repmat(ii,length(hd),1)];
       
    load(fullfile(dataset{ii},'Analysis','Isomap_WK_HdOnly'));
    tcReal = [tcReal hReal];
    tcWake = [tcWake hWk];
    tcRv = [tcRv hRev];
    iso_real = mapping;
    totCells = [totCells; length(ixEx)];
    
    load(fullfile(dataset{ii},'Analysis','Isomap_WK_HdOnly_Shuffle'));
    tcRealSh = [tcRealSh hReal];
    tcWakeSh = [tcWakeSh hWk];
    tcRvSh = [tcRvSh hRev];
    iso_sh = mapping;
    
    % calculate the Isomap radius in Wake 
    totS = size(iso_real,1);
    r = nan(totS,1);
    r_sh = nan(totS,1);
    
    for nS = 1:totS   
        pts = [0,0; iso_real(nS,1),iso_real(nS,2)];
        d = pdist(pts,'euclidean');
        r(nS) = d;

        pts = [0,0; iso_sh(nS,1),iso_sh(nS,2)];
        d = pdist(pts,'euclidean');
        r_sh(nS) = d;
    end    
    
   radius{ii} = r;
   radius_sh{ii} = r_sh;
end

%% Define Groups
totC = length(isGd);
% Good cells 
ixGd = find(isGd == 1);

% Ex cells
ixEx = find(isEx == 1);

% HD cells
ixHd = find(isHd == 1);

% FS cells
ixFs = find(isFs == 1);

%% Calculations 

% Calculate offset of each Isomap and correct tuning curves
recs = unique(whichRec);
totR = length(recs);

for nR = 1:totR    
    ixR = find(whichRec == nR); 
    ix = ixHd(find(ismember(ixHd,ixR)));
    tC = length(ix);
    ofW = nan(tC,1);
    ofR = nan(tC,1);
    for nC = 1:tC        
        cc = TCcrosscorr(tcReal(:,ix(nC)),tcWake(:,ix(nC))); % xcorr between two tuning curves
        [~,mx] = max(cc); % max point on the xcorr is the offset
        ofW(nC) = mx;
        cc = TCcrosscorr(tcReal(:,ix(nC)),tcRv(:,ix(nC))); % xcorr between two tuning curves
        [~,mx] = max(cc); % max point on the xcorr is the offset
        ofR(nC) = mx;
    end
    meanOf = round(rad2deg(CircularMean(deg2rad(ofW))),0);
    tcWake(:,ixR) = circshift(tcWake(:,ixR),meanOf);
    meanOf = round(rad2deg(CircularMean(deg2rad(ofR))),0);
    tcRv(:,ixR) = circshift(tcRv(:,ixR),meanOf);
end

for nR = 1:totR    
    ixR = find(whichRec == nR); 
    ix = ixHd(find(ismember(ixHd,ixR)));
    tC = length(ix);
    ofW = nan(tC,1);
    ofR = nan(tC,1);
    for nC = 1:tC        
        cc = TCcrosscorr(tcRealSh(:,ix(nC)),tcWakeSh(:,ix(nC))); % xcorr between two tuning curves
        [~,mx] = max(cc); % max point on the xcorr is the offset
        ofW(nC) = mx;
        cc = TCcrosscorr(tcRealSh(:,ix(nC)),tcRvSh(:,ix(nC))); % xcorr between two tuning curves
        [~,mx] = max(cc); % max point on the xcorr is the offset
        ofR(nC) = mx;
    end
    meanOf = round(rad2deg(CircularMean(deg2rad(ofW))),0);
    tcWakeSh(:,ixR) = circshift(tcWakeSh(:,ixR),meanOf);
    meanOf = round(rad2deg(CircularMean(deg2rad(ofR))),0);
    tcRvSh(:,ixR) = circshift(tcRvSh(:,ixR),meanOf);
end

% Calculate all correlations between real and iso turning curves
nBins = size(tcReal,1);
cc_RlWk = nan(nBins,totC);
cc_RlRv = nan(nBins,totC);

for nC = 1:totC 
    cc_RlWk(:,nC) = TCcrosscorr(tcReal(:,nC),tcWake(:,nC)); 
    cc_RlRv(:,nC) = TCcrosscorr(tcReal(:,nC),tcRv(:,nC)); 
end

cc_RlWkSh = nan(nBins,totC);
cc_RlRvSh = nan(nBins,totC);

for nC = 1:totC 
    cc_RlWkSh(:,nC) = TCcrosscorr(tcRealSh(:,nC),tcWakeSh(:,nC)); 
    cc_RlRvSh(:,nC) = TCcrosscorr(tcRealSh(:,nC),tcRvSh(:,nC)); 
end

% Make sure Isomap tuning curves now have zero offset
offWk = nan(totC,1);
offWkSh = nan(totC,1);

for nC = 1:totC          
        cc = cc_RlWk(:,nC);
        [~,mx] = max(cc); % max point on the xcorr is the offset
        offWk(nC) = mx;  
        cc = cc_RlWkSh(:,nC);
        [~,mx] = max(cc); % max point on the xcorr is the offset
        offWkSh(nC) = mx;
end
offWk = offWk*360./nBins;
offWkSh = offWkSh*360./nBins;

% calculate offset score
offScWk = offWk;
offScWk(find(offScWk >=180)) = offScWk(find(offScWk >=180)) - 360;
offScWk = abs(offScWk);

% get zero offset correlations
corrWk = cc_RlWk(1,:)';
corrRv = cc_RlRv(1,:)';
corrWkSh = cc_RlWkSh(1,:)';
corrRvSh = cc_RlRvSh(1,:)';


%%% Figures %%%

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];
wkc = [0.15 0.3 0.5];

fs = 11;
lw = 1;
    
 
 %% Figure 1: Zero-lag correlation between HD and FS cells

figure (1), clf
set(gcf,'Color','w')
tiledlayout(1,2)
         
nexttile
hold on
    vbl1 = corrWk(ixHd);
    vbl2 = corrWkSh(ixHd);
       
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.05;
    h2.BinWidth = 0.05;
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
    xlim([-1 1]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. HD cells'])
    xlabel(['Correlation (r)'])
    [p,h,stats] = signrank(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
nexttile
hold on
    vbl1 = corrWk(ixFs);
    vbl2 = corrWkSh(ixFs);
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.05;
    h2.BinWidth = 0.05;
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
    [p,h,stats] = signrank(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
 

%% Figure 2: Isomap radius
figure (2), clf
set(gcf,'Color','w')
tiledlayout(1,2)

totRecs = length(radius); 
meanR = nan(totR,1);
meanR_sh = nan(totR,1);

nexttile
hold on
    
    step = 0.02;
    edges = 0:step:4;
    bins = edges(1:end-1)+step./2;
    smTh = 3;
   
    for nR = 1:totRecs
        r = radius{nR} ./totCells(nR); 
        meanR(nR) = mean(r);
        r = r./meanR(nR); 
        n = histcounts(r,edges);
        n = n./sum(n);
        n = gaussFilt(n,smTh);
        ax = plot(bins,n);
        ax.LineWidth = 0.5;
        ax.Color = wkc;     
        
        r_sh = radius_sh{nR} ./totCells(nR); 
        meanR_sh(nR) = mean(r_sh);
        r_sh = r_sh./meanR_sh(nR); 
        n = histcounts(r_sh,edges);
        n = n./sum(n);
        n = gaussFilt(n,smTh);
        ax = plot(bins,n);
        ax.LineWidth = 0.5;
        ax.Color = [0.7 0.7 0.7]; 
    end
       
    ax = gca; 
    xlim([0 4]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. vectors'])
    xlabel(['Radius'])
    
nexttile
hold on
    
    vbl1 = meanR;
    vbl2 = meanR_sh;
    
    ax = plot([1 2],[vbl1,vbl2],'-','Color',[0.7 0.7 0.7]);
    ax = plot([0.8 1.2],[mean(vbl1),mean(vbl1)],'-','Color',[0 0 0]);
    ax.LineWidth = lw;
    ax = plot([1.8 2.2],[mean(vbl2),mean(vbl2)],'-','Color',[0 0 0]);
    ax.LineWidth = lw;
    
    ax = gca; 
    ax.XTick = [1 2];
    ax.XTickLabel = {'WAKE';'Shuffle'};
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    xlim([0.5 2.5]);
    ylim([0 0.3]);
    ylabel(['Mean radius (a.u.)'])   
    [p,h,stats] = signrank(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
      
return
    
