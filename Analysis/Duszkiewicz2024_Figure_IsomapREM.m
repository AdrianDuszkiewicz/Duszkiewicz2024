function Duszkiewicz2024_Figure_IsomapREM

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

rateRem = [];
tcReal = [];
tcWake = [];
tcRem = [];

tcRealSh = [];
tcWakeSh = [];
tcRemSh = [];

whichRec = [];
isHd = [];
isFs = [];
isGd = [];
totCells = [];

%% Load data
dataset = List2Cell('dataset_All_Rem.list');
noRec = length(dataset);
radius_wk = cell(noRec,1);
radius_rem = cell(noRec,1);

for ii = 1:length(dataset)

    %Load data
    fprintf(['Uploading recording ', dataset{ii},'\n'])
    load(fullfile(dataset{ii},'Data','CellTypes'));
    load(fullfile(dataset{ii},'Analysis','MeanFR'));
    
    isHd = [isHd; hd];
    isFs = [isFs; fs];
    isGd = [isGd; gd];
    whichRec = [whichRec; repmat(ii,length(hd),1)];
    rateRem = [rateRem; rateREM]; 
    
    load(fullfile(dataset{ii},'Analysis','Isomap_REMvsWK'));
    tcReal = [tcReal hReal];
    tcWake = [tcWake hWk];
    tcRem = [tcRem hRem];
    iso_wk = mapping(ixWk,:);
    iso_rem = mapping(ixRem,:);
    totCells = [totCells; length(ixHd)];
    
    load(fullfile(dataset{ii},'Analysis','Isomap_REMvsWK_Shuffle'));
    tcRealSh = [tcRealSh hReal];
    tcWakeSh = [tcWakeSh hWk];
    tcRemSh = [tcRemSh hRem];
        
    % calculate the Isomap radius in Wake and REM 
    totS = size(iso_wk,1);
    r_wk = nan(totS,1);
    r_rem = nan(totS,1);
    
    for nS = 1:totS   
        pts = [0,0; iso_wk(nS,1),iso_wk(nS,2)];
        d = pdist(pts,'euclidean');
        r_wk(nS) = d;

        pts = [0,0; iso_rem(nS,1),iso_rem(nS,2)];
        d = pdist(pts,'euclidean');
        r_rem(nS) = d;
    end    
    
   radius_wk{ii} = r_wk;
   radius_rem{ii} = r_rem;
       
end

%% Define Groups
totC = length(isGd);
% Good cells 
ixGd = find(isGd == 1 & rateRem > 0.5);

% HD cells
ixHd = find(isHd == 1 & rateRem > 0.5);

% FS cells
ixFs = find(isFs == 1 & rateRem > 10);

%% Calculations 

% Calculate offset of each Isomap and correct tuning curves
recs = unique(whichRec);
totR = length(recs);

for nR = 1:totR    
    ixR = find(whichRec == nR); 
    ix = ixHd(find(ismember(ixHd,ixR)));
    tC = length(ix);
    ofW = nan(tC,1);
    ofWsh = nan(tC,1);
    for nC = 1:tC        
        cc = TCcrosscorr(tcReal(:,ix(nC)),tcWake(:,ix(nC))); % xcorr between two tuning curves
        [~,mx] = max(cc); % max point on the xcorr is the offset
        ofW(nC) = mx;
        cc = TCcrosscorr(tcRealSh(:,ix(nC)),tcWakeSh(:,ix(nC))); % xcorr between two tuning curves
        [~,mx] = max(cc); % max point on the xcorr is the offset
        ofWsh(nC) = mx;
    end
    meanOf = round(mean(ofW),0);
    tcWake(:,ixR) = circshift(tcWake(:,ixR),meanOf);
    tcRem(:,ixR) = circshift(tcRem(:,ixR),meanOf);
    meanOf = round(mean(ofWsh),0);
    tcWakeSh(:,ixR) = circshift(tcWakeSh(:,ixR),meanOf);
    tcRemSh(:,ixR) = circshift(tcRemSh(:,ixR),meanOf);
end

% Calculate all correlations between real and iso turning curves
nBins = size(tcReal,1);
cc_RlWk = nan(nBins,totC);
cc_RlRem = nan(nBins,totC);
cc_WkRem = nan(nBins,totC);

for nC = 1:totC 
    cc_RlWk(:,nC) = TCcrosscorr(tcReal(:,nC),tcWake(:,nC));
    cc_RlRem(:,nC) = TCcrosscorr(tcReal(:,nC),tcRem(:,nC));
    cc_WkRem(:,nC) = TCcrosscorr(tcReal(:,nC),tcRem(:,nC));   
end

cc_RlWkSh = nan(nBins,totC);
cc_RlRemSh = nan(nBins,totC);
cc_WkRemSh = nan(nBins,totC);

for nC = 1:totC 
    cc_RlWkSh(:,nC) = TCcrosscorr(tcRealSh(:,nC),tcWakeSh(:,nC));
    cc_RlRemSh(:,nC) = TCcrosscorr(tcRealSh(:,nC),tcRemSh(:,nC));
    cc_WkRemSh(:,nC) = TCcrosscorr(tcWakeSh(:,nC),tcRemSh(:,nC));   
end

% Make sure Isomap tuning curves now have zero offset
offWk = [];
offRem = [];

for nR = 1:totR    
    ix = find(whichRec == nR);  
    tC = length(ix);
    ofW = nan(tC,1);
    ofR = nan(tC,1);
    for nC = 1:tC        
        cc = cc_RlWk(:,ix(nC));
        [~,mx] = max(cc); % max point on the xcorr is the offset
        ofW(nC) = mx;
        cc = cc_RlRem(:,ix(nC));
        [~,mx] = max(cc); % max point on the xcorr is the offset
        ofR(nC) = mx;
    end
    offWk = [offWk; ofW];
    offRem = [offRem; ofR];
end
offWk = offWk*360./nBins;
offRem = offRem*360./nBins;

offWkSh = [];
offRemSh = [];

for nR = 1:totR    
    ix = find(whichRec == nR);  
    tC = length(ix);
    ofW = nan(tC,1);
    ofR = nan(tC,1);
    for nC = 1:tC        
        cc = cc_RlWkSh(:,ix(nC));
        [~,mx] = max(cc); % max point on the xcorr is the offset
        ofW(nC) = mx;
        cc = cc_RlRemSh(:,ix(nC));
        [~,mx] = max(cc); % max point on the xcorr is the offset
        ofR(nC) = mx;
    end
    offWkSh = [offWkSh; ofW];
    offRemSh = [offRemSh; ofR];
end
offWkSh = offWkSh*360./nBins;
offRemSh = offRemSh*360./nBins;


% get zero correlations

corrWk = cc_RlWk(1,:)';
corrRem = cc_WkRem(1,:)';

corrWkSh = cc_RlWkSh(1,:)';
corrRemSh = cc_WkRemSh(1,:)';


%%% Figures %%%

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];
unc = [0.5 0.5 0.5];
adc = [1 0.4 1];

wkc = [0.15 0.3 0.5];
rmc = [0.3 0.5 0.2];

fsc_light = [1 0.4 0.4];
fsc_dark = [0.6 0.1 0.1];

fs = 11;
lw = 1;

  %% Figure 1: Zero-lag correlation between HD and FS cells

figure (1), clf
set(gcf,'Color','w')
tiledlayout(1,2)
         
nexttile
hold on
    vbl1 = corrRem(ixHd);
    vbl2 = corrRemSh(ixHd);
       
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
    [p,h,stats] = signrank(vbl1,vbl2)
    ax = title(p);
    ax.FontWeight = 'normal';
    
nexttile
hold on
    vbl1 = corrRem(ixFs);
    vbl2 = corrRemSh(ixFs);
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
    [p,h,stats] = signrank(vbl1,vbl2)
    ax = title(p);
    ax.FontWeight = 'normal';
    
%% Figure 2: Isomap radius
figure (2), clf
set(gcf,'Color','w')
tiledlayout(1,2)

totRecs = length(radius_wk); 
meanR_wk = nan(totR,1);
meanR_rem = nan(totR,1);

nexttile
hold on
    
    step = 0.05;
    edges = 0:step:3;
    bins = edges(1:end-1)+step./2;
    smTh = 2;
   
    for nR = 1:totRecs
        r_wk = radius_wk{nR}./totCells(nR); 
        meanR_wk(nR) = mean(r_wk);
        r_wk = r_wk./meanR_wk(nR); 
        n = histcounts(r_wk,edges);
        n = n./sum(n);
        n = gaussFilt(n,smTh);
        ax = plot(bins,n);
        ax.LineWidth = 0.5;
        ax.Color = wkc;     
        
        r_rem = radius_rem{nR}./totCells(nR); 
        meanR_rem(nR) = mean(r_rem);
        r_rem = r_rem./meanR_rem(nR); 
        n = histcounts(r_rem,edges);
        n = n./sum(n);
        n = gaussFilt(n,smTh);
        ax = plot(bins,n);
        ax.LineWidth = 0.5;
        ax.Color = rmc; 
    end
       
    ax = gca; 
    xlim([0 2.5]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. vectors'])
    xlabel(['Radius (norm.)'])
    
nexttile
hold on
    
    vbl1 = meanR_wk;
    vbl2 = meanR_rem;
    
    ax = plot([1 2],[vbl1,vbl2],'-','Color',[0.7 0.7 0.7]);
    ax = plot([0.8 1.2],[mean(vbl1),mean(vbl1)],'-','Color',[0 0 0]);
    ax.LineWidth = lw;
    ax = plot([1.8 2.2],[mean(vbl2),mean(vbl2)],'-','Color',[0 0 0]);
    ax.LineWidth = lw;
    
    ax = gca; 
    ax.XTick = [1 2];
    ax.XTickLabel = {'WAKE';'REM'};
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
