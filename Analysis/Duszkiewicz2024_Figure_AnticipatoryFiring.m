function Duszkiewicz2024_Figure_AnticipatoryFiring

% This script reporoduces panels from:
% Extended Data Figure 3

% Dependencies: 
% TStoolbox

%TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%% Parameters variavles

clear all 

tcAll_lag = [];
isHd = [];
isFs = [];
isGd = [];
isEx = [];
hdi_lag = [];
isAdn = [];
isPos = [];

%% Load data

dataset = List2Cell('dataset_All_FT.list');

for ii = 1:length(dataset)
    
    fprintf(['Uploading recording ', dataset{ii},'\n'])     
    load(fullfile(dataset{ii},'Data','CellTypes'));
    load(fullfile(dataset{ii},'Data','BrainArea'));
    load(fullfile(dataset{ii},'Analysis','HdTuning_AnticipCorr'));

        
    isHd = [isHd; hd];
    isFs = [isFs; fs];
    isGd = [isGd; gd];
    isEx = [isEx; ex]; 
    tcAll_lag = [tcAll_lag hAll]; % tuning curves
    hdi_lag = [hdi_lag hdInfo_lag];
    isPos = [isPos; pos];
    isAdn = [isAdn; adn];
    
end

totC = length(isHd);

%% define groups


% HD cells
ixHd = find(isHd);

% FS cells
ixFs = find(isFs); 

% Excitatory cells
ixEx = find(isEx);


%% pick tuning curves with lags of the highest HD info
[hdi,lagIx] = max(hdi_lag,[],1);

tcAll = nan(360,totC);
for nC = 1:totC
    tcAll(:,nC) = tcAll_lag(:,nC,lagIx(nC));
end

antInt = lagVec(lagIx);

%% Align tuning curves and calculate metrics
nbins = length(b);

% All tuning curves
totLag = size(tcAll_lag,3);
tcAllN_lag = nan(nbins,totC,totLag);
tcSide_lag = nan(nbins/2+1,totC,totLag);
tcWidth_lag = nan(totC,totLag);
tcWidth_ai = nan(totC,1); 

for nL = 1:totLag
    [~,muMax] = max(tcAll_lag(:,:,nL),[],1);     
    for nC = 1:totC    
      if ~isnan(muMax(nC))   
        tc = circshift(tcAll_lag(:,nC,nL),-muMax(nC)+1) ./ max(tcAll_lag(:,nC,nL));
        tcAllN_lag(:,nC,nL) = tc;
        tcSide_lag(:,nC,nL) = (tc(1:nbins/2+1) + flipud([tc(nbins/2+1:end,:);tc(1)]))./2;
        tcw = find(tcSide_lag(:,nC,nL) < 0.5,1);
        if ~isempty(tcw)
             tcWidth_lag(nC,nL) = tcw*2;
        end        
      end
    end
end

for nC = 1:totC
    tcWidth_ai(nC) = tcWidth_lag(nC,lagIx(nC));
end



%% Figures

lw = 1;
fs = 11;

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];
unc = [0.5 0.5 0.5];
adc = [1 0.4 1];

%% Figure 1: Anticipatory interval analysis

figure(1),clf
set(gcf, 'Color', 'w')

subplot(2,3,1)
hold on
    vbl = hdi_lag';
    vbl = vbl./vbl(:,51);
    vbl = vbl./max(vbl,[],2);
    ix1 = find(isHd & isPos);
    ix2 = find(isHd & isAdn);
    m = nanmean(vbl(ix1,:));
    s = sem(vbl(ix1,:));     
    ax = boundedline(lagVec,m,s,'cmap',hdc);
    m = mean(vbl(ix2,:));
    s = sem(vbl(ix2,:)); 
    ax = boundedline(lagVec,m,s,'cmap',adc);
    ax = gca;
    %ax.XTick = [1:maxComp];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    %ylim([0.5 2]);
    %xlim([-0.5 0.5]);
    xlabel('Lag (s)')
    ylabel('Fold change in HD info')   
    
subplot(2,3,2)
hold on
    vbl = hdi_lag';
    vbl = vbl./vbl(:,51);
    vbl = vbl./max(vbl,[],2);
    ix1 = find(isHd & isPos);
    ix2 = find(isHd & isAdn);
    m = nanmean(vbl(ix1,:));
    s = sem(vbl(ix1,:));     
    ax = boundedline(lagVec,m,s,'cmap',hdc);
    m = mean(vbl(ix2,:));
    s = sem(vbl(ix2,:)); 
    ax = boundedline(lagVec,m,s,'cmap',adc);
    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylim([0.95 1]);
    xlim([-0.05 0.1]);
    xlabel('Lag (s)')
    ylabel('Fold change in HD info')
    
subplot(2,3,3)
hold on
    vbl = antInt;
    ix1 = find(isHd & isAdn);
    ix2 = find(isHd & isPos);
    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.BinEdges = lagVec;
    h2.BinEdges = lagVec;
    h1.BinWidth = 0.01;
    h2.BinWidth = 0.01;
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
    xlim([-0.1 0.2]);
    ylim([0 0.3]);
    ax.TickLength = [0.03 0.025];
    %ax.XTick = [0 60 120];
    %ax.YTick = [0 0.25 0.5];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. cells'])
    xlabel(['Anticipatory interval (s)'])
    [p,r,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
    
subplot(2,3,4)
hold on
    vbl = tcWidth_lag;
    vbl = vbl./vbl(:,51);
    ix1 = find(isHd & isPos);
    ix2 = find(isHd & isAdn);
    m = nanmean(vbl(ix1,:));
    s = sem(vbl(ix1,:));     
    ax = boundedline(lagVec,m,s,'cmap',hdc);
    m = mean(vbl(ix2,:));
    s = sem(vbl(ix2,:)); 
    ax = boundedline(lagVec,m,s,'cmap',adc);
    ax = gca;
    %ax.XTick = [1:maxComp];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylim([0.5 3]);
    xlim([-0.5 0.5]);
    xlabel('Lag (s)')
    ylabel('Fold change in TC width')   
    
subplot(2,3,5)
hold on
    vbl = tcWidth_lag;
    vbl = vbl./vbl(:,51);
    ix1 = find(isHd & isPos);
    ix2 = find(isHd & isAdn);
    m = nanmean(vbl(ix1,:));
    s = sem(vbl(ix1,:));     
    ax = boundedline(lagVec,m,s,'cmap',hdc);
    m = mean(vbl(ix2,:));
    s = sem(vbl(ix2,:)); 
    ax = boundedline(lagVec,m,s,'cmap',adc);
    ax = gca;
    %ax.XTick = [1:maxComp];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylim([0.95 1.1]);
    xlim([-0.05 0.1]);
    xlabel('Lag (s)')
    ylabel('Fold change in TC width')
       
    
return   

