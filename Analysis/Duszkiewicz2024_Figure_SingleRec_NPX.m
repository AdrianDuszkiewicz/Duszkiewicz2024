
function Duszkiewicz2024_Figure_SingleRec_NPX

% This script reporoduces panels from:
% Figure 1
% Run in folder A8603-210602 to reproduce the manuscript figure

% Dependencies: 
% TStoolbox

%TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


%% Parameters and data
clear all 

smoothTC = 3; % smoothing of tuning curves (3 used in Duszkiewicz et al., 2024)

    %Load data

    load(fullfile('Data','WaveformFeatures'));
    load(fullfile('Data','SpikeData'));
    load(fullfile('Data','CellTypes'));
    load(fullfile('Analysis','MeanFR'));
    load(fullfile('Analysis','HdTuning_moveEp'));
    
    depth = depth - max(depth);
    cellDep = depth ./ 1000; 

    hdi = hdInfo(:,smoothTC+1);
    totC = length(fs);

    hAll = hAll(:,:,smoothTC+1);

%% define groups

% HD cells
ixHd = find(hd == 1);

% FS cells
ixFs = find(fs == 1);

% Ex cells
ixEx = find(ex == 1);


%% calculations

%Normalize tuning curves
[~,prefAng] = max(hAll,[],1);
hNormAll = nan(360,totC);
for nC = 1:totC
    m = max(hAll(:,nC));
    hNormAll(:,nC) = circshift(hAll(:,nC),180-prefAng(nC)) ./ m(1);
end


%% Figures

hdc = [128 0 255]./255;


%% Figure 1: Depth profile     

figure (1), clf
set(gcf,'color','w');

ix = ixEx;
ixC = [70 91 112 12];

vbl = hdi;
[~,sIx] = sort(cellDep(ixEx(ixC)));
ixC = ixC(flipud(sIx));

subplot(2,3,1)
hold on

    isCells = ex;
    ix1 = find(isCells); 
    % defne brackets
    binsize = 0.05;
    binMin = -2.5;
    binMax = 0;
    step = 0.05;
        
    v1 = [binMax:-step:binMin+binsize]';
    v2 = [binMax-binsize:-step:binMin]';
    brackets = [v1 v2];
    nBins = size(brackets,1);   

    % compute mean hdi for each bracket
    hdi_mean = nan(nBins,1);
    for nB = 1:nBins
         is = cellDep < brackets(nB,1) & cellDep >= brackets (nB,2) & isCells;
         if sum(is) > 0
            hdi_mean(nB,1) = mean(vbl(is));
         end
    end         
        
    % interpolate NaNs and smoothen (excluding NaNs)
    t = 1:numel(hdi_mean);
    hdi_mean(isnan(hdi_mean)) = interp1(t(~isnan(hdi_mean)), hdi_mean(~isnan(hdi_mean)), t(isnan(hdi_mean)));
    hdi_mean(~isnan(hdi_mean)) = gaussFilt(hdi_mean(~isnan(hdi_mean)),2);  

    bins = mean(brackets,2);
    ax = scatter(vbl(ix1),cellDep(ix1));
       ax.MarkerFaceColor = [0.6 0.6 0.6];
       ax.MarkerEdgeColor = 'none';
       ax.SizeData = 20;

    ax = scatter(vbl(ix(ixC(1))),cellDep(ix(ixC(1))));
       ax.MarkerFaceColor = 'k';
       ax.MarkerEdgeColor = 'none';
       ax.SizeData = 30;           
    ax = scatter(vbl(ix(ixC(2))),cellDep(ix(ixC(2))));
       ax.MarkerFaceColor = 'k';
       ax.MarkerEdgeColor = 'none';
       ax.SizeData = 30;  
    ax = scatter(vbl(ix(ixC(3))),cellDep(ix(ixC(3))));
       ax.MarkerFaceColor = 'k';
       ax.MarkerEdgeColor = 'none';
       ax.SizeData = 30;   
    ax = scatter(vbl(ix(ixC(4))),cellDep(ix(ixC(4))));
       ax.MarkerFaceColor = 'k';
       ax.MarkerEdgeColor = 'none';
       ax.SizeData = 30;   
           
       ax = plot(hdi_mean, bins);
       ax.LineWidth = 1.5;
       ax.Color = hdc;
           
       ax = gca;
        ax.LineWidth = 1;
        ax.FontSize = 11;
        ax.YTick = [-2 -1.5 -1 -0.5 0];
        ax.TickLength = [0.03 0.025];
        xlim([-0.1 2.5])
        ylim([-2.4 0])
        xlabel(['HD info (bits/spike)']);
        ylabel(['Dist. from dura (mm)']);
        axis(ax, 'square');
        
subplot(2,3,2)    
    ix = ixEx;
    tc = hAll(:,ix(ixC(1)));  
    tc = tc./max(tc);
    ax = polarplot(b,tc);
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
        
 subplot(2,3,3)   
    ix = ixEx;
    tc = hAll(:,ix(ixC(2)));  
    tc = tc./max(tc);
    ax = polarplot(b,tc);
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
    
subplot(2,3,4)
 hold on 
 
     [~,sIx] = sort(depth(ix1));
     tc = hNormAll(:,ix1(sIx));

     imagesc(tc');
         ax = gca;
         xlim([0 360])
         ylim([0 size(tc,2)])
         ax.Visible = 'off';
         colormap(summer)
         axis(ax, 'square');
     ax = axes('Position',ax.Position);
         ax.Color = 'none';
         ax.YDir = 'reverse';
         xlim([-180 180])
         ylim([0 size(tc,2)])
         ax.LineWidth = 1;
        ax.FontSize = 11;
        ax.XTick = [-180 0 180];
        ax.YTick = [1 size(tc,2)];
        xlabel(['Angle from peak (deg)']);
        ylabel(['Cell number']);
        axis(ax, 'square');
        box on
        
  subplot(2,3,5)   
    ix = ixEx;
    tc = hAll(:,ix(ixC(3))); 
    tc = tc./max(tc);
    ax = polarplot(b,tc);
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
        
  subplot(2,3,6)   
    ix = ixEx;
    tc = hAll(:,ix(ixC(4))); 
    tc = tc./max(tc);
    ax = polarplot(b,tc);
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
        
        
     
