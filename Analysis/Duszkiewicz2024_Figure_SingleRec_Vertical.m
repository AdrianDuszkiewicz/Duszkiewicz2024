% This script reporoduces panels from:
% Extended Data Figure 1
% Run in folder A6216-210623 to reproduce the manuscript figure 

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

thSmooth = 3; % smoothing of tuning curves (3 used in Duszkiewicz et al., 2024)

%Load data
load(fullfile('Data','Waveforms'));
load(fullfile('Data','CellTypes'));
load(fullfile('Data','CellDepth'));
load(fullfile('Analysis','HdTuning_moveEp'));

hdi = hdInfo(:,thSmooth+1);
totC = length(fs);
cellDep = cellDep ./ 1000; 
hAll = hAll(:,:,thSmooth+1);

%% define groups
% Good cells 
ixGd = find(gd);

% HD cells
ixHd = find(hd);

% FS cells
ixFs = find(fs);

% Ex cells
ixEx = find(ex);


%% Calculations

[~,prefAng] = max(hAll,[],1);


%% Figures

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];

fs = 11;
lw = 1;

% pick cells to display 

ixC_hd = [5 28 9]; 
[~,sIx] = sort(cellDep(ixEx(ixC_hd)));
ixC_hd = ixC_hd(flipud(sIx));

ixC_fs = [12 6 2];
[~,sIx] = sort(cellDep(ixFs(ixC_fs)));
ixC_fs = ixC_fs(flipud(sIx));

%% Figure 1: Depth profile and examples    

figure (1), clf
set(gcf,'color','w');

subplot(2,4,1)
hold on

    isCells = hd;
    vbl = hdi;           
    % defne brackets
    binsize = 0.05;
    binMin = -0.8;
    binMax = 0;
    step = 0.025;
        
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
       
    % interpolate NaNs and smoothen (excluding end NaNs)
    t = 1:numel(hdi_mean);
    hdi_mean(isnan(hdi_mean)) = interp1(t(~isnan(hdi_mean)), hdi_mean(~isnan(hdi_mean)), t(isnan(hdi_mean)));
    hdi_mean(~isnan(hdi_mean)) = gaussFilt(hdi_mean(~isnan(hdi_mean)),2);  
    
       % plot
       bins = mean(brackets,2);
       
       ax = plot([0.2 0.2],[-0.85 0.05],':k');
       ax.LineWidth = lw;
       
       ix = find(ex & ~hd);  
       ax = scatter(vbl(ix),cellDep(ix));
           ax.MarkerFaceColor = [0.7 0.7 0.7];
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 20;
           
       ix = ixHd;  
       ax = scatter(vbl(ix),cellDep(ix));
           ax.MarkerFaceColor = hdc;
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 20;
           
       ix = ixEx; 
       ax = scatter(vbl(ix(ixC_hd(1))),cellDep(ix(ixC_hd(1))));
           ax.MarkerFaceColor = 'k';
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 30;           
       ax = scatter(vbl(ix(ixC_hd(2))),cellDep(ix(ixC_hd(2))));
           ax.MarkerFaceColor = 'k';
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 30;  
       ax = scatter(vbl(ix(ixC_hd(3))),cellDep(ix(ixC_hd(3))));
           ax.MarkerFaceColor = 'k';
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 30;                      
                   
       ax = plot(hdi_mean, bins);
       ax.LineWidth = 1.5;
       ax.Color = hdc;
           
       ax = gca;
       xlim([-0.2 3])
       ylim([-0.85 0.05])
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.TickLength = [0.03 0.025];

        xlabel(['HD info (bits/spike)']);
        ylabel(['Dist. from top (mm)']);
        axis(ax, 'square');
        

subplot(2,4,2)    
    ix = ixEx;
    tc = hAll(:,ix(ixC_hd(1)));  
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
        
 subplot(2,4,3)   
    ix = ixEx;
    tc = hAll(:,ix(ixC_hd(2)));  
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
    
  subplot(2,4,4)   
    ix = ixEx;
    tc = hAll(:,ix(ixC_hd(3))); 
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

% 2nd row 
subplot(2,4,5)
hold on

    isCells = fs;
    vbl = hdi;
    
    % defne brackets
    binsize = 0.05;
    binMin = -0.8;
    binMax = 0;
    step = 0.025;
        
    v1 = [binMax:-step:binMin+binsize]';
    v2 = [binMax-binsize:-step:binMin]';
    brackets = [v1 v2];
    nBins = size(brackets,1);   
      
    % compute mean hdi for each bracket
    hdi_mean = nan(nBins,1);
    hdi_sem = nan(nBins,1);
    for nB = 1:nBins
         is = cellDep < brackets(nB,1) & cellDep >= brackets (nB,2) & isCells;
         if sum(is) > 0
            hdi_mean(nB,1) = mean(vbl(is));
            hdi_sem(nB,1) = std(vbl(is))./sqrt(numel(vbl(is)));
         end
    end        
        
    % interpolate NaNs and smoothen (excluding end NaNs)
    t = 1:numel(hdi_mean);
    hdi_mean(isnan(hdi_mean)) = interp1(t(~isnan(hdi_mean)), hdi_mean(~isnan(hdi_mean)), t(isnan(hdi_mean)));
    hdi_mean(~isnan(hdi_mean)) = gaussFilt(hdi_mean(~isnan(hdi_mean)),2);

    % plot        
       bins = mean(brackets,2);
       
       ix = ixFs;  
       ax = scatter(vbl(ix),cellDep(ix));
           ax.MarkerFaceColor = fsc;
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 20;
           
       ax = scatter(vbl(ix(ixC_fs(1))),cellDep(ix(ixC_fs(1))));
           ax.MarkerFaceColor = 'k';
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 30;           
       ax = scatter(vbl(ix(ixC_fs(2))),cellDep(ix(ixC_fs(2))));
           ax.MarkerFaceColor = 'k';
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 30;  
       ax = scatter(vbl(ix(ixC_fs(3))),cellDep(ix(ixC_fs(3))));
           ax.MarkerFaceColor = 'k';
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 30;  
                     
       ax = plot(hdi_mean, bins);
       ax.LineWidth = 1.5;
       ax.Color = fsc;
           
       ax = gca;
       xlim([-0.01 0.12])
       ylim([-0.85 0.05])
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.TickLength = [0.03 0.025];

        xlabel(['HD info (bits/spike)']);
        ylabel(['Dist. from top (mm)']);
        axis(ax, 'square');
        
subplot(2,4,6)    
    ix = ixFs;
    tc = hAll(:,ix(ixC_fs(1)));  
    tc = tc./max(tc);
    ax = polarplot(b,tc);
        ax.Color = fsc;
        ax.LineWidth = 2;
        ax = gca;
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'bottom';
        ax.RColor = [0.75 0.75 0.75];
        thetaticks(0:90:270);
        ax.RColorMode = 'manual';
        ax.RColor = [0.75 0.75 0.75];
        ax.ThetaColorMode = 'manual';
        ax.ThetaColor = [0.75 0.75 0.75];
        ax.GridAlpha = 1;
        ax.ThetaTickLabel = {'';'';'';''};
        rticks([]);  
        ax.LineWidth = 2;
        
 subplot(2,4,7)   
    ix = ixFs;
    tc = hAll(:,ix(ixC_fs(2)));  
    tc = tc./max(tc);
    ax = polarplot(b,tc);
        ax.Color = fsc;
        ax.LineWidth = 2;
        ax = gca;
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'bottom';
        ax.RColor = [0.75 0.75 0.75];
        thetaticks(0:90:270);
        ax.RColorMode = 'manual';
        ax.RColor = [0.75 0.75 0.75];
        ax.ThetaColorMode = 'manual';
        ax.ThetaColor = [0.75 0.75 0.75];
        ax.GridAlpha = 1;
        ax.ThetaTickLabel = {'';'';'';''};
        rticks([]);  
        ax.LineWidth = 2;
    
  subplot(2,4,8)   
    ix = ixFs;
    tc = hAll(:,ix(ixC_fs(3))); 
    tc = tc./max(tc);
    ax = polarplot(b,tc);
        ax.Color = fsc;
        ax.LineWidth = 2;
        ax = gca;
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'bottom';
        ax.RColor = [0.75 0.75 0.75];
        thetaticks(0:90:270);
        ax.RColorMode = 'manual';
        ax.RColor = [0.75 0.75 0.75];
        ax.ThetaColorMode = 'manual';
        ax.ThetaColor = [0.75 0.75 0.75];
        ax.GridAlpha = 1;
        ax.ThetaTickLabel = {'';'';'';''};
        rticks([]);  
        ax.LineWidth = 2;
                
%% Figure 2: Waveforms

figure (2), clf
set(gcf,'color','w');
tiledlayout(1,6);


step = 3;
yVec = [0:step:63*step]*-1;

for nC = 1:length(ixC_hd)
    
    nexttile
    hold on

    ix = ixEx;
    wf = meanWaveforms{ix(ixC_hd(nC))};
    wf = wf./max(wf,[],'all');
    [totRow,totCol] = size(wf);
    xVal = 1:totCol;

    for nR  = 1:totRow
        yVal = wf(nR,:) + yVec(nR);
        ax = plot(xVal,yVal);
        ax.Color = hdc;
        ax.LineWidth = 1.5;
    end

    ax = gca; 
    set(ax, 'visible', 'off')
    xlim([5 30])
    ylim([yVec(end)-step step])
end


for nC = 1:length(ixC_fs)
    
    nexttile
    hold on

    ix = ixFs;
    wf = meanWaveforms{ix(ixC_fs(nC))};
    wf = wf./max(wf,[],'all');
    [totRow,totCol] = size(wf);
    xVal = 1:totCol;

    for nR  = 1:totRow
        yVal = wf(nR,:) + yVec(nR);
        ax = plot(xVal,yVal);
        ax.Color = fsc;
        ax.LineWidth = 1.5;
    end

    ax = gca; 
    set(ax, 'visible', 'off')
    xlim([5 30])
    ylim([yVec(end)-step step])
end

    
%% Figure 3: All cells (coloured by identity)

figure (3), clf
set(gcf,'color','w');
tiledlayout('flow');

    ix = ixGd;
    b = deg2rad(0:359);
    [~,sIx] = sort(cellDep(ix));
    sIx = flipud(sIx);
    ix = ix(sIx);
    
   
    for nC = 1:length(ix)
        
        nexttile
        
        tc = hAll(:,ix(nC));
        ax = polarplot(b,tc);
        if ismember(ix(nC),ixFs)
            ax.Color = fsc;
        elseif ismember(ix(nC),ixHd)
                 ax.Color = hdc;   
        else
            ax.Color = [0.5 0.5 0.5];  
        end
        
        ax.LineWidth = 2;
        ax = gca;
        ax.ThetaDir = 'counterclockwise';
        ax.ThetaZeroLocation = 'left';
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
    end
    

    
    
    