function Duszkiewicz2024_Figure_SingleRec_IsomapWake

% This script reporoduces panels from:
% Figure 7, Extended Data Figure 10
% Run in folder A3706-200313 to reproduce the manuscript figure

% Dependencies: 
% TStoolbox
% CircStat toolbox

% TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%% Load data

clear all

load(fullfile('Data','BehavEpochs'));
load(fullfile('Data','CellTypes'));
load(fullfile('Data','Angle'));
load(fullfile('Data','SpikeData'));
load(fullfile('Analysis','Isomap_WK_HdOnly'));
iso_real = mapping;
load(fullfile('Analysis','Isomap_WK_HdOnly_Shuffle'));
iso_sh = mapping;

ep = wake1Ep;
binsize = median(diff(Range(Q)));

% define cell types

ixHd = find(hd);
ixEx = find(ex);
ixFs = find(fs);
    
%% Calculations 

% get real angle values
rAng = Restrict(ang,Q); % angle values in the middle of Q bins - ang value sampled at much higher freq.

% Get isomap angle values
isoAng = deg2rad(atan2d(iso_real(:,1),iso_real(:,2))); % compute the angle from centre for each point

% Determine whether the angle matches clockwise or counterclockwise:
dr = Data(rAng);
di = isoAng;

ix = find(isnan(dr));
dr(ix) = []; % remove nans for now
di(ix) = [];

err_for = nan(length(dr),1);
err_rev = nan(length(dr),1);

for nS = 1:length(dr) % first determine error between real and iso angle
    err_for(nS) = angdiff(dr(nS),di(nS)); %in one direction
    err_rev(nS) = angdiff(dr(nS),-di(nS)); % and in reverse
end
    
var_for = circ_var(err_for); 
var_rev = circ_var(err_rev); % calculate circular variance of both errors

if var_rev < var_for
    isoAng = -isoAng;
end

rQ = Range(Q);
isoAng = tsd(rQ,isoAng);
isoAngRev = tsd(rQ,flipud(Data(isoAng)));

% compute tuning curves for each cell 
totC = length(S); 

nBins = 360;
sdSmooth = 3;
B = 2*pi*(0:1:nBins-1)'/nBins;

tcReal = nan(nBins,totC);
tcWk = nan(nBins,totC);
tcRev = nan(nBins,totC);

h0_real = hist(mod(Data(rAng),2*pi),B); % histograms of angle distribution
h0_wk = hist(mod(Data(isoAng),2*pi),B); % histograms of angle distribution

% Real tuning curves
for nC = 1:totC
    Sr = Restrict(S{nC},ep);
    if ~isempty(Sr)
        angt = Restrict(rAng,Sr);
        h = hist(mod(Data(angt),2*pi),B);
        if sum(h) == 0
            h = h'; % dumb bug where 0-matrix is vertical
        end
        dt = 1./binsize;
        h = dt*h./h0_real;
        h(isnan(h))=0;
        h = gaussFiltAng(h,sdSmooth,1);
    else
        h = zeros(nBins,1);
    end
    tcReal(:,nC) = h;
end

% Isomap wake tuning curves
for nC = 1:totC
    Sr = Restrict(S{nC},ep);
    if ~isempty(Sr)
        angt = Restrict(isoAng,Sr);
        h = hist(mod(Data(angt),2*pi),B);
        if sum(h) == 0
            h = h';
        end
        dt = 1./binsize;
        h = dt*h./h0_wk;
        h(isnan(h)) = 0;
        h = gaussFiltAng(h,sdSmooth,1);
    else
        h = zeros(nBins,1);
    end
    tcWk(:,nC) = h;
end

% Isomap shuffled tuning curves
for nC = 1:totC
    Sr = Restrict(S{nC},ep);
    if ~isempty(Sr)
        angt = Restrict(isoAngRev,Sr);
        h = hist(mod(Data(angt),2*pi),B);
        if sum(h) == 0
            h = h';
        end
        dt = 1./binsize;
        h = dt*h./h0_wk;
        h(isnan(h)) = 0;
        h = gaussFiltAng(h,sdSmooth,1);
    else
        h = zeros(nBins,1);
    end
    tcRev(:,nC) = h;
end

% calculate tuning curve offset between Real and Wake 
offset = nan(totC,1);
ix = ixHd;
for nC = 1:length(ix)       
    cc = TCcrosscorr(tcReal(:,ix(nC)),tcWk(:,ix(nC))); % xcorr between two tuning curves
    [~,mx] = max(cc); % max point on the xcorr is the offset
    offset(nC) = mx;
end
     
offset_m = round(rad2deg(CircularMean(deg2rad(offset))),0);
tcWk_off = circshift(tcWk,offset_m);

% calculate tuning curve correlation between Real and Wake 

ccAll = nan(nBins,totC);
for nC = 1:totC       
    cc = TCcrosscorr(tcReal(:,nC),tcWk_off(:,nC)); % xcorr between two tuning curves
    ccAll(:,nC) = cc;
end


% calculate the Isomap radius in Wake 
totS = size(iso_real,1);
radius = nan(totS,1);

for nS = 1:totS   
    pts = [0,0; iso_real(nS,1),iso_real(nS,2)];
    d = pdist(pts,'euclidean');
    radius(nS) = d;
end    

totS = size(iso_sh,1);
radius_sh = nan(totS,1);

for nS = 1:totS   
    pts = [0,0; iso_sh(nS,1),iso_sh(nS,2)];
    d = pdist(pts,'euclidean');
    radius_sh(nS) = d;
end   

%%% Figures %%%

wkCol = [0.15 0.3 0.5];
remCol = [0.2 0.5 0.2];

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];
unc = [0.5 0.5 0.5];
adc = [1 0.4 1];

fsc_light = [1 0.4 0.4];
fsc_dark = [0.6 0.1 0.1];

hdc_light = [0.7 0.45 1];
hdc_dark = [0.4 0.15 0.8];


lw = 1;
fs = 11;

%% Figure 1: Isomap 

figure (1), clf
set(gcf, 'Color','w')
hold on

    % this part transforms angle into a color
    dTemp = Data(rAng);
    dTemp = mat2gray(dTemp)*256; %express angle as numbers 1:256
    colIx = floor(dTemp); % express as integers
    colIx(find(colIx == 0)) = 1; % get rid of 0
    col = hsv; % color map matrix (you can use diferent one)
    col = col(colIx,:);

    ax = scatter(iso_real(:,1),iso_real(:,2),80,col,'filled');
    ax.SizeData = 2;

    ax = gca; 
    ax.LineWidth = lw;
    ax.FontSize = fs;
    axis(ax, 'square');
    xlabel('Dimension 1')
    ylabel('Dimension 2')
    xlim([-20 20]);
    ylim([-20 20]);
   % set(gca, 'visible', 'off')
   
%% Figure 2: Isomap shuffle

figure (2), clf
set(gcf, 'Color','w')
hold on

    % this part transforms angle into a color
    dTemp = Data(rAng);
    dTemp = mat2gray(dTemp)*256; %express angle as numbers 1:256
    colIx = floor(dTemp); % express as integers
    colIx(find(colIx == 0)) = 1; % get rid of 0
    col = hsv; % color map matrix (you can use diferent one)
    col = col(colIx,:);

    ax = scatter(iso_sh(:,1),iso_sh(:,2),80,col,'filled');
    ax.SizeData = 2;

    ax = gca; 
    ax.LineWidth = lw;
    ax.FontSize = fs;
    axis(ax, 'square');
    xlabel('Dimension 1')
    ylabel('Dimension 2')
    xlim([-40 40]);
    ylim([-40 40]);
   % set(gca, 'visible', 'off')
    
%% Figure 3: Ring quality and examples 

figure (3), clf
set(gcf, 'Color','w')

nexttile
hold on

    vbl1 = radius; 
    vbl2 = radius_sh; 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.5;
    h2.BinWidth = 0.5;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = wkCol;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = wkCol;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    ylim([0 0.15]);
    ax.TickLength = [0.03 0.025];
    ax.YTick = [0 0.15];
    ax.FontSize = 11;
    ax.LineWidth = 1;
    axis(ax, 'square');
    xlabel(['Distance (a.u.)'])
    ylabel(['Probability'])

   
%% Figure 4: Examples FS

figure (4), clf
set(gcf, 'Color','w')
tiledlayout(2,2)

    ixC1 = 3;
    ixC2 = 18;
    ix = ixFs;

    b = deg2rad(0:1:360);
    
nexttile  
    tc1 = tcReal(:,ix(ixC1));
    tc2 = tcWk_off(:,ix(ixC1));  
    ax = polarplot(b,[tc1; tc1(1)]);
    hold on
        ax.Color = fsc_light;
        ax.LineWidth = 2;
    ax = polarplot(b,[tc2; tc2(1)]);
        ax.Color = fsc_dark;
        ax.LineWidth = 2;
        ax = gca;
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
        tcmx = max([tc1; tc2]);
        rlim([0 tcmx*1.1]);
               
nexttile
hold on
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
    tc1 = tcReal(:,ix(ixC2));
    tc2 = tcWk_off(:,ix(ixC2));  
    ax = polarplot(b,[tc1; tc1(1)]);
    hold on
        ax.Color = fsc_light;
        ax.LineWidth = 2;
    ax = polarplot(b,[tc2; tc2(1)]);
        ax.Color = fsc_dark;
        ax.LineWidth = 2;
        ax = gca;
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
        tcmx = max([tc1; tc2]);
        rlim([0 tcmx*1.1]);
               
nexttile
hold on
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
 
%% Figure 5: Examples HD

figure (5), clf
set(gcf, 'Color','w')
tiledlayout(2,2)

    ixC1 = 23; %3 9 14 10 6
    ixC2 = 8;
    ix = ixHd;

    b = deg2rad(0:1:360);
    
nexttile  
    tc1 = tcReal(:,ix(ixC1));
    tc2 = tcWk_off(:,ix(ixC1));  
    ax = polarplot(b,[tc1; tc1(1)]);
    hold on
        ax.Color = hdc_light;
        ax.LineWidth = 2;
    ax = polarplot(b,[tc2; tc2(1)]);
        ax.Color = hdc_dark;
        ax.LineWidth = 2;
        ax = gca;
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
        tcmx = max([tc1; tc2]);
        rlim([0 tcmx*1.1]);
               
nexttile
hold on
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
    tc1 = tcReal(:,ix(ixC2));
    tc2 = tcWk_off(:,ix(ixC2));  
    ax = polarplot(b,[tc1; tc1(1)]);
    hold on
        ax.Color = hdc_light;
        ax.LineWidth = 2;
    ax = polarplot(b,[tc2; tc2(1)]);
        ax.Color = hdc_dark;
        ax.LineWidth = 2;
        ax = gca;
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
        tcmx = max([tc1; tc2]);
        rlim([0 tcmx*1.1]);
               
nexttile
hold on
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
    

%% Figure 6: plot all HD and FS cells 
 
figure (6), clf
set(gcf,'color','w');
tiledlayout('flow');

    ix = ixHd;
    vbl = ccAll(1,ix)';
    [~,sIx] = sort(vbl);
    sIx = flipud(sIx);
    vbl = vbl(sIx);
    ix = ix(sIx);

   
    for nC = 1:length(ix)
        
    nexttile
      
        tc1 = tcReal(:,ix(nC));
        tc2 = tcWk_off(:,ix(nC));  
        tc1 = tc1./max(tc1);
        tc2 = tc2./max(tc2); 
    
        ax = polarplot(b,[tc1; tc1(1)]);
        ax.Color = hdc_light;       
        ax.LineWidth = 2;
        hold on
        
        ax = polarplot(b,[tc2; tc2(1)]);
        ax.Color = hdc_dark;       
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
        ax = title(round(vbl(nC),2));
        ax.FontWeight = 'normal';
    end
    
    ix = ixFs;
    vbl = ccAll(1,ix)';
    [~,sIx] = sort(vbl);
    sIx = flipud(sIx);
    vbl = vbl(sIx);
    ix = ix(sIx);

    for nC = 1:length(ix)
        
    nexttile
      
        tc1 = tcReal(:,ix(nC));
        tc2 = tcWk_off(:,ix(nC));  
        tc1 = tc1./max(tc1);
        tc2 = tc2./max(tc2); 
    
        ax = polarplot(b,[tc1; tc1(1)]);
        ax.Color = fsc_light;       
        ax.LineWidth = 2;
        hold on
        
        ax = polarplot(b,[tc2; tc2(1)]);
        ax.Color = fsc_dark;       
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
        ax = title(round(vbl(nC),2));
        ax.FontWeight = 'normal';
    end
  



    