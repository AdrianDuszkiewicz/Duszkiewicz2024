function Duszkiewicz2024_Figure_SingleRec_IsomapREM

% This script reporoduces panels from:
% Figure 7, Extended Data Figure 10
% run in folder A3706-200313

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

load(fullfile('Data','SpikeData'));
load(fullfile('Data','BehavEpochs'));
load(fullfile('Data','Angle'));
load(fullfile('Data','CellTypes'));
load(fullfile('Analysis','Isomap_REMvsWK'));
[~, foldername, ~] = fileparts(pwd);
load(fullfile('Sleep',[foldername '.SleepState.states.mat'] ));

binsize = median(diff(Range(Qrw)));

%% Cell types

ixHd = find(hd);
ixEx = find(ex);
ixFs = find(fs);
    
%% Calculations for tuning curves

% get real angle values
rQrw = Range(Qrw);
dQrw = Data(Qrw);
tsRem = rQrw(ixRem);

tsWk = rQrw(ixWk);
tsWk = ts(tsWk);

rAngWk = Restrict(ang,tsWk); % angle values in the middle of Q bins - ang value sampled at much higher freq.
r = Range(rAngWk);
r = r - binsize/2; % move the angle values to the beginning of Q bins, otherwise first half of the bin will have a wrong ang
d = Data(rAngWk);
rAngWk = tsd(r,d);

% Get isomap angle values
isoAng = deg2rad(atan2d(mapping(:,1),mapping(:,2))); % compute the angle from centre for each point

% Determine whether the angle matches clockwise or counterclockwise:
dr = Data(rAngWk);
di = isoAng(ixWk);

ix = find(isnan(d));
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

% ts for isoAng need to be at the beginning of Q bins because they apply to
% the whole Q bin. Otherwise Restrict will not work properly. 

isoAngWk = tsd(rQrw(ixWk)-binsize./2,isoAng(ixWk));
isoAngRem = tsd(rQrw(ixRem)-binsize./2,isoAng(ixRem));

% create intervalSets for each bin 
int = [rQrw - binsize./2 rQrw + binsize./2];
epRem = int(ixRem,:);
epRem = intervalSet(epRem(:,1),epRem(:,2));
epWk = int(ixWk,:);
epWk = intervalSet(epWk(:,1),epWk(:,2));

% compute tuning curves for each cell 
totC = length(S); 

nBins = 60;
sdSmooth = 1;
B = 2*pi*(0:1:nBins-1)'/nBins;

tcReal = nan(nBins,totC);
tcWk = nan(nBins,totC);
tcRem = nan(nBins,totC);

% Headdirectionfield crashes the script so we go manual
h0_real = hist(mod(Data(rAngWk),2*pi),B); % histograms of angle distribution
h0_wk = hist(mod(Data(isoAngWk),2*pi),B); % histograms of angle distribution
h0_rem = hist(mod(Data(isoAngRem),2*pi),B); % histograms of angle distribution

% Real tuning curves
for nC = 1:totC
    Sr = Restrict(S{nC},epWk);
    if ~isempty(Sr)
        angt = Restrict(rAngWk,Sr);
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
    Sr = Restrict(S{nC},epWk);
    if ~isempty(Sr)
        angt = Restrict(isoAngWk,Sr);
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

% Isomap REM tuning curves
for nC = 1:totC
    Sr = Restrict(S{nC},epRem);
    if ~isempty(Sr)
        angt = Restrict(isoAngRem,Sr);
        h = hist(mod(Data(angt),2*pi),B);
        if sum(h) == 0
            h = h';
        end
        dt = 1./binsize;
        h = dt*h./h0_rem;
        h(isnan(h)) = 0;
        h = gaussFiltAng(h,sdSmooth,1);
    else
        h = zeros(nBins,1);
    end
    tcRem(:,nC) = h;
end

% calculate tuning curve offset between Real and Wake 

offset = nan(totC,1);
ccAll = nan(nBins,totC);
 
for nC = 1:totC        
    cc = TCcrosscorr(tcRem(:,nC),tcWk(:,nC)); % xcorr between two tuning curves
    [~,mx] = max(cc); % max point on the xcorr is the offset
    offset(nC) = mx;
    ccAll(:,nC) = cc;
end
     
% calculate the Isomap radius in Wake and REM

totS = size(mapping,1);
radius = nan(totS,1);

for nS = 1:totS   
    pts = [0,0; mapping(nS,1),mapping(nS,2)];
    d = pdist(pts,'euclidean');
    radius(nS) = d;
end
    

%%% Figures %%%


wkc = [0.15 0.3 0.5];
rmc = [0.3 0.5 0.2];

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];
unc = [0.5 0.5 0.5];
adc = [1 0.4 1];

fsc_light = [1 0.4 0.4];
fsc_dark = [0.6 0.1 0.1];

hdc_light = [0.7 0.45 1];
hdc_dark = [0.4 0.15 0.8];

fs = 11;
lw = 1;

%% Figure 1: Isomap  


figure (1), clf
set(gcf, 'Color','w')
tiledlayout(1,1)

sd = 12;

nexttile
hold on

    ax = scatter(mapping(ixWk,1),mapping(ixWk,2));
    ax.MarkerFaceColor = wkc;
    ax.MarkerEdgeColor = 'none';
    ax.SizeData = sd;

    ax = scatter(mapping(ixRem,1),mapping(ixRem,2));
    ax.MarkerFaceColor = rmc;
    ax.MarkerEdgeColor = 'none';
    ax.SizeData = sd;
    ax = gca;
    axis(ax, 'square');
    xlim([-20 20])
    ylim([-20 20])
    %set(gca, 'visible', 'off')

%% Figure 2: Ring radius  

figure (2), clf
set(gcf, 'Color','w')
hold on

        vbl1 = radius(ixWk); 
        vbl2 = radius(ixRem); 
        m1 = median(vbl1);
        m2 = median(vbl2);
        h1 = histogram(vbl1);
        h2 = histogram(vbl2);
        %h1.BinWidth = 0.5;
        %h2.BinWidth = 0.5;
        h1.Normalization = 'probability';
        h2.Normalization = 'probability';
        h1.DisplayStyle = 'stairs';
        h2.DisplayStyle = 'stairs';
        h1.LineWidth = 1;
        h2.LineWidth = 1;
        h1.EdgeColor = wkc;
        h1.FaceColor = 'none';
        h2.EdgeColor = rmc;
        h2.FaceColor = 'none';
        mx1 = find(h1.BinEdges >= m1,1)-1;
        mx1 = h1.Values(mx1);
        mx2 = find(h2.BinEdges >= m2,1)-1;
        mx2 = h2.Values(mx2);
        ax = plot([m1 m1],[0 mx1],':');
        ax.LineWidth = 2;
        ax.Color = wkc;
        ax = plot([m2 m2],[0 mx2],':');
        ax.LineWidth = 2;
        ax.Color = rmc;
        ax = gca; 
        ylim([0 0.2]);
        ax.TickLength = [0.03 0.025];
        ax.YTick = [0 0.2];
        ax.FontSize = 11;
        ax.LineWidth = 1;
        axis(ax, 'square');
        xlabel(['Distance (a.u.)']);
        ylabel(['Prob.']);
        [p,h,stats] = signrank(vbl1,vbl2);
        ax = title(p);
        ax.FontWeight = 'normal';

%% Figure 3: Examples HD  

figure (3), clf
set(gcf, 'Color','w')


ix = ixHd; 
ixC1 = [23 8]; % HD cells in A3706

b = deg2rad(0:6:360);
totC = length(ixC1);

tiledlayout(totC,2);

for nC = 1:totC
      
nexttile  

    tc1 = tcWk(:,ix(ixC1(nC)));
    tc2 = tcRem(:,ix(ixC1(nC)));  
    tc1 = tc1./max(tc1);
    tc2 = tc2./max(tc2);    
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
    xVal = -180:6:180;
    cc = TCcrosscorr(tc1,tc2);
    cc = circshift(cc,31);
    cc = [cc(end); cc];
    ax = plot(xVal,cc);
    ax.Color = 'k';
    ax.LineWidth = 1;
    [~,maxIx] = max(cc);
    maxIx = maxIx - 29;
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
    
end

%% Figure 4: Examples FS 

figure (4), clf
set(gcf, 'Color','w')


ix = ixFs; 
ixC1 = [18 3]; % FS cells in A3706

b = deg2rad(0:6:360);
totC = length(ixC1);

tiledlayout(totC,2);

for nC = 1:totC
      
nexttile  

    tc1 = tcWk(:,ix(ixC1(nC)));
    tc2 = tcRem(:,ix(ixC1(nC)));  
    tc1 = tc1./max(tc1);
    tc2 = tc2./max(tc2);    
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
    xVal = -180:6:180;
    cc = TCcrosscorr(tc1,tc2);
    cc = circshift(cc,30);
    cc = [cc(end); cc];
    ax = plot(xVal,cc);
    ax.Color = 'k';
    ax.LineWidth = 1;
    [~,maxIx] = max(cc);
    maxIx = maxIx - 30;
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
    
end


%% Figure 5: plot all HD and FS cells 
 
figure (5), clf
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
      
        tc1 = tcWk(:,ix(nC));
        tc2 = tcRem(:,ix(nC));  
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
      
        tc1 = tcWk(:,ix(nC));
        tc2 = tcRem(:,ix(nC));  
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
    end
    