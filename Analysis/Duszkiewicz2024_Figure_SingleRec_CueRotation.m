function Duszkiewicz2024_Figure_SingleRec_CueRotation

% This script reporoduces panels from:
% Extended Data Figure 2
% Run in folder A3707-200318 to reproduce manuscript figure

% Dependencies: 
% TStoolbox

%TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


clear all

%% Collect data

load(fullfile('Data','CellTypes'));
load(fullfile('Analysis','TuningCurvesCue'));
load(fullfile('Analysis','BayesianDecoding_coarse')); 
load(fullfile('Analysis','BayesianDecoding_fine')); 
load(fullfile('Analysis','Sigmoids'));
load(fullfile('Data','CueEpochs'));

b = deg2rad(1:360);
       
[~,~,totR] = size(allCCcue);   
 %% Cell types
 
 ixHd = find(hd);
 ixFs = find(fs);
 
%% Calculate rotation 

[~,maxIx] = max(allCCcue,[],1);
maxIx = squeeze(maxIx);
maxIx(maxIx > 180) = maxIx(maxIx > 180) - 360;
     
%%% FIGURES %%%

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];

fs = 11;
lw = 1;

%% Figure 1: Rotation for FS and HD

figure (1), clf
set(gcf,'Color','w');

[~,maxIx] = max(allCCcue,[],1);
maxIx = squeeze(maxIx);
maxIx(maxIx > 180) = maxIx(maxIx > 180) - 360;

ixC1 = 5;
ixC3 = 2;

subplot(2,6,[1 2 7 8])
hold on

    xVec1 = 1:totR;

    ax = plot([0 16.5],[90 90],':');
    ax.LineWidth = 1.5;
    ax.Color = 'k';
    ax = plot([0 16.5],[-90 -90],':');
    ax.LineWidth = 1.5;
    ax.Color = 'k';
    
    for nR  = 1:totR
        ix = ixHd;
        ax = scatter(repmat(xVec1(nR)-0.2,length(ix),1),maxIx(ix,nR));
        ax.MarkerEdgeColor = hdc; 
        ax.SizeData = 5;
        ix = ixFs;
        ax = scatter(repmat(xVec1(nR)+0.2,length(ix),1),maxIx(ix,nR));
        ax.MarkerEdgeColor = fsc; 
        ax.SizeData = 5;
    end
    
    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.XTick = [1:totR];    
    ax.YTick = [-180 -90 0 90 180];
    ax.XTickLabel = {};
    ax.TickLength = [0.03 0.025];
    ylim([ -200 200])
    xlim([0 16.5])
    xlabel('Cue rotation epoch')
    ylabel('Rotation (deg)')
 

subplot(2,6,4)
    tc1 = meanTCcue(:,ixHd(ixC3),1);
    ax = polarplot(b,tc1);
    ax.Color = hdc;
    ax.LineWidth = 2;
    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top';
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
    
subplot(2,6,5)
    tc2 = meanTCcue(:,ixHd(ixC3),2);
    ax = polarplot(b,tc2);
    ax.Color = hdc;
    ax.LineWidth = 2;
    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top';
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
    
subplot(2,6,6)
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
    ax.XTickLabel = [{}];
    ax.TickLength = [0.03 0.025];
    xlim([-180 180])
    axis(ax, 'square');
    

subplot(2,6,10)
    tc1 = meanTCcue(:,ixFs(ixC1),1);
    ax = polarplot(b,tc1);
    ax.Color = fsc;
    ax.LineWidth = 2;
    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top';
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
    
subplot(2,6,11)
    tc2 = meanTCcue(:,ixFs(ixC1),2);
    ax = polarplot(b,tc2);
    ax.Color = fsc;
    ax.LineWidth = 2;
    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top';
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
    
subplot(2,6,12)
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
    ax.XTickLabel = [{}];
    ax.TickLength = [0.03 0.025];
    xlim([-180 180])
    axis(ax, 'square');
     

%% Figure 2: Rotation example

figure (2), clf
set(gcf,'Color','w');

t = Range(errBayes);
xVal = (t - t(1))./60;
              
subplot(1,2,1)
hold on

    yVal = rad2deg(Data(errBayes));

    ax = plot([0 xVal(end)],[90 90],':r');
    ax.LineWidth = 2;
    
    ax = plot([0 xVal(end)],[0 0],':r');
    ax.LineWidth = 2;
    
    ax = plot(xVal,yVal);
    ax.Color = [0.5 0.5 0.5];
     
    ax = gca;
    ax.LineWidth = lw;
    ax.FontSize = fs;
    ax.XColor = 'none';
    ax.YTick = [-180 -90 0 90 180];
    ax.TickLength = [0.03 0.025];
    box off
    ylim([-180 180]);
    
subplot(1,2,2)
hold on

    noEp = 5;
    ep = intervalSet(cueEp(noEp,1)-10,cueEp(noEp,2));
    err = Restrict(errBayesFine,ep);
    yVal = rad2deg(Data(err));
    xVal = Range(err) - cueEp(noEp,1);
    
    ax = plot(xVal,yVal);
    ax.Color = [0.5 0.5 0.5];
    ax.LineWidth = 2;
    

    ax = plot(sigTS(:,noEp)-cueEp(noEp,1),rad2deg(sigVal(:,noEp)));
    ax.Color = 'r';
    ax.LineWidth = 1;
    
    ax = gca;
    ax.LineWidth = lw;
    ax.FontSize = fs;
    
    ylim([-50 180]);
    xlim([-10 30]);
    ax.YTick = [-180 -90 0 90 180];


return 


