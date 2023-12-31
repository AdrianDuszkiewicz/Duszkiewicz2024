function Duszkiewicz2024_Figure_SingleRec_IsomapRaster

% This script reporoduces panels from:
% Figure 7, Extended Data Figure 10
% run in folder A3706-200313 to reproduce the manuscript figure

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

%% parameters

clear all 
smoothTC = 3; % smoothing of tuning curves (3 used in Duszkiewicz et al., 2024)

%% load data

load(fullfile('Data','SpikeData')); 
load(fullfile('Data','BehavEpochs')); 
load(fullfile('Data','CellTypes')); 
load(fullfile('Data','Angle'));
load(fullfile('Analysis','HdTuning_moveEp')); 
load(fullfile('Analysis','Isomap_WK_HdOnly'));
QisoWk = Q;
mappingWk = mapping;
load(fullfile('Analysis','Isomap_REMvsWK'));
QisoREM = Qrw;
mappingREM = mapping;
[~, foldername, ~] = fileparts(pwd);
load(fullfile('Sleep',[foldername '.EMGFromLFP.LFP.mat'] ));
load(fullfile('Sleep',[foldername '.SleepState.states.mat'] ));

% Get EMG
emgTs = EMGFromLFP.timestamps;
emgData = EMGFromLFP.data;
emgData = rescale(emgData);
emg = tsd(emgTs,emgData);

% get REM episodes
rem = SleepState.ints.REMstate;
remEp = intervalSet(rem(:,1),rem(:,2));

% get tuning curves
tcAll = hAll(:,:,smoothTC+1);
ixCells = find(hd);

%% Get Isomap angle 
% for WAKE
rAng = Restrict(ang,QisoWk); % angle values in the middle of Q bins - ang value sampled at much higher freq.
isoAngWk = deg2rad(atan2d(mappingWk(:,1),mappingWk(:,2))); % compute the angle from centre for each point

% Determine whether the angle matches clockwise or counterclockwise:
dr = Data(rAng);
di = isoAngWk;

ix = find(isnan(dr));
dr(ix) = []; % remove nans for now
di(ix) = [];

err_for = nan(length(dr),1);
err_rev = nan(length(dr),1);

for nS = 1:length(dr) % first determine error between real and iso angle
    err_for(nS) = angdiff(dr(nS),di(nS)); %in one direction
    err_rev(nS) = angdiff(dr(nS),-di(nS)); % and in reverse
end
    
var_for = circ_var(err_for); % requires CircStat toolbox (google it)
var_rev = circ_var(err_rev); % calculate circular variance of both errors

if var_rev < var_for
    isoAngWk = -isoAngWk;
end
isoAngWk = tsd(Range(QisoWk),isoAngWk);

% for REM
rAng = Restrict(ang,QisoREM); % angle values in the middle of Q bins - ang value sampled at much higher freq.
isoAngREM = deg2rad(atan2d(mappingREM(:,1),mappingREM(:,2))); % compute the angle from centre for each point

% Determine whether the angle matches clockwise or counterclockwise:
dr = Data(rAng);
di = isoAngREM;

ix = find(isnan(dr));
dr(ix) = []; % remove nans for now
di(ix) = [];

err_for = nan(length(dr),1);
err_rev = nan(length(dr),1);

for nS = 1:length(dr) % first determine error between real and iso angle
    err_for(nS) = angdiff(dr(nS),di(nS)); %in one direction
    err_rev(nS) = angdiff(dr(nS),-di(nS)); % and in reverse
end
    
var_for = circ_var(err_for); % requires CircStat toolbox (google it)
var_rev = circ_var(err_rev); % calculate circular variance of both errors

if var_rev < var_for
    isoAngREM = -isoAngREM;
end
isoAngREM = tsd(Range(QisoREM),isoAngREM);

%% Figures

lw = 1;
fs = 11;

%% Figure 1: Decoding WAKE

figure (1), clf;
set(gcf,'Color','w')

    epStart = 4253.3;
    epEnd = 4302.3;

    ep = intervalSet(epStart, epEnd); 
    decAngRes = Restrict(isoAngWk,ep);
    realAngRes = Restrict(ang,ep);

%angle
subplot(3,2,[1 2]);
hold on
    
    offset1 = -140;
    offset2 = 99; 

    yVal = rad2deg(wrapToPi(Data(realAngRes)+deg2rad(offset2)));
    ax = plot(Range(realAngRes),yVal);
    ax.Color = 'k';
    ax.LineWidth = 2;

    yVal = rad2deg(wrapToPi(Data(decAngRes)+deg2rad(offset1)));
    yVal = gaussFilt(yVal,0,0);
    ax = plot(Range(decAngRes),yVal);
    ax.Color = 'r';
    ax.LineWidth = 2;

    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YColor = 'k';
    ax.Color = 'none';
    ax.YTick = [-180 -90 0 90 180];
    ylabel ('Head direction (deg)');
    ax.XColor = 'none';
    box off
    ylim([-180 180])
    xlim([epStart epEnd])

% EMG 
subplot(3,2,[3 4]);

    emgRes = Restrict(emg,ep);  
    dt = gaussFilt(Data(emgRes),2);
    ax = bar(Range(emgRes),dt); 
    ax.FaceColor = [0.8 0.8 0.8];
    ax.EdgeColor = [0.8 0.8 0.8];
    ax.LineWidth = 2;
    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YColor = 'k';
    ax.Color = 'none';
    ax.YTick = [0 1];
    ylabel ('EMG (a.u.)');
    ax.XColor = 'none';
    box off
    xlim([epStart epEnd])
    ylim([0 1])

% raster 
subplot(3,2,[5 6]);
    
    axis xy
    ax = gca;
    ylim([-48 0])
    [~,mu] =  max(tcAll,[],1);
    mu = wrapTo360(mu(ixCells) + 180 + offset2);
    [~,sIx] = sort(mu);      
    mu = mu(sIx);
    mu = mu./360;
    colIx = round(mu*256);
    col = hsv;
    col = col(colIx,:);
    Sr = S(ixCells(sIx));

    for nC = 1:length(Sr)
        RasterPlot(Sr(nC),'TStart',epStart,'TEnd',epEnd,'YOffset',nC-1-100,'AxHandle',ax,'BarColor',{col(nC,:)},'FigureHandle',gcf)
        hold on
    end
    
    set(ax, 'Color', 'none')
    ax.YColor = 'none'; 
    ax.XColor = 'none';
    ax.XColor = [0.5 0.5 0.5];
    ax.Color = 'w';
    ax.FontSize = 10;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    ax.Visible = 'off';
    xlim([epStart epEnd])
    ylim([-105 -10])
      
%% Figure 2: Decoding SLEEP

figure (2), clf;
set(gcf,'Color','w')

    strt = rem(2,1);
    epStart = Start(sleep1Ep)+strt-20;
    epEnd = Start(sleep1Ep)+ strt+30;

    ep = intervalSet(epStart, epEnd); 
    decAngRes = Restrict(isoAngREM,ep);
    realAngRes = Restrict(ang,ep);
 
%angle
    subplot(3,2,[1 2]);
    offset1 = 28;
    offset2 = -11;
    
    hold on
    yVal = rad2deg(wrapToPi(Data(decAngRes)+deg2rad(offset1)));
    yVal = gaussFilt(yVal,0,0);
    ax = plot(Range(decAngRes),yVal);
    ax.Color = 'r';
    ax.LineWidth = 2;

    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YColor = 'k';
    ax.Color = 'none';
    ax.YTick = [-180 -90 0 90 180];
    ylabel ('Head direction (deg)');
    ax.XColor = 'none';
    box off
    ylim([-180 180])
    xlim([epStart epEnd])
            
% EMG 
subplot(3,2,[3 4]);

    emgRes = Restrict(emg,ep);  
    dt = gaussFilt(Data(emgRes),2);
    ax = bar(Range(emgRes),dt); 
    ax.FaceColor = [0.8 0.8 0.8];
    ax.EdgeColor = [0.8 0.8 0.8];
    ax.LineWidth = 2;
    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YColor = 'k';
    ax.Color = 'none';
    ax.YTick = [0 1];
    ylabel ('EMG (a.u.)');
    ax.XColor = 'none';
    box off
    xlim([epStart epEnd])
    ylim([0 1])

%raster 
subplot(3,2,[5 6]);
    
    axis xy
    ax = gca;
    ylim([-48 0])
    [~,mu] =  max(tcAll,[],1);
    mu = wrapTo360(mu(ixCells) + 180 + offset2);
    [~,sIx] = sort(mu);      
    mu = mu(sIx);
    mu = mu./360;
    colIx = round(mu*256);
    col = hsv;
    col = col(colIx,:);
    Sr = S(ixCells(sIx));

    for nC = 1:length(Sr)
        RasterPlot(Sr(nC),'TStart',epStart,'TEnd',epEnd,'YOffset',nC-1-100,'AxHandle',ax,'BarColor',{col(nC,:)},'FigureHandle',gcf)
        hold on
    end
    
    ax.YColor = 'none'; 
    ax.XColor = 'none';
    ax.Color = 'w';
    ax.FontSize = 10;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    xlim([epStart epEnd])
    ylim([-105 -10])

