function Duszkiewicz2024_Figure_SingleRec_LFP

% This script reporoduces panels from Figure:
% Extended Data Figure 9
% Run in folder A3706-200313 to reproduce manuscript figure
% 

% Dependencies: 
% TStoolbox
% Buzcode

% TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

clear all 

%% parameters

smoothTC = 3; % smoothing of tuning curves (3 used in Duszkiewicz et al., 2024)

%% load data

load(fullfile('Data','SpikeData')); 
load(fullfile('Data','BehavEpochs')); 
load(fullfile('Data','CellTypes')); 
load(fullfile('Analysis','HdTuning_moveEp'));

[~, foldername, ~] = fileparts(pwd);
load(fullfile('Sleep',[foldername '.EMGFromLFP.LFP.mat'] ));
load(fullfile('Sleep',[foldername '.SleepState.states.mat'] ));
load(fullfile('Sleep',[foldername '.SleepScoreLFP.LFP.mat'] ));
load(fullfile('Sleep',[foldername '.sessionInfo.mat'] ));

% Get EMG
emgTs = EMGFromLFP.timestamps;
emgData = EMGFromLFP.data;
emgData = rescale(emgData);
emg = tsd(emgTs,emgData);

% Get lFP
chan = SleepScoreLFP.THchanID+1; % Get the best theta channel from sleep scoring (+1 because of 0-indexing)
lfp = bz_GetLFP(chan,'downsample',5);
dt = lfp.data;
lfp = tsd(lfp.timestamps,dt);

% get REM episodes
rem = SleepState.ints.REMstate;

% tuning curves and cell type
tcAll = hAll(:,:,smoothTC+1);
ixCells = find(hd);
%% Figures

lw = 1;
fs = 11;

%% Figure 1: Decoding WAKE

figure (1), clf;
set(gcf,'Color','w');

    offset = 100; % circshift the raster plot so it fits in one window
    strt = 1779; % offset of the epoch

    epStart = Start(wake1Ep)+strt;
    epEnd = Start(wake1Ep)+ strt+5;

    ep = intervalSet(epStart, epEnd); % pick a nice interval after checking the whole epoch
    lfpRes = Restrict(lfp,ep);

%lfp
subplot(3,4,[1 2]);
    yVal = Data(lfpRes);
    ax = plot(Range(lfpRes),yVal);
    ax.Color = 'k';
    ax.LineWidth = 1;
       
    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YColor = 'k';
    ax.Color = 'none';
    ax.XColor = 'none';
    box off
    %ylim([-2000 4000])
    xlim([epStart epEnd])
    ylim([-2600 2600])

% EMG 
subplot(3,4,[5 6]);

    emgRes = Restrict(emg,ep);  
    dt = gaussFilt(Data(emgRes),2);
    ax = bar(Range(emgRes),dt); 
    ax.BarWidth = 1;
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
subplot(3,4,[9 10]);
    
    axis xy
    ax = gca;
    ylim([-48 0])
    [~,mu] =  max(tcAll,[],1);
     mu = wrapTo360(mu(ixCells) + offset);
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
    %ax.Visible = 'off';
    xlim([epStart epEnd])
    ylim([-105 -10])
        

% REM

strt = rem(2,1)+1;
% Tuned for A3706 
epStart = Start(sleep1Ep)+strt;
epEnd = Start(sleep1Ep)+ strt+5;

ep = intervalSet(epStart, epEnd); 
lfpRes = Restrict(lfp,ep);

%LFP
subplot(3,4,[3 4]);

    yVal = Data(lfpRes);
    ax = plot(Range(lfpRes),yVal);
    ax.Color = 'k';
    ax.LineWidth = 1;
       
    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YColor = 'k';
    ax.Color = 'none';
    ax.XColor = 'none';
    box off
    xlim([epStart epEnd])
    ylim([-1800 1800])
    
% EMG 
subplot(3,4,[7 8]);

    emgRes = Restrict(emg,ep);  
    dt = gaussFilt(Data(emgRes),2);
    ax = bar(Range(emgRes)+0.25,dt); 
    ax.BarWidth = 1;
    ax.FaceColor = [0.8 0.8 0.8];
    ax.EdgeColor = [0.8 0.8 0.8];
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
subplot(3,4,[11 12]);
    
    axis xy
    ax = gca;
    ylim([-48 0])
    [~,mu] =  max(tcAll,[],1);
    mu = wrapTo360(mu(ixCells) + offset);
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

%print('-painters','-dsvg','ExampleLFP');