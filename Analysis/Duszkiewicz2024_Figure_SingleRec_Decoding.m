function Duszkiewicz2024_Figure_SingleRec_Decoding

% This script reporoduces panels from:
% Extended Data Figure 1
% Run in folder A3705-200306 to reproduce the manuscript figure

% Dependencies: 
% TStoolbox

%TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


% input parameters

 binsize = 0.1; % bin size for decoding in seconds
 nBins = 360; % number of angular bins for tuning curves
 sdSmoothTC = 3; % smoothing widnow for tuning curves
 sdSmoothBD = 2; % smoothing window for spike trains
 velTh = 2;

%% Load data

load(fullfile('Data','SpikeData'),'S');
load(fullfile('Data','BehavEpochs'));
load(fullfile('Data','Angle'));
load(fullfile('Data','CellTypes'));
load(fullfile('Data','Velocity'));


% Restrict to the first half of wake1Ep (only tracked bit)
rAng = Restrict(ang,wake1Ep);
ep = Range(rAng);
ep = intervalSet(ep(1),ep(end));
rAng_rev = tsd(Range(rAng),flipud(Data(rAng)));

% restrict to velocity 

epVel = thresholdIntervals(vel,velTh,'Direction','Above');   
    
eps = regIntervals(ep,2);
eps{1} = intersect(eps{1},epVel); 
baseEp = eps{1};
decEp = eps{2};

%% Run decoding
totC = length(S);

% compute tuning curves
tcAll = nan(nBins,totC);
for nC = 1:totC 
    [h,b] = HeadDirectionField(S{nC},rAng,baseEp,nBins,sdSmoothTC);
    tcAll(:,nC) = h(1:end-1);    
end

b = b(1:end-1);

% bin spikes
Q = MakeQfromS(S,binsize);
Q = Restrict(Q,wake1Ep);
dQ = Data(Q);
dQ = gaussFilt(dQ,sdSmoothBD,0);

% Run decoder
ixCells = find(fs == 1);
[thetaEst,thetaP,thetaMat] = BayesReconstruction_1D(tcAll(:,ixCells),dQ(:,ixCells),b,binsize);
decAng_FS = tsd(Range(Q),thetaEst);
decPval_FS = tsd(Range(Q),thetaP);
decMat_FS = tsd(Range(Q),thetaMat);
er1 = Restrict(decAng_FS,decEp);
er2 = Restrict(rAng,er1);
errDec_FS = angdiff(Data(er1),Data(er2));
medianErrorFS = rad2deg(nanmedian(abs(errDec_FS)));

ixCells = find(hd == 1);
[thetaEst,thetaP,thetaMat] = BayesReconstruction_1D(tcAll(:,ixCells),dQ(:,ixCells),b,binsize);
decAng_HD = tsd(Range(Q),thetaEst);
decPval_HD = tsd(Range(Q),thetaP);
decMat_HD = tsd(Range(Q),thetaMat);
er1 = Restrict(decAng_HD,decEp);
er2 = Restrict(rAng,er1);
errDec_HD = angdiff(Data(er1),Data(er2));
medianErrorHD = rad2deg(nanmedian(abs(errDec_HD)));

% save decoder data

%SaveAnalysis(pwd,'DecodedAngle_FS',{decAng_FS; decPval_FS; decMat_FS; medianErrorFS; decAng_HD; decPval_HD; decMat_HD; medianErrorHD; medianErrorFS_sh},{'decAng_FS';'decPval_FS';'decMat_FS';'medianErrorFS'; 'decAng_HD';'decPval_HD';'decMat_HD';'medianErrorHD'; 'medianErrorFS_sh'});
    

%return 

%%% FIGURES %%% 

%% Figures

lw = 1;
fs = 11;

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];

%% Figure 1: Decoding (fragment) 

figure (1), clf;
set(gcf,'Color','w')


window = 15; % displayed decoding window in seconds
start = 560; % start of decoded epoch
epStart = Start(decEp)+start;
epEnd = Start(decEp)+start+window;

offset = 115; % angular offset of HD so it fits the window

ep = intervalSet(epStart, epEnd); 
decAngRes_FS = Restrict(decAng_FS,ep);
decAngRes_HD = Restrict(decAng_HD,ep);
realAngRes = Restrict(rAng,ep);


hold on
    
    yVal = rad2deg(wrapToPi(Data(realAngRes)+deg2rad(offset)));
    ax = plot(Range(realAngRes),yVal);
    ax.Color = 'k';
    ax.LineWidth = 2;

    yVal = rad2deg(wrapToPi(Data(decAngRes_FS)+deg2rad(offset)));
    yVal = gaussFilt(yVal,0,0);
    ax = plot(Range(decAngRes_FS),yVal);
    ax.Color = fsc;
    ax.LineWidth = 2;
    
    yVal = rad2deg(wrapToPi(Data(decAngRes_HD)+deg2rad(offset)));
    yVal = gaussFilt(yVal,0,0);
    ax = plot(Range(decAngRes_HD),yVal);
    ax.Color = hdc;
    ax.LineWidth = 2;

    ax = gca;
    ax.FontSize = 11;
    ax.LineWidth = 1;
    ax.YColor = 'k';
    ax.Color = 'none';
    ax.YTick = [-180 -90 0 90 180];
    ylabel ('Head direction (deg)');
    ax.XColor = 'none';
    box off
    ylim([-180 180])
    xlim([epStart epEnd])




