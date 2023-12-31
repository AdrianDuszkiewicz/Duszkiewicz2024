function Duszkiewicz2024_Figure_CueRotation

% This script reporoduces panels from:
% Figure 2 and Extended Data Figure 2

% Dependencies: 
% TStoolbox
% boundedline

%TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.


%% Parameters and variables

clear all

isHd = [];
isFs = [];
isEx = [];
whichRec = [];
isCw = [];
isCwR = [];

tcAll = [];
tcAllMean = [];
ccAll = [];
ccAll_sh = [];

%% Load data

dataset = List2Cell('dataset_cuerot.list');

for ii = 1:length(dataset) 
    
    fprintf(['Uploading recording ', dataset{ii},'\n'])
    load(fullfile(dataset{ii}, 'Data','CellTypes'));
    load(fullfile(dataset{ii},'Analysis','TuningCurvesCue'));
    load(fullfile(dataset{ii},'Data','RotDirection'));
    
    isHd = [isHd; hd];
    isFs = [isFs; fs];
    isEx = [isEx; ex];
    
    whichRec = [whichRec; repmat(ii,length(hd),1)]; % session ID    
    tcAll = [tcAll allTCcue]; % tuning curves
    tcAllMean = [tcAllMean meanTCcue];
    ccAll = [ccAll allCCcue]; % TC autocorrs
    ccAll_sh = [ccAll_sh allCCcue_rev];
    
    isCw = [isCw; isClockw];
    isCwR = [isCwR; repmat(isClockw',length(hd),1)];

end
       
[nbins,totC,totRnd] = size(ccAll);   

 %% Cell types
 
 ixHd = find(isHd);
 ixFs = find(isFs);
 
%% Calculations - rotation


% Degree of rotation (real)
[~,maxIx] = max(ccAll,[],1);
maxIx = squeeze(maxIx);
maxIx = maxIx-1;
maxIx(maxIx > 180) = maxIx(maxIx > 180) - 360;
maxIx = deg2rad(maxIx);

% mean rotation of HD cells (real)
totRec = length(unique(whichRec));
rotHD = nan(totRnd,totRec);

for nRec = 1:totRec
    ix = find(isHd & whichRec == nRec);
    rotHD(:,nRec) = CircularMean(maxIx(ix,:));
end
    
rotHD(rotHD > pi) = rotHD(rotHD > pi) - 2*pi;

% Calculate mean rotation diff between HD and FS (real)
rotDiff = nan(totC,1);
ix = ixFs;
for nC = 1:totC
    nRec = whichRec(nC);
    rc = maxIx(nC,:);
    rhd = rotHD(:,nRec)';
    rotDiff(nC) = CircularMean(angdiff(rc,rhd)');
end
rotDiff(rotDiff > pi) = rotDiff(rotDiff > pi) - 2*pi;

% Degree of rotation (rev)
[~,maxIx_sh] = max(ccAll_sh,[],1);
maxIx_sh = squeeze(maxIx_sh);
maxIx_sh = maxIx_sh-1;
maxIx_sh(maxIx_sh > 180) = maxIx_sh(maxIx_sh > 180) - 360;
maxIx_sh = deg2rad(maxIx_sh);
    
% mean rotation of HD cells (rev)
totRec = length(unique(whichRec));
rotHD_sh = nan(totRnd,totRec);

for nRec = 1:totRec
    ix = find(isHd & whichRec == nRec);
    rotHD_sh(:,nRec) = CircularMean(maxIx_sh(ix,:));
end
    
rotHD_sh(rotHD_sh > pi) = rotHD(rotHD_sh > pi) - 2*pi;

% Calculate mean rotation diff between HD and FS (rev)
rotDiff_sh = nan(totC,1);
ix = ixFs;
for nC = 1:totC
    nRec = whichRec(nC);
    rc = maxIx_sh(nC,:);
    rhd = rotHD_sh(:,nRec)';
    rotDiff_sh(nC) = CircularMean(angdiff(rc,rhd)');
end
rotDiff_sh(rotDiff_sh > pi) = rotDiff_sh(rotDiff_sh > pi) - 2*pi;

% calculate mean cw and ccw rotation for each cell

meanRot = nan(totC,2);

ix = find(isCwR(:,1) == 1);
meanRot(ix,1) = mean(maxIx(ix,1:2:end),2);
meanRot(ix,2) = mean(maxIx(ix,2:2:end),2);
ix = find(isCwR(:,1) == 0);
meanRot(ix,2) = mean(maxIx(ix,1:2:end),2);
meanRot(ix,1) = mean(maxIx(ix,2:2:end),2);

t = deg2rad([0:1:nbins-1]);
maxBasis = 180;
maxComp = 10;
maxComp_score = 3;

%% Calculations - Fourier 

% Fourier basis
basis = [];
for n = 1:maxBasis
    basis = [basis;cos(n*t) + 1i*sin(n*t)];
end    

vblFT = tcAllMean(:,:,1); 
    % FFT code
    ft   = basis * vblFT/length(t); %complex Fourier components
    ftN  = ft .* conj(ft);% ./ repmat(var(vblFT)/2,[size(ft,1) 1]);
    ftN = abs(ftN);
    % cut FT to N first components and normalize
    ftN = ftN(1:maxComp,:)./repmat(sum(ftN(1:maxComp,:),1),maxComp,1); 
    ft1 = ftN;
    symScore1 = (maxComp_score*max(ftN,[],1) - 1) ./ (maxComp_score - 1);
    
vblFT = tcAllMean(:,:,2); 
    % FFT code
    ft   = basis * vblFT/length(t); %complex Fourier components
    ftN  = ft .* conj(ft);% ./ repmat(var(vblFT)/2,[size(ft,1) 1]);
    ftN = abs(ftN);
    % cut FT to N first components and normalize
    ftN = ftN(1:maxComp,:)./repmat(sum(ftN(1:maxComp,:),1),maxComp,1); 
    ft2 = ftN;
    symScore2 = (maxComp_score*max(ftN,[],1) - 1) ./ (maxComp_score - 1);
    
    
%% Repeated measures ANOVA for FT spectra
% via https://uk.mathworks.com/matlabcentral/answers/124353-does-fitrm-ranova-support-within-subject-models-without-between-subject-factors

 % for two halves   
    vbl = [ft1(:,ixFs)' ft2(:,ixFs)'];
    grp = ones(size(ixFs));

    tbl = array2table([grp vbl]);
    tbl.Properties.VariableNames = {'CellType','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'}; % add variable names as appropriate
    Meas = table(categorical([ones(1,10) ones(1,10)+1])',categorical([1:10 1:10])','VariableNames',{'a','b'});

    rm = fitrm(tbl,'a1-b10~1','WithinDesign',Meas); % add variable names as appropriate
    anTbl12 = ranova(rm,'WithinModel','a*b-1')    
    

%%% FIGURES %%%

fsc = [0.9 0.2 0.2];
hdc = [0 0.3 0.8];

fsc_light = [1 0.4 0.4];
fsc_dark = [0.6 0.1 0.1];

fs = 11;
lw = 1;
%% Figure 1: Stats for FS 

figure (1), clf
set(gcf, 'Color', 'w');
tiledlayout(2,2);

binsize = 20; % bin size in degrees
edges = 0:binsize:360;
edges = deg2rad(edges - binsize/2);


nexttile
    ix = ixFs;
    vbl = meanRot(:,1);
    ax = polarhistogram(vbl(ix),edges);
    hold on
    %ax.Normalization = 'probability';
    ax.FaceColor = fsc;
    ax.EdgeColor = 'k';
    ax.FaceAlpha = 1;
    ix = ixFs;
    vbl = meanRot(:,2);
    ax = polarhistogram(vbl(ix),edges);
    %ax.Normalization = 'probability';
    ax.FaceColor = 'w';
    ax.EdgeColor = fsc;
    ax.LineWidth = 1;
    ax.FaceAlpha = 1;
    
    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top';
    ax.GridColor = 'k';
    thetaticks(0:90:360);
    %rticks([200]);
    ax.RColorMode = 'manual';
    ax.RColor = 'k';
    ax.ThetaColorMode = 'manual';
    ax.ThetaColor = 'k';
    ax.GridAlpha = 0.3;
    %ax.ThetaTickLabel = {'0' '90' '180' '270'};  
    ax.FontSize = fs;
    ax.LineWidth = lw;
    rlim([0 35])
    
nexttile
    ix = ixHd;
    vbl = meanRot(:,1);
    ax = polarhistogram(vbl(ix),edges);
    hold on
    %ax.Normalization = 'probability';
    ax.FaceColor = hdc;
    ax.EdgeColor = 'k';
    ax.FaceAlpha = 1;
    vbl = meanRot(:,2);
    ax = polarhistogram(vbl(ix),edges);
    %ax.Normalization = 'probability';
    ax.FaceColor = 'w';
    ax.EdgeColor = hdc;
    ax.LineWidth = 1;
    ax.FaceAlpha = 1;
    
    ax = gca;
    ax.ThetaDir = 'clockwise';
    ax.ThetaZeroLocation = 'top';
    ax.GridColor = 'k';
    thetaticks(0:90:360);
    %rticks([2000]);
    ax.RColorMode = 'manual';
    ax.RColor = 'k';
    ax.ThetaColorMode = 'manual';
    ax.ThetaColor = 'k';
    ax.GridAlpha = 0.3;
    %ax.ThetaTickLabel = {'0' '90' '180' '270'};  
    ax.FontSize = fs;
    ax.LineWidth = lw;
    rlim([0 250])

nexttile
hold on
    ix = ixFs;
    vbl1 = abs(rad2deg(rotDiff(ix)));
    vbl2 = abs(rad2deg(rotDiff_sh(ix)));
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 15;
    h2.BinWidth = 15;
    %h1.Normalization = 'probability';
    %h2.Normalization = 'probability'; 
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.FaceColor = 'none';
    h1.EdgeColor = fsc;
    h2.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h1.LineWidth = 1;
    h2.LineWidth = 1;
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
    ax.TickLength = [0.03 0.025];
    ax.XTick = [0 90 180];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    xlim([0 180]);
    axis(ax, 'square');
    xlabel(['Rot. diff. vs HD'])
    ylabel(['No. FS cells'])
    [p,h,stats] = signrank(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
       
nexttile
    ix = find(isFs);
    vbl = ft1;
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:maxComp,m(1:maxComp),s(1:maxComp),'cmap',fsc_light);
    vbl = ft2;
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:maxComp,m(1:maxComp),s(1:maxComp),'cmap',fsc_dark);

    xlabel('Fourier component')
    ylabel('Norm. power')   
    ax = gca;
%      ax.XScale = 'log';
%      ax.YScale = 'log';
    ax.XTick = [1:maxComp];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.YTick = [0:0.1:0.5];
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylim([0 0.5])
    ylim([0 0.5]);
    %xlim([0 11]);
       
    



