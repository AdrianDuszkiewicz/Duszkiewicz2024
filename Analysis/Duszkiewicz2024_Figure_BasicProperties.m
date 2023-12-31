function Duszkiewicz2024_Figure_BasicProperties

% This script reporoduces panels from: 
% Figures 1 and 2, Extended Data Figure 2


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

dataset = List2Cell('dataset_All.list');

cellSpk = [];
tpk = [];
frate = [];
meanWf = {};

whichRec = [];
hdi = [];
hdi_sh = [];
isHd = [];
isFs = [];
isGd = [];
isEx = [];

medianErrHD = [];
medianErrFS = [];
medianErrFS_sh = [];
totHD = [];
totFS = [];

errHistHD = [];
errHistFS = [];

tcAll = [];
tcAll_sh = [];

for ii = 1:length(dataset)

    %Load data
    fprintf(['Uploading recording ', dataset{ii},'\n'])
    load(fullfile(dataset{ii},'Data','Waveforms'));
    load(fullfile(dataset{ii},'Data','WaveformFeatures'));
    load(fullfile(dataset{ii},'Data','CellTypes'));
    load(fullfile(dataset{ii},'Analysis','MeanFR'));
    load(fullfile(dataset{ii},'Analysis','DecodedAngle_FS'));      
    load(fullfile(dataset{ii},'Analysis','HdTuning_moveEp'));
         
    totC = length(hd);
    
    cellSpk = [cellSpk; ~isnan(spkWidth)];
    tpk = [tpk; tr2pk];
    frate = [frate; rateS];
       
    % collect variables
    isHd = [isHd; hd];
    isFs = [isFs; fs];
    isGd = [isGd; gd];
    isEx = [isEx; ex];
    whichRec = [whichRec; repmat(ii,length(spkWidth),1)];
    hdi = [hdi; hdInfo(:,smoothTC+1)];
    hdi_sh = [hdi_sh; hdInfo_sh(:,smoothTC+1)];
    tcAll = [tcAll hAll(:,:,smoothTC+1)];
    tcAll_sh = [tcAll_sh hAll_sh(:,:,smoothTC+1)];
    
    medianErrHD = [medianErrHD; medianErrorHD];
    medianErrFS = [medianErrFS; medianErrorFS];
    medianErrFS_sh = [medianErrFS_sh; medianErrorFS_sh];
    totHD = [totHD; sum(hd)];
    totFS = [totFS; sum(fs)];
    
     % Decoding 

    edges = -180:5:180;
    bins = edges(1:end-1)+median(diff(edges))/2;
   
    counts = histcounts(rad2deg(Data(decErr_HD)),edges); %
    errHistHD = [errHistHD counts'];    
         
    counts = histcounts(rad2deg(Data(decErr_FS)),edges); %
    errHistFS = [errHistFS counts'];    
end

totC = length(isFs); % total cells

%% define groups

% Good cells 
ixGd = find(isGd == 1);
% HD cells
ixHd = find(isHd == 1 & isGd == 1);
% FS cells
ixFs = find(isFs == 1 & isGd == 1);
% Ex cells
ixEx = find(isEx == 1 & isGd == 1);

%% Calculations
    
% Get threshold for HD cells
vbl_ex = sort(hdi_sh(ixEx));
thr_ex = round(0.99*length(vbl_ex));
thr_ex = vbl_ex(thr_ex+1);
thr_ex = log10(round(thr_ex,1));

% calculate tuning curve width
tcAllSide = zeros(181,totC);
hwTc = nan(totC,1);

for nC = 1:totC    
 [~,muMax] = max(tcAll(:,nC));  
  muMax = muMax(1);   
  if ~isnan(muMax)    
    hn = circshift(tcAll(:,nC),-muMax+1);
    hn = [hn./max(hn); 1];  
    hn = mean([hn(1:181),flipud(hn(181:361))],2);
    tcAllSide(:,nC) = hn;     
    hw = find(hn < 0.5,1,'first') * 2;
    if ~isempty(hw)
        hwTc(nC) = hw;
    end
  end
end

% Calculate FR modulation 
minFR = min(tcAll,[],1);
maxFR = max(tcAll,[],1);
ampl = maxFR - minFR;
modIndex = ampl./maxFR;

% normalize tuning curves
tcAllN = tcAll./max(tcAll,[],1);
 
%% Figures

fsc = [255 175 100]./255;
fsc1 = [255 175 100]./255;
fsc2 = [155 128 0]./255;
hdc = [200 150 255]./255;
hdc1 = [200 150 255]./255;
hdc2 = []./255;
unc = [0.5 0.5 0.5];
exc = 'k';

lw = 1;
fs = 11;

%% Figure 1: Classification criteria and tuning
  
figure (1), clf
set(gcf, 'Color', 'w')
tiledlayout(2,4)
     

nexttile
hold on
    vbl1 = tpk(ixGd);
    h1 = histogram(vbl1);
    h1.BinWidth = 0.02;
    h1.DisplayStyle = 'stairs';
    h1.FaceColor = 'none';
    h1.EdgeColor = [0.5 0.5 0.5];
    h1.LineWidth = 1;
    
    tpk_th = 0.35;
    ax = plot([tpk_th tpk_th],[0 400],'--','Color','k');
    ax.LineWidth = lw;
              
    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    ax.XTick = [0:0.4:1.2];
    xlim([0 1.2]);
    axis(ax, 'square');
    xlabel(['Trough to peak (ms)']);
    ylabel(['No. cells'])
        
nexttile
hold on
       ax = scatter(tpk(ixGd),frate(ixGd));
           ax.MarkerFaceColor = unc;
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 5;
       
       ax = scatter(tpk(ixFs),frate(ixFs));
           ax.MarkerFaceColor = fsc;
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 5;
       
       ax = scatter(tpk(ixEx),frate(ixEx));
           ax.MarkerFaceColor = exc;
           ax.MarkerEdgeColor = 'none';
           ax.SizeData = 5;

        ax = gca;
        ax.XTick = [0:0.4:1.2];
        ax.YTick = [0:20:100];
        ax.TickLength = [0.03 0.025];
        ax.FontSize = fs;
        ax.LineWidth = lw;
        axis(ax, 'square');
        xlabel(['Trough to peak (ms)']);
        ylabel(['Firing rate (Hz)']);
        xlim([0 1.2]);
        ylim([0 100]);
        
        mx = max(tpk(ixGd));
        my = max(frate(ixGd));
        tpk_th = 0.35;
        frate_th = 10;
        
        ax = plot([tpk_th tpk_th],[0 100],'--','Color','k');
        ax.LineWidth = lw;
        ax = plot([0 1.2],[frate_th frate_th],'--','Color','k');
        ax.LineWidth = lw;

nexttile
hold on
    vbl1 = log10(hdi(ixGd));
    vbl2 = log10(hdi(ixHd));
    vbl3 = log10(hdi(ixFs));
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h3 = histogram(vbl3);
    h1.BinWidth = 0.1;
    h2.BinWidth = 0.1;
    h3.BinWidth = 0.1;
    h1.DisplayStyle = 'stairs';
    h1.FaceColor = 'none';
    h1.EdgeColor = [0.5 0.5 0.5];
    h1.LineWidth = 1;
    h2.DisplayStyle = 'stairs';
    h2.FaceColor = 'none';
    h2.EdgeColor = hdc;
    h2.LineWidth = 1;
    h3.DisplayStyle = 'stairs';
    h3.FaceColor = 'none';
    h3.EdgeColor = fsc;
    h3.LineWidth = 1;
    tpk_th = 0.35;
   % ax = plot([tpk_th tpk_th],[0 400],'--','Color','k');
    ax.LineWidth = lw;
              
    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    ax.XTick = [-3:1:1];
    xlim([-3 1]);
    ylim([0 300]);
    axis(ax, 'square');
    xlabel(['log10(HD info)']);
    ylabel(['No. cells'])

nexttile
hold on
    vbl1 = log10(hdi(ixEx));
    vbl2 = log10(hdi_sh(ixEx));
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.1;
    h2.BinWidth = 0.1;
    %h1.Normalization = 'probability';
    %h2.Normalization = 'probability'; 
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.FaceColor = 'none';
    h1.EdgeColor = hdc;
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
    ax.Color = hdc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = plot([thr_ex thr_ex],[0 400],'--');
    ax.LineWidth = lw;
    ax.Color = 'k';
    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    ax.XTick = [-3 -2 -1 0 1];
    xlim([-3 1]);
    axis(ax, 'square');
    xlabel(['log10(HD info)'])
    ylabel(['No. Ex cells'])
    mRef = m1;
    [p,h,stats] = signrank(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';

%2nd row
nexttile
hold on
    vbl1 = log10(hdi(ixFs));
    vbl2 = log10(hdi_sh(ixFs));
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
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
    ax = plot([mRef mRef],[0 100]);
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = gca;
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.XTick = [-4 -2 0];
    xlim([-5 1]);
    ylim([0 120]);
    axis(ax, 'square');
    xlabel(['log10(HD info)'])
    ylabel(['No. FS cells'])
    [p,h,stats] = signrank(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
      
nexttile
hold on
    vbl1 = modIndex(ixFs) * 100;
    m1 = median(vbl1);
    h1 = histogram(vbl1);
    h1.BinWidth = 5;
    h1.DisplayStyle = 'stairs';
    h1.FaceColor = 'none';
    h1.EdgeColor = fsc;
    h1.LineWidth = 1;
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = fsc;
    ax = gca;
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.XTick = [0 50 100];
    xlim([0 100]);
    axis(ax, 'square');
    xlabel(['% Firing rate mod.'])
    ylabel(['No. FS cells'])

 
%% Figure 2: Decoding accuracy 

figure(2),clf
set(gcf, 'Color', 'w');   

subplot(2,4,1)
hold on
        
    vbl1 = totHD;
    vbl2 = medianErrHD;
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 7;
    ax.MarkerFaceColor = hdc;
    ax.MarkerEdgeColor = 'none';
    
    vbl1 = totFS;
    vbl2 = medianErrFS;
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 7;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    
    vbl = mean(medianErrFS_sh);
    ax = plot([0 120],[vbl vbl],':');
    ax.Color = 'k';
    ax.LineWidth = 1;
    

    
    [r,p] = corr(vbl1,vbl2,'type','Pearson');
    
    ax = gca;  
    ax.LineWidth = lw;
    ax.FontSize = fs;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    xlabel({'# cells'})
    ylabel('Decoding error (deg)')
    xlim([0 120])
    %ylim([-0.2 0.2])
    
    
subplot(2,4,2)
hold on
    vbl1 = modIndex(ixFs)*100;
    m1 = median(vbl1);
    h1 = histogram(vbl1);
    h1.BinWidth = 5;
    h1.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h1.FaceColor = 'none';
    h1.EdgeColor = fsc;
    h1.LineWidth = 1;
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = fsc;
    ax = gca;
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    %ax.XTick = [0 0.5 1];
    %xlim([-5 1]);
    axis(ax, 'square');
    xlabel(['% modulation'])
    ylabel(['Prop. FS cells'])
    
    
subplot(2,4,3)
hold on
        
    vbl1 = totFS;
    vbl2 = medianErrFS;
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 20;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    
    vbl = mean(medianErrFS_sh);
    ax = plot([0 120],[vbl vbl],':');
    ax.Color = 'k';
    ax.LineWidth = 1;
    
    vbl = mean(medianErrHD);
    ax = plot([0 120],[vbl vbl],':');
    ax.Color = 'k';
    ax.LineWidth = 1;

    Fit = polyfit(vbl1,vbl2,1);
    xVal = [min(vbl1),max(vbl1)];
    ax = plot(xVal,polyval(Fit,xVal),'-');
    ax.LineWidth = 1;
    ax.Color = 'k';
         
    ax = gca;  
    ax.LineWidth = lw;
    ax.FontSize = fs;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    xlabel({'# cells'})
    ylabel('Decoding error (deg)')
    xlim([0 40])
    %ylim([-0.2 0.2])  
    
    [r,p] = corr(vbl1,vbl2,'type','Spearman');
    ax = title([r;p]);
    ax.FontWeight = 'normal';
    
subplot(2,4,4)
hold on 

    vbl = errHistHD./sum(errHistHD,1);
    m = mean(vbl,2);
    s = sem(vbl,2);
    xVal = bins;  
    boundedline(xVal,m(:,1),s(:,1),'cmap',hdc)
    
    vbl = errHistFS./sum(errHistFS,1);
    m = mean(vbl,2);
    s = sem(vbl,2);
    xVal = bins;   
    boundedline(xVal,m(:,1),s(:,1),'cmap',fsc)
      
    ax = gca; 
    xlim([-180 180]);
    %ylim([0 0.5]);
    ax.TickLength = [0.03 0.025];
    ax.XTick = [-180 0 180];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. bins'])
    xlabel(['Decoding error (deg)'])    
       
