function Duszkiewicz2024_Figure_GainSimulation

% This script reporoduces panels from:
% Extended Data Figure 8

% Dependencies: 


% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%% Parameters and data
clear all 

dataset = List2Cell('dataset_Opto.list');


isHd = [];
isEx = [];
isFs = [];
isGd = [];
group = [];

tcAll_on = [];
tcAll_off = [];

frate_on = [];
frate_off = [];
frate_base = [];
isPos = [];

for ii = 1:length(dataset)

    %Load data

    load(fullfile(dataset{ii},'Analysis','BehavEpochs'));
    load(fullfile(dataset{ii},'Analysis','CellTypes'));
    load(fullfile(dataset{ii},'Analysis','SpikeData'));
    load(fullfile(dataset{ii},'Analysis','BrainArea'));
    load(fullfile(dataset{ii},'Analysis','Angle'));
    load(fullfile(dataset{ii},'Analysis','Groups'));
    optoEp = csvread(fullfile(dataset{ii},'Analysis','Opto_TS.csv'));
    
    offStart = 4;
    offEnd = 0;
    ep = wakeOptoEp;
    angrange = Range(Restrict(ang,ep));
    angrange = intervalSet(angrange(1),angrange(end));
    ang = Restrict(ang,ep);
  
    % calculate epochs

    onEp    = intervalSet(optoEp(:,1), optoEp(:,2));
    offEp   = intervalSet(optoEp(:,1)-offStart,optoEp(:,1)-offEnd); 
    
    % get tuning curves for correlations    
    optoSm = 0;
    nbins = 60;
    for nC = 1:length(S)   
        [h1] = HeadDirectionField(S{nC},ang,onEp,nbins,optoSm);
        [h2] = HeadDirectionField(S{nC},ang,offEp,nbins,optoSm);    
        
        tcAll_on = [tcAll_on h1(1:end-1)];
        tcAll_off = [tcAll_off h2(1:end-1)];   
    end

    frate_on = [frate_on; Rate(S,onEp)];
    frate_off = [frate_off; Rate(S,offEp)];
    frate_base = [frate_base; Rate(S,wakeEp)];
    
     % collect variables
    isHd = [isHd; hd];
    isFs = [isFs; fs];
    isGd = [isGd; gd];
    isEx = [isEx; ex];
    group = [group; g];
    isPos = [isPos; pos];

end

totC = length(isFs);

%% define groups
% Good cells 

isGd = isGd == 1 & frate_on > 0.5 & frate_off > 0.5 & frate_base > 0.5;

% HD cells
isHd = isHd == 1 & isGd == 1;
ixHd = find(isHd);

% FS cells
isFs = isFs == 1 & isGd == 1;
ixFs = find(isFs);

% Ex cells
isEx = isEx == 1 & isGd == 1;
ixEx = find(isEx);


%% Compute additive and multiplicative factors

% Normalize tuning curves
tcAllN_off = tcAll_off./max(tcAll_off,[],1);
tcAllN_on = tcAll_on./max(tcAll_off,[],1);       

% get multiplication factors without normalization 
factorAdd_opto = nan(totC,1);
pAll_opto = nan(totC,2);

for nC = 1:totC    
    pAll_opto(nC,:) = polyfit(tcAll_off(:,nC), tcAll_on(:,nC),1);   % slope
    factorAdd_opto(nC) = polyval(pAll_opto(nC,:),0); % Y intercept
end

factorMult_opto = pAll_opto(:,1);

% recompute addition factors with normalization 
for nC = 1:totC    
    pAll_opto(nC,:) = polyfit(tcAllN_off(:,nC), tcAllN_on(:,nC),1);   % slope
    factorAdd_opto(nC) = polyval(pAll_opto(nC,:),0); % Y intercept
end

% compute tuning curve correlations 
[corrOpto,pvalTC] =corr(tcAllN_on,tcAllN_off);
corrOpto = diag(corrOpto); 
pvalTC = diag(pvalTC);
isGdGain = pvalTC < 0.05;

%% Simulate gain of tuning curves

% gain parameters
allGain_MagMult = [1 1.2 1];
allGain_SDmult = [0 0.33 0.33];

allGain_MagAdd = [0 0 0.2];
allGain_SDadd = [0 0.1 0.1];


% generate simiulated tuning curves
totSims = length(allGain_MagMult);
factorAdd_sim = nan(totC,totSims);
factorMult_sim = nan(totC,totSims);
corrSim = nan(totC,totSims);

for nS = 1:totSims
    
    gainMagMult = allGain_MagMult(nS);
    gainSDmult = allGain_SDmult(nS);

    gainMagAdd = allGain_MagAdd(nS);
    gainSDAdd = allGain_SDadd(nS);

    % multiplicative gain
    multi_on = normrnd(gainMagMult,gainSDmult,1,totC);
    
    % additive gain 
    gainMagAdd = normrnd(gainMagAdd,gainSDAdd,1,totC);
    add_on = gainMagAdd.*max(tcAll_off,[],1);   
    
    % create simulated opto and control group 
    tcAllSim_off = tcAll_off;
    tcAllSim_on = nan(size(tcAllSim_off,1),size(tcAllSim_off,2));
    ix = find(group == 1); 
    tcAllSim_on(:,ix) = (tcAll_off(:,ix)+add_on(:,ix)).*multi_on(:,ix); 
    ix = find(group == 2); 
    tcAllSim_on(:,ix) = tcAll_off(:,ix); 
    
    % noise parameters
    noiseMag = 1.8; 
    noiseSD = noiseMag.*repmat(noiseMag,totC,1);
 
    % add noise
    noise_off = nan(nbins,totC);
    noise_on = nan(nbins,totC);
    for nC = 1:totC
        noise_off(:,nC) = normrnd(0,noiseSD(nC),nbins,1);
        noise_on(:,nC) = normrnd(0,noiseSD(nC),nbins,1);
    end

    tcAllSim_off = tcAllSim_off + noise_off; 
    tcAllSim_on = tcAllSim_on + noise_on;

    % normalize tuning curves
    tcAllSimN_off = tcAllSim_off./max(tcAllSim_off,[],1);
    tcAllSimN_on = tcAllSim_on./max(tcAllSim_off,[],1);

    % get multi factors without normalization 
    pAll_sim = nan(totC,2);
    for nC = 1:totC    
        pAll_sim(nC,:) = polyfit(tcAllSim_off(:,nC), tcAllSim_on(:,nC),1);   % slope
    end
    factorMult_sim(:,nS) = pAll_sim(:,1);

    % recompute add factors with normalization 
    for nC = 1:totC    
        pAll_sim(nC,:) = polyfit(tcAllSimN_off(:,nC), tcAllSimN_on(:,nC),1);   % slope
        factorAdd_sim(nC,nS) = polyval(pAll_sim(nC,:),0); % Y intercept
    end

    % compute correlations for factor analysis 
    [r,~] =corr(tcAllSimN_on,tcAllSimN_off);
    corrSim(:,nS) = diag(r); 
 
end

%%% Figures %%%

fsc = [1 0.5 0];
unc = [0.5 0.5 0.5];

lw = 1;
fs = 11;

       
 %% Figure 1: Additive vs multiplicative factors (simulations) FS vs contol

figure (1), clf
set(gcf, 'Color', 'w');

tiledlayout(4,4);

nexttile
hold on
    vbl = factorAdd_opto;
    ix1 = isFs == 1 & group == 1 & isPos & isGdGain;
    ix2 = isFs == 1 & group == 2 & isPos;

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.1;
    h2.BinWidth = 0.1;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = fsc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
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
    xlim([-0.25 1]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. cells'])
    xlabel(['Additive factor'])
    [p,~,stats] = ranksum(vbl1,vbl2);
    z = stats.zval;
    ax = title([z;p]);
    if p >= 0.05
        ax.FontWeight = 'normal';
    end
    
    
nexttile
hold on
    vbl = factorMult_opto;
    ix1 = isFs == 1 & group == 1 & isPos & isGdGain;
    ix2 = isFs == 1 & group == 2 & isPos;

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = fsc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
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
    xlim([0 3]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. cells'])
    xlabel(['Multiplicative factor'])
    [p,~,stats] = ranksum(vbl1,vbl2);
    z = stats.zval;
    ax = title([z;p]);
    if p >= 0.05
        ax.FontWeight = 'normal';
    end
    
nexttile
hold on
    
    
    ix = isFs == 1 & group == 1 & isPos & isGdGain;
    vbl1 = factorMult_opto(ix);
    vbl2 = factorAdd_opto(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 7;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    ax = plot([0 3], [0 0], ':k');
    ax.LineWidth = 1.5;
    ax = plot([1 1], [-1 1], ':k');
    ax.LineWidth = 1.5;
    
    ix = isFs == 1 & group == 2 & isPos;
    vbl1 = factorMult_opto(ix);
    vbl2 = factorAdd_opto(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 7;
    ax.MarkerFaceColor = [0.5 0.5 0.5];
    ax.MarkerEdgeColor = 'none';
    ax = plot([0 3], [0 0], ':k');
    ax.LineWidth = 1.5;
    ax = plot([1 1], [-1 1], ':k');
    ax.LineWidth = 1.5;
    
    ax = gca;  
    ax.LineWidth = lw;
    ax.FontSize = fs;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    xlabel({'Multiplicative factor'})
    ylabel('Additive factor')
    xlim([0 3])
    ylim([-0.5 1])
    
nexttile
hold on
    vbl = corrOpto;
    ix1 = isFs == 1 & group == 1 & isPos & isGdGain;
    ix2 = isFs == 1 & group == 2 & isPos;

    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.1;
    h2.BinWidth = 0.1;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = fsc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
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
    xlim([-0.25 1]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. cells'])
    [p,~,stats] = ranksum(vbl1,vbl2);
    z = stats.zval;
    ax = title([z;p]);
    if p >= 0.05
        ax.FontWeight = 'normal';
    end
    

for nS = 1:totSims

    nexttile
    hold on
        vbl = factorAdd_sim(:,nS);
        ix1 = find(isFs & isPos & group == 1 & isGdGain);
        ix2 = find(isFs & isPos & group == 2);

        vbl1 = vbl(ix1); 
        vbl2 = vbl(ix2); 
        m1 = median(vbl1);
        m2 = median(vbl2);
        h1 = histogram(vbl1);
        h2 = histogram(vbl2);
        h1.BinWidth = 0.1;
        h2.BinWidth = 0.1;
        h1.Normalization = 'probability';
        h2.Normalization = 'probability';
        h1.DisplayStyle = 'stairs';
        h2.DisplayStyle = 'stairs';
        h1.LineWidth = 1;
        h2.LineWidth = 1;
        h1.EdgeColor = fsc;
        h1.FaceColor = 'none';
        h2.EdgeColor = [0.5 0.5 0.5];
        h2.FaceColor = 'none';
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
        xlim([-0.25 1]);
        ax.TickLength = [0.03 0.025];
        ax.FontSize = fs;
        ax.LineWidth = lw;
        axis(ax, 'square');
        ylabel(['Prop. cells'])
        xlabel(['Additive factor'])
        [p,~,stats] = ranksum(vbl1,vbl2);
        z = stats.zval;
        ax = title([z;p]);
        if p >= 0.05
          ax.FontWeight = 'normal';
        end
    
    nexttile
    hold on
        vbl = factorMult_sim(:,nS);
        ix1 = find(isFs & isPos & group == 1 & isGdGain);
        ix2 = find(isFs & isPos & group == 2);

        vbl1 = vbl(ix1); 
        vbl2 = vbl(ix2); 
        m1 = median(vbl1);
        m2 = median(vbl2);
        h1 = histogram(vbl1);
        h2 = histogram(vbl2);
        h1.BinWidth = 0.2;
        h2.BinWidth = 0.2;
        h1.Normalization = 'probability';
        h2.Normalization = 'probability';
        h1.DisplayStyle = 'stairs';
        h2.DisplayStyle = 'stairs';
        h1.LineWidth = 1;
        h2.LineWidth = 1;
        h1.EdgeColor = fsc;
        h1.FaceColor = 'none';
        h2.EdgeColor = [0.5 0.5 0.5];
        h2.FaceColor = 'none';
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
        xlim([0 3]);
        ax.TickLength = [0.03 0.025];
        ax.FontSize = fs;
        ax.LineWidth = lw;
        axis(ax, 'square');
        ylabel(['Prop. cells'])
        xlabel(['Multiplicative factor'])
        [p,~,stats] = ranksum(vbl1,vbl2);
        z = stats.zval;
        ax = title([z;p]);
        if p >= 0.05
          ax.FontWeight = 'normal';
        end

    nexttile
    hold on

        ix = isFs == 1 & group == 1 & isPos & isGdGain;
        vbl1 = factorMult_sim(ix,nS);
        vbl2 = factorAdd_sim(ix,nS);
        ax = scatter(vbl1,vbl2);
        ax.SizeData = 7;
        ax.MarkerFaceColor = fsc;
        ax.MarkerEdgeColor = 'none';
        ax = plot([0 3], [0 0], ':k');
        ax.LineWidth = 1.5;
        ax = plot([1 1], [-1 1], ':k');
        ax.LineWidth = 1.5;

        ix = isFs == 1 & group == 2 & isPos;
        vbl1 = factorMult_sim(ix);
        vbl2 = factorAdd_sim(ix);
        ax = scatter(vbl1,vbl2);
        ax.SizeData = 7;
        ax.MarkerFaceColor = [0.5 0.5 0.5];
        ax.MarkerEdgeColor = 'none';
        ax = plot([0 3], [0 0], ':k');
        ax.LineWidth = 1.5;
        ax = plot([1 1], [-1 1], ':k');
        ax.LineWidth = 1.5;

        ax = gca;  
        ax.LineWidth = lw;
        ax.FontSize = fs;
        ax.TickLength = [0.03 0.025];
        axis(ax, 'square');
        xlabel({'Multiplicative factor'})
        ylabel('Additive factor')
        xlim([0 3])
        ylim([-0.5 1])
    
    nexttile
    hold on
        vbl = corrSim(:,nS);
        ix1 = find(isFs & isPos & group == 1 & isGdGain);
        ix2 = find(isFs & isPos & group == 2);

        vbl1 = vbl(ix1); 
        vbl2 = vbl(ix2); 
        m1 = median(vbl1);
        m2 = median(vbl2);
        h1 = histogram(vbl1);
        h2 = histogram(vbl2);
        h1.BinWidth = 0.1;
        h2.BinWidth = 0.1;
        h1.Normalization = 'probability';
        h2.Normalization = 'probability';
        h1.DisplayStyle = 'stairs';
        h2.DisplayStyle = 'stairs';
        h1.LineWidth = 1;
        h2.LineWidth = 1;
        h1.EdgeColor = fsc;
        h1.FaceColor = 'none';
        h2.EdgeColor = [0.5 0.5 0.5];
        h2.FaceColor = 'none';
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
        xlim([-0.25 1]);
        ax.TickLength = [0.03 0.025];
        ax.FontSize = fs;
        ax.LineWidth = lw;
        axis(ax, 'square');
        ylabel(['Prop. cells'])
        xlabel(['Correlation (r)'])
        [p,~,stats] = ranksum(vbl1,vbl2);
        z = stats.zval;
        ax = title([z;p]);
        if p >= 0.05
          ax.FontWeight = 'normal';
        end
        
end