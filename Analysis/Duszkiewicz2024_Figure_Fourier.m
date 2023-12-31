function Duszkiewicz2024_Figure_Fourier

% This script reporoduces panels from:
% Figure 3, Extended Data Figures 3 and 4

% Dependencies: 
% TStoolbox
% CircStat
% Matlab Toolbox for DImensionality Reduction 

%TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%% Parameters

clear all 

smoothTC = 3; % smoothing of tuning curves (3 used in Duszkiewicz et al., 2024)
noComps = 10; % total number of Fourier components used
saveTC = 0; % 1 to save tuning curves
onlyGoodClusters = 0; % 1 to limit to good clusters

%% Load data 

dataset = List2Cell('dataset_All_FT.list');

% get animal ID
idRec = [];
for ii = 1:length(dataset)
    d = cell2mat(dataset(ii));
    idRec = [idRec; str2num(d(2:5))];
end

tcAll = [];
tcAllN = [];
tcAll_sh = [];

acAll = [];
acAll_sh = [];
whichRec = [];
isHd = [];
isFs = [];
isGd = [];
isEx = [];
cDep = [];
hdi = [];
isAdn = [];
isPos = [];
id = [];
cqSc = [];

for ii = 1:length(dataset)
    
    fprintf(['Uploading recording ', dataset{ii},'\n'])     
    load(fullfile(dataset{ii},'Data','CellTypes'));
    load(fullfile(dataset{ii},'Data','BrainArea'));
    load(fullfile(dataset{ii},'Analysis','HdTuning_moveEp'));
    load(fullfile(dataset{ii},'Analysis','TCautocorrs_xval_moveEp')); % for xval
    load(fullfile(dataset{ii},'Analysis','ClusterQualityScore')); % for xval
    
    isHd = [isHd; hd];
    isFs = [isFs; fs];
    isGd = [isGd; gd];
    isEx = [isEx; ex];
    whichRec = [whichRec; repmat(ii,length(hd),1)]; % session ID   
    tcAll = [tcAll hAll(:,:,smoothTC+1)]; % tuning curves
    tcAll_sh = [tcAll_sh hAll_sh(:,:,smoothTC+1)]; % tuning curves
    hdi = [hdi; hdInfo(:,smoothTC+1)];
    isPos = [isPos; pos];
    isAdn = [isAdn; adn];
    acAll = [acAll ac(:,:,smoothTC+1)]; % TC autocorrs
    acAll_sh = [acAll_sh ac_sh(:,:,smoothTC+1)]; % TC autocorrs
    id = [id; repmat(idRec(ii),length(hd),1)];
    cqSc = [cqSc; cqs];
    
end

totC = length(isHd);
nbins = length(b);

%% Align tuning curves and calculate metrics

tcAllN = nan(360,totC);
tcAllN_sh = nan(360,totC);
[~,muMax] = max(tcAll,[],1); 
[~,muMax_sh] = max(tcAll_sh,[],1);  

for nC = 1:totC    
  if ~isnan(muMax(nC))    
    tcAllN(:,nC) = circshift(tcAll(:,nC),-muMax(nC)+1);
    tcAllN_sh(:,nC) = circshift(tcAll_sh(:,nC),-muMax_sh(nC)+1);
  end
end

tcAllN = tcAllN./max(tcAllN,[],1);
tcAllN_sh = tcAllN_sh./max(tcAllN_sh,[],1);

% Compute tuning curve width

nBins = size(tcAll,1);
tcSide = (tcAllN(1:nBins/2+1,:) + flipud([tcAllN(nBins/2+1:end,:);tcAllN(1,:)]))./2; 
tcWidth = nan(totC,1); 

for nC = 1:totC
    tcw = find(tcSide(:,nC) < 0.5,1); 
    if ~isempty(tcw)
        tcWidth(nC) = tcw*2;
    end
end

% Define multipeaked cells

offFR = nan(totC,1);
rVec = nan(totC,1);

for nC = 1:totC  
    if ~isnan(tcWidth(nC))
        offFR(nC) = mean(tcAllN(tcWidth(nC):end-tcWidth(nC),nC),1);
    end
    rVec(nC) = circ_r(b,tcAll(:,nC),[],1);
end


%% define groups

if onlyGoodClusters == 1
    isHd = isHd & cqSc < 0.2;
    isFs = isFs & cqSc < 0.2;
    isGd = isGd & cqSc < 0.2;
end

% HD cells
isHd =  isHd & offFR < 0.05;
ixHd = find(isHd);

% nonHD cells
isNhd = isGd & ~isHd;
ixNhd = find(isNhd);

% FS cells
ixFs = find(isFs); 

% Excitatory cells
ixEx = find(isEx);

%% Compute FFT
   
t = deg2rad([0:1:nbins-1]);
maxBasis = 180;

% Fourier basis
basis = [];
for n = 1:maxBasis
    basis = [basis;cos(n*t) + 1i*sin(n*t)];
end    

vblFT = tcAll; 
vblFT_sh = tcAll_sh;

% FFT code
ft   = basis * vblFT/length(t); %complex Fourier components
ftN  = ft .* conj(ft);% ./ repmat(var(vblFT)/2,[size(ft,1) 1]);
phFt = angle(ft);
ft = abs(ft);
ftN = abs(ftN);

ft_sh   = basis * vblFT_sh/length(t); %complex Fourier components
ftN_sh  = ft_sh .* conj(ft_sh);% ./ repmat(var(vblFT_sh)/2,[size(ft_sh,1) 1]);
phFt_sh = angle(ft_sh);
ft_sh = abs(ft_sh);
ftN_sh = abs(ftN_sh);
    
% cut FT to N first components and normalize
ftN = ftN(1:noComps,:)./repmat(sum(ftN(1:noComps,:),1),noComps,1);
ftN_sh = ftN_sh(1:noComps,:)./repmat(sum(ftN_sh(1:noComps,:),1),noComps,1);   
    

% get highest component out of first three
[~,highComp] = max(ftN(1:3,:)',[],2);

% compute the euclidean distance from FS mean
ftDistFS = nan(totC,1);
meanFtFS = mean(ftN(:,ixFs),2);
meanFtPos = mean(ftN(:,find(isHd & isPos)),2); 
meanFtAdn = mean(ftN(:,find(isHd & isAdn)),2); 

for nC = 1:totC
    ftd = 0;
    for nComp = 1:noComps
        ftd = ftd + ((ftN(nComp,nC) - meanFtFS(nComp))^2);
    end
    ftDistFS(nC) = sqrt(ftd);
end

% compute the euclidean distance from own mean

ftDist = nan(totC,1); 

ix = ixFs;
for nC = 1:length(ix)
    ftd = 0;
    for nComp = 1:noComps
        ftd = ftd + ((ftN(nComp,ix(nC)) - meanFtFS(nComp))^2);
    end
    ftDist(ix(nC)) = sqrt(ftd);
end

ix = find(isHd & isPos);
for nC = 1:length(ix)
    ftd = 0;
    for nComp = 1:noComps
        ftd = ftd + ((ftN(nComp,ix(nC)) - meanFtPos(nComp))^2);
    end
    ftDist(ix(nC)) = sqrt(ftd);
end

ix = find(isHd & isAdn);
for nC = 1:length(ix)
    ftd = 0;
    for nComp = 1:noComps
        ftd = ftd + ((ftN(nComp,ix(nC)) - meanFtAdn(nComp))^2);
    end
    ftDist(ix(nC)) = sqrt(ftd);
end


% compute Kullback-Leibler Divergence

ftDistFsKL = nan(totC,1);

ftDistFsKL(ixHd) = KLDiv(repmat(meanFtFS',length(ixHd),1),ftN(:,ixHd)');

ftDistKL = nan(totC,1);
ftDistKL(ixFs) = KLDiv(repmat(meanFtFS',length(ixFs),1),ftN(:,ixFs)');
ftDistKL(find(isHd & isPos)) = KLDiv(repmat(meanFtPos',length(find(isHd & isPos)),1),ftN(:,find(isHd & isPos))');
ftDistKL(find(isHd & isAdn)) = KLDiv(repmat(meanFtAdn',length(find(isHd & isAdn)),1),ftN(:,find(isHd & isAdn))');

%% Run all ANOVAs

% Repeated measures ANOVA for FT spectra

    vbl = ftN';
    vbl = [vbl(ixFs,:); vbl(find(isHd & isPos),:); vbl(find(isHd & isAdn),:)];
    grp = [ones(size(ixFs)); ones(size(find(isHd & isPos)))+1; ones(size(find(isHd & isAdn)))+2];

    tbl = array2table([grp vbl]);
    tbl.Properties.VariableNames = {'CellType','c1','c2','c3','c4','c5','c6','c7','c8','c9','c10'}; % add variable names as appropriate
    Meas = table([1:size(vbl,2)]','VariableNames',{'Measurements'});
    rm = fitrm(tbl,'c1-c10 ~CellType','WithinDesign',Meas); % add variable names as appropriate
    anTbl = ranova(rm);
      

 % for shuffle  
    vbl = [ftN(:,ixFs)' ftN_sh(:,ixFs)'];
    grp = ones(size(ixFs));

    tbl = array2table([grp vbl]);
    tbl.Properties.VariableNames = {'CellType','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'}; % add variable names as appropriate
    Meas = table(categorical([ones(1,10) ones(1,10)+1])',categorical([1:10 1:10])','VariableNames',{'a','b'});

    rm = fitrm(tbl,'a1-b10~1','WithinDesign',Meas); % add variable names as appropriate
    anTblSh = ranova(rm,'WithinModel','a*b-1');
    
    %pairwise tests 
    vbl1 = ftN(:,ixFs)';
    vbl2 = ftN_sh(:,ixFs)';
    pvals = nan(nComp,1);
    zstat = nan(nComp,1);
    
    for nC = 1:nComp
         [p,~,stats] = signrank(vbl1(:,nC),vbl2(:,nC)); 
         pvals(nC) = p  * nComp;
         zstat(nC) = stats.zval ./nComp;
    end
       
% Independent ANOVA on Distance from mean FT

    ix1 = ixFs;
    ix2 = find(isHd & isPos);
    ix3 = find(isHd & isAdn);
    vbl = ftDist([ix1; ix2; ix3]);
    grp = [ones(size(ix1)); ones(size(ix2))+1; ones(size(ix3))+2];

    [p,tbl,stats] = anova1(vbl,grp,'off')

    pvals = nan(3,1);
    zstat = nan(3,1);
    
    [p,~,stats] = ranksum(ftDist(ix1),ftDist(ix2)); 
    pvals(1) = p * 3;
    zstat(1) = stats.zval ./3;
    [p,~,stats] = ranksum(ftDist(ix1),ftDist(ix3)); 
    pvals(2) = p * 3;
    zstat(2) = stats.zval ./3;
    [p,~,stats] = ranksum(ftDist(ix2),ftDist(ix3)); 
    pvals(3) = p * 3;
    zstat(3) = stats.zval ./3;
    
    
%%% FIGURES %%%

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];
unc = [0.5 0.5 0.5];
adc = [1 0.4 1];

lw = 1;
fs = 11;


%% Figure 1: Main figure panels

figure(1),clf
set(gcf, 'Color', 'w')

subplot(4,4,1)
hold on
    vbl = tcWidth;
    ix1 = find(isHd & isAdn);
    ix2 = find(isHd & isPos);
    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.BinWidth = 10;
    h2.BinWidth = 10;
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = adc;
    h1.FaceColor = 'none';
    h2.EdgeColor = hdc;
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = adc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = gca; 
    xlim([0 120]);
    ylim([0 0.5]);
    ax.TickLength = [0.03 0.025];
    ax.XTick = [0 60 120];
    ax.YTick = [0 0.25 0.5];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. cells'])
    xlabel(['Tuning curve width (deg)'])
    [p,r,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
 
    
subplot(4,4,2)
    labels = {'','',''};
    vbl = groupcounts(highComp(ixFs));
    ax = pie(vbl);
    ax(1).FaceColor = [1 0 0];
    ax(3).FaceColor = [0 1 0];
    ax(5).FaceColor = [0 0 1];


subplot(4,4,3)
hold on

    vbl = ftN(1:3,ixFs) ./ sum(ftN(1:3,ixFs),1);
    vbl = vbl';
    grp = ones(size(vbl)) + [0 1 2];
   % ax = boxplot(vbl,grp,'BoxStyle','outline','OutlierSize',0.1,'Widths',0.5,'Symbol','o','Colors',[fsc; hdc; adc],'Notch','on');
    ax = boxchart(vbl(:),'GroupByColor',grp(:));
    ax(1).BoxFaceColor = [1 0 0];
    ax(1).MarkerColor = 'none';
    ax(1).WhiskerLineColor = [1 0 0];
    ax(1).BoxWidth = 0.3;
    ax(2).BoxFaceColor = [0 1 0];
    ax(2).MarkerColor = 'none';
    ax(2).WhiskerLineColor = [0 1 0];
    ax(2).BoxWidth = 0.3;
    ax(3).BoxFaceColor = [0 0 1]; 
    ax(3).MarkerColor = 'none';
    ax(3).WhiskerLineColor = [0 0 1];
    ax(3).BoxWidth = 0.3;
    ax = gca;
    axis(ax, 'square');
    ax.XTick = [];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ylim([0 1])
            
% 2nd row
    
subplot(4,4,5)
hold on
    vbl = ftN;
    ix = find(isFs);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:noComps,m(1:noComps),s(1:noComps),'cmap',fsc);
    ix = find(isHd & isPos);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:noComps,m(1:noComps),s(1:noComps),'cmap',hdc);
    ix = find(isHd & isAdn);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:noComps,m(1:noComps),s(1:noComps),'cmap',adc);
    xlabel('Fourier component')
    ylabel('Norm. power')   
    ax = gca;
    ax.XTick = [1:noComps];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    %ylim([10^-4 1]);
    %xlim([0 11]);

      
subplot(4,4,6)
hold on
    vbl = cumsum(ftN);
    ix = find(isFs);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:noComps,m(1:noComps),s(1:noComps),'cmap',fsc);
    ix = find(isHd & isPos);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:noComps,m(1:noComps),s(1:noComps),'cmap',hdc);
    ix = find(isHd & isAdn);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:noComps,m(1:noComps),s(1:noComps),'cmap',adc);
    xlabel('Fourier component')
    ylabel('Cumsum(norm. power)')   
    ax = gca;
    ax.XTick = [1:noComps];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylim([0.2 1]);
    %xlim([0 11]);
     
subplot(4,4,7)
hold on
    vbl = ftDistFsKL;
    ix1 = find(isHd & isPos);
    ix2 = find(isHd & isAdn);
    m1 = median(vbl(ix1));
    m2 = median(vbl(ix2));
    h1 = histogram(vbl(ix1));
    h2 = histogram(vbl(ix2));
    h1.BinWidth = 0.1;
    h2.BinWidth = 0.1;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = hdc;
    h1.FaceColor = 'none';
    h2.EdgeColor = adc;
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = adc;
    ax = gca; 
    ylim([0 0.8])
    xlim([0 1.2]);
    ax.XTick = [0:0.5:1];
    ax.YTick = [0:0.2:0.8];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. HD cells'])
    xlabel(['KL dist. from FS (a.u.)']) 
    [p,r,stats] = ranksum(vbl(ix1),vbl(ix2));
    ax = title(p);
    ax.FontWeight = 'normal';
    
subplot(4,4,8)
hold on
  
    ix1 = ixFs;
    ix2 = find(isHd & isPos);
    ix3 = find(isHd & isAdn);
    vbl = ftDistKL([ix1; ix2; ix3]);
    grp = [ones(size(ix1)); ones(size(ix2))+1; ones(size(ix3))+2];
   % ax = boxplot(vbl,grp,'BoxStyle','outline','OutlierSize',0.1,'Widths',0.5,'Symbol','o','Colors',[fsc; hdc; adc],'Notch','on');
    ax = boxchart(vbl,'GroupByColor',grp);
    ax(1).BoxFaceColor = fsc;
    ax(1).MarkerColor = 'none';
    ax(1).WhiskerLineColor = fsc;
    ax(1).BoxWidth = 0.3;
    ax(2).BoxFaceColor = hdc;
    ax(2).MarkerColor = 'none';
    ax(2).WhiskerLineColor = hdc;
    ax(2).BoxWidth = 0.3;
    ax(3).BoxFaceColor = adc; 
    ax(3).MarkerColor = 'none';
    ax(3).WhiskerLineColor = adc;
    ax(3).BoxWidth = 0.3;
    ax = gca;
    axis(ax, 'square');
    ax.XTick = [];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ylim([0 1.5])
    

% 3rd row    
subplot(4,4,9)
hold on
    ix = find(isHd & isPos);
    vbl1 = mod(phFt(1,:)',pi)*2;
    vbl2 = mod(phFt(2,:)',2*pi);
    ax = scatter(rad2deg(vbl1(ix)),rad2deg(vbl2(ix)));
    ax.SizeData = 5;
    ax.MarkerFaceColor = hdc;
    ax.MarkerEdgeColor = 'none';
    
    ax = gca;
    %box on
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ax.YTick = [0 360];
    ax.XTick = [0 360]; 
    %ax.XTickLabel = {'-180' '0' '\pi'};  
    %ax.YTickLabel = {'-\pi' '0' '\pi'};  
    ylim([0 360])
    xlim([0 360])
    xlabel('1st comp. phase (deg)')
    ylabel('2nd comp. phase (deg)')
    [c,p] = circ_corrcc(vbl1(ix),vbl2(ix));
    ax = title({abs(c);p});
    ax.FontWeight = 'normal';
    
    
subplot(4,4,10)
hold on

    ix = find(isHd & isAdn);
    vbl1 = mod(phFt(1,:)',pi)*2;
    vbl2 = mod(phFt(2,:)',2*pi);
    ax = scatter(rad2deg(vbl1(ix)),rad2deg(vbl2(ix)));
    ax.SizeData = 5;
    ax.MarkerFaceColor = adc;
    ax.MarkerEdgeColor = 'none';
    
    ax = gca;
    %box on
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ax.YTick = [0 360];
    ax.XTick = [0 360]; 
    %ax.XTickLabel = {'-180' '0' '\pi'};  
    %ax.YTickLabel = {'-\pi' '0' '\pi'};  
    ylim([0 360])
    xlim([0 360])
    xlabel('1st comp. phase (deg)')
    ylabel('2nd comp. phase (deg)')
    [c,p] = circ_corrcc(vbl1(ix),vbl2(ix));
    ax = title({abs(c);p});
    ax.FontWeight = 'normal';

           
subplot(4,4,11)
hold on

     ix = ixFs;
    vbl1 = mod(phFt(1,:)',pi)*2;
    vbl2 = mod(phFt(2,:)',2*pi);
    ax = scatter(rad2deg(vbl1(ix)),rad2deg(vbl2(ix)));
    ax.SizeData = 5;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    
    ax = gca;
    %box on
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ax.YTick = [0 360];
    ax.XTick = [0 360]; 
    %ax.XTickLabel = {'-180' '0' '\pi'};  
    %ax.YTickLabel = {'-\pi' '0' '\pi'};  
    ylim([0 360])
    xlim([0 360])
    xlabel('1st comp. phase (deg)')
    ylabel('2nd comp. phase (deg)')
    [c,p] = circ_corrcc(vbl1(ix),vbl2(ix));
    ax = title({abs(c);p});
    ax.FontWeight = 'normal';

subplot(4,4,12)
hold on
    vbl = ftN;
    ix = find(isFs);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:noComps,m(1:noComps),s(1:noComps),'cmap',fsc);
    vbl = ftN_sh;
    ix = find(isFs);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:noComps,m(1:noComps),s(1:noComps),'cmap',[0.5 0.5 0.5]);
    xlabel('Fourier component')
    ylabel('Norm. power')   
    ax = gca;
%      ax.XScale = 'log';
%      ax.YScale = 'log';
    ax.XTick = [1:noComps];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    %ylim([10^-4 1]);
    %xlim([0 11]);
    [p1,~,stats1] = signrank(ftN(1,ix),ftN_sh(1,ix));
    [p2,~,stats2] = signrank(ftN(2,ix),ftN_sh(2,ix));
    [p3,~,stats3] = signrank(ftN(3,ix),ftN_sh(3,ix));
    [p4,~,stats4] = signrank(ftN(4,ix),ftN_sh(4,ix));
       
%4th row

subplot(4,4,13)
hold on
    ix = find(isHd & isPos);
    vbl1 = mod(phFt(1,:)',2/3*pi)*3;
    vbl2 = mod(phFt(3,:)',2*pi);
    ax = scatter(rad2deg(vbl1(ix)),rad2deg(vbl2(ix)));
    ax.SizeData = 5;
    ax.MarkerFaceColor = hdc;
    ax.MarkerEdgeColor = 'none';
    
    ax = gca;
    %box on
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ax.YTick = [0 360];
    ax.XTick = [0 360]; 
    %ax.XTickLabel = {'-180' '0' '\pi'};  
    %ax.YTickLabel = {'-\pi' '0' '\pi'};  
    ylim([0 360])
    xlim([0 360])
    xlabel('1st comp. phase (deg)')
    ylabel('3rd comp. phase (deg)')
    [c,p] = circ_corrcc(vbl1(ix),vbl2(ix));
    ax = title({abs(c);p});
    ax.FontWeight = 'normal';
    
subplot(4,4,14)
hold on
    ix = find(isHd & isAdn);
    vbl1 = mod(phFt(1,:)',2/3*pi)*3;
    vbl2 = mod(phFt(3,:)',2*pi);
    ax = scatter(rad2deg(vbl1(ix)),rad2deg(vbl2(ix)));
    ax.SizeData = 5;
    ax.MarkerFaceColor = adc;
    ax.MarkerEdgeColor = 'none';
    ax = gca;
    %box on
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ax.YTick = [0 360];
    ax.XTick = [0 360]; 
    %ax.XTickLabel = {'-180' '0' '\pi'};  
    %ax.YTickLabel = {'-\pi' '0' '\pi'};  
    ylim([0 360])
    xlim([0 360])
    xlabel('1st comp. phase (deg)')
    ylabel('3rd comp. phase (deg)') 
    [c,p] = circ_corrcc(vbl1(ix),vbl2(ix));
    ax = title({abs(c);p});
    ax.FontWeight = 'normal';
    
subplot(4,4,15)
hold on
    ix = ixFs;
    vbl1 = mod(phFt(1,:)',2/3*pi)*3;
    vbl2 = mod(phFt(3,:)',2*pi);
    ax = scatter(rad2deg(vbl1(ix)),rad2deg(vbl2(ix)));
    ax.SizeData = 5;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    ax = gca;
    %box on
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ax.YTick = [0 360];
    ax.XTick = [0 360]; 
    %ax.XTickLabel = {'-180' '0' '\pi'};  
    %ax.YTickLabel = {'-\pi' '0' '\pi'};  
    ylim([0 360])
    xlim([0 360])
    xlabel('1nd comp. phase (deg)')
    [c,p] = circ_corrcc(vbl1(ix),vbl2(ix));
    ax = title({abs(c);p});
    ax.FontWeight = 'normal';
    

    
%% Figure 2: Extended Data: various

figure(2),clf
set(gcf, 'Color', 'w');


% 1st row
subplot(2,4,1)
hold on
    ix = [ixFs; ixHd]; 
    ax = scatter(sum(ft(:,ix).*conj(ft(:,ix))),var(vblFT(:,ix)));
    ax.MarkerEdgeColor = [0.5 0.5 0.5];
    ax.SizeData = 40;   
    
    plot([0 200],[0 400],'k--')
    ax = gca;
    ax.FontSize = fs;
    ax.LineWidth = lw;
    xlim([0 200])
    ylim([0 400])
    ax.TickLength = [0.03 0.025];
    ax.XTick = [0:100:400];
    ax.YTick = [0:100:400];
    axis square
    ylabel('var(r)')
    xlabel('\Sigma c_n')
    
  
subplot(2,4,2)
hold on

    ids = unique(id);
    totID = length(ids);

    for nID = 1:totID
        ix = find(id == ids(nID) & isHd & isPos);
        if ~isempty(ix)
            ax = plot([1:noComps]',mean(ftN(:,ix),2));
            ax.Color = hdc; 
        end
    end
          
    ax = gca;
    ax.XTick = [1:noComps];
    ax.FontSize = 11;
    ax.LineWidth = 1;
    ax.TickLength = [0.03 0.025];
    xlabel('Fourier component')
    ylabel('Norm. power') 
    axis(ax, 'square');
    ylim([0 0.7]); 
    xlim([1 noComps])
  
subplot(2,4,3)
hold on
  
    ids = unique(id);
    totID = length(ids);

    for nID = 1:totID
        ix = find(id == ids(nID) & isFs & isPos);
        if ~isempty(ix)
            ax = plot([1:noComps]',mean(ftN(:,ix),2));
            ax.Color = fsc; 
        end
    end
          
    ax = gca;
    ax.XTick = [1:noComps];
    ax.FontSize = 11;
    ax.LineWidth = 1;
    ax.TickLength = [0.03 0.025];
    xlabel('Fourier component')
    ylabel('Norm. power') 
    axis(ax, 'square');
    ylim([0 0.7]); 
    xlim([1 noComps])
 
subplot(2,4,4)
hold on

    ids = unique(id);
    totID = length(ids);

    for nID = 1:totID
        ix = find(id == ids(nID) & isHd & isAdn);
        if ~isempty(ix)
            ax = plot([1:noComps]',mean(ftN(:,ix),2));
            ax.Color = adc; 
        end
    end
          
    ax = gca;
    ax.XTick = [1:noComps];
    ax.FontSize = 11;
    ax.LineWidth = 1;
    ax.TickLength = [0.03 0.025];
    xlabel('Fourier component')
    ylabel('Norm. power') 
    axis(ax, 'square');
    ylim([0 0.7]); 
    xlim([1 noComps])
    
 
%2nd row
subplot(2,4,5)
hold on
    hdTh = log10(0.2);
    vbl = log10(hdi);
    ix1 = find(isHd & isAdn);
    ix2 = find(isNhd & isAdn);
    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = adc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    ax = plot([hdTh hdTh],[0 50],'--');
    ax.LineWidth = 1;
    ax.Color = 'k';
    ax = gca; 
    xlim([-3 1]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. cells'])
    xlabel(['log10(HD info)'])
    
   
subplot(2,4,6)
hold on

    vbl = offFR;
    ix1 = find(isHd & isPos);
    ix2 = find(isHd & isAdn);
    vbl1 = vbl(ix1); 
    vbl2 = vbl(ix2); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.BinWidth = 0.025;
    h2.BinWidth = 0.025;
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = hdc;
    h1.FaceColor = 'none';
    h2.EdgeColor = adc;
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = adc;
    ax = gca; 
    xlim([0 0.5]);
    ylim([0 0.6]);
    ax.TickLength = [0.03 0.025];
    ax.XTick = [0 0.2 0.4];
    %ax.YTick = [0 0.25 0.5];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. cells']) 
    [p,r,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    

subplot(2,4,7)
hold on
    vbl = acAll;
    xVal = 1:360;
    ix = find(isFs & highComp == 1);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(xVal,m,s,'cmap',[1 0 0]);
    ix = find(isFs & highComp == 2);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(xVal,m,s,'cmap',[0 1 0]);
    ix = find(isFs & highComp == 3);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(xVal,m,s,'cmap',[0 0 1]);
    ix = find(isHd & isPos);
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(xVal,m,s,'cmap',[0.5 0.5 0.5]);
    
    xlabel('Offset (deg)')
    ylabel('Correlation (r)')   
    ax = gca;
    ax.XTick = [0 60 120 180 240 300 360];
    xtickangle(-45);
    ax.FontSize = fs;
    ax.LineWidth = lw;
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    xlim([0 360]);
   
 %% Figure 3: Extended Data: Autocorrelation

figure(3),clf
set(gcf, 'Color', 'w');
             
subplot(3,3,1)
hold on
    ix = find(isFs & highComp == 1);
    vbl1 = acAll(1,:);
    vbl2 = acAll_sh(1,:);
    m1 = median(vbl1(ix));
    m2 = median(vbl2(ix));
    h1 = histogram(vbl1(ix));
    h2 = histogram(vbl2(ix));
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = [1 0 0];
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [1 0 0];
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
   % ylim([0 0.8])
    xlim([-1 1]);
    %ax.XTick = [-1 0 1];
   % ax.YTick = [0:0.2:0.8];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. FS cells'])
    xlabel(['Correlation (r)']) 
    [p,r,stats] = ranksum(vbl1(ix),vbl2(ix));
    z = stats.zval;
    ax = title([z,p]);
    ax.FontWeight = 'normal';
    
subplot(3,3,2)
hold on
    ix = find(isFs & highComp == 1);
    vbl1 = acAll(181,:);
    vbl2 = acAll_sh(181,:);
    m1 = median(vbl1(ix));
    m2 = median(vbl2(ix));
    h1 = histogram(vbl1(ix));
    h2 = histogram(vbl2(ix));
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = [0 1 0];
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [0 1 0];
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
   % ylim([0 0.8])
    xlim([-1 1]);
    %ax.XTick = [-1 0 1];
   % ax.YTick = [0:0.2:0.8];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. FS cells'])
    xlabel(['Correlation (r)']) 
    [p,r,stats] = ranksum(vbl1(ix),vbl2(ix));
    z = stats.zval;
    ax = title([z,p]);
    ax.FontWeight = 'normal';

subplot(3,3,3)
hold on
    ix = find(isFs & highComp == 1);
    vbl1 = acAll(121,:);
    vbl2 = acAll_sh(121,:);
    m1 = median(vbl1(ix));
    m2 = median(vbl2(ix));
    h1 = histogram(vbl1(ix));
    h2 = histogram(vbl2(ix));
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = [0 0 1];
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [0 0 1];
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
   % ylim([0 0.8])
    xlim([-1 1]);
    %ax.XTick = [-1 0 1];
   % ax.YTick = [0:0.2:0.8];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. FS cells'])
    xlabel(['Correlation (r)']) 
    [p,r,stats] = ranksum(vbl1(ix),vbl2(ix));
    z = stats.zval;
    ax = title([z,p]);
    ax.FontWeight = 'normal';
    
%2nd row
subplot(3,3,4)
hold on
    ix = find(isFs & highComp == 2);
    vbl1 = acAll(1,:);
    vbl2 = acAll_sh(1,:);
    m1 = median(vbl1(ix));
    m2 = median(vbl2(ix));
    h1 = histogram(vbl1(ix));
    h2 = histogram(vbl2(ix));
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = [1 0 0];
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [1 0 0];
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
   % ylim([0 0.8])
    xlim([-1 1]);
    %ax.XTick = [-1 0 1];
   % ax.YTick = [0:0.2:0.8];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. FS cells'])
    xlabel(['Correlation (r)']) 
    [p,r,stats] = ranksum(vbl1(ix),vbl2(ix));
    z = stats.zval;
    ax = title([z,p]);
    ax.FontWeight = 'normal';
    
subplot(3,3,5)
hold on
    ix = find(isFs & highComp == 2);
    vbl1 = acAll(181,:);
    vbl2 = acAll_sh(181,:);
    m1 = median(vbl1(ix));
    m2 = median(vbl2(ix));
    h1 = histogram(vbl1(ix));
    h2 = histogram(vbl2(ix));
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = [0 1 0];
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [0 1 0];
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
   % ylim([0 0.8])
    xlim([-1 1]);
    %ax.XTick = [-1 0 1];
   % ax.YTick = [0:0.2:0.8];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. FS cells'])
    xlabel(['Correlation (r)']) 
    [p,r,stats] = ranksum(vbl1(ix),vbl2(ix));
    z = stats.zval;
    ax = title([z,p]);
    ax.FontWeight = 'normal';

subplot(3,3,6)
hold on
    ix = find(isFs & highComp == 2);
    vbl1 = acAll(121,:);
    vbl2 = acAll_sh(121,:);
    m1 = median(vbl1(ix));
    m2 = median(vbl2(ix));
    h1 = histogram(vbl1(ix));
    h2 = histogram(vbl2(ix));
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = [0 0 1];
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [0 0 1];
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
   % ylim([0 0.8])
    xlim([-1 1]);
    %ax.XTick = [-1 0 1];
   % ax.YTick = [0:0.2:0.8];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. FS cells'])
    xlabel(['Correlation (r)']) 
    [p,r,stats] = ranksum(vbl1(ix),vbl2(ix));
    z = stats.zval;
    ax = title([z,p]);
    ax.FontWeight = 'normal';
    
%3rd row
subplot(3,3,7)
hold on
    ix = find(isFs & highComp == 3);
    vbl1 = acAll(1,:);
    vbl2 = acAll_sh(1,:);
    m1 = median(vbl1(ix));
    m2 = median(vbl2(ix));
    h1 = histogram(vbl1(ix));
    h2 = histogram(vbl2(ix));
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = [1 0 0];
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [1 0 0];
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
   % ylim([0 0.8])
    xlim([-1 1]);
    %ax.XTick = [-1 0 1];
   % ax.YTick = [0:0.2:0.8];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. FS cells'])
    xlabel(['Correlation (r)']) 
    [p,r,stats] = ranksum(vbl1(ix),vbl2(ix));
    z = stats.zval;
    ax = title([z,p]);
    ax.FontWeight = 'normal';
    
subplot(3,3,8)
hold on
    ix = find(isFs & highComp == 3);
    vbl1 = acAll(181,:);
    vbl2 = acAll_sh(181,:);
    m1 = median(vbl1(ix));
    m2 = median(vbl2(ix));
    h1 = histogram(vbl1(ix));
    h2 = histogram(vbl2(ix));
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = [0 1 0];
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [0 1 0];
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
   % ylim([0 0.8])
    xlim([-1 1]);
    %ax.XTick = [-1 0 1];
   % ax.YTick = [0:0.2:0.8];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. FS cells'])
    xlabel(['Correlation (r)']) 
    [p,r,stats] = ranksum(vbl1(ix),vbl2(ix));
    z = stats.zval;
    ax = title([z,p]);
    ax.FontWeight = 'normal';

subplot(3,3,9)
hold on
    ix = find(isFs & highComp == 3);
    vbl1 = acAll(121,:);
    vbl2 = acAll_sh(121,:);
    m1 = median(vbl1(ix));
    m2 = median(vbl2(ix));
    h1 = histogram(vbl1(ix));
    h2 = histogram(vbl2(ix));
    h1.BinWidth = 0.2;
    h2.BinWidth = 0.2;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = [0 0 1];
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = [0 0 1];
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
   % ylim([0 0.8])
    xlim([-1 1]);
    %ax.XTick = [-1 0 1];
   % ax.YTick = [0:0.2:0.8];
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['Prop. FS cells'])
    xlabel(['Correlation (r)']) 
    [p,r,stats] = ranksum(vbl1(ix),vbl2(ix));
    z = stats.zval;
    ax = title([z,p]);
    ax.FontWeight = 'normal';

    
     
 %% Figure 4: Isomap of autocorrelograms (PoSub-FS, PoSub-HD, ADN-HD)

figure(4),clf
set(gcf, 'Color', 'w');

ix1 = ixFs;
ix2 = find(isHd & isPos); 
ix3 = find(isHd & isAdn);
ix = [ix1; ix2; ix3];

ix1 = 1:length(ix1);
ix2 = [1:length(ix2)]+length(ix1);
ix3 = length([ix1 ix2])+[1:length(ix3)];

% Compute autocorr embeddings
vbl = acAll;
mapped_data = compute_mapping(vbl(:,ix)', 'Isomap');

% calculate distance between mid of FS and the rest
midFs = mean(mapped_data(ix1,:),1);
md = mapped_data - repmat(midFs,length(ix),1);
dist = sqrt(sum(md.^2,2));
col = ftN(1:3,:)./repmat(sum(ftN(1:3,:),1),[3 1]);

intPar = 1; 
maxBin = 12;
binSize = 0.5;

bgc = [0.9 0.9 0.9];

subplot(3,3,1)
hold on
    ax = scatter(mapped_data(ix1,1),mapped_data(ix1,2),30,col(1:3,ixFs)','.');
    ax.SizeData = 300;
    ax.LineWidth = 1;
    ax = gca;
    ax = gca;
    axis(ax, 'square');
    set(gca, 'visible', 'off')
    xlim([-10 10])
    ylim([-10 10])  
    
subplot(3,3,2)
hold on
    ax = scatter(mapped_data(ix1,1),mapped_data(ix1,2),30,bgc,'.');
    ax.SizeData = 300;
    ax = gca;
    set(ax, 'visible', 'off')
    axis(ax, 'square');
    xlim([-10 10])
    ylim([-10 10]) 
    
    ax = axes('Position',ax.Position);
    [N,x,y] = hist2d(mapped_data(ix2,1),mapped_data(ix2,2),-maxBin:binSize:maxBin,-maxBin:binSize:maxBin); % for real
    N = N./max(N,[],'all');
    N = interp2(N,intPar);
    N = flipud(N);
    if intPar > 0
       x = -maxBin:(binSize/intPar):maxBin;
    end
    cmap = customcolormap([0 1], [hdc; 1 1 1],256);
    a = zeros(size(N)); 
    a(find(N > 0.1)) = 0.9;
    %a = a + 0.5;
    imagesc(x,x,N,'AlphaData',a);
    ax = gca;
    colormap(ax,cmap);
    set(ax, 'visible', 'off')
    axis(ax, 'square');
    xlim([-10 10])
    ylim([-10 10]) 
    
    ax = axes('Position',ax.Position);
    ax = scatter(midFs(1),midFs(2),'x');
    ax.MarkerFaceColor = 'r';
    ax.MarkerEdgeColor = 'r';
    ax.SizeData = 100;
    ax = gca;
    set(ax, 'visible', 'off')
    axis(ax, 'square');
    xlim([-10 10])
    ylim([-10 10]) 
       
    ax = gca;
    ax.LineWidth = 1;
    axis(ax, 'square');
    set(gca, 'visible', 'off')
    xlim([-10 10])
    ylim([-10 10]) 
    
    pos = ax.Position;
    ax = axes('Position',pos);
    ax = gca;
    set(ax, 'visible', 'off')
    axis(ax, 'square');
    colormap(ax,cmap);
    ax = colorbar;
    ax.Position = [pos(1)+0.23 pos(2)+0.1 0.015 0.15];

subplot(3,3,3)
hold on

    ax = scatter(mapped_data(ix1,1),mapped_data(ix1,2),30,bgc,'.');
    ax.SizeData = 300;
    ax = gca;
    set(ax, 'visible', 'off')
    axis(ax, 'square');
    xlim([-10 10])
    ylim([-10 10]) 
    
    ax = axes('Position',ax.Position);
    [N,x,y] = hist2d(mapped_data(ix3,1),mapped_data(ix3,2),-12:0.5:12,-12:0.5:12); % for real
    N = N./max(N,[],'all');
    N = flipud(N);
    N = interp2(N,intPar);
    if intPar > 0
       x = -maxBin:(binSize/intPar):maxBin;
    end
    cmap = customcolormap([0 1], [adc; 1 1 1],256);
    a = zeros(size(N)); 
    %a = a+ 0.5;
    a(find(N > 0.1)) = 1;
    ax = imagesc(ax,x,y,N,'AlphaData',a);
    ax = gca;
    colormap(ax,cmap);
    set(ax, 'visible', 'off')
    axis(ax, 'square');
    xlim([-10 10])
    ylim([-10 10]) 
    
    ax = axes('Position',ax.Position);
    ax = scatter(midFs(1),midFs(2),'x');
    ax.MarkerFaceColor = 'r';
    ax.MarkerEdgeColor = 'r';
    ax.SizeData = 100;
    ax = gca;
    set(ax, 'visible', 'off')
    axis(ax, 'square');
    xlim([-10 10])
    ylim([-10 10]) 
    
    ax.LineWidth = 1;
    ax = gca;
    ax = gca;
    axis(ax, 'square');
    set(gca, 'visible', 'off')
    xlim([-10 10])
    ylim([-10 10])   
    
    pos = ax.Position;
    ax = axes('Position',pos);
    ax = gca;
    set(ax, 'visible', 'off')
    axis(ax, 'square');
    colormap(ax,cmap);
    ax = colorbar;
    ax.Position = [pos(1)+0.23 pos(2)+0.1 0.015 0.15];
       
        
subplot(3,3,4)
hold on

    vbl1 = dist(ix2); 
    vbl2 = dist(ix3); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 1;
    h2.BinWidth = 1;
    h1.Normalization = 'probability';
    h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = hdc;
    h1.FaceColor = 'none';
    h2.EdgeColor = adc;
    h2.FaceColor = 'none';
    mx1 = find(h1.BinEdges >= m1,1)-1;
    mx1 = h1.Values(mx1);
    mx2 = find(h2.BinEdges >= m2,1)-1;
    mx2 = h2.Values(mx2);
    ax = plot([m1 m1],[0 mx1],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = plot([m2 m2],[0 mx2],':');
    ax.LineWidth = 2;
    ax.Color = adc;
    ax = gca; 
    ylim([0 0.3]);
    ax.TickLength = [0.03 0.025];
    ax.YTick = [0 0.1 0.2 0.3];
    ax.FontSize = 11;
    ax.LineWidth = 1;
    axis(ax, 'square');
    xlabel(['Distance (a.u.)'])
    ylabel(['Probability'])
    [p,r,stats] = ranksum(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
subplot(3,3,5)
hold on
    ax = scatter(mapped_data(ix1,1),mapped_data(ix1,2),30,[0.85 0.85 0.85],'.');
    ax.SizeData = 300;
    
    ax = scatter(mapped_data(ix2,1),mapped_data(ix2,2),30,hdc,'.');
    ax.SizeData = 5;
    
    ax = scatter(midFs(1),midFs(2),'x');
    ax.MarkerFaceColor = 'r';
    ax.MarkerEdgeColor = 'r';
    ax.SizeData = 100;
    
    ax = gca;
    colormap(ax,cmap);
    set(ax, 'visible', 'off')
    axis(ax, 'square');
    xlim([-10 10])
    ylim([-10 10]) 

subplot(3,3,6)
hold on
    ax = scatter(mapped_data(ix1,1),mapped_data(ix1,2),30,[0.85 0.85 0.85],'.');
    ax.SizeData = 300;
    
    ax = scatter(mapped_data(ix3,1),mapped_data(ix3,2),30,adc,'.');
    ax.SizeData = 30;
    
    ax = scatter(midFs(1),midFs(2),'x');
    ax.MarkerFaceColor = 'r';
    ax.MarkerEdgeColor = 'r';
    ax.SizeData = 100;
    
    ax = gca;
    colormap(ax,cmap);
    set(ax, 'visible', 'off')
    axis(ax, 'square');
    xlim([-10 10])
    ylim([-10 10]) 

     
       
%% Figure 5 & 6: Isomap of PoSub-FS cell autocorrs

figure(5),clf
set(gcf, 'Color', 'w');

ix = ixFs; 

% Compute autocorr embeddings
vbl = acAll;
mapped_data = compute_mapping(vbl(:,ix)', 'Isomap');
[N,x,y] = hist2d(mapped_data(:,1),mapped_data(:,2),-24:1:24,-24:1:24); % for real
dist = sqrt(sum(mapped_data.^2,2));
isoAng = atan2(mapped_data(:,2),mapped_data(:,1));
isoAng = wrapToPi(isoAng+pi/2);
ftGain = ftN(1:3,:)./mean(ftN(1:3,ix),2);
col = ftGain ./ repmat(sum(ftGain(1:3,:)),[3 1]);

% compute autocorr embeddings for shuffle 
vbl_sh = acAll_sh;
mapped_data_sh = compute_mapping(vbl_sh(:,ix)', 'Isomap'); 
[Njit,x,y] = hist2d(mapped_data_sh(:,1),mapped_data_sh(:,2),-24:1:24,-24:1:24); % for jittered
dist_sh = sqrt(sum(mapped_data_sh.^2,2));
isoAng_sh = atan2(mapped_data_sh(:,2),mapped_data_sh(:,1));
isoAng_sh = wrapToPi(isoAng_sh+pi/2);
col_sh = ftN_sh(1:3,:)./mean(ftN_sh(1:3,ix),2);
col_sh = col_sh ./ repmat(sum(col_sh(1:3,:)),[3 1]);


ixC = [265 136 134 366 284 300 21 357 160]; % edge

totCC = length(ixC);

subplot(2,2,1)
hold on
    ax = scatter(mapped_data(:,1),mapped_data(:,2),30,(col(1:3,ixFs))','.');
    ax.LineWidth = 1;
    ax.SizeData = 400;
    sd = 150;
    
    for nC = 1:totCC
        ix = ixC(nC);
        ax = scatter(mapped_data(ixC,1),mapped_data(ixC,2),sd,'k');
        ax.LineWidth = 2;
    end
   
    ax = gca;
    axis(ax, 'square');
    set(gca, 'visible', 'off');
    xlim([min(mapped_data(:,1)) max(mapped_data(:,1))])
    ylim([min(mapped_data(:,2)) max(mapped_data(:,2))]) 
    
subplot(2,2,2)
hold on

    vbl1 = rad2deg(isoAng);
    edges = [-180:15:180];
    xVal = edges(1:end-1)+diff(edges)./2;
    
    vbl_m = nan(length(edges)-1,1);
    vbl_s = nan(length(edges)-1,1);
    
    vbl2 = col(1,ixFs);
    for nB = 1:length(edges)-1
        ix2 = find(vbl1 >= edges(nB) & vbl1 < edges (nB+1));
        vbl_m(nB) = mean(vbl2(ix2)');
        vbl_s(nB) = sem(vbl2(ix2)');      
    end    
    ax = boundedline(xVal,vbl_m,vbl_s,'cmap',[1 0 0]);
 
    vbl2 = col(2,ixFs);
    for nB = 1:length(edges)-1
        ix2 = find(vbl1 >= edges(nB) & vbl1 < edges (nB+1));
        vbl_m(nB) = mean(vbl2(ix2)');
        vbl_s(nB) = sem(vbl2(ix2)');      
    end    
    ax = boundedline(xVal,vbl_m,vbl_s,'cmap',[0 1 0]);
    
    vbl2 = col(3,ixFs);
    for nB = 1:length(edges)-1
        ix2 = find(vbl1 >= edges(nB) & vbl1 < edges (nB+1));
        vbl_m(nB) = mean(vbl2(ix2)');
        vbl_s(nB) = sem(vbl2(ix2)');      
    end    
    ax = boundedline(xVal,vbl_m,vbl_s,'cmap',[0 0 1]);
    
    ax = gca;
    xlim([-180 180])
    ylim([0 0.8])
    ax.TickLength = [0.03 0.025];
    ax.XTick = [-180 0 180];
    xlabel('Angle from centre (deg)')
    ylabel('Power of component (a.u.)')
    axis(ax, 'square'); 
    ax.LineWidth = lw;
    ax.FontSize = fs;

subplot(2,2,3)
hold on
    ax = scatter(mapped_data_sh(:,1),mapped_data_sh(:,2),30,(col_sh(1:3,ixFs))','.');
    ax.LineWidth = 1;
    ax.SizeData = 400;   
    ax = gca;
    axis(ax, 'square');
    set(gca, 'visible', 'off')
    xlim([min(mapped_data(:,1)) max(mapped_data(:,1))])
    ylim([min(mapped_data(:,2)) max(mapped_data(:,2))]) 
      
subplot(2,2,4)
hold on

    vbl1 = rad2deg(isoAng_sh);
    edges = [-180:15:180];
    xVal = edges(1:end-1)+diff(edges)./2;
    
    vbl_m = nan(length(edges)-1,1);
    vbl_s = nan(length(edges)-1,1);
    
    vbl2 = col(1,ixFs);
    for nB = 1:length(edges)-1
        ix2 = find(vbl1 >= edges(nB) & vbl1 < edges (nB+1));
        vbl_m(nB) = mean(vbl2(ix2)');
        vbl_s(nB) = sem(vbl2(ix2)');      
    end    
    ax = boundedline(xVal,vbl_m,vbl_s,'cmap',[1 0 0]);
 
    vbl2 = col(2,ixFs);
    for nB = 1:length(edges)-1
        ix2 = find(vbl1 >= edges(nB) & vbl1 < edges (nB+1));
        vbl_m(nB) = mean(vbl2(ix2)');
        vbl_s(nB) = sem(vbl2(ix2)');      
    end    
    ax = boundedline(xVal,vbl_m,vbl_s,'cmap',[0 1 0]);
    
    vbl2 = col(3,ixFs);
    for nB = 1:length(edges)-1
        ix2 = find(vbl1 >= edges(nB) & vbl1 < edges (nB+1));
        vbl_m(nB) = mean(vbl2(ix2)');
        vbl_s(nB) = sem(vbl2(ix2)');      
    end    
    ax = boundedline(xVal,vbl_m,vbl_s,'cmap',[0 0 1]);
    
    ax = gca;
    xlim([-180 180])
    ylim([0 0.8])
    ax.TickLength = [0.03 0.025];
    ax.XTick = [-180 0 180];
    xlabel('Angle from centre (deg)')
    ylabel('Power of component (a.u.)')
    axis(ax, 'square'); 
    ax.LineWidth = lw;
    ax.FontSize = fs;  

figure(6),clf
set(gcf, 'Color', 'w');
  
noCols = totCC;

for nC = 1:totCC
    
    subplot(4,noCols,nC)         

    tc = tcAllN(:,ixFs(ixC(nC)));  
    tc = gaussFiltAng(tc,3,0);

        ax = polarplot(b,tc);
        hold on
            ax.Color = col(:,ixFs(ixC(nC)))';
            ax.LineWidth = 2;
            ax = gca;
            ax.ThetaDir = 'clockwise';
            ax.ThetaZeroLocation = 'top';
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
            title(ixC(nC))
            
    subplot(4,noCols,noCols+nC) 
    hold on
            
        xVal = -179:180;
        ac = acAll(:,ixFs(ixC(nC)));
        ac = circshift(ac,179);
        ac = gaussFiltAng(ac,3,0);
        ax = plot(xVal,ac);
        ax.Color = 'k';
        ax.LineWidth = 1;
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
        ylim([-1 1])
            
    subplot(4,noCols,2*noCols+nC) 
    hold on
            
            xVal = 1:3;
            vbl = col(1:3,ixFs(ixC(nC)));
            fc = [1 0 0; 0 1 0; 0 0 1];
            for nB = 1:3
                ax = bar(xVal(nB),vbl(nB));
                ax.FaceColor = fc(nB,:);
                ax.EdgeColor = 'none';
            end
            ax = gca;
            ax.XTick = [1:3];
            ax.FontSize = fs;
            ax.LineWidth = lw;
            ax.YTick = [0 1];
            ax.TickLength = [0.03 0.025];
            axis(ax, 'square');  
            xlim([0.5 3.5])
            ylim([0 1])
            xlabel('Fourier comp.')
            ylabel('Power (a.u.)')
end             


%% Figure 7 & 8: Isomap of PoSub-FS cell autocorrs (examples from one mouse)

figure(7),clf
set(gcf, 'Color', 'w');

ix = ixFs;
recNo = 25;
ixC = find(isFs & whichRec == recNo);
ixC = find(ismember(ix,ixC));
totCC = length(ixC);

subplot(1,1,1)
hold on
    ax = scatter(mapped_data(:,1),mapped_data(:,2),30,[0.7 0.7 0.7],'.');
    ax.LineWidth = 1;
    ax.SizeData = 400;
    sd = 150;
    dx = 0.1; dy = 0.1; % displacement so the text does not overlay the data points

    for nC = 1:totCC
        ix = ixC(nC);
        ax = scatter(mapped_data(ix,1),mapped_data(ix,2),400,col(:,ixFs(ix))','.');
        ax = scatter(mapped_data(ix,1),mapped_data(ix,2),sd,'k');
        ax.LineWidth = 2;
        text(mapped_data(ix,1)+dx, mapped_data(ix,2)+dy, num2str(ixC(nC)));
    end
   
    ax = gca;
    axis(ax, 'square');
    set(gca, 'visible', 'off');
    xlim([min(mapped_data(:,1)) max(mapped_data(:,1))])
    ylim([min(mapped_data(:,2)) max(mapped_data(:,2))]) 
   
figure(8),clf
set(gcf, 'Color', 'w');
  
noCols = totCC;

for nC = 1:totCC
    
    subplot(4,noCols,nC)         

    tc = tcAll(:,ixFs(ixC(nC)));  
    tc = gaussFiltAng(tc,3,0);

        ax = polarplot(b,tc);
        hold on
            ax.Color = col(:,ixFs(ixC(nC)))';
            ax.LineWidth = 2;
            ax = gca;
            ax.ThetaDir = 'clockwise';
            ax.ThetaZeroLocation = 'top';
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
            title(ixC(nC))
            
    subplot(4,noCols,noCols+nC) 
    hold on
            
        xVal = -179:180;
        ac = acAll(:,ixFs(ixC(nC)));
        ac = circshift(ac,179);
        ac = gaussFiltAng(ac,3,0);
        ax = plot(xVal,ac);
        ax.Color = 'k';
        ax.LineWidth = 1;
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
        ylim([-1 1])
            
    subplot(4,noCols,2*noCols+nC) 
    hold on
            
            xVal = 1:3;
            vbl = col(1:3,ixFs(ixC(nC)));
            fc = [1 0 0; 0 1 0; 0 0 1];
            for nB = 1:3
                ax = bar(xVal(nB),vbl(nB));
                ax.FaceColor = fc(nB,:);
                ax.EdgeColor = 'none';
            end
            ax = gca;
            ax.XTick = [1:3];
            ax.FontSize = fs;
            ax.LineWidth = lw;
            ax.YTick = [0 1];
            ax.TickLength = [0.03 0.025];
            axis(ax, 'square');  
            xlim([0.5 3.5])
            ylim([0 1])
            xlabel('Fourier comp.')
            ylabel('Power (a.u.)')
end             

%% Figure 9: Fourier examples: PoSub

% pick cells
ix =  find(isHd & isPos);
ixC =ix([1017 26,15]);

figure(9),clf
set(gcf, 'Color', 'w');

noCols = 3;
noRows = 4;

for nC = 1:length(ixC)
              
    subplot(noCols,noRows,noRows*(nC-1)+1)         

    tc = tcAll(:,ixC(nC));  
    tc = gaussFiltAng(tc,3,0);

        ax = polarplot(b,tc);
        hold on
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
            
 subplot(noCols,noRows,[noRows*(nC-1)+2 noRows*(nC-1)+3])  
    hold on
            
            tc = rescale(tc); 
            tc = [tc; tc];
            xVal = [rad2deg(b); rad2deg(b)+360];
            ax = plot(xVal,tc);
            ax.Color = [0.5 0.5 0.5];
            ax.LineWidth = 2;
            ax = gca;
            ax.XTick = [0 360 720];
            ax.XTickLabel = {'0/360','0/360','0/360'};
            ax.FontSize = fs;
            ax.LineWidth = lw;
            ax.YTick = [0 1];
            ax.TickLength = [0.03 0.025]; 
            xlim([0 720])
            ylim([0 1])
            xlabel('Heading angle (deg)')
            ylabel('Tuning mplitude (a.u.)')
            ax2 = axes('Position',ax.Position);
            hold on
            % shift and scale sine waves
            sw = nan(nbins,3);
            phShift = round(rad2deg(phFt(1:3,ixC(nC)))./[1 2 3]');
            %phShift = [0 0 0];
            for nS = 1:3
                sw(:,nS) = circshift(sin(b*nS+pi/2),phShift(nS)) * ftN(nS,ixC(nC));
            end
            yVal = [sw; sw];
            % plot sine waves
            ax = plot(ax2,xVal,yVal(:,1));
            ax.Color = 'r';
            ax.LineWidth = lw;
            ax = plot(ax2,xVal,yVal(:,2));
            ax.Color = 'g';
            ax.LineWidth = lw;
            ax = plot(ax2,xVal,yVal(:,3));
            ax.Color = 'b';
            ax.LineWidth = lw;
            ax = gca;
            ax.XColor = 'none';
            ax.Color = 'none';
            ax.YAxisLocation = 'right';
            ax.FontSize = fs;
            ax.LineWidth = lw;
            ax.YTick = [-1 0 1];
            ax.TickLength = [0.03 0.025]; 
            xlim([0 720])
            ylim([-1 1])
            %axis(ax, 'square');  
            
    subplot(noCols,noRows,noRows*(nC-1)+4)  
    hold on
            xVal = 1:10;
            ax = plot(xVal,meanFtFS);
            ax.Color = 'k';
            ax.LineWidth = lw;
            vbl = ftN(:,ixC(nC));
            ax = stairs(xVal-0.5,vbl);
            ax.Color = hdc;
            ax.LineWidth = lw;

            ax = gca;
            ax.XTick = [1 10];
            ax.FontSize = fs;
            ax.LineWidth = lw;
            ax.YTick = [0 1];
            ax.TickLength = [0.03 0.025]; 
            axis(ax, 'square');  
            xlim([0.5 10.5])
            ylim([0 1])
            xlabel('Fourier comp.')
            ylabel('Norm power (a.u.)')
    

end  

%% Figure 10: Fourier examples: ADN-HD cells

% pick cells
ix =  find(isHd & isAdn);
ixC =ix([3 4 5]);

noRows = 4;
noCols = 3;

figure(10),clf
set(gcf, 'Color', 'w');

for nC = 1:length(ixC)
              
    subplot(noCols,noRows,noRows*(nC-1)+1)         

    tc = tcAll(:,ixC(nC));  
    tc = gaussFiltAng(tc,3,0);

        ax = polarplot(b,tc);
        hold on
            ax.Color = adc;
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
            
 subplot(noCols,noRows,[noRows*(nC-1)+2 noRows*(nC-1)+3])  
    hold on
            
            tc = rescale(tc); 
            tc = [tc; tc];
            xVal = [rad2deg(b); rad2deg(b)+360];
            ax = plot(xVal,tc);
            ax.Color = [0.5 0.5 0.5];
            ax.LineWidth = 2;
            ax = gca;
            ax.XTick = [0 360 720];
            ax.XTickLabel = {'0/360','0/360','0/360'};
            ax.FontSize = fs;
            ax.LineWidth = lw;
            ax.YTick = [0 1];
            ax.TickLength = [0.03 0.025]; 
            xlim([0 720])
            ylim([0 1])
            xlabel('Heading angle (deg)')
            ylabel('Tuning mplitude (a.u.)')
            ax2 = axes('Position',ax.Position);
            hold on
            % shift and scale sine waves
            sw = nan(nbins,3);
            phShift = round(rad2deg(phFt(1:3,ixC(nC)))./[1 2 3]');
            %phShift = [0 0 0];
            for nS = 1:3
                sw(:,nS) = circshift(sin(b*nS+pi/2),phShift(nS)) * ftN(nS,ixC(nC));
            end
            yVal = [sw; sw];
            % plot sine waves
            ax = plot(ax2,xVal,yVal(:,1));
            ax.Color = 'r';
            ax.LineWidth = lw;
            ax = plot(ax2,xVal,yVal(:,2));
            ax.Color = 'g';
            ax.LineWidth = lw;
            ax = plot(ax2,xVal,yVal(:,3));
            ax.Color = 'b';
            ax.LineWidth = lw;
            ax = gca;
            ax.XColor = 'none';
            ax.Color = 'none';
            ax.YAxisLocation = 'right';
            ax.FontSize = fs;
            ax.LineWidth = lw;
            ax.YTick = [-1 0 1];
            ax.TickLength = [0.03 0.025]; 
            xlim([0 720])
            ylim([-1 1])
            %axis(ax, 'square');  
            
    subplot(noCols,noRows,noRows*(nC-1)+4)  
    hold on
            xVal = 1:10;
            ax = plot(xVal,meanFtFS);
            ax.Color = 'k';
            ax.LineWidth = lw;
            vbl = ftN(:,ixC(nC));
            ax = stairs(xVal-0.5,vbl);
            ax.Color = adc;
            ax.LineWidth = lw;

            ax = gca;
            ax.XTick = [1 10];
            ax.FontSize = fs;
            ax.LineWidth = lw;
            ax.YTick = [0 1];
            ax.TickLength = [0.03 0.025]; 
            axis(ax, 'square');  
            xlim([0.5 10.5])
            ylim([0 1])
            xlabel('Fourier comp.')
            ylabel('Norm power (a.u.)')
    

          
end

%% Figure 11: Examples FS with FT

% pick cells
ix =  ixFs;
ixC = ix([35 109 39 357]); %209, 219 for mixed

noCols = length(ixC);
noRows = 6;

figure(11),clf
set(gcf, 'Color', 'w');

for nC = 1:length(ixC)
              
    subplot(noCols,noRows,noRows*(nC-1)+1)         

    tc = tcAll(:,ixC(nC));  
    tc = gaussFiltAng(tc,3,0);

        ax = polarplot(b,tc);
        hold on
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
                         
    subplot(noCols,noRows,[noRows*(nC-1)+2 noRows*(nC-1)+3])  
    hold on
            
            tc = rescale(tc); 
            tc = [tc; tc];
            xVal = [rad2deg(b); rad2deg(b)+360];
            ax = plot(xVal,tc);
            ax.Color = [0.5 0.5 0.5];
            ax.LineWidth = 2;
            ax = gca;
            ax.XTick = [0 360 720];
            ax.XTickLabel = {'0/360','0/360','0/360'};
            ax.FontSize = fs;
            ax.LineWidth = lw;
            ax.YTick = [0 1];
            ax.TickLength = [0.03 0.025]; 
            xlim([0 720])
            ylim([0 1])
            xlabel('Heading angle (deg)')
            ylabel('Tuning mplitude (a.u.)')
            ax2 = axes('Position',ax.Position);
            hold on
            % shift and scale sine waves
            sw = nan(nbins,3);
            phShift = round(rad2deg(phFt(1:3,ixC(nC)))./[1 2 3]');
            %phShift = [0 0 0];
            for nS = 1:3
                sw(:,nS) = circshift(sin(b*nS+pi/2),phShift(nS)) * ftN(nS,ixC(nC));
            end
            yVal = [sw; sw];
            % plot sine waves
            ax = plot(ax2,xVal,yVal(:,1));
            ax.Color = 'r';
            ax.LineWidth = lw;
            ax = plot(ax2,xVal,yVal(:,2));
            ax.Color = 'g';
            ax.LineWidth = lw;
            ax = plot(ax2,xVal,yVal(:,3));
            ax.Color = 'b';
            ax.LineWidth = lw;
            ax = gca;
            ax.XColor = 'none';
            ax.Color = 'none';
            ax.YAxisLocation = 'right';
            ax.FontSize = fs;
            ax.LineWidth = lw;
            ax.YTick = [-1 0 1];
            ax.TickLength = [0.03 0.025]; 
            xlim([0 720])
            ylim([-1 1])
            %axis(ax, 'square');  

            
  
    subplot(noCols,noRows,noRows*(nC-1)+4)  
    hold on
            
            xVal = 1:10;
            vbl = ftN(:,ixC(nC));
            ax = plot(xVal,meanFtFS);
            ax.Color = 'k';
            ax.LineWidth = lw;
            ax = stairs(xVal-0.5,vbl);
            ax.Color = fsc;
            ax.LineWidth = lw;
            ax = gca;
            ax.XTick = [1 10];
            ax.FontSize = fs;
            ax.LineWidth = lw;
           % ax.YTick = [0 1];
            ax.TickLength = [0.03 0.025]; 
            axis(ax, 'square');  
            xlim([0.5 10.5])
            ymin = min([round(0.8*min(vbl),1); round(0.8*min(cumsum(meanFtFS)),1)]);
            ylim([0 1])
            ax.YTick = [0 1];
            xlabel('Fourier comp.')
            ylabel('Power (a.u.)')   
            
    subplot(noCols,noRows,noRows*(nC-1)+5)  
    hold on
            
            xVal = 1:10;
            vbl = cumsum(ftN(:,ixC(nC)));
  
            ax = stairs(xVal-0.5,vbl);
            ax.Color = fsc;
            ax.LineWidth = lw;
            ax = plot(xVal,cumsum(meanFtFS));
            ax.Color = 'k';
            ax.LineWidth = lw;
            ax = gca;
            ax.XTick = [1 10];
            ax.FontSize = fs;
            ax.LineWidth = lw;
           % ax.YTick = [0 1];
            ax.TickLength = [0.03 0.025]; 
            axis(ax, 'square');  
            xlim([0.5 10.5])
            ymin = min([round(0.8*min(vbl),1); round(0.8*min(cumsum(meanFtFS)),1)]);
            ylim([0 1])
            ax.YTick = [0 1];
            xlabel('Fourier comp.')
            ylabel('Cumsum power (a.u.)')   
            
              
end

%% Figure 12: Examples multipeak

% pick cells
ix =  find(isHd & isPos);
ixC =ix([77 82 14 155]);

figure(12),clf
set(gcf, 'Color', 'w');
set(gcf,'Position');
tiledlayout(2,4)

for nC = 1:length(ixC)
                 
nexttile

    tc = tcAll(:,ixC(nC)); 
    [mx,maxIx] = max(tc);
    tc = tc./mx;
    tcw = tcWidth(ixC(nC));
    vec = zeros(360,1);
    vec(1:2*tcw+1) = 1;
    vec = circshift(vec,maxIx-tcw);
    vec = deg2rad(find(vec == 0));
    
        ax = polarhistogram(0);
        hold on
            ax.BinEdges = deg2rad(0:360);   
            ax.Data = vec;
            ax.FaceColor = [0.8 0.8 0.8];
            ax.EdgeColor = 'none';

        ax = polarplot(b,tc);
        hold on
            ax.Color = hdc;
            ax.LineWidth = 2;
            ax = gca;
            ax.ThetaDir = 'clockwise';
            ax.ThetaZeroLocation = 'top';
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
            title(round(offFR(ixC(nC)),2));  
            
end  

% pick cells
ix =  find(isHd & isAdn);
ixC =ix([65 91 27 36]);


for nC = 1:length(ixC)
                 
nexttile

    tc = tcAll(:,ixC(nC)); 
    [mx,maxIx] = max(tc);
    tc = tc./mx;
    tcw = tcWidth(ixC(nC));
    vec = zeros(360,1);
    vec(1:2*tcw+1) = 1;
    vec = circshift(vec,maxIx-tcw);
    vec = deg2rad(find(vec == 0));
    
        ax = polarhistogram(0);
        hold on
            ax.BinEdges = deg2rad(0:360);   
            ax.Data = vec;
            ax.FaceColor = [0.8 0.8 0.8];
            ax.EdgeColor = 'none';

        ax = polarplot(b,tc);
        hold on
            ax.Color = adc;
            ax.LineWidth = 2;
            ax = gca;
            ax.ThetaDir = 'clockwise';
            ax.ThetaZeroLocation = 'top';
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
            title(round(offFR(ixC(nC)),2));
end



    
    
