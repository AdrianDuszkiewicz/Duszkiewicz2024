function Duszkiewicz2024_Figure_Stability

% This script reporoduces panels from:
% Figure 2, Extended Data Figure 2

% Dependencies: 
% TStoolbox
% CircStat

%TODO: script dependencies, links to external functions

% Copyright (C) 2023 by Adrian Duszkiewicz and Adrien Peyrache
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%% Parameters and variables
clear all

% cell properties
goodSpk = [];
spw = [];
tpk = [];
hpw = [];
frateS = [];
frateT = [];

hdi = [];
hdi1 = [];
hdi2 = [];
hdiT = [];

isHd = [];
isFs = [];
isGd = [];
isTri = [];
isTriRec = [];
whichRec = [];

% for HD-FS symmetry
tcAll1 = [];
tcAll2 = [];
tcAllS = [];
tcAllT = [];
tcAll1_sh = [];      
tcAll2_sh = []; 
tcAllS_sh = [];
tcAllT_sh = [];
 
realign = 1; % 1 if Triangle HD tuning curves should be rotated to compensate for ensamble realignment 
smoothTC = 3;  % smoothing of tuning curves (3 used in Duszkiewicz et al., 2024)
%% Load data

dataset = List2Cell('dataset_All.list');

 
 for ii = 1:length(dataset)

    fprintf(['Uploading recording number ', num2str(ii),'\n'])
    load(fullfile(dataset{ii},'Data','CellTypes'));
    load(fullfile(dataset{ii},'Analysis','MeanFR'));
    load(fullfile(dataset{ii},'Analysis','HdTuning_moveEp'));
    load(fullfile(dataset{ii},'Analysis','HdTuning_xval_moveEp'));

    
    frateS = [frateS; rateS];
    frateT = [frateT; rateT];
    hdi = [hdi; hdInfo(:,smoothTC)];
    hdi1 = [hdi1; hdInfo1(:,smoothTC)];
    hdi2 = [hdi2; hdInfo2(:,smoothTC)];
    hdiT = [hdiT; hdInfoT(:,smoothTC)];
    whichRec = [whichRec; repmat(ii,length(hdInfo),1)];
        
    isHd = [isHd; hd];
    isFs = [isFs; fs];
    isGd = [isGd; gd];
    isTri = [isTri; repmat(tri,length(hd),1)];
    isTriRec = [isTriRec; tri(1)];

    % get tuning curves 
    tcAll1 = [tcAll1 hAll1(:,:,smoothTC)];      
    tcAll2 = [tcAll2 hAll2(:,:,smoothTC)];
    tcAllS = [tcAllS hAll(:,:,smoothTC)]; 
    tcAllT = [tcAllT hAllT(:,:,smoothTC)];
    
    tcAll1_sh = [tcAll1_sh hAll1_sh(:,:,smoothTC)];      
    tcAll2_sh = [tcAll2_sh hAll2_sh(:,:,smoothTC)];  
    tcAllS_sh = [tcAllS_sh hAll_sh(:,:,smoothTC)];  
    tcAllT_sh = [tcAllT_sh hAllT_sh(:,:,smoothTC)];  
      
end

totC = size(isFs,1);
nbins = size(tcAll1,1);

%% Define groups

ixHdT = find(isHd & frateS > 0.5 & frateT > 0.5); % typo in the manuscript
ixFsT = find(isFs & frateS > 10 & frateT > 10);

% Good cells 
ixGd = find(isGd == 1);

% HD cells
ixHd = find(isHd == 1);

% FS cells
ixFs = find(isFs == 1);
      
%% Get correlations 

% Compute all correlations (real) 
ccHH = nan(360,totC);
for nC = 1:totC        
    tc1 = tcAll1(:,nC);
    tc2 = tcAll2(:,nC);
    cc = TCcrosscorr(tc1,tc2); 
    ccHH(:,nC) = cc;
end                  
        
ccST = nan(360,totC);
for nC = 1:totC        
    tc1 = tcAllS(:,nC);
    tc2 = tcAllT(:,nC);
    cc = TCcrosscorr(tc1,tc2); 
    ccST(:,nC) = cc;
end          
  
[~,ccMaxHH] = max(ccHH,[],1);
    ccMaxHH = deg2rad(ccMaxHH');

[~,ccMaxST] = max(ccST,[],1);
ccMaxST = deg2rad(ccMaxST'-1);

% Compute all correlations (reversed) 
ccHH_sh = nan(360,totC);
for nC = 1:totC        
    tc1 = tcAll1_sh(:,nC);
    tc2 = tcAll2_sh(:,nC);
    cc = TCcrosscorr(tc1,tc2); 
    ccHH_sh(:,nC) = cc;
end                  
        
ccST_sh = nan(360,totC);
for nC = 1:totC        
    tc1 = tcAllS_sh(:,nC);
    tc2 = tcAllT_sh(:,nC);
    cc = TCcrosscorr(tc1,tc2); 
    ccST_sh(:,nC) = cc;
end          
  
[~,ccMaxHH_sh] = max(ccHH_sh,[],1);
    ccMaxHH_sh = deg2rad(ccMaxHH_sh');

[~,ccMaxST_sh] = max(ccST_sh,[],1);
ccMaxST_sh = deg2rad(ccMaxST_sh'-1);

% Circshift all triangle tuning curves based on the average HD cell shift
recs = unique(whichRec);
totRecs = length(recs);
ccMaxFs = nan(totRecs,1);
ccMaxHd = nan(totRecs,1);
ccMaxFs_sh = nan(totRecs,1);
ccMaxHd_sh = nan(totRecs,1);
off = [];
off_sh = [];
     
for nR = 1:totRecs   
    if isTriRec(nR) == 1
        ixR = find(whichRec == nR);   

        ixC = ixHd;
        ix = ixC(ismember(ixC,ixR)); 

        cm = CircularMean(ccMaxST(ix)); 
        ccMaxHd(nR) = cm;
        offset = round(rad2deg(cm)); 
        off = [off; offset];

        cm = CircularMean(ccMaxST_sh(ix)); 
        ccMaxHd_sh(nR) = cm;
        offset_sh = round(rad2deg(cm)); 
        off_sh = [off_sh; offset_sh];

        ixC = ixFs;
        ix = ixC(ismember(ixC,ixR));        
        cm = CircularMean(ccMaxST(ix)); 
        ccMaxFs(nR) = cm;  
        cm = CircularMean(ccMaxST_sh(ix)); 
        ccMaxFs_sh(nR) = cm;   

        if realign == 1
            tcAllT(:,ixR) = circshift(tcAllT(:,ixR),offset,1);
            tcAllT_sh(:,ixR) = circshift(tcAllT_sh(:,ixR),offset_sh,1);
        end 
    end
end
    
 % compute all Triangle correlations again
 if realign == 1
     
     ccST = nan(360,totC);
     ccST_sh = nan(360,totC);

    for nC = 1:totC        
        tc1 = tcAllS(:,nC);
        tc2 = tcAllT(:,nC);
        cc = TCcrosscorr(tc1,tc2); 
        ccST(:,nC) = cc;
        
        tc1 = tcAllS_sh(:,nC);
        tc2 = tcAllT_sh(:,nC);
        cc = TCcrosscorr(tc1,tc2); 
        ccST_sh(:,nC) = cc;
        
    end   
    
    [~,ccMaxST] = max(ccST,[],1);
    ccMaxST = deg2rad(ccMaxST'-1);
    
    [~,ccMaxST_sh] = max(ccST_sh,[],1);
    ccMaxST_sh = deg2rad(ccMaxST_sh'-1);
    
 end
      
%Get relevant correlations from xcorrs
  
corrHH = ccHH(1,:);
corrST = ccST(1,:);

zeroHH_rev = ccHH_sh(1,:);
zeroST_rev = ccST_sh(1,:);
%% compute Fourier components 
 
t = deg2rad([0:1:nbins-1]);
maxBasis = 180;
maxComp = 10;
maxComp_score = 3;
    
% Fourier basis
basis = [];
for n = 1:maxBasis
    basis = [basis;cos(n*t) + 1i*sin(n*t)];
end    

vblFT = tcAll1; 
    % FFT code
    ft   = basis * vblFT/length(t); %complex Fourier components
    ftN  = ft .* conj(ft);% ./ repmat(var(vblFT)/2,[size(ft,1) 1]);
    ftN = abs(ftN);
    % cut FT to N first components and normalize
    ftN = ftN(1:maxComp,:)./repmat(sum(ftN(1:maxComp,:),1),maxComp,1); 
    ft1 = ftN;

vblFT = tcAll2; 
    % FFT code
    ft   = basis * vblFT/length(t); %complex Fourier components
    ftN  = ft .* conj(ft);% ./ repmat(var(vblFT)/2,[size(ft,1) 1]);
    ftN = abs(ftN);
    % cut FT to N first components and normalize
    ftN = ftN(1:maxComp,:)./repmat(sum(ftN(1:maxComp,:),1),maxComp,1); 
    ft2 = ftN;
    
vblFT = tcAllS; 
    % FFT code
    ft   = basis * vblFT/length(t); %complex Fourier components
    ftN  = ft .* conj(ft);% ./ repmat(var(vblFT)/2,[size(ft,1) 1]);
    ftN = abs(ftN);
    % cut FT to N first components and normalize
    ftN = ftN(1:maxComp,:)./repmat(sum(ftN(1:maxComp,:),1),maxComp,1); 
    ftS = ftN;
    
vblFT = tcAllT; 
    % FFT code
    ft   = basis * vblFT/length(t); %complex Fourier components
    ftN  = ft .* conj(ft);% ./ repmat(var(vblFT)/2,[size(ft,1) 1]);
    ftN = abs(ftN);
    % cut FT to N first components and normalize
    ftN = ftN(1:maxComp,:)./repmat(sum(ftN(1:maxComp,:),1),maxComp,1); 
    ftT = ftN;

%% Repeated measures ANOVA for FT spectra
% via https://uk.mathworks.com/matlabcentral/answers/124353-does-fitrm-ranova-support-within-subject-models-without-between-subject-factors

 % for two halves   
    vbl = [ft1(:,ixFs)' ft2(:,ixFs)'];
    grp = ones(size(ixFs));

    tbl = array2table([grp vbl]);
    tbl.Properties.VariableNames = {'CellType','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'}; % add variable names as appropriate
    Meas = table(categorical([ones(1,10) ones(1,10)+1])',categorical([1:10 1:10])','VariableNames',{'a','b'});

    rm = fitrm(tbl,'a1-b10~1','WithinDesign',Meas); % add variable names as appropriate
    anTbl12 = ranova(rm,'WithinModel','a*b-1');

 % for two boxes 
    vbl = [ftS(:,ixFsT)' ftT(:,ixFsT)'];
    grp = ones(size(ixFsT));

    tbl = array2table([grp vbl]);
    tbl.Properties.VariableNames = {'CellType','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'}; % add variable names as appropriate
    Meas = table(categorical([ones(1,10) ones(1,10)+1])',categorical([1:10 1:10])','VariableNames',{'a','b'});

    rm = fitrm(tbl,'a1-b10~1','WithinDesign',Meas); % add variable names as appropriate
    anTblST = ranova(rm,'WithinModel','a*b-1');

%% %%% FIGURES %%% %%

fsc = [1 0.5 0];
hdc = [0.6 0.35 1];
unc = [0.5 0.5 0.5];
adc = [1 0.4 1];

fsc_light = [1 0.4 0.4];
fsc_dark = [0.6 0.1 0.1];

lw = 1;
fs = 11;

          
%% Figure 1: Histogram of zero lag correlation

figure (1), clf
set(gcf,'Color','w')
    
ixC1 = 34; % 15 23 34
       
subplot(3,4,1)
    hold on
    ix = ixHdT;
    vbl1 = corrST(ix); 
    vbl2 = zeroST_rev(ix); 
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.05;
    h2.BinWidth = 0.05;
    %h1.Normalization = 'probability';
    %h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = hdc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
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
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    xlim([-1 1]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. FS cells'])
    xlabel(['Correlation (r)'])
    [p,h,stats] = signrank(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    mRef = m1;
                  
subplot(3,4,2)
hold on
    ix = ixFsT;
    vbl1 = corrST(ix); 
    vbl2 = zeroST_rev(ix); 
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
    ax = plot([mRef mRef],[0 80],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = gca; 
    xlim([-1 1]);
    ylim([0 80]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. FS cells'])
    xlabel(['Correlation (r)'])
    [p,h,stats] = signrank(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
    
 subplot(3,4,3)   
    ix = ixFsT;
    tc1 = tcAllS(:,ix(ixC1));  
    ax = polarplot([b],[tc1; tc1(1)]);
    hold on
        ax.Color = fsc_light;
        ax.LineWidth = 2;
    tc2 = tcAllT(:,ix(ixC1));  
    ax = polarplot([b],[tc2; tc2(1)]);
        ax.Color = fsc_dark;
        ax.LineWidth = 2;
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
    
subplot(3,4,4) 
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
    
%2nd row      
subplot(3,4,5)
hold on
    vbl1 = corrHH(ixHd); 
    vbl2 = zeroHH_rev(ixHd);
    m1 = median(vbl1);
    m2 = median(vbl2);
    h1 = histogram(vbl1);
    h2 = histogram(vbl2);
    h1.BinWidth = 0.05;
    h2.BinWidth = 0.05;
    %h1.Normalization = 'probability';
    %h2.Normalization = 'probability';
    h1.DisplayStyle = 'stairs';
    h2.DisplayStyle = 'stairs';
    h1.LineWidth = 1;
    h2.LineWidth = 1;
    h1.EdgeColor = hdc;
    h1.FaceColor = 'none';
    h2.EdgeColor = [0.5 0.5 0.5];
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
    ax.Color = [0.5 0.5 0.5];
    ax = gca; 
    xlim([-1 1]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. HD cells'])
    xlabel(['Correlation (r)'])
    mRef = m1;
    [p,h,stats] = signrank(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
        
subplot(3,4,6)
hold on
    vbl1 = corrHH(ixFs); 
    vbl2 = zeroHH_rev(ixFs); 
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
    ax = plot([mRef mRef],[0 150],':');
    ax.LineWidth = 2;
    ax.Color = hdc;
    ax = gca; 
    xlim([-1 1]);
    ylim([0 150]);
    ax.TickLength = [0.03 0.025];
    ax.FontSize = fs;
    ax.LineWidth = lw;
    axis(ax, 'square');
    ylabel(['No. HD cells'])
    xlabel(['Correlation (r)'])
    [p,h,stats] = signrank(vbl1,vbl2);
    ax = title(p);
    ax.FontWeight = 'normal';
   
subplot(3,4,7)   
    ix = ixFsT;
    tc1 = tcAll1(:,ix(ixC1));  
    ax = polarplot(b,[tc1; tc1(1)]);
    hold on
        ax.Color = fsc_light;
        ax.LineWidth = 2;
    tc2 = tcAll2(:,ix(ixC1));  
    ax = polarplot([b],[tc2; tc2(1)]);
        ax.Color = fsc_dark;
        ax.LineWidth = 2;
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
    
subplot(3,4,8) 
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
 
% 3rd row
subplot(3,4,9) 
hold on
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
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylim([0 0.5])
    %ylim([10^-4 1]);
    %xlim([0 11]);
    
subplot(3,4,10)
hold on
    ix = ixFsT;
    vbl = ftS;
    m = mean(vbl(:,ix),2);
    s = sem(vbl(:,ix),2);
    ax = boundedline(1:maxComp,m(1:maxComp),s(1:maxComp),'cmap',fsc_light);
    vbl = ftT;
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
    ax.TickLength = [0.03 0.025];
    axis(ax, 'square');
    ylim([0 0.5])
    %ylim([10^-4 1]);
    %xlim([0 11]);
    

    
%% Figure 2: Histogram of max representational shift (sanity check)

figure (2), clf
set(gcf,'Color','w')

polVec = deg2rad([0:6:360]);
        
    subplot(2,2,1)
        vbl = deg2rad(ccMaxHH(ixHd));
        ax = polarhistogram(vbl,polVec);
        ax.FaceColor = hdc;
        ax = gca;
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'top';
        ax.GridColor = 'k';
        thetaticks(0:90:360);
        ax.RColorMode = 'manual';
        ax.RColor = 'k';
        ax.ThetaColorMode = 'manual';
        ax.ThetaColor = 'k';
        ax.GridAlpha = 0.3;
        ax.ThetaTickLabel = {'0' '90' '180' '270'};
        rticks([]);   
        ax.FontSize = 15;
        ax.LineWidth = 1.5;
        ax = title('Two halves of square');
        ax.FontWeight = 'normal';
        
    subplot(2,2,2)
        vbl = deg2rad(ccMaxST(ixHd));
        ax = polarhistogram(vbl,polVec);
        ax.FaceColor = hdc;
        ax = gca;
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'top';
        ax.GridColor = 'k';
        thetaticks(0:90:360);
        ax.RColorMode = 'manual';
        ax.RColor = 'k';
        ax.ThetaColorMode = 'manual';
        ax.ThetaColor = 'k';
        ax.GridAlpha = 0.3;
        ax.ThetaTickLabel = {'0' '90' '180' '270'};
        rticks([]);   
        ax.FontSize = 15;
        ax.LineWidth = 1.5;
        ax = title('Square vs Triangle');
        ax.FontWeight = 'normal';
        
    subplot(2,2,3)
        vbl = deg2rad(ccMaxHH(ixFs));
        ax = polarhistogram(vbl,polVec);
        ax.FaceColor = fsc;
        ax = gca;
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'top';
        ax.GridColor = 'k';
        thetaticks(0:90:360);
        ax.RColorMode = 'manual';
        ax.RColor = 'k';
        ax.ThetaColorMode = 'manual';
        ax.ThetaColor = 'k';
        ax.GridAlpha = 0.3;
        ax.ThetaTickLabel = {'0' '90' '180' '270'};
        rticks([]);   
        ax.FontSize = 15;
        ax.LineWidth = 1.5;
        ax = title('Two halves of square');
        ax.FontWeight = 'normal';
          
    subplot(2,2,4)
        vbl = deg2rad(ccMaxST(ixFs));
        ax = polarhistogram(vbl,polVec);
        ax.FaceColor = fsc;
        ax = gca;
        ax.ThetaDir = 'clockwise';
        ax.ThetaZeroLocation = 'top';
        ax.GridColor = 'k';
        thetaticks(0:90:360);
        ax.RColorMode = 'manual';
        ax.RColor = 'k';
        ax.ThetaColorMode = 'manual';
        ax.ThetaColor = 'k';
        ax.GridAlpha = 0.3;
        ax.ThetaTickLabel = {'0' '90' '180' '270'};
        rticks([]);   
        ax.FontSize = 15;
        ax.LineWidth = 1.5;
        ax = title('Square vs Triangle');
        ax.FontWeight = 'normal';
     
