

%% parameters and data
binsize = 0.1; % Bin size in seconds
rng = 10; %Range of CCG in seconds
sdSmooth1 = 1; % smoothing kernel of MUA, ideally std = binsize
sdSmooth2 = 1; % smoothing kernel of spike trains, ideally std ~ 150 ms

load(fullfile('Data','SpikeData'),'S');
load(fullfile('Data','BehavEpochs'));
[~, foldername, ~] = fileparts(pwd);
load(fullfile('Sleep',[foldername '.SleepState.states.mat'] ));

% get REM episodes
rem = SleepState.ints.REMstate;
remDur = sum(rem(:,2) - rem(:,1));
epRem = intervalSet(rem(:,1),rem(:,2));
epWake = wake1Ep;

%% Bin and smoothen spike trains
Q = MakeQfromS(S,binsize);
dQ = Data(Q);
rQ = Range(Q);

dQm = sum(dQ,2); % get MUA for regressor A
dQm = gaussFilt(dQm,sdSmooth1,0);
dQm = zscore(dQm); % regressor needs to be z-scored
Qm = tsd(rQ,dQm);

% Get spikes for reference cell and regressor B
dQ = gaussFilt(dQ,sdSmooth2,0); 
dQs = zscore(dQ); % regressor needs to be z-scored
Qs = tsd(rQ,dQs);

%% Compute GLM with jitter for REM

fprintf('Computing spike train GLMs for REM \n');

totC = length(S);
nbins = (1 ./ binsize) * rng;
pairs = nchoosek(1:totC,2);
totP = size(pairs,1);
beta_rem = nan(nbins*2+1,totP,3);
beta_wake = nan(nbins*2+1,totP,3);
warn_wake = zeros(nbins*2+1,totP); % vector for catching glmfit warnings
warn_rem = zeros(nbins*2+1,totP); % vector for catching glmfit warnings

lagVec = [-rng:binsize:rng]; % ts vector for jitter (if used - not recommended)
stp = rng./binsize;
shiftVec = -stp:1:stp; % step vector for circshift 
totL = length(lagVec); % total number of jitters


Qs_rem = Data(Restrict(Qs,epRem));
Qm_rem = Data(Restrict(Qm,epRem)); 


for nL  = 1:totL
    fprintf(['Computing REM GLMs for time lag ' num2str(nL) ' of ' num2str(totL) '...\n']);    
    dQjit = circshift(dQ,shiftVec(nL)); % add time lag to pop vec ts's 
    Qlag = tsd(rQ,dQjit); % pack back into tsd
    Qlag = Data(Restrict(Qlag,epRem)); % restrict to REM 
       
    for nP = 1:totP
        c = Qlag(:,pairs(nP,1)); % pick cell 1
        X = [Qm_rem Qs_rem(:,pairs(nP,2))]; % pack both regressors
        lastwarn('', ''); % clear last warning
        b = glmfit(X,c,'poisson'); % run GLM
        [warnMsg, ~] = lastwarn(); % catch last warning
        if ~isempty(warnMsg)
            warn_rem(nL,nP) = 1;
        end
        beta_rem(nL,nP,:) = b;
    end
end
        
% compute GLM for WAKE 
fprintf('Computing spike train GLMs for WAKE \n');
Qm_wake = Data(Restrict(Qm,epWake)); 
Qs_wake = Data(Restrict(Qs,epWake));

for nL  = 1:totL
    fprintf(['Computing WAKE GLMs for time lag ' num2str(nL) ' of ' num2str(totL) '...\n']);    
    dQjit = circshift(dQ,shiftVec(nL)); % add time lag to pop vec ts's 
    Qlag = tsd(rQ,dQjit); % pack back into tsd
    Qlag = Data(Restrict(Qlag,epWake)); % restrict to WAKE 
       
    for nP = 1:totP
        c = Qlag(:,pairs(nP,1)); % pick cell 1
        X = [Qm_wake Qs_wake(:,pairs(nP,2))]; % pack both regressors
            lastwarn('', ''); % clear last warning
        b = glmfit(X,c,'poisson'); % run GLM
        [warnMsg, ~] = lastwarn(); % catch last warning
        if ~isempty(warnMsg)
            warn_wake(nL,nP) = 1;
        end
        beta_wake(nL,nP,:) = b;
    end
end       
        
SaveAnalysis(pwd,'GLM_RemWake',{pairs, beta_rem, beta_wake, lagVec, warn_wake, warn_rem},{'pairs', 'beta_rem','beta_wake', 'bins', 'warn_wake', 'warn_rem'});

return

%%% Sanity check figures %%% 

%% Figure 1: Display HD cell GLM ccg's in WAKE and REM 


fsc = [0.9 0.2 0.2];
hdc = [0 0.3 0.8];
unc = [0.5 0.5 0.5];

% load additional data
    load(fullfile('Analysis','GLM_RemWake'));
    load(fullfile('Analysis','TCxcorrs'));
    load(fullfile('Analysis','CellTypes'));
    load(fullfile('Analysis','BehavEpochs'));
    load(fullfile('Analysis','SpikeData'));
    
    [~, foldername, ~] = fileparts(pwd);
    load(fullfile('Sleep',[foldername '.SleepState.states.mat'] ));

    % get REM episodes
    rem = SleepState.ints.REMstate;
    remDur = sum(rem(:,2) - rem(:,1));
    fprintf(['REM duration: ' num2str(remDur) '\n']);
    epRem = intervalSet(rem(:,1),rem(:,2));
    epWake = wakeEp; % CHANGE DEPENDING ON REC!!!
      
% Get cell types
    frate_rem = Rate(S,epRem);
    frate_wake = Rate(S,epWake);
    totC = length(gd);

    isGd = gd == 1 & frate_rem > 1 & frate_wake > 1;
    ixEx = find(ex == 1 & isGd == 1);
    ixHd = find(hd == 1 & isGd == 1);
    ixFs = find(fs == 1 & isGd == 1);

% Get angular differences between cells 
    [~,diffAng] = max(ccgTC,[],1);
    diffAng(diffAng >=180) = diffAng(diffAng >=180) - 360;
    diffAng = abs(diffAng'); 
 
 
% find GLMs with warnings 

[~,colW] = find(warn_wake == 1);
[~,colR] = find(warn_rem == 1);
ixWarn = unique([colW; colR]);

% Plot figure
figure (1), clf
set(gcf,'color','w');

% find HD-HD cell pairs and sort according to ang diff
ix = find(ismember(pairs(:,1),ixHd) & ismember(pairs(:,2),ixHd));
[~,sIx] = sort(diffAng(ix));
ix = ix(sIx); 

cmap = redblue(256);
colormap(cmap); 
sdGauss1 = 4;
sdGauss2 = 4;
    
subplot(2,6,[1 2])
    vbl = beta_wake(:,:,3);
    vbl(:,ixWarn) = nan;
    vbl = vbl(:,ix);
    [row,col] = find(isnan(vbl));
    vbl(:,col) = [];

    
    vbl = gaussFilt(vbl,sdGauss1,sdGauss2);
    vbl = zscore(vbl);
    ax = imagesc(vbl');
    ax = gca;
    ax.XTick = [1 101 201];
    ax.XTickLabel = {'-10', '0','10'};
    xlabel('Time lag (s)')
    ylabel('Cell pair number')
    ax.FontSize = 15;
    ax.LineWidth = 2;
    title('Wake')

subplot(2,6,[3 4])
    vbl = beta_rem(:,:,3);
    vbl(:,ixWarn) = nan;
    vbl = vbl(:,ix);
    [row,col] = find(isnan(vbl));
    vbl(:,col) = [];

    vbl = gaussFilt(vbl,sdGauss1,sdGauss2);
    vbl = zscore(vbl);
    ax = imagesc(vbl');
    ax = gca;
    ax.XTick = [1 101 201];
    ax.XTickLabel = {'-10', '0','10'};
    xlabel('Time lag (s)')
    ylabel('Cell pair number')
    ax.FontSize = 15;
    ax.LineWidth = 2;
    title('REM')
    
subplot(2,6,5)   
   vbl = diffAng(ix);
   %vbl(ix_nan) = [];  
   ax = plot(vbl,1:length(vbl));   
   ax.LineWidth = 2;
   ax.Color = 'k'
   ax = gca;
   ax.YDir = 'reverse';
   ax.XTick = [0 60 120 180]
   ylim([0 length(vbl)])
   xlim([0 180])
   xlabel('Angular offset (deg)')
    ylabel('Cell pair number')
    ax.FontSize = 15;
    ax.LineWidth = 2;
    title('Offset')   

% find FS-FS cell pairs and sort according to ang diff
ix = find(ismember(pairs(:,1),ixFs) & ismember(pairs(:,2),ixFs));
[~,sIx] = sort(diffAng(ix));
ix = ix(sIx); 

cmap = redblue(256);
colormap(cmap); 
sdGauss1 = 1;
    
subplot(2,6,[7 8])
    vbl = beta_wake(:,ix,3);
    %vbl(:,ix_nan) = [];
    vbl = gaussFilt(vbl,sdGauss1);
    vbl = zscore(vbl);
    ax = imagesc(vbl');
    ax = gca;
    ax.XTick = [1 101 201];
    ax.XTickLabel = {'-10', '0','10'};
    xlabel('Time lag (s)')
    ylabel('Cell pair number')
    ax.FontSize = 15;
    ax.LineWidth = 2;
    title('Wake')

subplot(2,6,[9 10])
    vbl = beta_rem(:,ix,3);
    %vbl(:,ix_nan) = [];
    vbl = gaussFilt(vbl,sdGauss1);
    vbl = zscore(vbl);
    ax = imagesc(vbl');
    ax = gca;
    ax.XTick = [1 101 201];
    ax.XTickLabel = {'-10', '0','10'};
    xlabel('Time lag (s)')
    ylabel('Cell pair number')
    ax.FontSize = 15;
    ax.LineWidth = 2;
    title('REM')
    
subplot(2,6,11)   
   vbl = diffAng(ix);
   %vbl(ix_nan) = [];  
   ax = plot(vbl,1:length(vbl));   
   ax.LineWidth = 2;
   ax.Color = 'k'
   ax = gca;
   ax.YDir = 'reverse';
   ax.XTick = [0 60 120 180]
   ylim([0 length(vbl)])
   xlim([0 180])
   xlabel('Angular offset (deg)')
    ylabel('Cell pair number')
    ax.FontSize = 15;
    ax.LineWidth = 2;
    title('Offset')   
    
   
%% Figure 2: Correlations

figure (2), clf
set(gcf,'color','w');
tiledlayout(3,2);

nexttile
hold on
    ix = find(ismember(pairs(:,1),ixHd) & ismember(pairs(:,2),ixHd));
    mid = round(size(beta_wake,1)./2);
    vbl2 = beta_wake(mid,ix,3);
    vbl1 = diffAng(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 10;
    ax.MarkerFaceColor = hdc;
    ax.MarkerEdgeColor = 'none';
    ax = gca;  
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.XTick = [0 60 120 180];
    ylabel('Beta0 (wake)')
    xlabel('Angular offset (deg)')
    [r,p] = corr(vbl1,vbl2','type','Spearman');
    title([r; p]);
    xlim([0 180]);
    ylim([-1 1]);
    
nexttile
hold on
    ix = find(ismember(pairs(:,1),ixHd) & ismember(pairs(:,2),ixHd));
    mid = round(size(beta_rem,1)./2);
    vbl2 = beta_rem(mid,ix,3);
    vbl1 = diffAng(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 10;
    ax.MarkerFaceColor = hdc;
    ax.MarkerEdgeColor = 'none';
    ax = gca;  
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.XTick = [0 60 120 180];
    ylabel('Beta0 (REM)')
    xlabel('Angular offset (deg)')
    [r,p] = corr(vbl1,vbl2','type','Spearman');
    title([r; p]);
    xlim([0 180]);
    ylim([-1 1]);
    
nexttile
hold on
    ix = find(ismember(pairs(:,1),ixFs) & ismember(pairs(:,2),ixFs));
    mid = round(size(beta_wake,1)./2);
    vbl2 = beta_wake(mid,ix,3);
    vbl1 = diffAng(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 15;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    ax = gca;  
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.XTick = [0 60 120 180];
    ylabel('Beta0 (wake)')
    xlabel('Angular offset (deg)')
    [r,p] = corr(vbl1,vbl2','type','Spearman');
    title([r; p]);
    xlim([0 180]);
    ylim([-0.2 0.2]);
    
nexttile
hold on
    ix = find(ismember(pairs(:,1),ixFs) & ismember(pairs(:,2),ixFs));
    mid = round(size(beta_rem,1)./2);
    vbl2 = beta_rem(mid,ix,3);
    vbl1 = diffAng(ix);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 15;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    ax = gca;  
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ax.XTick = [0 60 120 180];
    ylabel('Beta0 (REM)')
    xlabel('Angular offset (deg)')
    [r,p] = corr(vbl1,vbl2','type','Spearman');
    title([r; p]);
    xlim([0 180]);
    ylim([-0.2 0.2]);
    
nexttile
hold on
    ix = find(ismember(pairs(:,1),ixHd) & ismember(pairs(:,2),ixHd));
    mid = round(size(beta_wake,1)./2);
    vbl1 = beta_wake(mid,ix,3);
    vbl2 = beta_rem(mid,ix,3);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 10;
    ax.MarkerFaceColor = hdc;
    ax.MarkerEdgeColor = 'none';
    ax = plot([-1 1], [-1 1],'--k');
    ax.LineWidth = 2;
    ax = gca;  
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ylabel('Beta0 (REM)')
    xlabel('Beta0 (wake)')
    [r,p] = corr(vbl1',vbl2','type','Spearman');
    title([r; p]);
    xlim([-1 1]);
    ylim([-1 1]);
    
nexttile
hold on
    ix = find(ismember(pairs(:,1),ixFs) & ismember(pairs(:,2),ixFs));
    mid = round(size(beta_wake,1)./2);
    vbl1 = beta_wake(mid,ix,3);
    vbl2 = beta_rem(mid,ix,3);
    ax = scatter(vbl1,vbl2);
    ax.SizeData = 10;
    ax.MarkerFaceColor = fsc;
    ax.MarkerEdgeColor = 'none';
    ax = plot([-0.2 0.2], [-0.2 0.2],'--k');
    ax.LineWidth = 2;
    ax = gca;  
    ax.LineWidth = 2;
    ax.FontSize = 15;
    ylabel('Beta0 (REM)')
    xlabel('Beta0 (wake)')
    [r,p] = corr(vbl1',vbl2','type','Spearman');
    title([r; p]);
