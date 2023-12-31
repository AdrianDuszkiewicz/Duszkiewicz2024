

%% Load data


binsize = 0.2; % bin size for population vectors, needs to be large
thSmooth = 2; % smoothing factor of population vectors
iso = 1; % set to 1 if need to re-run Isomap

load(fullfile('Data','SpikeData'));
load(fullfile('Data','BehavEpochs'));
[~, foldername, ~] = fileparts(pwd);
load(fullfile('Sleep',[foldername '.SleepState.states.mat'] ));
load(fullfile('Data','Angle'));
load(fullfile('Data','CellTypes'));


%% Calculations for Isomap

% get REM-like episodes
rem = SleepState.ints.REMstate;
remDur = sum(rem(:,2) - rem(:,1));
remEp = intervalSet(rem(:,1),rem(:,2));

% bin spike trains
Q = MakeQfromS(S,binsize); 
dQ = Data(Q);
rQ = Range(Q);

% calculate firing rates from population vectors
Qwk = Restrict(Q,wake1Ep);
dQwk = Data(Qwk);
rateWk = sum(dQwk,1)' ./size(dQwk,1);
rateWk = rateWk ./ binsize;

Qrem = Restrict(Q,remEp);
dQrem = Data(Qrem);
rateRem = sum(dQrem,1)' ./size(dQrem,1);
rateRem = rateRem ./ binsize;

% define cell types
goodC = rateWk > 0.5 & rateRem > 0.5;
ixGood = find(goodC == 1);

ixHd = find(hd == 1 & goodC == 1);

ixEx = find(ex == 1 & goodC == 1);

ixFs = find(fs == 1 & goodC == 1);
    
if iso == 0    
    load(fullfile('Analysis','Isomap_REMvsWK_Shuffle'));
else

    % process spike trains
    dQ = gaussFilt(dQ,thSmooth,0);
    dQ = sqrt(dQ);
    Q = tsd(rQ,dQ);

    % calculate REM/Wake pop vectors again and match vector count
    Qrem = Restrict(Q,remEp);
    dQrem = Data(Qrem);
    rQrem = Range(Qrem);
    noVrEp = size(dQrem,1);

    Qwk = Restrict(Q,wake1Ep);
    dQwk = Data(Qwk);
    rQwk = Range(Qwk);
    noVrWk = size(dQwk,1);

    
    rd = rand(noVrWk,1); % pick random bins from Wake
    [~,sIx] = sort(rd);
    dQwk = dQwk(sIx,:);
    rQwk = rQwk(sIx,:);
    dQwk = dQwk(1:noVrEp,:);  
    rQwk = rQwk(1:noVrEp,:);
    
    % Jitter population vectors cell-wise
    [l,totC] = size(dQwk);
    randVecWk = rand(totC,1);
    randVecWk = round(randVecWk * l,0);
    
    for nC = 1:totC
        dQwk(:,nC) = circshift(dQwk(:,nC),randVecWk(nC)); % circshift each cell's spike train by a random amount
    end
    
    [l,totC] = size(dQrem);
    randVecRem = rand(totC,1);
    randVecRem = round(randVecRem * l,0);
    
    for nC = 1:totC
        dQrem(:,nC) = circshift(dQrem(:,nC),randVecRem(nC)); % circshift each cell's spike train by a random amount
    end

    dq = [dQwk; dQrem];
    rq = [rQwk; rQrem];

    idIso = zeros(noVrEp*2,1);
    idIso(1:noVrEp) = 1;
    idIso(noVrEp+1:end) = 2;

    % sort timestamps and id vector
    [rq_rw,sIx] = sort(rq);
    dq_rw = dq(sIx,:);
    Qrw = tsd(rq_rw,dq_rw);

    idIso = idIso(sIx);
    ixWk = find(idIso == 1);
    ixRem = find(idIso == 2);


    % run Isomap mapping
    ixCells = ixHd;
    dq_res = dq_rw(:,ixCells);
    mapping = compute_mapping(dq_res, 'Isomap',3); % run the Isomap script
    
end


%% Calculations for tuning curves


% get real angle values
rQrw = Range(Qrw);
dQrw = Data(Qrw);
tsRem = rQrw(ixRem);

%rQrw = rQrw - binsize./2; %move ts from mid to start of Q bins

tsWk = rQrw(ixWk);
tsWk = ts(tsWk);

rAngWk = Restrict(ang,tsWk); % angle values in the middle of Q bins - ang value sampled at much higher freq.
r = Range(rAngWk);
r = r - binsize/2; % move the angle values to the beginning of Q bins, otherwise first half of the bin will have a wrong ang
d = Data(rAngWk);
rAngWk = tsd(r,d);


% Get isomap angle values
isoAng = deg2rad(atan2d(mapping(:,1),mapping(:,2))); % compute the angle from centre for each point

% Determine whether the angle matches clockwise or counterclockwise:
dr = Data(rAngWk);
di = isoAng(ixWk);

ix = find(isnan(d));
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
    isoAng = -isoAng;
end

% ts for isoAng need to be at the beginning of Q bins because they apply to
% the whole Q bin. Otherwise Restrict will not work properly. 

isoAngWk = tsd(rQrw(ixWk)-binsize./2,isoAng(ixWk));
isoAngRem = tsd(rQrw(ixRem)-binsize./2,isoAng(ixRem));

% create intervalSets for each bin 
int = [rQrw - binsize./2 rQrw + binsize./2];
epRem = int(ixRem,:);
epRem = intervalSet(epRem(:,1),epRem(:,2));
epWk = int(ixWk,:);
epWk = intervalSet(epWk(:,1),epWk(:,2));

% compute tuning curves for each cell 
totC = length(S); 

nBins = 60;
sdSmooth = 1;
B = 2*pi*(0:1:nBins-1)'/nBins;

tcReal = nan(nBins,totC);
tcWk = nan(nBins,totC);
tcRem = nan(nBins,totC);

% Headdirectionfield crashes the script
h0_real = hist(mod(Data(rAngWk),2*pi),B); % histograms of angle distribution
h0_wk = hist(mod(Data(isoAngWk),2*pi),B); % histograms of angle distribution
h0_rem = hist(mod(Data(isoAngRem),2*pi),B); % histograms of angle distribution

% Real tuning curves
for nC = 1:totC
    Sr = Restrict(S{nC},epWk);
    if ~isempty(Sr)
        angt = Restrict(rAngWk,Sr);
        h = hist(mod(Data(angt),2*pi),B);
        if sum(h) == 0
            h = h'; % dumb bug where 0-matrix is vertical
        end
        dt = 1./binsize;
        h = dt*h./h0_real;
        h(isnan(h))=0;
        h = gaussFiltAng(h,sdSmooth,1);
    else
        h = zeros(nBins,1);
    end
    tcReal(:,nC) = h;
end

% Isomap wake tuning curves
for nC = 1:totC
    Sr = Restrict(S{nC},epWk);
    if ~isempty(Sr)
        angt = Restrict(isoAngWk,Sr);
        h = hist(mod(Data(angt),2*pi),B);
        if sum(h) == 0
            h = h';
        end
        dt = 1./binsize;
        h = dt*h./h0_wk;
        h(isnan(h)) = 0;
        h = gaussFiltAng(h,sdSmooth,1);
    else
        h = zeros(nBins,1);
    end
    tcWk(:,nC) = h;
end

% Isomap REM tuning curves
for nC = 1:totC
    Sr = Restrict(S{nC},epRem);
    if ~isempty(Sr)
        angt = Restrict(isoAngRem,Sr);
        h = hist(mod(Data(angt),2*pi),B);
        if sum(h) == 0
            h = h';
        end
        dt = 1./binsize;
        h = dt*h./h0_rem;
        h(isnan(h)) = 0;
        h = gaussFiltAng(h,sdSmooth,1);
    else
        h = zeros(nBins,1);
    end
    tcRem(:,nC) = h;
end

% calculate tuning curve offset between Real and Wake 

offset = nan(totC,1);
ccAll = nan(nBins,totC);
 
for nC = 1:totC        
    cc = TCcrosscorr(tcReal(:,nC),tcWk(:,nC)); % xcorr between two tuning curves
    [~,mx] = max(cc); % max point on the xcorr is the offset
    offset(nC) = mx;
    ccAll(:,nC) = cc;
end
     
% calculate the Isomap radius in Wake and REM

totS = size(mapping,1);
radius = nan(totS,1);

for nS = 1:totS   
    pts = [0,0; mapping(nS,1),mapping(nS,2)];
    d = pdist(pts,'euclidean');
    radius(nS) = d;
end
    

% Save your work
if iso == 1
    SaveAnalysis(pwd,'Isomap_REMvsWK_Shuffle',{Qrw; mapping; ixHd; ixWk; ixRem; tcReal; tcWk; tcRem; B; randVecWk; randVecRem},{'Qrw'; 'mapping'; 'ixHd'; 'ixWk'; 'ixRem'; 'hReal'; 'hWk'; 'hRem'; 'B'; 'randVecWk'; 'randVecRem'});
end

return
%%% Figures %%%

wkCol = [0.15 0.3 0.5];
remCol = [0.2 0.5 0.2];

%% Figure 1: Isomap and ring quality 

figure (1), clf
set(gcf, 'Color','w')

subplot(1,2,1)
ax = scatter(mapping(ixWk,1),mapping(ixWk,2));
ax.MarkerFaceColor = wkCol;
ax.MarkerEdgeColor = 'none';
ax.SizeData = 40;
hold on

ax = scatter(mapping(ixRem,1),mapping(ixRem,2));
ax.MarkerFaceColor = remCol;
ax.MarkerEdgeColor = 'none';
ax.SizeData = 40;
set(gca, 'visible', 'off')

% ax = scatter(0,0);
% ax.MarkerFaceColor = 'r';
% ax.MarkerEdgeColor = 'none';
% ax.SizeData = 80;

ax = gca; 
ax.LineWidth = 1.5;
ax.FontSize = 20;
xlabel('Dimension 1')
ylabel('Dimension 2')

subplot(1,2,2)
hold on
rd = radius;
[h,p1] = kstest(rd(ixWk));
[h,p2] = kstest(rd(ixRem));
h1 = histogram(rd(ixWk));
h1.Normalization = 'probability';
h1.BinWidth = 1;
h1.FaceColor = wkCol;
h1.EdgeColor = 'none';
h2 = histogram(rd(ixRem));
h2.Normalization = 'probability';
h2.BinWidth = 1;
h2.FaceColor = remCol;
h2.EdgeColor = 'none';
% ax = title([p1; p2]);
% ax.FontSize = 15;
% ax.FontWeight = 'normal';
ax = gca; 
ax.LineWidth = 2;
ax.FontSize = 20;
xlabel('Distance to center')
ylabel('Prop. bins')

ax = sgtitle({[foldername]; ['REM duration: ' num2str(remDur) ' seconds']});
ax.FontSize = 20;

%% Figure 2: Examples of tuning curves

maxC = 16;
b = [B; B(1)];

figure (2), clf
set(gcf, 'Color','w')
tiledlayout(maxC./2,6)
ix = ixFs;
for nC = 1:maxC
          
    nexttile
    tc = [tcReal(:,ix(nC)); tcReal(1,ix(nC))];
    ax = polarplot(b,tc);
    ax.LineWidth = 2;
    ax.Color = 'k';
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
    ax.LineWidth = 1;
    title(round(max(tc),0))
    
    nexttile
    tc = [tcWk(:,ix(nC)); tcWk(1,ix(nC))];
    ax = polarplot(b,tc);
    ax.LineWidth = 2;
    ax.Color = wkCol;
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
    ax.LineWidth = 1;
    title(round(max(tc),0))
    
    nexttile
    tc = [tcRem(:,ix(nC)); tcRem(1,ix(nC))];
    ax = polarplot(b,tc);
    ax.LineWidth = 2;
    ax.Color = remCol;
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
    ax.LineWidth = 1;
    title(round(max(tc),0))
        
end
    
%% Figure 3: Offset histogram 
 
figure (3), clf
set(gcf, 'Color','w')

    polVec = deg2rad([0:12:360]);
    vbl = deg2rad(offset(ixHd)*(360./nBins));
    subplot(1,2,1)
    ax = polarhistogram(vbl,polVec);
    ax.FaceColor = 'b';
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
    %ax.ThetaTickLabel = {'0' '90' '180' '270'};
    rticks([]);   
    ax.FontSize = 15;
    ax.LineWidth = 1.5;
    title('HD cell offset')