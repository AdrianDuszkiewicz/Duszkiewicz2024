

%% Load data


binsize = 0.2; % bin size for population vectors, needs to be large
thSmooth = 2; % smoothing factor of population vectors
iso = 0; % set to 1 if need to re-run Isomap

velTh = 2; % velocity threshold (2 = default)
mergeTh = 0; % threshold for merging run epochs (0 = no merging)


[~, foldername, ~] = fileparts(pwd);
load(fullfile('Data','SpikeData'));
load(fullfile('Data','BehavEpochs'));
load(fullfile('Data','CellTypes'));
load(fullfile('Data','Angle'));
load(fullfile('Data','Velocity'));

ep = wake1Ep;
  

% define cell types
frate = Rate(S,ep);
ixHd = find(hd == 1);
ixEx = find(ex == 1);
ixCells = ixHd;
    
%% Calculations for Isomap
   
if iso == 0    
    load(fullfile('Analysis','Isomap_WK'));
else

    % calculate firing rates from population vectors
    Q = MakeQfromS(S,binsize); 
    
    Q = Restrict(Q,ep);
    dQ = Data(Q);
    rQ = Range(Q);
    frate = sum(dQ,1)' ./size(dQ,1);
    frate = frate ./ binsize;

    % process spike trains
    dQ = gaussFilt(dQ,thSmooth,0);
    dQ = sqrt(dQ);
    Q = tsd(rQ,dQ);

    % run Isomap mapping
    dq_iso = dQ(:,ixCells);
    mapping = compute_mapping(dq_iso, 'Isomap',3); % run the Isomap script
    
end


%% Calculations for tuning curves

% get real angle values
rAng = Restrict(ang,Q); % angle values in the middle of Q bins - ang value sampled at much higher freq.
% r = Range(rAng);
% r = r - binsize/2; % move the angle values to the beginning of Q bins, otherwise first half of the bin will have a wrong ang
% d = Data(rAng);
% rAng = tsd(r,d);

% Get isomap angle values
isoAng = deg2rad(atan2d(mapping(:,1),mapping(:,2))); % compute the angle from centre for each point

% Determine whether the angle matches clockwise or counterclockwise:
dr = Data(rAng);
di = isoAng;

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
    isoAng = -isoAng;
end

rQ = Range(Q);
isoAng = tsd(rQ,isoAng);
isoAngRev = tsd(rQ,flipud(Data(isoAng)));

% compute tuning curves for each cell 
totC = length(S); 

nBins = 360;
sdSmooth = 3;
B = 2*pi*(0:1:nBins-1)'/nBins;

tcReal = nan(nBins,totC);
tcWk = nan(nBins,totC);
tcRev = nan(nBins,totC);

h0_real = hist(mod(Data(rAng),2*pi),B); % histograms of angle distribution
h0_wk = hist(mod(Data(isoAng),2*pi),B); % histograms of angle distribution

% Real tuning curves
for nC = 1:totC
    Sr = Restrict(S{nC},ep);
    if ~isempty(Sr)
        angt = Restrict(rAng,Sr);
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
    Sr = Restrict(S{nC},ep);
    if ~isempty(Sr)
        angt = Restrict(isoAng,Sr);
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

% Isomap shuffled tuning curves
for nC = 1:totC
    Sr = Restrict(S{nC},ep);
    if ~isempty(Sr)
        angt = Restrict(isoAngRev,Sr);
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
    tcRev(:,nC) = h;
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
     
% calculate the Isomap radius in Wake 
totS = size(mapping,1);
radius = nan(totS,1);

for nS = 1:totS   
    pts = [0,0; mapping(nS,1),mapping(nS,2)];
    d = pdist(pts,'euclidean');
    radius(nS) = d;
end    

% Save your work
%if iso == 1
    SaveAnalysis(pwd,'Isomap_WK_HdOnly',{Q; mapping; ixCells; tcReal; tcWk; tcRev; B},{'Q'; 'mapping'; 'ixEx'; 'hReal'; 'hWk'; 'hRev'; 'B'});
%end


%return

%%% Figures %%%

wkCol = [0.15 0.3 0.5];
remCol = [0.2 0.5 0.2];

%% Figure 1: Isomap and ring quality 

figure (1), clf
set(gcf, 'Color','w')

subplot(1,2,1)

% this part transforms angle into a color
dTemp = Data(rAng);
dTemp = mat2gray(dTemp)*256; %express angle as numbers 1:256
colIx = floor(dTemp); % express as integers
colIx(find(colIx == 0)) = 1; % get rid of 0
col = hsv; % color map matrix (you can use diferent one)
col = col(colIx,:);

ax = scatter(mapping(:,1),mapping(:,2),80,col,'filled');
%ax.MarkerFaceColor = col;
%ax.MarkerEdgeColor = 'none';
ax.SizeData = 5;

ax = gca; 
ax.LineWidth = 1.5;
ax.FontSize = 20;
axis(ax, 'square');
xlabel('Dimension 1')
ylabel('Dimension 2')
set(gca, 'visible', 'off')

subplot(1,2,2)
hold on
rd = radius;
h = histogram(rd);
h.BinWidth = 0.002;
h.FaceColor = wkCol;
ax.FontSize = 15;
ax.FontWeight = 'normal';
ax = gca; 
ax.LineWidth = 1.5;
ax.FontSize = 20;
xlabel('Distance to center')
ylabel('Counts')

%% Figure 2: Examples of tuning curves

maxC = 16;
b = [B; B(1)];

figure (2), clf
set(gcf, 'Color','w')
tiledlayout(maxC./2,6)
ix = ixHd;
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
    
        
end
    
%% Figure 3: Offset histogram 
 
figure (3), clf
set(gcf, 'Color','w')
    
    ix = ixCells;
    polVec = deg2rad([0:12:360]);
    vbl = deg2rad(offset(ixCells)*(360./nBins));
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