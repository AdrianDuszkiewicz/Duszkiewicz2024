
%% Load data

% input parameters

 binsize = 0.1; % bin size for decoding in seconds
 nBins = 360; % number of angular bins for tuning curves
 sdSmoothTC = 3; % smoothing widnow for tuning curves
 sdSmoothBD = 2; % smoothing window for spike trains
 velTh = 2; % velocity threshold

% load session files

load(fullfile('Data','SpikeData'),'S');
load(fullfile('Data','BehavEpochs'));
load(fullfile('Data','Angle'));
load(fullfile('Analysis','CellTypes'));
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
%eps{2} = intersect(eps{2},epVel); 
%% Define Groups
totC = length(S);

% compute tuning curves
tcAll = nan(nBins,totC);
for nC = 1:totC 
    [h,b] = HeadDirectionField(S{nC},rAng,baseEp,nBins,sdSmoothTC);
    tcAll(:,nC) = h(1:end-1);    
end

tcAll_sh = nan(nBins,totC);
for nC = 1:totC 
    [h,b] = HeadDirectionField(S{nC},rAng_rev,baseEp,nBins,sdSmoothTC);
    tcAll_sh(:,nC) = h(1:end-1);    
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
er1 = Restrict(decAng_FS,eps{2});
er2 = Restrict(rAng,er1);
errDec_FS = angdiff(Data(er1),Data(er2));
medianErrorFS = rad2deg(nanmedian(abs(errDec_FS)));
decErr_FS = tsd(Range(er2),errDec_FS);


[thetaEst,thetaP,thetaMat] = BayesReconstruction_1D(tcAll_sh(:,ixCells),dQ(:,ixCells),b,binsize); %shuffle
decAng_FS_sh = tsd(Range(Q),thetaEst);
decPval_FS_sh = tsd(Range(Q),thetaP);
decMat_FS_sh = tsd(Range(Q),thetaMat);
er1 = Restrict(decAng_FS_sh,eps{2});
er2 = Restrict(rAng_rev,er1);
errDec_FS_sh = angdiff(Data(er1),Data(er2));
medianErrorFS_sh = rad2deg(nanmedian(abs(errDec_FS_sh)));
decErr_FS_sh = tsd(Range(er2),errDec_FS_sh);

ixCells = find(hd == 1);
[thetaEst,thetaP,thetaMat] = BayesReconstruction_1D(tcAll(:,ixCells),dQ(:,ixCells),b,binsize);
decAng_HD = tsd(Range(Q),thetaEst);
decPval_HD = tsd(Range(Q),thetaP);
decMat_HD = tsd(Range(Q),thetaMat);
er1 = Restrict(decAng_HD,eps{2});
er2 = Restrict(rAng,er1);
errDec_HD = angdiff(Data(er1),Data(er2));
medianErrorHD = rad2deg(nanmedian(abs(errDec_HD)));
decErr_HD = tsd(Range(er2),errDec_HD);

% save decoder data

SaveAnalysis(pwd,'DecodedAngle_FS',{decAng_FS; decPval_FS; decMat_FS; medianErrorFS; decAng_HD; decPval_HD; decMat_HD; medianErrorHD; medianErrorFS_sh; decErr_FS; decErr_FS_sh; decErr_HD},{'decAng_FS';'decPval_FS';'decMat_FS';'medianErrorFS'; 'decAng_HD';'decPval_HD';'decMat_HD';'medianErrorHD'; 'medianErrorFS_sh'; 'decErr_FS'; 'decErr_FS_sh'; 'decErr_HD'});
    



