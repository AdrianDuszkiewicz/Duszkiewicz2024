

clear all

%% Load data

    load(fullfile('Data','SpikeData'),'S');
    load(fullfile('Data','BehavEpochs'),'wakeEp');
    load(fullfile('Data','Angle'));   
    load(fullfile('Data','CellTypes'));   
    led1 = load(fullfile('Data','LED1_TS.csv'));
    baseEp = intervalSet(led1(1,1),led1(2,1)); % we can take the duration of baseline epoch directly from led1 (600+ sec)
    
%% Define cell types

    ixEx = find(ex); % Let's not be arbitrary and take all excitatory cells
    totC = length(ixEx); 
    S = S(ixEx); 
     
%% Compute tuning curves

    tcAll    = nan(360,totC); % always prefill a matrix when possible (for speed)
    nbins = 360; 
    smTh = 3;
    for nC = 1:totC
        [h,b]   =  HeadDirectionField(S{nC},ang,baseEp,nbins,smTh);
        tcAll(:,nC)   = h(1:end-1);
    end
    b = b(1:end-1); % the last value is the same as the first
    
%% Run Bayesian decoder
    %binSize = 0.05; 
    binSize = 0.2;
    sdGauss = 2;
     

    Q          = MakeQfromS(S,binSize);
    Q          = Restrict(Q,wakeEp);
    dQ         = Data(Q); % gives you a vector of firing rates for each time bin
    rQ         = Range(Q);
    dQ         = gaussFilt(dQ,sdGauss,0,0);

    [angBayes,pBayes] = BayesReconstruction_1D(tcAll,dQ,b,binSize);
    angBayes     = mod(angBayes,2*pi);
    angBayes     = tsd(rQ,angBayes);
    pBayes       = tsd(rQ,pBayes);

    angReal         = Restrict(ang,rQ);
    angReal         = tsd(Range(angReal),mod(Data(angReal),2*pi));

    errBayes      = angdiff(Data(angReal), Data(angBayes)); 

    errBayes      = tsd(rQ,errBayes);

%% Save analysis

%SaveAnalysis(pwd,'BayesianDecoding_fine_test',{angReal; angBayes; errBayes; pBayes},{'angRealFine'; 'angBayesFine'; 'errBayesFine'; 'pBayesFine'});
SaveAnalysis(pwd,'BayesianDecoding_coarse_test',{angReal; angBayes; errBayes; pBayes},{'angReal'; 'angBayes'; 'errBayes'; 'pBayes'});



