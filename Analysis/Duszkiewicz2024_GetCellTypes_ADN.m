
function PoSubPaper_GetCellTypesADN


    sdSmooth = 3;
   
    %% load data
    load(fullfile('Analysis','MeanFR'));
    load(fullfile('Analysis','HdTuning_moveEp'));
    hdi = hdInfo(:,sdSmooth+1);

    frate = rateS;
    gd = frate > 0.5;
    %% define groups

    % FS cells
    fs = zeros(length(frate),1);

    % Excitatory cells
    ex = gd == 1;
    
    % HD cells
    hd = ex == 1 & hdi >= 0.2;
    
    % Non-HD cells
    nhd = ex == 1 & hdi < 0.2;
    
    %% 
    SaveAnalysis(pwd,'CellTypes',{gd; hd; fs; ex; nhd},{'gd'; 'hd'; 'fs'; 'ex'; 'nhd'});
    
