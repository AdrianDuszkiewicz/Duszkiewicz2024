



load(fullfile('Analysis','SynCrossCorrs'));
binsize = median(diff(bins));

fprintf('Finding synaptic connections \n');
[synLat,synStrZ,synStrR,bounds,ccg_conv] = FindSynapse_Exc(ccg,'synWin',[0 2],'alpha',0.01,'bins',binsize,'convwin',41,'excSigWin',0.4);

    

SaveAnalysis(pwd,'SynPairs_P01',{pairs, ccg, bins, synLat,synStrZ,synStrR,bounds,ccg_conv},{'pairs', 'ccg', 'bins', 'lat', 'strZ', 'strR','bounds','ccg_conv'});


