function [spkWidth,tr2pk,halfPkW,maxIx,wave] = Make_WaveFormFeatures(waveforms,varargin)

%  Extract spike width, trough to peak and peak half-width from waveform(s)
%  The program assumes spikes are negative. Exclude positive spikes for
%  now.
%
%  USAGE
%
%    Make_WaveFormFeatures(waveforms,<Fq>)
%
%    INPUT
%    waveforms      a vector or matrix of average waveforms (if recorded on a
%                   polytrode). Matrix should be Channels X Samples
%                   OR a cell array of matrices for different neurons
%    Fq (optional)  sampling frequency (default 20,000Hz)
%
%    OUTPUT
%    spkWidth       total spike width (inverse of peak frequency in a
%                   wavelet transform)
%    tr2pk          trough to peak
%    maxIx          channel index of max spike amplitude
% Copyright (C) 2017 Adrien Peyrache

%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

Fq = 20000;

if ~isempty(varargin)
   Fq = varargin{1};
end

spkWidth = [];
tr2pk = [];
maxIx = [];
wave = [];
halfPkW  = [];

if isa(waveforms,'cell')%call itself for each item of the cell array
    for c=1:length(waveforms)
        [spkW,t2p,halfW,mxI,wv]   = Make_WaveFormFeatures(waveforms{c},Fq);
        spkWidth            = [spkWidth;spkW];
        tr2pk               = [tr2pk;t2p];
        maxIx               = [maxIx;mxI];
        halfPkW             = [halfPkW;halfW];
        
    end
    
else
        %if more than one channel, select the one with highest amplitude
        %Restrict to the 1sst 32 samples as sometimes late samples
        %contribute to the power

        if min(size(waveforms))>1
            l               = min(32,size(waveforms,2));
            [~,mxIx]        = max(sum(waveforms(:,1:l).^2,2));
            maxIx(end+1)    = mxIx;
            w               = waveforms(mxIx,:);
        else
            w               = waveforms;
        end

        wu = w-mean(w);
        wu = resample(wu,10,1);
        Fq = Fq*10;
        
        nSamples = length(wu);
        t = [0:1/(Fq):(nSamples-1)/(Fq)]*1000;

        %Positive or negative spike?
        baseLine = mean(wu(1:7));
        [~,absPk] = max(abs(wu));

        if wu(absPk) > baseLine
            fprintf('Positive spike, skipping\n')
            spkWidth    = NaN;
            tr2pk       = NaN;
            halfPkW     = NaN;
            wave        = w;
        else
            
            % Trough is supposed to be sample #17, but we never know...
            [~,minPos] = min(wu); %and should be the same as absPk...

            % Where is the following peak? Restrict to the 1st 240 samples
            l = min(240,length(wu)-minPos);
            [maxVal,maxPos] = max(wu(minPos+1:minPos+l));

            maxPos = maxPos+minPos;
            p2v = t(maxPos)-t(minPos);

            [waveLET,f] = cwt(wu,Fq,'VoicesPerOctave',48);
            waveLET = waveLET(f>500 & f<3000,:);
            f = f(f>500 & f<3000);

            %Where is the max power?
            [maxPow,maxF] = max(waveLET);

            %Which frequency does it correspond to?
            [~,fIx] = max(maxPow);

            maxF = maxF(fIx);

            %Spike width is the inverse of the peak frequency
            spkW = 1000/f(maxF);

            %Half Peak Width
            baseLine = mean(wu(end-10:end));
            baseLine = (maxVal-baseLine)/2;
 
            ix1 = find(wu(minPos:maxPos)>=baseLine);

            if isempty(ix1)
                halfW = NaN;
                figure(1),clf
                plot(wu)
                keyboard
            else

                ix1 = minPos+ix1(1);

                if maxPos ~= length(wu)
                    ix2 = find(wu(maxPos:end)>=baseLine);
                    ix2 = ix2(end)+maxPos;
                else
                    ix2 = length(wu);
                end

                halfW = t(ix2-1)-t(ix1);

            end

            spkWidth    = spkW;
            tr2pk       = p2v;
            halfPkW     = halfW;
            wave        = w;

        end


end

SaveAnalysis(pwd,'WaveformFeatures',{spkWidth tr2pk halfPkW wave maxIx},{'spkWidth' 'tr2pk' 'halfPkW' 'wave' 'maxIx'});