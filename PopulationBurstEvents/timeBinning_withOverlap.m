function [binnedFiring, binCenters, noFiringUnits] = timeBinning_withOverlap(PBEs, spikes, binDur, binOverlapRatio, fileInfo)

% This function is intended to bin the spikes trains of pyramidal units
% (idnetifies by qclus) within each event. 

% INPUTS
% binOverlapRatio: the overlap between the successive bins

overlap  = binOverlapRatio * binDur;
stepSize = binDur - overlap; 


nPBEs = size(PBEs, 1);
startTs = PBEs(:,1); 
endTs   = PBEs(:,2); 


nUnits = numel(spikes); % total number of units

% okUnits  = fileInfo.stablePyr;
okUnits  = fileInfo.okUnits;
nOKUnits = numel(okUnits); % total number of units which are pyramidal and stable


binnedFiring  = cell(nPBEs, 2);
noFiringUnits = zeros(nPBEs, 1);


for pbe = 1:nPBEs

    PBElen = endTs(pbe)- startTs(pbe); 
    
    binStarts = 0:stepSize:(PBElen - binDur);% for coarse time bins in ms
    if isempty(binStarts)
        binStarts = 0;
    end
    nBins     = length(binStarts);
    
    
    binnedFiring{pbe, 2} = zeros(nUnits, nBins);
    
    for iUnit = 1:nOKUnits
        
        unit = okUnits(iUnit);
        
        spikeTimes = spikes(unit).time;
        spikeTimes = spikeTimes(spikeTimes >= startTs(pbe) & spikeTimes < endTs(pbe));
        
        % % %
        spikeTimes = spikeTimes - startTs(pbe);
        % % %
        
        spikeTrainSegSize = 200; % Due to the long recording, the spike trains were divided into mutiple parts and a summation is calaculated at the end

        nSegments = ceil(numel(spikeTimes)/spikeTrainSegSize);
        
        
        segCounts = zeros(nSegments, nBins);
        for jj = 1:nSegments 

           currentSpikeTimes = spikeTimes((jj-1)*spikeTrainSegSize+1: min(numel(spikeTimes), jj*spikeTrainSegSize));

           firstBin = find(binStarts < currentSpikeTimes(1), 1, 'last');
           lastBin  = find(binStarts > currentSpikeTimes(end), 1, 'first');

           currBins = firstBin:lastBin;

           respSpikeTimes  = (repmat(currentSpikeTimes, [1 numel(currBins)]) - repmat(binStarts(currBins), [numel(currentSpikeTimes) 1])); % spike times in respect to start times of the bins
           ifWithinBin     = respSpikeTimes >= 0 & respSpikeTimes < binDur; % a row corresponding to each spike

           segCounts(jj, currBins) = sum(ifWithinBin, 1);
           
        end

        binnedFiring{pbe, 2}(unit, :) = sum(segCounts, 1); 
     
    end

end



