function [binnedSpikeCounts, noFiringUnits, PBElen] = finalBinningResult(PBEs, spikes, binDur, fileInfo)

%% binnig the spikes within an event and qualifying the event based on number of active units and length


%%% binning the spikes within each event

[binnedSpikeCounts, noFiringUnits, PBElen] = timeBinning(PBEs, spikes, binDur, fileInfo);


% Remove the flanking zero-firing periods(silent bins) from the start
% and end of each event. 

% [binnedSpikeCounts, PBElen, startTs, endTs] = removeSideSilentBins(binnedSpikeCounts_allPBEs, PBEs(:, 1), binDur);

% secondaryPBEs = [startTs endTs PBEs(:, 3:4)]; 


end
