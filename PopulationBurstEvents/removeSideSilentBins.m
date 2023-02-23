function [eventsBinnedfiring2, pbeLen, newStartTs, newEndTs] = removeSideSilentBins(eventsBinnedfiring, startTs, binDur)

% This function removes the silent (non-firing) bins from the start and
% end of each event.

nPBEs = size(eventsBinnedfiring, 1); %% original pbes

eventsBinnedfiring2 = cell(size(eventsBinnedfiring));
pbeLen = zeros(nPBEs, 1);

newStartTs = zeros(nPBEs, 1);
newEndTs   = zeros(nPBEs, 1);


% Remove the surrounding silent bins from each event

for pbe = 1:nPBEs
    
    pbeData = eventsBinnedfiring{pbe, 2}; % working with the coarse time bins (e.g., 20 ms for PBEs)
    
    spikeCounts = sum(pbeData, 1); % number of spikes within each bin
    
    firstBin = find(spikeCounts > 0, 1, 'first');
    lastBin  = find(spikeCounts > 0, 1, 'last');
    
    
    if ~isempty(firstBin) && ~isempty(lastBin) 
        
        pbeLen(pbe) = lastBin - firstBin + 1; %% length of the truncated event
        
        % create new truncated PBEs
        
        eventsBinnedfiring2{pbe, 2} = pbeData(:, firstBin:lastBin);
        eventsBinnedfiring2{pbe, 1} = eventsBinnedfiring{pbe, 1}(:, floor((firstBin-1)*binDur*1000+1) : floor(lastBin*binDur*1000));
        
        
        % update the start and end time of the PBEs
        newStartTs(pbe) = startTs(pbe) + (firstBin-1) * binDur;
        newEndTs(pbe)   = startTs(pbe) + lastBin * binDur;

    end
    
    
end


end
