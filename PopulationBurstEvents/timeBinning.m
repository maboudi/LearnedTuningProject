function [Binnedfiring, nFiringUnits, PBElen] = timeBinning(PBEs, spikes, binDur, fileInfo)

% This function is intended to bin the spikes trains of pyramidal units
% (idnetified by qclus) within each event. 


binSizes = [0.001 binDur]; % binning for two time resolutions (1 ms and 20 ms durations)

nPBEs = numel(PBEs);
startTs = [PBEs.startT];
endTs   = [PBEs.endT];

nUnits = numel(spikes); % total number of units

% okUnits  = fileInfo.stablePyr;
okUnits  = fileInfo.okUnits;
nOKUnits = numel(okUnits); % total number of units which are pyramidal and stable


Binnedfiring  = cell(nPBEs, 2);
nFiringUnits  = zeros(nPBEs, 1);
PBElen        = zeros(nPBEs, 1);

for pbe = 1:nPBEs
    
    for ib = 1:2

        binSize = binSizes(ib);
        
        nBins = ceil((endTs(pbe) - startTs(pbe))/binSize); 
        binEdges = startTs(pbe):binSize:(startTs(pbe)+nBins*binSize);
        
        if ib == 2
           PBElen(pbe) = nBins;
        end
        
%         binEdges = startTs(pbe):binSize:endTs(pbe);
%         nBins    = length(binEdges)-1;
        
        Binnedfiring{pbe, ib} = zeros(nUnits, nBins);

        for iUnit = 1:nOKUnits
            
            unit = okUnits(iUnit);
                        
            spikeTimes = spikes(unit).time;
            spikeTimes = spikeTimes(spikeTimes >= startTs(pbe) & spikeTimes <= endTs(pbe));
            
            if ~isempty(spikeTimes)
                spikeCounts      = histc(spikeTimes, binEdges);
                spikeCounts(end) = [];
                
                Binnedfiring{pbe, ib}(unit, :) = spikeCounts; 
            else
                Binnedfiring{pbe, ib}(unit, :) = zeros(1, nBins); 
            end
            
        end
    end
    
    nFiringUnits(pbe) = numel(find(sum(Binnedfiring{pbe, 1}, 2)));
    
end
    

end