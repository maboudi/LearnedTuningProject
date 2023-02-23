function spikeStruct = spikeBehaviorAnalysis3(spikes, speed, laps, thetaPeriods, qual2consider, fileinfo)



%% position and speed data

xpos = fileinfo.xyt(:, 1);
ypos = fileinfo.xyt(:, 2);
tpos = fileinfo.xyt(:, 3);

linearPos = fileinfo.xyt2(:, 1);


%% spike data


% initializing the new structure for spike data
spikeStruct = struct('t', [], 'unit', [], 'qclu', [], 'x', [], 'y', [], 'lap', [], 'speed', [], 'linearPos', [], 'theta', []);

% , 'cluster', [], 'shank', , 'theta', []

% For now we are not considering the stability of unit as a inclusion
% criterion

switch qual2consider
    case 'pyr'
        spikesIncluded = spikes.PyrIDs;
    case 'int'
        spikesIncluded = spikes.IntIDs;
    case 'all'
        spikesIncluded = sort([spikes.PyrIDs; spikes.IntIDs]);
end

spikeIdxIncluded = find(ismember(spikes.SpikeIDs, spikesIncluded));

spikeStruct.t = spikes.SpikeTimes(spikeIdxIncluded);
spikeunitIDs = spikes.SpikeIDs(spikeIdxIncluded);
spikeStruct.qclu = 1*ismember(spikeunitIDs, spikes.PyrIDs) + 9*ismember(spikeunitIDs, spikes.IntIDs);

spikeStruct.unit = zeros(length(spikeIdxIncluded), 1);
spikeStruct.lap = zeros(length(spikeIdxIncluded), 1);
spikeStruct.theta = zeros(length(spikeIdxIncluded), 1);

for ii = 1:length(spikesIncluded)
    
    unit = spikesIncluded(ii);
    
    unitSpikesIdx = find(spikeunitIDs == unit);
    unitSpikeTimes = spikeStruct.t(unitSpikesIdx);
    
    spikeStruct.unit(unitSpikesIdx) = ii;
    numSpikes = length(unitSpikesIdx);
    
    
    currSpikeslap = zeros(numSpikes, 1);
    
    for jj = 1: numSpikes
        
        % lap
        temp = find(unitSpikeTimes(jj) > laps(:, 1) & unitSpikeTimes(jj) <= laps(:, 2));
        
        if ~isempty(temp)
            currSpikeslap(jj) = laps(temp, 3);
        end 
        
        % check if spikes occurred during theta
        
        if ~isempty(find(unitSpikeTimes(jj) >= thetaPeriods(:,1) & unitSpikeTimes(jj) <= thetaPeriods(:,2)))
            spikeStruct.theta(unitSpikesIdx(jj)) = 1;
        end
        
    end
    
    spikeStruct.lap(unitSpikesIdx) = currSpikeslap;

end

spikeStruct.x         = interp1(tpos, xpos, spikeStruct.t);
spikeStruct.y         = interp1(tpos, ypos, spikeStruct.t);
spikeStruct.speed     = interp1(speed.t, speed.v, spikeStruct.t);
spikeStruct.linearPos = interp1(tpos, linearPos, spikeStruct.t);


end
