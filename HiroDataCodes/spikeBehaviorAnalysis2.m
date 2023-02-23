function [spikeStruct, acceptedUnits] = spikeBehaviorAnalysis2(spikes, speed, laps, qual2consider, fileinfo)

% [spikeStruct, acceptedUnits] = spikeBehaviorAnalysis2(spikes, speed, laps, thetaPeriods, qual2consider, fileinfo)


%% position and speed data

xpos = fileinfo.xyt(:, 1);
ypos = fileinfo.xyt(:, 2);
tpos = fileinfo.xyt(:, 3);

linearPos = fileinfo.xyt2(:, 1);


%% spike data

% initializing the new structure for spike data
spikeStruct = struct('t', [], 'unit', [], 'cluster', [], 'shank', [], 'qclu', [], 'x', [], 'y', [], 'lap', [], 'speed', [], 'linearPos', []);


%%% quality number <= 3 --> pyramidal, quality number = 9 --> interneurons
if strcmp(qual2consider, 'pyramidal')
    spikes = spikes([spikes.quality] <= 3 & [spikes.StablePrePostFiner] == 1);

elseif strcmp(qual2consider, 'interneuron')
    spikes = spikes([spikes.quality] == 9 & [spikes.StablePrePostFiner] == 1);

elseif strcmp(qual2consider, 'all')
    spikes = spikes([spikes.StablePrePostFiner] == 1);

end



nUnits = numel(spikes);

for unit = 1:nUnits
    
    fprintf(1,'currently procceing the unit# %d (out of %d)\n',  unit, nUnits)
    
    nSpikes = length(spikes(unit).time);
    
    spikeStruct.t = [spikeStruct.t; spikes(unit).time]; 
    
    spikeStruct.unit    = [spikeStruct.unit; unit*ones(nSpikes, 1)];
    spikeStruct.shank   = [spikeStruct.shank; spikes(unit).id(2) * ones(nSpikes, 1)];
    spikeStruct.cluster = [spikeStruct.cluster; spikes(unit).id(2) * ones(nSpikes, 1)];
    spikeStruct.qclu    = [spikeStruct.qclu; spikes(unit).quality * ones(nSpikes, 1)];
    
    spikeStruct.x         = [spikeStruct.x; interp1(tpos, xpos, spikes(unit).time)];
    spikeStruct.y         = [spikeStruct.y; interp1(tpos, ypos, spikes(unit).time)]; 
    spikeStruct.linearPos = [spikeStruct.linearPos; interp1(tpos, linearPos, spikes(unit).time)];
    
    spikeStruct.speed     = [spikeStruct.speed; interp1(speed.t, speed.v, spikes(unit).time)];
 
    currSpikeslap      = zeros(nSpikes, 1);
    
    for jj = 1: nSpikes
        temp = find(spikes(unit).time(jj) > laps(:, 1) & spikes(unit).time(jj) <= laps(:, 2));
        if ~isempty(temp)
            currSpikeslap(jj) = laps(temp, 3);
        end 
    end
        
    spikeStruct.lap   = [spikeStruct.lap; currSpikeslap];
end


% sorting the spike time stamps

[~, spikeSortIdx] = sort(spikeStruct.t, 'ascend');
spikeStruct = structfun(@(x)(x(spikeSortIdx)), spikeStruct,'UniformOutput',false);


end
