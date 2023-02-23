function [spikeStruct, acceptedUnits] = spikeBehaviorAnalysis_BG(spikes, speed, laps, thetaPeriods, qual2consider, fileinfo)



%% position and speed data

xpos = fileinfo.xyt(:, 1);
ypos = fileinfo.xyt(:, 2);
tpos = fileinfo.xyt(:, 3);

linearPos = fileinfo.xyt2(:, 1);


%% spike data

% initializing the new structure for spike data
spikeStruct = struct('t', [], 'unit', [], 'cluster', [], 'shank', [], 'qclu', [], 'x', [], 'y', [], 'lap', [], 'speed', [], 'linearPos', []);



% including units that are stable and have specific qualities 

isokunit = zeros(length(spikes), 1);

for unit = 1 : length(spikes)
    
    if strcmp(qual2consider, 'pyramidal')
        
        %%% quality number < 4 --> pyramidal, quality number = 9 --> interneurons
        if spikes(unit).q < 4 %&& spikes(unit).isStable 
            isokunit(unit) = 1;
        end
        
    elseif strcmp(qual2consider, 'interneuron')
        
        if spikes(unit).q == 8 %&& spikes(unit).isStable
            isokunit(unit) = 1;
        end
        
    elseif strcmp(qual2consider, 'all')
        
%         if  spikes(unit).isStable
%             isokunit(unit) = 1;
%         end

        isokunit = ones(length(spikes), 1);
    end
    
end


acceptedUnits = find(isokunit); %% index of pyramidal units

nUnits = length(acceptedUnits); %% number of pyramidal units

 
for ii = 1 : nUnits
    
    unit = acceptedUnits(ii); %% ii-th unit
    
    fprintf(1,'currently procceing the unit# %d (%d / %d)\n',  unit, ii, nUnits)
    
    nSpikes = length(spikes(unit).t);
    currSpikes = spikes(unit).t';  
    
    spikeStruct.t = [spikeStruct.t; currSpikes]; 
    
    spikeStruct.unit    = [spikeStruct.unit; ii*ones(nSpikes, 1)];
    spikeStruct.shank   = [spikeStruct.shank; spikes(unit).shank * ones(nSpikes, 1)];
    spikeStruct.cluster = [spikeStruct.cluster; spikes(unit).id * ones(nSpikes, 1)];
    spikeStruct.qclu    = [spikeStruct.qclu; spikes(unit).q * ones(nSpikes, 1)];
    
    spikeStruct.x         = [spikeStruct.x; interp1(tpos, xpos, spikes(unit).t)'];
    spikeStruct.y         = [spikeStruct.y; interp1(tpos, ypos, spikes(unit).t)']; 
    spikeStruct.linearPos = [spikeStruct.linearPos; interp1(tpos, linearPos, spikes(unit).t)'];
    
    spikeStruct.speed     = [spikeStruct.speed; interp1(speed.t, speed.v, spikes(unit).t)'];
 
    currSpikeslap      = zeros(nSpikes, 1);
%     currSpikesThetaIdx = zeros(nSpikes, 1);
    
    for jj = 1: nSpikes
           
        temp = find(spikes(unit).t(jj) > laps(:, 1) & spikes(unit).t(jj) <= laps(:, 2));
        
        if ~isempty(temp)
            currSpikeslap(jj) = laps(temp, 3);
        end 
        
        % check if spikes ocurred during theta
        
%         if ~isempty(find(currSpikes(jj) >= thetaPeriods(:,1) & currSpikes(jj) <= thetaPeriods(:,2)))
%             currSpikesThetaIdx(jj) = 1;
%         end
        
    end
        
%     spikeStruct.theta = [spikeStruct.theta; currSpikesThetaIdx];
    spikeStruct.lap   = [spikeStruct.lap; currSpikeslap];

end


% sorting the spike time stamps

[~, spikeSortIdx] = sort(spikeStruct.t, 'ascend');
spikeStruct = structfun(@(x)(x(spikeSortIdx)), spikeStruct,'UniformOutput',false);


end
