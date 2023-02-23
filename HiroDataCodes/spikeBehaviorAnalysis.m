function [newSpikes, units]  = spikeBehaviorAnalysis(spikes, laps, ripples, speed, qual2consider, fileinfo)

%%% qual2consider: signifies the units considered to process.
%%% It can be assigned three values as 'pyramidal',
%%% 'interneuron', and 'all'


%% position and speed data

xpos = fileinfo.xyt(:, 1);
ypos = fileinfo.xyt(:, 2);
tpos = fileinfo.xyt(:, 3);

linearPos = fileinfo.xyt2(:, 1);


velocity  = speed.v;
tvelocity = speed.t;


%% spike data


% initializing the new structure for spike data
newSpikes = struct('t', [], 'unit', [], 'cluster', [], 'shank', [], 'qclu', [], 'lap', [], 'x', [], 'y', [], 'linearPos', [], 'speed', [], 'theta', []);

%%% signifying units to be processed

isokunit = zeros(length(spikes), 1);
for unit = 1 : length(spikes)
    
    if strcmp(qual2consider, 'pyramidal')
        
        %%% quality number < 4 --> pyramidal, quality number = 9 --> interneurons
        if spikes(1, unit).quality < 4 && spikes(1, unit).isStable 
            isokunit(unit) = 1;
        end
        
    elseif strcmp(qual2consider, 'interneuron')
        
        if spikes(1, unit).quality == 9 && spikes(1, unit).isStable 
            isokunit(unit) = 1;
        end
        
    elseif strcmp(qual2consider, 'all')
        
        if  spikes(1, unit).isStable 
            isokunit(unit) = 1;
        end
    end
    
end


units = find(isokunit); %% index of pyramidal units
numUnits = length(units); %% number of pyramidal units


for ii = 1 : numUnits
    
    unit = units(ii); %% ii-th unit
    
    fprintf(1,'currently procceing the unit# %d (%d / %d)\n',  unit, ii, numUnits)
    
    nSpikes = length(spikes(unit).time);
    currSpikes = spikes(unit).time;  
    
    newSpikes.t = [newSpikes.t; currSpikes']; 
    
    newSpikes.unit = [newSpikes.unit; ii*ones(nSpikes, 1)];
    newSpikes.shank = [newSpikes.shank; spikes(unit).id(2) * ones(nSpikes, 1)];
    newSpikes.cluster = [newSpikes.cluster; spikes(unit).id(2) * ones(nSpikes, 1)];
    
    newSpikes.qclu = [newSpikes.qclu; spikes(unit).quality * ones(nSpikes, 1)];
    
    newSpikes.x = [newSpikes.x; interp1(tpos, xpos, spikes(unit).time)'];
    newSpikes.y = [newSpikes.y; interp1(tpos, ypos, spikes(unit).time)']; 
    newSpikes.speed = [newSpikes.speed; interp1(tvelocity, velocity, spikes(unit).time)'];
    newSpikes.linearPos = [newSpikes.linearPos; interp1(tpos, linearPos, spikes(unit).time)'];
    
    
    currSpikes_lap = zeros(nSpikes, 1);

    for jj = 1: nSpikes
           
        temp = find(spikes(unit).time(jj) > laps(:, 1) & spikes(unit).time(jj) <= laps(:, 2));
        
        if ~isempty(temp)
            currSpikes_lap(jj) = laps(temp, 3);
        end 
    end
    
    newSpikes.lap = [newSpikes.lap; currSpikes_lap];
  
end


[~, timeSortInd] = sort(newSpikes.t, 'ascend');

newSpikes.t = newSpikes.t(timeSortInd);


newSpikes.unit = newSpikes.unit(timeSortInd);
newSpikes.shank = newSpikes.shank(timeSortInd);
newSpikes.cluster = newSpikes.cluster(timeSortInd);

newSpikes.qclu = newSpikes.qclu(timeSortInd);

newSpikes.x = newSpikes.x(timeSortInd);
newSpikes.y = newSpikes.y(timeSortInd);
newSpikes.linearPos = newSpikes.linearPos(timeSortInd);
newSpikes.speed = newSpikes.speed(timeSortInd);
newSpikes.lap = newSpikes.lap(timeSortInd);



%%% indexing the spikes based on whether they occur during a ripple

runSpikeTimes = newSpikes.t(newSpikes.lap > 0);
ripples = ripples(ripples(:, 1) > runSpikeTimes(1) & ripples(:, 1) < runSpikeTimes(end), :);


nSpikes = length(newSpikes.t);
newSpikes.ripple = zeros(nSpikes, 1);

for ii = 1 : nSpikes
    temp = find(newSpikes.t(ii) >= ripples(:, 1) & newSpikes.t(ii) <= ripples(:, 2));

    if ~isempty(temp)
        newSpikes.ripple(ii) = 1;
    end
end


end
