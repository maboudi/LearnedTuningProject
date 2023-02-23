function spikes = spatialTuning_1D_V2(spikes, selectlaps, thetaPeriods, turningPeriods, runSpeedThresh, direction, posBinSize, fileInfo)



nUnits = numel(spikes);

% including only the spikes meeting certain criteria

if strcmp(direction, 'LR')
    desiredMod = 0; % even traversals
elseif strcmp(direction, 'RL')
    desiredMod = 1; % odd traversals
elseif strcmp(direction, 'uni')
    desiredMod = [0 1];
end

behavior = fileInfo.behavior;
speed    = fileInfo.speed;


% finding position samples pertaining to the maze period and certain direction of travel 

speedatPos   = interp1(speed.t, speed.v, fileInfo.xyt(:, 3));

positionIdx  = find(fileInfo.xyt(:, 3) > behavior.time(2,1) & fileInfo.xyt(:, 3) < behavior.time(2,2) ... % within the run period
                    & speedatPos > runSpeedThresh ... % position bins when animal's velocity is higher than threshold
                    & fileInfo.linearPos(:, 3) > 0 & ismember(mod(fileInfo.linearPos(:, 3), 2), desiredMod) ... % limiting to the travels in specific direction
                    & ismember(fileInfo.linearPos(:, 3), selectlaps(:,3))); % limiting to a selected number of laps


                
% include only the position bins within the theta peridos

if ~isempty(thetaPeriods)
    thetaPositionInd = []; 
    for ii = 1: size(thetaPeriods, 1)
        thetaPositionInd = [thetaPositionInd; find(fileInfo.xyt(:, 3) >= thetaPeriods(ii, 1) & fileInfo.xyt(:, 3) <= thetaPeriods(ii, 2))];
    end     
    positionIdx = intersect(positionIdx, thetaPositionInd);
end


if ~isempty(turningPeriods)
    turningPositionInd = []; % should be excluded
    for ii = 1: size(turningPeriods, 1)
        turningPositionInd = [turningPositionInd; find(fileInfo.xyt(:, 3) >= turningPeriods(ii, 1) & fileInfo.xyt(:, 3) <= turningPeriods(ii, 2))];
    end
    positionIdx = setdiff(positionIdx, turningPositionInd);
end


linearPos = fileInfo.linearPos(:, 1);                    
directionLinearPos = nan(size(linearPos));          
directionLinearPos(positionIdx) = linearPos(positionIdx);




% defining the position bins

nPosBins = floor((max(linearPos) - min(linearPos))/posBinSize);
posBinEdges = min(linearPos): posBinSize: max(linearPos); % center of the position bins
posSamplingPeriod = median(diff(fileInfo.xyt(:, 3)))/fileInfo.timeUnit; % 1/sampling frequency 

% spatialTunings = zeros(nUnits, nPosBins);
% peakPosBin     = zeros(nUnits, 1);

for unit = 1: nUnits
    

    if size(spikes(unit).time, 1) ~= size(spikes(unit).lap, 1)
        spikes(unit).lap = spikes(unit).lap';
    end


    spikeInd = find(spikes(unit).time >= behavior.time(2,1) & spikes(unit).time < behavior.time(2,2) ... % within the RUN period
                    & spikes(unit).speed > runSpeedThresh ... % only the spikes happening when the velocity of the animal is higher than the threshold (usually 10 cm/sec)          
                    & spikes(unit).lap > 0 & ismember(mod(spikes(unit).lap, 2), desiredMod) ...  % limiting to the travels in specific direction
                    & ismember(spikes(unit).lap, selectlaps(:,3))); % limiting to a selected number of laps
    %                 & spikeStruct.theta == 1); % only spike during high theta power
    
    % % excluding the spikes that occurr when the animal stopped or turned
    % around in the middle of the track

    turningSpikes = [];
    for ii = 1: size(turningPeriods, 1)
        turningSpikes = [turningSpikes; find(spikes(unit).time >= turningPeriods(ii, 1) & spikes(unit).time <= turningPeriods(ii, 2))]; 
    end
    
    spikeInd = setdiff(spikeInd, turningSpikes);
    
    
    
    spikePositions = spikes(unit).linearPos(spikeInd);
    
    posBinnSpikeCnts      = histc(spikePositions, posBinEdges);
    posBinnSpikeCnts(end) = [];
    
    posBinDwelltime       = histc(directionLinearPos, posBinEdges) * posSamplingPeriod;
    posBinDwelltime(end)  = [];
    
    
    
    if size(posBinnSpikeCnts, 2) > 1
        posBinnSpikeCnts = posBinnSpikeCnts';
    end
    
    unsmoothed_tunings = posBinnSpikeCnts ./ posBinDwelltime;

    unsmoothed_tunings(isnan(unsmoothed_tunings)) = 0;
    unsmoothed_tunings(isinf(unsmoothed_tunings)) = 0;
    
    
    win = gausswindow(2,5);
    
    spikes(unit).CVspatialTuning.(direction) = conv(unsmoothed_tunings, win, 'same');
    [~, spikes(unit).CVpeakPosBin.(direction)] = max(spikes(unit).CVspatialTuning.(direction));
    
end

end