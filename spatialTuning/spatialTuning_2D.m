function placeTunings = spatialTuning_2D(spikeStruct, qclus, fileinfo, behavior, thetaPeriods, turningPeriods, speed, direction, posBinSize, runSpeedThresh, timeUnit, FileBase)

% EXPLAIN THE FUNCTION


totoalNumofUnits = max(spikeStruct.unit);

% including only the spikes meeting certain criteria

if strcmp(direction, 'LR')
    desiredMod = 0; % even traversals
elseif strcmp(direction, 'RL')
    desiredMod = 1; % odd traversals
elseif strcmp(direction, 'uni')
    desiredMod = [0 1];
end


spikeInd = find(spikeStruct.t >= behavior.time(2,1) & spikeStruct.t < behavior.time(2,2) ... % within the RUN period
                & ismember(spikeStruct.qclu, qclus) ... % only stable pyramidal units are included
                & spikeStruct.speed > runSpeedThresh ... % only the spikes happening when the velocity of the animal is higher than the threshold (usually 10 cm/sec)
                & spikeStruct.theta == 1 ... % only the spikes during high theta
                & spikeStruct.lap > 0 & ismember(mod(spikeStruct.lap, 2), desiredMod)); % limiting to the travels in specific direction

% % excluding that spikes that occurr when the animal stopped or turned
% around in the middle of the track
turningSpikes = [];
for ii = 1: size(turningPeriods, 1); turningSpikes = [turningSpikes; find(spikeStruct.t >= turningPeriods(ii, 1) & spikeStruct.t <= turningPeriods(ii, 2))]; end     

spikeInd = setdiff(spikeInd, turningSpikes);


            
spikePositions = [spikeStruct.y(spikeInd) spikeStruct.x(spikeInd)];
spikeUnit      = spikeStruct.unit(spikeInd);

activeUnits    = unique(spikeUnit); 


% finding position samples pertaining to the maze period and certain direction of travel 

speedatPos = interp1(speed.t, speed.v, fileinfo.xyt(:, 3));

positionIdx    = find(fileinfo.xyt(:, 3) > behavior.time(2,1) & fileinfo.xyt(:, 3) < behavior.time(2,2) ... % within the run period
                        &  speedatPos > runSpeedThresh ... % position bins when animal's velocity is higher than threshold
                        &  fileinfo.xyt2(:, 3) > 0 & ismember(mod(fileinfo.xyt2(:, 3), 2), desiredMod)); % limiting to the travels in specific direction
              
% % include only the position bins within the theta peridos

thetaPositionInd = []; 
for ii = 1: size(thetaPeriods, 1); thetaPositionInd = [thetaPositionInd; find(fileinfo.xyt(:, 3) >= thetaPeriods(ii, 1) & fileinfo.xyt(:, 3) <= thetaPeriods(ii, 2))]; end     

positionIdx = intersect(positionIdx, thetaPositionInd);


turningPositionInd = []; % should be excluded
for ii = 1: size(turningPeriods, 1); turningPositionInd = [turningPositionInd; find(fileinfo.xyt(:, 3) >= turningPeriods(ii, 1) & fileinfo.xyt(:, 3) <= turningPeriods(ii, 2))]; end     

positionIdx = setdiff(positionIdx, turningPositionInd);

                    

xpos = fileinfo.xyt(:, 1);                    
directionXPos = nan(size(xpos));          
directionXPos(positionIdx) = xpos(positionIdx);


ypos = fileinfo.xyt(:, 2);                    
directionYPos = nan(size(ypos));          
directionYPos(positionIdx) = ypos(positionIdx);



% defining the position bins

runIdx = find(fileinfo.xyt(:,3) > behavior.time(2,1) & fileinfo.xyt(:,3)< behavior.time(2,2));


% x pos bins
nXPosBins = floor((max(xpos(runIdx)) - min(xpos(runIdx)))/posBinSize);
xposBins = min(xpos(runIdx)): posBinSize: max(xpos(runIdx)); 

% y pos bins
nYPosBins = floor((max(ypos(runIdx)) - min(ypos(runIdx)))/posBinSize);
yposBins = min(ypos(runIdx)): posBinSize: max(ypos(runIdx)); 


% 2D postion bins
posBinEdges    = cell(1,2);
posBinEdges{1} = yposBins;
posBinEdges{2} = xposBins;


% for calcualting the dwelling times in each position bin
posSamplingPeriod = median(diff(fileinfo.xyt(:, 3)))/timeUnit; % 1/sampling frequency



placeTunings = zeros(nYPosBins+1, nXPosBins+1, totoalNumofUnits); % which dimensions of the tuning matrix correspond to x or y positions
 
% posBinDwell          = hist3([directionYPos directionXPos], 'Edges', posBinEdges);
% posBinDwell          = posBinDwell*posSamplingPeriod;

% totalNumSpikes = zeros(length(activeUnits), 1);

for ii = 1: length(activeUnits)
    
    unit               = activeUnits(ii);
    unitSpikePositions = spikePositions(spikeUnit == unit, :);
    totalNumSpikes(ii) = length(unitSpikePositions);
    
    
    posBinnSpikes        = hist3(unitSpikePositions, 'Edges', posBinEdges);
    
    posBinDwell          = hist3([directionYPos directionXPos], 'Edges', posBinEdges);
    posBinDwell          = posBinDwell*posSamplingPeriod;
    
    
    prim_tunings = posBinnSpikes ./ posBinDwell; % unsmoothed tunings

    prim_tunings(isnan(prim_tunings)) = 0;
    prim_tunings(isinf(prim_tunings)) = 0;
    

    win = gausswindow(3,3);
    win2D = win' * win;
    
    placeTunings(: , :, unit) = conv2(prim_tunings, win2D, 'same');
end

figure; imagesc(posBinDwell); colormap('hot'); set(gca , 'YDir', 'normal')

plotThePlaceFields(placeTunings, spikePositions, spikeUnit, activeUnits, fileinfo, behavior, direction, posBinSize, FileBase)


end
