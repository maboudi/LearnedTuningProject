function [runBinnedfiring, posbinIdx, linposcenters] = binRunData_1D(spikes, activePeriods, runSpeedThresh, turningPeriods, binDur, binOverlapRatio, posBinSize, fileInfo)


% This function is intended to divide the RUN neuronal activity to activity bouts within time bins with
% duration determined by binDur. Also, it calculates index of the position
% bin corresponding to each time bin.


% linearized track positions are used here (its difference with binRunData)

% (we might get rid of this later)


behavior = fileInfo.behavior;
speed    = fileInfo.speed;

linearPos = fileInfo.linearPos(:,1);
tpos      = fileInfo.linearPos(:, 2)*1/fileInfo.lfpSampleRate;

nanIdx = find(isnan(linearPos));

linearPos_woNans = linearPos;
linearPos_woNans(isnan(linearPos)) = [];

tpos_woNans = tpos;
tpos_woNans(isnan(linearPos)) = [];

pos_at_nans = interp1(tpos_woNans, linearPos_woNans, tpos(nanIdx));

linearPos(nanIdx) = pos_at_nans;




%% binning the run data

runBinnedfiring2 = timeBinning_withOverlap(activePeriods , spikes, binDur, binOverlapRatio, fileInfo); % the theta periods
runBinnedfiring2 = runBinnedfiring2(:, 2); % we need only data with coarse time bins


% Calculate the track positions (continuous not discretized) corresponding to each time bin
% exclude the positon bins with speed lower than the runSpeedThresh
% and also the position bins within the turning periods or on the platforms
% in the middle of a lap


periods   = activePeriods; % * 1/fileInfo.lfpSampleRate;

runBinPos = cell(size(periods, 1), 1);

stepSize = (1-binOverlapRatio)*binDur;

for ii = 1 : size(periods, 1)
    
    binStarts = periods(ii, 1):stepSize:(periods(ii, 2)-binDur);
    binCenters = binStarts + binDur/2;
    
    runBinPos{ii}(:,1)    = interp1(tpos, linearPos, binCenters); %% position at the center of each time bin
    
    runBinPos{ii}(:,2)    = 1:length(runBinPos{ii}(:,1)); % indexing, before excluding bins 
    
    
    % exclude bins with animal speed lower then threshold (or the time bins during the turning periods)
    
    isTuringAtBin = nan(length(binCenters), 1);
    
    if ~isempty(turningPeriods)

        for ibin = 1:length(binCenters)

            aa = find(turningPeriods(:, 1) < binCenters(ibin) & turningPeriods(:, 2) > binCenters(ibin), 1, 'first');

            if ~isempty(aa)
                isTuringAtBin(ibin) = 1;
            end

        end
    end
    
    
    if ~isempty(runSpeedThresh)
        
        speedAtBin = interp1(speed.t/fileInfo.lfpSampleRate, speed.v, binStarts+binDur/2);  
        
        ifLowSpeed_or_Turning = find((speedAtBin < runSpeedThresh) | isTuringAtBin' == 1); 
    
        runBinPos{ii}(ifLowSpeed_or_Turning, :)        = [];
        runBinnedfiring2{ii}(:, ifLowSpeed_or_Turning) = [];
    end
    
end


%% Indexing the continuous positions based on discretized position bins

% Defining the position bins; the same we did in calculating the place
% fields


mazePeriodIdx = fileInfo.linearPos(:, 2) > behavior.time(2,1) & fileInfo.linearPos(:, 2) < behavior.time(2,2);
linPos        = fileInfo.linearPos(mazePeriodIdx, 1);


linposBins      = min(linPos): posBinSize: max(linPos); % center of the position bins
linposBins(end) = max(linPos);

linposcenters      = linposBins + posBinSize/2;
linposcenters(end) = [];



posbinIdx2 = cell(size(periods, 1), 1);
for ii = 1 : size(periods, 1)
    
    [~, posbinIdx2{ii}(:,1)] = histc(runBinPos{ii}(:, 1), linposBins);
    posbinIdx2{ii}(:,2)      = runBinPos{ii}(:, 2);

end


runBinnedfiring = cell2mat(runBinnedfiring2'); % concatenating all the run binned spikes
posbinIdx       = cell2mat(posbinIdx2);


end
