function peakDecodedPosition = PBEsLocationCoding(eventsBinnedfiring, spatialTunings, activeUnits, posBinSize, binDur, subfolder)

% EXPLAIN THE FUNCTION


noEvents = size(eventsBinnedfiring, 1);
noPositionBins = size(spatialTunings, 2);


if isempty(activeUnits)
   activeUnits = 1:size(spatialTunings, 1); 
end

%%%%%%%%%%%% Initiation %%%%%%%%%%%%%


%%% posterior probability matrices
%%% the matrix comprised of posterior probability of the track positions (as rows) within each time bins of the event (as columns)


posteriorProbMatrix = cell(noEvents,1); 

%%% 1st array -->  leftward(RL) trajectory only
%%% 2nd array --> rightward(LR) trajectory only
%%% 3rd array -->  summed over the two directions
%%% we can use either the sum or the two trajectories separately for replay
%%% detection (here we have used just the sum)


for evt = 1: noEvents
    
    currEvent = eventsBinnedfiring{evt, 2}(activeUnits, :); %%% using just the active units
    
     % exclude the time bins with no unit firing
    
    sumFiring = sum(currEvent, 1);
    currEvent = currEvent(:, find(sumFiring));
    
    
    %%% calculate the posterior probabilty matrix without any normalization 
    
    postprob = baysDecoder(currEvent, spatialTunings(activeUnits, :), binDur); 

    posteriorProbMatrix{evt} = postprob ./ repmat(sum(postprob, 1), [noPositionBins, 1]);

end


cnct_posteriorProbMatrix = cell2mat(posteriorProbMatrix');

[~, peakDecodedPosition] = max(cnct_posteriorProbMatrix, [], 1);
peakDecodedPosition = peakDecodedPosition*posBinSize;



end