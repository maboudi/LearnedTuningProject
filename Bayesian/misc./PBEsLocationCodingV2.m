function peakDecodedPosition = PBEsLocationCodingV2(eventsBinnedfiring, spatialTunings_LR, spatialTunings_RL, activeUnits, posBinSize, binDur, subfolder)

% EXPLAIN THE FUNCTION


nEvents = size(eventsBinnedfiring, 1);
nPositionBins = size(spatialTunings_LR, 2);


if isempty(activeUnits)
   activeUnits = 1:size(spatialTunings_LR, 1); 
end

%%%%%%%%%%%% Initiation %%%%%%%%%%%%%


%%% posterior probability matrices
%%% the matrix comprised of posterior probability of the track positions (as rows) within each time bins of the event (as columns)


posteriorProbMatrix = cell(nEvents,1); 

%%% 1st array -->  leftward(RL) trajectory only
%%% 2nd array --> rightward(LR) trajectory only
%%% 3rd array -->  summed over the two directions
%%% we can use either the sum or the two trajectories separately for replay
%%% detection (here we have used just the sum)


for pbe = 1: nEvents
    
    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); %%% using just the active units
    
    
    %remove all bins with no firings
    
    sumFiring = sum(currEvent, 1);
    currEvent = currEvent(:, find(sumFiring));
    
    
    %%% Bayesian Sequential Score

    %%% calculate the posterior probabilty matrix without any normalization 
    
    postprobLR = baysDecoder(currEvent, spatialTunings_LR(activeUnits, :), binDur); 
%     postprob_LR = postprob_LR + eps*eps;
    
    postprobRL = baysDecoder(currEvent, spatialTunings_RL(activeUnits, :), binDur); 
%     postprob_RL = postprob_RL + eps*eps;


    posteriorProbMatrix{pbe} = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [nPositionBins, 1]);

end


cnct_posteriorProbMatrix = cell2mat(posteriorProbMatrix');

[~, peakDecodedPosition] = max(cnct_posteriorProbMatrix, [], 1);
peakDecodedPosition = peakDecodedPosition*posBinSize;


end