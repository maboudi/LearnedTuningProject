function [positionProb2D, peakPosBin, distFromPrevBin, PBETrajsegments] = BayesianPosReconstrcution2D_v2(eventsBinnedfiring, tuning, binDur, posBinSize, frame2framedisThresh, noStepsThresh, start2endThresh, plotTraj)

% This function is made to calcualte the virtual trajectory through the
% environment for each PBE

%%% inputs

% eventBinnedfiring is the 20 milisecond binned number of spikes for each
% PBE

% Tuning is a 3-dimensional matrix consisting of the 2D rate tuning of different neurons 

% frame2framedisThresh: the maximum tolerable distance to consider to
% successive bins decoding adjacent track positions. 

% noStepsThresh: minimum number of steps to consider a segment (of
% successive bins decoding adjacent track positions) as a (virtual)
% trajectory event

% start2endThresh: minimum distance on the track (open field, 2D env) covered to consider a segment as a trajectory event 


%%% outputs

% positionProb2D: for each PBE is a 3d matrix with the 3rd dimension
% corresponding to the individual time bins and the first and 2nd dimensions
% corresponding to the decoding likelihhood of each 2d position pixel

% peakPosBin: within each PBE, x and y coordinates of peak position pixel for each time bin 
% distFromPrevBin: distance between the successive time bins

% PBETrajsegments: trajectory segments within each PBE. There can be more
% than one within each PBE. Each trajectory is characterized by start bin,
% end bin, number of steps, its comparison with the determined threshold, coverd length of the env and its comparison with the determined threshold

tuning = tuning + 0.001;

noPBEs = size(eventsBinnedfiring, 1);

[noYbins, noXbins, noUnits] = size(tuning);

tuning1D = reshape(tuning, [noYbins*noXbins noUnits 1])';
noPositionBins = size(tuning1D, 2);

positionProb2D = cell(noPBEs, 1);
peakPosBin = cell(noPBEs, 1);
distFromPrevBin = cell(noPBEs, 1);
PBETrajsegments = cell(noPBEs, 1);

for pbe = 1: noPBEs
    
    currEvent = eventsBinnedfiring{pbe, 2}; %%% using just the active units
    noTimeBins = size(currEvent, 2);

    % likelihood of each position bin given the population firing within
    % the bin and the neuronal tuning across the position bins
    
    posLikelihoods = baysDecoder(currEvent, tuning1D, binDur); 
    
    %%% normalizing the probability distribution across the positions
    
    posLikelihoods_normalized = posLikelihoods./repmat(sum(posLikelihoods, 1), [noPositionBins, 1]);
    
    positionProb2D{pbe} = reshape(posLikelihoods_normalized, [noYbins noXbins noTimeBins]);
    
    peakPosBin{pbe} = zeros(noTimeBins, 2);
    distFromPrevBin{pbe} = zeros(noTimeBins, 1);
    
    for bin = 1:noTimeBins
        
        currentBinMap = positionProb2D{pbe}(:,:, bin);
        [r, c] = find(currentBinMap == max(currentBinMap(:)));
        
        try
            peakPosBin{pbe}(bin, :) = [r c]; % y and x in first and second dimesions
        catch
             peakPosBin{pbe}(bin, :) = [inf inf];
        end

        if bin > 1
           distFromPrevBin{pbe}(bin) = norm(peakPosBin{pbe}(bin, :) - peakPosBin{pbe}(bin-1, :))*posBinSize;
        end
        
    end
    
    
    distFromPrevBin{pbe} = [distFromPrevBin{pbe} (1:noTimeBins)']; % added the second column to keep track of the bins' indices (later we'll need to calculate the number of steps)
    
    distFromPrevBin{pbe}(distFromPrevBin{pbe}(:,1) < 1, :) = []; % removing all the stationary segments
    distFromPrevBin{pbe} = [0 1; distFromPrevBin{pbe}]; % adding the first bin again (was removed by the operation above)
    
    
    isDistLessThresh = distFromPrevBin{pbe}(2:end, 1) <= frame2framedisThresh;
    isDistLessThresh = [0; isDistLessThresh]; % add zero for the first bin
    
    % find the consequitive time bins with bin to bin distance less than
    % the threshold
    temp = abs(diff([0; isDistLessThresh])); 
    
    startOrEndBins = find(temp > 0)-1; % decrement by 1 
    
    if mod(length(startOrEndBins), 2) == 1
        startOrEndBins = [startOrEndBins; size(distFromPrevBin{pbe}, 1)];
    end

    numSegs = length(startOrEndBins)/2;

    PBETrajsegments{pbe} = zeros(numSegs, 6); % start and end of each (potentially) trajectory segments, without filtering based on the length yet
    for seg = 1: numSegs
        
        segStart  = startOrEndBins((seg-1)*2+1);
        segEnd    = startOrEndBins(seg*2);
        segNSteps = segEnd - segStart;
        
        segDistanceCovered = norm(peakPosBin{pbe}(distFromPrevBin{pbe}(segEnd, 2), :) - peakPosBin{pbe}(distFromPrevBin{pbe}(segStart, 2), :))*posBinSize;
        
        PBETrajsegments{pbe}(seg, :) = [distFromPrevBin{pbe}(segStart, 2) distFromPrevBin{pbe}(segEnd, 2) ...
                                        segNSteps segNSteps>=noStepsThresh ... % trajectory's number of steps
                                        segDistanceCovered segDistanceCovered>=start2endThresh]; % trajectory's covered distance
                                    
                                    
    end
    
    PBETrajsegments{pbe}(PBETrajsegments{pbe}(:, 5) < 1, :) = [];
   
end



% statistical tests
% Two shuffling methods were used by Pfeiffer and Foster 2013: shuffling the cell identities and shuffling the place tuniings in both x and y directions




end