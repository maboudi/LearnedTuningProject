function [weightedCorr, jumpDistance, coveredLength] = replayevaluation_v2(postPr)

% This function evaluates PBEs (frames) regarding to criteria: Weihgted
% correlation and jump distance (3 different methods of measurement)


nPosBins = size(postPr, 1);


%% weighted correlation

weightedCorr = calWeightedCorr(postPr);



%% jump distance

[~, peakPos] = max(postPr, [], 1);


% eliminating non-firing time bins
sumProb                = sum(postPr, 1);
peakPos(sumProb < 0.1) = [];


% calculating jump distances



if numel(peakPos) < 2
   
   maxJump     = 1;
   medianJump  = 1;
   normMaxJump = 1;
   
else

    jumpDistances = abs(diff(peakPos));
    
    maxJump       = max(jumpDistances)/nPosBins; % maximum jump distance as a ratio of the track length
    medianJump    = median(jumpDistances)/nPosBins; % median jump distance as a ratio of the track length


    % normalized maximum jump distance (data as a percentile of shuffle distribution:within-PBE time bin shuffle (ref: Dragoi2019 Science))

    nShuffles      = 500; 
    shuffleMaxJump = zeros(nShuffles, 1);

    for ii = 1: nShuffles

        shufflePeakPos      = peakPos(randperm(length(peakPos)));
        shuffleJumpDistance = abs(diff(shufflePeakPos));


        shuffleMaxJump(ii)  = max(shuffleJumpDistance)/nPosBins;

    end

    normMaxJump  = length(find(shuffleMaxJump < maxJump))/nShuffles;

end


jumpDistance = [maxJump normMaxJump medianJump];


%% covered length of the track

% coveredLength = (max(peakPos) - min(peakPos))/nPosBins;
coveredLength = nan;

end