function [posteriorProbMatrix, posteriorProbMatrix_nn, weightedCorr, weightedDistCorr, BDseqscore, weightedCorr_null, weightedDistCorr_null] = BDreplayDetect_weightedDist(eventsBinnedfiring, tuningRL, tuningLR, activeUnits, binDur)



% This function is intended to calculate the bayesian replay score for each
% event. The significance of the replay scores is evaluated through
% comparing it with a distribution of scores resulted from shuffle
% datasets, i.e. column cycle shuffle, time bin shuffle, etc. 


%%% for now we skip the replay type-we get back to it later ..


% first replace all of zeros in tuning with a small value (at the level of smallest non-zero)

tuningRL(~tuningRL) = 1e-4;
tuningLR(~tuningLR) = 1e-4;


nEvents       = size(eventsBinnedfiring, 1);
nPositionBins = size(tuningRL, 2);

% noshufflingMethods = 2; %% number of methods for shuffling (time swap and column cycle shuffle here)
nShuffles = 500; %% number of shuffles in each of the methods 


if isempty(activeUnits)
   activeUnits = 1:size(tuningRL, 1); 
end

%%%%%%%%%%%% Initiation %%%%%%%%%%%%%


%%% posterior probability matrices
%%% the matrix comprised of posterior probability of the track positions (as rows) within each time bins of the event (as columns)


posteriorProbMatrix    = cell(nEvents, 1); 
posteriorProbMatrix_nn = cell(nEvents, 1);

%%% 1st array -->  leftward(RL) trajectory only
%%% 2nd array --> rightward(LR) trajectory only
%%% 3rd array -->  summed over the two directions
%%% we can use either the sum or the two trajectories separately for replay
%%% detection (here we have used just the sum)



weightedCorr     = zeros(nEvents, 1);
weightedDistCorr = zeros(nEvents, 1);

baseStruct  = struct('weightedCorr', [], 'weightedDistCorr', []);
baseStruct.weightedCorr     = struct('prctilescore', zeros(nEvents, 1), 'zscore', zeros(nEvents, 1));
baseStruct.weightedDistCorr = struct('prctilescore', zeros(nEvents, 1), 'zscore', zeros(nEvents, 1));

BDseqscore  = struct('wPBEtimeswap', baseStruct, 'unitIDshuffle', baseStruct);

fprintf('\n Calculating replay scores of actuall PBEs')

tic

% fitLineDictionary = generateLineDictionary(nPositionBins, 1:30);

for pbe = 1:nEvents


    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); %%% using just the active units
    nTimebins = size(currEvent, 2);
    % exclude the time bins with no unit firing
    
    sumFiring = sum(currEvent, 1);
    idx       = find(sumFiring == 0);
    
    %%% Bayesian Sequential Score

    %%% calculate the posterior probabilty matrix without any normalization 
    
    postprobRL = baysDecoder(currEvent, tuningRL(activeUnits, :), binDur); 
    postprobLR = baysDecoder(currEvent, tuningLR(activeUnits, :), binDur);
    
    posteriorProbMatrix_nn{pbe} = postprobRL + postprobLR;
    posteriorProbMatrix_nn{pbe}(:, idx) = [];
    %%% Summed probability distribution, every column is normalized to one
    %%% (this is the same method as used in Davidson et al. 2009)
    
    posteriorProbMatrix{pbe} = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [nPositionBins, 1]);
    posteriorProbMatrix{pbe}(:, idx) = [];
    
    
    % weighted correlation
    weightedCorr(pbe) = calWeightedCorr(posteriorProbMatrix{pbe});
    
    % weighted distance correlation
    weightedDistCorr(pbe) = calWeightedDistCorr(posteriorProbMatrix{pbe});
    

    
    if mod(pbe, 50) == 0
        fprintf('.')
    end

end

toc



%% calculating shuffle distributions 

weightedCorr_null     = zeros(nEvents, nShuffles, 2);
weightedDistCorr_null = zeros(nEvents, nShuffles, 2);


% %% 1) within-PBE time swap
% 
% 
% 
% fprintf('\n Measuring statistical significances by comparing against within-PBE time swap')
% 
% for pbe = 1: nEvents
%     
%     currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); 
%     
%     currPosteriorProbMatrix = posteriorProbMatrix(pbe);
%     
%     
%     weightedCorr_null_pbe     = zeros(nShuffles, 1);
%     weightedDistCorr_null_pbe = zeros(nShuffles, 1);
%     
%     parfor sn = 1: nShuffles
%         
%         posteriorProbMatrix_shuffle = genTimeSwap(currPosteriorProbMatrix);
%         posteriorProbMatrix_shuffle = posteriorProbMatrix_shuffle{1};
% 
%         % weighted correlation
%         weightedCorr_null_pbe(sn) = calWeightedCorr(posteriorProbMatrix_shuffle);
% 
%         % weighted distance correlation
%         weightedDistCorr_null_pbe(sn) = calWeightedDistCorr(posteriorProbMatrix_shuffle);
% 
%     end
%     
%     weightedCorr_null(pbe, :, 1)     = permute(weightedCorr_null_pbe, [2, 1]);
%     weightedDistCorr_null(pbe, :, 1) = permute(weightedDistCorr_null_pbe, [2, 1]);
%     
%     % sequence scores
%     
%     BDseqscore.wPBEtimeswap.weightedCorr.prctilescore(pbe) = length(find(abs(weightedCorr_null(pbe, :, 1)) < abs(weightedCorr(pbe))))/nShuffles * 100;
%     BDseqscore.wPBEtimeswap.weightedCorr.zscore(pbe)       = (abs(weightedCorr(pbe)) - mean(abs(weightedCorr_null(pbe, :, 1)))) / std(abs(weightedCorr_null(pbe, :, 1)));
% 
%     BDseqscore.wPBEtimeswap.weightedDistCorr.prctilescore(pbe) = length(find(abs(weightedDistCorr_null(pbe, :, 1)) < abs(weightedDistCorr(pbe))))/nShuffles * 100;
%     BDseqscore.wPBEtimeswap.weightedDistCorr.zscore(pbe)       = (abs(weightedDistCorr(pbe)) - mean(abs(weightedDistCorr_null(pbe, :, 1)))) / std(abs(weightedDistCorr_null(pbe, :, 1)));
%     
%     
%     if mod(pbe, 50) == 0
%        fprintf('\n event #%d: pscore wc %.1f wdc %.1f  and zscore wc %.2f wdc %.2f', ...
%            pbe, BDseqscore.wPBEtimeswap.weightedCorr.prctilescore(pbe, 1), BDseqscore.wPBEtimeswap.weightedDistCorr.prctilescore(pbe, 1), BDseqscore.wPBEtimeswap.weightedCorr.zscore(pbe, 1), BDseqscore.wPBEtimeswap.weightedDistCorr.zscore(pbe, 1))
%     end
%     
% end
% 
% 
% %% 2) unit ID shuffle
% % There is much redundancy between the procedures for the two shuffle
% % methods, can I make it one?
% 
% fprintf('\n Measuring statistical significances by comparing against unit ID shuffle')
% % shuffleInclusionFRthresh = 1;
% for pbe = 1:nEvents
%     
%     currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); %%% using just the active units
%     nTimebins = size(currEvent, 2);
%     
%     pbeActUnits = find(sum(currEvent, 2));
%     
%     
%     % calculate new spatial tunings with shuffled cell IDs
%         
%     [tuningRLshuffles, tuningLRshuffles] = tuningsUnitIDshuffle_normalize(tuningRL, tuningLR, pbeActUnits, nShuffles, 'same-peak'); % the shuffle type can be 'regular' or 'same-peak'
%     
%     sumFiring = sum(currEvent, 1);
%     idx       = find(sumFiring == 0);
%         
%     weightedCorr_null_pbe     = zeros(nShuffles, 1);
%     weightedDistCorr_null_pbe = zeros(nShuffles, 1);
%     
%     
%     parfor sn = 1:nShuffles % shuffle index
% 
%         %%% calculate the posterior probabilty matrix without any normalization 
% 
%         postprobRL = baysDecoder(currEvent, tuningRLshuffles(activeUnits, :, sn), binDur); 
%         postprobLR = baysDecoder(currEvent, tuningLRshuffles(activeUnits, :, sn), binDur);
% 
%         posteriorProbMatrix_PFshuffle = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [nPositionBins, 1]);
%         posteriorProbMatrix_PFshuffle(:, idx) = [];
%         
%         weightedCorr_null_pbe(sn)     = calWeightedCorr(posteriorProbMatrix_PFshuffle);
%         weightedDistCorr_null_pbe(sn) = calWeightedDistCorr(posteriorProbMatrix_PFshuffle);
%         
%     end
%         
%     weightedCorr_null(pbe, :, 2)     = permute(weightedCorr_null_pbe, [2, 1]);
%     weightedDistCorr_null(pbe, :, 2) = permute(weightedDistCorr_null_pbe, [2, 1]);
%     
%     
%     % sequence scores
%     
%     BDseqscore.unitIDshuffle.weightedCorr.prctilescore(pbe) = length(find(abs(weightedCorr_null(pbe, :, 2)) < abs(weightedCorr(pbe))))/nShuffles * 100;
%     BDseqscore.unitIDshuffle.weightedCorr.zscore(pbe)       = (abs(weightedCorr(pbe)) - mean(abs(weightedCorr_null(pbe, :, 2)))) / std(abs(weightedCorr_null(pbe, :, 2)));
% 
%     BDseqscore.unitIDshuffle.weightedDistCorr.prctilescore(pbe) = length(find(abs(weightedDistCorr_null(pbe, :, 2)) < abs(weightedDistCorr(pbe))))/nShuffles * 100;
%     BDseqscore.unitIDshuffle.weightedDistCorr.zscore(pbe)       = (abs(weightedDistCorr(pbe)) - mean(abs(weightedDistCorr_null(pbe, :, 2)))) / std(abs(weightedDistCorr_null(pbe, :, 2)));
%     
%     
%     if mod(pbe, 50) == 0
%        fprintf('\n event #%d: pscore wc %.1f wdc %.1f  and zscore wc %.2f wdc %.2f', ...
%            pbe, BDseqscore.unitIDshuffle.weightedCorr.prctilescore(pbe, 1), BDseqscore.unitIDshuffle.weightedDistCorr.prctilescore(pbe, 1), BDseqscore.unitIDshuffle.weightedCorr.zscore(pbe, 1), BDseqscore.unitIDshuffle.weightedDistCorr.zscore(pbe, 1))
%     end
%    
% end
% 
% 



end
