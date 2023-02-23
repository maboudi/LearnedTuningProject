function [posteriorProbMatrix, replayScore, weightedCorr, jumpDist, BDseqscore, sigMatrix, onlyLineElements, begPosition, endPosition, replayScore_null, weightedCorr_null, coveredLen, coveredLen_null, jumpDist_null] = BDreplayDetect_1D_vFeb20(eventsBinnedfiring, tuning, activeUnits, binDur)


% This function is intended to calculate the bayesian replay score for each
% event. The significance of the replay scores is evaluated through
% comparing it with a distribution of scores resulted from shuffle
% datasets, i.e. column cycle shuffle, time bin shuffle, etc. 


%%% for now we skip the replay type-we get back to it later ..


% first replace all of zeros in tuning with a small value (at the level of smallest non-zero)

tuning(~tuning) = 1e-4;


nEvents       = size(eventsBinnedfiring, 1);
nPositionBins = size(tuning, 2);

% noshufflingMethods = 2; %% number of methods for shuffling (time swap and column cycle shuffle here)
nShuffles = 500; %% number of shuffles in each of the methods 


if isempty(activeUnits)
   activeUnits = 1:size(tuning, 1); 
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



% %%% best fitting line
% coveredLen       = zeros(nEvents, 3); %% RL, LR, and both together (?)
% rho         = zeros(nEvents,3);  

begPosition = zeros(nEvents,2);   %% y intercept; the third dimention is comprised of RL,LR, and both, respectively
endPosition = zeros(nEvents,2);  


%%% replay score by Bayesian (goodness of fit of the best fitting line)

replayScore  = zeros(nEvents, 1);
coveredLen   = zeros(nEvents, 1);

weightedCorr = zeros(nEvents, 1);

onlyLineElements = cell(nEvents, 1);


temp = arrayfun(@(x) struct('maxJump', [], 'normMaxJump', [], 'medianJump', []), 1:nEvents, 'UniformOutput',0);
jumpDist = horzcat(temp{:});% looks dumb but can be extended to two templates (directions)
clear temp


baseStruct  = struct('weightedCorr', [], 'replayScore', []);
baseStruct.weightedCorr = struct('prctilescore', zeros(nEvents, 1), 'zscore', zeros(nEvents, 1));
baseStruct.replayScore  = struct('prctilescore', zeros(nEvents, 1), 'zscore', zeros(nEvents, 1));

BDseqscore  = struct('wPBEtimeswap', baseStruct, 'unitIDshuffle', baseStruct);
% %%% replay type
% replaydirection = zeros(nEvents,1); %% if above zero, current event is replaying rightward (LR) trajectory template and if less than zero, leftward (RL) template
% replayOrder = zeros(nEvents, 1); %% if above zero is a forward replay and if less than zero is a reverse replay
% 


fprintf('\n Calculating replay scores of actuall PBEs')

tic

% lineDictionary = generateLineDictionary(nPositionBins, 1:30);

for pbe =  1: nEvents


    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); %%% using just the active units
    nTimebins = size(currEvent, 2);
    % exclude the time bins with no unit firing
    
    sumFiring = sum(currEvent, 1);
    idx       = find(sumFiring == 0);
    
    
    
    %%% Bayesian Sequential Score

    %%% calculate the posterior probabilty matrix without any normalization 
    
    postprob = baysDecoder(currEvent, tuning(activeUnits, :), binDur); 

  
    
    posteriorProbMatrix{pbe, 1} = postprob./repmat(sum(postprob, 1), [nPositionBins, 1]);
    posteriorProbMatrix{pbe, 1}(:, idx) = 0;

    %%% calculate the best fitting line and the replay score and replay
    %%% type for the current event
   
    fitLineDictionary = generateLineDictionary(nPositionBins, nTimebins);
    [replayScore(pbe), weightedCorr(pbe), pbeJumpDistRL, onlyLineElements{pbe}, begPosition(pbe,:), endPosition(pbe,:)] = replayevaluation_polar2(posteriorProbMatrix{pbe, 1}, fitLineDictionary);

    if endPosition(pbe, 2) >= begPosition(pbe, 2) % positive slope
        coveredLen(pbe) = abs(min(endPosition(pbe, 2), nPositionBins) - max(begPosition(pbe, 2), 1))/nPositionBins;
    else
        coveredLen(pbe) = abs(min(begPosition(pbe, 2), nPositionBins) - min(endPosition(pbe, 2), 1))/nPositionBins;
    end 
    
    % concatenate jump distance for the two directions
    f = fieldnames(jumpDist);
    for ii = 1: length(f)
%         jumpDist(pbe).(f{ii}) = [pbeJumpDistRL(ii) pbeJumpDistLR(ii)];
          jumpDist(pbe).(f{ii}) = pbeJumpDistRL(ii);
    end
    
    if mod(pbe, 50) == 0
        fprintf('.')
    end

end

toc


%% calculating shuffle distributions

replayScore_null  = zeros(nEvents, nShuffles, 2); % third dimension corresponds to different methods of shuffle
coveredLen_null   = zeros(nEvents, nShuffles, 2);

weightedCorr_null = zeros(nEvents, nShuffles, 2);


temp = arrayfun(@(x) struct('maxJump', zeros(nShuffles, 2), 'normMaxJump', zeros(nShuffles, 2), 'medianJump', zeros(nShuffles, 2)), 1:nEvents, 'UniformOutput',0);
jumpDist_null = horzcat(temp{:}); 
clear temp


%% 1) within-PBE time swap

fprintf('\n Measuring statistical significances by comparing against within-PBE time swap')


for pbe = 1: nEvents
   
    
    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); 
    nTimeBins = size(currEvent, 2);
    
    fitLineDictionary = generateLineDictionary(nPositionBins, nTimeBins);
    
    currPosteriorProbMatrix = posteriorProbMatrix(pbe);
%     posteriorProbMatrixLR = posteriorProbMatrix(pbe, 2);
    

    replayScore_null_pbe = zeros(nShuffles, 1); % for individual PBEs
    coveredLen_null_pbe  = zeros(nShuffles, 1);
    
    weightedCorr_null_pbe = zeros(nShuffles, 1);

    jumpDist_null_pbe    = zeros(nShuffles, 3); % each column corresponds to one way of calculating the jump distance

    parfor sn = 1: nShuffles
        
        posteriorProbMatrix_shuffle = genTimeSwap(currPosteriorProbMatrix);
        posteriorProbMatrix_shuffle = posteriorProbMatrix_shuffle{1};

        [replayScore_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_shuffle, fitLineDictionary);
        
        if endPosition_null(2) >= begPosition_null(2) % positive slope
            coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPositionBins) - max(begPosition_null(2), 1))/nPositionBins;
        else
            coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPositionBins) - min(endPosition_null(2), 1))/nPositionBins;
        end
    end
    
    
    replayScore_null(pbe, :, 1)  = permute(replayScore_null_pbe, [2, 1]); % concatenating two directions and permuting the dimensions 
    coveredLen_null(pbe, :, 1)   = permute(coveredLen_null_pbe, [2, 1]);
    
    weightedCorr_null(pbe, :, 1) = permute(weightedCorr_null_pbe, [2, 1]);
    
    
    f = fieldnames(jumpDist_null);
    for ii = 1: length(f)
        jumpDist_null(pbe).(f{ii})(:, 1) = jumpDist_null_pbe(:, ii);
    end


    % sequence scores
    
    BDseqscore.wPBEtimeswap.replayScore.prctilescore(pbe)  = length(find(replayScore_null(pbe, :, 1) < replayScore(pbe)))/nShuffles*100;
    BDseqscore.wPBEtimeswap.replayScore.zscore(pbe)        = (replayScore(pbe) - nanmean(replayScore_null(pbe, :, 1)))/std(replayScore_null(pbe, :, 1));

    BDseqscore.wPBEtimeswap.weightedCorr.prctilescore(pbe) = length(find(abs(weightedCorr_null(pbe, :, 1)) < abs(weightedCorr(pbe))))/nShuffles * 100;
    BDseqscore.wPBEtimeswap.weightedCorr.zscore(pbe)       = (abs(weightedCorr(pbe)) - nanmean(abs(weightedCorr_null(pbe, :, 1)))) / std(abs(weightedCorr_null(pbe, :, 1)));

    
    
    if mod(pbe, 50) == 0
       fprintf('\n event #%d - %d: pscore rs %.1f wc %.1f  and zscore rs %.2f wc %.2f', ...
           pbe-49, pbe, nanmean(BDseqscore.wPBEtimeswap.replayScore.prctilescore(pbe-49:pbe, 1)), nanmean(BDseqscore.wPBEtimeswap.weightedCorr.prctilescore(pbe-49:pbe, 1)), nanmean(BDseqscore.wPBEtimeswap.replayScore.zscore(pbe-49:pbe, 1)), nanmean(BDseqscore.wPBEtimeswap.weightedCorr.zscore(pbe-49:pbe, 1)))
    end
    
end


%% 2) unit ID shuffle
% There is much redundancy between the procedures for the two shuffle
% methods, can I make it one?

fprintf('\n Measuring statistical significances by comparing against unit ID shuffle')

for pbe = 1:nEvents
    
    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); %%% using just the active units
    nTimebins = size(currEvent, 2);
    
    pbeActUnits = find(sum(currEvent, 2));
    
    
    % calculate new spatial tunings with shuffled cell IDs
        
    tuningshuffles = tuningsUnitIDshuffle_normalize(tuning, [], pbeActUnits, nShuffles, 'same-peak'); % the shuffle type can be 'regular' or 'same-peak'
    
    sumFiring = sum(currEvent, 1);
    idx       = find(sumFiring == 0);
    
    replayScore_null_pbe = zeros(nShuffles, 1); % for individual PBEs
    coveredLen_null_pbe  = zeros(nShuffles, 1);
    
    weightedCorr_null_pbe = zeros(nShuffles, 1);
    jumpDist_null_pbe     = zeros(nShuffles, 3);
    
%     [m n] = size(posteriorProbMatrix{1000});
%     temp = zeros(m,n,50);
%     posteriorProbMatrix_PFshuffle = cell(nShuffles, 1);
    parfor sn = 1:nShuffles % shuffle index

        %%% calculate the posterior probabilty matrix without any normalization 

        postprob = baysDecoder(currEvent, tuningshuffles(activeUnits, :, sn), binDur); 

        posteriorProbMatrix_PFshuffle = postprob./repmat(sum(postprob, 1), [nPositionBins, 1]);
        posteriorProbMatrix_PFshuffle(:, idx) = 0;

        
        fitLineDictionary = generateLineDictionary(nPositionBins, nTimebins);
        [replayScore_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), ~, begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_PFshuffle,fitLineDictionary);
        
        if endPosition_null(2) >= begPosition_null(2) % positive slope
            coveredLen_null_pbe(sn) = abs(min(endPosition_null(2), nPositionBins) - max(begPosition_null(2), 1))/nPositionBins;
        else
            coveredLen_null_pbe(sn) = abs(min(begPosition_null(2), nPositionBins) - min(endPosition_null(2), 1))/nPositionBins;
        end
    end
    
    
    replayScore_null(pbe, :, 2)  = permute(replayScore_null_pbe, [2, 1]); % concatenating two directions and permuting the dimensions 
    coveredLen_null(pbe, :, 2)   = permute(coveredLen_null_pbe, [2, 1]);
    
    weightedCorr_null(pbe, :, 2) = permute(weightedCorr_null_pbe, [2, 1]);
    
    
    f = fieldnames(jumpDist_null);
    for ii = 1: length(f)
        jumpDist_null(pbe).(f{ii})(:, 2) = jumpDist_null_pbe(:, ii);
    end


    % sequence scores
    
    BDseqscore.unitIDshuffle.replayScore.prctilescore(pbe)  = length(find(replayScore_null(pbe, :, 2) < replayScore(pbe)))/nShuffles*100;
    BDseqscore.unitIDshuffle.replayScore.zscore(pbe)        = (replayScore(pbe) - nanmean(replayScore_null(pbe, :, 2)))/std(replayScore_null(pbe, :, 2));

    BDseqscore.unitIDshuffle.weightedCorr.prctilescore(pbe) = length(find(abs(weightedCorr_null(pbe, :, 2)) < abs(weightedCorr(pbe))))/nShuffles * 100;
    BDseqscore.unitIDshuffle.weightedCorr.zscore(pbe)       = (abs(weightedCorr(pbe)) - nanmean(abs(weightedCorr_null(pbe, :, 2)))) / std(abs(weightedCorr_null(pbe, :, 2)));

    
    
    if mod(pbe, 50) == 0
%        fprintf('\n event #%d: pscore rs %.1f wc %.1f  and zscore rs %.2f wc %.2f', ...
%            pbe, max(BDseqscore.replayScore.prctilescore(pbe, :)), max(BDseqscore.weightedCorr.prctilescore(pbe, :)), max(BDseqscore.replayScore.zscore(pbe, :)), max(BDseqscore.weightedCorr.zscore(pbe, :)))
       fprintf('\n event #%d - %d: pscore rs %.1f wc %.1f  and zscore rs %.2f wc %.2f', ...
           pbe-49, pbe, nanmean(BDseqscore.unitIDshuffle.replayScore.prctilescore(pbe-49:pbe, 1)), nanmean(BDseqscore.unitIDshuffle.weightedCorr.prctilescore(pbe-49:pbe, 1)), nanmean(BDseqscore.unitIDshuffle.replayScore.zscore(pbe-49:pbe, 1)), nanmean(BDseqscore.unitIDshuffle.weightedCorr.zscore(pbe-49:pbe, 1)))
    end
   
end



% 
% figure; 
% for ii = 1:100
%     subplot(10 ,10, ii)
%     imagesc(posteriorProbMatrix{ii})
%     colormap('hot'); set(gca, 'YDir', 'normal', 'xtick', [], 'ytick', [])
%     ylabel(['pbe ' num2str(ii)], 'fontsize', 10)
%     title(sprintf('rs_ts = %.0f wc_ts = %.0f\nrs_ui = %.0f wc_ui = %.0f', BDseqscore.wPBEtimeswap.replayScore.prctilescore(ii, 1), BDseqscore.wPBEtimeswap.weightedCorr.prctilescore(ii, 1), ...
%         BDseqscore.unitIDshuffle.replayScore.prctilescore(ii, 1), BDseqscore.unitIDshuffle.weightedCorr.prctilescore(ii, 1)), 'fontweight', 'normal', 'fontsize', 6, 'interpreter', 'none')
%     
%     hold on
%     
%     line([1 size(posteriorProbMatrix{ii}, 2)], [begPosition(ii, 2) endPosition(ii, 2)], 'color', 'c', 'linestyle', ':')
% end
% 
% 
% figure;
% set(gcf, 'position', [100 100 650 950])
% for ii = 1:nShuffles
%    
%     subplot(20,10 ,ii)
%     imagesc(posteriorProbMatrix_PFshuffle{ii})
%     colormap('hot'); set(gca, 'YDir', 'normal', 'xtick', [], 'ytick', [])
%     text(0.5, nPositionBins*1.1, ['shuffle' num2str(ii)], 'fontsize', 4, 'fontweight', 'normal', 'horizontalAlignment', 'left')
% end
% print(gcf, 'pbe56_UIshuffles', '-dpdf')

%% Significance matrices (WC vs JD   and    RS vs SLP)

shuffleMethodNames = {'wPBEtimeswap';'unitIDshuffle'};
sigMatrix = struct('wPBEtimeswap', [], 'unitIDshuffle', []);

for shuffleMethod = 1:2


    %%% absolute weighted correlation and jump distance


    maxJump     = zeros(nEvents, 1);
    normMaxJump = zeros(nEvents, 1);
    medianJump  = zeros(nEvents, 1);

    for pbe = 1:nEvents    
        maxJump(pbe)     = jumpDist(pbe).maxJump;
        normMaxJump(pbe) = jumpDist(pbe).normMaxJump;
        medianJump(pbe)  = jumpDist(pbe).medianJump;
    end

    maxJump_surr     = zeros(nEvents, nShuffles);
    normMaxJump_surr = zeros(nEvents, nShuffles);
    medianJump_Surr  = zeros(nEvents, nShuffles);

    for pbe = 1: nEvents 

      maxJump_surr(pbe, :)     = jumpDist_null(pbe).maxJump(:, shuffleMethod);
      normMaxJump_surr(pbe, :) = jumpDist_null(pbe).normMaxJump(:, shuffleMethod);
      medianJump_Surr(pbe, :)  = jumpDist_null(pbe).medianJump(:, shuffleMethod);

    end

    % calculating the matrices of PBE counts passing pairs of threshold on the two
    % metrics

    wcthresholds = 0:0.1:0.9;
    jdthresholds = 0.1:0.1:1;

    nwcthreshs   = length(wcthresholds);
    njdthreshs   = length(jdthresholds);



    % DATA

    % using max jump distance
    twoDMatrix1 = calRScoreJdistMatrix(weightedCorr, maxJump, wcthresholds, jdthresholds);

    % using normalized max jump distance
    twoDMatrix2 = calRScoreJdistMatrix(weightedCorr, normMaxJump, wcthresholds, jdthresholds);

    % using median jump distance
    twoDMatrix3 = calRScoreJdistMatrix(weightedCorr, medianJump, wcthresholds, jdthresholds);



    % SURROGATES

    twoDMatrix_surr1 = zeros(nwcthreshs, njdthreshs, nShuffles); % max jump distance
    twoDMatrix_surr2 = zeros(nwcthreshs, njdthreshs, nShuffles); % normalized max jump distance
    twoDMatrix_surr3 = zeros(nwcthreshs, njdthreshs, nShuffles); % median Jump distance

    for sn = 1: nShuffles

        % max jump distance
        twoDMatrix_surr1(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null(:, sn, shuffleMethod), maxJump_surr(:, sn), wcthresholds, jdthresholds);

        % normalized max jump distance
        twoDMatrix_surr2(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null(:, sn, shuffleMethod), normMaxJump_surr(:, sn), wcthresholds, jdthresholds);

        % median Jump distance
        twoDMatrix_surr3(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null(:, sn, shuffleMethod), medianJump_Surr(:, sn), wcthresholds, jdthresholds);

    end


    % SIGNIFICANCE MATRIX

    sigMatrix.(shuffleMethodNames{shuffleMethod}) = struct('wc_maxJump', zeros(nwcthreshs, njdthreshs), 'wc_normMaxJump', zeros(nwcthreshs, njdthreshs), 'wc_medianJump', zeros(nwcthreshs, njdthreshs), 'rs_slp', []);


    for ii = 1: nwcthreshs
       for jj = 1: njdthreshs

           nSurrGreater = length(find(twoDMatrix_surr1(ii, jj, :) >= twoDMatrix1(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
           sigMatrix.(shuffleMethodNames{shuffleMethod}).wc_maxJump(ii, jj) = (nSurrGreater+1)/(nShuffles+1); % Monte Carlo P-value

           nSurrGreater = length(find(twoDMatrix_surr2(ii, jj, :) >= twoDMatrix2(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
           sigMatrix.(shuffleMethodNames{shuffleMethod}).wc_normMaxJump(ii, jj) = (nSurrGreater+1)/(nShuffles+1);

           nSurrGreater = length(find(twoDMatrix_surr3(ii, jj, :) >= twoDMatrix3(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
           sigMatrix.(shuffleMethodNames{shuffleMethod}).wc_medianJump(ii, jj) = (nSurrGreater+1)/(nShuffles+1);


       end
    end


    %% significance matrix2: replay score and coveredLen


%     coveredLen(coveredLen > 1) = 1;
%     coveredLen_null(coveredLen_null > 1) = 1;

    %%% significance matrix

    rsthresholds = 0:0.1:0.9;
    coveredLenthresholds =  0:0.1:0.9; 

    nrsthreshs   = length(rsthresholds);
    ncoveredLenthreshs   = length(coveredLenthresholds);


    twoDMatrix = calRScoreJdistMatrix2(replayScore, coveredLen, rsthresholds, coveredLenthresholds);


    % SURROGATES

    twoDMatrix_surr = zeros(nrsthreshs, ncoveredLenthreshs, nShuffles); % max jump distance

    for sn = 1: nShuffles   
        twoDMatrix_surr(:, :, sn) = calRScoreJdistMatrix2(replayScore_null(:, sn, shuffleMethod), coveredLen_null(:, sn, shuffleMethod), rsthresholds, coveredLenthresholds);
    end


    sigMatrix.(shuffleMethodNames{shuffleMethod}).rs_slp = zeros(nrsthreshs, ncoveredLenthreshs);
    for ii = 1: nrsthreshs
       for jj = 1: ncoveredLenthreshs

           nSurrGreater             = length(find(twoDMatrix_surr(ii, jj, :) >= twoDMatrix(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
           sigMatrix.(shuffleMethodNames{shuffleMethod}).rs_slp(ii, jj) = (nSurrGreater+1)/(nShuffles+1); % Monte Carlo P-value

       end
    end



end



end


function count = calRScoreJdistMatrix(replayScore, jumpDistance, wcthresholds, jdthresholds)

nwcthreshs = length(wcthresholds);
njdthreshs = length(jdthresholds);

count = zeros(nwcthreshs, njdthreshs);

for ii = 1:nwcthreshs
    
    wcthresh = wcthresholds(ii);
    
    for jj = 1: njdthreshs
        
        jdthresh = jdthresholds(jj);
        
        count(ii, jj) = length(find(abs(replayScore) > wcthresh & jumpDistance < jdthresh)); % number of PBEs meeting passing the two thresholds
        
    end
end

end
 

function count = calRScoreJdistMatrix2(replayScore, coveredLen, rsthresholds, coveredLenthresholds)

nrsthreshs = length(rsthresholds);
ncoveredLenthreshs = length(coveredLenthresholds);

count = zeros(nrsthreshs, ncoveredLenthreshs);

for ii = 1:nrsthreshs
    
    rsthresh = rsthresholds(ii);
    
    for jj = 1: ncoveredLenthreshs
        
        coveredLenthresh = coveredLenthresholds(jj);
        
        count(ii, jj) = length(find(replayScore > rsthresh & coveredLen > coveredLenthresh)); % number of PBEs meeting passing the two thresholds
        
    end
end

end




