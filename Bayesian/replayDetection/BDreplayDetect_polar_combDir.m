function [posteriorProbMatrix, replayScore, weightedCorr, jumpDist, BDseqscore, sigMatrix, onFitLineElements, begPosition, endPosition, replayScore_null, weightedCorr_null] = BDreplayDetect_polar_combDir(eventsBinnedfiring, tuningRL, tuningLR, activeUnits, binDur)



% This function is intended to calculate the bayesian replay score for each
% event. The significance of the replay scores is evaluated through
% comparing it with a distribution of scores resulted from shuffle
% datasets, i.e. column cycle shuffle, time bin shuffle, etc. 


%%% for now we skip the replay type-we get back to it later ..


% first replace all of zeros in tuning with a small value (at the level of smallest non-zero)

% tuningRL(~tuningRL) = 1e-4;
% tuningLR(~tuningLR) = 1e-4;


nEvents       = size(eventsBinnedfiring, 1);
nPositionBins = size(tuningRL, 2);

% noshufflingMethods = 2; %% number of methods for shuffling (time swap and column cycle shuffle here)
nShuffles = 200; %% number of shuffles in each of the methods 


if isempty(activeUnits)
   activeUnits = 1:size(tuningRL, 1); 
end

%%%%%%%%%%%% Initiation %%%%%%%%%%%%%


%%% posterior probability matrices
%%% the matrix comprised of posterior probability of the track positions (as rows) within each time bins of the event (as columns)


posteriorProbMatrix = cell(nEvents, 1); 

%%% 1st array -->  leftward(RL) trajectory only
%%% 2nd array --> rightward(LR) trajectory only
%%% 3rd array -->  summed over the two directions
%%% we can use either the sum or the two trajectories separately for replay
%%% detection (here we have used just the sum)



% %%% best fitting line
% slope       = zeros(nEvents, 3); %% RL, LR, and both together (?)
% rho         = zeros(nEvents,3);  

begPosition = zeros(nEvents,2);   %% y intercept; the third dimension is comprised of RL,LR, and both, respectively
endPosition = zeros(nEvents,2);  

%%% replay score by Bayesian (goodness of fit of the best fitting line)

replayScore  = zeros(nEvents, 1);
slope        = zeros(nEvents, 1);

weightedCorr = zeros(nEvents, 1);

onFitLineElements = cell(nEvents, 1);

jumpDist    = struct('maxJump', cell(nEvents, 1), 'normMaxJump', cell(nEvents, 1), 'medianJump', cell(nEvents, 1));

BDseqscore  = struct('weightedCorr', [], 'replayScore', []);
BDseqscore.weightedCorr = struct('prctilescore', zeros(nEvents, 1), 'zscore', zeros(nEvents, 1));
BDseqscore.replayScore  = struct('prctilescore', zeros(nEvents, 1), 'zscore', zeros(nEvents, 1));


% %%% replay type
% replaydirection = zeros(nEvents,1); %% if above zero, current event is replaying rightward (LR) trajectory template and if less than zero, leftward (RL) template
% replayOrder = zeros(nEvents, 1); %% if above zero is a forward replay and if less than zero is a reverse replay
% 


fprintf('\n Calculating replay scores of actuall PBEs')

tic

% fitLineDictionary = generateLineDictionary(nPositionBins, 1:30);

for pbe =  1: nEvents


    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); %%% using just the active units
    nTimebins = size(currEvent, 2);
    % exclude the time bins with no unit firing
    
    sumFiring = sum(currEvent, 1);
    idx       = find(sumFiring == 0);
    
    %%% Bayesian Sequential Score

    %%% calculate the posterior probabilty matrix without any normalization 
    
    postprobRL = baysDecoder(currEvent, tuningRL(activeUnits, :), binDur); 
    postprobLR = baysDecoder(currEvent, tuningLR(activeUnits, :), binDur);

 
    %%% Summed probability distribution, every column is normalized to one
    %%% (this is the same method as used in Davidson et al. 2009)
    
    posteriorProbMatrix{pbe} = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [nPositionBins, 1]);
    posteriorProbMatrix{pbe}(:, idx) = 0;

    %%% calculate the best fitting line and the replay score and replay
    %%% type for the current event
    

    fitLineDictionary = generateLineDictionary(nPositionBins, nTimebins);
    [replayScore(pbe), weightedCorr(pbe), pbeJumpDist, onFitLineElements{pbe}, begPosition(pbe,:), endPosition(pbe,:)] = replayevaluation_polar2(posteriorProbMatrix{pbe, 1}, fitLineDictionary);
    slope(pbe) = (endPosition(pbe, 2) - begPosition(pbe, 2))/(endPosition(pbe, 1) - begPosition(pbe, 1));

   
    
    % concatenate jump distance for the two directions
    f = fieldnames(jumpDist);
    for ii = 1: length(f)
%         jumpDist(pbe).(f{ii}) = [pbeJumpDistRL(ii) pbeJumpDistLR(ii)];
          jumpDist(pbe).(f{ii}) = pbeJumpDist(ii);
    end
    
    if mod(pbe, 50) == 0
        fprintf('.')
    end

end

toc



%% calculating shuffle distributions (within-PBE time swap)

replayScore_null  = zeros(nEvents, nShuffles);
slope_null        = zeros(nEvents, nShuffles);

weightedCorr_null = zeros(nEvents, nShuffles);

jumpDist_null    = struct('maxJump', cell(nEvents, 1), 'normMaxJump', cell(nEvents, 1), 'medianJump', cell(nEvents, 1)); 


fprintf('\n Measuring statistical significances using shuffled PBEs')
tic
for pbe = 1: nEvents
   
    
    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); 
    nTimeBins = size(currEvent, 2);
    
    fitLineDictionary = generateLineDictionary(nPositionBins, nTimeBins);
    
    currPosteriorProbMatrix = posteriorProbMatrix(pbe);
%     posteriorProbMatrixLR = posteriorProbMatrix(pbe, 2);
    

    replayScore_null_pbe = zeros(nShuffles, 1); % for individual PBEs
    slope_null_pbe       = zeros(nShuffles, 1);
    
    weightedCorr_null_pbe = zeros(nShuffles, 1);

    jumpDist_nullRL_pbe    = zeros(nShuffles, 3);

    parfor sn = 1: nShuffles
        
        posteriorProbMatrix_shuffle = genTimeSwap(currPosteriorProbMatrix);
        posteriorProbMatrix_shuffle = posteriorProbMatrix_shuffle{1};

        [replayScore_null_pbe(sn), weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :), begPosition_null, endPosition_null] = replayevaluation_polar2(posteriorProbMatrix_shuffle, fitLineDictionary);
        slope_null_pbe(sn) = (endPosition_null(2) - begPosition_null(2))/(endPosition_null(1) - begPosition_null(1));

    end
    
    
    replayScore_null(pbe, :)  = permute(replayScore_null_pbe, [2, 1]); % concatenating two directions and permuting the dimensions 
    slope_null(pbe, :)        = permute(sloep_null_pbe, [2, 1]);
    
    weightedCorr_null(pbe, :) = permute(weightedCorr_null_pbe, [2, 1]);
    
    
    f = fieldnames(jumpDist_null);
    for ii = 1: length(f)
        jumpDist_null(pbe).(f{ii}) = jumpDist_null_pbe(:, ii);
    end


    % sequence scores
    
        
        BDseqscore.replayScore.prctilescore(pbe)  = length(find(replayScore_null(pbe, :) < replayScore(pbe)))/nShuffles*100;
        BDseqscore.replayScore.zscore(pbe)        = (replayScore(pbe) - mean(replayScore_null(pbe, :)))/std(replayScore_null(pbe, :));
        
        BDseqscore.weightedCorr.prctilescore(pbe) = length(find(abs(weightedCorr_null(pbe, :)) < abs(weightedCorr(pbe))))/nShuffles * 100;
        BDseqscore.weightedCorr.zscore(pbe)       = (abs(weightedCorr(pbe)) - mean(abs(weightedCorr_null(pbe, :)))) / std(abs(weightedCorr_null(pbe, :)));

    
    
    if mod(pbe, 50) == 0
%        fprintf('\n event #%d: pscore rs %.1f wc %.1f  and zscore rs %.2f wc %.2f', ...
%            pbe, max(BDseqscore.replayScore.prctilescore(pbe, :)), max(BDseqscore.weightedCorr.prctilescore(pbe, :)), max(BDseqscore.replayScore.zscore(pbe, :)), max(BDseqscore.weightedCorr.zscore(pbe, :)))
       fprintf('\n event #%d: pscore rs %.1f wc %.1f  and zscore rs %.2f wc %.2f', ...
           pbe, BDseqscore.replayScore.prctilescore(pbe, 1), BDseqscore.weightedCorr.prctilescore(pbe, 1), BDseqscore.replayScore.zscore(pbe, 1), BDseqscore.weightedCorr.zscore(pbe, 1))
    end
    
end

toc

%% significance matrix1: absolute weighted correlation and jump distance


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
  maxJump_surr(pbe, :)     = jumpDist_null(pbe).maxJump;
  normMaxJump_surr(pbe, :) = jumpDist_null(pbe).normMaxJump;
  medianJump_Surr(pbe, :)  = jumpDist_null(pbe).medianJump;
end




%%% significance matrix

wcthresholds = 0:0.1:0.9;
jdthresholds = 0.1:0.1:1;

nwcthreshs   = length(wcthresholds);
njdthreshs   = length(jdthresholds);

% calculating the matrices of PBE counts passing the two thresholds


% DATA

% using max jump distance
twoDMatrix1 = calRScoreJdistMatrix(weightedCorr, maxJumpP, wcthresholds, jdthresholds);

% using normalized max jump distance
twoDMatrix2 = calRScoreJdistMatrix(weightedCorr, normMaxJumpP, wcthresholds, jdthresholds);

% using median jump distance
twoDMatrix3 = calRScoreJdistMatrix(weightedCorr, medianJumpP, wcthresholds, jdthresholds);



% SURROGATES

twoDMatrix_surr1 = zeros(nwcthreshs, njdthreshs, nShuffles); % max jump distance
twoDMatrix_surr2 = zeros(nwcthreshs, njdthreshs, nShuffles); % normalized max jump distance
twoDMatrix_surr3 = zeros(nwcthreshs, njdthreshs, nShuffles); % median Jump distance

for sn = 1: nShuffles
    
    % max jump distance
    twoDMatrix_surr1(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null(:, sn), maxJump_surr(:, sn), wcthresholds, jdthresholds);
    
    % normalized max jump distance
    twoDMatrix_surr2(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null(:, sn), normMaxJump_surr(:, sn), wcthresholds, jdthresholds);
    
    % median Jump distance
    twoDMatrix_surr3(:, :, sn) = calRScoreJdistMatrix(weightedCorr_null(:, sn), medianJump_Surr(:, sn), wcthresholds, jdthresholds);
    
end


% SIGNIFICANCE MATRIX

sigMatrix = struct('wc_maxJump', zeros(nwcthreshs, njdthreshs), 'wc_normMaxJump', zeros(nwcthreshs, njdthreshs), 'wc_medianJump', zeros(nwcthreshs, njdthreshs), 'rs_slp', []);


for ii = 1: nwcthreshs
   for jj = 1: njdthreshs
      
       nSurrGreater                     = length(find(twoDMatrix_surr1(ii, jj, :) >= twoDMatrix1(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
       sigMatrix.wc_maxJump(ii, jj)     = (nSurrGreater+1)/(nShuffles+1); % Monte Carlo P-value
       
       nSurrGreater                     = length(find(twoDMatrix_surr2(ii, jj, :) >= twoDMatrix2(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
       sigMatrix.wc_normMaxJump(ii, jj) = (nSurrGreater+1)/(nShuffles+1);
       
       nSurrGreater                     = length(find(twoDMatrix_surr3(ii, jj, :) >= twoDMatrix3(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
       sigMatrix.wc_medianJump(ii, jj)  = (nSurrGreater+1)/(nShuffles+1);
       
       
       
   end
end


%% significance matrix2: replay score and slope


% find the ranges (the replay score and slope are not within 0 to 1 range)

replayScorePooled = [replayScore replayScore_null];
slopePooled       = [slope slope_null];


%%% significance matrix

rsthresholds = linspace(min(replayScorePooled(:)), max(replayScorePooled(:)), 10);
slpthresholds = linspace(min(slopePooled(:)), max(slopePooled(:)), 10);

nrsthreshs   = length(rsthresholds);
nslpthreshs   = length(slpthresholds);


twoDMatrix = calRScoreJdistMatrix(replayScore, slope, rsthresholds, slpthresholds);


% SURROGATES

twoDMatrix_surr = zeros(nrsthreshs, nslpthreshs, nShuffles); % max jump distance

for sn = 1: nShuffles   
    twoDMatrix_surr(:, :, sn) = calRScoreJdistMatrix(replayScore_null(:, sn), slope_null(:, sn), rsthresholds, slpthresholds);
end


sigMatrix.rs_slp = zeros(nrsthreshs, nslpthreshs);
for ii = 1: nrsthreshs
   for jj = 1: nslpthreshs
      
       nSurrGreater             = length(find(twoDMatrix_surr(ii, jj, :) >= twoDMatrix(ii, jj))); % number of surrogates with number of PBEs meeting both critera greater than actual dataset
       sigMatrix.rs_slp(ii, jj) = (nSurrGreater+1)/(nShuffles+1); % Monte Carlo P-value

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
        
        count(ii, jj) = length(find(abs(replayScore > wcthresh & jumpDistance < jdthresh))); % number of PBEs meeting passing the two thresholds
        
    end
end

end
    
