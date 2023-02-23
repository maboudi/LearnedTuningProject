function [posteriorProbMatrix, weightedCorr, jumpDist, BDseqscore, sigMatrix, coveredTrackLen, weightedCorr_null] = BDreplayDetect_v2(eventsBinnedfiring, tuningRL, tuningLR, activeUnits, binDur)

% This function enables evlaution of PBEs using two features, absolute weighted
% correlation and jump distance between decoded positions in consequitive time bins 

% This function returns:

% weightedCorr      : absolute weighted correlation 
% jumpDist          : jump distance (maximum / normalized maximum/ median jump distance)
% BDseqscore        : sequence score -- weighted correlation as a percentile of the
% shuffle distribution (column cycle or within-PBE-time bin shuffle)
% sigMatrix         : significance matrix under parameter space of absolute
% weighted correlation and maximum jump distance (SilvaFosterNatNeurosci2015) 
% coveredTrackLen: covered length of the track by the range of the
% decoded positions with a PBE calculated as ratios of the track length


% first replace all of zeros in tuning with a small value (at the level of smallest non-zero)
tuningRL(~tuningRL) = 1e-4;
tuningLR(~tuningLR) = 1e-4;


nEvents        = size(eventsBinnedfiring, 1);
noPositionBins = size(tuningRL, 2);

% noshufflingMethods = 2; %% number of methods for shuffling (time swap and column cycle shuffle here)
nShuffles      = 500; %% number of shuffles in each of the methods 


if isempty(activeUnits)
   activeUnits = 1:size(tuningRL, 1); 
end

%%%%%%%%%%%% Initiation %%%%%%%%%%%%%


%%% posterior probability matrices
%%% the matrix comprised of posterior probability of the track positions (as rows) within each time bins of the event (as columns)


posteriorProbMatrix = cell(nEvents,3); 

%%% 1st array -->  leftward(RL) trajectory only
%%% 2nd array --> rightward(LR) trajectory only
%%% 3rd array -->  summed over the two directions
%%% we can use either the sum or the two trajectories separately for replay
%%% detection (here we have used just the sum)


weightedCorr    = zeros(nEvents, 2); 
jumpDist        = struct('maxJump', cell(nEvents, 1), 'normMaxJump', cell(nEvents, 1), 'medianJump', cell(nEvents, 1));
coveredTrackLen = zeros(nEvents, 2);

BDseqscore      = struct('prctilescore', zeros(nEvents, 2), 'zscore', zeros(nEvents, 2)); % significance: data as a percentile of shuffle distribution or zscore


for pbe =  1: nEvents
    
    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); %%% using just the active units
    
    % exclude the time bins with no unit firing
    
    sumFiring = sum(currEvent, 1);
    idx       = find(sumFiring == 0);
    
    
    %%% Bayesian Sequential Score

    %%% calculate the posterior probabilty matrix without any normalization 
    
    postprobRL = baysDecoder(currEvent, tuningRL(activeUnits, :), binDur); 
    postprobLR = baysDecoder(currEvent, tuningLR(activeUnits, :), binDur);

    
    postprobRL = postprobRL + eps*eps;
    postprobLR = postprobLR + eps*eps;
    
    
    %%% normalize posterior probabilty matrices of each direction, in case we are interested in
    %%% doing replay evaluation for each direction independently (as done in Grosmark&Buzsaki, Science2015)
        
    posteriorProbMatrix{pbe, 1} = postprobRL ./ repmat(sum(postprobRL, 1), [noPositionBins, 1]); %% normalized posterior probabilty matrix for each direction
    posteriorProbMatrix{pbe, 2} = postprobLR ./ repmat(sum(postprobLR, 1), [noPositionBins, 1]); 

  
    posteriorProbMatrix{pbe, 1}(:, idx) = eps*eps; 
    posteriorProbMatrix{pbe, 2}(:, idx) = eps*eps;

    
    % replay evaluation

    [weightedCorr(pbe, 1), pbeJumpDistRL, coveredTrackLen(pbe, 1)] = replayevaluation_v2(posteriorProbMatrix{pbe, 1});

    [weightedCorr(pbe, 2), pbeJumpDistLR, coveredTrackLen(pbe, 2)] = replayevaluation_v2(posteriorProbMatrix{pbe, 2});
    
    
    % concatenate jump distance for the two directions
    f = fieldnames(jumpDist);
    for ii = 1: length(f)
        jumpDist(pbe).(f{ii}) = [pbeJumpDistRL(ii) pbeJumpDistLR(ii)];
    end
    
    
end



% calculating shuffle distributions

weightedCorr_null   = zeros(nEvents, 2,  nShuffles);
jumpDist_null       = struct('maxJump', cell(nEvents, 1), 'normMaxJump', cell(nEvents, 1), 'medianJump', cell(nEvents, 1)); 


for pbe = 1: nEvents
    
    
    currEvent  = eventsBinnedfiring{pbe, 2}(activeUnits, :); 
    nTimeBins = size(currEvent, 2);
    
    
    posteriorProbMatrixRL   = posteriorProbMatrix(pbe, 1);
    posteriorProbMatrixLR   = posteriorProbMatrix(pbe, 2);
    
    
    weightedCorr_nullRL_pbe = zeros(nShuffles, 1); % for individual PBEs
    weightedCorr_nullLR_pbe = zeros(nShuffles, 1);

    jumpDist_nullRL_pbe     = zeros(nShuffles, 3);
    jumpDist_nullLR_pbe     = zeros(nShuffles, 3);

    
    
    % method one: column cycle shuffle
    parfor sn = 1 : nShuffles 
        
        randShift = randi(noPositionBins-1, 1, nTimeBins);

        posteriorProbMatrixRL_shuffle = column_cycle_shuffle(posteriorProbMatrixRL, randShift);
        posteriorProbMatrixLR_shuffle = column_cycle_shuffle(posteriorProbMatrixLR, randShift);        
        
        
        [weightedCorr_nullRL_pbe(sn), jumpDist_nullRL_pbe(sn, :)] = replayevaluation_v2(posteriorProbMatrixRL_shuffle);    
        [weightedCorr_nullLR_pbe(sn), jumpDist_nullLR_pbe(sn, :)] = replayevaluation_v2(posteriorProbMatrixLR_shuffle); 
            
    end
    
     
%     % method two: time swap
%     parfor sn = 1: nShuffles
%          
% 
%         posteriorProbMatrixRL_shuffle = genTimeSwap(posteriorProbMatrixRL);
%         posteriorProbMatrixRL_shuffle = posteriorProbMatrixRL_shuffle{1};
% 
%         posteriorProbMatrixLR_shuffle = genTimeSwap(posteriorProbMatrixLR);
%         posteriorProbMatrixLR_shuffle = posteriorProbMatrixLR_shuffle{1};
% 
%         [weightedCorr_nullRL_pbe(sn), jumpDist_nullRL_pbe(sn, :)] = replayevaluation_v2(posteriorProbMatrixRL_shuffle);    
%         [weightedCorr_nullLR_pbe(sn), jumpDist_nullLR_pbe(sn, :)] = replayevaluation_v2(posteriorProbMatrixLR_shuffle);
% 
% 
%     end
    
    weightedCorr_null(pbe, :, :) = permute([weightedCorr_nullRL_pbe weightedCorr_nullLR_pbe], [3, 2, 1]); % concatenating two directions and permuting the dimensions 
    
    f = fieldnames(jumpDist_null);
    for ii = 1: length(f)
        jumpDist_null(pbe).(f{ii}) = [jumpDist_nullRL_pbe(:, ii) jumpDist_nullLR_pbe(:, ii)];
    end
    

    % BD sequence score for each PBE calcualted as a percentile of the shuffle distribution
    
    for dir = 1:2
        
        % percentile score
        BDseqscore.prctilescore(pbe, dir) = length(find(abs(weightedCorr_null(pbe, dir, :)) < abs(weightedCorr(pbe, dir))))/nShuffles * 100;
        
        % z score
%         BDseqscore.zscore(pbe, dir)       = (abs(weightedCorr(pbe, dir)) - mean(abs(weightedCorr_null(pbe, dir, :)))) / std(abs(weightedCorr_null(pbe, dir, :)));
        BDseqscore.zscore(pbe, dir)       = abs((weightedCorr(pbe, dir) - mean(weightedCorr_null(pbe, dir, :))) / std(weightedCorr_null(pbe, dir, :)));

        
    end


    if mod(pbe, 100) == 0
        fprintf('replay scores of event #%d are pscore %.1f  and zscore %.2f \n', pbe, max(BDseqscore.prctilescore(pbe, :)), max(BDseqscore.zscore(pbe, :)))
    end
    
end

%% significance matrix: two features together, absolute weighted correlation and jump distance


%%% weighted correlation

weightedCorrP = reshape(weightedCorr, [nEvents*2, 1]); % pooling over the directions


% surrogate distributions

weightedCorr_surrP = reshape(weightedCorr_null, [nEvents*2, nShuffles, 1]); % pooling over the directions


%%% jump distance

maxJumpP     = zeros(nEvents*2, 1);
normMaxJumpP = zeros(nEvents*2, 1);
medianJumpP  = zeros(nEvents*2, 1);

for pbe = 1:nEvents
    
    maxJumpP((pbe-1)*2+1 : pbe*2)     = jumpDist(pbe).maxJump(:);
    normMaxJumpP((pbe-1)*2+1 : pbe*2) = jumpDist(pbe).normMaxJump(:);
    medianJumpP((pbe-1)*2+1 : pbe*2)  = jumpDist(pbe).medianJump(:);

end
    

%%% surrogate distributions

maxJump_surrP     = zeros(nEvents*2, nShuffles);
normMaxJump_surrP = zeros(nEvents*2, nShuffles);
medianJump_SurrP  = zeros(nEvents*2, nShuffles);

for sn = 1:nShuffles
    
   for pbe = 1: nEvents
       
      maxJump_surrP((pbe-1)*2+1 : pbe*2, sn)     = jumpDist_null(pbe).maxJump(sn, :);
      normMaxJump_surrP((pbe-1)*2+1 : pbe*2, sn) = jumpDist_null(pbe).normMaxJump(sn, :);
      medianJump_SurrP((pbe-1)*2+1 : pbe*2, sn)  = jumpDist_null(pbe).medianJump(sn, :);
      
   end
    
end



%%% significance matrix

wcthresholds = 0:0.1:0.9;
jdthresholds = 0.1:0.1:1;

nwcthreshs   = length(wcthresholds);
njdthreshs   = length(jdthresholds);

% calculating the matrices of PBE counts passing the two thresholds


% DATA

% using max jump distance
twoDMatrix1 = calWCorrJdistMatrix(weightedCorrP, maxJumpP, wcthresholds, jdthresholds);

% using normalized max jump distance
twoDMatrix2 = calWCorrJdistMatrix(weightedCorrP, normMaxJumpP, wcthresholds, jdthresholds);

% using median jump distance
twoDMatrix3 = calWCorrJdistMatrix(weightedCorrP, medianJumpP, wcthresholds, jdthresholds);



% SURROGATES

twoDMatrix_surr1 = zeros(nwcthreshs, njdthreshs, nShuffles); % max jump distance
twoDMatrix_surr2 = zeros(nwcthreshs, njdthreshs, nShuffles); % normalized max jump distance
twoDMatrix_surr3 = zeros(nwcthreshs, njdthreshs, nShuffles); % median Jump distance

for sn = 1: nShuffles
    
    % max jump distance
    twoDMatrix_surr1(:, :, sn) = calWCorrJdistMatrix(weightedCorr_surrP(:, sn), maxJump_surrP(:, sn), wcthresholds, jdthresholds);
    
    % normalized max jump distance
    twoDMatrix_surr2(:, :, sn) = calWCorrJdistMatrix(weightedCorr_surrP(:, sn), normMaxJump_surrP(:, sn), wcthresholds, jdthresholds);
    
    % median Jump distance
    twoDMatrix_surr3(:, :, sn) = calWCorrJdistMatrix(weightedCorr_surrP(:, sn), medianJump_SurrP(:, sn), wcthresholds, jdthresholds);
    
end


% SIGNIFICANCE MATRIX

sigMatrix = struct('wc_maxJump', zeros(nwcthreshs, njdthreshs), 'wc_normMaxJump', zeros(nwcthreshs, njdthreshs), 'wc_medianJump', zeros(nwcthreshs, njdthreshs));


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



% %%
% subfolder = fullfile(FileBase, dataType);
% mkdir(subfolder)
% 
% 
% save(fullfile(subfolder, 'BDoutputs'),  'posteriorProbMatrix', 'weightedCorr', 'jumpDist', 'coveredTrackLen', 'BDseqscore', 'sigMatrix') 

end



function count = calWCorrJdistMatrix(weightedCorr, jumpDistance, wcthresholds, jdthresholds)

nwcthreshs = length(wcthresholds);
njdthreshs = length(jdthresholds);

count = zeros(nwcthreshs, njdthreshs);

for ii = 1:nwcthreshs
    
    wcthresh = wcthresholds(ii);
    
    for jj = 1: njdthreshs
        
        jdthresh = jdthresholds(jj);
        
        count(ii, jj) = length(find(abs(weightedCorr > wcthresh & jumpDistance < jdthresh))); % number of PBEs meeting passing the two thresholds
        
    end
end

end
    
