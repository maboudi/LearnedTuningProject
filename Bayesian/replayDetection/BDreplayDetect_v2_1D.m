function [posteriorProbMatrix, weightedCorr, jumpDist, BDseqscore, sigMatrix, coveredTrackLen] = BDreplayDetect_v2_1D(eventsBinnedfiring, tuning, activeUnits, binDur, FileBase, dataType)

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


nPBEs        = size(eventsBinnedfiring, 1);
noPositionBins = size(tuning, 2);

% noshufflingMethods = 2; %% number of methods for shuffling (time swap and column cycle shuffle here)
nShuffles      = 500; %% number of shuffles in each of the methods 


if isempty(activeUnits)
   activeUnits = 1:size(tuning, 1); 
end

%%%%%%%%%%%% Initiation %%%%%%%%%%%%%


%%% posterior probability matrices
%%% the matrix comprised of posterior probability of the track positions (as rows) within each time bins of the event (as columns)


posteriorProbMatrix = cell(nPBEs,1); 

%%% 1st array -->  leftward(RL) trajectory only
%%% 2nd array --> rightward(LR) trajectory only
%%% 3rd array -->  summed over the two directions
%%% we can use either the sum or the two trajectories separately for replay
%%% detection (here we have used just the sum)


weightedCorr    = zeros(nPBEs, 1); 
jumpDist        = struct('maxJump', cell(nPBEs, 1), 'normMaxJump', cell(nPBEs, 1), 'medianJump', cell(nPBEs, 1));
coveredTrackLen = zeros(nPBEs, 1);

BDseqscore      = struct('prctilescore', zeros(nPBEs, 1), 'zscore', zeros(nPBEs, 1)); % significance: data as a percentile of shuffle distribution or zscore


for pbe =  1: nPBEs
    
    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); %%% using just the active units
    
    % exclude the time bins with no unit firing
    
    sumFiring = sum(currEvent, 1);
    idx       = find(sumFiring == 0);
    
    
    %%% Bayesian Sequential Score

    %%% calculate the posterior probabilty matrix without any normalization 
    
    postprob = baysDecoder(currEvent, tuning(activeUnits, :), binDur)+ eps*eps; 
    
    
    %%% normalize posterior probabilty matrices of each direction, in case we are interested in
    %%% doing replay evaluation for each direction independently (as done in Grosmark&Buzsaki, Science2015)
        
    posteriorProbMatrix{pbe} = postprob ./ repmat(sum(postprob, 1), [noPositionBins, 1]); %% normalized posterior probabilty matrix for each direction
 
    posteriorProbMatrix{pbe}(:, idx) = eps*eps; 
 
    
    % replay evaluation

    [weightedCorr(pbe), pbeJumpDist, coveredTrackLen(pbe)] = replayevaluation_v2(posteriorProbMatrix{pbe});

    % concatenate jump distance for the two directions
    f = fieldnames(jumpDist);
    for ii = 1: length(f)
        jumpDist(pbe).(f{ii}) = pbeJumpDist(ii);
    end
    
    
end



% calculating shuffle distributions

weightedCorr_null   = zeros(nPBEs,  nShuffles);
jumpDist_null       = struct('maxJump', cell(nPBEs, 1), 'normMaxJump', cell(nPBEs, 1), 'medianJump', cell(nPBEs, 1)); 


for pbe = 1: nPBEs
    
    
    currEvent  = eventsBinnedfiring{pbe, 2}(activeUnits, :); 
    nTimeBins = size(currEvent, 2);
    
    currPPM = posteriorProbMatrix(pbe);
    
    weightedCorr_null_pbe = zeros(nShuffles, 1); % for individual PBEs
    
    jumpDist_null_pbe     = zeros(nShuffles, 3);
    
    
    
    % method one: column cycle shuffle
    parfor sn = 1 : nShuffles 
        
        randShift = randi(noPositionBins-1, 1, nTimeBins);

        posteriorProbMatrix_shuffle = column_cycle_shuffle(currPPM, randShift);
                
        [weightedCorr_null_pbe(sn), jumpDist_null_pbe(sn, :)] = replayevaluation_v2(posteriorProbMatrix_shuffle);    
            
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
    
    weightedCorr_null(pbe, :) = weightedCorr_null_pbe; % concatenating two directions and permuting the dimensions 
    
    f = fieldnames(jumpDist_null);
    for ii = 1: length(f)
        jumpDist_null(pbe).(f{ii}) = jumpDist_null_pbe(:, ii);
    end
    

    % BD sequence score for each PBE calcualted as a percentile of the shuffle distribution
    
    for dir = 1:2
        
        % percentile score
        BDseqscore.prctilescore(pbe) = length(find(abs(weightedCorr_null(pbe, :)) < abs(weightedCorr(pbe))))/nShuffles * 100;
        
        % z score
        BDseqscore.zscore(pbe)       = (abs(weightedCorr(pbe)) - mean(abs(weightedCorr_null(pbe, :)))) / std(abs(weightedCorr_null(pbe, :)));
        
    end


    if mod(pbe, 50) == 0
        fprintf('replay scores of event #%d are pscore %.1f  and zscore %.2f \n', pbe, BDseqscore.prctilescore(pbe), BDseqscore.zscore(pbe))
    end
    
end

%% significance matrix: two features together, absolute weighted correlation and jump distance



%%% jump distance

maxJumpP     = zeros(nPBEs, 1);
normMaxJumpP = zeros(nPBEs, 1);
medianJumpP  = zeros(nPBEs, 1);

for pbe = 1:nPBEs
    
    maxJumpP(pbe)     = jumpDist(pbe).maxJump;
    normMaxJumpP(pbe) = jumpDist(pbe).normMaxJump;
    medianJumpP(pbe)  = jumpDist(pbe).medianJump;

end
    

%%% surrogate distributions

maxJump_surrP     = zeros(nPBEs, nShuffles);
normMaxJump_surrP = zeros(nPBEs, nShuffles);
medianJump_SurrP  = zeros(nPBEs, nShuffles);

for sn = 1:nShuffles
    
   for pbe = 1: nPBEs
       
      maxJump_surrP(pbe, sn)     = jumpDist_null(pbe).maxJump(sn);
      normMaxJump_surrP(pbe, sn) = jumpDist_null(pbe).normMaxJump(sn);
      medianJump_SurrP(pbe, sn)  = jumpDist_null(pbe).medianJump(sn);
      
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
twoDMatrix1 = calWCorrJdistMatrix(weightedCorr, maxJumpP, wcthresholds, jdthresholds);

% using normalized max jump distance
twoDMatrix2 = calWCorrJdistMatrix(weightedCorr, normMaxJumpP, wcthresholds, jdthresholds);

% using median jump distance
twoDMatrix3 = calWCorrJdistMatrix(weightedCorr, medianJumpP, wcthresholds, jdthresholds);



% SURROGATES

twoDMatrix_surr1 = zeros(nwcthreshs, njdthreshs, nShuffles); % max jump distance
twoDMatrix_surr2 = zeros(nwcthreshs, njdthreshs, nShuffles); % normalized max jump distance
twoDMatrix_surr3 = zeros(nwcthreshs, njdthreshs, nShuffles); % median Jump distance

for sn = 1: nShuffles
    
    % max jump distance
    twoDMatrix_surr1(:, :, sn) = calWCorrJdistMatrix(weightedCorr_null(:, sn), maxJump_surrP(:, sn), wcthresholds, jdthresholds);
    
    % normalized max jump distance
    twoDMatrix_surr2(:, :, sn) = calWCorrJdistMatrix(weightedCorr_null(:, sn), normMaxJump_surrP(:, sn), wcthresholds, jdthresholds);
    
    % median Jump distance
    twoDMatrix_surr3(:, :, sn) = calWCorrJdistMatrix(weightedCorr_null(:, sn), medianJump_SurrP(:, sn), wcthresholds, jdthresholds);
    
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



%%
subfolder = fullfile(FileBase, dataType);
mkdir(subfolder)


save(fullfile(subfolder, 'BDoutputs'),  'posteriorProbMatrix', 'weightedCorr', 'jumpDist', 'coveredTrackLen', 'BDseqscore', 'sigMatrix') 

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
    
