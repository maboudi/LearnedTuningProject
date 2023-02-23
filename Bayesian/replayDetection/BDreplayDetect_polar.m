function [posteriorProbMatrix, replayScore, weightedCorr, jumpDist, BDseqscore, sigMatrix, onlyLineElements, begPosition, endPosition, replayScore_null, weightedCorr_null] = BDreplayDetect_polar(eventsBinnedfiring, tuningRL, tuningLR, activeUnits, binDur)



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
nShuffles = 200; %% number of shuffles in each of the methods 


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



% %%% best fitting line
% slope       = zeros(nEvents, 3); %% RL, LR, and both together (?)
% rho         = zeros(nEvents,3);  

begPosition = zeros(nEvents,2,3);   %% y intercept; the third dimention is comprised of RL,LR, and both, respectively
endPosition = zeros(nEvents,2,3);  

%%% replay score by Bayesian (goodness of fit of the best fitting line)

replayScore  = zeros(nEvents, 2);
weightedCorr = zeros(nEvents, 2);

onlyLineElements = cell(nEvents, 2);

jumpDist    = struct('maxJump', cell(nEvents, 1), 'normMaxJump', cell(nEvents, 1), 'medianJump', cell(nEvents, 1));

BDseqscore  = struct('weightedCorr', [], 'replayScore', []);
BDseqscore.weightedCorr = struct('prctilescore', zeros(nEvents, 2), 'zscore', zeros(nEvents, 2));
BDseqscore.replayScore  = struct('prctilescore', zeros(nEvents, 2), 'zscore', zeros(nEvents, 2));


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
    
    postprobRL = baysDecoder(currEvent, tuningRL(activeUnits, :), binDur); 
    postprobLR = baysDecoder(currEvent, tuningLR(activeUnits, :), binDur);

    
%     postprobRL = postprobRL + eps*eps; % what was this step exactly??
%     postprobLR = postprobLR + eps*eps;
    
    
    %%% Summed probability distribution, every column is normalized to one
    %%% (this is the same method as used in Davidson et al. 2009)
    
%     posteriorProbMatrix{pbe, 3} = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [nPositionBins, 1]);
%     posteriorProbMatrix{pbe, 3}(:, idx) = 0;

    %%% calculate the best fitting line and the replay score and replay
    %%% type for the current event
                       
%     [replayScore(pbe, 3), slope(pbe, 3), rho(pbe, 3), begPosition(pbe,:,3), endPosition(pbe,:,3)] = replayevaluation(posteriorProbMatrix{pbe, 3});
    


    %%% calculating replay order(Forward or reverse) and replay direction
    %%% (LR or RL)    
    
%     [replaydirection(pbe), replayOrder(pbe)] = replaytype(postprobRL, postprobLR, slope(pbe, 3), rho(pbe, 3)); 
    
    
    %%% normalize posterior probabilty matrices of each direction in case we are interested in
    %%% doing replay evaluation for each direction independently (as used by Grosmark&Buzsaki, Science2015)
  
    posteriorProbMatrix{pbe, 1} = postprobRL ./ repmat(sum(postprobRL, 1), [nPositionBins, 1]); %% normalized posterior probabilty matrix for each direction
    posteriorProbMatrix{pbe, 2} = postprobLR ./ repmat(sum(postprobLR, 1), [nPositionBins, 1]); 

  
    posteriorProbMatrix{pbe, 1}(:, idx) = 0; 
    posteriorProbMatrix{pbe, 2}(:, idx) = 0;
    
    
    % replay evaluation
    
%     
%     
%     [replayScore(pbe, 1), weightedCorr, pbeJumpDistRL, begPosition(pbe,:, 1), endPosition(pbe,:, 1)] = replayevaluation_polar(posteriorProbMatrix{pbe, 1});
%     
%     [replayScore(pbe, 2), weightedCorr, pbeJumpDistLR, begPosition(pbe,:, 2), endPosition(pbe,:, 2)] = replayevaluation_polar(posteriorProbMatrix{pbe, 2});
    
    
    lineDictionary = generateLineDictionary(nPositionBins, nTimebins);
    
    [replayScore(pbe, 1), weightedCorr(pbe, 1), pbeJumpDistRL, onlyLineElements{pbe, 1}, begPosition(pbe,:, 1), endPosition(pbe,:, 1)] = replayevaluation_polar2(posteriorProbMatrix{pbe, 1}, lineDictionary);
    
    [replayScore(pbe, 2), weightedCorr(pbe, 2), pbeJumpDistLR, onlyLineElements{pbe, 2}, begPosition(pbe,:, 2), endPosition(pbe,:, 2)] = replayevaluation_polar2(posteriorProbMatrix{pbe, 2}, lineDictionary);
    
    
    % concatenate jump distance for the two directions
    f = fieldnames(jumpDist);
    for ii = 1: length(f)
        jumpDist(pbe).(f{ii}) = [pbeJumpDistRL(ii) pbeJumpDistLR(ii)];
    end
    
    if mod(pbe, 50) == 0
        fprintf('.')
    end

end

toc



% calculating shuffle distributions
%%
replayScore_null  = zeros(nEvents, 2,  nShuffles);
weightedCorr_null = zeros(nEvents, 2,  nShuffles);

jumpDist_null    = struct('maxJump', cell(nEvents, 1), 'normMaxJump', cell(nEvents, 1), 'medianJump', cell(nEvents, 1)); 


fprintf('\n Measuring statistical significances using shuffled PBEs')
tic
for pbe = 1: nEvents
    
    
    
    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); 
    nTimeBins = size(currEvent, 2);
    
    lineDictionary = generateLineDictionary(nPositionBins, nTimeBins);
    
    posteriorProbMatrixRL = posteriorProbMatrix(pbe, 1);
    posteriorProbMatrixLR = posteriorProbMatrix(pbe, 2);
    

    replayScore_nullRL_pbe = zeros(nShuffles, 1); % for individual PBEs
    replayScore_nullLR_pbe = zeros(nShuffles, 1);
    
    weightedCorr_nullRL_pbe = zeros(nShuffles, 1);
    weightedCorr_nullLR_pbe = zeros(nShuffles, 1);

    jumpDist_nullRL_pbe    = zeros(nShuffles, 3);
    jumpDist_nullLR_pbe    = zeros(nShuffles, 3);
    
%     
%     % method one: column cycle shuffle
%     parfor sn = 1 : nShuffles 
%         
%         randShift = randi(nPositionBins-1, 1, nTimeBins);
% 
%         posteriorProbMatrixRL_shuffle = column_cycle_shuffle(posteriorProbMatrixRL, randShift);
%         posteriorProbMatrixLR_shuffle = column_cycle_shuffle(posteriorProbMatrixLR, randShift);        
%         
%         [replayScore_nullRL_pbe(sn), jumpDist_nullRL_pbe(sn, :)] = replayevaluation_polar(posteriorProbMatrixRL_shuffle);
%         [replayScore_nullLR_pbe(sn), jumpDist_nullLR_pbe(sn, :)] = replayevaluation_polar(posteriorProbMatrixLR_shuffle);
%         
%     end
%     
     
    % method two: time swap
    parfor sn = 1: nShuffles
        

        posteriorProbMatrixRL_shuffle = genTimeSwap(posteriorProbMatrixRL);
        posteriorProbMatrixRL_shuffle = posteriorProbMatrixRL_shuffle{1};

        posteriorProbMatrixLR_shuffle = genTimeSwap(posteriorProbMatrixLR);
        posteriorProbMatrixLR_shuffle = posteriorProbMatrixLR_shuffle{1};

        [replayScore_nullRL_pbe(sn), weightedCorr_nullRL_pbe(sn), jumpDist_nullRL_pbe(sn, :)] = replayevaluation_polar2(posteriorProbMatrixRL_shuffle, lineDictionary);
        [replayScore_nullLR_pbe(sn), weightedCorr_nullLR_pbe(sn), jumpDist_nullLR_pbe(sn, :)] = replayevaluation_polar2(posteriorProbMatrixLR_shuffle, lineDictionary);

% 
%         [replayScore_nullRL_pbe(sn), weightedCorr_nullRL_pbe(sn), jumpDist_nullRL_pbe(sn, :)] = replayevaluation_polar(posteriorProbMatrixRL_shuffle);
%         [replayScore_nullLR_pbe(sn), weightedCorr_nullLR_pbe(sn), jumpDist_nullLR_pbe(sn, :)] = replayevaluation_polar(posteriorProbMatrixLR_shuffle);
    end
    
    
    replayScore_null(pbe, :, :)  = permute([replayScore_nullRL_pbe replayScore_nullLR_pbe], [3, 2, 1]); % concatenating two directions and permuting the dimensions 
    weightedCorr_null(pbe, :, :) = permute([weightedCorr_nullRL_pbe weightedCorr_nullLR_pbe], [3, 2, 1]);
    
    f = fieldnames(jumpDist_null);
    for ii = 1: length(f)
        jumpDist_null(pbe).(f{ii}) = [jumpDist_nullRL_pbe(:, ii) jumpDist_nullLR_pbe(:, ii)];
    end


    % sequence scores
    
    
    for dir = 1:2
        

        BDseqscore.replayScore.prctilescore(pbe, dir)  = length(find(replayScore_null(pbe, dir, :) < replayScore(pbe, dir)))/nShuffles*100;
        BDseqscore.replayScore.zscore(pbe, dir)        = (replayScore(pbe, dir) - mean(replayScore_null(pbe, dir, :)))/std(replayScore_null(pbe, dir, :));
        
        BDseqscore.weightedCorr.prctilescore(pbe, dir) = length(find(abs(weightedCorr_null(pbe, dir, :)) < abs(weightedCorr(pbe, dir))))/nShuffles * 100;
        BDseqscore.weightedCorr.zscore(pbe, dir)       = abs((weightedCorr(pbe, dir) - mean(weightedCorr_null(pbe, dir, :))) / std(weightedCorr_null(pbe, dir, :)));

        
    end
    
    
    if mod(pbe, 50) == 0
       fprintf('\n event #%d: pscore rs %.1f wc %.1f  and zscore rs %.2f wc %.2f', ...
           pbe, max(BDseqscore.replayScore.prctilescore(pbe, :)), max(BDseqscore.weightedCorr.prctilescore(pbe, :)), max(BDseqscore.replayScore.zscore(pbe, :)), max(BDseqscore.weightedCorr.zscore(pbe, :)))
    end
    
end

toc

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
twoDMatrix1 = calRScoreJdistMatrix(weightedCorrP, maxJumpP, wcthresholds, jdthresholds);

% using normalized max jump distance
twoDMatrix2 = calRScoreJdistMatrix(weightedCorrP, normMaxJumpP, wcthresholds, jdthresholds);

% using median jump distance
twoDMatrix3 = calRScoreJdistMatrix(weightedCorrP, medianJumpP, wcthresholds, jdthresholds);



% SURROGATES

twoDMatrix_surr1 = zeros(nwcthreshs, njdthreshs, nShuffles); % max jump distance
twoDMatrix_surr2 = zeros(nwcthreshs, njdthreshs, nShuffles); % normalized max jump distance
twoDMatrix_surr3 = zeros(nwcthreshs, njdthreshs, nShuffles); % median Jump distance

for sn = 1: nShuffles
    
    % max jump distance
    twoDMatrix_surr1(:, :, sn) = calRScoreJdistMatrix(weightedCorr_surrP(:, sn), maxJump_surrP(:, sn), wcthresholds, jdthresholds);
    
    % normalized max jump distance
    twoDMatrix_surr2(:, :, sn) = calRScoreJdistMatrix(weightedCorr_surrP(:, sn), normMaxJump_surrP(:, sn), wcthresholds, jdthresholds);
    
    % median Jump distance
    twoDMatrix_surr3(:, :, sn) = calRScoreJdistMatrix(weightedCorr_surrP(:, sn), medianJump_SurrP(:, sn), wcthresholds, jdthresholds);
    
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

% 
% 
% subfolder = fullfile(FileBase, dataType);
% mkdir(subfolder)
% 
% save(fullfile(subfolder, 'BDresults'), 'BDprctile', 'posteriorProbMatrix', 'replayScore', 'begPosition', 'endPosition', 'slope') % , 'pval' , 'replayScore_null'

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
    
