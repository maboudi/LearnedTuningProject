function [posteriorProbMatrix, replayScore, weightedCorr, BDseqscore] = BDreplayDetect_PFshuffle(eventsBinnedfiring, tuningRL, tuningLR, tuningRLshuffles, tuningLRshuffles, activeUnits, binDur)



% This function is intended to calculate the bayesian replay score for each
% event. The significance of the replay scores is evaluated through
% comparing it with a distribution of scores resulted from shuffle
% datasets, i.e. column cycle shuffle, time bin shuffle, etc. 


%%% for now we skip the replay type-we get back to it later ..


% first replace all of zeros in tuning with a small value (at the level of smallest non-zero)

% % tuningRL(~tuningRL) = 1e-4;
% % tuningLR(~tuningLR) = 1e-4;


nEvents       = size(eventsBinnedfiring, 1);
nPositionBins = size(tuningRL, 2);


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



%%% replay score by Bayesian (goodness of fit of the best fitting line)

replayScore  = zeros(nEvents, 1);
weightedCorr = zeros(nEvents, 1);


BDseqscore  = struct('weightedCorr', [], 'replayScore', []);
BDseqscore.weightedCorr = struct('prctilescore', zeros(nEvents, 1), 'zscore', zeros(nEvents, 1));
BDseqscore.replayScore  = struct('prctilescore', zeros(nEvents, 1), 'zscore', zeros(nEvents, 1));


% %%% replay type
% replaydirection = zeros(nEvents,1); %% if above zero, current event is replaying rightward (LR) trajectory template and if less than zero, leftward (RL) template
% replayOrder = zeros(nEvents, 1); %% if above zero is a forward replay and if less than zero is a reverse replay
% 


fprintf('\n Calculating replay scores of actuall PBEs')



%fitLineDictionary = generateLineDictionary(nPositionBins, 1:30);

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
    
    posteriorProbMatrix{pbe, 1} = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [nPositionBins, 1]);
    posteriorProbMatrix{pbe, 1}(:, idx) = 0;

    %%% calculate the best fitting line and the replay score and replay
    %%% type for the current event
                       
%     [replayScore(pbe, 3), slope(pbe, 3), rho(pbe, 3), begPosition(pbe,:,3), endPosition(pbe,:,3)] = replayevaluation(posteriorProbMatrix{pbe, 3});
    
    fitLineDictionary = generateLineDictionary(nPositionBins, nTimebins);
    [replayScore(pbe, 1), weightedCorr(pbe, 1)] = replayevaluation_polar2(posteriorProbMatrix{pbe, 1},fitLineDictionary);

    
    
    %%% calculating replay order(Forward or reverse) and replay direction
    %%% (LR or RL)    
    
%     [replaydirection(pbe), replayOrder(pbe)] = replaytype(postprobRL, postprobLR, slope(pbe, 3), rho(pbe, 3)); 
    
    
    %%% normalize posterior probabilty matrices of each direction in case we are interested in
    %%% doing replay evaluation for each direction independently (as used by Grosmark&Buzsaki, Science2015)
    %%%%%%
%     posteriorProbMatrix{pbe, 1} = postprobRL ./ repmat(sum(postprobRL, 1), [nPositionBins, 1]); %% normalized posterior probabilty matrix for each direction
%     posteriorProbMatrix{pbe, 2} = postprobLR ./ repmat(sum(postprobLR, 1), [nPositionBins, 1]); 
% 
%   
%     posteriorProbMatrix{pbe, 1}(:, idx) = 0; 
%     posteriorProbMatrix{pbe, 2}(:, idx) = 0;
    %%%%%%
    
    % replay evaluation
    
%     
%     
%     [replayScore(pbe, 1), weightedCorr, pbeJumpDistRL, begPosition(pbe,:, 1), endPosition(pbe,:, 1)] = replayevaluation_polar(posteriorProbMatrix{pbe, 1});
%     
%     [replayScore(pbe, 2), weightedCorr, pbeJumpDistLR, begPosition(pbe,:, 2), endPosition(pbe,:, 2)] = replayevaluation_polar(posteriorProbMatrix{pbe, 2});
    
%     
%    fitLineDictionary = generateLineDictionary(nPositionBins, nTimebins);
%     
%     [replayScore(pbe, 1), weightedCorr(pbe, 1), pbeJumpDistRL, onlyLineElements{pbe, 1}, begPosition(pbe,:, 1), endPosition(pbe,:, 1)] = replayevaluation_polar2(posteriorProbMatrix{pbe, 1},fitLineDictionary);
%     
%     [replayScore(pbe, 2), weightedCorr(pbe, 2), pbeJumpDistLR, onlyLineElements{pbe, 2}, begPosition(pbe,:, 2), endPosition(pbe,:, 2)] = replayevaluation_polar2(posteriorProbMatrix{pbe, 2},fitLineDictionary);
    
    
    % concatenate jump distance for the two directions
    
%     
%     f = fieldnames(jumpDist);
%     for ii = 1: length(f)
% %         jumpDist(pbe).(f{ii}) = [pbeJumpDistRL(ii) pbeJumpDistLR(ii)];
%           jumpDist(pbe).(f{ii}) = [pbeJumpDistRL(ii)];
%     end
%     
    if mod(pbe, 50) == 0
        fprintf('.')
    end

end




%% Scores of place field unit ID shuffle


numPFshuffles = size(tuningRLshuffles, 3);

[replayScore_PFshuffle, weightedCorr_PFshuffle] = deal(zeros(nEvents, numPFshuffles));

for pbe = 1: nEvents
    
    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); %%% using just the active units
    nTimebins = size(currEvent, 2);
    % exclude the time bins with no unit firing
    
    sumFiring = sum(currEvent, 1);
    idx       = find(sumFiring == 0);
    
    
    [replayScore_PFshuffle_pbe, weightedCorr_PFshuffle_pbe] = deal(zeros(numPFshuffles, 1));
    
    parfor ss = 1:numPFshuffles % shuffle index

        %%% calculate the posterior probabilty matrix without any normalization 

        postprobRL = baysDecoder(currEvent, tuningRLshuffles(activeUnits, :, ss), binDur); 
        postprobLR = baysDecoder(currEvent, tuningLRshuffles(activeUnits, :, ss), binDur);

        posteriorProbMatrix_PFshuffle = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [nPositionBins, 1]);
        posteriorProbMatrix_PFshuffle(:, idx) = 0;

       fitLineDictionary = generateLineDictionary(nPositionBins, nTimebins);
        [replayScore_PFshuffle_pbe(ss), weightedCorr_PFshuffle_pbe(ss)] = replayevaluation_polar2(posteriorProbMatrix_PFshuffle,fitLineDictionary);

    end
    
    
    replayScore_PFshuffle(pbe, :)  = replayScore_PFshuffle_pbe;
    weightedCorr_PFshuffle(pbe, :) = weightedCorr_PFshuffle_pbe;
    

    BDseqscore.replayScore.prctilescore(pbe)  = length(find(replayScore_PFshuffle_pbe < replayScore(pbe)))/numPFshuffles * 100;
    BDseqscore.replayScore.zscore(pbe)        = (replayScore(pbe) - mean(replayScore_PFshuffle_pbe))/std(replayScore_PFshuffle_pbe);

    BDseqscore.weightedCorr.prctilescore(pbe) = length(find(abs(weightedCorr_PFshuffle_pbe) < abs(weightedCorr(pbe))))/numPFshuffles * 100;
    BDseqscore.weightedCorr.zscore(pbe)       = (abs(weightedCorr(pbe)) - mean(abs(weightedCorr_PFshuffle_pbe))) / std(abs(weightedCorr_PFshuffle_pbe));

    
    if mod(pbe, 50) == 0
       fprintf('\n event #%d: pscore rs %.1f wc %.1f  and zscore rs %.2f wc %.2f', ...
           pbe, BDseqscore.replayScore.prctilescore(pbe), BDseqscore.weightedCorr.prctilescore(pbe), BDseqscore.replayScore.zscore(pbe), BDseqscore.weightedCorr.zscore(pbe))
    end
    

end


end