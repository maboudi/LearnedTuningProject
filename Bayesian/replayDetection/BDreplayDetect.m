function [PosteriorProbMatrix, BDseqscore, begPosition, endPosition, slope, replaydirection, replayOrder] = BDreplayDetect(eventsBinnedfiring, tuningRL, tuningLR, activeUnits, binDur, FileBase, dataType)

% BDseqscore instead of BDprctile

% This function is intended to calculate the bayesian replay score for each
% event. The significance of the replay scores is evaluated through
% comparing it with a distribution of scores resulted from shuffle
% datasets, i.e. column cycle shuffle, time bin shuffle, etc. 


%%% for now we skip the replay type-we get back to it later ..


% first replace all of zeros in tuning with a small value (at the level of smallest non-zero)
tuningRL(~tuningRL) = 1e-4;
tuningLR(~tuningLR) = 1e-4;


noEvents       = size(eventsBinnedfiring, 1);
noPositionBins = size(tuningRL, 2);

% noshufflingMethods = 2; %% number of methods for shuffling (time swap and column cycle shuffle here)
nShuffles = 500; %% number of shuffles in each of the methods 


if isempty(activeUnits)
   activeUnits = 1:size(tuningRL, 1); 
end

%%%%%%%%%%%% Initiation %%%%%%%%%%%%%


%%% posterior probability matrices
%%% the matrix comprised of posterior probability of the track positions (as rows) within each time bins of the event (as columns)


PosteriorProbMatrix = cell(noEvents,3); 

%%% 1st array -->  leftward(RL) trajectory only
%%% 2nd array --> rightward(LR) trajectory only
%%% 3rd array -->  summed over the two directions
%%% we can use either the sum or the two trajectories separately for replay
%%% detection (here we have used just the sum)



%%% best fitting line
slope       = zeros(noEvents, 3); %% RL, LR, and both together (?)
rho         = zeros(noEvents,3);  

begPosition = zeros(noEvents,2,3);   %% y intercept; the third dimention is comprised of RL,LR, and both, respectively
endPosition = zeros(noEvents,2,3);  

%%% replay score by Bayesian (goodness of fit of the best fitting line)
maxgoodnessofFit = zeros(noEvents,3);


%%% replay type
replaydirection = zeros(noEvents,1); %% if above zero, current event is replaying rightward (LR) trajectory template and if less than zero, leftward (RL) template
replayOrder = zeros(noEvents, 1); %% if above zero is a forward replay and if less than zero is a reverse replay



%%% significance evaluation
BDseqscore = struct('prctilescore', zeros(nEvents, 2), 'zscore', zeros(nEvents, 2)); % significance: data as a percentile of shuffle distribution or zscore


for pbe =  1: noEvents
    
    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); %%% using just the active units
    
    % exclude the time bins with no unit firing
    
    sumFiring = sum(currEvent, 1);
    idx       = find(sumFiring == 0);
    
    %%% Bayesian Sequential Score

    %%% calculate the posterior probabilty matrix without any normalization 
    
    postprobRL = baysDecoder(currEvent, tuningRL(activeUnits, :), binDur); 
    postprobLR = baysDecoder(currEvent, tuningLR(activeUnits, :), binDur);

    
    postprobRL = postprobRL + eps*eps; % what was this step exactly??
    postprobLR = postprobLR + eps*eps;
    
    
    %%% Summed probability distribution, every column is normalized to one
    %%% (this is the same method as used in Davidson et al. 2009)
    
%     PosteriorProbMatrix{pbe, 3} = (postprobRL + postprobLR)./repmat(sum(postprobRL + postprobLR, 1), [noPositionBins, 1]);
%     PosteriorProbMatrix{pbe, 3}(:, idx) = 0;

    %%% calculate the best fitting line and the replay score and replay
    %%% type for the current event
                       
%     [maxgoodnessofFit(pbe, 3), slope(pbe, 3), rho(pbe, 3), begPosition(pbe,:,3), endPosition(pbe,:,3)] = replayevaluation(PosteriorProbMatrix{pbe, 3});
    


    %%% calculating replay order(Forward or reverse) and replay direction
    %%% (LR or RL)    
    
%     [replaydirection(pbe), replayOrder(pbe)] = replaytype(postprobRL, postprobLR, slope(pbe, 3), rho(pbe, 3)); 
    
    
    %%% normalize posterior probabilty matrices of each direction in case we are interested in
    %%% doing replay evaluation for each direction independently (as used by Grosmark&Buzsaki, Science2015)
        
    PosteriorProbMatrix{pbe, 1} = postprobRL ./ repmat(sum(postprobRL, 1), [noPositionBins, 1]); %% normalized posterior probabilty matrix for each direction
    PosteriorProbMatrix{pbe, 2} = postprobLR ./ repmat(sum(postprobLR, 1), [noPositionBins, 1]); 

  
    PosteriorProbMatrix{pbe, 1}(:, idx) = eps*eps; 
    PosteriorProbMatrix{pbe, 2}(:, idx) = eps*eps;
    
    
    [maxgoodnessofFit(pbe, 1), pbeJumpDistRL, slope(pbe, 1), rho(pbe, 1), begPosition(pbe,:, 1), endPosition(pbe,:, 1)] = replayevaluation_polar(PosteriorProbMatrix{pbe, 1});
    [maxgoodnessofFit(pbe, 2), pbeJumpDistLR, slope(pbe, 2), rho(pbe, 2), begPosition(pbe,:, 2), endPosition(pbe,:, 2)] = replayevaluation_polar(PosteriorProbMatrix{pbe, 2});
    
    
    
    % concatenate jump distance for the two directions
    f = fieldnames(jumpDist);
    for ii = 1: length(f)
        jumpDist(pbe).(f{ii}) = [pbeJumpDistRL(ii) pbeJumpDistLR(ii)];
    end
    
end


%%% calculating the replay score for the shuffled data (with two
%%% shuffling methods)

for pbe = 1: noEvents
    
    
    currEvent = eventsBinnedfiring{pbe, 2}(activeUnits, :); 
    notimeBins = size(currEvent, 2);
    

    maxgoodnessofFit_nullRL = zeros(nShuffles, 1);
    maxgoodnessofFit_nullLR = zeros(nShuffles, 1);

    PosteriorProbMatrixRL = PosteriorProbMatrix(pbe, 1);
    PosteriorProbMatrixLR = PosteriorProbMatrix(pbe, 2);
    
    
    % method one: column cycle shuffle
    parfor sn = 1 : nShuffles 
        
        randShift = randi(noPositionBins-1, 1, notimeBins);

        PosteriorProbMatrixRL_shuffle = column_cycle_shuffle(PosteriorProbMatrixRL, randShift);
        PosteriorProbMatrixLR_shuffle = column_cycle_shuffle(PosteriorProbMatrixLR, randShift);        
        
        maxgoodnessofFit_nullRL(sn) = replayevaluation_polar(PosteriorProbMatrixRL_shuffle, 1);
        maxgoodnessofFit_nullLR(sn) = replayevaluation_polar(PosteriorProbMatrixLR_shuffle, 1);
        
    end
    
     
%     % method two: time swap
%     parfor sn = 1: nShuffles
%          
% 
%         PosteriorProbMatrixRL_shuffle = genTimeSwap(PosteriorProbMatrixRL);
%         PosteriorProbMatrixRL_shuffle = PosteriorProbMatrixRL_shuffle{1};
% 
%         PosteriorProbMatrixLR_shuffle = genTimeSwap(PosteriorProbMatrixLR);
%         PosteriorProbMatrixLR_shuffle = PosteriorProbMatrixLR_shuffle{1};
% 
%         maxgoodnessofFit_nullRL(sn) = replayevaluation(PosteriorProbMatrixRL_shuffle, 1);
%         maxgoodnessofFit_nullLR(sn) = replayevaluation(PosteriorProbMatrixLR_shuffle, 1);
% 
% 
%     end



    % sequence score for each PBE calcualted as a percentile of the shuffle distribution
    
    for template = 1:2
    
        BDseqscore.prctilescore(pbe, dir) = length(find(maxgoodnessofFit_nullRL < maxgoodnessofFit(pbe, 1)))/nShuffles*100;
        
        BDseqscore.zscore(pbe, dir)       = (maxgoodnessofFit(pbe, 1) - mean(maxgoodnessofFit_nullRL))/std(maxgoodnessofFit_nullRL);
        
    end
    

    if mod(pbe, 100) == 0
       fprintf('replay scores of event #%d are pscore %.1f  and zscore %.2f \n', pbe, max(BDseqscore.prctilescore(pbe, :)), max(BDseqscore.zscore(pbe, :)))
    end
    
end
% 
% 
% subfolder = fullfile(FileBase, dataType);
% mkdir(subfolder)
% 
% save(fullfile(subfolder, 'BDresults'), 'BDprctile', 'PosteriorProbMatrix', 'maxgoodnessofFit', 'begPosition', 'endPosition', 'slope') % , 'pval' , 'maxgoodnessofFit_null'

end

