function [BDprctile, PosteriorProbMatrix, begPosition, endPosition, slope, replaydirection, replayOrder] = BDreplayDetect_grosmark(eventsBinnedfiring, tuning, activeUnits, binDur, FileBase, dataType)

% This function is intended to calculate the bayesian replay score for each
% event. The significance of the replay scores is evaluated through
% comparing it with a distribution of scores resulted from shuffle
% datasets, i.e. column cycle shuffle, time bin shuffle, etc. 



noEvents = size(eventsBinnedfiring, 1);
noPositionBins = size(tuning, 2);

noshufflingMethods = 1; %% number of methods for shuffling (time swap and column cycle shuffle here)
noShuffle = 500; %% number of shuffles in each of the methods 


if isempty(activeUnits)
   activeUnits = 1:size(tuning, 1); 
end

%%%%%%%%%%%% Initiation %%%%%%%%%%%%%


%%% posterior probability matrices
%%% the matrix comprised of posterior probability of the track positions (as rows) within each time bins of the event (as columns)


PosteriorProbMatrix = cell(noEvents,1); 

%%% 1st array -->  leftward(RL) trajectory only
%%% 2nd array --> rightward(LR) trajectory only
%%% 3rd array -->  summed over the two directions
%%% we can use either the sum or the two trajectories separately for replay
%%% detection (here we have used just the sum)



%%% best fitting line
slope = zeros(noEvents, 1); %% RL, LR, and both
rho = zeros(noEvents,1);  
begPosition = zeros(noEvents,2);   %% y intercept; the third dimention is comprised of RL,LR, and both, respectively
endPosition = zeros(noEvents,2);  

%%% replay score by Bayesian (goodness of fit of the best fitting line)
replayScore = zeros(noEvents,1);


%%% replay type
replaydirection = zeros(noEvents,1); %% if above zero, current event is replaying rightward (LR) trajectory template and if less than zero, leftward (RL) template
replayOrder = zeros(noEvents, 1); %% if above zero is a forward replay and if less than zero is a reverse replay

%%% the matrix for shuffled data 
% 
% replayScore_null = cell(1,3);
% 
% for ii = 1:3
%     replayScore_null{ii} = zeros(noEvents, noShuffle, noshufflingMethods);
% end
% 

%%% significance evaluation
BDprctile = zeros(noEvents, noshufflingMethods); 
% pval = zeros(noEvents, noshufflingMethods, 3);


for evt = 1: noEvents
    
    currEvent = eventsBinnedfiring{evt, 2}(activeUnits, :); %%% using just the active units
    
    
    % exclude the time bins with no unit firing
    
    sumFiring = sum(currEvent, 1);
    idx = find(sumFiring == 0);
    
    %%% Bayesian Sequential Score

    %%% calculate the posterior probabilty matrix without any normalization 
    
    postprob = baysDecoder(currEvent, tuning(activeUnits, :), binDur); 
    

    
    postprob = postprob + eps*eps;
    
    
    
    %%% normalize posterior probabilty matrices of each direction in case we are interested in
    %%% doing replay evaluation for each direction independently (as used by Grosmark and Buzsaki 2015 science paper)
        
    PosteriorProbMatrix{evt} = postprob ./ repmat(sum(postprob, 1), [noPositionBins, 1]); %% normalized posterior probabilty matrix for each direction

  
    PosteriorProbMatrix{evt}(:, idx) = eps*eps; 
    
    
    [replayScore(evt), slope(evt), rho(evt), begPosition(evt,:), endPosition(evt,:)] = replayevaluation(PosteriorProbMatrix{evt}, 1);
    
    
end





%%% calculating the replay score for the shuffled data (with two
%%% shuffling methods)

for evt = 1 : noEvents
    
    
        method = 1;
        

        replayScore_null = zeros(noShuffle, 1);
        
        currPosteriorProbMatrix = PosteriorProbMatrix(evt);
        
        parfor sn = 1 : noShuffle 
            
%             if method == 2 %%% column cycle shuffle
% 
% %                 % generating random position shift for each time bin
% %                 
% %                 randShift = randi(noPositionBins-1, 1, notimeBins);
%                 
% %                 PosteriorProbMatrix_shuffle = column_cycle_shuffle(PosteriorProbMatrix{evt, 3}, randShift(sn, :));
%                 
%                 %%% each direction separately
%                 PosteriorProbMatrixRL_shuffle = column_cycle_shuffle(PosteriorProbMatrix{evt, 1}, randShift(sn, :));
%                 PosteriorProbMatrixLR_shuffle = column_cycle_shuffle(PosteriorProbMatrix{evt, 2}, randShift(sn, :));
            
%             elseif method == 1
                
%                 PosteriorProbMatrix_shuffle = timebinshuffle(PosteriorProbMatrix(:,3), evt);
%                 PosteriorProbMatrixRL_shuffle = timebinshuffle(PosteriorProbMatrix(:,1), evt);
%                 PosteriorProbMatrixLR_shuffle = timebinshuffle(PosteriorProbMatrix(:,2), evt);
%            
%                 PosteriorProbMatrix_shuffle = genTimeSwap(PosteriorProbMatrix(evt, 3));
                PosteriorProbMatrix_shuffle = genTimeSwap(currPosteriorProbMatrix);
                PosteriorProbMatrix_shuffle = PosteriorProbMatrix_shuffle{1};
                
                
%             end
            
            
%             replayScore_null{3}(evt,sn,method) = replayevaluation(PosteriorProbMatrix_shuffle); 
            

            %%% each direction separately
%             replayScore_null{1}(evt,sn,method) = replayevaluation(PosteriorProbMatrixRL_shuffle, 1);
%             replayScore_null{2}(evt,sn,method) = replayevaluation(PosteriorProbMatrixLR_shuffle, 1);
              replayScore_null(sn) = replayevaluation(PosteriorProbMatrix_shuffle, 1);
            
        end

        %%% replayScore(z) or zrScore is defined as how distanct the
        %%% replay score of raw data is from median of the shuffle
        %%% distribution normalized to the distance of 95 percentile of the
        %%% distribution to the median. We do the same for each direction independently as well
        
        BDprctile(evt, method) = length(find(replayScore_null < replayScore(evt)))/noShuffle;
        
%         for ii = 1:2 %% direction
%             
%             BDprctile(evt, method, ii) = length(find(replayScore_null{ii}(evt,:,method) < replayScore(evt, ii)))/noShuffle;
%             
% %             [BDprctile(evt, method, ii), pval(evt, method, ii)] = rawScoreCompare2null(replayScore(evt, ii), replayScore_null{ii}(evt,:,method)); 
%         end
%         
        
%     end

    fprintf('replay scores of event #%d are %.2f \n', evt, BDprctile(evt,1))
    
end


subfolder = fullfile(FileBase, dataType);
mkdir(subfolder)

save(fullfile(subfolder, 'BDresults'), 'BDprctile', 'PosteriorProbMatrix', 'replayScore', 'begPosition', 'endPosition', 'slope') % , 'pval' , 'replayScore_null'

end

