
function sequenceScores = sequenceScorevstime_CV(data, pbeEpochInd, numofStates, fileinfo)



noEpochs = max(pbeEpochInd);

noEvents = size(data, 1);
wholeEvts = 1 : noEvents;

noActiveUnits = size(data{1,1}, 1);

sequenceScores = cell(noEpochs, 1);
% PBELikelihood = cell(noEpochs, 1);

% initialize model parameters

noShuffle = 500;

prior0 = normalise(rand(numofStates,1));
transmat0 = mk_stochastic(rand(numofStates,numofStates));
lambda0 = rand(noActiveUnits, numofStates);


for epoch = 1: noEpochs
    
    
    % train the model on a single epoch
    
    trainPBEset = find(pbeEpochInd == epoch);
    testPBEset = setdiff(wholeEvts, trainPBEset);
    
    trainData = data(trainPBEset, 2);
    testData = data(testPBEset, :);
    
    [~, ~, curr_transmat, curr_lambda] = phmm_em(trainData, prior0, transmat0, lambda0, 'max_iter', 200, 'verbose', 1);
    
    
    
    
    % calcualte the likelihood of the test events (the train events are excluded)
    
%     PBELikelihood{epoch} = nan(noEvents, 1);
%     
%     for evt = 1: length(testData)
%         
%         currEvent = testData{evt};
%         B = poisson_prob(currEvent, curr_lambda,1);
%         prior = 1/numofStates * ones(numofStates,1); %%% a uniform prior probability distribution over the states
% 
%         [~, ~, ~,  PBELikelihood{epoch}(testPBEset(evt)), ~] = fwdback(prior, curr_transmat, B);
%     end


    sequenceScores{epoch} = nan(noEvents, 1);
    
    sequenceScores{epoch}(testPBEset) = modelCongruence(testData, curr_transmat, curr_lambda, noShuffle, fileinfo, [], [], 0);
    
end


end

