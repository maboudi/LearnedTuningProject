function [transmats, lambdas, testEvts] = trainCVmodels(data, numofStates, nFolds)

% [transmat, lambda, testEvts] = trainCVmodels(data, tuningRL, tuningLR, positionBins, numofStates, nFolds, binDur, fileinfo, longORshort)

%% train HMMs using k-1 folds of data (K models totally)

nEvents = size(data, 1);
nActiveUnits = size(data{1,1}, 1);

foldSize = floor(nEvents/nFolds); %% the number of events within each fold

wholeEvts = 1 : nEvents;

transmats = zeros(numofStates, numofStates, nFolds);
lambdas = zeros(nActiveUnits, numofStates, nFolds);

testEvts = cell(1, nFolds);
trainEvts = cell(1, nFolds);

for fold = 1 : nFolds
    
    fold
    %%% assiginging the events to train and test sets
    
    if fold ~= nFolds
       testEvts{fold} = wholeEvts((fold-1)*foldSize+1 : fold*foldSize); 
    else
       testEvts{fold} = wholeEvts((fold-1)*foldSize+1 : nEvents); %% in case nFolds is not an divisor of nEvents
    end
       
    trainEvts{fold} = wholeEvts(~ismember(wholeEvts, testEvts{fold})); 
    

    
    %%% making the training dataset
    
    trainData = cell(1, length(trainEvts{fold}));
    
    ntrainEvts = 0;
    for event = trainEvts{fold}
        ntrainEvts = ntrainEvts + 1;
        trainData{ntrainEvts} = data{event, 2}; 
    end
   
    
    %%% training HMM

    % initialize the model parameters
    prior0 = normalise(rand(numofStates,1));
    transmat0 = mk_stochastic(rand(numofStates,numofStates));
    lambda0 = rand(nActiveUnits, numofStates);
    
    

    [~, curr_prior, curr_transmat, curr_lambda] = phmm_em(trainData, prior0, transmat0, lambda0, 'max_iter', 200, 'verbose', 1);
    
%     [~, startIdx] = max(curr_prior);
    transmats(:, :, fold) = curr_transmat;
    lambdas(:, :, fold)   = curr_lambda;
%     [transmats(:, :, fold), lambdas(:, :, fold)] = sortStates(curr_transmat, curr_lambda, curr_prior);
    
%     transmats(:, :, fold) = curr_transmat(sortIdx, sortIdx);
%     lambdas(:, :, fold) = curr_lambda(:, sortIdx);
    
end

% 
% currDir = pwd;
% FileBase = [currDir '/' fileinfo.name  '/HMM/sequenceDetection/cross-validation/' period];
% mkdir(FileBase)
% 
% save([FileBase '/crossValidationHMMs.mat'], 'transmats', 'lambdas', 'testEvts');

end

