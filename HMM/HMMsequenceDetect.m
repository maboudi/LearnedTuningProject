function  [HMMScoreZ, HMMScoreZ_poisson, HMMScoreZ_tisw] = ...
                            HMMsequenceDetect(data, tuningRL, tuningLR, positionBins, numofStates, numofFolds, binDur, fileinfo, longORshort)

% compute_poisson = 0;
% if nargout > 2
%     compute_poisson = 1;
% end


%% train HMMs using k-1 folds of data (K models totally)

noEvents = size(data, 1);
noActiveUnits = size(data{1,1}, 1);

foldSize = floor(noEvents/numofFolds); %% the number of events within each fold

wholeEvts = 1 : noEvents;

transmat = zeros(numofStates, numofStates, numofFolds);
lambda = zeros(noActiveUnits, numofStates, numofFolds);

testEvts = cell(1, numofFolds);
trainEvts = cell(1, numofFolds);

for fold = 1 : numofFolds
    
    
    %%% assiginging the events to train and test sets
    
    if fold ~= numofFolds
       testEvts{fold} = wholeEvts((fold-1)*foldSize+1 : fold*foldSize); 
    else
       testEvts{fold} = wholeEvts((fold-1)*foldSize+1 : noEvents); %% in case numofFolds is not an divisor of noEvents
    end
       
    trainEvts{fold} = wholeEvts(~ismember(wholeEvts, testEvts{fold})); 
    

    
    %%% making the training dataset
    
    trainData = cell(1, length(trainEvts{fold}));
    
    notrainEvts = 0;
    for event = trainEvts{fold}
        notrainEvts = notrainEvts + 1;
        trainData{notrainEvts} = data{event, 2}; 
    end
   
    
    %%% training HMM

    % initialize the model parameters
    prior0 = normalise(rand(numofStates,1));
    transmat0 = mk_stochastic(rand(numofStates,numofStates));
    lambda0 = rand(noActiveUnits, numofStates);
    
    

    [LL, prior, transmat(:, :, fold), lambda(:, :, fold)] = phmm_em(trainData, prior0, transmat0, lambda0, 'max_iter', 20);
    
    
    
    %% visulization of the model parameters
    
    currDir = pwd;
    FileBase = [currDir '/' fileinfo.name '/data part' num2str(longORshort) '/HMM/sequenceDetection/trainedModels/fold #' num2str(fold)];
    mkdir(FileBase)

    HMMParametersVisualiztion(lambda(:, :, fold), transmat(:, :, fold), tuningRL, tuningLR, binDur, positionBins, FileBase);
    
    
    fold_transmat = transmat(:, :, fold);
    fold_lambda = lambda(:, :, fold);
    save([FileBase '/foldoutputs.mat'], 'fold_transmat', 'fold_lambda');
    
    close all

    
end



%% recruiting the HMM model on test data (the remaining fold)

%%% we do sequence detection in raw and shuffled (or simulated) dataset to
%%% make a comparison (Like as we did with BD replay detection, we aim to
%%% make a comparison between raw and shuffle data to see how much the
%%% number of detected sequences in the two datasets differs from each other)


noShuffle = 1000;

currDir = pwd;
FileBase = [currDir '/' fileinfo.name '/data part' num2str(longORshort) '/HMM/sequenceDetection/measure congruence/'];
mkdir(FileBase)

%%% raw Data

testDataType = 'raw';
[HMMScoreZ, maxCongruenceLength] = measureHMMcongruence(data, testEvts, transmat, lambda, noShuffle, testDataType, FileBase, binDur);


%%% poisson simulated data

% if compute_poisson
testDataType = 'Poisson'; 
[HMMScoreZ_poisson, maxCongruenceLength_poisson] = measureHMMcongruence(data, testEvts, transmat, lambda, noShuffle, testDataType, FileBase, binDur);
% end



end
