function [trainLikelihood, CVLikelihood] = HMMaRUNquality(runBouts, spikeStruct, qclus, binDur, nStates, fileinfo, ifshuffle)


% time-binning the events (here the theta bouts/high speed bouts instead of PBEs)

runBoutsBinnedfiring  = timeBinning(runBouts, spikeStruct, qclus, binDur, fileinfo);
runBoutsBinnedfiring  = removeSideSilentBins(runBoutsBinnedfiring, runBouts(:, 1), binDur, fileinfo);


% randomizing the temporal order of the events
rndgen   = randperm(size(runBoutsBinnedfiring, 1)); 
data     = runBoutsBinnedfiring(rndgen, :); 


if ifshuffle
    
    data = timeswap(data, binDur);
end


nActiveUnits = size(data{1,2}, 1);



% training likelihood

[transmat0, lambda0, prior0] = initialModel(nStates, nActiveUnits);
trainLikelihood              = phmm_em(runBoutsBinnedfiring(:,2), prior0, transmat0, lambda0, 'verbose', 1, 'max_iter', 200);



%% cross-validation


% train HMMs using k-1 folds of data 


nEvents      = size(data, 1);
nFolds       = 5;


foldSize     = floor(nEvents/nFolds); %% the number of events within each fold

wholeEvts    = 1:nEvents;

CVLikelihood = zeros(1, nEvents); %% log likelihood of the raw data


testEvts  = cell(1, nFolds);
trainEvts = cell(1, nFolds);


for fold = 1 : nFolds
    
    %%% assiginging the events to train and test sets
    
    if fold ~= nFolds
       testEvts{fold} = wholeEvts((fold-1)*foldSize+1 : fold*foldSize); 
    else
       testEvts{fold} = wholeEvts((fold-1)*foldSize+1 : nEvents); %% in case nFolds is not an divisor of nEvents
    end
       
    trainEvts{fold}   = wholeEvts(~ismember(wholeEvts, testEvts{fold})); 
    

    
    %%% building the training dataset
    
    trainData  = cell(1, length(trainEvts{fold}));
    
    ntrainEvts = 0;
    for event = trainEvts{fold}
        ntrainEvts            = ntrainEvts + 1;
        trainData{ntrainEvts} = data{event, 2}; 
    end
   
    
    %%% training HMM

    % initialize the model parameters
    prior0    = normalise(rand(nStates,1));
    transmat0 = mk_stochastic(rand(nStates,nStates));
    lambda0   = rand(nActiveUnits, nStates);
    
    
    [~, ~, curr_transmat, curr_lambda] = phmm_em(trainData, prior0, transmat0, lambda0, 'max_iter', 300, 'verbose', 0);
    
    
    for evt = testEvts{fold}
        
        currEvent = data{evt, 2};
        
        B         = poisson_prob(currEvent, curr_lambda,1);
        prior     = 1/nStates * ones(nStates,1); %%% a uniform prior probability distribution over the states
        
        %%% observation likelihood for the raw data
        [~, ~, ~,  CVLikelihood(evt), ~] = fwdback(prior, curr_transmat, B);
    
    end
    

end




end